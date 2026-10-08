"""Process-isolated access to the parameterized COMSOL helium/SS304 pipe.

No Java VM is started in the web server. One worker owns one COMSOL client and
reuses its loaded template across calls. All public inputs and outputs use SI.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import queue
import shutil
import subprocess
import sys
import threading
import time
from typing import Callable

from models.config import PROJECT_ROOT as WORKSPACE, DATA_DIR, COMSOL_TEMPLATE as TEMPLATE, COMSOL_EXECUTABLE
PROTOCOL = "FLUID_EDA_JSON "
SCHEMA_VERSION = 1


class ComsolError(RuntimeError):
    """COMSOL startup, numerical solution, or protocol failed; no fallback used."""


class ComsolCancelled(ComsolError):
    pass


def normalize_parameters(params, p_a, p_b, t_a, t_b):
    """Validate geometry, material temperatures and positive absolute pressures."""
    p = dict(params or {})

    def first(*keys, default):
        return next((p[key] for key in keys if key in p), default)

    data = {
        "length_m": first("length_m", "L", default=0.288),
        "outer_diameter_m": first("outer_diameter_m", "Do", default=0.00635),
        "wall_thickness_m": first("wall_thickness_m", "wall_w", default=0.000762),
        "pressure_a": p_a, "pressure_b": p_b,
        "temperature_a": t_a, "temperature_b": t_b,
        "axial_elements": first("axial_elements", default=320),
        "radial_elements": first("radial_elements", default=16),
        "profile_points": first("profile_points", default=31),
    }
    for name, value in data.items():
        if isinstance(value, bool):
            raise ValueError(f"{name} must be a number, not a boolean.")
        try:
            data[name] = float(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{name} must be a finite number.") from exc
        if not math.isfinite(data[name]):
            raise ValueError(f"{name} must be finite.")
    for name in ("length_m", "outer_diameter_m", "wall_thickness_m"):
        if data[name] <= 0:
            raise ValueError(f"{name} must be positive.")
    if 2 * data["wall_thickness_m"] >= data["outer_diameter_m"]:
        raise ValueError("Wall thickness must be less than half the outer diameter.")
    for name in ("pressure_a", "pressure_b"):
        if data[name] <= 0:
            raise ValueError(f"{name} must be greater than zero Pa absolute.")
    for name in ("temperature_a", "temperature_b"):
        if data[name] <= 0:
            raise ValueError(f"{name} must be greater than zero K absolute.")
    for name, low, high in (("axial_elements", 40, 1280),
                            ("radial_elements", 4, 64),
                            ("profile_points", 3, 501)):
        if not data[name].is_integer() or not low <= data[name] <= high:
            raise ValueError(f"{name} must be an integer in [{low}, {high}].")
        data[name] = int(data[name])
    boundary = p.get("wall_boundary", "adiabatic")
    if boundary not in ("adiabatic", "prescribed"):
        raise ValueError("wall_boundary must be 'adiabatic' or 'prescribed'.")
    data["wall_boundary"] = boundary
    data["wall_profile"] = []
    if boundary == "prescribed":
        sections = p.get("wall_profile")
        if not isinstance(sections, list) or not 1 <= len(sections) <= 100:
            raise ValueError("A prescribed wall needs 1–100 wall_profile sections.")
        previous = data["temperature_a"]
        for section in sections:
            if not isinstance(section, dict):
                raise ValueError("Each wall_profile section must be an object.")
            try:
                row = {key: float(section[key]) for key in ("length_m", "T_start_K", "T_end_K")}
            except (KeyError, ValueError, TypeError) as exc:
                raise ValueError("A wall section requires length_m, T_start_K, and T_end_K.") from exc
            if not all(math.isfinite(value) for value in row.values()) or row["length_m"] <= 0:
                raise ValueError("Wall section length must be positive and values finite.")
            if not all(row[key] > 0 for key in ("T_start_K", "T_end_K")):
                raise ValueError("Wall section temperatures must be greater than zero K absolute.")
            if not math.isclose(previous, row["T_start_K"], rel_tol=0, abs_tol=1e-6):
                raise ValueError("The wall profile must be continuous and start at temperature_a.")
            previous = row["T_end_K"]
            data["wall_profile"].append(row)
        if not math.isclose(previous, data["temperature_b"], rel_tol=0, abs_tol=1e-6):
            raise ValueError("The wall profile must end at temperature_b.")
        if not math.isclose(sum(row["length_m"] for row in data["wall_profile"]), data["length_m"], rel_tol=1e-8, abs_tol=1e-10):
            raise ValueError("The wall profile section lengths must sum to the pipe length.")
    return data


def _sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def model_fingerprint(template=TEMPLATE):
    """Cache identity includes template, adapter, initialization, and properties."""
    paths = [Path(template), Path(__file__), Path(__file__).with_name("comsol_worker.py"),
             WORKSPACE / "models/comsol_pipe/steel_properties.py",
             WORKSPACE / "models/comsol_pipe/helium_properties.py",
             WORKSPACE / "models/one_k_pot/he4_viscosity.py",
             WORKSPACE / "models/one_k_pot/properties.py"]
    return {"schema": SCHEMA_VERSION, "comsol": "6.1",
            "files": {str(path.relative_to(WORKSPACE)) if path.is_relative_to(WORKSPACE)
                      else str(path): _sha256_file(path) for path in paths}}


def cache_key(inputs, fingerprint):
    return hashlib.sha256(json.dumps({"inputs": inputs, "model": fingerprint},
                                    sort_keys=True, separators=(",", ":"),
                                    allow_nan=False).encode()).hexdigest()


class ComsolPipeBackend:
    """Lazy, persistent worker, usable as a context manager for one network job.

    ``progress`` receives plain strings. ``cancel_event`` follows threading.Event.
    A timeout or cancellation terminates the worker, including its in-process JVM.
    Calls on one backend are serialized; independent jobs use separate backends.
    """

    def __init__(self, cache_dir=None, *, template=TEMPLATE, executable=COMSOL_EXECUTABLE,
                 cores=2, timeout=900, progress: Callable[[str], None] | None = None,
                 cancel_event=None):
        self.template = Path(template).resolve()
        self.executable = Path(executable).resolve()
        self.cache_dir = Path(cache_dir or DATA_DIR / "comsol_cache")
        self.cores = max(1, int(cores))
        self.timeout = float(timeout)
        self.progress = progress
        self.cancel_event = cancel_event
        self._process = None
        self._queue = queue.Queue()
        self._lock = threading.Lock()
        self._fingerprint = None
        self._counter = 0
        self._log_tail = []

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()

    def _report(self, message):
        if self.progress:
            self.progress(str(message))

    def _check_cancel(self):
        if self.cancel_event is not None and self.cancel_event.is_set():
            self.close()
            raise ComsolCancelled("COMSOL calculation was cancelled.")

    def _reader(self, process, target_queue):
        for line in iter(process.stdout.readline, ""):
            if PROTOCOL in line:
                try:
                    target_queue.put(json.loads(line.split(PROTOCOL, 1)[1]))
                    continue
                except json.JSONDecodeError:
                    pass
            self._log_tail.append(line.rstrip())
            self._log_tail = self._log_tail[-60:]
        target_queue.put({"type": "eof"})

    def _wait(self, expected_type, request_id, deadline):
        while True:
            self._check_cancel()
            if time.monotonic() > deadline:
                self.close()
                raise ComsolError(f"COMSOL exceeded the {self.timeout:g} s evaluation timeout.")
            try:
                message = self._queue.get(timeout=0.2)
            except queue.Empty:
                continue
            kind = message.get("type")
            if kind == "progress":
                self._report(message.get("message", ""))
            elif kind == "error":
                raise ComsolError(message.get("error", "COMSOL worker failed."))
            elif kind == "eof":
                self.close()
                raise ComsolError("COMSOL worker exited unexpectedly. " + "\n".join(self._log_tail[-12:]))
            elif kind == expected_type and message.get("id") == request_id:
                return message

    def _start(self, deadline):
        if self._process is not None and self._process.poll() is None:
            return
        if not self.template.is_file():
            raise ComsolError(f"COMSOL template does not exist: {self.template}. Run python -m models.comsol_pipe.build_model --build-only first.")
        if not self.executable.is_file():
            raise ComsolError(f"COMSOL installation does not exist: {self.executable}")
        self._report("Starting isolated COMSOL 6.1 worker and loading pipe template…")
        self._queue = queue.Queue()
        command = [sys.executable, "-u", "-m", "models.fluid_eda.comsol_worker",
                   "--template", str(self.template), "--executable", str(self.executable),
                   "--cores", str(self.cores)]
        env = dict(os.environ, PYTHONIOENCODING="utf-8", PYTHONUTF8="1")
        self._process = subprocess.Popen(command, cwd=WORKSPACE, env=env,
                                         stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                         stderr=subprocess.STDOUT, text=True, encoding="utf-8",
                                         errors="replace", bufsize=1,
                                         creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
        threading.Thread(target=self._reader, args=(self._process, self._queue), daemon=True).start()
        self._wait("ready", None, deadline)

    def evaluate_pipe(self, params, p_a, p_b, t_a, t_b, *, save_path=None, profiles=True):
        inputs = normalize_parameters(params, p_a, p_b, t_a, t_b)
        if save_path is not None:
            save_path = Path(save_path).resolve()
            if save_path == self.template or save_path == TEMPLATE.resolve():
                raise ValueError("The original COMSOL template cannot be overwritten.")
            if save_path.suffix.lower() != ".mph":
                raise ValueError("COMSOL artifacts must have the .mph extension.")
        with self._lock:
            self._check_cancel()
            if self._fingerprint is None:
                if not self.template.is_file():
                    raise ComsolError(f"COMSOL template does not exist: {self.template}. Run python -m models.comsol_pipe.build_model --build-only first.")
                self._fingerprint = model_fingerprint(self.template)
            key = cache_key(inputs, self._fingerprint)
            cache_path = self.cache_dir / f"{key}.json"
            if cache_path.is_file():
                try:
                    cached = json.loads(cache_path.read_text(encoding="utf-8"))
                except (OSError, ValueError):
                    cached = {}  # An incomplete cache is recomputed by COMSOL.
                if cached.get("cache_key") == key and cached.get("inputs") == inputs:
                    reusable = save_path is None
                    if save_path is not None and (cached.get("artifact_cache") or cached.get("artifact")):
                        artifact = Path(cached.get("artifact_cache") or cached["artifact"]).resolve()
                        if artifact.is_file():
                            expected_hash = cached.get("artifact_sha256")
                            actual_hash = _sha256_file(artifact)
                            if expected_hash is None or expected_hash == actual_hash:
                                self._check_cancel()
                                try:
                                    save_path.parent.mkdir(parents=True, exist_ok=True)
                                    if artifact != save_path:
                                        shutil.copy2(artifact, save_path)
                                except OSError as exc:
                                    raise ComsolError(f"Cannot copy the cached COMSOL artifact to {save_path}: {exc}") from exc
                                cached["artifact"] = str(save_path)
                                cached["artifact_sha256"] = actual_hash
                                reusable = True
                    if reusable:
                        cached["cache_hit"] = True
                        self._report("COMSOL pipe: reusing exact matching solved result.")
                        return cached
            deadline = time.monotonic() + self.timeout
            self._start(deadline)
            self._counter += 1
            request_id = self._counter
            payload = {"id": request_id, "inputs": inputs,
                       "save_path": str(save_path) if save_path else None}
            try:
                self._process.stdin.write(json.dumps(payload, allow_nan=False) + "\n")
                self._process.stdin.flush()
            except (OSError, BrokenPipeError) as exc:
                self.close()
                raise ComsolError("Could not send the pipe evaluation to COMSOL.") from exc
            result = self._wait("result", request_id, deadline)["result"]
            result.update(cache_key=key, cache_hit=False, inputs=inputs,
                          model_fingerprint=self._fingerprint)
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            if result.get("artifact") and Path(result["artifact"]).is_file():
                result["artifact_sha256"] = _sha256_file(result["artifact"])
                # A network job replaces pipe_1.mph on each correction round.
                # Keep a separate immutable artifact for each exact cache key.
                cached_artifact = (self.cache_dir / f"{key}.mph").resolve()
                if Path(result["artifact"]).resolve() != cached_artifact:
                    artifact_temporary = cached_artifact.with_suffix(f".{os.getpid()}.{threading.get_ident()}.mph.tmp")
                    shutil.copy2(result["artifact"], artifact_temporary)
                    artifact_temporary.replace(cached_artifact)
                result["artifact_cache"] = str(cached_artifact)
            temporary = cache_path.with_suffix(f".{os.getpid()}.{threading.get_ident()}.tmp")
            temporary.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
            temporary.replace(cache_path)
            return result

    def close(self):
        process, self._process = self._process, None
        if process is not None:
            if process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait(timeout=5)
            for stream in (process.stdin, process.stdout):
                if stream:
                    stream.close()


def evaluate_pipe(params, p_a, p_b, t_a, t_b, *, cache_dir=None, save_path=None, profiles=True, **options):
    """Convenience call; use ComsolPipeBackend for repeated network evaluations."""
    with ComsolPipeBackend(cache_dir, **options) as backend:
        return backend.evaluate_pipe(params, p_a, p_b, t_a, t_b, save_path=save_path, profiles=profiles)
