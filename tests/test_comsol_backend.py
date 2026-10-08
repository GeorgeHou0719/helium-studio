"""Boundary, initialization, cache, and process lifecycle checks without a license.

COMSOL execution is mocked here; a licensed integration check is separate.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import threading
import unittest
from unittest.mock import patch

import numpy as np

from models.fluid_eda.comsol_backend import (
    ComsolCancelled, ComsolError, ComsolPipeBackend, cache_key,
    normalize_parameters,
)
from models.fluid_eda.comsol_worker import initial_profiles, wall_mesh_positions, checked_energy_balance


class PipeInputTests(unittest.TestCase):
    def test_temperatures_outside_old_window_reach_initialization_unchanged(self):
        for ta, tb in [(.9, 310), (310, .9), (.9932677279113648, .9932677279113648)]:
            with self.subTest(ta=ta, tb=tb):
                data = normalize_parameters({}, 15, 10, ta, tb)
                _, temperature, _, _ = initial_profiles(data)
                np.testing.assert_allclose(temperature[[0, -1]], [ta, tb])
        params = {'length_m': .5, 'wall_boundary': 'prescribed',
                  'wall_profile': [{'length_m': .5, 'T_start_K': .9, 'T_end_K': 310}]}
        data = normalize_parameters(params, 15, 10, .9, 310)
        _, temperature, _, _ = initial_profiles(data)
        np.testing.assert_allclose(temperature[[0, -1]], [.9, 310])

    def test_invalid_absolute_temperatures_are_still_rejected(self):
        for value in [0, -1, float('nan'), float('inf'), True]:
            for ta, tb in [(value, 293), (1, value)]:
                with self.assertRaises(ValueError):
                    normalize_parameters({}, 15, 10, ta, tb)
            params = {'wall_boundary': 'prescribed', 'wall_profile': [
                {'length_m': .288, 'T_start_K': 1, 'T_end_K': value}]}
            with self.assertRaises(ValueError):
                normalize_parameters(params, 15, 10, 1, 293)

    def test_solved_temperature_extrema_no_longer_use_property_window(self):
        from types import SimpleNamespace
        from unittest.mock import Mock
        from models.fluid_eda.comsol_worker import PipeSolver
        class ReachedSampling(Exception):
            pass
        solver = PipeSolver.__new__(PipeSolver)
        data = normalize_parameters({}, 15, 10, 1, 293)
        for minimum, maximum in [(.9932677279113648, 293.15), (1, 310)]:
            values = {'comp1.mingas(T)': minimum, 'comp1.maxgas(T)': maximum,
                      'comp1.mingas(p)': 9.74262}
            solver.model = SimpleNamespace(evaluate=lambda expression: values.get(expression, 0.0))
            solver._sample = Mock(side_effect=ReachedSampling)
            with self.assertRaises(ReachedSampling):
                solver.extract(data)
        for invalid in [0, -1, float('nan'), float('inf')]:
            values['comp1.mingas(T)'] = invalid
            with self.assertRaises(RuntimeError):
                solver.extract(data)

    def test_positive_pressure_has_no_adapter_floor_or_ceiling(self):
        for pa, pb in [(15.57, .00248221), (.00248221, 15.57),
                       (1e-12, 300000), (1e7, 1e-12)]:
            with self.subTest(pa=pa, pb=pb):
                data = normalize_parameters({}, pa, pb, 1, 293.15)
                self.assertEqual(data['pressure_a'], pa)
                self.assertEqual(data['pressure_b'], pb)

    def test_nonpositive_or_nonfinite_absolute_pressures_are_rejected(self):
        for value in (0, -1, float('inf'), float('-inf'), float('nan'), True):
            for port in ('a', 'b'):
                with self.subTest(value=value, port=port):
                    with self.assertRaisesRegex(ValueError, 'pressure_' + port):
                        normalize_parameters({}, value if port=='a' else 15,
                                             value if port=='b' else .5, 1, 293.15)

    def test_energy_acceptance_uses_outward_heat_and_advective_signs(self):
        values = {"heat_a": -.2, "heat_b": .01, "heat_environment_W": .29,
                  "flow_energy_a": .11, "flow_energy_b": .01}
        residual, scale = checked_energy_balance(values)
        self.assertAlmostEqual(residual, 0)
        self.assertEqual(scale, .29)
        with self.assertRaisesRegex(RuntimeError, "energy conservation failed"):
            checked_energy_balance({**values, "heat_environment_W": .30})
        near_zero = {key: 0.0 for key in values}
        checked_energy_balance({**near_zero, "heat_a": 9e-9})
        with self.assertRaises(RuntimeError):
            checked_energy_balance({**near_zero, "heat_a": 2e-8})

    def test_prescribed_wall_profile_continuity_length_and_initialization(self):
        profile = [{"length_m": .3, "T_start_K": 293.15, "T_end_K": 35},
                   {"length_m": .3, "T_start_K": 35, "T_end_K": 35},
                   {"length_m": .1, "T_start_K": 35, "T_end_K": 4},
                   {"length_m": .1, "T_start_K": 4, "T_end_K": 4},
                   {"length_m": .1, "T_start_K": 4, "T_end_K": 2.8}]
        params = {"length_m": .9, "wall_boundary": "prescribed", "wall_profile": profile}
        inputs = normalize_parameters(params, 50000.005, 50000, 293.15, 2.8)
        x, temperature, _, _ = initial_profiles(inputs)
        self.assertGreater(float(np.min(np.diff(x))), 1e-12)
        np.testing.assert_allclose(np.interp(np.array([0, .3, .45, .6, .7, .8, .9])/.9, x, temperature),
                                   [293.15, 35, 35, 35, 4, 4, 2.8], atol=1e-10)
        for changed in ({**params, "length_m": 1}, {**params, "wall_profile": []},
                        {**params, "wall_profile": [{"length_m": .9, "T_start_K": 295, "T_end_K": 2.8}]}):
            with self.assertRaises(ValueError):
                normalize_parameters(changed, 50000.005, 50000, 293.15, 2.8)
        self.assertNotEqual(cache_key(inputs, {}), cache_key(normalize_parameters({"length_m": .9}, 50000.005, 50000, 293.15, 2.8), {}))
        mesh_positions = wall_mesh_positions(inputs)
        for transition in np.array([0, .3, .6, .7, .8, .9]) / .9:
            self.assertLess(float(np.min(abs(mesh_positions - transition))), 1e-12)
        np.testing.assert_allclose(mesh_positions, 1 - mesh_positions[::-1], atol=1e-12)
        self.assertTrue(np.all(np.diff(mesh_positions) > 1e-12))

    def test_rejects_invalid_geometry_states_and_mesh(self):
        for params, pa, pb, ta, tb in [
            ({"wall_w": .004}, 15, .5, 1, 293),
            ({"L": -1}, 15, .5, 1, 293),
            ({}, 0, .5, 1, 293),
            ({}, 15, .5, 0, 293),
            ({}, 15, .5, 1, -1),
            ({}, float("nan"), .5, 1, 293),
            ({"axial_elements": 60.5}, 15, .5, 1, 293),
            ({"radial_elements": True}, 15, .5, 1, 293),
        ]:
            with self.subTest(params=params, pa=pa, ta=ta, tb=tb):
                with self.assertRaises(ValueError):
                    normalize_parameters(params, pa, pb, ta, tb)

    def test_initialization_respects_heating_cooling_isothermal_and_reverse_flow(self):
        for pa, pb, ta, tb in [(15, .5, 1, 293.15), (.5, 15, 293.15, 1),
                               (1000, 900, 293.15, 293.15), (1000, 1000, 4, 4)]:
            with self.subTest(pa=pa, pb=pb, ta=ta, tb=tb):
                data = normalize_parameters({}, pa, pb, ta, tb)
                x, temperature, pressure, mass = initial_profiles(data)
                np.testing.assert_allclose(x[[0, -1]], [0, 1])
                np.testing.assert_allclose(temperature[[0, -1]], [ta, tb])
                np.testing.assert_allclose(pressure[[0, -1]], [pa, pb])
                self.assertTrue(np.all(np.diff(x) > 0))
                self.assertTrue(np.all(pressure > 0))
                self.assertEqual(np.sign(mass), np.sign(pa - pb))

    def test_cache_identity_covers_state_geometry_and_fingerprint(self):
        a = normalize_parameters({}, 15, .5, 1, 293.15)
        b = normalize_parameters({"L": .3}, 15, .5, 1, 293.15)
        self.assertNotEqual(cache_key(a, {"template": "a"}), cache_key(b, {"template": "a"}))
        self.assertNotEqual(cache_key(a, {"template": "a"}), cache_key(a, {"template": "b"}))


FAKE_WORKER = r'''
import json,os,sys,time
from pathlib import Path
prefix='FLUID_EDA_JSON '
print('ordinary COMSOL initialization log',flush=True)
print(prefix+json.dumps({'type':'ready','id':None}),flush=True)
for line in sys.stdin:
    req=json.loads(line)
    p=req['inputs']
    print(prefix+json.dumps({'type':'progress','id':req['id'],'message':'working'}),flush=True)
    if p['pressure_a']==666:
        print(prefix+json.dumps({'type':'error','id':req['id'],'error':'license unavailable'}),flush=True)
        continue
    if p['pressure_a']==777:
        time.sleep(30)
    result={'backend':'COMSOL protocol test','mass_flow':p['pressure_a']-p['pressure_b'],'pid':os.getpid(),'profile':[],'diagnostics':{}}
    if req.get('save_path'):
        output=Path(req['save_path'])
        output.parent.mkdir(parents=True,exist_ok=True)
        output.write_bytes(b'computed fake MPH artifact')
        result['artifact']=str(output)
    print(prefix+json.dumps({'type':'result','id':req['id'],'result':result}),flush=True)
'''


class WorkerLifecycleTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.template = self.root / "template.mph"
        self.template.write_bytes(b"test placeholder")
        self.actual_popen = subprocess.Popen

    def tearDown(self):
        self.tmp.cleanup()

    def fake_popen(self, command, **kwargs):
        return self.actual_popen([sys.executable, "-u", "-c", FAKE_WORKER], **kwargs)

    def backend(self, **kwargs):
        return ComsolPipeBackend(self.root / "cache", template=self.template,
                                 executable=Path(sys.executable), **kwargs)

    def test_reuses_one_process_then_cache_without_startup(self):
        with patch("models.fluid_eda.comsol_backend.subprocess.Popen", self.fake_popen):
            with self.backend() as backend:
                first = backend.evaluate_pipe({}, 15, .5, 1, 293.15)
                second = backend.evaluate_pipe({}, 16, .5, 1, 293.15)
                self.assertEqual(first["pid"], second["pid"])
                self.assertFalse(first["cache_hit"])
            with self.backend() as backend:
                cached = backend.evaluate_pipe({}, 15, .5, 1, 293.15)
                self.assertTrue(cached["cache_hit"])
                self.assertIsNone(backend._process)
                self.assertEqual(first["mass_flow"], cached["mass_flow"])

    def test_worker_errors_are_not_cached_or_silently_replaced(self):
        with patch("models.fluid_eda.comsol_backend.subprocess.Popen", self.fake_popen):
            with self.backend() as backend:
                with self.assertRaisesRegex(ComsolError, "license unavailable"):
                    backend.evaluate_pipe({}, 666, 500, 293, 293)
                self.assertFalse((self.root / "cache").exists())

    def test_exact_cache_copies_artifact_without_starting_comsol(self):
        original = self.root / "job_one" / "pipe.mph"
        copied = self.root / "job_two" / "pipe.mph"
        with patch("models.fluid_eda.comsol_backend.subprocess.Popen", self.fake_popen):
            with self.backend() as backend:
                first = backend.evaluate_pipe({}, 15, .5, 1, 293.15, save_path=original)
            # Job output paths are overwritten by later network corrections;
            # their old exact-state cache artifacts must remain reusable.
            original.write_bytes(b"different later-round solution")
            with self.backend() as backend:
                cached = backend.evaluate_pipe({}, 15, .5, 1, 293.15, save_path=copied)
                self.assertTrue(cached["cache_hit"])
                self.assertIsNone(backend._process)
                self.assertEqual(cached["artifact"], str(copied.resolve()))
                self.assertEqual(cached["artifact_sha256"], first["artifact_sha256"])
                self.assertEqual(copied.read_bytes(), b"computed fake MPH artifact")
                same = backend.evaluate_pipe({}, 15, .5, 1, 293.15, save_path=original)
                self.assertTrue(same["cache_hit"])

    def test_missing_or_changed_artifact_requires_an_actual_worker_evaluation(self):
        for change in ("missing", "changed"):
            original = self.root / f"{change}.mph"
            copied = self.root / f"{change}_new.mph"
            pa = 15 if change == "missing" else 16
            with patch("models.fluid_eda.comsol_backend.subprocess.Popen", self.fake_popen):
                with self.backend() as backend:
                    first = backend.evaluate_pipe({}, pa, .5, 1, 293.15, save_path=original)
                cached_artifact = Path(first["artifact_cache"])
                if change == "missing":
                    cached_artifact.unlink()
                else:
                    cached_artifact.write_bytes(b"changed unrelated model")
                with self.backend() as backend:
                    recomputed = backend.evaluate_pipe({}, pa, .5, 1, 293.15, save_path=copied)
                    self.assertFalse(recomputed["cache_hit"])
                    self.assertIsNotNone(backend._process)
                    self.assertEqual(copied.read_bytes(), b"computed fake MPH artifact")

    def test_timeout_and_cancel_terminate_worker(self):
        with patch("models.fluid_eda.comsol_backend.subprocess.Popen", self.fake_popen):
            with self.backend(timeout=.5) as backend:
                with self.assertRaisesRegex(ComsolError, "timeout"):
                    backend.evaluate_pipe({}, 777, 500, 293, 293)
                self.assertIsNone(backend._process)
            cancel = threading.Event()
            with self.backend(cancel_event=cancel) as backend:
                backend.evaluate_pipe({}, 15, .5, 1, 293)
                process = backend._process
                cancel.set()
                with self.assertRaises(ComsolCancelled):
                    backend.evaluate_pipe({}, 15, .5, 1, 293)
                self.assertIsNotNone(process.poll())

    def test_refuses_overwriting_template(self):
        with self.backend() as backend:
            with self.assertRaisesRegex(ValueError, "overwritten"):
                backend.evaluate_pipe({}, 15, .5, 1, 293, save_path=self.template)


if __name__ == "__main__":
    unittest.main()
