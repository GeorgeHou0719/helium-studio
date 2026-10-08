"""Steady, mass-conserving helium networks assembled from schematic ports.

Wires are ideal connections (one pressure and one prescribed node temperature).
Component-to-component conduction is excluded. Advected gas energy is retained
at each node, including the incoming enthalpy of an adiabatic JT stream.
Temperature anchors are physical baths: pipe end temperatures, reservoirs, and
pump discharge temperatures. The pipe evaluator solves the distributed thermal
problem. An adiabatic ideal-gas restriction transports its inlet enthalpy to the
receiving bath; it does not cool gas to the receiving bath temperature.

Pressure unknowns use logarithms, and Kirchhoff mass balances are solved with
SciPy least_squares. A fast pipe evaluator can be replaced by a COMSOL-corrected
constitutive law without changing topology or mass balance.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from numbers import Real
from typing import Callable

import numpy as np
from scipy.optimize import brentq, least_squares

from models.comsol_pipe.helium_properties import (
    CP_J_KG_K, SPECIFIC_GAS_CONSTANT_J_KG_K,
)
from models.continuous_he4.properties import saturation_pressure
from .safe_expression import ScalarExpression, ExpressionError
from .pumps import (PumpModelError, validate_forward_backflow,
                    evaluate_forward_backflow, check_operating_limits)


R = SPECIFIC_GAS_CONSTANT_J_KG_K
CP = CP_J_KG_K
PORTS = {"source": ("a", "b"), "junction": ("a",),
         "pipe": ("a", "b"), "pump": ("a", "b"), "jt": ("a", "b")}
FLOW_TYPES = {"pipe", "pump", "jt"}


def _phase_screen(edge):
    """Flag supersaturated gas states; do not invent a two-phase solution."""
    if edge.get("type") != "pipe":
        return
    prefix = "The single-phase gas assumption is invalid"
    edge["warnings"] = [w for w in edge.get("warnings", []) if not w.startswith(prefix)]
    ratios = []
    for point in edge.get("profile", []):
        temperature = point.get("temperature", point.get("temperature_K"))
        pressure = point.get("pressure", point.get("pressure_Pa"))
        if temperature is None or pressure is None:
            continue
        if 1.0 <= temperature <= 4.3 and pressure > 0:
            ratios.append((float(pressure) / float(saturation_pressure(temperature)), float(temperature)))
    maximum, at_temperature = max(ratios, default=(0.0, None))
    diagnostics = edge.setdefault("diagnostics", {})
    diagnostics["gas_phase_screen_passed"] = bool(maximum <= 1 + 1e-6)
    diagnostics["max_pressure_over_saturation"] = maximum if ratios else None
    if maximum > 1 + 1e-6:
        edge["warnings"].append(
            f"{prefix}: T={at_temperature:.4g} K at P/Psat={maximum:.4g}; "
            "condensation may occur, The current result excludes real-fluid and two-phase behavior."
        )


def _advected_energy(edge, ta, tb):
    """Signed a-to-b energy transport, using COMSOL total flux when present."""
    mass = edge["mass_flow"]
    diagnostics = edge.get("diagnostics", {})
    if edge["type"] == "pipe" and all(key in diagnostics for key in ("flow_energy_a", "flow_energy_b")):
        qa, qb = float(diagnostics["flow_energy_a"]), float(diagnostics["flow_energy_b"])
        if not math.isfinite(qa) or not math.isfinite(qb):
            raise NetworkError(f"{edge['id']}: COMSOL convective energy flow is not finite.", [edge["id"]], "component_evaluation")
        return qa, qb
    return mass * CP * ta, mass * CP * tb


THERMAL_NODE_MODEL = "no_intercomponent_conduction"
THERMAL_CLOSURE = "prescribed port states + advected fluid energy; no intercomponent conduction"


def _fluid_energy_duty(node, energy_out, supply):
    """Keep advected energy; pipe-end conduction is isolated from this node.

    Prescribed port temperatures remain boundary states, not an implicit sink
    for an incoming JT stream. Pot phase/enthalpy accounting is applied to the
    final network result separately using actual incoming and outgoing streams.
    """
    net_out = energy_out - supply * CP * node["temperature"]
    return {
        "thermal_model": THERMAL_NODE_MODEL,
        "fluid_cooling_power_W": net_out,
        "fluid_heat_into_node_W": -net_out,
        "cooling_load_W": -net_out,
        "external_heat_W": net_out,
    }


class NetworkError(ValueError):
    def __init__(self, message, component_ids=(), code="invalid_network", *,
                 stage=None, hint=None, context=None, details=None):
        super().__init__(message)
        self.message = message
        self.component_ids = list(dict.fromkeys(component_ids))
        self.code = code
        self.stage = stage
        self.hint = hint
        self.context = context
        self.details = details

    def as_dict(self):
        return {"message": self.message, "component_ids": self.component_ids,
                "code": self.code,
                **{key: getattr(self, key) for key in ("stage", "hint", "context", "details")
                   if getattr(self, key) is not None}}


def _number(value, name, cid, *, positive=False, minimum=None, maximum=None):
    if isinstance(value, bool) or not isinstance(value, Real):
        raise NetworkError(f"{cid}: {name} must be a finite number.", [cid], "invalid_parameter")
    value = float(value)
    if not math.isfinite(value) or (positive and value <= 0):
        raise NetworkError(f"{cid}: {name} must be finite and positive.", [cid], "invalid_parameter")
    if minimum is not None and value < minimum or maximum is not None and value > maximum:
        raise NetworkError(f"{cid}: {name} Outside the supported range [{minimum}, {maximum}].", [cid], "parameter_range")
    return value


def _temperature(value, name, cid):
    return _number(value, name, cid, positive=True)


@dataclass
class CompiledNetwork:
    components: dict
    nodes: dict
    edges: list
    port_nodes: dict
    settings: dict


def validate_network(circuit: dict) -> CompiledNetwork:
    """Validate topology, laws, dimensions, and thermal/pressure anchors.

    Returns the compiled port-to-node graph; failures include component IDs for
    schematic highlighting. Geometry coordinates are intentionally irrelevant.
    """
    if not isinstance(circuit, dict):
        raise NetworkError("The circuit must be a JSON object.")
    raw_components, wires = circuit.get("components"), circuit.get("wires")
    if not isinstance(raw_components, list) or not raw_components:
        raise NetworkError("Add components first.", code="empty_network")
    if len(raw_components) > 128 or not isinstance(wires, list) or len(wires) > 512:
        raise NetworkError("The circuit supports at most 128 components and 512 wires.", code="network_size")
    settings = circuit.get("settings", {})
    if not isinstance(settings, dict):
        raise NetworkError("settings must be an object.")
    if not isinstance(settings.get("initial_pressures_Pa", {}), dict):
        raise NetworkError("initial_pressures_Pa must be a node-pressure object.", code="invalid_initial_state")
    components, parent, connections = {}, {}, {}

    def find(port):
        while parent[port] != port:
            parent[port] = parent[parent[port]]
            port = parent[port]
        return port

    def union(a, b):
        parent[find(b)] = find(a)

    for component in raw_components:
        if not isinstance(component, dict):
            raise NetworkError("Each component must be an object.")
        cid, kind = component.get("id"), component.get("type")
        if not isinstance(cid, str) or not cid or len(cid) > 80 or ":" in cid:
            raise NetworkError("Component ID must have a length of 1–80 characters, and must not contain a colon.")
        if cid in components:
            raise NetworkError(f"Duplicate component ID: {cid}.", [cid], "duplicate_id")
        if kind not in PORTS:
            raise NetworkError(f"{cid}: Unsupported component type {kind}.", [cid], "unsupported_component")
        params = component.get("params", {})
        if not isinstance(params, dict):
            raise NetworkError(f"{cid}: params must be an object.", [cid])
        components[cid] = dict(component, params=dict(params))
        for port in PORTS[kind]:
            key = f"{cid}:{port}"
            parent[key], connections[key] = key, 0
        if kind == "source":
            union(f"{cid}:a", f"{cid}:b")
    seen_wires, wire_ids = set(), set()
    for wire in wires:
        if not isinstance(wire, dict):
            raise NetworkError("Each wire must be an object.")
        wid = wire.get("id")
        if not isinstance(wid, str) or not wid or wid in wire_ids:
            raise NetworkError("Each wire requires a unique nonempty string ID.", code="duplicate_wire_id")
        wire_ids.add(wid)
        endpoints = []
        for side in ("from", "to"):
            endpoint = wire.get(side)
            if not isinstance(endpoint, dict):
                raise NetworkError(f"Wire {wid} missing {side} Port.", code="invalid_wire")
            cid, port = endpoint.get("component"), endpoint.get("port")
            key = f"{cid}:{port}"
            if key not in parent:
                raise NetworkError(f"Wire {wid} references a nonexistent port {key}.", [cid] if isinstance(cid, str) else [], "invalid_wire")
            endpoints.append(key)
        pair = tuple(sorted(endpoints))
        if pair[0] == pair[1] or pair in seen_wires:
            raise NetworkError(f"Wire {wid} Duplicated or connected to the same port.", [k.split(":")[0] for k in pair], "invalid_wire")
        seen_wires.add(pair)
        for port in endpoints:
            connections[port] += 1
        union(*endpoints)

    for cid, component in components.items():
        kind = component["type"]
        if kind in FLOW_TYPES:
            open_ports = [p for p in PORTS[kind] if not connections[f"{cid}:{p}"]]
            if open_ports:
                raise NetworkError(f"{cid}: Port {', '.join(open_ports)} Not connected.", [cid], "open_port")
        elif sum(connections[f"{cid}:{p}"] for p in PORTS[kind]) == 0:
            raise NetworkError(f"{cid}: The component is not connected.", [cid], "open_port")

    roots, nodes, port_nodes = {}, {}, {}
    for port in parent:
        root = find(port)
        if root not in roots:
            node_id = f"n{len(roots) + 1}"
            roots[root] = node_id
            nodes[node_id] = {"id": node_id, "ports": [], "pressure": None,
                              "temperature": None, "source_ids": [], "temperature_anchors": []}
        nid = roots[root]
        cid, p = port.split(":")
        nodes[nid]["ports"].append({"component": cid, "port": p})
        port_nodes[port] = nid

    def anchor_temperature(nid, value, cid, label):
        temperature = _temperature(value, label, cid)
        previous = nodes[nid]["temperature"]
        if previous is not None and not math.isclose(previous, temperature, rel_tol=1e-8, abs_tol=1e-6):
            owners = [a["component"] for a in nodes[nid]["temperature_anchors"]]
            raise NetworkError(f"{cid}: Conflicting temperature boundaries on the same connected node ({previous:g} K and {temperature:g} K).", owners + [cid], "conflicting_temperature")
        nodes[nid]["temperature"] = temperature
        nodes[nid]["temperature_anchors"].append({"component": cid, "parameter": label})

    edges = []
    for cid, component in components.items():
        kind, params = component["type"], component["params"]
        a = port_nodes[f"{cid}:a"]
        b = port_nodes.get(f"{cid}:b")
        if kind in FLOW_TYPES:
            if a == b:
                raise NetworkError(f"{cid}: Both ends are shorted by ideal connections.", [cid], "shorted_component")
            edges.append({"id": cid, "type": kind, "node_a": a, "node_b": b, "params": params})
        if kind == "source":
            if params.get("thermal_mode", "fixed_temperature") != "fixed_temperature":
                raise NetworkError(f"{cid}: This pressure source requires a specified fixed temperature.", [cid], "unsupported_thermal_mode")
            temperature = _temperature(params.get("temperature_K"), "temperature_K", cid)
            _number(params.get("heat_load_W", 0.0), "heat_load_W", cid, minimum=0)
            anchor_temperature(a, temperature, cid, "temperature_K")
            mode = params.get("mode", "constant")
            if mode == "constant":
                pressure = _number(params.get("pressure_Pa"), "pressure_Pa", cid, positive=True)
            elif mode == "saturation":
                try:
                    pressure = float(saturation_pressure(temperature))
                except ValueError as error:
                    raise NetworkError(f"{cid}: {error}", [cid], "property_range") from error
            elif mode == "expression":
                try:
                    pressure = ScalarExpression(params.get("expression"), ("T",))(T=temperature)
                except ExpressionError as error:
                    raise NetworkError(f"{cid}: {error}", [cid], "invalid_expression") from error
                _number(pressure, "P(T)", cid, positive=True)
            else:
                raise NetworkError(f"{cid}: Unknown pressure-source mode {mode}.", [cid], "unsupported_law")
            previous = nodes[a]["pressure"]
            if previous is not None and not math.isclose(previous, pressure, rel_tol=1e-8, abs_tol=1e-8):
                raise NetworkError(f"{cid}: Fixed-pressure sources on the same node have conflicting pressures.", nodes[a]["source_ids"] + [cid], "conflicting_pressure")
            if not 1e-9 <= pressure <= 1e9:
                raise NetworkError(f"{cid}: Pressure must be within the numerical range 1e-9–1e9 Pa.", [cid], "parameter_range")
            nodes[a]["pressure"] = pressure
            nodes[a]["source_ids"].append(cid)
        elif kind == "junction":
            if params.get("temperature_K") is not None:
                anchor_temperature(a, params["temperature_K"], cid, "temperature_K")
        elif kind == "pipe":
            length = _number(params.get("length_m"), "length_m", cid, positive=True)
            diameter = _number(params.get("outer_diameter_m"), "outer_diameter_m", cid, positive=True)
            wall = _number(params.get("wall_thickness_m"), "wall_thickness_m", cid, positive=True)
            if 2 * wall >= diameter:
                raise NetworkError(f"{cid}: Outer diameter must exceed twice the wall thickness.", [cid], "invalid_geometry")
            if length < 1e-6 or length > 1e4 or diameter < 1e-7 or diameter > 10:
                raise NetworkError(f"{cid}: Geometry is outside the supported numerical range.", [cid], "invalid_geometry")
            for endpoint, nid in (("a", a), ("b", b)):
                value = params.get(f"T_{endpoint}_K")
                if value is not None:
                    anchor_temperature(nid, value, cid, f"T_{endpoint}_K")
        elif kind == "pump":
            mode = params.get("mode", "constant_speed")
            if mode not in ("constant_speed", "expression", "forward_backflow"):
                raise NetworkError(f"{cid}: Unknown pump-curve mode.", [cid], "unsupported_law")
            if mode == "forward_backflow":
                try:
                    params.update(validate_forward_backflow(params))
                except PumpModelError as error:
                    raise NetworkError(f"{cid}: {error}", [cid], "invalid_parameter",
                                       stage="Pump parameter validation", context=error.context) from error
            else:
                params["speed_m3_s"] = _number(params.get("speed_m3_s", 0.004), "speed_m3_s", cid, positive=True)
            if mode == "expression":
                try:
                    component["compiled_expression"] = ScalarExpression(params.get("expression"), ("Pin", "Pout", "Tin", "S", "R"))
                except ExpressionError as error:
                    raise NetworkError(f"{cid}: {error}", [cid], "invalid_expression") from error
            if params.get("outlet_temperature_K") is not None:
                anchor_temperature(b, params["outlet_temperature_K"], cid, "outlet_temperature_K")
        elif kind == "jt":
            law = params.get("law", "orifice")
            if law not in ("orifice", "capillary"):
                raise NetworkError(f"{cid}: Unknown J–T flow-resistance model.", [cid], "unsupported_law")
            _number(params.get("diameter_m"), "diameter_m", cid, minimum=1e-9, maximum=1.0)
            if law == "capillary":
                _number(params.get("length_m"), "length_m", cid, positive=True)
                if params.get("impedance_mode", "geometry") not in ("geometry", "measured"):
                    raise NetworkError(f"{cid}: Unknown capillary impedance mode.", [cid], "unsupported_law")
                if params.get("impedance_mode") == "measured":
                    _number(params.get("impedance_m_inv3"), "impedance_m_inv3", cid, positive=True)
            else:
                _number(params.get("discharge_coefficient", 0.7), "discharge_coefficient", cid, positive=True, maximum=1.0)

    if not edges:
        raise NetworkError("The circuit requires at least one pipe, pump or restriction.", code="no_flow_components")
    fixed_nodes = [nid for nid, node in nodes.items() if node["pressure"] is not None]
    if not fixed_nodes:
        raise NetworkError("Add at least one pressure source as an absolute-pressure reference.", code="missing_pressure_reference")
    adjacency = {nid: set() for nid in nodes}
    for edge in edges:
        a, b = edge["node_a"], edge["node_b"]
        adjacency[a].add(b)
        adjacency[b].add(a)
    visited, queue = set(), [fixed_nodes[0]]
    while queue:
        nid = queue.pop()
        if nid in visited:
            continue
        visited.add(nid)
        queue.extend(adjacency[nid] - visited)
    if visited != set(nodes):
        ids = [p["component"] for nid in nodes.keys() - visited for p in nodes[nid]["ports"]]
        raise NetworkError("A region is disconnected from the main circuit; Connect it before running a simulation.", ids, "disconnected_network")

    # A massless, adiabatic ideal-gas JT element preserves T. This only fills
    # unanchored nodes; different prescribed baths remain distinct heat sinks.
    changed = True
    while changed:
        changed = False
        for edge in edges:
            if edge["type"] != "jt":
                continue
            a, b = nodes[edge["node_a"]], nodes[edge["node_b"]]
            if (a["temperature"] is None) != (b["temperature"] is None):
                blank, known = (a, b) if a["temperature"] is None else (b, a)
                blank["temperature"] = known["temperature"]
                blank["temperature_inferred"] = "ideal-gas adiabatic restriction"
                changed = True
    for nid, node in nodes.items():
        if node["temperature"] is None:
            ids = [p["component"] for p in node["ports"]]
            raise NetworkError(f"Node {nid} Missing temperature boundary; Set the connected pipe-end temperature, pump outlet or node temperature.", ids, "missing_temperature_anchor")
    for edge in edges:
        if edge["type"] == "pipe":
            params, cid = edge["params"], edge["id"]
            boundary = params.get("wall_boundary", "adiabatic")
            if boundary not in ("adiabatic", "prescribed"):
                raise NetworkError(f"{cid}: Unknown outer-wall thermal boundary.", [cid], "invalid_wall_boundary")
            if boundary == "prescribed":
                segments = params.get("wall_profile")
                if not isinstance(segments, list) or not 1 <= len(segments) <= 32:
                    raise NetworkError(f"{cid}: Outer-wall temperature requires 1–32 temperature segments.", [cid], "invalid_wall_profile")
                length, previous = 0.0, nodes[edge["node_a"]]["temperature"]
                for segment in segments:
                    if not isinstance(segment, dict):
                        raise NetworkError(f"{cid}: Each wall-temperature segment must be an object.", [cid], "invalid_wall_profile")
                    length += _number(segment.get("length_m"), "segment.length_m", cid, minimum=1e-6, maximum=100)
                    start = _temperature(segment.get("T_start_K"), "segment.T_start_K", cid)
                    end = _temperature(segment.get("T_end_K"), "segment.T_end_K", cid)
                    if abs(previous - start) > 1e-6:
                        raise NetworkError(f"{cid}: Wall-temperature segments must be continuous, and the first segment must match the inlet temperature.", [cid], "invalid_wall_profile")
                    previous = end
                if abs(previous - nodes[edge["node_b"]]["temperature"]) > 1e-6:
                    raise NetworkError(f"{cid}: The final outer-wall temperature must match the outlet temperature.", [cid], "invalid_wall_profile")
                if not math.isclose(length, params["length_m"], rel_tol=1e-6, abs_tol=1e-8):
                    raise NetworkError(f"{cid}: Outer-wall segment lengths must sum to the pipe length.", [cid], "invalid_wall_profile")
    return CompiledNetwork(components, nodes, edges, port_nodes, settings)


def solve_network(circuit: dict, pipe_evaluator: Callable | None = None,
                  progress: Callable | None = None) -> dict:
    """Solve node pressure and obtain thermal profiles for every component.

    ``pipe_evaluator(params, p_a, p_b, T_a, T_b, profiles=False)`` must return
    mass_flow_kg_s (positive a to b); with profiles=True, also a profile and
    outward conductive heat_a_W / heat_b_W into the endpoint baths.
    """
    graph = validate_network(circuit)
    if pipe_evaluator is None:
        from .thermal import reduced_pipe
        pipe_evaluator = reduced_pipe
    if any(edge["type"] == "jt" for edge in graph.edges):
        from .jt import evaluate_restriction

    nodes, edges = graph.nodes, graph.edges
    unknown = [nid for nid, node in nodes.items() if node["pressure"] is None]
    fixed = {nid: node["pressure"] for nid, node in nodes.items() if node["pressure"] is not None}
    baseline = math.exp(float(np.mean(np.log(list(fixed.values())))))
    initial = {nid: baseline for nid in unknown}
    for edge in edges:
        if edge["type"] == "pump":
            if edge["node_a"] in initial:
                initial[edge["node_a"]] = baseline * 0.5
            if edge["node_b"] in initial:
                initial[edge["node_b"]] = baseline * 10.0
    # A pump connected directly to a pressure reservoir through a suction
    # pipe has a useful scalar initial estimate. It is only an initial guess;
    # the full graph still determines every final flow and pressure.
    seed_flows = []
    for pump in [e for e in edges if e["type"] == "pump"]:
        inlet = pump["node_a"]
        if inlet not in initial:
            continue
        for pipe in edges:
            if pipe["type"] != "pipe" or pipe["node_b"] != inlet or pipe["node_a"] not in fixed:
                continue
            pa = fixed[pipe["node_a"]]
            ta, tb = nodes[pipe["node_a"]]["temperature"], nodes[inlet]["temperature"]
            speed = pump["params"].get("speed_m3_s", 0.004)
            pump_expression = graph.components[pump["id"]].get("compiled_expression")
            pout_guess = initial.get(pump["node_b"], fixed.get(pump["node_b"], baseline))
            def pump_seed_flow(pin):
                if pump["params"].get("mode") == "forward_backflow":
                    return evaluate_forward_backflow(pump["params"], pin, pout_guess, tb,
                                                     nodes[pump["node_b"]]["temperature"],
                                                     diagnostics=False)["mass_flow_kg_s"]
                if pump_expression is not None:
                    return pump_expression(Pin=pin, Pout=pout_guess, Tin=tb, S=speed, R=R)
                return pin * speed / (R * tb)
            def suction_balance(pb):
                row = pipe_evaluator(pipe["params"], pa, pb, ta, tb, profiles=False)
                return float(row.get("mass_flow_kg_s", row.get("mass_flow"))) - pump_seed_flow(pb)
            try:
                guess = brentq(suction_balance, max(1e-8, pa * 1e-8), pa, xtol=1e-11)
                initial[inlet] = guess
                seed_flows.append(max(pump_seed_flow(guess), 0.0))
            except (ValueError, ArithmeticError):
                pass
    if seed_flows:
        target_flow = float(sum(seed_flows))
        known_seeds = dict(fixed)
        known_seed_flows = {}
        # The orifice pressure scale can be orders of magnitude above pot
        # pressure. Seed it by matching a representative pump throughput.
        groups = {}
        for jt in [e for e in edges if e["type"] == "jt"]:
            a, b = jt["node_a"], jt["node_b"]
            if b in fixed and a in initial:
                groups.setdefault((a, b), []).append(jt)
        # Parallel restrictions share a pressure drop; seed their SUMMED flow.
        # Separate return branches receive a conductance-weighted throughput
        # estimate. The full network solve subsequently corrects every split.
        weights = {}
        for (a, b), group in groups.items():
            probe = min(1e9, max(fixed[b] * 100, 100.0))
            weights[a, b] = sum(abs(evaluate_restriction(
                jt["params"], probe, fixed[b], nodes[a]["temperature"], nodes[b]["temperature"]
            )["mass_flow_kg_s"]) for jt in group)
        weight_sum = sum(weights.values())
        for (a, b), group in groups.items():
            ta, tb = nodes[a]["temperature"], nodes[b]["temperature"]
            group_target = target_flow * weights[a, b] / weight_sum if weight_sum else target_flow / len(groups)
            def jt_balance(pa):
                return sum(float(evaluate_restriction(jt["params"], pa, fixed[b], ta, tb)["mass_flow_kg_s"])
                           for jt in group) - group_target
            try:
                upper = max(fixed[b] * 2, 1.0)
                while upper < 1e9 and jt_balance(upper) < 0:
                    upper *= 2
                guessed = brentq(jt_balance, fixed[b], min(upper, 1e9), xtol=1e-9)
                known_seeds[a] = initial[a] = guessed
                known_seed_flows[a] = known_seed_flows.get(a, 0.0) + group_target
            except (ValueError, ArithmeticError):
                pass
        # Work upstream through passive pipes from known discharge/restriction
        # pressure guesses. Avoid a huge artificial trial gradient that would
        # exceed the thermal model's permitted trial mass flow.
        for _ in range(len(edges)):
            changed = False
            for pipe in [e for e in edges if e["type"] == "pipe"]:
                a, b = pipe["node_a"], pipe["node_b"]
                if b not in known_seeds or a in known_seeds or a not in initial:
                    continue
                pb = known_seeds[b]
                ta, tb = nodes[a]["temperature"], nodes[b]["temperature"]
                inlet_count = sum(e["type"] == "pipe" and e["node_b"] == b for e in edges)
                branch_target = known_seed_flows.get(b, target_flow) / max(inlet_count, 1)
                def pipe_balance(pa):
                    row = pipe_evaluator(pipe["params"], pa, pb, ta, tb, profiles=False)
                    return float(row.get("mass_flow_kg_s", row.get("mass_flow"))) - branch_target
                try:
                    delta = max(pb * 1e-7, 1e-7)
                    upper = pb + delta
                    while upper < 1e9 and pipe_balance(upper) < 0:
                        delta *= 2
                        upper = pb + delta
                    guess = brentq(pipe_balance, pb, min(upper, 1e9), xtol=1e-10)
                    known_seeds[a] = initial[a] = guess
                    known_seed_flows[a] = known_seed_flows.get(a, 0.0) + branch_target
                    changed = True
                except (ValueError, ArithmeticError):
                    pass
            if not changed:
                break
    physical_initial = dict(initial)
    for nid in unknown:
        raw = graph.settings.get("initial_pressures_Pa", {}).get(nid)
        if raw is not None:
            initial[nid] = _number(raw, "initial_pressure", nid, positive=True)

    def pressures(x):
        return {**fixed, **{nid: float(math.exp(value)) for nid, value in zip(unknown, x)}}

    calls = 0

    def edge_results(p, profiles=False):
        nonlocal calls
        calls += 1
        results = []
        for edge in edges:
            cid, kind, params = edge["id"], edge["type"], edge["params"]
            a, b = edge["node_a"], edge["node_b"]
            ta, tb = nodes[a]["temperature"], nodes[b]["temperature"]
            try:
                if kind == "pipe":
                    result = pipe_evaluator(params, p[a], p[b], ta, tb, profiles=profiles)
                elif kind == "jt":
                    result = evaluate_restriction(params, p[a], p[b], ta, tb)
                elif params.get("mode") == "forward_backflow":
                    result = evaluate_forward_backflow(params, p[a], p[b], ta, tb, diagnostics=profiles)
                else:
                    speed = params["speed_m3_s"]
                    if params.get("mode", "constant_speed") == "expression":
                        mdot = graph.components[cid]["compiled_expression"](Pin=p[a], Pout=p[b], Tin=ta, S=speed, R=R)
                    else:
                        mdot = p[a] * speed / (R * ta)
                    result = {"pump_mode": params.get("mode", "constant_speed"),
                              "mass_flow": mdot, "mass_flow_kg_s": mdot,
                              "temperature_in": ta, "temperature_out": tb,
                              "outlet_fluid_temperature_K": tb,
                              "enthalpy_change_W": mdot * CP * (tb - ta),
                              "pumping_speed_m3_s": speed, "warnings": []}
                mdot = float(result.get("mass_flow_kg_s", result.get("mass_flow")))
                if not math.isfinite(mdot):
                    raise ValueError("Flow rate is not finite.")
            except (ValueError, ArithmeticError, TypeError) as error:
                if isinstance(error, NetworkError):
                    raise
                raise NetworkError(f"{cid}: {error}", [cid], "component_evaluation") from error
            result = dict(result, id=cid, label=graph.components[cid].get("label", cid),
                          type=kind, node_a=a, node_b=b, mass_flow=mdot,
                          mass_flow_kg_s=mdot, pressure_in=p[a], pressure_out=p[b],
                          temperature_in=ta, temperature_out=tb)
            if kind == "jt":
                result["outlet_fluid_temperature_K"] = ta if mdot >= 0 else tb
                result["temperature_out"] = result["outlet_fluid_temperature_K"]
                result["receiving_node_temperature_K"] = tb if mdot >= 0 else ta
                if mdot < 0 and result.get("profile"):
                    old_profile = result["profile"]
                    length = float(old_profile[-1]["x"])
                    result["profile"] = [dict(row, x=length - float(row["x"]))
                                         for row in reversed(old_profile)]
            if profiles and not result.get("profile") and result.get("pump_mode") != "forward_backflow":
                tout = result.get("outlet_fluid_temperature_K", tb)
                result["profile"] = [
                    {"x": 0.0, "pressure": p[a], "temperature": ta,
                     "wall_temperature": None, "velocity": None},
                    {"x": 1.0, "pressure": p[b], "temperature": tout,
                     "wall_temperature": None, "velocity": None},
                ]
            results.append(result)
        return results

    def balances(results):
        balance, total = {nid: 0.0 for nid in nodes}, {nid: 0.0 for nid in nodes}
        for edge in results:
            a, b, mdot = edge["node_a"], edge["node_b"], edge["mass_flow"]
            balance[a] += mdot
            balance[b] -= mdot
            total[a] += abs(mdot)
            total[b] += abs(mdot)
        return balance, total

    flow_floor = 1e-16
    initial_p = {**fixed, **initial}
    initial_results = edge_results(initial_p)
    _, initial_total = balances(initial_results)
    characteristic_flow = max([abs(e["mass_flow"]) for e in initial_results] or [flow_floor])
    residual_scales = {nid: max(initial_total[nid], characteristic_flow * 0.1, flow_floor)
                       for nid in nodes}

    def residual(x):
        balance, total = balances(edge_results(pressures(x)))
        # Fixed characteristic scaling preserves the KCL Jacobian. Dividing
        # by instantaneous absolute flows would flatten residuals whenever all
        # trial flows point out of a node, and hinders high-pressure tees.
        return np.array([balance[nid] / residual_scales[nid] for nid in unknown])

    def relative_error(x):
        balance, total = balances(edge_results(pressures(x)))
        return max([abs(balance[nid]) / max(total[nid], flow_floor) for nid in unknown] or [0.0])

    if progress:
        progress("Solving node pressures and mass conservation")
    if unknown:
        x0 = np.clip(np.log([initial[nid] for nid in unknown]), math.log(1e-9) + 1e-8, math.log(1e9) - 1e-8)
        if initial != physical_initial:
            physical_x0 = np.clip(np.log([physical_initial[nid] for nid in unknown]), math.log(1e-9) + 1e-8, math.log(1e9) - 1e-8)
            # COMSOL corrections alter branch conductances. A stale warm start
            # can be much worse than the newly rebuilt hydraulic estimate.
            if relative_error(physical_x0) < relative_error(x0):
                x0 = physical_x0
        solution = least_squares(residual, x0, bounds=(math.log(1e-9), math.log(1e9)),
                                 xtol=2e-10, ftol=2e-10, gtol=2e-10,
                                 max_nfev=350, diff_step=2e-7)
        norm = relative_error(solution.x)
        # A second scale in the initial pressure improves convergence for a
        # very high-compression closed pump loop, without accepting a bad root.
        if norm > 2e-6:
            second = np.array(x0)
            for edge in edges:
                if edge["type"] == "pump" and edge["node_b"] in unknown:
                    second[unknown.index(edge["node_b"])] = min(math.log(baseline * 1000), math.log(1e9) - 1e-6)
            try:
                retry = least_squares(residual, second, bounds=(math.log(1e-9), math.log(1e9)),
                                      xtol=2e-10, ftol=2e-10, gtol=2e-10,
                                      max_nfev=350, diff_step=2e-7)
                retry_norm = relative_error(retry.x)
                if retry_norm < norm:
                    solution, norm = retry, retry_norm
            except NetworkError:
                pass
        if norm > 2e-5:
            raise NetworkError(f"Mass conservation did not converge (maximum relative residual {norm:.3g}); Check the pump curve, boundary conditions and restriction dimensions.", [e["id"] for e in edges], "nonconvergence")
        final_p, nfev = pressures(solution.x), int(solution.nfev)
    else:
        final_p, nfev = fixed, 0
    if progress:
        progress("Calculating pipe temperatures, pressure profiles and energy balance")
    final_edges = edge_results(final_p, profiles=True)
    balance, total = balances(final_edges)
    max_relative = max([abs(balance[nid]) / max(total[nid], flow_floor) for nid in unknown] or [0.0])
    if max_relative > 2e-5:
        raise NetworkError("Detailed profiles are inconsistent with the solved flow resistance, Mass-conservation check failed.", [e["id"] for e in edges], "inconsistent_evaluator")

    # Node energy includes flowing gas only. Pipe-end conduction stays in
    # pipe diagnostics and never contributes to another component or pot.
    energy_out = {nid: 0.0 for nid in nodes}
    warnings = []
    for edge in final_edges:
        _phase_screen(edge)
        a, b, mdot = edge["node_a"], edge["node_b"], edge["mass_flow"]
        ta, tb = nodes[a]["temperature"], nodes[b]["temperature"]
        if edge.get("pump_mode") == "forward_backflow":
            try:
                check_operating_limits(graph.components[edge["id"]]["params"], final_p[a], final_p[b], edge)
            except PumpModelError as error:
                raise NetworkError(f"{edge['id']}: {error}", [edge["id"]], "pump_operating_limit",
                                   stage="Pump operating point", context=error.context,
                                   hint="Check this pump's operating point, forevacuum arrangement and entered helium limits.") from error
        elif edge["type"] == "pump" and mdot < -1e-16:
            raise NetworkError(f"{edge['id']}: The pump curve predicts reverse flow at the solution; Adjust the curve or reverse the pump orientation.", [edge["id"]], "reverse_pump_flow")
        if edge["type"] == "jt":
            upstream_T = ta if mdot >= 0 else tb
            energy_out[a] += mdot * CP * upstream_T
            energy_out[b] -= mdot * CP * upstream_T
        else:
            fa, fb = _advected_energy(edge, ta, tb)
            energy_out[a] += fa
            energy_out[b] -= fb
        for warning in edge.get("warnings", []):
            warnings.append(f"{edge['label']}: {warning}")

    node_results = []
    component_results = {e["id"]: e for e in final_edges}
    for nid, node in nodes.items():
        # Fixed-pressure reservoirs may provide/remove net mass at their own
        # temperature. Subtract that supply enthalpy from external thermal duty.
        supply = balance[nid] if node["source_ids"] else 0.0
        duty = _fluid_energy_duty(node, energy_out[nid], supply)
        result = dict(node, pressure=final_p[nid], pressure_Pa=final_p[nid],
                      temperature_K=node["temperature"], fixed_pressure=nid in fixed,
                      mass_imbalance_kg_s=balance[nid], source_mass_flow_kg_s=supply,
                      **duty)
        node_results.append(result)
        for cid in node["source_ids"]:
            component_results[cid] = {
                "id": cid, "label": graph.components[cid].get("label", cid),
                "type": "source", "node": nid, "pressure": final_p[nid],
                "temperature": node["temperature"], "pressure_in": final_p[nid],
                "pressure_out": final_p[nid], "temperature_in": node["temperature"],
                "temperature_out": node["temperature"],
                "mass_flow": supply if len(node["source_ids"]) == 1 else None,
                **duty,
                "profile": [], "warnings": [],
            }
        if len(node["source_ids"]) > 1:
            warnings.append(f"Node {nid} multiple equal-pressure sources share the total flow; Individual source supply rates cannot be determined uniquely.")
        for port in node["ports"]:
            cid = port["component"]
            if graph.components[cid]["type"] == "junction":
                component_results[cid] = {"id": cid, "type": "junction", "node": nid,
                                          "pressure": final_p[nid], "temperature": node["temperature"],
                                          **duty,
                                          "profile": [], "warnings": []}
    warnings.append("Intercomponent conduction is excluded from node energy balances. Gas flow energy is retained; prescribed port temperatures do not erase the incoming JT enthalpy.")
    diagnostics = {
        "converged": True, "unknown_pressures": len(unknown), "node_count": len(nodes),
        "edge_count": len(edges), "function_evaluations": nfev, "component_evaluations": calls,
        "max_relative_mass_residual": max_relative,
        "max_mass_residual_kg_s": max([abs(balance[nid]) for nid in unknown] or [0.0]),
        "net_source_mass_flow_kg_s": float(sum(balance[nid] for nid in fixed)),
        "pressure_min_Pa": min(final_p.values()), "pressure_max_Pa": max(final_p.values()),
        "thermal_node_model": THERMAL_NODE_MODEL,
        "thermal_closure": THERMAL_CLOSURE,
    }
    return {"nodes": node_results, "components": component_results, "edges": final_edges,
            "port_nodes": graph.port_nodes, "diagnostics": diagnostics,
            "warnings": list(dict.fromkeys(warnings)), "engine": graph.settings.get("engine", "reduced")}


def refresh_result(result: dict, max_relative_mass_residual: float = 0.005,
                   raise_on_failure: bool = False) -> dict:
    """Recompute balances after replacing ``result['edges']`` with COMSOL rows.

    The input result is updated and returned. Component rows, reservoir thermal
    fluid energy and diagnostics are refreshed from actual edge flow energies.
    End-conduction fields remain pipe diagnostics, excluded from node balances.
    """
    nodes = {row["id"]: row for row in result["nodes"]}
    balance = {nid: 0.0 for nid in nodes}
    magnitude = {nid: 0.0 for nid in nodes}
    energy = {nid: 0.0 for nid in nodes}
    environment = 0.0
    for edge in result["edges"]:
        _phase_screen(edge)
        cid, a, b = edge["id"], edge["node_a"], edge["node_b"]
        mass = float(edge.get("mass_flow_kg_s", edge.get("mass_flow")))
        if not math.isfinite(mass):
            raise NetworkError(f"{cid}: Nonfinite flow rate.", [cid], "component_evaluation")
        edge["mass_flow"] = edge["mass_flow_kg_s"] = mass
        balance[a] += mass
        balance[b] -= mass
        magnitude[a] += abs(mass)
        magnitude[b] += abs(mass)
        ta, tb = nodes[a]["temperature"], nodes[b]["temperature"]
        if edge["type"] == "jt":
            upstream = ta if mass >= 0 else tb
            energy[a] += mass * CP * upstream
            energy[b] -= mass * CP * upstream
            edge["outlet_fluid_temperature_K"] = upstream
        else:
            fa, fb = _advected_energy(edge, ta, tb)
            energy[a] += fa
            energy[b] -= fb
        environment += float(edge.get("heat_environment_W", 0))
        result["components"][cid] = edge
    unknown, fixed = [], []
    for nid, node in nodes.items():
        sources = node.get("source_ids", [])
        (fixed if sources else unknown).append(nid)
        supply = balance[nid] if sources else 0.0
        duty = _fluid_energy_duty(node, energy[nid], supply)
        node.update(mass_imbalance_kg_s=balance[nid], source_mass_flow_kg_s=supply,
                    **duty)
        for cid in sources:
            result["components"][cid].update(mass_flow=supply if len(sources) == 1 else None,
                                              **duty)
        for port in node.get("ports", []):
            cid = port["component"]
            if result["components"].get(cid, {}).get("type") == "junction":
                result["components"][cid].update(**duty)
    relative = max([abs(balance[nid]) / max(magnitude[nid], 1e-16) for nid in unknown] or [0.0])
    result["diagnostics"].update(
        converged=relative <= max_relative_mass_residual,
        max_relative_mass_residual=relative,
        max_mass_residual_kg_s=max([abs(balance[nid]) for nid in unknown] or [0.0]),
        net_source_mass_flow_kg_s=float(sum(balance[nid] for nid in fixed)),
        pipe_heat_to_environment_W=environment,
        thermal_node_model=THERMAL_NODE_MODEL,
        thermal_closure=THERMAL_CLOSURE,
    )
    prefixes = tuple(f"{edge.get('label', edge['id'])}: " for edge in result["edges"])
    warnings = [warning for warning in result.get("warnings", [])
                if not warning.startswith(prefixes)]
    for edge in result["edges"]:
        warnings.extend(f"{edge.get('label', edge['id'])}: {warning}"
                        for warning in edge.get("warnings", []))
    result["warnings"] = list(dict.fromkeys(warnings))
    if raise_on_failure and relative > max_relative_mass_residual:
        raise NetworkError(f"COMSOL Branch flows do not yet satisfy node mass conservation ({relative:.3g}).", [e["id"] for e in result["edges"]], "comsol_mass_balance")
    return result
