"""Analytical graph checks, energy signs, invalid circuits, and expression safety."""
import copy
import math

import pytest

from models.fluid_eda.network import NetworkError, R, CP, solve_network, validate_network, refresh_result
from models.fluid_eda.safe_expression import ExpressionError, ScalarExpression


def component(cid, kind, **params):
    return {"id": cid, "type": kind, "label": cid, "params": params}


def source(cid, pressure, temperature=293.15):
    return component(cid, "source", mode="constant", pressure_Pa=pressure,
                     temperature_K=temperature)


def pipe(cid, ta=293.15, tb=293.15, coefficient=1e-12, **extra):
    return component(cid, "pipe", length_m=0.288, outer_diameter_m=0.00635,
                     wall_thickness_m=0.000762, T_a_K=ta, T_b_K=tb,
                     coefficient=coefficient, **extra)


def wire(n, a, b):
    ca, pa = a.split(":")
    cb, pb = b.split(":")
    return {"id": f"w{n}", "from": {"component": ca, "port": pa},
            "to": {"component": cb, "port": pb}}


def circuit(components, connections):
    return {"components": components,
            "wires": [wire(n, *pair) for n, pair in enumerate(connections)],
            "settings": {"engine": "reduced"}}


def analytical_pipe(params, pa, pb, ta, tb, profiles=False):
    mdot = params["coefficient"] * (pa * pa - pb * pb)
    result = {"mass_flow_kg_s": mdot, "heat_a_W": params.get("heat_a_W", 0),
              "heat_b_W": params.get("heat_b_W", 0), "warnings": []}
    if profiles:
        result["profile"] = [
            {"x": 0, "temperature": ta, "pressure": pa},
            {"x": params["length_m"], "temperature": tb, "pressure": pb},
        ]
    return result


def simple_circuit():
    return circuit([source("high", 100), pipe("pipe"), source("low", 10)],
                   [("high:a", "pipe:a"), ("pipe:b", "low:a")])


def test_two_fixed_reservoirs_have_analytic_flow_and_source_conservation():
    result = solve_network(simple_circuit(), analytical_pipe)
    expected = 1e-12 * (100**2 - 10**2)
    assert result["components"]["pipe"]["mass_flow"] == pytest.approx(expected)
    assert result["components"]["high"]["mass_flow"] == pytest.approx(expected)
    assert result["components"]["low"]["mass_flow"] == pytest.approx(-expected)
    assert result["diagnostics"]["net_source_mass_flow_kg_s"] == pytest.approx(0, abs=1e-24)


def test_series_resistance_matches_pressure_squared_analytic_solution():
    graph = circuit([source("high", 100), pipe("p1", coefficient=1e-12),
                     pipe("p2", coefficient=3e-12), source("low", 10)],
                    [("high:a", "p1:a"), ("p1:b", "p2:a"), ("p2:b", "low:a")])
    result = solve_network(graph, analytical_pipe)
    expected_pmid = math.sqrt((1e-12 * 100**2 + 3e-12 * 10**2) / 4e-12)
    assert result["components"]["p1"]["pressure_out"] == pytest.approx(expected_pmid, rel=1e-7)
    assert result["components"]["p1"]["mass_flow"] == pytest.approx(result["components"]["p2"]["mass_flow"], rel=1e-7)


def test_parallel_branches_and_tee_conserve_mass():
    graph = circuit([source("high", 100), pipe("feed", coefficient=2e-12),
                     component("tee", "junction"), pipe("left", coefficient=1e-12),
                     pipe("right", coefficient=3e-12), source("low", 10)],
                    [("high:a", "feed:a"), ("feed:b", "tee:a"),
                     ("tee:a", "left:a"), ("tee:a", "right:a"),
                     ("left:b", "low:a"), ("right:b", "low:b")])
    result = solve_network(graph, analytical_pipe)
    components = result["components"]
    assert components["right"]["mass_flow"] / components["left"]["mass_flow"] == pytest.approx(3)
    assert components["feed"]["mass_flow"] == pytest.approx(components["left"]["mass_flow"] + components["right"]["mass_flow"], rel=1e-7)
    expected = math.sqrt((2 * 100**2 + 4 * 10**2) / 6)
    assert components["tee"]["pressure"] == pytest.approx(expected, rel=1e-7)


def test_closed_loop_constant_speed_pump_matches_analytic_suction():
    graph = circuit([source("pot", 15, 1), pipe("suction", 1, 293.15),
                     component("pump", "pump", mode="constant_speed", speed_m3_s=.004,
                               outlet_temperature_K=293.15), pipe("return", 293.15, 1)],
                    [("pot:a", "suction:a"), ("suction:b", "pump:a"),
                     ("pump:b", "return:a"), ("return:b", "pot:b")])
    result = solve_network(graph, analytical_pipe)
    speed_factor = .004 / (R * 293.15)
    c = 1e-12
    # Stable positive quadratic root of c*p^2 + k*p - c*Ppot^2 = 0.
    expected = 2 * c * 15**2 / (speed_factor + math.sqrt(speed_factor**2 + 4 * c*c * 15**2))
    assert result["components"]["pump"]["pressure_in"] == pytest.approx(expected, rel=1e-6)
    flows = [e["mass_flow"] for e in result["edges"]]
    assert max(flows) - min(flows) < 1e-6 * max(flows)
    assert abs(result["components"]["pot"]["mass_flow"]) < 1e-16


def test_user_pump_expression_reproduces_constant_speed():
    graph = circuit([source("a", 10), source("b", 100),
                     component("pump", "pump", mode="expression", speed_m3_s=.004,
                               expression="Pin*S/(R*Tin)", outlet_temperature_K=293.15)],
                    [("a:a", "pump:a"), ("pump:b", "b:a")])
    result = solve_network(graph, analytical_pipe)
    assert result["components"]["pump"]["mass_flow"] == pytest.approx(10 * .004 / (R * 293.15))


def test_pipe_end_conduction_is_excluded_from_source_fluid_energy():
    graph = simple_circuit()
    graph["components"][1]["params"].update(heat_a_W=0.1368, heat_b_W=-0.1368)
    result = solve_network(graph, analytical_pipe)
    assert result["components"]["high"]["cooling_load_W"] == pytest.approx(0, abs=1e-15)
    assert result["components"]["low"]["external_heat_W"] == pytest.approx(0, abs=1e-15)


def test_orientation_reverses_signed_flow_without_changing_solution():
    graph = simple_circuit()
    graph["wires"] = [wire(0, "high:a", "pipe:b"), wire(1, "pipe:a", "low:a")]
    result = solve_network(graph, analytical_pipe)
    assert result["components"]["pipe"]["mass_flow"] == pytest.approx(-1e-12*(10000-100))


def test_expression_source_is_evaluated_at_its_temperature():
    graph = simple_circuit()
    graph["components"][0]["params"].update(mode="expression", expression="T/2")
    result = solve_network(graph, analytical_pipe)
    assert result["components"]["high"]["pressure"] == pytest.approx(293.15 / 2)


@pytest.mark.parametrize("mutation,code", [
    (lambda c: c["wires"].pop(), "open_port"),
    (lambda c: c["components"][1]["params"].update(wall_thickness_m=.004), "invalid_geometry"),
    (lambda c: c["components"][1]["params"].update(T_a_K=4), "conflicting_temperature"),
    (lambda c: c["components"][1].update(type="unknown"), "unsupported_component"),
    (lambda c: c["components"][0]["params"].update(pressure_Pa=float("nan")), "invalid_parameter"),
    (lambda c: c["components"][0]["params"].update(pressure_Pa=True), "invalid_parameter"),
])
def test_invalid_networks_are_rejected_with_actionable_codes(mutation, code):
    graph = simple_circuit()
    mutation(graph)
    with pytest.raises(NetworkError) as error:
        validate_network(graph)
    assert error.value.code == code
    assert error.value.component_ids


def test_pressure_sources_on_same_net_cannot_disagree():
    graph = simple_circuit()
    graph["components"].append(source("extra", 99))
    graph["wires"].append(wire(99, "extra:a", "high:a"))
    with pytest.raises(NetworkError, match="conflicting pressures"):
        validate_network(graph)


def test_disconnected_network_rejected():
    graph = simple_circuit()
    other = circuit([source("otherhigh", 100), pipe("otherpipe"), source("otherlow", 10)],
                    [("otherhigh:a", "otherpipe:a"), ("otherpipe:b", "otherlow:a")])
    graph["components"].extend(other["components"])
    for n, w in enumerate(other["wires"]):
        w["id"] = f"other{n}"
    graph["wires"].extend(other["wires"])
    with pytest.raises(NetworkError) as error:
        validate_network(graph)
    assert error.value.code == "disconnected_network"


@pytest.mark.parametrize("expression", [
    "__import__('os').system('whoami')", "(1).__class__", "Pin[0]",
    "[x for x in (1,2)]", "lambda: 1", "open('file')", "True", "1e999",
])
def test_untrusted_expression_cannot_access_python_objects(expression):
    with pytest.raises(ExpressionError):
        ScalarExpression(expression, ("Pin",))


@pytest.mark.parametrize("expression", ["1/0", "sqrt(-1)", "2**1000", "(-1)**0.5", "exp(1000)"])
def test_undefined_expression_fails_cleanly(expression):
    with pytest.raises(ExpressionError):
        ScalarExpression(expression)()


def test_safe_expression_mathematical_functions_and_si_variables():
    expr = ScalarExpression("max(0, Pin - Pout)*S/(R*Tin)", ("Pin", "Pout", "S", "R", "Tin"))
    assert expr(Pin=100, Pout=20, S=.004, R=R, Tin=293.15) == pytest.approx(80*.004/(R*293.15))


def test_validation_does_not_mutate_input():
    graph = simple_circuit()
    before = copy.deepcopy(graph)
    validate_network(graph)
    assert graph == before


def test_real_continuous_loop_conserves_mass_and_jt_does_not_invent_cooling():
    from models.fluid_eda.examples import continuous_pump_example
    result = solve_network(continuous_pump_example())
    flows = [edge["mass_flow"] for edge in result["edges"]]
    assert min(flows) > 0
    assert max(flows) - min(flows) < 1e-6 * max(flows)
    jt = result["components"]["jt1"]
    assert jt["temperature_in"] == pytest.approx(4)
    assert jt["outlet_fluid_temperature_K"] == pytest.approx(4)
    assert jt["receiving_node_temperature_K"] == pytest.approx(1)
    pipe = result["components"]["pipe2"]
    assert pipe["profile"][-1]["temperature"] == pytest.approx(4)
    assert pipe["profile"][-1]["x"] == pytest.approx(1.5)
    expected_load = jt["mass_flow"] * CP * (4 - 1)
    assert result["components"]["pot"]["cooling_load_W"] == pytest.approx(expected_load, rel=1e-8)


def test_refresh_result_uses_actual_replaced_pipe_fluxes_and_exposes_nonconvergence():
    graph = circuit([source("high", 100), pipe("p1"), pipe("p2"), source("low", 10)],
                    [("high:a", "p1:a"), ("p1:b", "p2:a"), ("p2:b", "low:a")])
    result = solve_network(graph, analytical_pipe)
    result["edges"][0].update(mass_flow_kg_s=result["edges"][0]["mass_flow"] * 1.1,
                                heat_a_W=.2, heat_environment_W=.4,
                                warnings=["new real-pipe diagnostic"])
    refreshed = refresh_result(result)
    assert not refreshed["diagnostics"]["converged"]
    assert refreshed["diagnostics"]["max_relative_mass_residual"] == pytest.approx(.1 / 2.1)
    assert refreshed["components"]["high"]["cooling_load_W"] == pytest.approx(0, abs=1e-15)
    assert refreshed["diagnostics"]["pipe_heat_to_environment_W"] == pytest.approx(.4)
    assert any("new real-pipe diagnostic" in w for w in refreshed["warnings"])
    with pytest.raises(NetworkError) as error:
        refresh_result(result, raise_on_failure=True)
    assert error.value.code == "comsol_mass_balance"


@pytest.mark.parametrize("mutation", [
    lambda p: p["wall_profile"][0].update(length_m=.2),
    lambda p: p["wall_profile"][1].update(T_start_K=40),
    lambda p: p["wall_profile"][-1].update(T_end_K=3),
    lambda p: p.update(wall_profile=[]),
])
def test_prescribed_wall_profile_must_be_continuous_and_match_geometry(mutation):
    from models.fluid_eda.examples import continuous_pump_example
    graph = continuous_pump_example()
    params = next(c["params"] for c in graph["components"] if c["id"] == "pipe2")
    mutation(params)
    with pytest.raises(NetworkError) as error:
        validate_network(graph)
    assert error.value.code == "invalid_wall_profile"


def test_source_reservoir_removes_the_sensible_heat_of_an_incoming_jt_stream():
    graph = circuit([source("warm", 100, 20), source("cold", 10, 1),
                     component("jt", "jt", law="orifice", diameter_m=1e-5, discharge_coefficient=.7)],
                    [("warm:a", "jt:a"), ("jt:b", "cold:a")])
    result = solve_network(graph)
    mass = result["components"]["jt"]["mass_flow"]
    assert result["components"]["cold"]["cooling_load_W"] == pytest.approx(mass * CP * 19)
    assert result["components"]["warm"]["cooling_load_W"] == pytest.approx(0, abs=1e-15)


def test_cryogenic_pipe_condensation_is_flagged_without_fabricated_phase_solution():
    graph = circuit([source("high", 100, 1), pipe("pipe", 1, 1), source("low", 10, 1)],
                    [("high:a", "pipe:a"), ("pipe:b", "low:a")])
    result = solve_network(graph, analytical_pipe)
    row = result["components"]["pipe"]
    assert row["diagnostics"]["gas_phase_screen_passed"] is False
    assert row["diagnostics"]["max_pressure_over_saturation"] > 6
    assert any("The single-phase gas assumption is invalid" in warning for warning in result["warnings"])


def test_missing_pipe_temperature_boundary_requires_explicit_user_closure():
    graph = circuit([source("high", 100), pipe("p1", tb=None), pipe("p2", ta=None), source("low", 10)],
                    [("high:a", "p1:a"), ("p1:b", "p2:a"), ("p2:b", "low:a")])
    with pytest.raises(NetworkError) as error:
        validate_network(graph)
    assert error.value.code == "missing_temperature_anchor"


def test_junction_temperature_provides_explicit_missing_pipe_boundary():
    graph = circuit([source("high", 100), pipe("p1", tb=None),
                     component("j", "junction", temperature_K=100),
                     pipe("p2", ta=None), source("low", 10)],
                    [("high:a", "p1:a"), ("p1:b", "j:a"),
                     ("j:a", "p2:a"), ("p2:b", "low:a")])
    result = solve_network(graph, analytical_pipe)
    assert result["components"]["p1"]["temperature_out"] == pytest.approx(100)
    assert result["components"]["p2"]["temperature_in"] == pytest.approx(100)


def test_comsol_corrected_network_refreshes_hydraulic_seed_instead_of_stale_warm_start():
    from models.fluid_eda.examples import continuous_pump_example
    from models.fluid_eda.thermal import reduced_pipe
    graph = continuous_pump_example()
    first = solve_network(graph)
    graph["settings"]["initial_pressures_Pa"] = {n["id"]: n["pressure"] for n in first["nodes"]}

    def corrected(params, pa, pb, ta, tb, profiles=False):
        result = reduced_pipe(params, pa, pb, ta, tb, profiles=profiles)
        result["mass_flow"] *= 1.08
        result["mass_flow_kg_s"] *= 1.08
        return result

    second = solve_network(graph, pipe_evaluator=corrected)
    assert second["diagnostics"]["max_relative_mass_residual"] < 1e-6
    assert second["components"]["pump1"]["mass_flow"] > first["components"]["pump1"]["mass_flow"]


def test_refresh_uses_comsol_integrated_enthalpy_plus_kinetic_energy():
    result = solve_network(simple_circuit(), analytical_pipe)
    edge = result["edges"][0]
    base = edge["mass_flow"] * CP * 293.15
    edge["diagnostics"].update(flow_energy_a=base+.01, flow_energy_b=base+.02)
    refreshed = refresh_result(result)
    assert refreshed["components"]["high"]["external_heat_W"] == pytest.approx(.01)
    assert refreshed["components"]["low"]["cooling_load_W"] == pytest.approx(.02)


@pytest.mark.parametrize("split_pipe", [False, True])
@pytest.mark.parametrize("diameter", [1e-5, 1.6e-5, 2e-6])
def test_high_pressure_return_with_parallel_jt_branches_conserves_every_node(split_pipe, diameter):
    from models.fluid_eda.examples import continuous_pump_example
    graph = continuous_pump_example()
    # This regression exercises orifice area scaling, independent of UI defaults.
    graph["components"][1]["params"].update(length_m=.288,outer_diameter_m=.00635,wall_thickness_m=.000762)
    graph["components"][2]["params"].update(mode="constant_speed",speed_m3_s=.004)
    graph["components"][3]["params"].update(length_m=.9,outer_diameter_m=.00635,wall_thickness_m=.000762,T_b_K=2.8,
        wall_profile=[{"length_m":length,"T_start_K":start,"T_end_K":end} for length,start,end in
                      [(0.3,293.15,35),(0.3,35,35),(0.1,35,4),(0.1,4,4),(0.1,4,2.8)]])
    graph["components"][-1]["params"].update(law="orifice",diameter_m=1e-5)
    jt = copy.deepcopy(graph["components"][-1])
    jt["id"] = "jt2"
    jt["params"]["diameter_m"] = diameter
    graph["components"].append(jt)
    if split_pipe:
        branch = copy.deepcopy(graph["components"][3])
        branch["id"] = "pipe3"
        graph["components"].append(branch)
        graph["wires"].extend([wire(100, "pump1:b", "pipe3:a"), wire(101, "pipe3:b", "jt2:a")])
    else:
        graph["wires"].append(wire(100, "pipe2:b", "jt2:a"))
    graph["wires"].append(wire(102, "jt2:b", "pot:a"))
    result = solve_network(graph)
    rows = result["components"]
    assert result["diagnostics"]["max_relative_mass_residual"] < 1e-6
    assert rows["pump1"]["mass_flow"] == pytest.approx(rows["jt1"]["mass_flow"] + rows["jt2"]["mass_flow"], rel=1e-6)
    if not split_pipe:
        assert rows["jt2"]["mass_flow"] / rows["jt1"]["mass_flow"] == pytest.approx((diameter/1e-5)**2, rel=1e-6)
    else:
        assert rows["pipe2"]["mass_flow"] == pytest.approx(rows["jt1"]["mass_flow"], rel=1e-6)
        assert rows["pipe3"]["mass_flow"] == pytest.approx(rows["jt2"]["mass_flow"], rel=1e-6)


@pytest.mark.parametrize("expression", ["Pin*S/(R*Tin)",
                                         "Pin*S/(R*Tin)*(1-Pout/100000)",
                                         "Pin*S/(R*Tin)/(1+Pout/1000)"])
def test_closed_loop_pump_expression_with_discharge_pressure_dependence(expression):
    from models.fluid_eda.examples import continuous_pump_example
    graph = continuous_pump_example()
    graph["components"][2]["params"].update(mode="expression", expression=expression)
    result = solve_network(graph)
    pump = result["components"]["pump1"]
    expected = ScalarExpression(expression, ("Pin", "Pout", "Tin", "S", "R"))(
        Pin=pump["pressure_in"], Pout=pump["pressure_out"], Tin=pump["temperature_in"], S=.004, R=R)
    assert pump["mass_flow"] == pytest.approx(expected)
    assert result["diagnostics"]["max_relative_mass_residual"] < 1e-6


@pytest.mark.parametrize('ta,tb', [(.9, 310), (310, .9)])
@pytest.mark.parametrize('prescribed', [False, True])
def test_temperature_window_removed_from_network_and_preview(ta, tb, prescribed):
    from models.fluid_eda.thermal import reduced_pipe
    extra = {'wall_boundary': 'prescribed', 'wall_profile': [
        {'length_m': .288, 'T_start_K': ta, 'T_end_K': tb}]} if prescribed else {}
    tube = pipe('pipe', ta, tb, **extra)
    graph = circuit([source('high', 100, ta), tube, source('low', 100, tb)],
                    [('high:a', 'pipe:a'), ('pipe:b', 'low:a')])
    validate_network(graph)
    result = reduced_pipe(tube['params'], 100, 100, ta, tb, profiles=True)
    assert result['profile'][0]['temperature'] == pytest.approx(ta)
    assert result['profile'][-1]['temperature'] == pytest.approx(tb)


@pytest.mark.parametrize("parallel", [False, True])
@pytest.mark.parametrize("reverse", [False, True])
def test_intercomponent_conduction_does_not_enter_fluid_energy_balance(parallel, reverse):
    """Solid-end heat changes never enter a fluid/pot balance in either path."""
    high, low = (10, 100) if reverse else (100, 10)
    components = [source("a", high, 293.15), pipe("p1", 293.15, 20, imposed_heat=1.),
                  component("j", "junction", temperature_K=20), pipe("p2", 20, 4),
                  source("b", low, 4)]
    pairs = [("a:a", "p1:a"), ("p1:b", "j:a"), ("j:a", "p2:a"), ("p2:b", "b:a")]
    if parallel:
        components.append(pipe("p3", 20, 4))
        pairs.extend([("j:a", "p3:a"), ("p3:b", "b:a")])
    graph = circuit(components, pairs)

    def energy_conserving_pipe(params, pa, pb, ta, tb, profiles=False):
        row = analytical_pipe(params, pa, pb, ta, tb, profiles)
        q = params.get("imposed_heat", 0.)
        row.update(heat_a_W=-q, heat_b_W=q-row["mass_flow_kg_s"]*CP*(tb-ta))
        return row

    before = solve_network(graph, energy_conserving_pipe)
    graph["components"][1]["params"]["imposed_heat"] += 10
    after = solve_network(graph, energy_conserving_pipe)
    assert after["components"]["j"]["fluid_cooling_power_W"] == before["components"]["j"]["fluid_cooling_power_W"]
    assert after["components"]["a"]["fluid_cooling_power_W"] == before["components"]["a"]["fluid_cooling_power_W"]
    for cid in ["p2", "b"] + (["p3"] if parallel else []):
        assert after["components"][cid] == before["components"][cid]
    assert after["components"]["p2"]["temperature_in"] == 20
    for node in after["nodes"]:
        assert node["thermal_model"] == "no_intercomponent_conduction"
        assert "reservoir_cooling_W" not in node
        assert node["external_heat_W"] == node["fluid_cooling_power_W"]

    # COMSOL refresh must also exclude end conduction while keeping flow energy.
    updated = copy.deepcopy(before)
    updated["edges"][0]["heat_a_W"] -= 10
    updated["edges"][0]["heat_b_W"] += 10
    refresh_result(updated)
    for old, new in zip(after["nodes"], updated["nodes"]):
        for key in ("fluid_cooling_power_W", "fluid_heat_into_node_W"):
            assert new[key] == pytest.approx(old[key])
    assert updated["components"]["p2"] == before["components"]["p2"]


def test_jt_incoming_energy_is_not_erased_by_receiving_node_temperature():
    graph = circuit([source("warm", 100, 20),
                     component("jt", "jt", law="orifice", diameter_m=1e-5, discharge_coefficient=.7),
                     component("j", "junction", temperature_K=1), pipe("p", 1, 1), source("cold", 10, 1)],
                    [("warm:a", "jt:a"), ("jt:b", "j:a"), ("j:a", "p:a"), ("p:b", "cold:a")])
    result = solve_network(graph, analytical_pipe)
    jt, tube, junction = [result["components"][cid] for cid in ("jt", "p", "j")]
    assert jt["outlet_fluid_temperature_K"] == 20
    assert tube["temperature_in"] == 1
    assert "reservoir_conduction_W" not in junction
    assert junction["fluid_heat_into_node_W"] == pytest.approx(jt["mass_flow"]*CP*19, rel=1e-6)
    assert junction["fluid_heat_into_node_W"] > 0
    assert result["components"]["cold"]["fluid_heat_into_node_W"] == pytest.approx(0, abs=1e-15)
    refresh_result(result)
    assert result["components"]["j"]["fluid_heat_into_node_W"] == pytest.approx(junction["fluid_heat_into_node_W"])
