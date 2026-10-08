"""Independent pump limits, signed network conservation and input validation."""
import copy
import json
import math

import pytest

from models.fluid_eda.network import NetworkError, R, solve_network, validate_network, refresh_result
from models.fluid_eda.pumps import evaluate_forward_backflow, validate_forward_backflow
from models.fluid_eda.simulation import run_simulation
from tests.test_network import analytical_pipe, circuit, component, pipe, source


def pump_circuit(pin=1, pout=100, **params):
    return circuit([source("inlet", pin), source("outlet", pout),
                    component("pump", "pump", mode="forward_backflow", outlet_temperature_K=293.15, **params)],
                   [("inlet:a", "pump:a"), ("pump:b", "outlet:a")])


def test_default_matches_hi_pace_and_does_not_mutate_input():
    graph = pump_circuit()
    before = copy.deepcopy(graph)
    result = solve_network(graph)
    pump = result["components"]["pump"]
    assert pump["mass_flow"] == pytest.approx(.240 * (1 - 100 / 1e8) / (R * 293.15))
    assert pump["diagnostics"]["backward_speed_m3_s"] == pytest.approx(2.4e-9)
    assert pump["diagnostics"]["forward_speed_m3_s"] == .240
    assert graph == before
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize("ratio, sign", [(0.5, 1), (1, 0), (2, -1)])
def test_compression_equilibrium_and_signed_backflow(ratio, sign):
    # At K0=1000 and fixed 100 Pa outlet, zero flow requires 0.1 Pa inlet.
    result = solve_network(pump_circuit(pin=.1, pout=100 * ratio, compression_ratio=1000))
    row = result["components"]["pump"]
    d = row["diagnostics"]
    expected = .240 * .1 * (1-ratio) / (R * 293.15)
    assert row["mass_flow"] == pytest.approx(expected, abs=1e-24)
    assert d["forward_mass_flow_kg_s"] - d["backflow_mass_flow_kg_s"] == pytest.approx(expected, rel=1e-14, abs=1e-22)
    assert d["flow_direction"] == {1: "forward", 0: "balanced", -1: "backflow"}[sign]
    assert result["components"]["inlet"]["mass_flow"] == pytest.approx(expected, abs=1e-24)
    assert result["components"]["outlet"]["mass_flow"] == pytest.approx(-expected, abs=1e-24)
    assert result["diagnostics"]["net_source_mass_flow_kg_s"] == pytest.approx(0, abs=1e-24)
    assert refresh_result(result)["diagnostics"]["net_source_mass_flow_kg_s"] == pytest.approx(0, abs=1e-24)


def test_unknown_pressure_network_with_net_reverse_flow_matches_quadratic():
    # source -- passive pipe -- pump -- high-pressure source; net B->A flow.
    c, alpha = 1e-9, .01/(R*293.15)
    graph = circuit([source("low", 1), pipe("tube", coefficient=c),
                     component("pump", "pump", mode="forward_backflow", forward_speed_m3_s=.01,
                               compression_ratio=10, outlet_temperature_K=293.15), source("high", 100)],
                    [("low:a", "tube:a"), ("tube:b", "pump:a"), ("pump:b", "high:a")])
    result = solve_network(graph, analytical_pipe)
    # c*(1-p^2) = alpha*(p-10)
    expected = 2*(c+10*alpha)/(alpha+math.sqrt(alpha*alpha+4*c*(c+10*alpha)))
    assert result["components"]["pump"]["pressure_in"] == pytest.approx(expected, rel=1e-7)
    assert result["components"]["pump"]["mass_flow"] < 0
    assert result["components"]["pump"]["mass_flow"] == pytest.approx(result["components"]["tube"]["mass_flow"], rel=1e-7)
    assert result["diagnostics"]["max_relative_mass_residual"] < 1e-7


def test_low_backpressure_reduces_to_legacy_constant_speed():
    forward = pump_circuit(pin=1, pout=1, compression_ratio=1e20, forward_speed_m3_s=.004)
    legacy = copy.deepcopy(forward)
    legacy["components"][2]["params"].update(mode="constant_speed", speed_m3_s=.004)
    assert solve_network(forward)["components"]["pump"]["mass_flow"] == pytest.approx(
        solve_network(legacy)["components"]["pump"]["mass_flow"], rel=1e-12)


def test_k_equal_one_becomes_passive_pressure_conductance():
    result = solve_network(pump_circuit(pin=1, pout=2, compression_ratio=1))
    assert result["components"]["pump"]["mass_flow"] == pytest.approx(-.240/(R*293.15))


def test_log_log_curves_override_nominal_and_endpoints_warn():
    p = validate_forward_backflow({"forward_speed_m3_s": 9,
                                  "forward_speed_curve": [[1,.1],[100,.001]],
                                  "compression_curve": [[10,1e6],[1000,1e2]]})
    row = evaluate_forward_backflow(p, 10, 100, 293.15, 293.15)
    assert row["diagnostics"]["forward_speed_m3_s"] == pytest.approx(.01)
    assert row["diagnostics"]["compression_ratio"] == pytest.approx(1e4)
    assert not row["diagnostics"]["curve_range_exceeded"]
    row = evaluate_forward_backflow(p, .1, 2000, 293.15, 293.15)
    assert row["diagnostics"]["forward_speed_m3_s"] == pytest.approx(.1)
    assert row["diagnostics"]["compression_ratio"] == pytest.approx(1e2)
    assert row["diagnostics"]["curve_range_exceeded"]
    assert sum("nearest endpoint" in w for w in row["warnings"]) == 2


@pytest.mark.parametrize("params", [
    {"forward_speed_m3_s": 0}, {"forward_speed_m3_s": True},
    {"forward_speed_m3_s": float("inf")}, {"compression_ratio": .5},
    {"compression_ratio": "1e8"}, {"compression_ratio": float("nan")},
    {"forward_speed_curve": None}, {"forward_speed_curve": [[1,.1]]},
    {"forward_speed_curve": [[1,.1],[1,.2]]}, {"forward_speed_curve": [[2,.1],[1,.2]]},
    {"forward_speed_curve": [[1,.1],[2,0]]}, {"forward_speed_curve": [[1,.1],[2,.2,3]]},
    {"compression_curve": [[1,10],[2,.2]]}, {"compression_curve": [[0,10],[2,20]]},
    {"compression_curve": [[1,10],[2,True]]}, {"max_outlet_pressure_Pa": -1},
    {"max_inlet_pressure_Pa": float("nan")}, {"max_throughput_Pa_m3_s": 0},
])
def test_invalid_parameters_name_the_pump(params):
    with pytest.raises(NetworkError) as error:
        validate_network(pump_circuit(**params))
    assert error.value.component_ids == ["pump"]
    assert error.value.context["field"]
    assert error.value.stage == "Pump parameter validation"


@pytest.mark.parametrize("key,limit,port", [
    ("max_inlet_pressure_Pa", .5, "a"),
    ("max_outlet_pressure_Pa", 50, "b"),
    ("max_throughput_Pa_m3_s", .1, None),
])
def test_operating_limit_errors_name_parameter_and_port(key, limit, port):
    with pytest.raises(NetworkError) as error:
        solve_network(pump_circuit(**{key:limit}))
    assert error.value.code == "pump_operating_limit"
    assert error.value.component_ids == ["pump"]
    assert error.value.context["field"] == key
    assert error.value.context.get("port") == port
    assert error.value.context["value"] > limit


def test_limit_check_accepts_boundary_and_ignores_solver_trials():
    graph = circuit([source("source", 1), pipe("tube", coefficient=1e-9),
                     component("pump", "pump", mode="forward_backflow", outlet_temperature_K=293.15,
                               max_inlet_pressure_Pa=.1), source("outlet", 10)],
                    [("source:a", "tube:a"), ("tube:b", "pump:a"), ("pump:b", "outlet:a")])
    # Initial guess is outside the limit; its converged inlet is ~0.00254 Pa.
    result = solve_network(graph, analytical_pipe)
    assert result["components"]["pump"]["pressure_in"] < .1
    solve_network(pump_circuit(max_inlet_pressure_Pa=1, max_outlet_pressure_Pa=100))


def test_reverse_legacy_expression_behavior_is_preserved():
    graph = pump_circuit()
    graph["components"][2]["params"].update(mode="expression", expression="-1e-8", speed_m3_s=.004)
    with pytest.raises(NetworkError, match="reverse flow") as error:
        solve_network(graph)
    assert error.value.code == "reverse_pump_flow"


def test_new_mode_ignores_dormant_legacy_parameters_and_reports_own_warnings():
    graph = pump_circuit(speed_m3_s=-1, expression="unused invalid expression")
    result = run_simulation(graph)
    assert result["status"] == "complete"
    assert not any("constant-speed pump" in w for w in result["warnings"])
    assert any("not a validated" in w for w in result["warnings"])


def test_port_temperature_difference_is_explicit_inlet_reference_assumption():
    p = validate_forward_backflow({})
    row = evaluate_forward_backflow(p, 1, 100, 200, 300)
    assert row["mass_flow"] == pytest.approx(.24*(1-100/1e8)/(R*200))
    assert any("temperatures differ" in w for w in row["warnings"])
    assert row["diagnostics"]["throughput_reference_temperature_K"] == 200


def test_continuous_loop_with_hi_pace_defaults_conserves_mass():
    from models.fluid_eda.examples import continuous_pump_example
    graph = continuous_pump_example()
    graph["settings"]["engine"] = "reduced"
    graph["components"][2]["params"]["mode"] = "forward_backflow"
    result = run_simulation(graph)
    flows = [row["mass_flow"] for row in result["edges"]]
    assert min(flows) > 0
    assert max(flows) - min(flows) < 1e-6 * max(flows)
    assert result["summary"]["mass_balance_relative"] < 1e-6
    assert abs(result["components"]["pot"]["mass_flow"]) < 1e-8 * max(flows)
