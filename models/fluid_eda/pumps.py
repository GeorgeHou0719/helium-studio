"""Signed forward/backflow pump law, separate from the legacy pump modes.

Q = Sf(Pa) * (Pa - Pb / K0(Pb)), referenced to the inlet temperature.
The HiPace starting values are nominal He data, not a measured 2-D map.
See PHYSICS.md for temperature assumptions, sources and calibration limits.
"""
from __future__ import annotations

import math
from numbers import Real

import numpy as np

from models.comsol_pipe.helium_properties import (
    CP_J_KG_K as CP, SPECIFIC_GAS_CONSTANT_J_KG_K as R,
)


HIPACE_SOURCE = "https://www.pfeiffervacuum.com/ca/en/product/hipace-300-m.html"
FORWARD_BACKFLOW_DEFAULTS = {
    "forward_speed_m3_s": 0.240,
    "compression_ratio": 1e8,
    "forward_speed_curve": [],       # [inlet pressure Pa, forward speed m3/s]
    "compression_curve": [],         # [outlet pressure Pa, zero-flow compression]
    "max_inlet_pressure_Pa": None,
    "max_outlet_pressure_Pa": None,
    "max_throughput_Pa_m3_s": None,
}


class PumpModelError(ValueError):
    def __init__(self, message, field, *, value=None, maximum=None, port=None):
        super().__init__(message)
        self.context = {"field": field}
        for key, item in (("value", value), ("maximum", maximum), ("port", port)):
            if item is not None:
                self.context[key] = item


def _positive(value, field, *, minimum=None, maximum=1e20):
    if (isinstance(value, bool) or not isinstance(value, Real)
            or not math.isfinite(value) or value <= 0
            or (minimum is not None and value < minimum) or value > maximum):
        rule = f"at least {minimum:g}" if minimum is not None else "positive"
        raise PumpModelError(f"{field} must be finite, {rule}, and at most {maximum:g}.", field)
    return float(value)


def validate_forward_backflow(params):
    """Return normalized SI parameters without changing the saved project."""
    values = {key: params.get(key, default) for key, default in FORWARD_BACKFLOW_DEFAULTS.items()}
    values["forward_speed_m3_s"] = _positive(values["forward_speed_m3_s"], "forward_speed_m3_s", minimum=1e-15, maximum=1e6)
    values["compression_ratio"] = _positive(values["compression_ratio"], "compression_ratio", minimum=1)
    for key in ("forward_speed_curve", "compression_curve"):
        rows = values[key]
        if not isinstance(rows, list) or (len(rows) != 0 and not 2 <= len(rows) <= 128):
            raise PumpModelError(f"{key} must be empty or contain 2–128 [pressure, value] rows.", key)
        normalized, previous = [], 0.0
        for index, row in enumerate(rows):
            label = f"{key} row {index + 1}"
            if not isinstance(row, (list, tuple)) or len(row) != 2:
                raise PumpModelError(f"{label} requires exactly two numbers.", key)
            pressure = _positive(row[0], label + " pressure", minimum=1e-9, maximum=1e9)
            value = _positive(row[1], label + " value", minimum=1 if key == "compression_curve" else 1e-15,
                              maximum=1e20 if key == "compression_curve" else 1e6)
            if pressure <= previous:
                raise PumpModelError(f"{label}: pressures must be strictly increasing without duplicates.", key)
            normalized.append([pressure, value])
            previous = pressure
        values[key] = normalized
    for key in ("max_inlet_pressure_Pa", "max_outlet_pressure_Pa", "max_throughput_Pa_m3_s"):
        if values[key] is not None:
            values[key] = _positive(values[key], key, maximum=1e9)
    return values


def _curve(rows, pressure, fallback):
    if not rows:
        return fallback, False
    # Positive power-law interpolation; endpoint holds keep solver trials finite.
    # Any accepted point outside a table is explicitly flagged in the results.
    x, y = np.asarray(rows, dtype=float).T
    value = float(np.exp(np.interp(math.log(pressure), np.log(x), np.log(y))))
    return value, bool(pressure < x[0] or pressure > x[-1])


def evaluate_forward_backflow(params, pa, pb, ta, tb, *, diagnostics=True):
    """Evaluate validated parameters. A/B identities never swap on reversal.

    Both directional throughputs use the same inlet-temperature reference.
    Different bath temperatures do not constitute a calibrated nonisothermal
    turbopump model; an explicit warning accompanies that approximation.
    """
    speed, speed_outside = _curve(params["forward_speed_curve"], pa, params["forward_speed_m3_s"])
    compression, compression_outside = _curve(params["compression_curve"], pb, params["compression_ratio"])
    backward_speed = speed / compression
    q_forward, q_backward = speed * pa, backward_speed * pb
    throughput = speed * (pa - pb / compression)
    mass = throughput / (R * ta)
    result = {
        "pump_mode": "forward_backflow", "mass_flow": mass, "mass_flow_kg_s": mass,
        "outlet_fluid_temperature_K": tb if mass >= 0 else ta,
        "enthalpy_change_W": mass * CP * (tb - ta),
        "pumping_speed_m3_s": throughput / pa, "warnings": [],
    }
    if not diagnostics:
        return result
    direction = "balanced" if abs(throughput) <= 1e-12 * max(q_forward, q_backward) else "forward" if mass > 0 else "backflow"
    result["diagnostics"] = {
        "forward_speed_m3_s": speed, "backward_speed_m3_s": backward_speed,
        "compression_ratio": compression, "pressure_ratio": pb / pa,
        "forward_mass_flow_kg_s": q_forward / (R * ta),
        "backflow_mass_flow_kg_s": q_backward / (R * ta),
        "net_throughput_Pa_m3_s": throughput,
        "backflow_fraction": q_backward / q_forward,
        "zero_flow_inlet_pressure_Pa": pb / compression,
        "flow_direction": direction, "throughput_reference_temperature_K": ta,
        "curve_range_exceeded": speed_outside or compression_outside,
        "model_validated": False,
    }
    warnings = result["warnings"]
    warnings.append("Forward–backflow is a reduced helium model, not a validated two-pressure performance map. "
                    "HiPace starting values are 240 L/s and K0=1e8; confirm the pump revision and calibrate the curves. "
                    "Compression work, pump cooling and speed reduction under load are not calculated.")
    if not params["forward_speed_curve"] or not params["compression_curve"]:
        warnings.append("Untabulated coefficients remain constant. Nominal speed and compression do not establish the loaded operating envelope.")
    if speed_outside:
        warnings.append("Inlet pressure is outside the forward-speed table; the nearest endpoint value was held.")
    if compression_outside:
        warnings.append("Outlet pressure is outside the compression table; the nearest endpoint value was held.")
    if not math.isclose(ta, tb, rel_tol=1e-6):
        warnings.append("Pump port temperatures differ. Both directional throughputs use port A temperature as their reference; "
                        "nonisothermal pump performance has not been calibrated.")
    if direction == "backflow":
        warnings.append("Backflow exceeds forward pumping: net flow is B → A. The signed flow is retained in mass conservation.")
    if any(params[key] is None for key in ("max_inlet_pressure_Pa", "max_outlet_pressure_Pa", "max_throughput_Pa_m3_s")):
        warnings.append("The helium operating envelope is not fully specified. Supply verified limits and model the forevacuum arrangement.")
    return result


def check_operating_limits(params, pa, pb, result):
    """Check only converged operating points, never nonlinear solver trials."""
    for field, value, port, unit, label in (
        ("max_inlet_pressure_Pa", pa, "a", "Pa", "Inlet pressure"),
        ("max_outlet_pressure_Pa", pb, "b", "Pa", "Outlet pressure"),
        ("max_throughput_Pa_m3_s", abs(result["diagnostics"]["net_throughput_Pa_m3_s"]), None, "Pa m³/s", "Absolute net throughput"),
    ):
        limit = params[field]
        if limit is not None and value > limit * (1 + 1e-9):
            raise PumpModelError(f"{label} {value:.6g} {unit} exceeds the entered limit of {limit:.6g} {unit}.",
                                 field, value=value, maximum=limit, port=port)
