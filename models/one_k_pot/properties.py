"""Bounded He-4 properties for a *reduced*, unlicensed 1-K-pot model.

Saturation P, latent heat and liquid enthalpy: Donnelly & Barenghi (1998),
DOI 10.1063/1.556028, Tables 17.2, 16.2 and 7.5 respectively. The liquid-H
coefficients were checked visually against printed p. 1242 (PDF p. 26).
Primary author mirror: https://pages.uoregon.edu/rjd/vapor8.htm
Primary reference: https://doi.org/10.1063/1.556028 (PDF not bundled).

This is NOT HEPAK or a general fluid EOS. For the compressed feed, we use
h(T,P) ~= h_sat_liquid(T) + (P-Psat(T))/145 kg/m3. It ignores thermal
expansion, compressibility and the pressure shift of the lambda transition.
Benchmark this liquid correction before design predictions. The active model
prescribes the post-HEX5 liquid state; no return-vapor warming EOS is used.

Public quantities use K, Pa, J/kg and kg/m3. No property extrapolation.
"""

from __future__ import annotations

import numpy as np
from scipy.interpolate import BSpline, PchipInterpolator
from scipy.optimize import brentq

MOLAR_MASS_KG_MOL = 0.0040026  # DB98's molar conversion convention.
MIN_T_K = 0.8
MAX_T_K = 4.3
MAX_FEED_PRESSURE_PA = 200_000.0  # Deliberately narrow approximation scope.
PRESSURE_CORRECTION_DENSITY_KGM3 = 145.0
LIQUID_REFERENCE_DENSITY_KGM3 = 125.4075  # DB98 Table 1.2: SVP at 4.20 K.
PROPERTY_SOURCE = "Donnelly & Barenghi (1998), DOI 10.1063/1.556028"

# P in Pa, not log(P); coefficients previously checked in he4_simplified.
_P_KNOTS = [0.5] * 4 + [
    .52, .53, .56, .6, .65, .7, .75, .8, .85, .92, 1, 1.1, 1.25,
    1.5, 1.7, 1.85, 2, 2.2, 2.5, 2.7, 3, 3.3, 3.7, 4.05, 4.5,
] + [5.1] * 4
_P_COEFF = [
    .00183797, .0025, .00376, .00624, .01236, .03306, .09339, .25856,
    .6269, 1.37025, 2.93305, 6.37318, 15.0042, 38.0570, 116.050, 362.700,
    990.900, 1886., 3194., 5569.5, 9279.1, 15370., 23480., 37355., 57050.,
    87170., 132825., 179650., 211567.,
]
# Latent heat in J/mol; divide by molar mass once at the public boundary.
_L_KNOTS = [0.] * 4 + [
    1.268, 1.980630, 2.1768, 2.1768, 2.1768, 3.7265,
    4.3789, 5.106127, 5.178168,
] + [5.195767] * 4
_L_COEFF = [
    59.82983, 68.59642, 82.78576, 94.70215, 92.50350, 90.73236,
    93.68727, 97.81899, 80.05877, 64.58454, 29.47641, 15.94783, .0003110047,
]
# Liquid H in J/mol, relative to its ground-state enthalpy. This is the
# publication's convenience spline, not an exact interpolant of Table 7.6.
_H_KNOTS = [0.] * 4 + [
    .1, .2, .3, .37, .5, .61, .743950, .87, 1.027550, 1.25, 1.5,
    1.75, 2., 2.13, 2.190120, 2.232940, 2.8542, 4.,
] + [4.9] * 4
_H_COEFF = [
    0., 0., -1.65e-6, 1.801e-5, 1.185e-4, 4.10855e-4, 1.0332e-3,
    2.710e-3, 6.09e-3, 1.784e-2, 6.032740e-2, .24005, .8646372,
    2.4454, 5.55436, 8.8, 12.2013, 15.581, 20.04, 32.231, 47.0464, 58.4893,
]
_pressure = BSpline(_P_KNOTS, _P_COEFF, 3, extrapolate=False)
_latent = BSpline(_L_KNOTS, _L_COEFF, 3, extrapolate=False)
_liquid_h = BSpline(_H_KNOTS, _H_COEFF, 3, extrapolate=False)
_pressure_derivative = _pressure.derivative()
_liquid_h_derivative = _liquid_h.derivative()


def _numeric(value, name):
    try:
        result = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be numeric.") from error
    if np.any(~np.isfinite(result)):
        raise ValueError(f"{name} must be finite.")
    return result


def _temperature(value):
    result = _numeric(value, "Temperature")
    if np.any(result < MIN_T_K) or np.any(result > MAX_T_K):
        raise ValueError(f"Properties support {MIN_T_K:g}-{MAX_T_K:g} K only.")
    return result


def _feed_pressure(value):
    result = _numeric(value, "Feed pressure")
    if np.any(result <= 0) or np.any(result > MAX_FEED_PRESSURE_PA):
        raise ValueError(
            f"Reduced feed properties require 0 < P <= {MAX_FEED_PRESSURE_PA:g} Pa."
        )
    return result


def _scalar_or_array(value):
    result = np.asarray(value)
    return float(result) if result.ndim == 0 else result


def saturation_pressure(T):
    """Saturated He-4 pressure [Pa], with an explicit 0.8-4.3 K guard."""
    return _scalar_or_array(_pressure(_temperature(T)))


def saturation_temperature(P):
    """Exact bounded inversion of the saturation-pressure spline [K]."""
    pressure = _numeric(P, "Saturation pressure")
    p_min, p_max = _pressure(MIN_T_K), _pressure(MAX_T_K)
    if np.any(pressure < p_min) or np.any(pressure > p_max):
        raise ValueError(
            f"Saturation pressure must be {p_min:.6g}-{p_max:.6g} Pa."
        )
    result = np.array([
        brentq(lambda T: float(_pressure(T)) - p, MIN_T_K, MAX_T_K, xtol=1e-12)
        for p in pressure.ravel()
    ]).reshape(pressure.shape)
    return _scalar_or_array(result)


def latent_heat(T):
    """Saturation vaporization enthalpy [J/kg], not J/mol."""
    return _scalar_or_array(_latent(_temperature(T)) / MOLAR_MASS_KG_MOL)


def saturated_liquid_enthalpy(T):
    """DB98 Table 7.5 saturated-liquid H spline [J/kg]."""
    return _scalar_or_array(_liquid_h(_temperature(T)) / MOLAR_MASS_KG_MOL)


def saturated_vapor_enthalpy(T):
    """Saturated-vapor H [J/kg], on the same reference as liquid H."""
    return _scalar_or_array(saturated_liquid_enthalpy(T) + latent_heat(T))


def liquid_enthalpy(T, P):
    """Approximate compressed-liquid H [J/kg]; reject a vapor-side state.

    Constant-density pressure correction is an engineering closure, not a
    source-provided compressed-liquid property. Pressure in Pa is absolute.
    """
    temperature, pressure = np.broadcast_arrays(_temperature(T), _feed_pressure(P))
    p_sat = _pressure(temperature)
    if np.any(pressure < p_sat - 1e-7):
        raise ValueError("Liquid feed pressure is below saturation at its temperature.")
    result = (_liquid_h(temperature) / MOLAR_MASS_KG_MOL
              + (pressure - p_sat) / PRESSURE_CORRECTION_DENSITY_KGM3)
    return _scalar_or_array(result)


def _pressure_free_h(T):
    return (_liquid_h(T) / MOLAR_MASS_KG_MOL
            - _pressure(T) / PRESSURE_CORRECTION_DENSITY_KGM3)


# With fixed correction density, h(T,P) = g(T) + P/rho. One monotone
# inverse suffices at every pressure. PCHIP supplies a vectorized initial
# estimate; two exact-spline Newton updates remove the interpolation error.
_inverse_grid = np.unique(np.r_[np.linspace(MIN_T_K, MAX_T_K, 3501),
                              [t for t in _H_KNOTS if MIN_T_K <= t <= MAX_T_K],
                              2.1768])
_inverse_h = _pressure_free_h(_inverse_grid)
if not np.all(np.diff(_inverse_h) > 0):
    raise RuntimeError("Reduced liquid enthalpy is not monotone; inversion is unsafe.")
_inverse_liquid = PchipInterpolator(_inverse_h, _inverse_grid, extrapolate=False)


def liquid_temperature(h, P):
    """Vectorized inversion of liquid_enthalpy [K], with domain/phase checks."""
    enthalpy, pressure = np.broadcast_arrays(_numeric(h, "Liquid enthalpy"),
                                            _feed_pressure(P))
    target = enthalpy - pressure / PRESSURE_CORRECTION_DENSITY_KGM3
    # Permit only floating-point endpoint roundoff, never physical extrapolation.
    tolerance = 1e-8  # J/kg
    if np.any(target < _inverse_h[0] - tolerance) or np.any(target > _inverse_h[-1] + tolerance):
        raise ValueError("Liquid enthalpy lies outside the 0.8-4.3 K property domain.")
    target = np.clip(target, _inverse_h[0], _inverse_h[-1])
    temperature = np.asarray(_inverse_liquid(target))
    for _ in range(2):
        gradient = (_liquid_h_derivative(temperature) / MOLAR_MASS_KG_MOL
                    - _pressure_derivative(temperature) / PRESSURE_CORRECTION_DENSITY_KGM3)
        temperature = temperature - (_pressure_free_h(temperature) - target) / gradient
    temperature = np.clip(temperature, MIN_T_K, MAX_T_K)
    if np.any(pressure < _pressure(temperature) - 1e-7):
        raise ValueError("Requested liquid enthalpy implies a vapor-side state at this pressure.")
    return _scalar_or_array(temperature)
