"""NIST dilute He4 GAS viscosity, extended through the warm pumping line.

Hurly & Moldover (2000), J. Res. NIST 105, 667-688, Table A1, pp. 682-683.
https://nvlpubs.nist.gov/nistpubs/jres/105/5/j55hur.pdf
The 1-5 K interpolation is retained exactly from he4_simplified; the same
published viscosity and derivative columns extend it to 350 K. The 0.5-1 K
power extension is user-authorized, not published NIST data. No liquid or
finite-density viscosity is represented here.
"""
import math
from numbers import Real

import numpy as np
from scipy.interpolate import CubicHermiteSpline, PPoly

VISCOSITY_SOURCE = "Hurly-Moldover 2000 NIST Table A1 (4He zero-density gas)"
VISCOSITY_SOURCE_URL = "https://nvlpubs.nist.gov/nistpubs/jres/105/5/j55hur.pdf"
VISCOSITY_MIN_K, VISCOSITY_TABLE_MIN_K, VISCOSITY_MAX_K = 0.5, 1.0, 350.0

# Temperature K, viscosity microPa s, temperature derivative microPa s/K.
NIST_TABLE = (
    (1.0, .3279, .05296), (1.2, .3395, .06980),
    (1.4, .3573, .1096), (1.6, .3834, .1504),
    (1.8, .4171, .1860), (2.0, .4573, .2149),
    (2.25, .5146, .2419), (2.5, .5775, .2600),
    (2.75, .6440, .2705), (3.0, .7123, .2749),
    (3.5, .8492, .2708), (4.0, .9815, .2573),
    (4.5, 1.106, .2407), (5.0, 1.222, .2244),
    (6.0, 1.433, .1974), (7.0, 1.620, .1783),
    (8.0, 1.791, .1649), (9.0, 1.951, .1551),
    (10.0, 2.102, .1474), (11.0, 2.246, .1412),
    (12.0, 2.385, .1359), (14.0, 2.648, .1274),
    (16.0, 2.895, .1206), (18.0, 3.131, .1151),
    (20.0, 3.356, .1104), (22.0, 3.573, .1064),
    (23.0, 3.678, .1046), (24.0, 3.782, .1029),
    (25.0, 3.884, .1013),
    # Continuation of Table A1, printed page 683 (PDF page 17).
    (26.0, 3.985, .09980), (28.0, 4.181, .09706),
    (30.0, 4.373, .09460), (35.0, 4.833, .08941),
    (40.0, 5.269, .08523), (45.0, 5.686, .08176),
    (50.0, 6.087, .07882), (60.0, 6.851, .07406),
    (70.0, 7.572, .07035), (80.0, 8.260, .06734),
    (90.0, 8.920, .06483), (100.0, 9.558, .06270),
    (120.0, 10.775, .05923), (140.0, 11.932, .05650),
    (160.0, 13.039, .05428), (180.0, 14.105, .05242),
    (200.0, 15.137, .05083), (225.0, 16.386, .04914),
    (250.0, 17.596, .04769), (275.0, 18.772, .04644),
    (300.0, 19.919, .04534), (325.0, 21.040, .04436),
    (350.0, 22.138, .04348),
)
_INTERPOLANT = CubicHermiteSpline(
    [r[0] for r in NIST_TABLE], [r[1] * 1e-6 for r in NIST_TABLE],
    [r[2] * 1e-6 for r in NIST_TABLE], extrapolate=False,
)
# Preserve even last-bit rounding at the old 5 K endpoint. Cubic Hermite
# intervals are local, but appending a knot changes the endpoint's selected
# interval from the left to the right. Use the original spline on 1-5 K.
_COLD_INTERPOLANT = CubicHermiteSpline(
    [r[0] for r in NIST_TABLE[:14]], [r[1] * 1e-6 for r in NIST_TABLE[:14]],
    [r[2] * 1e-6 for r in NIST_TABLE[:14]], extrapolate=False,
)
EXTRAPOLATION_EXPONENT = NIST_TABLE[0][0] * NIST_TABLE[0][2] / NIST_TABLE[0][1]

# On interval j, the spline uses z = T - T_j. Multiply its descending-power
# coefficients by T = z + T_j to obtain the quartic mu(T)*T. PPoly integrates
# each polynomial analytically and sums all crossed intervals, without an
# assumed average temperature or numerical quadrature grid.
_weighted_coefficients = np.zeros((5, len(NIST_TABLE) - 1))
_weighted_coefficients[:-1] += _INTERPOLANT.c
_weighted_coefficients[1:] += _INTERPOLANT.c * _INTERPOLANT.x[:-1]
_WEIGHTED_INTERPOLANT = PPoly(
    _weighted_coefficients, _INTERPOLANT.x, extrapolate=False,
)


def _validate_temperature(T):
    if isinstance(T, bool) or not isinstance(T, Real) or not math.isfinite(T):
        raise ValueError("Viscosity temperature must be a finite number.")
    if not VISCOSITY_MIN_K <= T <= VISCOSITY_MAX_K:
        raise ValueError("Gas viscosity supports 0.5-350 K; 0.5-1 K is extrapolated.")


def gas_viscosity(T):
    """Dynamic GAS viscosity [Pa s]; Hermite 1-350 K, power below 1 K."""
    _validate_temperature(T)
    if T < VISCOSITY_TABLE_MIN_K:
        return NIST_TABLE[0][1] * 1e-6 * (T / VISCOSITY_TABLE_MIN_K)**EXTRAPOLATION_EXPONENT
    if T <= 5.0:
        return float(_COLD_INTERPOLANT(T))
    return float(_INTERPOLANT(T))


def viscosity_temperature_integral(T_low, T_high):
    """Return signed integral of mu(T)*T dT [Pa s K^2] on 0.5-350 K.

    Equal limits give zero. The part below 1 K is the analytic integral of
    the existing power extension; the remainder is the exact integral of
    the piecewise Hermite polynomials, not an extrapolation to room temperature.
    For a linear profile the mean mu*T is this integral / (T_high-T_low).
    """
    _validate_temperature(T_low)
    _validate_temperature(T_high)
    if T_low == T_high:
        return 0.0
    if T_high < T_low:
        return -viscosity_temperature_integral(T_high, T_low)
    integral = 0.0
    if T_low < VISCOSITY_TABLE_MIN_K:
        upper = min(T_high, VISCOSITY_TABLE_MIN_K)
        exponent = EXTRAPOLATION_EXPONENT + 2.0
        coefficient = (NIST_TABLE[0][1] * 1e-6 /
                       VISCOSITY_TABLE_MIN_K**EXTRAPOLATION_EXPONENT)
        # expm1/log1p also preserve accuracy for nearly equal cold limits.
        power_difference = T_low**exponent * math.expm1(
            exponent * math.log1p((upper - T_low) / T_low))
        integral += coefficient * power_difference / exponent
    if T_high > VISCOSITY_TABLE_MIN_K:
        integral += float(_WEIGHTED_INTERPOLANT.integrate(
            max(T_low, VISCOSITY_TABLE_MIN_K), T_high))
    return integral
