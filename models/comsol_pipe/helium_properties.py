"""He-4 GAS properties for the steady COMSOL pipe model (SI units).

Viscosity is imported, unchanged, from the existing ``one_k_pot`` model.
Density and cp use that model's ideal-gas constants and the earlier cryostat
model's monatomic cp closure.  The primary conductivity is the companion
zero-density He-4 table of Hurly & Moldover (2000), Table A1, pp. 682-683:
https://nvlpubs.nist.gov/nistpubs/jres/105/5/j55hur.pdf
The optional ``legacy`` conductivity calls the previous cryostat Python
model's engineering anchor interpolation instead; it is not NIST data.

These dilute-gas material coefficients DO NOT validate no-slip continuum
transport in a rarefied tube.  Density/cp are ideal-gas approximations;
phase checks and Knudsen screens are independent of the material table.
No liquid properties, phase change, slip or wall accommodation are included.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.interpolate import CubicHermiteSpline

from models.one_k_pot.he4_viscosity import (
    NIST_TABLE as NIST_VISCOSITY_TABLE,
    gas_viscosity,
)
from models.one_k_pot.properties import saturation_pressure, saturation_temperature
M_HE4 = 0.004002602  # kg/mol; same constant as the original model.
R = 8.314462618  # J/(mol K)


MOLAR_MASS_KG_MOL = M_HE4
SPECIFIC_GAS_CONSTANT_J_KG_K = R / M_HE4
CP_J_KG_K = 2.5 * R / M_HE4
CV_J_KG_K = 1.5 * R / M_HE4
GAMMA = 5.0 / 3.0
MIN_T_K, MAX_T_K = 1.0, 350.0
NIST_SOURCE_URL = "https://nvlpubs.nist.gov/nistpubs/jres/105/5/j55hur.pdf"
RAREFACTION_SOURCE_URL = (
    "https://doc.comsol.com/6.3/doc/com.comsol.help.mfl/"
    "mfl_ug_modeling.05.18.html"
)

# T [K], k [mW/(m K)], dk/dT [mW/(m K^2)]. Table A1, He-4, not He-3.
NIST_CONDUCTIVITY_TABLE = (
    (1.0, 2.624, .4512), (1.2, 2.713, .5030),
    (1.4, 2.839, .7759), (1.6, 3.026, 1.097),
    (1.8, 3.276, 1.395), (2.0, 3.581, 1.642),
    (2.25, 4.022, 1.872), (2.5, 4.510, 2.023),
    (2.75, 5.028, 2.108), (3.0, 5.560, 2.143),
    (3.5, 6.628, 2.110), (4.0, 7.658, 2.005),
    (4.5, 8.630, 1.878), (5.0, 9.537, 1.754),
    (6.0, 11.180, 1.547), (7.0, 12.650, 1.400),
    (8.0, 14.000, 1.295), (9.0, 15.250, 1.218),
    (10.0, 16.440, 1.157), (11.0, 17.570, 1.108),
    (12.0, 18.660, 1.066), (14.0, 20.720, .9990),
    (16.0, 22.660, .9458), (18.0, 24.510, .9021),
    (20.0, 26.270, .8654), (22.0, 27.970, .8338),
    (23.0, 28.800, .8196), (24.0, 29.610, .8063),
    (25.0, 30.410, .7938), (26.0, 31.200, .7821),
    (28.0, 32.740, .7606), (30.0, 34.240, .7412),
    (35.0, 37.840, .7004), (40.0, 41.260, .6675),
    (45.0, 44.530, .6402), (50.0, 47.670, .6171),
    (60.0, 53.650, .5798), (70.0, 59.290, .5507),
    (80.0, 64.680, .5270), (90.0, 69.846, .5073),
    (100.0, 74.833, .4906), (120.0, 84.360, .4633),
    (140.0, 93.405, .4419), (160.0, 102.063, .4245),
    (180.0, 110.403, .4099), (200.0, 118.474, .3975),
    (225.0, 128.241, .3842), (250.0, 137.701, .3729),
    (275.0, 146.897, .3630), (300.0, 155.863, .3544),
    (325.0, 164.625, .3467), (350.0, 173.205, .3398),
)
_K_SPLINE = CubicHermiteSpline(
    [row[0] for row in NIST_CONDUCTIVITY_TABLE],
    [row[1] * 1e-3 for row in NIST_CONDUCTIVITY_TABLE],
    [row[2] * 1e-3 for row in NIST_CONDUCTIVITY_TABLE],
    extrapolate=False,
)


def _temperature(temperature_K, maximum_K=MAX_T_K):
    result = np.asarray(temperature_K, dtype=float)
    if np.any(~np.isfinite(result)) or np.any(result < MIN_T_K) or np.any(result > maximum_K):
        raise ValueError(f"Gas properties require {MIN_T_K:g} <= T <= {maximum_K:g} K.")
    return result


def _positive(value, name):
    result = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(result)) or np.any(result <= 0):
        raise ValueError(f"{name} must be finite and positive.")
    return result


def _result(value):
    result = np.asarray(value)
    return float(result) if result.ndim == 0 else result


def gas_thermal_conductivity(temperature_K, source="nist"):
    """k [W/(m K)], NIST Hermite 1-350 K; optional prior-model legacy 1-300 K.

    The legacy option calls its original log-log interpolation unchanged, but
    this wrapper rejects out-of-range temperatures instead of silently clipping.
    """
    if source == "nist":
        return _result(_K_SPLINE(_temperature(temperature_K)))
    if source == "legacy":
        from .legacy_properties import kappa_gas

        return _result(kappa_gas(_temperature(temperature_K, 300.0)))
    raise ValueError("Conductivity source must be 'nist' or 'legacy'.")


def gas_density(temperature_K, pressure_Pa):
    """Ideal He-4 gas density [kg/m^3]; pressure is ABSOLUTE."""
    return _result(_positive(pressure_Pa, "Absolute pressure") /
                   (SPECIFIC_GAS_CONSTANT_J_KG_K * _temperature(temperature_K)))


def gas_heat_capacity(temperature_K):
    """Ideal monatomic constant-pressure heat capacity [J/(kg K)]."""
    return _result(np.full_like(_temperature(temperature_K), CP_J_KG_K))


def mean_free_path(temperature_K, pressure_Pa):
    """Viscosity-based kinetic-theory regime estimate [m], reused from return_line."""
    temperatures, pressures = np.broadcast_arrays(
        _temperature(temperature_K), _positive(pressure_Pa, "Absolute pressure")
    )
    viscosities = np.array([gas_viscosity(float(t)) for t in temperatures.ravel()])
    viscosities = viscosities.reshape(temperatures.shape)
    return _result(viscosities / pressures * np.sqrt(
        math.pi * SPECIFIC_GAS_CONSTANT_J_KG_K * temperatures / 2.0
    ))


def knudsen_number(temperature_K, pressure_Pa, inner_diameter_m):
    """Kn based on tube bore; independent screen, not a rarefaction correction."""
    return _result(mean_free_path(temperature_K, pressure_Pa) /
                   _positive(inner_diameter_m, "Inner diameter"))


def conductivity_piecewise_rows(argument="t"):
    """Exact NIST k Hermite intervals for a COMSOL Piecewise function.

    Set COMSOL argunit='K', fununit='W/(m*K)', and extrapolation='none'.
    Rows are lower bound, upper bound, polynomial in numerical argument t.
    No smoothing should be added: the published-derivative Hermite function
    is already C1 continuous.  Function arguments use their stated K unit.
    """
    rows = []
    for j, (lo, hi) in enumerate(zip(_K_SPLINE.x[:-1], _K_SPLINE.x[1:])):
        z = f"({argument}-({lo:.17g}))"
        a, b, c, d = _K_SPLINE.c[:, j]
        expression = f"(({a:.17g}*{z}+({b:.17g}))*{z}+({c:.17g}))*{z}+({d:.17g})"
        rows.append([f"{lo:.17g}", f"{hi:.17g}", expression])
    return rows


def transport_table(temperatures_K=None, conductivity_source="nist"):
    """Rows [T/K, viscosity/(Pa*s), k/(W/(m*K))] for COMSOL import.

    Default sampling subdivides every NIST interval, preserving all original
    knots. COMSOL should use exact Piecewise polynomials where possible.
    """
    if temperatures_K is None:
        knots = np.array([row[0] for row in NIST_CONDUCTIVITY_TABLE])
        if conductivity_source == "legacy":
            knots = knots[knots <= 300]
        temperatures_K = np.unique(np.concatenate([
            np.linspace(a, b, 9) for a, b in zip(knots[:-1], knots[1:])
        ]))
    temperatures = np.atleast_1d(_temperature(temperatures_K))
    viscosities = np.array([gas_viscosity(float(t)) for t in temperatures])
    return np.column_stack((temperatures, viscosities,
                            gas_thermal_conductivity(temperatures, conductivity_source)))


def boundary_diagnostics(T1=1.0, P1=15.0, T2=293.15, P2=0.5,
                         inner_diameter_m=0.004826):
    """Default endpoint regime and cold-phase screens, all SI units."""
    result = {
        "inner_diameter_m": float(inner_diameter_m),
        "cp_J_kg_K": CP_J_KG_K,
        "inlet_mean_free_path_m": mean_free_path(T1, P1),
        "outlet_mean_free_path_m": mean_free_path(T2, P2),
        "inlet_Kn": knudsen_number(T1, P1, inner_diameter_m),
        "outlet_Kn": knudsen_number(T2, P2, inner_diameter_m),
        "inlet_density_kg_m3": gas_density(T1, P1),
        "outlet_density_kg_m3": gas_density(T2, P2),
    }
    if 0.8 <= T1 <= 4.3:
        psat = saturation_pressure(T1)
        result.update(inlet_saturation_pressure_Pa=psat,
                      inlet_pressure_over_saturation=P1 / psat,
                      inlet_condensation_warning=bool(P1 > psat))
    result["continuum_warning"] = bool(max(result["inlet_Kn"], result["outlet_Kn"]) >= .01)
    return result


PROPERTY_PROVENANCE = {
    "viscosity": {
        "implementation": "models.one_k_pot.he4_viscosity.gas_viscosity (unchanged import)",
        "source": "Hurly & Moldover 2000, Table A1, pp. 682-683; zero-density He-4 gas",
        "url": NIST_SOURCE_URL,
        "table_range_K": [1.0, 350.0],
    },
    "conductivity_nist": {
        "source": "Hurly & Moldover 2000, Table A1, pp. 682-683; zero-density He-4 gas",
        "url": NIST_SOURCE_URL,
        "interpolation": "Cubic Hermite with published temperature derivatives",
        "implemented_range_K": [1.0, 350.0],
    },
    "conductivity_legacy": {
        "implementation": "models.comsol_pipe.legacy_properties.kappa_gas (preserved interpolation)",
        "source": "Prior Python model engineering anchors; no primary citation in that code",
        "implemented_range_K": [1.0, 300.0],
    },
    "density_and_cp": {
        "density": "rho = P_abs*M/(R*T); inherited from one_k_pot.return_line",
        "cp": "cp = (5/2)*R/M; inherited monatomic closure from cryostat.heat_leak",
        "R_J_mol_K": R,
        "M_kg_mol": M_HE4,
        "limitations": "Ideal gas, not HEPAK or a real-fluid equation of state",
    },
    "phase_screen": {
        "implementation": "models.one_k_pot.properties.saturation_pressure (unchanged import)",
        "source": "Donnelly & Barenghi 1998, Table 17.2",
        "url": "https://doi.org/10.1063/1.556028",
        "author_table_url": "https://pages.uoregon.edu/rjd/vapor18.htm",
    },
    "flow_validity": {
        "mean_free_path": "lambda = mu/P * sqrt(pi*R*T/(2*M)); inherited return_line estimate",
        "source": RAREFACTION_SOURCE_URL,
        "limitation": "Kn>=0.01 invalidates a no-slip continuum interpretation; Kn>0.1 requires kinetic treatment",
    },
}


if __name__ == "__main__":
    import json

    print(json.dumps({"provenance": PROPERTY_PROVENANCE,
                      "default_diagnostics": boundary_diagnostics()}, indent=2))
