"""Traceable AISI 304 stainless-steel properties for the steady pipe model.

Conductivity is the exact NIST Cryogenic Materials Property Database fit for
UNS S30400. NIST gives an equation range of 1--300 K, a data range of 4--300 K,
and a 2% curve-fit error relative to those data. The 1--4 K part is therefore
within the published equation range but below the underlying data range; 2%
must not be represented as its uncertainty.

This stationary, undeformed solid needs only conductivity in its heat-balance
equation. Density and the room-temperature heat-capacity constant below are
provided only to satisfy material-input fields if COMSOL requires them. The
constant Cp is NOT a cryogenic heat-capacity model and must not be used in a
transient study. A separately bounded NIST Cp function is included for 4--300 K.
"""

from __future__ import annotations

import numpy as np


NIST_SOURCE_URL = (
    "https://trc.nist.gov/cryogenics/materials/304Stainless/304Stainless_rev.htm"
)
DENSITY_SOURCE_URL = (
    "https://www.outokumpu.com/en/products/product-ranges/-/media/files/"
    "products/core/outokumpu-core-range-datasheet.pdf"
)

# a, b, ..., i in log10(y) = sum(c_i * log10(T/K)**i).
CONDUCTIVITY_COEFFICIENTS = (
    -1.4087, 1.3982, 0.2543, -0.6260, 0.2334, 0.4256,
    -0.4658, 0.1650, -0.0199,
)
HEAT_CAPACITY_COEFFICIENTS = (
    22.0061, -127.5528, 303.647, -381.0098, 274.0328, -112.9212,
    24.7593, -2.239153, 0.0,
)
DENSITY_KG_M3 = 7900.0
ROOM_TEMPERATURE_CP_J_KG_K = 500.0

PROPERTY_METADATA = {
    "material": "AISI 304 stainless steel (UNS S30400)",
    "conductivity": {
        "source": "NIST Cryogenic Materials Property Database, 304 Stainless",
        "source_url": NIST_SOURCE_URL,
        "unit": "W/(m K)",
        "equation": "log10(k) = sum(c_i * log10(T/K)**i), i=0..8",
        "coefficients": list(CONDUCTIVITY_COEFFICIENTS),
        "equation_range_K": [1.0, 300.0],
        "data_range_K": [4.0, 300.0],
        "curve_fit_error_percent_relative_to_data": 2.0,
        "low_temperature_note": (
            "1--4 K is inside NIST's published equation range but below its "
            "data range. NIST's 2% fit error is not a validated uncertainty "
            "in that interval. No additional extrapolation is introduced."
        ),
    },
    "density": {
        "value_kg_m3": DENSITY_KG_M3,
        "source": "Outokumpu Core range datasheet, physical properties Table 7",
        "source_url": DENSITY_SOURCE_URL,
        "source_value_kg_dm3": 7.9,
        "assumption": (
            "Constant room-temperature value; no thermal contraction modeled. "
            "It does not affect stationary conduction in an immobile solid."
        ),
    },
    "stationary_material_cp_placeholder": {
        "value_J_kg_K": ROOM_TEMPERATURE_CP_J_KG_K,
        "source": "Outokumpu Core range datasheet, Table 7, value at 20 degC",
        "source_url": DENSITY_SOURCE_URL,
        "limitation": (
            "Material-field placeholder only, absent from the stationary solid "
            "energy equation. Not a cryogenic Cp model; do not use transiently."
        ),
    },
    "heat_capacity_fit": {
        "source_url": NIST_SOURCE_URL,
        "unit": "J/(kg K)",
        "data_and_equation_range_K": [4.0, 300.0],
        "curve_fit_error_percent_relative_to_data": 5.0,
    },
}


def _fit(temperature_K, coefficients, lower_K: float):
    """Evaluate a published fit while refusing values outside its range."""
    temperature = np.asarray(temperature_K, dtype=float)
    if np.any(~np.isfinite(temperature)) or np.any(
        (temperature < lower_K) | (temperature > 300.0)
    ):
        raise ValueError(f"Temperature must be finite and in [{lower_K}, 300] K.")
    value = np.power(
        10.0, np.polynomial.polynomial.polyval(np.log10(temperature), coefficients)
    )
    return float(value) if value.ndim == 0 else value


def thermal_conductivity_w_m_k(temperature_K):
    """NIST 304 conductivity in W/(m K), 1 <= T <= 300 K."""
    return _fit(temperature_K, CONDUCTIVITY_COEFFICIENTS, 1.0)


def heat_capacity_j_kg_k(temperature_K):
    """NIST 304 Cp in J/(kg K), 4 <= T <= 300 K; no low-T extrapolation."""
    return _fit(temperature_K, HEAT_CAPACITY_COEFFICIENTS, 4.0)


def conductivity_table(n: int = 241) -> np.ndarray:
    """Return log-spaced [T/K, k/(W/m/K)] rows including 1, 4, 293.15, 300 K."""
    if isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < 2:
        raise ValueError("n must be an integer of at least 2.")
    temperature = np.unique(np.r_[np.geomspace(1.0, 300.0, n), 4.0, 293.15])
    return np.column_stack((temperature, thermal_conductivity_w_m_k(temperature)))


def comsol_conductivity_expression(temperature_argument: str = "T/1[K]") -> str:
    """Exact dimensioned COMSOL expression; argument must be dimensionless T/K.

    Range validation is the caller's responsibility. This form deliberately
    does not clamp or extend the fit beyond its published 1--300 K range.
    """
    log_temperature = f"log10({temperature_argument})"
    polynomial = f"({CONDUCTIVITY_COEFFICIENTS[-1]:.8g})"
    for coefficient in reversed(CONDUCTIVITY_COEFFICIENTS[:-1]):
        polynomial = f"(({coefficient:.8g})+({log_temperature})*{polynomial})"
    return f"10^({polynomial})*1[W/(m*K)]"
