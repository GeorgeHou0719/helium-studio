"""Hydraulic iteration seed for COMSOL coupling; never a final thermal result.

Reuse COMSOL's prescribed-wall/conductive initial temperature guess to estimate
C = pi*d^4 / (256*R*integral(mu(T)*T dz)), mdot=C*(Pa^2-Pb^2).
The temperature guess does not solve an advection problem. Actual COMSOL flow
updates C in the coupling loop, and only actual COMSOL fields/heat fluxes are
accepted into final results. This is an initialization law, not a new pipe model.
"""
from functools import lru_cache
import math
import numpy as np
from .comsol_backend import normalize_parameters
from .comsol_worker import initial_temperature_profile
from models.comsol_pipe.helium_properties import gas_viscosity, SPECIFIC_GAS_CONSTANT_J_KG_K as R


@lru_cache(maxsize=256)
def _conductance(length, outer, wall, ta, tb, boundary, sections):
    inputs = {'length_m': length, 'temperature_a': ta, 'temperature_b': tb,
              'wall_boundary': boundary,
              'wall_profile': [dict(zip(('length_m', 'T_start_K', 'T_end_K'), row)) for row in sections]}
    x, temperature = initial_temperature_profile(inputs)
    integrand = np.array([gas_viscosity(float(t)) * t for t in temperature])
    resistance = length * float(np.trapezoid(integrand, x))
    coefficient = math.pi * (outer - 2 * wall)**4 / (256 * R * resistance)
    if not math.isfinite(coefficient) or coefficient <= 0:
        raise ValueError('Hydraulic initialization requires a finite positive pipe conductance.')
    return coefficient


def hydraulic_pipe_seed(params, pa, pb, ta, tb, *, profiles=False):
    data = normalize_parameters(params, pa, pb, ta, tb)
    sections = tuple(tuple(row[k] for k in ('length_m', 'T_start_K', 'T_end_K'))
                     for row in data['wall_profile'])
    coefficient = _conductance(data['length_m'], data['outer_diameter_m'],
                               data['wall_thickness_m'], data['temperature_a'],
                               data['temperature_b'], data['wall_boundary'], sections)
    mass = coefficient * (data['pressure_a'] - data['pressure_b']) * (data['pressure_a'] + data['pressure_b'])
    if not math.isfinite(mass):
        raise ValueError('Hydraulic initialization returned a nonfinite mass flow.')
    # Intentionally no thermal profiles or heat duties: this row is internal to
    # network iterations and must be replaced with a real COMSOL result.
    return {'mass_flow': mass, 'mass_flow_kg_s': mass,
            'backend': 'hydraulic_iteration_seed', 'warnings': []}
