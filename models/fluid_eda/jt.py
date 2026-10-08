"""Gas restriction flow laws and explicit ideal-gas isenthalpic closure.

Short orifice: compressible reservoir-to-orifice law with choking, derived
from NASA Glenn mass-flow equations; Cd is an empirical coefficient.
Long capillary: isothermal, fully developed compressible Poiseuille law.
Neither is a universal JT resistance. Both require continuum single gas.
Ideal h=cp*T gives Tout=Tin across a thermally insulated throttle: no claimed
real-helium liquefaction or JT cooling. See PHYSICS.md for source and bounds.
"""
import math
from models.comsol_pipe.helium_properties import gas_viscosity, SPECIFIC_GAS_CONSTANT_J_KG_K as R

GAMMA=5/3
CRITICAL_RATIO=(2/(GAMMA+1))**(GAMMA/(GAMMA-1))


def evaluate_restriction(params,p_a,p_b,t_a,t_b):
    vals=[float(p_a),float(p_b),float(t_a),float(t_b)]
    if any(not math.isfinite(v) or v<=0 for v in vals): raise ValueError('Restrictions require positive absolute pressures and temperatures')
    pa,pb,ta,tb=vals
    forward=pa>=pb; hi,lo=(pa,pb) if forward else (pb,pa); tin=ta if forward else tb
    d=float(params.get('diameter_m',1e-5)); length=float(params.get('length_m',.01)); cd=float(params.get('discharge_coefficient',.7))
    if not math.isfinite(d) or not 1e-8<=d<=.1: raise ValueError('Restriction diameter must be within 1e-8–0.1 m range')
    ratio=lo/hi;law=params.get('law','orifice'); choked=False; impedance=None
    if law=='orifice':
        if not math.isfinite(cd) or not 0<cd<=1: raise ValueError('Discharge coefficient Cd must be within (0,1] range')
        if ratio<=CRITICAL_RATIO:
            factor=math.sqrt(GAMMA)*(2/(GAMMA+1))**((GAMMA+1)/(2*(GAMMA-1))); choked=True
        else:
            factor=math.sqrt(max(0,2*GAMMA/(GAMMA-1)*(ratio**(2/GAMMA)-ratio**((GAMMA+1)/GAMMA))))
        magnitude=cd*math.pi*d*d/4*hi/math.sqrt(R*tin)*factor
    elif law=='capillary':
        if not math.isfinite(length) or not 1e-7<=length<=100: raise ValueError('Capillary length must be within 1e-7–100 m range')
        mode=params.get('impedance_mode','geometry')
        if mode=='geometry':
            impedance=128*length/(math.pi*d**4)
        elif mode=='measured':
            raw=params.get('impedance_m_inv3')
            if isinstance(raw,bool) or not isinstance(raw,(int,float)) or not math.isfinite(raw) or raw<=0:
                raise ValueError('Measured capillary impedance z must be finite and positive (m^-3)')
            impedance=float(raw)
        else:
            raise ValueError('Capillary impedance mode must be geometry or measured')
        # z = delta_P / (mu * volume_flow) for incompressible Poiseuille flow.
        # Isothermal ideal-gas integration gives mdot = (pa^2-pb^2)/(2*mu*R*T*z).
        magnitude=(hi-lo)*(hi+lo)/(2*gas_viscosity(tin)*R*tin*impedance)
    else: raise ValueError('Restriction flow law must be orifice or capillary')
    mass=magnitude if forward else -magnitude
    kn_hi=gas_viscosity(tin)/hi*math.sqrt(math.pi*R*tin/2)/d
    kn_lo=gas_viscosity(tin)/lo*math.sqrt(math.pi*R*tin/2)/d
    warnings=['Ideal-gas isenthalpic throttling: outlet temperature equals upstream temperature; Real-helium J–T cooling and two-phase flash are not calculated']
    if max(kn_hi,kn_lo)>=.01: warnings.append(f'Rarefied gas in restriction, Kn estimated up to {max(kn_hi,kn_lo):.3g}; Continuum flow resistance requires calibration')
    if law=='capillary' and ratio<.1: warnings.append('The capillary pressure ratio is large; The Poiseuille model excludes inertia and choking; separate validation is required')
    if law=='capillary' and params.get('impedance_mode')=='measured':
        warnings.append('Measured impedance is used as a fixed coefficient in the isothermal gas model; its transfer to cryogenic or two-phase conditions is not validated')
    return {'mass_flow':mass,'mass_flow_kg_s':mass,'pressure_in':pa,'pressure_out':pb,
            'temperature_in':ta,'temperature_out':tin,'outlet_fluid_temperature_K':tin,
            'upstream_temperature_K':tin,'heat_a':0.,'heat_b':0.,'heat_a_W':0.,'heat_b_W':0.,
            'profile':[{'x':0.,'pressure':pa,'temperature':tin},{'x':length if law=='capillary' else d,'pressure':pb,'temperature':tin}],
            'warnings':warnings,'diagnostics':{'choked':choked,'pressure_ratio':ratio,'critical_ratio':CRITICAL_RATIO,'max_Kn':max(kn_hi,kn_lo),
                **({'impedance_m_inv3':impedance,'impedance_mode':params.get('impedance_mode','geometry')} if law=='capillary' else {})},
            'material':params.get('material',''),
            'law':law,'thermodynamic_model':'ideal_gas_isenthalpic','backend':'analytic'}
