"""Adiabatic HEM capillary: friction, acceleration, kinetic energy and choking.

P is the integration coordinate. h + G**2*v**2/2 = h0. Available impedance
is the integral from Pexit to Pentry of
 (1 + G**2 * dv/dP|h0)/(mdot * mu * v * C_f), C_f = f_D*Re/64.
Integration stops at the sonic condition (numerator zero), independently of
lower pot pressure. z_equiv must equal measured/geometric Poiseuille impedance.
The outlet jet equilibrates in the pot at h0 after kinetic-energy recovery.
This is a homogeneous equilibrium approximation, not a He II two-fluid model.
"""
import math
import numpy as np
from scipy.optimize import brentq
from .two_phase_properties import Helium, PhaseDomainError


def darcy_factor(reynolds, relative_roughness=0.):
    """Churchill (1977) all-regime Darcy factor; 64/Re in laminar limit."""
    re=max(float(reynolds),1e-12)
    if re < 100: return 64/re
    a=(2.457*math.log(1/((7/re)**.9+.27*relative_roughness)))**16
    b=(37530/re)**16
    return 8*((8/re)**12+1/(a+b)**1.5)**(1/12)


def impedance(params):
    d=float(params['diameter_m']); length=float(params['length_m'])
    if d<=0 or length<=0:
        raise ValueError('Capillary diameter and length must be positive.')
    if params.get('law','capillary') != 'capillary':
        raise ValueError('The two-phase mode currently requires a capillary restriction.')
    z=float(params['impedance_m_inv3']) if params.get('impedance_mode')=='measured' else 128*length/(math.pi*d**4)
    if not math.isfinite(z) or z<=0:
        raise ValueError('Capillary impedance must be finite and positive.')
    return d,length,z


def path_at_flow(params, pressure_in, enthalpy_in, pressure_back, mass_flow, *, cells=160, helium=None):
    """Return impedance traversed before back pressure, choking or property limit.

    Inlet enthalpy is the upstream reservoir/stagnation enthalpy. This includes
    entrance acceleration without an invented entrance-loss coefficient.
    """
    he=helium or Helium();d,length,target=impedance(params)
    if not 0 < pressure_back < pressure_in or mass_flow<=0:
        raise ValueError('Forward capillary flow requires Pin > Pback > 0 and positive mass flow.')
    area=math.pi*d*d/4; flux=mass_flow/area
    upstream=he.ph(pressure_in,enthalpy_in)
    def entrance_energy(p):
        state=he.ps(float(p),upstream.entropy)
        return state.enthalpy+.5*(flux*state.volume)**2-enthalpy_in
    p_high=pressure_in;previous=entrance_energy(p_high)
    entrance_pressure=None
    for p_low in np.geomspace(pressure_in,pressure_back,max(cells,160)+1)[1:]:
        try:current=entrance_energy(p_low)
        except PhaseDomainError:
            return dict(impedance=0.,choked=False,domain_limited=True,profile=[],mass_flow=mass_flow)
        if current<=0:
            entrance_pressure=brentq(entrance_energy,float(p_low),p_high,xtol=1e-7)
            break
        if current>previous*(1+1e-7):break
        p_high=float(p_low);previous=current
    if entrance_pressure is None:
        # Requested flux exceeds the lossless entrance/nozzle capacity or
        # consumes more pressure than is available before the pipe begins.
        return dict(impedance=0.,choked=True,domain_limited=False,profile=[],mass_flow=mass_flow)
    points=[]; blocked=False
    # Avoid calling the normal-liquid EOS deep into He II. A trial that reaches
    # its boundary before choking remains explicitly incomplete.
    for p in np.geomspace(entrance_pressure, pressure_back, cells+1):
        def energy(h):
            state=he.ph(float(p),float(h))
            return h + .5*(flux*state.volume)**2 - enthalpy_in
        try:
            upper=enthalpy_in
            high=he.ph(float(p),upper)
            delta=max(.55*(flux*high.volume)**2,1.)
            lower=upper-delta
            for _ in range(32):
                try:
                    low_value=energy(lower)
                except PhaseDomainError:
                    # At low pressure the lowest gas state is Tmin; at higher
                    # pressure this is the normal-liquid boundary.
                    lower=he.pt(float(p),he.minimum_temperature).enthalpy+1e-5
                    low_value=energy(lower)
                    if low_value>0:
                        raise PhaseDomainError('Kinetic-energy trajectory leaves the normal-fluid domain.')
                    break
                if low_value<=0:break
                delta*=2;lower=upper-delta
            h=brentq(energy,lower,upper,xtol=1e-6)
            state=he.ph(float(p),h)
        except PhaseDomainError:
            blocked=True
            break
        points.append(state)
        if len(points)>=3:
            # A sign change means the last interval has passed the critical
            # point. Enough information is now available for interpolation.
            v1,v2=points[-2].volume,points[-1].volume
            if 1+flux*flux*(v2-v1)/(points[-1].pressure-points[-2].pressure)<0:
                break
    if len(points)<3:
        return dict(impedance=0.,choked=True,domain_limited=blocked,profile=[],
                    mass_flow=mass_flow,critical_pressure_Pa=pressure_in)
    p=np.array([s.pressure for s in points]);v=np.array([s.volume for s in points])
    dv=np.gradient(v,p,edge_order=2)
    numerator=1+flux*flux*dv
    cross=np.flatnonzero(numerator<=0)
    z_values=[0.]; used=[points[0]]; choked=bool(len(cross))
    first=int(cross[0]) if choked else len(points)
    if first==0:
        return dict(impedance=0.,choked=True,domain_limited=False,profile=[],
                    mass_flow=mass_flow,critical_pressure_Pa=pressure_in)
    def integrand(i):
        re=4*mass_flow/(math.pi*d*points[i].viscosity)
        multiplier=darcy_factor(re,float(params.get('roughness_m',0.))/d)*re/64
        return numerator[i]/(mass_flow*points[i].viscosity*v[i]*multiplier)
    for i in range(1,first):
        z_values.append(z_values[-1]+.5*(integrand(i-1)+integrand(i))*(p[i-1]-p[i]))
        used.append(points[i])
    if choked:
        i=first
        fraction=numerator[i-1]/(numerator[i-1]-numerator[i])
        pc=float(p[i-1]+fraction*(p[i]-p[i-1]))
        z_values.append(z_values[-1]+.5*integrand(i-1)*(p[i-1]-pc))
        ha,hb=sorted([used[-1].enthalpy,points[i].enthalpy])
        hc=brentq(lambda h:h+.5*(flux*he.ph(pc,h).volume)**2-enthalpy_in,ha,hb,xtol=1e-7)
        used.append(he.ph(pc,hc))
    else:pc=float(p[-1])
    z_final=float(z_values[-1])
    profile=[s.point(length*z/target,mass_flow,area) for s,z in zip(used,z_values)]
    re=max(4*mass_flow/(math.pi*d*s.viscosity) for s in used)
    return dict(impedance=z_final,choked=choked,domain_limited=blocked and not choked,
                critical_pressure_Pa=pc if choked else None,exit_pressure_Pa=pc,
                entrance_pressure_Pa=entrance_pressure,upstream_entropy_J_kg_K=upstream.entropy,
                mass_flow=mass_flow,profile=profile,max_Re=re,
                total_enthalpy_J_kg=enthalpy_in,
                relative_impedance_residual=z_final/target-1)


def solve_capillary(params,pressure_in,enthalpy_in,pressure_back,*,cells=240):
    he=Helium();d,length,z=impedance(params)
    upstream=he.ph(pressure_in,enthalpy_in)
    initial=(pressure_in-pressure_back)/(z*upstream.viscosity*upstream.volume)
    def trial(m):
        return path_at_flow(params,pressure_in,enthalpy_in,pressure_back,m,cells=cells,helium=he)
    high=max(initial,1e-15)
    for _ in range(40):
        r=trial(high)
        if r['impedance']<z:break
        high*=2
    low=high/2
    for _ in range(60):
        r=trial(low)
        if r['impedance']>z:break
        if r['domain_limited']:
            raise PhaseDomainError('The capillary enters He II before a flow solution or choking point. Public equilibrium data cannot supply its superfluid hydraulic resistance.')
        low/=2
    else:raise ValueError('Could not bracket the capillary flow.')
    m=brentq(lambda m:trial(m)['impedance']/z-1,low,high,xtol=1e-16,rtol=1e-7)
    r=trial(m)
    if r['domain_limited']:
        raise PhaseDomainError('The capillary reaches He II before choking. A superfluid flow model is required.')
    r['mass_flow_kg_s']=m
    return r
