"""Fast 1-D conjugate-equilibrium pipe preview, distinct from COMSOL.

Cross sections share T for gas and wall. K(T)=As*ks+Ag*kg, and stationary
conduction/advection gives K(T) T' = mdot*cp*T + C. This conserves total
enthalpy plus conductive flux. Hydraulic closure is compressible, creeping,
no-slip Poiseuille, using the solved temperature profile. No inertia/phase
change/rarefaction or Joule-Thomson cooling is manufactured in this preview.
"""
from __future__ import annotations
from functools import lru_cache
import math
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import brentq
from models.comsol_pipe.helium_properties import (
    gas_viscosity, gas_thermal_conductivity, SPECIFIC_GAS_CONSTANT_J_KG_K as R,
    CP_J_KG_K as CP,
)
from models.comsol_pipe.steel_properties import thermal_conductivity_w_m_k

_TEMPS=np.unique(np.r_[np.geomspace(1,12,900),np.linspace(12,300,1800)])
_MU=np.array([gas_viscosity(float(t)) for t in _TEMPS])
_KG=gas_thermal_conductivity(_TEMPS)
_KS=thermal_conductivity_w_m_k(_TEMPS)
_GX,_GW=np.polynomial.legendre.leggauss(80)


def _number(value,name,lo,hi):
    if isinstance(value,bool): raise ValueError(f'{name} must be a number')
    v=float(value)
    if not math.isfinite(v) or not lo<=v<=hi: raise ValueError(f'{name} must be within {lo:g}–{hi:g} range')
    return v


def _temperature(value, name):
    if isinstance(value, bool):
        raise ValueError(f'{name} must be a number')
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f'{name} must be finite and greater than zero K absolute')
    return value


def dimensions(params):
    length=_number(params.get('length_m',.288),'Pipe length (m)',1e-5,100)
    outer=_number(params.get('outer_diameter_m',.00635),'Outer diameter (m)',1e-5,1)
    wall=_number(params.get('wall_thickness_m',.000762),'Wall thickness (m)',1e-7,.49*outer)
    return length,outer,wall


def wall_segments(params,ta,tb):
    """Validate a continuous, prescribed outer-wall temperature schedule."""
    length,_,_=dimensions(params)
    profile=params.get('wall_profile',[])
    if not isinstance(profile,list) or not 1<=len(profile)<=32:
        raise ValueError('A prescribed outer-wall temperature requires 1–32 temperature segments')
    knots=[0.]; temperatures=[]; previous=None
    for segment in profile:
        le=_number(segment['length_m'],'Temperature segment length',1e-6,100)
        start=_temperature(segment['T_start_K'],'Outer-wall temperature')
        end=_temperature(segment['T_end_K'],'Outer-wall temperature')
        if previous is not None and abs(start-previous)>1e-6:
            raise ValueError('Adjacent wall-temperature segments must be continuous')
        if not temperatures: temperatures.append(start)
        knots.append(knots[-1]+le);temperatures.append(end);previous=end
    if abs(knots[-1]-length)>max(1e-8,length*1e-6):
        raise ValueError('Outer-wall segment lengths must sum to the pipe length')
    if abs(temperatures[0]-ta)>1e-6 or abs(temperatures[-1]-tb)>1e-6:
        raise ValueError('The first and last wall temperatures must match the pipe end temperatures')
    return np.array(knots),np.array(temperatures)


def _prescribed(params,pa,pb,ta,tb):
    length,outer,wall=dimensions(params);di=outer-2*wall
    ag=math.pi*di**2/4;ass=math.pi*(outer**2-di**2)/4
    knots,temperatures=wall_segments(params,ta,tb)
    # Include each slope transition exactly; enough samples for mu(T)*T.
    x=np.unique(np.r_[np.linspace(0,length,401),knots])
    temp=np.interp(x,knots,temperatures)
    cumulative=cumulative_trapezoid(np.interp(temp,_TEMPS,_MU)*temp,x,initial=0)
    mass=math.pi*di**4/(256*R)*(pa-pb)*(pa+pb)/cumulative[-1]
    pressure=np.sqrt(np.maximum(1e-30,pa*pa-(pa*pa-pb*pb)*cumulative/cumulative[-1]))
    pressure[0]=pa;pressure[-1]=pb
    qa=float(_k(ta,ag,ass)*(temperatures[1]-ta)/(knots[1]-knots[0]))
    qb=-float(_k(tb,ag,ass)*(tb-temperatures[-2])/(knots[-1]-knots[-2]))
    return mass,x,temp,pressure,qa,qb


def _k(t,ag,ass):
    return ag*np.interp(t,_TEMPS,_KG)+ass*np.interp(t,_TEMPS,_KS)


def _temperature_profile(length,di,do,ta,tb,mass,n=101):
    """Solve the first-integral boundary value problem by a scalar root."""
    x=np.linspace(0,length,n)
    ag=math.pi*di**2/4; ass=math.pi*(do**2-di**2)/4
    if abs(tb-ta)<1e-10:
        return x,np.full(n,ta),0.,0.
    direction=1 if tb>ta else -1
    low=min(ta,tb); high=max(ta,tb); a=mass*direction*CP
    t=(high-low)/2*_GX+(high+low)/2
    weights=(high-low)/2*_GW
    kval=_k(t,ag,ass)
    q0=float(np.dot(weights,kval))/length
    lower=-min(a*low,a*high)
    eps=max(q0*1e-13,1e-18)
    def error(c): return float(np.dot(weights,kval/(a*t+c)))-length
    upper=lower+q0+abs(a)*(high-low)+eps
    # Extremely advective trials cannot be resolved by this gas/wall equilibrium
    # preview. Reject explicitly instead of silently inventing a profile.
    if error(lower+eps)<0:
        raise ValueError('Advection is too strong for the quick pipe model; Reduce the flow or use COMSOL coupled simulation')
    while error(upper)>0: upper=lower+2*(upper-lower)
    c=brentq(error,lower+eps,upper,xtol=max(q0*1e-12,1e-20),rtol=1e-12)
    # Integral in temperature space, with extra endpoint resolution.
    u=np.linspace(0,1,801)
    tt=low+(high-low)*(1-np.cos(np.pi*u))/2
    xx=cumulative_trapezoid(_k(tt,ag,ass)/(a*tt+c),tt,initial=0)
    xx*=length/xx[-1]
    temp=np.interp(x,xx,tt)
    if direction<0: temp=temp[::-1]
    flux_a=mass*CP*ta+direction*c
    flux_b=-(mass*CP*tb+direction*c)
    return x,temp,float(flux_a),float(flux_b)


@lru_cache(maxsize=2048)
def _solve(length,outer,wall,pa,pb,ta,tb):
    di=outer-2*wall
    coefficient=math.pi*di**4/(256*R)
    drive=(pa-pb)*(pa+pb)
    # Start from the nonlinear, zero-flow conductive temperature profile.
    x,temp,_,_=_temperature_profile(length,di,outer,ta,tb,0.)
    resistance=float(np.trapezoid(np.interp(temp,_TEMPS,_MU)*temp,x))
    mass=coefficient*drive/resistance
    for _ in range(12):
        x,temp,qa,qb=_temperature_profile(length,di,outer,ta,tb,mass)
        resistance=float(np.trapezoid(np.interp(temp,_TEMPS,_MU)*temp,x))
        new=coefficient*drive/resistance
        if abs(new-mass)<=1e-7*max(abs(new),1e-16):
            mass=new
            break
        mass=.7*new+.3*mass
    x,temp,qa,qb=_temperature_profile(length,di,outer,ta,tb,mass)
    weighted=np.interp(temp,_TEMPS,_MU)*temp
    cum=cumulative_trapezoid(weighted,x,initial=0)
    pressure=np.sqrt(np.maximum(1e-30,pa*pa-(pa*pa-pb*pb)*cum/cum[-1]))
    pressure[0]=pa; pressure[-1]=pb
    return mass,x,temp,pressure,qa,qb


def reduced_pipe(params,p_a,p_b,t_a,t_b,*,profiles=False):
    length,outer,wall=dimensions(params)
    pa=_number(p_a,'End pressure (Pa)',1e-8,1e9); pb=_number(p_b,'End pressure (Pa)',1e-8,1e9)
    ta=_temperature(t_a,'End temperature (K)'); tb=_temperature(t_b,'End temperature (K)')
    boundary=params.get('wall_boundary','adiabatic')
    if boundary not in ('adiabatic','prescribed'): raise ValueError('Unknown outer-wall boundary type')
    mass,x,temp,pressure,qa,qb=(_prescribed(params,pa,pb,ta,tb) if boundary=='prescribed'
                               else _solve(length,outer,wall,pa,pb,ta,tb))
    environment=-(qa+qb+mass*CP*(tb-ta)) if boundary=='prescribed' else 0.
    output={'mass_flow':mass,'mass_flow_kg_s':mass,'pressure_in':pa,'pressure_out':pb,
            'temperature_in':ta,'temperature_out':tb,'heat_a':qa,'heat_b':qb,
            'heat_a_W':qa,'heat_b_W':qb,'heat_environment_W':float(environment),'backend':'reduced',
            'warnings':[],'profile':[]}
    if profiles:
        di=outer-2*wall; area=math.pi*di**2/4
        mu=np.interp(temp,_TEMPS,_MU);rho=pressure/(R*temp)
        velocity=mass/(rho*area)
        kn=mu/pressure*np.sqrt(math.pi*R*temp/2)/di
        re=abs(mass)*di/(area*mu)
        mach=abs(velocity)/np.sqrt(5/3*R*temp)
        output['profile']=[{'x':float(z),'x_m':float(z),'temperature':float(t),'temperature_K':float(t),
                            'wall_temperature':float(t),'wall_temperature_K':float(t),'outer_wall_temperature':float(t),'pressure':float(p),'pressure_Pa':float(p),
                            'velocity':float(v),'velocity_m_s':float(v),'Kn':float(k)} for z,t,p,v,k in zip(x,temp,pressure,velocity,kn)]
        output['diagnostics']={'max_Kn':float(max(kn)),'max_Re':float(max(re)), 'max_Mach':float(max(mach)),
                               'energy_residual_W':float(qa+qb+environment+mass*CP*(tb-ta))}
        if max(kn)>=.01: output['warnings'].append(f'Rarefied gas: Maximum Kn={max(kn):.3g}, Use continuum results with caution')
        if max(mach)>=.3: output['warnings'].append(f'Maximum Mach={max(mach):.3g}, Outside the low-speed pipe-model range')
        if max(re)>=2300: output['warnings'].append(f'Maximum Re={max(re):.3g}, The laminar-flow assumption is invalid')
    return output
