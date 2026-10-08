"""One-dimensional pipe closures for the experimental continuous-pot mode.

Suction: dilute vapor, separate gas/steel temperatures, axial steel conduction,
adiabatic exterior, finite Nu=3.66 gas-wall exchange, no intercomponent conduction.
Precooler: real-fluid P,h marching with finite wall exchange and latent heat.
The prescribed exterior is an ideal local thermal boundary. Its axial steel
conduction does not alter the imposed temperature profile. Fluid axial conduction
and axial kinetic energy in the large pipes are neglected (reported assumptions).
"""
from functools import lru_cache
import math
import numpy as np
from scipy.integrate import solve_bvp, cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from .thermal import dimensions, wall_segments
from .two_phase_properties import Helium, PhaseDomainError
from .two_phase_capillary import darcy_factor
from .pot_energy import feed_state
from models.continuous_he4 import properties as he4
from models.comsol_pipe.helium_properties import gas_viscosity, gas_thermal_conductivity, SPECIFIC_GAS_CONSTANT_J_KG_K as R
from models.comsol_pipe.steel_properties import thermal_conductivity_w_m_k as steel_k


@lru_cache(maxsize=24)
def gas_table(pot_temperature, pressure):
    # A dilute-vapor enthalpy table, matching the saturated vapor at the pot.
    # A small reference transition joins DB98 below lambda to normal-fluid HEOS.
    ts=np.unique(np.r_[np.linspace(pot_temperature,2.18,140) if pot_temperature<2.18 else [],np.geomspace(max(2.18,pot_temperature),500,350)])
    hs=np.array([feed_state(float(t),pressure*.999)[0] for t in ts])
    hs[0]=float(he4.saturated_vapor_enthalpy(pot_temperature))
    if np.any(np.diff(hs)<=0):
        raise ValueError('The joined dilute-vapor enthalpy table is not monotonic.')
    enthalpy=PchipInterpolator(ts,hs)
    return enthalpy,enthalpy.derivative()


def suction_pipe(params,mass,pressure_pot,temperature_pot,*,rtol=2e-4,initial_guess=None):
    le,outer,wall=dimensions(params);di=outer-2*wall
    if params.get('wall_boundary','adiabatic')!='adiabatic':
        raise ValueError('The continuous-pot mode requires an adiabatic suction-pipe exterior.')
    ta,tb=float(params['T_a_K']),float(params['T_b_K'])
    if not 1<=temperature_pot<4.3 or min(ta,tb)<=0 or max(ta,tb)>300 or mass<=0:
        raise PhaseDomainError('This 1D suction closure uses public gas transport from 1 K and the 304 steel fit through 300 K. Supply additional material data for colder/hotter wall boundaries; this is not a COMSOL temperature restriction.')
    if abs(ta-temperature_pot)>1e-6:
        raise ValueError('The suction inlet wall must be anchored at the pot trial temperature.')
    area=math.pi*di**2/4;solid=math.pi*(outer**2-di**2)/4
    enthalpy,cp=gas_table(temperature_pot,pressure_pot)
    factor=float(params.get('heat_transfer_multiplier',1.))
    if not math.isfinite(factor) or factor<=0:
        raise ValueError('Heat-transfer multiplier must be finite and positive.')
    # Kirchhoff transform of steel temperature avoids dividing the Newton
    # system by very small cryogenic conductivity. Gas enthalpy is the energy
    # unknown, so cp kinks do not introduce spurious energy sources.
    grid_t=np.unique(np.r_[np.geomspace(1,20,600),np.linspace(20,300,1000)])
    kw=cumulative_trapezoid(steel_k(grid_t)*solid,grid_t,initial=0)
    forward_k=PchipInterpolator(grid_t,kw)
    inverse_k=PchipInterpolator(kw,grid_t)
    ht=np.unique(np.r_[np.linspace(temperature_pot,2.18,300) if temperature_pot<2.18 else [],np.geomspace(max(2.18,temperature_pot),500,500)])
    inverse_h=PchipInterpolator(enthalpy(ht),ht)
    hscale=1e6
    def ode(x,y):
        tg=inverse_h(y[0]*hscale);tw=inverse_k(y[1])
        g=3.66*math.pi*gas_thermal_conductivity(np.clip(tg,1,350))*factor
        q=g*(tw-tg)
        return np.vstack([q/(mass*hscale),-y[2],-q])
    def bc(a,b):return np.array([a[0]-float(enthalpy(temperature_pot))/hscale,a[1]-float(forward_k(ta)),b[1]-float(forward_k(tb))])
    x=np.unique(np.r_[0.,np.geomspace(le*1e-7,le,100),np.linspace(0,le,100)])
    potential=float(forward_k(ta))+(float(forward_k(tb))-float(forward_k(ta)))*x/le
    temp=inverse_k(potential)
    initial=np.vstack([enthalpy(temp)/hscale,potential,np.full_like(x,-(float(forward_k(tb))-float(forward_k(ta)))/le)])
    if initial_guess is not None:
        x=np.asarray(initial_guess['x']);initial=np.asarray(initial_guess['y'])
    solution=solve_bvp(ode,bc,x,initial,tol=rtol,max_nodes=8000)
    if not solution.success:
        raise ValueError(f'The 1D suction energy solve did not converge at {mass:g} kg/s: '+solution.message)
    x=np.unique(np.r_[solution.x,np.linspace(0,le,241)])
    hg,potential,q=solution.sol(x)
    tg=inverse_h(hg*hscale);tw=inverse_k(potential)
    if min(tg)<temperature_pot-1e-5 or max(tg)>max(ta,tb)+1e-3 or min(tw)<=0:
        raise ValueError('The suction solution violates its thermal boundaries.')
    mu=np.array([gas_viscosity(float(t)) for t in tg])
    re=4*mass/(math.pi*di*mu)
    multiplier=np.array([darcy_factor(v)*v/64 for v in re])
    integral=cumulative_trapezoid(mu*tg*multiplier,x,initial=0)
    p2=pressure_pot**2-256*R*mass/(math.pi*di**4)*integral
    if p2[-1]<=0:
        raise ValueError('Requested circulation exceeds the suction-pipe pressure-drop capacity.')
    p=np.sqrt(p2);rho=p/(R*tg);velocity=mass/(area*rho)
    kn=mu/p*np.sqrt(math.pi*R*tg/2)/di
    mach=velocity/np.sqrt((5/3)*R*tg)
    if max(mach)>.3:
        raise ValueError('Suction-pipe Mach number exceeds 0.3; this dilute-vapor pipe closure needs compressible inertia.')
    h=hg*hscale
    total=mass*h+q
    residual=float((max(total)-min(total))/max(abs(mass*(h[-1]-h[0])),1e-12))
    if residual>5e-3:
        raise ValueError('The suction-pipe fluid/steel energy balance did not meet tolerance.')
    indices=np.unique(np.linspace(0,len(x)-1,180).astype(int))
    profile=[dict(x=float(x[i]),pressure=float(p[i]),temperature=float(tg[i]),wall_temperature=float(tw[i]),enthalpy_J_kg=float(h[i]),velocity=float(velocity[i]),phase='gas',vapor_quality=None) for i in indices]
    return dict(_bvp_guess=dict(x=solution.x.tolist(),y=solution.y.tolist()),mass_flow=mass,mass_flow_kg_s=mass,pressure_in=pressure_pot,pressure_out=float(p[-1]),
                temperature_in=temperature_pot,temperature_out=float(tg[-1]),
                enthalpy_in_J_kg=float(h[0]),enthalpy_out_J_kg=float(h[-1]),profile=profile,
                diagnostics=dict(max_Re=float(max(re)),max_Kn=float(max(kn)),max_Mach=float(max(mach)),relative_energy_residual=residual),
                warnings=['Dilute-vapor suction: finite gas/steel exchange and axial steel conduction; gas axial conduction and axial inertia are neglected.']+([f'Suction Kn reaches {max(kn):.3g}; continuum no-slip friction and heat exchange require rarefaction/accommodation calibration.'] if max(kn)>=.01 else []),backend='one_dimensional_conjugate_gas')


def precooling_pipe(params,mass,pressure_in,temperature_in,*,max_step_fraction=1/100):
    he=Helium();le,outer,wall=dimensions(params);di=outer-2*wall;area=math.pi*di*di/4
    if params.get('wall_boundary')!='prescribed':
        raise ValueError('The two-phase precooler requires its prescribed outer-wall temperature profile.')
    knots,temps=wall_segments(params,float(params['T_a_K']),float(params['T_b_K']))
    if min(temps)<he.minimum_temperature:
        raise PhaseDomainError('The precooler wall enters the He II range. Normal-fluid heat transfer must not be extrapolated into superfluid liquid.')
    factor=float(params.get('heat_transfer_multiplier',1.))
    if not math.isfinite(factor) or factor<=0:
        raise ValueError('Heat-transfer multiplier must be finite and positive.')
    inlet=he.pt(pressure_in,temperature_in)
    def evaluate(x,p,h):
        state=he.ph(float(p),float(h));tw=float(np.interp(x,knots,temps))
        re=4*mass/(math.pi*di*state.viscosity)
        if re>2300:
            raise ValueError('Precooler Reynolds number exceeds 2300. Its current laminar heat-transfer closure is outside its range.')
        # Series radial wall/film resistances per unit length. Nu=3.66 is an
        # explicit fully-developed homogeneous estimate, not a condensation fit.
        g=1/(1/(3.66*math.pi*state.conductivity*factor)+math.log(outer/di)/(2*math.pi*float(steel_k((state.temperature+tw)/2))))
        heat=g*(tw-state.temperature)
        dp=-128*state.viscosity*mass*state.volume/(math.pi*di**4)
        return state,tw,heat,dp,re
    # Conservative implicit finite volumes: the outlet enthalpy is unknown,
    # including inside the saturation interval. A monotone backward step avoids
    # crossing the prescribed wall temperature on a coarse thermal grid.
    xs=[0.];ys=[[pressure_in,inlet.enthalpy,0.]]
    p,h,qtotal=pressure_in,inlet.enthalpy,0.
    for a,b in zip(knots[:-1],knots[1:]):
        count=max(1,int(math.ceil((b-a)/(le*max_step_fraction))))
        grid=np.linspace(a,b,count+1)
        for xa,xb in zip(grid[:-1],grid[1:]):
            dx=float(xb-xa);xm=float((xa+xb)/2)
            before=he.ph(p,h)
            tw=float(np.interp(xm,knots,temps))
            try:
                wall_enthalpy=he.pt(p,tw).enthalpy
            except PhaseDomainError:
                if not he.minimum_saturation_pressure<=p<he.critical_pressure:raise
                saturated_liquid,saturated_vapor=he.pq(p,0),he.pq(p,1)
                if abs(saturated_liquid.temperature-tw)>1e-5:raise
                # A wall exactly at saturation does not prescribe a quality.
                # Bracket toward the nearest equilibrium enthalpy, preserving
                # an existing mixture when no temperature difference exists.
                wall_enthalpy=min(max(h,saturated_liquid.enthalpy),saturated_vapor.enthalpy)
            def balance(hout):
                state,twall,heat,dp,re=evaluate(xm,p,hout)
                return mass*(hout-h)-heat*dx
            lower,upper=sorted([h,wall_enthalpy])
            fl,fu=balance(lower),balance(upper)
            if fl*fu>0:
                if min(abs(fl),abs(fu))>1e-7:raise ValueError(f'Precooler energy step is not bracketed at x={xm:g} m: {fl:g}, {fu:g} W')
                hout=lower if abs(fl)<abs(fu) else upper
            else:
                hout=brentq(balance,lower,upper,xtol=1e-6) if upper-lower>1e-7 else h
            state,twall,heat,dp,re=evaluate(xm,p,hout)
            pnext=p+dp*dx
            if pnext<=0:raise ValueError('Precooler pressure reaches zero; requested flow is too large.')
            qtotal+=heat*dx
            p,h=float(pnext),float(hout)
            xs.append(float(xb));ys.append([p,h,qtotal])
    value=ys[-1]
    profile=[];maxre=0.;maxmach=0.
    for x,y in zip(xs,ys):
        state,tw,heat,dp,re=evaluate(float(x),float(y[0]),float(y[1]))
        point=state.point(x,mass,area);point.update(wall_temperature=tw,outer_wall_temperature=tw)
        profile.append(point);maxre=max(maxre,re)
    outlet=he.ph(float(value[0]),float(value[1]))
    # Enthalpy is integrated against an independently accumulated wall heat.
    error=mass*(outlet.enthalpy-inlet.enthalpy)-float(value[2])
    energy_relative=abs(error)/max(abs(float(value[2])),1e-12)
    if energy_relative>1e-5:raise ValueError('Precooler energy balance failed.')
    return dict(mass_flow=mass,mass_flow_kg_s=mass,pressure_in=pressure_in,pressure_out=outlet.pressure,
                temperature_in=temperature_in,temperature_out=outlet.temperature,
                enthalpy_in_J_kg=inlet.enthalpy,enthalpy_out_J_kg=outlet.enthalpy,
                outlet_state=outlet,profile=profile,
                diagnostics=dict(max_Re=maxre,relative_energy_residual=energy_relative,
                                 heat_transfer_multiplier=factor),
                warnings=['Homogeneous equilibrium precooler: Nu=3.66 and mixture conductivity estimate finite heat exchange. The heat-transfer multiplier needs sensitivity analysis/calibration; slip, stratification and film condensation are not resolved.'],backend='one_dimensional_real_helium')
