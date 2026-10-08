"""Coupled experimental HEM continuous one-k-pot solver.

Initial topology scope: one saturated pot, suction pipe, pump, prescribed-wall
precooler, capillary, back to pot. This solves a common mass flow and discharge
pressure; it does not postprocess a gas-only hydraulic solution into liquid.
Other EDA topologies retain the existing solvers and receive an explicit scope
error if this experimental mode is selected.
"""
from copy import deepcopy
from functools import lru_cache
import math
import numpy as np
from scipy.optimize import least_squares
from .network import validate_network, NetworkError
from .pumps import evaluate_forward_backflow, check_operating_limits
from .two_phase_properties import PhaseDomainError
from .two_phase_pipes import suction_pipe, precooling_pipe
from .two_phase_capillary import path_at_flow, impedance
from .pot_energy import pot_balance, flash_at_pot, reference_shift
import CoolProp
from models.continuous_he4 import properties as he4
from models.comsol_pipe.helium_properties import SPECIFIC_GAS_CONSTANT_J_KG_K as R


MODEL='HEM_continuous_pot_public_properties_v1'


def _topology(project):
    graph=validate_network(project)
    pots=[c for c in graph.components.values() if c['type']=='source' and c['params'].get('mode')=='saturation']
    counts={kind:sum(e['type']==kind for e in graph.edges) for kind in ['pipe','pump','jt']}
    if len(pots)!=1 or counts!={'pipe':2,'pump':1,'jt':1} or len(graph.components)!=5:
        raise NetworkError('Two-phase continuous-pot mode currently supports one saturated pot, two pipes, one pump and one capillary in a single loop.',list(graph.components),'two_phase_topology')
    pot=pots[0];pn=graph.port_nodes[pot['id']+':a']
    order=[];node=pn
    for _ in range(4):
        candidates=[e for e in graph.edges if e['node_a']==node and e['id'] not in [r['id'] for r in order]]
        if len(candidates)!=1:break
        edge=candidates[0];order.append(edge);node=edge['node_b']
    if len(order)!=4 or node!=pn or [e['type'] for e in order]!=['pipe','pump','pipe','jt']:
        raise NetworkError('Connect the two-phase loop A to B: pot -> suction pipe -> pump -> precooler -> capillary -> pot. Component mirroring does not change port identities.',list(graph.components),'two_phase_topology')
    if order[-1]['params'].get('law')!='capillary':
        raise NetworkError('Two-phase continuous-pot mode currently requires the capillary J-T law.',[order[-1]['id']],'two_phase_scope')
    return graph,pot,order


def solve_two_phase_loop(project,*,report=None,cancel_event=None,pipe_cells=480,capillary_cells=640):
    graph,pot,order=_topology(project)
    suction,pump,precool,jt=order
    tp=float(pot['params']['temperature_K']);pp=float(he4.saturation_pressure(tp))
    tout=float(pump['params']['outlet_temperature_K'])
    for edge in (suction,precool):
        multiplier=edge['params'].get('heat_transfer_multiplier',1.)
        if isinstance(multiplier,bool) or not isinstance(multiplier,(int,float)) or not math.isfinite(multiplier) or not 0<multiplier<=1000:
            raise NetworkError('Heat-transfer multiplier must be numeric and within (0, 1000].',[edge['id']],'invalid_parameter')
    rough=jt['params'].get('roughness_m',0.)
    if isinstance(rough,bool) or not isinstance(rough,(int,float)) or not math.isfinite(rough) or not 0<=rough<=.1*jt['params']['diameter_m']:
        raise NetworkError('Capillary roughness must be numeric, nonnegative and at most 10% of its bore.',[jt['id']],'invalid_parameter')
    history=[]
    def progress(message):
        if cancel_event is not None and cancel_event.is_set():
            from .comsol_backend import ComsolCancelled
            raise ComsolCancelled('Calculation cancelled')
        if report:report(min(.9,.1+.015*len(history)),message)
    def component_call(edge,fn,*args,**kwargs):
        try:return fn(*args,**kwargs)
        except NetworkError:raise
        except (ValueError,ArithmeticError) as error:
            label=graph.components[edge['id']].get('label',edge['id'])
            raise NetworkError(f'{label}: {error}',[edge['id']],
                               'two_phase_property_range' if isinstance(error,PhaseDomainError) else 'two_phase_component',
                               stage='Two-phase component evaluation') from error
    suction_guesses={}
    @lru_cache(maxsize=128)
    def suction_at(mass):
        closest=min(suction_guesses,key=lambda m:abs(math.log(m/mass))) if suction_guesses else None
        value=component_call(suction,suction_pipe,suction['params'],mass,pp,tp,initial_guess=suction_guesses.get(closest))
        suction_guesses[mass]=value.pop('_bvp_guess')
        return value
    def evaluate(values):
        mass,pout=np.exp(values)
        progress('Solving real-helium circulation, precooling and capillary flashing')
        s=suction_at(float(mass));pin=s['pressure_out'];tin=s['temperature_out']
        par=pump['params'];mode=par.get('mode','constant_speed')
        if mode=='forward_backflow':
            pu=evaluate_forward_backflow(par,pin,pout,tin,tout)
        elif mode=='constant_speed':
            pu=dict(mass_flow=par['speed_m3_s']*pin/(R*tin),pump_mode=mode,warnings=[])
        else:
            formula=graph.components[pump['id']]['compiled_expression']
            pu=dict(mass_flow=formula(Pin=pin,Pout=pout,Tin=tin,S=par['speed_m3_s'],R=R),pump_mode=mode,warnings=[])
        pre=component_call(precool,precooling_pipe,precool['params'],mass,pout,tout,max_step_fraction=1/pipe_cells)
        cap=component_call(jt,path_at_flow,jt['params'],pre['pressure_out'],pre['enthalpy_out_J_kg'],pp,mass,cells=capillary_cells)
        if cap['domain_limited']:
            raise NetworkError('The capillary enters the He II domain before choking. Public equilibrium data alone cannot close superfluid hydraulic resistance.',[jt['id']],'two_phase_property_range')
        residual=[(pu['mass_flow']-mass)/max(mass,1e-16),math.log(max(cap['impedance'],1.)/impedance(jt['params'])[2])]
        history.append(dict(mass_flow_kg_s=float(mass),pump_outlet_pressure_Pa=float(pout),residual=[float(r) for r in residual]))
        return residual,(s,pu,pre,cap)
    # Bound trial circulation using actual suction closure. This is a numerical
    # bracket, not a physical flow clamp on the accepted solution.
    par=pump['params'];speed=float(par.get('forward_speed_m3_s',par.get('speed_m3_s',.24)))
    nominal=max(speed*pp/(R*tout),1e-12)
    upper=nominal*1.8
    for _ in range(50):
        try:suction_at(upper);break
        except NetworkError as error:
            if error.code=='two_phase_property_range':raise
            upper*=.8
    else:raise NetworkError('Cannot find a positive admissible suction flow.',[suction['id']],'two_phase_initialization')
    initial=np.log([min(nominal*.7,upper*.8),max(pp*2,1e5)])
    solved=least_squares(lambda y:evaluate(y)[0],initial,
                         bounds=(np.log([1e-12,max(pp*1.001,100.)]),np.log([upper,1e7])),
                         diff_step=2e-4,xtol=2e-8,ftol=2e-8,gtol=2e-8,max_nfev=90)
    residual,(s,pu,pre,cap)=evaluate(solved.x)
    if max(abs(v) for v in residual)>2e-3:
        raise NetworkError(f'Two-phase circulation did not converge; maximum constitutive residual is {max(abs(v) for v in residual):.3g}.',[e['id'] for e in order],'two_phase_nonconvergence')
    mass,pout=map(float,np.exp(solved.x))
    if s['diagnostics']['max_Re']>2300:
        raise NetworkError('Suction-pipe Reynolds number exceeds 2300. Its laminar heat-transfer closure needs a turbulent extension.',[suction['id']],'two_phase_scope')
    if pu['pump_mode']=='forward_backflow':
        component_call(pump,check_operating_limits,pump['params'],s['pressure_out'],pout,pu)
    feed=pre.pop('outlet_state')
    flash=flash_at_pot(feed.enthalpy,tp,pp)
    # The capillary exit and the fully mixed downstream pot are separate states.
    # Kinetic energy of a choked jet is recovered in the pot control volume.
    incoming=[dict(component_id=jt['id'],mass_flow_kg_s=mass,enthalpy_J_kg=feed.enthalpy,
                   upstream_phase=feed.phase,method='Coupled real-helium precooler enthalpy',state_at_pot_pressure=flash)]
    balance=pot_balance(tp,incoming,mass,float(pot['params'].get('heat_load_W',0.)))
    balance.update(model=MODEL,notes=['Homogeneous equilibrium engineering model; heat-transfer and friction closures require calibration. He II is used for pot equilibrium only.'],property_sources={'pot':'Donnelly & Barenghi (1998)','feed':'CoolProp HEOS, DB98 reference alignment','CoolProp_version':CoolProp.__version__,'reference_shift_J_kg':reference_shift()})
    profile=cap.pop('profile')
    # Normalize physical positions only within the accepted constitutive residual.
    if profile:
        scale=jt['params']['length_m']/profile[-1]['x']
        for point in profile:point['x']*=scale
    jrow=dict(mass_flow=mass,mass_flow_kg_s=mass,pressure_in=feed.pressure,pressure_out=pp,
              temperature_in=feed.temperature,temperature_out=flash['temperature_K'],
              enthalpy_in_J_kg=feed.enthalpy,enthalpy_out_J_kg=feed.enthalpy,
              outlet_fluid_temperature_K=flash['temperature_K'],flash_at_pot=flash,
              profile=profile,diagnostics=cap,thermodynamic_model='HEM_adiabatic_total_enthalpy',
              warnings=['HEM assumes equal phase velocity and equilibrium. Churchill friction uses smooth walls unless roughness is supplied; measured z sets effective hydraulic length and is not a two-phase calibration.',
                        'The axial profile ends at the capillary exit. If choked, this pressure exceeds pot pressure; the pot flash includes downstream jet kinetic-energy recovery.'])
    pu.update(mass_flow=mass,mass_flow_kg_s=mass,pressure_in=s['pressure_out'],pressure_out=pout,
              temperature_in=s['temperature_out'],temperature_out=tout,
              profile=[dict(x=0.,pressure=s['pressure_out'],temperature=s['temperature_out']),dict(x=1.,pressure=pout,temperature=tout)])
    rows=[s,pu,pre,jrow];edges=[];components={};warnings=[]
    for edge,row in zip(order,rows):
        row.update(id=edge['id'],type=edge['type'],label=graph.components[edge['id']].get('label',edge['id']),node_a=edge['node_a'],node_b=edge['node_b'])
        components[edge['id']]=row;edges.append(row)
        warnings.extend(row['label']+': '+w for w in row.get('warnings',[]))
    nodes=[]
    node_values=[(order[0]['node_a'],pp,tp),(order[0]['node_b'],s['pressure_out'],s['temperature_out']),
                 (order[1]['node_b'],pout,tout),(order[2]['node_b'],feed.pressure,feed.temperature)]
    for nid,p,t in node_values:
        original=graph.nodes[nid]
        nodes.append(dict(original,pressure=float(p),temperature=float(t),net_mass_supply_kg_s=0.,mass_balance_kg_s=0.,external_heat_W=balance['net_cooling_power_W'] if nid==order[0]['node_a'] else 0.,cooling_load_W=-balance['net_cooling_power_W'] if nid==order[0]['node_a'] else 0.,thermal_model='advected_enthalpy_no_intercomponent_conduction'))
    components[pot['id']]=dict(id=pot['id'],type='source',label=pot.get('label',pot['id']),node=order[0]['node_a'],pressure=pp,temperature=tp,
                               mass_flow=0.,net_mass_supply_kg_s=0.,pot_energy=balance)
    result=dict(status='complete',engine='two_phase',model=MODEL,components=components,edges=edges,nodes=nodes,
                warnings=warnings,assumptions=[
                    'One saturated pot -> adiabatic suction pipe -> pump -> prescribed-wall precooler -> capillary -> pot; forward circulation only.',
                    'No intercomponent solid conduction; incoming/outgoing helium enthalpy is retained. The same mass flow passes through every component.',
                    'HEM: local vapor-liquid thermal equilibrium and equal phase velocity. No stratification, phase slip, metastable flashing, superfluid film or He II channel-flow law.',
                    'Precooler Nu=3.66 with homogeneous mixture transport and a configurable heat-transfer multiplier. Finite gas-wall exchange determines outlet temperature and vapor quality.',
                    'Axial gas conduction and kinetic energy in the large pipes are neglected. Axial steel conduction is retained in the adiabatic suction pipe; the precooler outer wall is prescribed.',
                    'Capillary total enthalpy includes kinetic energy. Friction uses the Churchill all-regime law; the liquid/mixture laws are engineering approximations requiring validation.'
                ],diagnostics=dict(pot_energy_model=MODEL,max_relative_mass_residual=0.,max_constitutive_residual=max(abs(v) for v in residual),
                                    history=history,pipe_cells=pipe_cells,capillary_cells=capillary_cells),
                summary=dict(engine='two_phase',mass_balance_relative=0.,pressure_min_Pa=min(n['pressure'] for n in nodes),
                             pressure_max_Pa=max(n['pressure'] for n in nodes),max_Re=max(r.get('diagnostics',{}).get('max_Re',0.) for r in rows),max_Kn=s['diagnostics']['max_Kn']))
    progress('Two-phase circulation and pot enthalpy balance complete')
    return result
