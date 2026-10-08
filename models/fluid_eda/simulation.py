"""Job-level network/COMSOL coupling, with explicit convergence auditing."""
from __future__ import annotations
from copy import deepcopy
import math
from pathlib import Path
import threading
from .network import solve_network, validate_network, refresh_result, NetworkError
from .hydraulic_seed import hydraulic_pipe_seed
from .comsol_backend import ComsolPipeBackend, ComsolCancelled, normalize_parameters
from .error_reporting import pipe_error
from .pot_energy import apply_pot_balances
from models.config import DATA_DIR


def _summary(result,engine,project):
    apply_pot_balances(result,project)
    result["diagnostics"]["pot_energy_model"] = "DB98_pot_CoolProp_feed"
    obsolete = "The current ideal-helium model does not calculate real J"
    result["assumptions"] = [a for a in result.get("assumptions",[]) if not a.startswith(obsolete)]
    result["assumptions"].append("No intercomponent conduction. Pot cooling retains incoming feed enthalpy and outgoing saturated-vapor enthalpy; pipe/restriction hydraulics remain gas approximations.")
    edges=result['edges'];diags=[e.get('diagnostics',{}) for e in edges]
    result['engine']=engine;result['status']='complete'
    result['summary']={
        'engine':engine,
        'mass_balance_relative':result['diagnostics']['max_relative_mass_residual'],
        'max_Kn':max([d.get('max_Kn',d.get('max_knudsen',0)) for d in diags] or [0]),
        'max_Re':max([d.get('max_Re',d.get('max_reynolds',0)) for d in diags] or [0]),
        'total_external_heat_W':sum(n['external_heat_W'] for n in result['nodes']),
        'pressure_min_Pa':min(n['pressure'] for n in result['nodes']),
        'pressure_max_Pa':max(n['pressure'] for n in result['nodes']),
    }
    if engine=='reduced':
        result['warnings'].append('Quick preview assumes equal gas and wall temperatures in each cross section and neglects axial inertia. These pipe results were not computed by COMSOL.')
    if any(e['type']=='pump' and e.get('pump_mode','constant_speed')=='constant_speed' for e in edges):
        result['warnings'].append('The constant-speed pump is an ideal prescribed-flow component. Ultimate pressure, maximum compression ratio, power limits, compression work and pump cooling power are not modeled. Calibrate against measured pump curves.')
    if any(e.get('pump_mode')=='expression' for e in edges):
        result['warnings'].append('Custom pump flow expressions require calibration. Pump work and cooling power are not calculated.')
    return result


def run_simulation(project,*,progress=None,cancel_event=None,output_dir=None,backend_factory=ComsolPipeBackend):
    """Use COMSOL only at converged network trials, never inside finite differences.

    A hydraulic-only seed bypasses the quick thermal solver for every topology.
    Pipe conductances are corrected from real 2-D solves, then the entire
    network is re-solved. The accepted result contains the actual final COMSOL
    branches; its KCL error is audited independently (default <=0.5%).
    """
    circuit=deepcopy(project);cancel_event=cancel_event or threading.Event()
    def report(fraction,message):
        if cancel_event.is_set(): raise ComsolCancelled('Calculation cancelled')
        if progress: progress(fraction,message)
    engine=circuit.get('settings',{}).get('engine','comsol')
    if engine=='two_phase':
        from .two_phase_loop import solve_two_phase_loop
        report(.02,'Checking continuous-pot two-phase topology')
        result=solve_two_phase_loop(circuit,report=report,cancel_event=cancel_event)
        report(1,'Two-phase simulation complete')
        return result
    if engine not in ('reduced','comsol'): raise NetworkError('Unknown pipe-solver mode')
    report(.02,'Checking circuit topology and thermal boundaries')
    if engine=='reduced':
        result=solve_network(circuit,progress=lambda m:report(.2,m))
        result['assumptions']=project.get('assumptions',[])
        report(1,'Quick preview complete')
        return _summary(result,engine,circuit)

    # Stable component identity is carried only in a private, copied parameter.
    for component in circuit['components']:
        component.setdefault('params',{})['_component_id']=component['id']
    pipes=[c for c in circuit['components'] if c['type']=='pipe']
    factors={c['id']:1. for c in pipes}
    output_dir=Path(output_dir or DATA_DIR/'manual')
    output_dir.mkdir(parents=True,exist_ok=True)
    history=[];tolerance=.005;maximum_rounds=8
    def corrected(params,pa,pb,ta,tb,profiles=False):
        try:
            value=hydraulic_pipe_seed(params,pa,pb,ta,tb,profiles=profiles)
        except (ValueError, TypeError) as error:
            component=graph.components[params['_component_id']]
            raise pipe_error(error,component,pa,pb,ta,tb,validation=True) from error
        factor=factors[params['_component_id']]
        value['mass_flow']*=factor;value['mass_flow_kg_s']*=factor
        return value
    current_progress=.05
    graph=validate_network(circuit)
    with backend_factory(progress=lambda m:report(current_progress,m),cancel_event=cancel_event) as backend:
        def evaluate(component,index,pa,pb,ta,tb):
            params=component['params']
            try:
                normalize_parameters(params,pa,pb,ta,tb)
            except (ValueError, TypeError) as error:
                raise pipe_error(error,component,pa,pb,ta,tb,validation=True) from error
            try:
                # Never construct artifact paths from user-provided component IDs.
                return backend.evaluate_pipe(params,pa,pb,ta,tb,save_path=output_dir/f'pipe_{index+1}.mph')
            except ComsolCancelled:
                raise
            except Exception as error:
                raise pipe_error(error,component,pa,pb,ta,tb) from error

        if all(node['pressure'] is not None for node in graph.nodes.values()):
            # Reservoirs already determine every pipe's boundary state. No
            # preview temperature solution or conductance correction is needed.
            actual={}
            edges={edge['id']:edge for edge in graph.edges}
            for index,component in enumerate(pipes):
                current_progress=.08+.85*index/max(len(pipes),1)
                report(current_progress,f'COMSOL fixed-pressure solve · {component.get("label",component["id"])}')
                edge=edges[component['id']]
                a,b=graph.nodes[edge['node_a']],graph.nodes[edge['node_b']]
                actual[component['id']]=evaluate(component,index,a['pressure'],b['pressure'],a['temperature'],b['temperature'])
            def solved_pipe(params,pa,pb,ta,tb,profiles=False):
                return deepcopy(actual[params['_component_id']])
            result=solve_network(circuit,pipe_evaluator=solved_pipe,progress=lambda m:report(.96,m))
            result=refresh_result(result,raise_on_failure=True)
            result['comsol_coupling']={'converged':True,'relative_tolerance':tolerance,
                                      'iterations':[],'pipe_count':len(pipes),'mode':'fixed_pressures'}
            result['assumptions']=project.get('assumptions',[])
            report(1,'COMSOL fixed-pressure simulation complete')
            return _summary(result,engine,circuit)

        for iteration in range(maximum_rounds):
            result=solve_network(circuit,pipe_evaluator=corrected,progress=lambda m:report(current_progress,m))
            actual={};changes=[]
            for index,component in enumerate(pipes):
                current_progress=.08+.85*(iteration+(index/max(len(pipes),1)))/maximum_rounds
                cid=component['id'];edge=result['components'][cid]
                report(current_progress,f'COMSOL Correction {iteration+1} / {maximum_rounds} · {component.get("label",cid)}')
                params=component['params'];pa=edge['pressure_in'];pb=edge['pressure_out']
                ta=edge['temperature_in'];tb=edge['temperature_out']
                solved=evaluate(component,index,pa,pb,ta,tb)
                actual[cid]={**edge,**solved,'id':cid,'node_a':edge['node_a'],'node_b':edge['node_b'],
                             'type':'pipe','label':edge['label']}
                raw=hydraulic_pipe_seed(params,pa,pb,ta,tb)['mass_flow']
                mass=solved['mass_flow']
                if abs(raw)>1e-20:
                    correction=mass/raw
                    if not math.isfinite(correction) or correction<=0:
                        raise NetworkError(f'{cid}: COMSOL and the hydraulic starting estimate predict opposite flow directions',[cid],'comsol_correction')
                    factors[cid]=correction
                elif abs(mass)>1e-17:
                    raise NetworkError(f'{cid}: Cannot correct flow resistance near zero pressure difference',[cid],'comsol_correction')
                changes.append(abs(mass-edge['mass_flow'])/max(abs(mass),abs(edge['mass_flow']),1e-18))
            for cid,edge in actual.items():result['components'][cid]=edge
            result['edges']=[actual.get(e['id'],e) for e in result['edges']]
            # Recompute conservation and heat duties from true COMSOL branches.
            result=refresh_result(result)
            mismatch=result['diagnostics']['max_relative_mass_residual']
            history.append({'iteration':iteration+1,'max_relative_mass_residual':mismatch,
                            'max_pipe_change':max(changes or [0]),'conductance_factors':dict(factors)})
            if mismatch<=tolerance and max(changes or [0])<=2*tolerance:
                result['comsol_coupling']={'converged':True,'relative_tolerance':tolerance,
                                           'iterations':history,'pipe_count':len(pipes),
                                           'mode':'hydraulic_conductance_correction'}
                result['assumptions']=project.get('assumptions',[])
                report(1,'COMSOL coupled simulation complete; actual branch mass conservation has been verified')
                return _summary(result,engine,circuit)
            circuit.setdefault('settings',{})['initial_pressures_Pa']={n['id']:n['pressure'] for n in result['nodes']}
        raise NetworkError(f'COMSOL coupled simulation did not meet the mass-conservation threshold {tolerance:.1%}, Final residual {mismatch:.3%}.',
                           [c['id'] for c in pipes],'comsol_nonconvergence')
