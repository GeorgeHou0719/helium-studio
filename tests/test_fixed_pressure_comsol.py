"""Fixed-pressure COMSOL circuits bypass the preview thermal solver."""
from copy import deepcopy
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from models.fluid_eda.simulation import run_simulation
from models.fluid_eda.network import NetworkError, CP
from models.fluid_eda.comsol_backend import ComsolError, ComsolCancelled


def circuit():
    return {'settings':{'engine':'comsol'}, 'components':[
        {'id':'cold','type':'source','params':{'mode':'saturation','temperature_K':1}},
        {'id':'pipe','label':'PIPE-01','type':'pipe','params':{'length_m':.5,
         'outer_diameter_m':.0635,'wall_thickness_m':.00762,'T_a_K':1,'T_b_K':4}},
        {'id':'warm','type':'source','params':{'mode':'constant','pressure_Pa':.5,'temperature_K':4}}],
        'wires':[{'id':'a','from':{'component':'cold','port':'b'},'to':{'component':'pipe','port':'a'}},
                 {'id':'b','from':{'component':'pipe','port':'b'},'to':{'component':'warm','port':'a'}}]}


class FixedPressureTests(unittest.TestCase):
    def test_direct_comsol_preserves_states_energy_and_signed_flow_without_preview(self):
        for reverse in (False,True):
            project=circuit()
            if reverse:project['components'][2]['params']['pressure_Pa']=20
            original=deepcopy(project);calls=[]
            class Backend:
                def __init__(self,**kwargs):pass
                def __enter__(self):return self
                def __exit__(self,*args):pass
                def evaluate_pipe(self,params,pa,pb,ta,tb,**kwargs):
                    calls.append((pa,pb,ta,tb,kwargs['save_path']))
                    mass=1e-6 if pa>pb else -1e-6
                    fa,fb=mass*CP*ta+.01,mass*CP*tb+.02
                    return {'backend':'recording COMSOL test double','mass_flow':mass,
                        'mass_flow_kg_s':mass,'heat_a_W':.3,'heat_b_W':-.3-(fb-fa),
                        'diagnostics':{'flow_energy_a':fa,'flow_energy_b':fb},
                        'profile':[{'x':0,'pressure':pa,'temperature':ta},
                                   {'x':.5,'pressure':pb,'temperature':tb}], 'warnings':[]}
            with tempfile.TemporaryDirectory() as path, patch('models.fluid_eda.simulation.hydraulic_pipe_seed',side_effect=AssertionError('Preview must not run')), patch('models.fluid_eda.thermal.reduced_pipe',side_effect=AssertionError('Preview must not run')):
                result=run_simulation(project,backend_factory=Backend,output_dir=path)
            self.assertEqual(len(calls),1)
            pa,pb,ta,tb,artifact=calls[0]
            self.assertAlmostEqual(pa,15.56999892255892)
            self.assertEqual((pb,ta,tb),(20 if reverse else .5,1,4))
            self.assertEqual(Path(artifact).name,'pipe_1.mph')
            self.assertEqual(result['comsol_coupling']['mode'],'fixed_pressures')
            self.assertEqual(result['comsol_coupling']['iterations'],[])
            self.assertEqual(result['diagnostics']['function_evaluations'],0)
            self.assertAlmostEqual(result['components']['cold']['cooling_load_W'],-.01)
            self.assertAlmostEqual(result['diagnostics']['net_source_mass_flow_kg_s'],0)
            self.assertEqual(result['components']['pipe']['mass_flow'], -1e-6 if reverse else 1e-6)
            self.assertEqual(project,original)

    def test_errors_and_cancellation_retain_pipe_context(self):
        calls=[]
        class Backend:
            error=ComsolError('actual COMSOL solver failure')
            def __init__(self,**kwargs):pass
            def __enter__(self):return self
            def __exit__(self,*args):pass
            def evaluate_pipe(self,*args,**kwargs):calls.append(args);raise self.error
        with tempfile.TemporaryDirectory() as path:
            with self.assertRaises(NetworkError) as caught:
                run_simulation(circuit(),backend_factory=Backend,output_dir=path)
            self.assertEqual(caught.exception.code,'comsol_pipe_failure')
            self.assertEqual(caught.exception.component_ids,['pipe'])
            self.assertEqual(caught.exception.context['temperature_b'],4)
            Backend.error=ComsolCancelled('cancelled')
            with self.assertRaises(ComsolCancelled):
                run_simulation(circuit(),backend_factory=Backend,output_dir=path)
            calls.clear();project=circuit();project['components'][1]['params']['axial_elements']=39
            with self.assertRaises(NetworkError) as caught:
                run_simulation(project,backend_factory=Backend,output_dir=path)
            self.assertEqual(caught.exception.stage,'COMSOL input validation')
            self.assertEqual(calls,[])

    def test_quick_mode_is_not_silently_promoted_to_comsol(self):
        project=circuit();project['settings']['engine']='reduced'
        with self.assertRaisesRegex(NetworkError,'Advection is too strong'):
            run_simulation(project,backend_factory=lambda **kw: self.fail('COMSOL must not start'))


if __name__=='__main__':unittest.main()
