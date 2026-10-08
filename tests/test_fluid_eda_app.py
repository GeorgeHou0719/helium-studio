import json
from pathlib import Path
import tempfile
import threading
import time
import unittest
from copy import deepcopy
from models.fluid_eda.app import create_app
from models.fluid_eda.examples import continuous_pump_example
from models.fluid_eda.simulation import run_simulation
from models.fluid_eda.comsol_backend import ComsolCancelled, ComsolError
from models.fluid_eda.network import NetworkError
from models.fluid_eda.error_reporting import describe_error
from models.fluid_eda.thermal import reduced_pipe


class CouplingTests(unittest.TestCase):
    def test_low_pressure_reaches_backend_without_clamping(self):
        project=continuous_pump_example()
        # Preserve the narrow suction-line regression independently of the example defaults.
        project['components'][1]['params'].update(length_m=.288,outer_diameter_m=.00635,wall_thickness_m=.000762)
        project['components'][2]['params']['mode']='forward_backflow'
        captured=[]
        class RecordingBackend:
            def __init__(self,**kw):pass
            def __enter__(self):return self
            def __exit__(self,*a):pass
            def evaluate_pipe(self,params,pa,pb,ta,tb,**kw):
                captured.append((params['_component_id'],pa,pb))
                return reduced_pipe(params,pa,pb,ta,tb,profiles=True)
        with tempfile.TemporaryDirectory() as path:
            result=run_simulation(project,output_dir=path,backend_factory=RecordingBackend)
        suction=next(row for row in captured if row[0]=='pipe1')
        self.assertGreater(suction[2],0)
        self.assertLess(suction[2],.003)
        self.assertEqual(result['components']['pipe1']['pressure_out'],suction[2])
        self.assertTrue(result['comsol_coupling']['converged'])

    def test_backend_failure_and_cancellation_keep_their_identity(self):
        failure=ComsolError('Nonlinear solver failed on pipe 2\nDetailed solver trace')
        class FailingBackend:
            def __init__(self,**kw):pass
            def __enter__(self):return self
            def __exit__(self,*a):pass
            def evaluate_pipe(self,params,pa,pb,ta,tb,**kw):
                if params['_component_id']=='pipe2':raise failure
                return reduced_pipe(params,pa,pb,ta,tb,profiles=True)
        project=continuous_pump_example()
        with tempfile.TemporaryDirectory() as path:
            with self.assertRaises(NetworkError) as caught:
                run_simulation(project,output_dir=path,backend_factory=FailingBackend)
            detail=describe_error(caught.exception,project)
            self.assertEqual(detail['component_ids'],['pipe2'])
            self.assertEqual(detail['components'][0]['label'],'Staged precooling pipe')
            self.assertEqual(detail['stage'],'COMSOL pipe solve')
            self.assertIn('Detailed solver trace',detail['details'])
            failure=ComsolCancelled('User cancelled')
            with self.assertRaises(ComsolCancelled):
                run_simulation(project,output_dir=path,backend_factory=FailingBackend)

    def test_horizontal_flip_does_not_change_network_physics(self):
        project=continuous_pump_example();project['settings']['engine']='reduced'
        original=deepcopy(project)
        before=run_simulation(project)
        for c in project['components']:c['flipX']=True
        after=run_simulation(project)
        self.assertEqual(before['components'],after['components'])
        self.assertEqual(project['wires'],original['wires'])

    def test_comsol_correction_reaches_true_network_balance(self):
        class StandIn:
            calls=0
            def __init__(self,**kw):pass
            def __enter__(self):return self
            def __exit__(self,*a):pass
            def evaluate_pipe(self,params,pa,pb,ta,tb,**kw):
                StandIn.calls+=1
                r=reduced_pipe(params,pa,pb,ta,tb,profiles=True)
                r['mass_flow']*=1.08;r['mass_flow_kg_s']*=1.08;r['backend']='test_double'
                return r
        with tempfile.TemporaryDirectory() as path:
            r=run_simulation(continuous_pump_example(),output_dir=path,backend_factory=StandIn)
        self.assertGreaterEqual(StandIn.calls,4)
        self.assertLess(r['summary']['mass_balance_relative'],.005)
        self.assertTrue(r['comsol_coupling']['converged'])
        self.assertEqual(r['components']['pipe1']['backend'],'test_double')

    def test_cancel_before_comsol_start(self):
        cancel=threading.Event();cancel.set()
        with self.assertRaises(ComsolCancelled):run_simulation(continuous_pump_example(),cancel_event=cancel)

    def test_forward_backflow_survives_comsol_correction(self):
        class StandIn:
            def __init__(self,**kw):pass
            def __enter__(self):return self
            def __exit__(self,*a):pass
            def evaluate_pipe(self,params,pa,pb,ta,tb,**kw):
                row=reduced_pipe(params,pa,pb,ta,tb,profiles=True)
                row['mass_flow']*=1.08;row['mass_flow_kg_s']*=1.08
                return row
        project={'settings':{'engine':'comsol'},'components':[
            {'id':'a','type':'source','params':{'pressure_Pa':10,'temperature_K':293.15}},
            {'id':'b','type':'source','params':{'pressure_Pa':100,'temperature_K':293.15}},
            {'id':'pump','type':'pump','params':{'mode':'forward_backflow','outlet_temperature_K':293.15}},
            {'id':'pipe','type':'pipe','params':{'length_m':.288,'outer_diameter_m':.00635,
             'wall_thickness_m':.000762,'T_a_K':293.15,'T_b_K':293.15}}],
            'wires':[{'id':'w1','from':{'component':'a','port':'a'},'to':{'component':'pump','port':'a'}},
                     {'id':'w2','from':{'component':'pump','port':'b'},'to':{'component':'pipe','port':'a'}},
                     {'id':'w3','from':{'component':'pipe','port':'b'},'to':{'component':'b','port':'a'}}]}
        with tempfile.TemporaryDirectory() as path:
            result=run_simulation(project,output_dir=path,backend_factory=StandIn)
        self.assertTrue(result['comsol_coupling']['converged'])
        self.assertGreaterEqual(len(result['comsol_coupling']['iterations']),2)
        self.assertEqual(result['components']['pump']['pump_mode'],'forward_backflow')
        self.assertGreater(result['components']['pump']['diagnostics']['backflow_mass_flow_kg_s'],0)
        self.assertLess(result['summary']['mass_balance_relative'],.005)


class APITests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.app=create_app(results_dir=self.temp.name)
        self.client=self.app.test_client()
    def tearDown(self):
        self.app.extensions['eda']['executor'].shutdown(wait=True,cancel_futures=True)
        self.temp.cleanup()
    def test_example_to_completed_preview_and_export(self):
        project=self.client.get('/api/example').get_json();project['settings']['engine']='reduced'
        response=self.client.post('/api/simulate',json=project)
        self.assertEqual(response.status_code,202)
        job=response.get_json()['job_id']
        deadline=time.monotonic()+10
        while time.monotonic()<deadline:
            data=self.client.get('/api/jobs/'+job).get_json()
            if data['status'] in ('complete','error'):break
            time.sleep(.02)
        self.assertEqual(data['status'],'complete',data)
        self.assertLess(data['result']['summary']['mass_balance_relative'],1e-5)
        file=self.client.get(f'/api/jobs/{job}/files/results.json')
        self.assertEqual(file.status_code,200)
        self.assertEqual(json.loads(file.data)['engine'],'reduced')
        file.close()
        self.assertEqual(self.client.get(f'/api/jobs/{job}/files/../../app.py').status_code,404)
    def test_invalid_topology_and_unsafe_expression(self):
        self.assertEqual(self.client.post('/api/simulate',json={}).status_code,400)
        p=continuous_pump_example();p['components'][2]['params'].update(mode='expression',expression='__import__("os").getcwd()')
        response=self.client.post('/api/simulate',json=p)
        self.assertEqual(response.status_code,400)
        self.assertIn('pump1',response.get_json()['error']['component_ids'])
        self.assertEqual(response.get_json()['error']['components'][0]['label'],'Circulation pump')
        self.assertEqual(response.get_json()['error']['stage'],'Circuit validation')

    def test_async_invalid_mesh_error_retains_component_context(self):
        project=continuous_pump_example();project['components'][1]['params']['axial_elements']=39
        response=self.client.post('/api/simulate',json=project)
        self.assertEqual(response.status_code,202)
        job=response.get_json()['job_id'];deadline=time.monotonic()+10
        while time.monotonic()<deadline:
            data=self.client.get('/api/jobs/'+job).get_json()
            if data['status'] in ('complete','error'):break
            time.sleep(.02)
        self.assertEqual(data['status'],'error')
        self.assertEqual(data['error']['component_ids'],['pipe1'])
        self.assertIn('axial_elements',data['error']['message'])
        self.assertEqual(data['error']['components'][0]['label'],'Suction pipe')
        self.assertEqual(data['error']['stage'],'COMSOL input validation')
        self.app.extensions['eda']['jobs'][job]['future'].result(timeout=5)
        recorded=json.loads((Path(self.temp.name)/job/'job.json').read_text(encoding='utf-8'))
        self.assertEqual(recorded['error'],data['error'])

    def test_error_metadata_handles_multiple_devices_and_unattributed_errors(self):
        project=continuous_pump_example()
        detail=describe_error(NetworkError('Conflicting boundaries',['pipe1','pump1']),project)
        self.assertEqual([c['label'] for c in detail['components']],['Suction pipe','Circulation pump'])
        self.assertEqual(describe_error(RuntimeError('Server failure'),project)['components'],[])
        self.assertEqual(describe_error(ValueError('Bad JSON'),None)['components'],[])
    def test_cross_site_and_host_protection(self):
        self.assertEqual(self.client.get('/api/status',headers={'Host':'evil.test'}).status_code,403)
        self.assertEqual(self.client.post('/api/simulate',json={},headers={'Origin':'https://evil.test'}).status_code,403)
        self.assertEqual(self.client.get('/api/jobs/missing').status_code,404)
    def test_job_cancel_is_cooperative(self):
        def slow(project,*,cancel_event,**kw):
            cancel_event.wait(5)
            if cancel_event.is_set():raise ComsolCancelled('cancelled')
            raise RuntimeError('cancel request was not received')
        other=create_app(results_dir=Path(self.temp.name)/'cancel',runner=slow);client=other.test_client()
        try:
            job=client.post('/api/simulate',json=continuous_pump_example()).get_json()['job_id']
            self.assertEqual(client.post('/api/jobs/'+job+'/cancel').status_code,200)
            for _ in range(100):
                data=client.get('/api/jobs/'+job).get_json()
                if data['status']=='cancelled':break
                time.sleep(.01)
            self.assertEqual(data['status'],'cancelled')
        finally:other.extensions['eda']['executor'].shutdown(wait=True,cancel_futures=True)

    def test_forward_backflow_pump_api_preserves_signed_result_and_error_identity(self):
        project={'settings':{'engine':'reduced'},'components':[
            {'id':'a','type':'source','params':{'pressure_Pa':.1,'temperature_K':293.15}},
            {'id':'b','type':'source','params':{'pressure_Pa':200,'temperature_K':293.15}},
            {'id':'pump','type':'pump','label':'HiPace test pump','params':
             {'mode':'forward_backflow','compression_ratio':1000,'outlet_temperature_K':293.15}}],
            'wires':[{'id':'w1','from':{'component':'a','port':'a'},'to':{'component':'pump','port':'a'}},
                     {'id':'w2','from':{'component':'pump','port':'b'},'to':{'component':'b','port':'a'}}]}
        def run():
            response=self.client.post('/api/simulate',json=project)
            self.assertEqual(response.status_code,202)
            jid=response.get_json()['job_id']
            self.app.extensions['eda']['jobs'][jid]['future'].result(timeout=10)
            return self.client.get('/api/jobs/'+jid).get_json()
        data=run()
        self.assertEqual(data['status'],'complete',data)
        pump=data['result']['components']['pump']
        self.assertLess(pump['mass_flow'],0)
        self.assertEqual(pump['diagnostics']['flow_direction'],'backflow')
        self.assertGreater(pump['diagnostics']['backflow_mass_flow_kg_s'],pump['diagnostics']['forward_mass_flow_kg_s'])
        project['components'][2]['params']['max_outlet_pressure_Pa']=100
        data=run()
        self.assertEqual(data['status'],'error',data)
        self.assertEqual(data['error']['code'],'pump_operating_limit')
        self.assertEqual(data['error']['context']['port'],'b')
        self.assertEqual(data['error']['context']['value'],200)
        self.assertEqual(data['error']['components'][0]['label'],'HiPace test pump')
        self.assertIn('HiPace test pump',data['error']['message'])
        project['components'][2]['params']['compression_curve']=[[1,100],[1,10]]
        response=self.client.post('/api/simulate',json=project)
        self.assertEqual(response.status_code,400)
        self.assertEqual(response.get_json()['error']['component_ids'],['pump'])

if __name__=='__main__':unittest.main()
