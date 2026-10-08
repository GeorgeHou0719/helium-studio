"""General COMSOL coupling must never invoke the quick thermal BVP."""
from copy import deepcopy
import math
import tempfile
import unittest
from unittest.mock import patch
from models.fluid_eda.hydraulic_seed import hydraulic_pipe_seed
from models.fluid_eda.simulation import run_simulation
from models.fluid_eda.network import CP
from models.fluid_eda.examples import continuous_pump_example
from models.comsol_pipe.helium_properties import gas_viscosity, SPECIFIC_GAS_CONSTANT_J_KG_K as R


def source(cid, pressure, temperature=20):
    return {'id':cid,'type':'source','params':{'mode':'constant','pressure_Pa':pressure,'temperature_K':temperature}}

def pipe(cid):
    return {'id':cid,'type':'pipe','params':{'length_m':.3,'outer_diameter_m':.00635,
            'wall_thickness_m':.000762,'T_a_K':20,'T_b_K':20}}

def wires(pairs):
    return [{'id':str(i),'from':dict(zip(('component','port'),a.split(':'))),
             'to':dict(zip(('component','port'),b.split(':')))} for i,(a,b) in enumerate(pairs)]

class AnalyticalBackend:
    def __init__(self,**kw):self.calls=[]
    def __enter__(self):return self
    def __exit__(self,*args):pass
    def evaluate_pipe(self,params,pa,pb,ta,tb,**kw):
        cid=params['_component_id'];self.calls.append((cid,pa,pb))
        c={'p1':2e-10,'p2':3e-10,'p3':4e-10,'pipe1':1e-8,'pipe2':1e-10}[cid]
        mass=c*(pa-pb)*(pa+pb)
        return {'backend':'analytical test backend','mass_flow':mass,'mass_flow_kg_s':mass,
                'heat_a_W':.03,'heat_b_W':-.03-mass*CP*(tb-ta),
                'profile':[{'x':0.,'temperature':ta,'pressure':pa},
                           {'x':params['length_m'],'temperature':tb,'pressure':pb}],
                'diagnostics':{},'warnings':[]}


class HydraulicSeedTests(unittest.TestCase):
    def test_isothermal_analytic_limit_zero_flow_and_reversal(self):
        params=pipe('x')['params'];d=.00635-2*.000762
        exact=math.pi*d**4*(100**2-10**2)/(256*R*.3*gas_viscosity(20)*20)
        self.assertAlmostEqual(hydraulic_pipe_seed(params,100,10,20,20)['mass_flow']/exact,1,places=12)
        self.assertAlmostEqual(hydraulic_pipe_seed(params,10,100,20,20)['mass_flow']/-exact,1,places=12)
        self.assertEqual(hydraulic_pipe_seed(params,10,10,20,20)['mass_flow'],0)

    def test_strong_advection_and_prescribed_wall_need_no_thermal_solve(self):
        params={'length_m':.5,'outer_diameter_m':.0635,'wall_thickness_m':.00762}
        with patch('models.fluid_eda.thermal._temperature_profile',side_effect=AssertionError('Quick thermal solve called')):
            for extra in [{},{'wall_boundary':'prescribed','wall_profile':[{'length_m':.5,'T_start_K':1,'T_end_K':4}]}]:
                row=hydraulic_pipe_seed({**params,**extra},15.57,.5,1,4,profiles=True)
                self.assertTrue(math.isfinite(row['mass_flow']))
                self.assertGreater(row['mass_flow'],0)
                self.assertNotIn('profile',row)
                self.assertNotIn('heat_a_W',row)

    def test_unknown_pressure_series_parallel_and_reverse_flow(self):
        for parallel in (False,True):
            for reverse in (False,True):
                with self.subTest(parallel=parallel,reverse=reverse):
                    pa,pb=(10,100) if reverse else (100,10)
                    components=[source('a',pa),source('b',pb),{'id':'j','type':'junction','params':{'temperature_K':20}},pipe('p1'),pipe('p2')]
                    pairs=[('a:a','p1:a'),('p1:b','j:a'),('j:a','p2:a'),('p2:b','b:a')]
                    if parallel:components.append(pipe('p3'));pairs.extend([('j:a','p3:a'),('p3:b','b:a')])
                    project={'settings':{'engine':'comsol'},'components':components,'wires':wires(pairs)}
                    original=deepcopy(project);backend=AnalyticalBackend()
                    with tempfile.TemporaryDirectory() as path, patch('models.fluid_eda.thermal.reduced_pipe',side_effect=AssertionError('Preview called')),patch('models.fluid_eda.thermal._temperature_profile',side_effect=AssertionError('Thermal BVP called')):
                        result=run_simulation(project,backend_factory=lambda **kw:backend,output_dir=path)
                    c1,c2=2e-10,7e-10 if parallel else 3e-10
                    expected_p=math.sqrt((c1*pa**2+c2*pb**2)/(c1+c2))
                    self.assertAlmostEqual(result['components']['p1']['pressure_out']/expected_p,1,places=6)
                    self.assertLess(result['summary']['mass_balance_relative'],.005)
                    self.assertGreaterEqual(len(result['comsol_coupling']['iterations']),2)
                    self.assertEqual(project,original)
                    self.assertTrue(all(e['backend']=='analytical test backend' for e in result['edges']))
                    self.assertEqual(math.copysign(1,result['components']['p1']['mass_flow']),-1 if reverse else 1)

    def test_high_advection_pump_jt_loop_completes_without_preview(self):
        project=continuous_pump_example()
        project['components'][0]['params']['temperature_K']=1.1
        project['components'][1]['params'].update(length_m=.3,T_a_K=1.1)
        project['components'][2]['params']['forward_speed_m3_s']=.1
        params=project['components'][3]['params']
        params.update(length_m=1.8,T_b_K=2.8)
        params['wall_profile'].append({'length_m':.3,'T_start_K':4.,'T_end_K':2.8})
        with tempfile.TemporaryDirectory() as path, patch('models.fluid_eda.thermal.reduced_pipe',side_effect=AssertionError('Preview called')),patch('models.fluid_eda.thermal._temperature_profile',side_effect=AssertionError('Thermal BVP called')):
            result=run_simulation(project,backend_factory=AnalyticalBackend,output_dir=path)
        self.assertTrue(result['comsol_coupling']['converged'])
        self.assertLess(result['summary']['mass_balance_relative'],.005)
        self.assertEqual(result['components']['pipe1']['backend'],'analytical test backend')
        self.assertEqual(result['components']['pipe2']['backend'],'analytical test backend')

if __name__=='__main__':unittest.main()
