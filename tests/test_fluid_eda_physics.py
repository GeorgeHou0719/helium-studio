"""Independent limits and conservation checks for the preview/restriction."""
import math
import unittest
import numpy as np
from models.fluid_eda.thermal import reduced_pipe,wall_segments
from models.fluid_eda.jt import evaluate_restriction,CRITICAL_RATIO
from models.fluid_eda.examples import continuous_pump_example
from models.comsol_pipe.helium_properties import gas_viscosity, SPECIFIC_GAS_CONSTANT_J_KG_K as R, CP_J_KG_K as CP


class PipePhysicsTests(unittest.TestCase):
    def test_isothermal_compressible_poiseuille(self):
        p={'length_m':.7,'outer_diameter_m':.008,'wall_thickness_m':.001}
        r=reduced_pipe(p,110,100,30,30,profiles=True)
        exact=math.pi*.006**4*(110**2-100**2)/(256*gas_viscosity(30)*R*30*.7)
        # The preview interpolates its dense viscosity table (2 ppm allowance).
        self.assertLess(abs(r['mass_flow']/exact-1),2e-6)
        np.testing.assert_allclose([v['temperature'] for v in r['profile']],30)
        self.assertEqual(r['heat_environment_W'],0)

    def test_reversal_and_energy(self):
        a=reduced_pipe({},15,.5,1,293.15,profiles=True)
        b=reduced_pipe({},.5,15,293.15,1,profiles=True)
        self.assertAlmostEqual(a['mass_flow']/-b['mass_flow'],1,places=10)
        np.testing.assert_allclose([v['temperature'] for v in a['profile']],
                                   [v['temperature'] for v in b['profile']][::-1],atol=1e-9)
        self.assertLess(abs(a['heat_a']+a['heat_b']+a['mass_flow']*CP*(293.15-1)),1e-12)
        self.assertAlmostEqual(a['heat_a'],b['heat_b'],places=10)

    def test_no_flow_conduction(self):
        r=reduced_pipe({},5,5,1,293.15,profiles=True)
        self.assertEqual(r['mass_flow'],0)
        self.assertGreater(r['heat_a'],.13)
        self.assertLess(r['heat_a'],.15)
        self.assertAlmostEqual(r['heat_a'],-r['heat_b'],places=12)
        self.assertTrue(all(v['pressure']==5 for v in r['profile']))

    def test_prescribed_profile_exact_knots_and_global_energy(self):
        p=continuous_pump_example()['components'][3]['params']
        r=reduced_pipe(p,1865,1864.9,293.15,4.,profiles=True)
        z=np.array([v['x'] for v in r['profile']]);t=np.array([v['temperature'] for v in r['profile']])
        np.testing.assert_allclose(np.interp([0,.3,.6,.9,1.2,1.5],z,t),[293.15,293.15,35,35,4,4],atol=1e-10)
        residual=r['heat_a']+r['heat_b']+r['heat_environment_W']+r['mass_flow']*CP*(4.-293.15)
        self.assertLess(abs(residual),1e-12)
        self.assertGreater(r['heat_environment_W'],0)
        p['length_m']=.8
        with self.assertRaises(ValueError):wall_segments(p,293.15,4.)


class RestrictionTests(unittest.TestCase):
    def test_choking_independent_downstream_pressure(self):
        a=evaluate_restriction({},1000,100,30,30)
        b=evaluate_restriction({},1000,10,30,30)
        self.assertEqual(a['mass_flow'],b['mass_flow'])
        self.assertTrue(a['diagnostics']['choked'])
        self.assertAlmostEqual(CRITICAL_RATIO,.4871392896,places=8)

    def test_choked_to_subcritical_continuity_and_zero_flow(self):
        lo=evaluate_restriction({},1000,1000*CRITICAL_RATIO*(1-1e-6),30,30)
        hi=evaluate_restriction({},1000,1000*CRITICAL_RATIO*(1+1e-6),30,30)
        self.assertAlmostEqual(lo['mass_flow']/hi['mass_flow'],1,places=10)
        self.assertEqual(evaluate_restriction({},10,10,3,3)['mass_flow'],0)

    def test_reverse_uses_upstream_temperature(self):
        a=evaluate_restriction({},100,5,4,1)
        b=evaluate_restriction({},5,100,1,4)
        self.assertEqual(a['mass_flow'],-b['mass_flow'])
        self.assertEqual(b['outlet_fluid_temperature_K'],4)
        self.assertEqual(b['profile'][0]['pressure'],5)

    def test_long_capillary_analytic_limit(self):
        p={'law':'capillary','diameter_m':.0001,'length_m':.05}
        out=evaluate_restriction(p,100,99,10,10)
        analytic=math.pi*.0001**4*(100**2-99**2)/(256*gas_viscosity(10)*R*10*.05)
        self.assertAlmostEqual(out['mass_flow']/analytic,1,places=12)

    def test_measured_capillary_impedance_controls_flow_and_reversal(self):
        from models.fluid_eda.examples import DEFAULTS
        p=dict(DEFAULTS['jt'])
        self.assertEqual((p['material'],p['diameter_m'],p['length_m']),('CuNi',.0001,.1))
        self.assertEqual(p['impedance_m_inv3'],3.2e16)
        result=evaluate_restriction(p,100,90,10,5)
        expected=(100**2-90**2)/(2*gas_viscosity(10)*R*10*3.2e16)
        self.assertAlmostEqual(result['mass_flow']/expected,1,places=12)
        self.assertEqual(result['mass_flow'],-evaluate_restriction(p,90,100,5,10)['mass_flow'])
        self.assertEqual(evaluate_restriction(p,100,100,10,10)['mass_flow'],0)
        self.assertAlmostEqual(evaluate_restriction(dict(p,impedance_m_inv3=6.4e16),100,90,10,5)['mass_flow']/expected,.5)
        # Measured impedance overrides nominal dimensions, without double counting.
        self.assertEqual(evaluate_restriction(dict(p,length_m=.2),100,90,10,5)['mass_flow'],result['mass_flow'])
        geometric=evaluate_restriction(dict(p,impedance_mode='geometry'),100,90,10,5)
        self.assertAlmostEqual(result['mass_flow']/geometric['mass_flow'],128*.1/(math.pi*.0001**4)/3.2e16)

    def test_measured_impedance_validation(self):
        from models.fluid_eda.examples import DEFAULTS, continuous_pump_example
        from models.fluid_eda.network import validate_network, NetworkError
        for bad in [0,-1,float('nan'),float('inf'),True,None,'3.2e16']:
            with self.subTest(value=bad):
                p=dict(DEFAULTS['jt'],impedance_m_inv3=bad)
                with self.assertRaises(ValueError):evaluate_restriction(p,100,90,10,10)
                project=continuous_pump_example();project['components'][-1]['params']=p
                with self.assertRaises(NetworkError) as caught:validate_network(project)
                self.assertIn('jt1',caught.exception.component_ids)

if __name__=='__main__':unittest.main()
