"""Physics limits and integration tests for the public-data HEM mode."""
from copy import deepcopy
import json
import math
import unittest
from models.fluid_eda.two_phase_properties import Helium, PhaseDomainError
from models.fluid_eda.two_phase_capillary import solve_capillary, darcy_factor, impedance
from models.fluid_eda.two_phase_pipes import precooling_pipe, suction_pipe
from models.fluid_eda.two_phase_loop import solve_two_phase_loop
from models.fluid_eda.simulation import run_simulation
from models.fluid_eda.examples import continuous_pump_example, DEFAULTS
from models.fluid_eda.network import NetworkError
from models.fluid_eda.pot_energy import reference_shift
from models.continuous_he4 import properties as he4


class TwoPhasePropertiesTests(unittest.TestCase):
    def test_nist_saturation_check_values_and_mass_quality(self):
        he=Helium()
        # NIST IR8474 Table 4 at 2.8 K: Psat=17.562 kPa,
        # liquid/vapor molar densities=35.697/0.84784 mol/dm3.
        import CoolProp.CoolProp as CP
        p=CP.PropsSI('P','T',2.8,'Q',0,'Helium')
        self.assertAlmostEqual(p/17562,1,delta=5e-5)
        a,b=he.pq(p,0),he.pq(p,1)
        self.assertAlmostEqual(a.density/(35.697*4.002602),1,delta=3e-5)
        self.assertAlmostEqual(b.density/(.84784*4.002602),1,delta=3e-5)
        x=.3;mix=he.ph(p,(1-x)*a.enthalpy+x*b.enthalpy)
        self.assertAlmostEqual(mix.quality,x,places=8)
        self.assertAlmostEqual(mix.volume,(1-x)*a.volume+x*b.volume,places=10)
        self.assertGreater(mix.void_fraction,mix.quality)

    def test_heii_transport_is_not_normal_liquid_extrapolation(self):
        with self.assertRaises(PhaseDomainError):Helium().pt(100000,1.)
        with self.assertRaises(PhaseDomainError):Helium().pq(15.57,.5)

    def test_churchill_laminar_limit(self):
        for re in [1,100,1000]:self.assertAlmostEqual(darcy_factor(re),64/re,places=8)


class CapillaryTests(unittest.TestCase):
    def setUp(self):
        self.params=deepcopy(DEFAULTS['jt']);self.he=Helium()

    def test_small_pressure_drop_matches_poiseuille(self):
        pin,pout=120000.,119990.;state=self.he.pt(pin,4.)
        r=solve_capillary(self.params,pin,state.enthalpy,pout,cells=120)
        expected=(pin-pout)*state.density/(state.viscosity*impedance(self.params)[2])
        self.assertAlmostEqual(r['mass_flow']/expected,1,delta=.005)
        self.assertFalse(r['choked'])

    def test_flash_acceleration_total_energy_and_choked_back_pressure(self):
        state=self.he.pt(120000,4.)
        a=solve_capillary(self.params,120000,state.enthalpy,15.57,cells=640)
        b=solve_capillary(self.params,120000,state.enthalpy,100.,cells=640)
        self.assertTrue(a['choked']);self.assertGreater(a['exit_pressure_Pa'],100)
        self.assertAlmostEqual(a['mass_flow']/b['mass_flow'],1,delta=.001)
        self.assertTrue(any(p['phase']=='two_phase' for p in a['profile']))
        for point in a['profile']:
            self.assertAlmostEqual(point['enthalpy_J_kg']+.5*point['velocity']**2,state.enthalpy,delta=.002)
            self.assertGreaterEqual(point['entropy_J_kg_K'],state.entropy-1e-5)
        self.assertAlmostEqual(a['profile'][-1]['x'],self.params['length_m'],delta=1e-6)

    def test_mesh_refinement(self):
        h=self.he.pt(34000,2.8).enthalpy
        a=solve_capillary(self.params,34000,h,15.57,cells=320)
        b=solve_capillary(self.params,34000,h,15.57,cells=640)
        self.assertAlmostEqual(a['mass_flow']/b['mass_flow'],1,delta=.001)


class CoupledPhaseTests(unittest.TestCase):
    def test_precooler_condenses_and_conserves_latent_energy(self):
        params=deepcopy(continuous_pump_example()['components'][3]['params'])
        params['heat_transfer_multiplier']=10.
        r=precooling_pipe(params,4e-6,120000.,293.15,max_step_fraction=1/240)
        self.assertTrue(any(p['phase']=='two_phase' for p in r['profile']))
        self.assertLess(r['enthalpy_out_J_kg'],r['enthalpy_in_J_kg'])
        self.assertLess(r['diagnostics']['relative_energy_residual'],1e-6)
        self.assertGreaterEqual(r['temperature_out'],4.)

    def test_wall_at_saturation_does_not_impose_liquid_quality(self):
        import CoolProp.CoolProp as CP
        pressure=CP.PropsSI('P','T',4.,'Q',0,'Helium')
        params=deepcopy(continuous_pump_example()['components'][3]['params'])
        params['heat_transfer_multiplier']=10.
        result=precooling_pipe(params,4e-6,pressure,293.15,max_step_fraction=1/240)
        self.assertGreaterEqual(result['temperature_out'],4.-1e-5)
        self.assertGreaterEqual(result['enthalpy_out_J_kg'],Helium().pq(result['pressure_out'],1).enthalpy-.1)

    def test_suction_steel_and_advected_energy_are_conserved(self):
        p=continuous_pump_example()['components'][1]['params']
        r=suction_pipe(p,4e-6,float(he4.saturation_pressure(1)),1.)
        self.assertLess(r['temperature_out'],p['T_b_K'])
        self.assertLess(r['diagnostics']['relative_energy_residual'],1e-5)
        self.assertLess(r['pressure_out'],r['pressure_in'])

    def test_actual_two_phase_loop_not_a_gas_flow_postprocessor(self):
        p=continuous_pump_example();p['settings']['engine']='two_phase'
        p['components'][3]['params']['heat_transfer_multiplier']=10.
        r=run_simulation(p)
        self.assertEqual(r['engine'],'two_phase')
        flows=[e['mass_flow'] for e in r['edges']]
        self.assertTrue(all(m==flows[0] for m in flows))
        self.assertTrue(any(q['phase']=='two_phase' for q in r['components']['pipe2']['profile']))
        e=r['components']['pot']['pot_energy'];self.assertTrue(e['available'])
        self.assertGreater(e['net_cooling_power_W'],0)
        self.assertEqual(e['inventory_rate_kg_s'],0)
        self.assertEqual(e['net_cooling_power_W'],-e['energy_residual_W'])
        self.assertAlmostEqual(e['net_cooling_power_W'],flows[0]*(he4.saturated_vapor_enthalpy(1)-r['components']['jt1']['enthalpy_in_J_kg']))
        self.assertLess(r['diagnostics']['max_constitutive_residual'],.002)
        self.assertTrue(r['components']['jt1']['diagnostics']['choked'])
        json.dumps(r,allow_nan=False)
        self.assertEqual(p['components'][3]['params']['T_b_K'],4.)

    def test_unsupported_topology_has_component_ids(self):
        p=continuous_pump_example();p['components'][-1]['params']['law']='orifice'
        with self.assertRaises(NetworkError) as cm:solve_two_phase_loop(p)
        self.assertIn('jt1',cm.exception.component_ids)

    def test_geometry_mode_ignores_stale_measured_impedance(self):
        p=continuous_pump_example();p['settings']['engine']='two_phase'
        jt=p['components'][-1]['params'];jt['impedance_mode']='geometry'
        jt['impedance_m_inv3']=1e3  # Must not become the network target.
        r=solve_two_phase_loop(p,pipe_cells=160,capillary_cells=320)
        actual=r['components']['jt1']['diagnostics']['impedance']
        geometric=128*jt['length_m']/(math.pi*jt['diameter_m']**4)
        self.assertAlmostEqual(actual/geometric,1,delta=.002)

    def test_invalid_phase_parameters_are_attributed(self):
        for cid,field,value in [('pipe2','heat_transfer_multiplier',0),('jt1','roughness_m',-1)]:
            p=continuous_pump_example()
            next(c for c in p['components'] if c['id']==cid)['params'][field]=value
            with self.assertRaises(NetworkError) as cm:solve_two_phase_loop(p)
            self.assertIn(cid,cm.exception.component_ids)

    def test_cancel_before_numerical_work(self):
        import threading
        from models.fluid_eda.comsol_backend import ComsolCancelled
        event=threading.Event();event.set()
        p=continuous_pump_example();p['settings']['engine']='two_phase'
        with self.assertRaises(ComsolCancelled):run_simulation(p,cancel_event=event)
