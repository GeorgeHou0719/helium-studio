"""Independent pot first-law limits, common reference, flash and source checks."""
from copy import deepcopy
import math
import unittest
from unittest.mock import patch
from CoolProp.CoolProp import PropsSI
from models.continuous_he4 import properties as he4
from models.fluid_eda.pot_energy import pot_balance, feed_state, flash_at_pot, reference_shift, apply_pot_balances


def fixture(reverse=False, temperature=4., pressure=34188.813831):
    project={'components':[{'id':'pot','type':'source','params':{'mode':'saturation','temperature_K':1.,'heat_load_W':.001}}]}
    result={'nodes':[{'id':'potnode','temperature':1.,'pressure':float(he4.saturation_pressure(1.)),'source_ids':['pot']},
                     {'id':'feednode','temperature':temperature,'pressure':pressure,'source_ids':[]},
                     {'id':'warm','temperature':293.15,'pressure':10.,'source_ids':[]}],
            'components':{'pot':{'id':'pot','type':'source','node':'potnode'}},'warnings':[]}
    jt={'id':'jt','type':'jt','node_a':'potnode' if reverse else 'feednode',
        'node_b':'feednode' if reverse else 'potnode','mass_flow':-2e-6 if reverse else 2e-6,
        'warnings':[],'profile':[{'x':0,'temperature':temperature},{'x':1,'temperature':temperature}]}
    pipe={'id':'suction','type':'pipe','node_a':'potnode','node_b':'warm','mass_flow':2e-6,
          'heat_a_W':11.33,'heat_b_W':-14.,'warnings':[]}
    result['edges']=[jt,pipe];result['components'].update(jt=jt,suction=pipe)
    return result,project


class PotEnergyTests(unittest.TestCase):
    def test_saturated_liquid_feed_gives_latent_cooling_and_load_balance(self):
        m=2e-6;hl=float(he4.saturated_liquid_enthalpy(1));latent=float(he4.latent_heat(1))
        streams=[{'mass_flow_kg_s':m,'enthalpy_J_kg':hl}]
        r=pot_balance(1,streams,m,m*latent)
        self.assertAlmostEqual(r['net_cooling_power_W'],m*latent)
        self.assertAlmostEqual(r['energy_residual_W'],0)
        self.assertEqual(r['inventory_rate_kg_s'],0)
        self.assertFalse(r['solid_conduction_included'])

    def test_vapor_feed_zero_cooling_hot_feed_heats_and_load_sign(self):
        m=1e-6;hg=float(he4.saturated_vapor_enthalpy(1))
        r=pot_balance(1,[{'mass_flow_kg_s':m,'enthalpy_J_kg':hg}],m)
        self.assertAlmostEqual(r['net_cooling_power_W'],0)
        r=pot_balance(1,[{'mass_flow_kg_s':m,'enthalpy_J_kg':hg+1000}],m,.002)
        self.assertAlmostEqual(r['net_cooling_power_W'],-.001)
        self.assertAlmostEqual(r['energy_residual_W'],.003)

    def test_inventory_is_reported_and_enthalpy_reference_cancels(self):
        streams=[{'mass_flow_kg_s':2e-6,'enthalpy_J_kg':9000.}]
        before=pot_balance(1,streams,1.9e-6)
        offset=54321.;hl=float(he4.saturated_liquid_enthalpy(1));hg=float(he4.saturated_vapor_enthalpy(1))
        with patch.object(he4,'saturated_liquid_enthalpy',return_value=hl+offset),patch.object(he4,'saturated_vapor_enthalpy',return_value=hg+offset):
            after=pot_balance(1,[{'mass_flow_kg_s':2e-6,'enthalpy_J_kg':9000+offset}],1.9e-6)
        self.assertAlmostEqual(before['net_cooling_power_W'],after['net_cooling_power_W'])
        self.assertAlmostEqual(before['inventory_rate_kg_s'],1e-7)
        self.assertFalse(pot_balance(1,[],1e-6)['circulating'])

    def test_common_reference_and_overlap_with_independent_saturation_data(self):
        self.assertAlmostEqual(PropsSI('Hmass','T',4.2,'Q',0,'Helium')+reference_shift(),he4.saturated_liquid_enthalpy(4.2))
        # Independent sources are not identical: record/allow their actual
        # overlap discrepancy rather than assuming exact thermodynamic identity.
        for t in [2.18,2.5,2.8,3.,4.,4.2]:
            for q,fn in [(0,he4.saturated_liquid_enthalpy),(1,he4.saturated_vapor_enthalpy)]:
                self.assertLess(abs(PropsSI('Hmass','T',t,'Q',q,'Helium')+reference_shift()-fn(t)),150.)

    def test_actual_gas_and_liquid_feed_states_are_not_interchanged(self):
        hg,pg,_=feed_state(4.,34188.813831)
        hl,pl,_=feed_state(4.,120000.)
        self.assertEqual(pg,'gas');self.assertEqual(pl,'liquid')
        self.assertGreater(hg,he4.saturated_vapor_enthalpy(1))
        self.assertLess(hl,he4.saturated_vapor_enthalpy(1))
        self.assertEqual(flash_at_pot(hg,1,float(he4.saturation_pressure(1)))['phase'],'superheated_vapor')
        flash=flash_at_pot(hl,1,float(he4.saturation_pressure(1)))
        self.assertEqual(flash['phase'],'two_phase');self.assertTrue(0<flash['vapor_quality']<1)
        self.assertAlmostEqual(flash['vapor_quality']*he4.saturated_vapor_enthalpy(1)+(1-flash['vapor_quality'])*he4.saturated_liquid_enthalpy(1),hl)

    def test_conduction_is_excluded_jt_retains_enthalpy_and_reverse_orientation(self):
        for reverse in (False,True):
            with self.subTest(reverse=reverse):
                result,project=fixture(reverse)
                original=deepcopy(result)
                apply_pot_balances(result,project)
                e=result['components']['pot']['pot_energy'];jt=result['components']['jt']
                self.assertTrue(e['available']);self.assertLess(e['net_cooling_power_W'],0)
                self.assertEqual(jt['enthalpy_in_J_kg'],jt['enthalpy_out_J_kg'])
                self.assertGreater(jt['outlet_fluid_temperature_K'],1)
                self.assertEqual(jt['temperature_in' if reverse else 'temperature_out'],jt['outlet_fluid_temperature_K'])
                self.assertEqual(jt['mass_flow'],original['components']['jt']['mass_flow'])
                original['components']['suction']['heat_a_W']=99999.
                apply_pot_balances(original,project)
                self.assertEqual(original['components']['pot']['pot_energy'],e)

    def test_liquid_feed_is_conditional_on_gas_hydraulics(self):
        result,project=fixture(temperature=2.8)
        apply_pot_balances(result,project)
        e=result['components']['pot']['pot_energy']
        self.assertGreater(e['net_cooling_power_W'],0)
        self.assertTrue(any('gas hydraulics' in s for s in e['notes']))
        self.assertEqual(result['components']['jt']['temperature_out'],1)

    def test_missing_feed_phase_and_unsupported_properties_are_visible(self):
        with self.assertRaises(ValueError):feed_state(1.,float(he4.saturation_pressure(1.)))
        result,project=fixture(temperature=1600.)
        apply_pot_balances(result,project)
        self.assertFalse(result['components']['pot']['pot_energy']['available'])
        self.assertTrue(any('unavailable' in s for s in result['warnings']))
        for bad in [float('nan'),-1.,float('inf')]:
            with self.subTest(bad=bad),self.assertRaises(ValueError):pot_balance(1,[],bad)

if __name__=='__main__':unittest.main()
