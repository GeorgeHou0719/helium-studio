"""Normal He-4 EOS and homogeneous mixture properties in SI.

Thermodynamics: CoolProp HEOS (NIST normal-fluid EOS), aligned to DB98.
Mixture closure: no velocity slip; McAdams harmonic dynamic viscosity.
Conductivity uses void-fraction parallel mixing, an engineering heat-transfer
closure, not a validated helium condensation correlation. No He II transport
is supplied. Pot He II equilibrium is evaluated separately from DB98.
"""
from dataclasses import dataclass
import math
import CoolProp.CoolProp as CP
from .pot_energy import reference_shift, NORMAL_FLUID_MIN_T_K


class PhaseDomainError(ValueError):
    pass


@dataclass(frozen=True)
class State:
    pressure: float
    temperature: float
    enthalpy: float
    density: float
    viscosity: float
    conductivity: float
    quality: float | None
    phase: str
    void_fraction: float | None = None
    entropy: float | None = None

    @property
    def volume(self):
        return 1.0 / self.density

    def point(self, x, mass=0., area=1.):
        return dict(x=float(x), pressure=self.pressure, temperature=self.temperature,
                    enthalpy_J_kg=self.enthalpy, density_kg_m3=self.density,
                    vapor_quality=self.quality, void_fraction=self.void_fraction,
                    phase=self.phase, entropy_J_kg_K=self.entropy, velocity=float(mass / (area*self.density)))


class Helium:
    def __init__(self):
        self.eos = CP.AbstractState("HEOS", "Helium")
        self.liquid = CP.AbstractState("HEOS", "Helium")
        self.vapor = CP.AbstractState("HEOS", "Helium")
        self.shift = reference_shift()
        self.minimum_temperature = NORMAL_FLUID_MIN_T_K + 1e-5
        self.critical_pressure = self.eos.p_critical()
        self.minimum_saturation_pressure = CP.PropsSI('P','T',self.minimum_temperature,'Q',0,'Helium')

    def _read(self):
        s = self.eos
        t, p, h, rho = s.T(), s.p(), s.hmass()+self.shift, s.rhomass()
        if not (self.minimum_temperature-1e-7 <= t <= 1500 and 0 < p <= 1e9):
            raise PhaseDomainError('Normal-fluid transport reached the He II/unsupported property domain. A superfluid flow law is required; no normal-liquid extrapolation was used.')
        q = s.Q()
        if 0 <= q <= 1:
            self.liquid.update(CP.PQ_INPUTS,p,0)
            self.vapor.update(CP.PQ_INPUTS,p,1)
            a,b = self.liquid,self.vapor
            volume = (1-q)/a.rhomass() + q/b.rhomass()
            alpha = (q/b.rhomass())/volume
            mu = 1/((1-q)/a.viscosity() + q/b.viscosity())
            k = (1-alpha)*a.conductivity()+alpha*b.conductivity()
            return State(p,t,h,1/volume,mu,k,float(q),'two_phase',alpha,float(s.smass()))
        phase = 'liquid' if s.phase() == CP.iphase_liquid else 'gas' if s.phase() in (CP.iphase_gas,CP.iphase_supercritical_gas) else 'dense_fluid' if s.phase() == CP.iphase_supercritical_liquid else 'supercritical'
        return State(p,t,h,rho,s.viscosity(),s.conductivity(),None,phase,entropy=float(s.smass()))

    def pt(self, p, t):
        if not (math.isfinite(p) and math.isfinite(t) and p > 0 and t >= self.minimum_temperature):
            raise PhaseDomainError('Normal-fluid P,T state is outside its domain; He II liquid transport is not implemented.')
        try:
            self.eos.update(CP.PT_INPUTS,p,t)
            return self._read()
        except ValueError as error:
            raise PhaseDomainError(str(error)) from error

    def ph(self, p, h):
        try:
            self.eos.update(CP.HmassP_INPUTS,h-self.shift,p)
            return self._read()
        except ValueError as error:
            raise PhaseDomainError('Cannot resolve this P,h state with the normal-helium EOS. It may enter He II: '+str(error)) from error

    def pq(self,p,q):
        if not self.minimum_saturation_pressure <= p < self.critical_pressure:
            raise PhaseDomainError('Normal-fluid saturation pressure is outside the supported range.')
        self.eos.update(CP.PQ_INPUTS,p,q)
        return self._read()

    def ps(self,p,entropy):
        try:
            self.eos.update(CP.PSmass_INPUTS,p,entropy)
            return self._read()
        except ValueError as error:
            raise PhaseDomainError('Isentropic entrance trajectory leaves the normal-fluid domain: '+str(error)) from error
