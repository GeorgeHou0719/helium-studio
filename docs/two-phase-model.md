# Continuous-pot two-phase mode

Select **Two-phase pot (1D)** in Simulation model. Existing COMSOL
and quick-preview modes remain available and unchanged. This additional mode
uses public property data and a Python 1D solver; it does not turn the existing
COMSOL gas template into a two-phase model.

Status: experimental engineering model, with numerical verification. No claim
of experimental validation of helium microcapillary flow or condensation.

## Supported circuit and boundaries

The first implementation accepts the current five-component loop, connected
from A to B: saturated pot -> adiabatic suction pipe -> pump -> prescribed-wall
precooling pipe -> capillary -> pot. Schematic mirroring/rotation still works.
Other graphs and reverse net circulation are outside this mode's scope; errors
identify the relevant components. Legacy modes retain their existing topology
and pump support. Pump forward/backflow, constant-speed and expression laws are
reused, including configured operating-envelope checks for forward/backflow.

Pot temperature is a trial value; pressure follows its saturation curve.
Outgoing pot vapor is saturated. Pipe endpoint temperatures are steel boundary
conditions in this mode. The fluid outlet temperature is computed from finite
heat exchange rather than forced to the wall temperature. Pump discharge
fluid temperature remains prescribed. No solid conduction is passed between
components or added to pot duty. Within the suction pipe, axial steel conduction
and separate fluid/wall temperatures are retained. No wall-cooling duty is shown.

All four flow elements share one solved mass flow. The other loop unknown is
pump discharge pressure. Pump flow and capillary impedance residuals are solved
together; precooler phase change changes the hydraulic operating point. This
is not a liquid-fraction calculation appended to a gas-only network solution.

## Thermodynamics and sources

Normal helium uses CoolProp HEOS with a constant enthalpy reference shift to
DB98 saturated liquid at 4.2 K. Pot liquid enthalpy, vapor pressure and latent
heat reuse Donnelly and Barenghi (1998). Reference alignment does not establish
independent accuracy of the joined sources.

Within the saturation dome, pressure and specific enthalpy define vapor mass
fraction x = (h - h_l)/(h_g - h_l). Mixture specific volume is
v = (1-x)/rho_l + x/rho_g. Vapor volume fraction is x/(rho_g*v).
Out-of-dome states remain liquid, vapor or supercritical fluid; x is not clipped
into [0,1] to manufacture a two-phase state.

The normal-fluid EOS is not a superfluid EOS. Normal-liquid transport is never
extrapolated into He II. If a capillary reaches the unsupported region before
choking, the circuit returns a component-specific property/model error. If it
chokes while still in the supported normal-fluid domain, downstream equilibrium
at the 1 K pot is evaluated using DB98, without integrating fictitious He II
hydraulics inside the capillary.

The new suction closure currently uses the existing public gas transport table
from 1 K, and the existing 304 steel fit through 300 K. Its colder/hotter boundary
extension needs additional property work. This is a property scope of this
additional Python model, not a reinstated temperature restriction on COMSOL.

Primary/reference sources:

- [Donnelly and Barenghi, 1998](https://doi.org/10.1063/1.556028), saturation
  quantities; implemented property splines; publication PDF is not bundled.
- [NIST IR 8474, 2023](https://nvlpubs.nist.gov/nistpubs/ir/2023/NIST.IR.8474.pdf),
  normal-fluid EOS, Table 4 check values, and explicit exclusion of He II.
- [CoolProp helium documentation](https://coolprop.org/fluid_properties/fluids/Helium.html),
  implementation and transport references. The Python requirement is CoolProp >=8,<9.
- [Hurly and Moldover, 2000](https://nvlpubs.nist.gov/nistpubs/jres/105/5/j55hur.pdf),
  dilute-gas transport, reused from existing project functions.
- [NBS Technical Note 633, 1973](https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote633.pdf),
  homogeneous equilibrium critical-flow framework and its limitations for cryogens.
- [Churchill, 1977](https://files.engineering.com/files/85c0f3a6-a102-4a22-9d35-f15858c0dd2b/CEM_-_Friction-factor_equation_(1977).pdf),
  all-regime pipe friction; implementation uses the Darcy-factor convention.
- [NASA CR-202662, Arauz, 1996/1997](https://ntrs.nasa.gov/api/citations/19970005355/downloads/19970005355.pdf),
  sections 4.4-4.5 explain pressure/enthalpy mixture states and discuss mixture
  viscosity closures. This is model background, not experimental validation of
  the present helium capillary.
- [COMSOL, Internal Film Resistance](https://doc.comsol.com/6.3/doc/com.comsol.help.pipe/pipe_ug_heattransfer.06.08.html),
  Nu=3.66 for fully developed laminar circular-tube flow. Extending this to
  condensation through homogeneous mixture transport is our engineering
  assumption, not a source-provided helium condensation correlation.

## Capillary momentum and energy

Let G = mass flow / bore area, and v = 1/rho. Static enthalpy inside the
capillary satisfies h + (G*v)^2/2 = h0. The entrance is treated as lossless and
isentropic: it consumes pressure to accelerate fluid from the upstream pipe
into the small bore. An extra empirical entrance-loss coefficient is not assumed.

The distributed momentum equation is

    dP/dz + G^2 dv/dz = -f_D G^2 v/(2 D).

The laminar Darcy factor tends to 64/Re. Churchill's all-regime expression
handles higher Reynolds numbers instead of extending a laminar law without
checking it. Default roughness is zero, an explicit smooth-wall assumption.
Mixture viscosity uses harmonic mass-fraction mixing (McAdams closure).

In measured-impedance mode, Leff = z_measured*pi*D^4/128. Thus the measured
coefficient preserves its low-speed Poiseuille meaning and sets effective
hydraulic length. It does not supply a measured two-phase friction law. In
geometry mode the entered measured value is ignored and actual length is used.

The pressure-coordinate integration ends at the receiving pressure or when
1 + G^2*(dv/dP at fixed h0) = 0. That is the HEM critical-flow condition. When
choked, lowering pot pressure does not materially increase flow. The axial
profile ends at capillary exit pressure, which can exceed pot pressure.
Downstream expansion/mixing and recovery of jet kinetic energy are included in
the pot control volume: incoming energy per unit mass is h0, not the lower
static exit enthalpy. Additional jet heat loss and metastable flashing are omitted.

## Pipes and heat-transfer uncertainty

The adiabatic suction closure solves a gas/steel conjugate boundary-value
problem. Its conserved flux is mass_flow*h_gas + axial_steel_heat_flow. A
Kirchhoff transform of steel conductivity makes the cryogenic numerical solve
well conditioned. Vapor is dilute; below the normal EOS range its enthalpy is
approximated by DB98 saturated-vapor h(T). Axial gas conduction and gas kinetic
energy are neglected. Mach > 0.3 is outside this approximation; rarefaction is
reported using Kn rather than silently treated as validated continuum flow.

The precooler marches pressure and enthalpy with conservative implicit cells:

    mass_flow * dh/dz = U_per_length * (T_wall - T_fluid).

Film and radial-steel thermal resistances are in series. The outer steel wall
profile is prescribed; its external thermal boundary absorbs the required
heat. The approximation takes the imposed profile as the steel axial thermal
boundary and does not solve a separate 2D wall field. Latent heat remains in h,
so a constant saturation temperature can accompany changing liquid fraction.

The starting film estimate is Nu=3.66 times the editable heat-transfer multiplier.
Mixture conductivity uses void-fraction parallel mixing. These are explicit
engineering closures. Default multiplier 1 is not measured. Slip, stratified
flow, film condensation, boiling correlations, entrance development and He II
film/normal-superfluid counterflow are not resolved. Both pipe heat-transfer closures are restricted to laminar flow; a component
error identifies an accepted operating point with Re > 2300.

## Pot balance and equilibrium

    Q_net = m_out*(h_g - h_l) - sum(m_in*(h_in - h_l))
    heat_imbalance = applied_heat_load - Q_net

The loop imposes m_in = m_out. Positive Q_net is available refrigeration;
negative Q_net heats the pot. With continuous circulation, a trial temperature
is an equilibrium candidate when heat_imbalance is near zero and constitutive
residuals are small. Stability and experimental correctness are separate checks.

Results show fluid temperature, pressure, velocity, phase, enthalpy, vapor mass
fraction, vapor volume fraction, capillary choking and pot net cooling. Single-
phase quality is left undefined. Technical explanations are under question-mark
help. Superheated pot-entry gas below the normal EOS range may have unresolved
entry temperature while its incoming enthalpy and energy balance remain known.

## Verification and remaining calibration

`tests/test_two_phase.py` checks NIST saturation density values, phase enthalpy
and volume identities, rejection of He II transport, the laminar friction limit,
small-pressure-drop Poiseuille flow, choking/back-pressure independence, total
enthalpy and entropy at the entrance/along the capillary, capillary mesh
refinement, condensation, pipe energy conservation, coupled loop mass/energy
balance, component identity in errors and cancellation.

Run `python -m pytest` from the repository root. Run histories, result JSON,
solver logs and local audit artifacts are intentionally excluded from this
repository. Heat-transfer multipliers are calibration inputs, not measured
coefficients established by these numerical checks.
Calibration priorities: real precooler heat transfer, helium capillary mass flow
versus pressures/enthalpy, phase slip and nonequilibrium flashing, low-pressure
suction accommodation/slip, loaded pump operating envelope, and He II channel
flow if the calculated restriction trajectory enters that regime.
