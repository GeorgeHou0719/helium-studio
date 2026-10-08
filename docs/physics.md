# Physics and numerical scope

Helium Studio separates numerical conservation, empirical closures, property
validity, and experimental validation. All equations below use SI units and
absolute pressure. The two-phase extensions have a
[separate detailed specification](two-phase-model.md).

## Gas-network conservation

Connected ports share a pressure node. Unknown positive pressures are solved
using logarithmic variables. Each non-source node satisfies
`sum(signed branch mass flow) = 0`. Pressure sources can supply or absorb mass.
Pipe elements are solved independently; entrance interactions between neighboring
pipe segments are omitted. Temperature boundaries are prescribed rather than
obtained from finite-capacity refrigerator stages.

COMSOL network coupling begins with a hydraulic pressure estimate, then corrects
branch conductance using actual COMSOL evaluations. The final network mass
residual threshold is 0.5%; individual pipe conservation is checked separately.
COMSOL startup does not require the quick thermal solver to converge.

## Pipes

The COMSOL template is steady, 2D axisymmetric, compressible ideal helium with
laminar no-slip flow, gas/steel conjugate heat transfer, and axial steel
conduction. End pressure and gas/steel end temperatures are prescribed. The outer
wall is insulated or follows a continuous piecewise temperature profile.

The quick model assumes equal gas and wall temperature across each section and
neglects axial inertia. For an insulated exterior its heat equation is

```text
K(T) = Ag * kg(T) + As * ks(T)
K(T) * dT/dz = mdot * cp * T + C
mdot = pi * Di^4 * (Pa^2 - Pb^2) / (256 * RHe * integral(mu(T)*T dz))
```

For a prescribed wall, quick preview assumes gas temperature follows the wall
profile. It cannot represent finite gas/wall temperature differences. Pressure
squared varies with accumulated `mu*T` resistance; density is `P/(RHe*T)`.

The new two-phase mode has a different suction closure: separate gas/wall
temperatures, gas enthalpy advection, and axial steel conduction, with calculated
gas outlet temperature. Its precooler conserves real-fluid enthalpy while
resolving finite heat transfer and possible condensation.

## Pumps

Constant-speed pumping uses inlet-volume speed `S`:

```text
mdot = Pin * S / (RHe * Tin)
```

The custom expression returns mass flow in kg/s and can use `Pin`, `Pout`,
`Tin`, `S`, and `R`. The legacy expression model rejects negative flow.

The additional forward–backflow law is

```text
Qforward = Sf(Pin) * Pin
Qback = Sf(Pin) * Pout / K0(Pout)
mdot = (Qforward - Qback) / (RHe * Tin)
```

`Q` is inlet-temperature-referenced throughput in Pa·m³/s. `K0` is zero-flow
compression; with constant coefficients, zero net flow occurs at `Pout/Pin=K0`.
The net flow is signed. Both throughput terms use `Tin`; this is a reduced
isothermal-reference law, not a calibrated nonisothermal pump model.

Optional speed and compression tables use log-pressure/log-value interpolation.
The nearest endpoint is held outside the table, with a result warning. Each
table replaces its nominal constant. User-supplied helium inlet/outlet-pressure
and throughput limits are checked at the accepted operating point. Unknown
limits remain unset. Pump discharge temperature is prescribed; shaft work and
heat rejection are not predicted.

The HiPace 300 M preset uses 240 L/s and `K0=1e8`. Published documentation also
reports 215 L/s and compression greater than 1e8 for a TM700 configuration.
Confirm the hardware revision. A maximum or nominal compression ratio and
speed do not establish a complete two-pressure performance map or valid
foreline-pressure envelope. The separable `Sf(Pin), K0(Pout)` representation
and operating limits require measured data.

## Restrictions

The legacy capillary law follows fully developed, isothermal ideal-gas
Poiseuille flow:

```text
z_geometry = 128 * L / (pi * D^4)
mdot = (Pu^2 - Pd^2) / (2 * mu(Tu) * RHe * Tu * z)
```

Measured mode uses the entered `z` directly; geometry mode derives it from
length and bore. The default CuNi restriction has 0.10 mm bore, 100 mm length,
and `z=3.2e16 m^-3`. This is a supplied literature-based coefficient, not a
measurement of a newly built device. CuNi is descriptive metadata; no capillary
wall heat-transfer model is inferred from it.

The alternative legacy orifice law uses the standard isentropic ideal-gas nozzle
mass-flow relation with `gamma=5/3`, empirical discharge coefficient, and choking
at `Pd/Pu = (2/(gamma+1))^(gamma/(gamma-1))`. Neither gas-only hydraulic law
models liquid/two-phase resistance. Pot enthalpy postprocessing does not change
that hydraulic limitation.

The two-phase capillary instead conserves total enthalpy, includes acceleration
in momentum balance, and detects homogeneous-equilibrium critical flow. Its
measured impedance sets an effective hydraulic length, not a measured two-phase
friction correlation. See [the capillary equations](two-phase-model.md#capillary-momentum-and-energy).

## Pot energy and validity

All modes retain incoming feed enthalpy and outgoing saturated-vapor enthalpy in
the pot balance, excluding intercomponent solid conduction. See the
[pot-result guide](user-guide.md#interpreting-pot-results). Real-feed properties
are postprocessed for legacy gas hydraulics; the new mode couples phase changes
into the loop itself.

Gas viscosity/conductivity use public dilute-helium data. SS304 conductivity
uses the NIST fit: equation range 1–300 K, underlying data 4–300 K. The quoted
fit error over the data range is not a validated uncertainty at 1 K.
Source-specific ranges remain in force even where blanket adapter windows
have been removed.

Knudsen, Reynolds, and Mach diagnostics screen the approximations. Continuum
no-slip friction/heat transfer need correction or calibration as rarefaction
increases. No phase slip, metastable flashing, superfluid film flow, or He II
counterflow is resolved. In the two-phase mode, both pipe heat-transfer
closures require laminar accepted operating points. Numerical convergence and
conservation do not establish experimental accuracy.
