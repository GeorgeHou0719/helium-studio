# User guide

## Edit a circuit

1. Click or drag a component from the palette onto the canvas.
2. Click one port and then another to connect them. Use junctions for branches.
3. Select a component to edit its parameters in the inspector.
4. Select a solver and run the simulation. Select individual components to view
   results and pipe profiles.
5. Export the circuit JSON to keep a portable copy; download results separately.

Wires merge ports into nodes without adding pressure loss, thermal resistance,
or volume. The gas-network engines permit branches. The two-phase engine currently
requires the five-component example loop connected in the forward A→B direction.

## Canvas controls

| Action | Control |
| --- | --- |
| Zoom | Mouse wheel or zoom buttons |
| Pan | Space + drag or middle mouse button |
| Rotate selected component | `R` or toolbar |
| Flip selected component horizontally | `X`, toolbar, or inspector |
| Delete selected component | `Delete` |
| Undo / redo | `Ctrl+Z` / `Ctrl+Shift+Z` |

Mirroring and rotation affect appearance; port A/B identities and physical
orientation remain unchanged. Text fields keep normal editing behavior.
Short explanations and model details are available through the `?` controls.

## Units and pump inputs

The UI displays millimetres, litres/second, kelvin, pascals, and cooling power in
milliwatts. Circuit JSON and Python APIs use SI: metres, m³/s, kelvin, pascals,
kg/s, and watts. All pressures are absolute.

Pumping speed is volume flow at the inlet state, not standard-volume flow or
mass flow. The forward–backflow law accepts nominal speed and compression plus
optional pressure-dependent tables. Its displayed forward/backward contributions
are nonnegative; their difference is signed. Negative net flow is possible in
the gas-network engines. The two-phase loop currently accepts forward circulation.

## Interpreting pot results

The pot temperature is a user trial value, and pot pressure follows the
saturation curve. The application does not automatically search for an
equilibrium temperature or integrate changing liquid inventory.

The pot energy balance is

```text
Qcool = mout * (hg - hl) - sum(min * (hin - hl))
heat imbalance = applied heat load - Qcool
inventory rate = sum(min) - mout
```

Here `hg` and `hl` are saturated vapor/liquid specific enthalpies at the pot
temperature, and `hin` is incoming feed energy per unit mass. J–T jet kinetic
energy recovery is included in the two-phase feed's total incoming energy.

- Positive **Net cooling power** means refrigeration is available; a negative
  value means the incoming circulation heats the pot.
- Positive **Heat imbalance** means net heating; negative means net cooling.
- Nonzero inventory rate means the liquid inventory is not steady. In the
  two-phase loop, a common mass flow enforces equal circulating inflow/outflow.
- A general fixed-pressure source can supply or absorb mass. Its **Net mass
  supply** is therefore not evidence that all branches have identical flows.

A trial temperature is an equilibrium candidate when the heat imbalance and
inventory rate are near zero and the numerical residuals are acceptable. Choose
tolerances appropriate to the intended heat load and measurement uncertainty.
Convergence alone does not establish stability or experimental accuracy. A
stable thermal root normally has heating just below it and cooling just above it,
with other boundaries held fixed and positive effective heat capacity.

Solid conduction is excluded between components, as an explicit thermal-isolation
assumption. Gas enthalpy is retained. Each pipe can still conduct heat internally;
within-pipe heat exchange changes its fluid outlet state. There is no displayed
wall-cooling contribution added to the pot capacity.

## Errors and model limits

An error report names the affected component and links to its parameters.
Expand technical details for boundary values and solver messages. A successful
numerical solve can still carry warnings about Knudsen number, Reynolds number,
Mach number, phase assumptions, or pump curve extrapolation.

COMSOL startup or solve failures do not automatically switch engines. Quick
preview can fail in strong-advection regimes. Select the appropriate engine and
check the physical assumptions rather than treating a different solver as proof
of model validity. Unsupported two-phase topologies and He II transport states
produce explicit errors.

Changing the schematic while a calculation is running can leave a result tied
to the submitted version. The report identifies this situation; rerun to update.
Drafts in browser local storage are specific to the browser and origin, including
the port number. Use JSON export before changing browsers or ports.
