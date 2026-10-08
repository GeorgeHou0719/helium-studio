# Helium Studio

[![Tests](https://github.com/GeorgeHou0719/helium-studio/actions/workflows/tests.yml/badge.svg)](https://github.com/GeorgeHou0719/helium-studio/actions/workflows/tests.yml)

A local schematic editor and steady-state simulator for helium fluid circuits.
Place pipes, pumps, pressure sources, J–T restrictions, and junctions on a canvas;
connect their ports; inspect pressure, temperature, flow, and pot cooling results.

The application runs in your browser with a Python backend. The editor needs no
Node.js build step. Python solvers run without COMSOL; the separate COMSOL engine
requires a licensed installation.

## Features

- Circuit-style canvas with drag-and-drop components, wiring, rotation, horizontal
  flipping, undo/redo, local draft saving, and JSON import/export.
- Stainless-steel pipes with insulated or prescribed-temperature outer walls.
- Constant-speed, custom-expression, and forward–backflow pump laws.
- Constant-pressure, temperature-dependent, and saturated-vapor pressure sources.
- Capillary and orifice restrictions, with measured capillary impedance support.
- Component-specific errors, axial profiles, regime warnings, and downloadable results.
- A continuous one-kelvin-pot example and an experimental two-phase loop solver.

## Quick start

Requires **Python 3.11 or newer**. Start from a fresh checkout:

```sh
git clone https://github.com/GeorgeHou0719/helium-studio.git
cd helium-studio
python -m venv .venv
```

Activate the environment on Windows PowerShell:

```powershell
.\.venv\Scripts\Activate.ps1
```

Or on macOS/Linux:

```sh
source .venv/bin/activate
```

Install and launch:

```sh
python -m pip install -e .
python -m models.fluid_eda.app
```

Open **http://127.0.0.1:5004**. Stop the server with `Ctrl+C`.
You can also use `helium-studio --port 5005`. On Windows,
`./start_fluid_eda.ps1` creates the local environment when needed and starts the app.

For a first run **without COMSOL**, click **Example circuit**, select
**Two-phase pot (1D)** in the left sidebar, and click **Run simulation**.
Alternatively, import [continuous_pump_two_phase.json](examples/continuous_pump_two_phase.json).
The example button retains the COMSOL engine as its saved default; explicitly
select the Python engine when COMSOL has not been configured.

## Choose a solver

| UI mode | Engine | Scope |
| --- | --- | --- |
| Quick preview | Python / SciPy | Gas networks; reduced pipe model with equal gas/wall cross-sectional temperature. |
| COMSOL coupled simulation | Python network + COMSOL 6.1 pipe solves | Gas networks; axisymmetric laminar flow and conjugate gas/steel heat transfer. |
| Two-phase pot (1D) | Python / SciPy / CoolProp | One saturated pot → suction pipe → pump → precooler → capillary → pot, with forward circulation. |

**Two-phase pot (1D) does not call COMSOL.** Its precooler and capillary use
homogeneous equilibrium approximations. Its suction tube carries single-phase
vapor. Other topologies use the legacy gas-network engines. Selecting COMSOL
never silently falls back to a Python preview.

See [COMSOL setup](docs/comsol.md), [physics and assumptions](docs/physics.md), and
the [two-phase model](docs/two-phase-model.md).

## Continuous-pot example

```text
Saturated vapor pot → Suction pipe → Pump → Precooling pipe → J–T capillary
         ↑_______________________________________________________|
```

| Component | Supplied example settings |
| --- | --- |
| Pot | Saturated helium-4 at 1 K; applied heat load initially 0 mW. |
| Suction pipe | SS304; 500 mm long; 63.5 mm OD; 7.62 mm wall; exterior insulated; cold/warm boundaries 1 / 293.15 K. |
| Pump | Forward–backflow; nominal forward speed 240 L/s; zero-flow compression ratio 10⁸; discharge temperature 293.15 K. |
| Precooler | SS304; 63.5 mm OD; 7.62 mm wall; five 300 mm outer-wall segments: 293.15 K, 293.15→35 K, 35 K, 35→4 K, and 4 K. |
| J–T | CuNi capillary; 0.10 mm bore; 100 mm length; measured impedance 3.2 × 10¹⁶ m⁻³. |

In the two-phase mode, pipe fluid outlet temperatures are calculated from finite
heat transfer; a 4 K wall does **not** guarantee a 4 K fluid outlet. The HiPace
nominal pump coefficients are starting assumptions, not a validated operating map.
Successful convergence does not establish that the example refrigerates at 1 K.

Trial the pot temperature and inspect **Net cooling power**, **Applied heat load**,
**Heat imbalance**, and mass balance. A steady inventory and near-zero heat
imbalance are necessary equilibrium conditions. See the
[user guide](docs/user-guide.md#interpreting-pot-results) for signs and units.

## Development

```sh
python -m pip install -e ".[dev]"
python -m pytest
python -m build
```

The tests cover network conservation, pump laws, API errors/cancellation,
COMSOL adapter behavior through mocks, helium property checks, capillary choking,
and two-phase loop energy balance. They do not require a COMSOL license.
Real COMSOL execution is a separate integration check. If Node.js is installed,
`node --check models/fluid_eda/static/eda.js` checks frontend syntax.

## Repository layout

```text
models/
  fluid_eda/        Editor, API, network solvers, COMSOL adapter, static assets
  comsol_pipe/      COMSOL template generator and gas/steel property functions
  continuous_he4/   Saturated-helium thermodynamic properties
  one_k_pot/        Shared helium viscosity and saturation-property functions
  config.py        Runtime and COMSOL paths
examples/          Importable circuit definitions
tests/             Automated regression tests
docs/              Usage, equations, configuration, and source references
var/               Generated jobs, caches, and COMSOL files (ignored by Git)
```

The `models` package name is retained to preserve the existing solver imports.
Only the property modules needed by this application are included from the
earlier modeling packages.

## Runtime data and limitations

The server listens on loopback and accepts local requests. Drafts stay in browser
local storage; exported JSON files are the portable way to save/share a design.
Jobs, logs, caches, generated plots, and `.mph` files are written locally and
excluded from Git. Set `HELIUM_STUDIO_DATA_DIR` to choose another writable runtime
directory. Environment variables must be set before starting the process.

These are research/engineering models with explicit property and flow-regime
limits. They do not resolve superfluid film transport, phase slip, or all rarefied
gas effects. Two-phase heat-transfer/friction closures and pump maps require
experimental calibration. See [sources](docs/sources.md) and
[contributing](CONTRIBUTING.md) for the numerical verification scope.
