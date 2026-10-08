# COMSOL setup

The COMSOL engine solves each pipe using a parameterized axisymmetric helium/SS304
model and couples the resulting mass flow back into the Python pressure network.
The template remains a single-phase gas model. The two-phase engine is separate.

## Requirements

- COMSOL Multiphysics **6.1**, with licenses for the Laminar Flow,
  Heat Transfer in Solids and Fluids, and Nonisothermal Flow interfaces used by
  the template.
- Python dependencies from the `comsol` extra. COMSOL itself and its license are
  not included or installed by pip.
- MPh must be able to discover the COMSOL 6.1 installation.

From the activated environment:

```sh
python -m pip install -e ".[comsol]"
```

Set the executable path for your installation. PowerShell example:

```powershell
$env:COMSOL_EXECUTABLE = 'E:\COMSOL61\COMSOL61\Multiphysics\bin\win64\comsol.exe'
python -m models.comsol_pipe.build_model --build-only
python -m models.fluid_eda.app
```

On a compatible Linux installation:

```sh
export COMSOL_EXECUTABLE=/usr/local/comsol61/multiphysics/bin/comsol
python -m models.comsol_pipe.build_model --build-only
python -m models.fluid_eda.app
```

The builder generates `var/comsol/helium_ss304_pipe.mph` with geometry, mesh,
materials, boundary selections, and solver configuration. `--build-only` avoids
solving the standalone demonstration boundary conditions. The circuit adapter
sets each pipe's geometry and boundary states before solving.

No precomputed `.mph` solution is required or distributed. The default builder
mesh is 640 axial / 32 radial gas elements; the EDA adapter applies the circuit
mesh settings and refines around wall-profile transitions.

## Configuration

| Variable | Meaning | Default |
| --- | --- | --- |
| `COMSOL_EXECUTABLE` | Full path to the COMSOL executable | `comsol` on PATH, then a conventional platform installation path. |
| `HELIUM_STUDIO_DATA_DIR` | Writable jobs/cache/generated-model directory | `var/` beside the source checkout. |
| `HELIUM_STUDIO_COMSOL_TEMPLATE` | Full path to the generated template | `<data directory>/comsol/helium_ss304_pipe.mph`. |

Set these in the shell before both building and launching. The same template
setting is shared by the builder and adapter. A nonstandard installation may also
need configuration supported by [MPh](https://mph.readthedocs.io/).
The application does not load a `.env` file automatically.

## Execution and diagnostics

A separate worker owns the COMSOL Java client. The web server remains responsive
and can cancel the worker. Template, adapter, and material-function fingerprints
identify cached pipe evaluations. Regenerating the template invalidates its cache.
Runtime artifacts and logs remain under the local data directory.

The status endpoint reports whether the executable and template files exist;
this is not a license checkout or solver-convergence test. Component-specific
errors identify the pipe and boundary values that failed. COMSOL errors are not
replaced by a successful Python preview.

The adapter requires finite positive absolute pressures and temperatures, valid
geometry, and mass/energy conservation. It does not impose a blanket 1–300 K
acceptance window. Removing that window does not extend the validity of the
individual material fits or guarantee convergence.
