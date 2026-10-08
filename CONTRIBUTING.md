# Contributing

Use Python 3.11+ and a local virtual environment. Install with
`python -m pip install -e ".[dev]"`, then run `python -m pytest`.
The optional `comsol` extra is needed only for real COMSOL model construction
and execution. Browser assets are plain HTML/CSS/JavaScript.

## Changes and verification

- Keep internal quantities in SI and identify absolute pressures explicitly.
- Preserve component IDs in errors and A/B identities when changing schematic
  rendering, rotation, or mirroring.
- Document assumptions separately from sourced data and fitted parameters.
- Physics changes need conservation checks and a meaningful limiting case or
  independent property reference. Numerical tests are not experimental validation.
- Test COMSOL adapter/protocol changes without a license through mocks, and
  state separately whether a real licensed solve was performed.
- Update the relevant documentation when equations, boundary conditions, source
  ranges, or topology support change.

The normal test suite includes API integration, analytical pump/network checks,
real Python two-phase solves, conservation, and selected refinement checks.
GitHub Actions runs it on Windows and Linux and builds the source distribution
and wheel. CI does not contain COMSOL or a license.

## Runtime files

Put runtime outputs in `var/` (or the configured data directory). Do not commit
job histories, solver logs, caches, generated `.mph` files, virtual environments,
local credentials, or private experimental circuits. Share a minimal circuit
definition and reproduction steps when reporting a bug.

The standalone checkout intentionally preserves the `models.*` import layout.
Only shared property functions needed by the EDA are included; importing unrelated
research packages would break a clean installation.

## Licensing

A project-wide redistribution license has not been selected. Public visibility
does not by itself grant an open-source license. Third-party packages and
referenced scientific data retain their respective licenses and terms.
