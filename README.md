# Solid-gas thermosolver

## What this repository answers

This repository calculates whether a reacting gas can reduce rutile TiO2 at equilibrium. It also distributes a measured oxygen-vacancy inventory among surface, first-subsurface, and bulk rutile sites.

[Open the browser solver](https://robin-yk.github.io/Solid-gas-thermosolver/)

## Main results

An equilibrated CO2/H2 = 1:1 feed does not reduce rutile from 500 to 1500 °C. The nearest reduced phase remains 78.44 to 27.30 kJ mol-O⁻¹ above stability over this range.

| Gas charge | Equilibrium result |
| --- | --- |
| CO2/H2 = 1:1 | Rutile remains stable from 500 to 1500 °C |
| Pure H2 | Reaches the TiO2/Ti10O19 boundary |
| Sealed N2 | Produces only trace reduction at the highest temperatures |

Detailed tables and stability margins are in [docs/results.md](docs/results.md).

## Two workspaces

### Gas-solid equilibrium

Active-set Gibbs minimization determines the stable gas composition and Ti-O phase assemblage at fixed temperature, pressure, volume, and elemental inventory. Condensed-phase stability is evaluated from KKT reduced costs. The model is a closed, finite-inventory equilibrium calculation.

### Defect statistical mechanics

Given a measured oxygen-deficiency inventory, a canonical lattice model distributes vacancies among rutile (110) surface, first-subsurface, and bulk sites. The calculation uses particle geometry, site energetics, site exclusion, and surface interactions. It reports equilibrium occupancy, not a diffusion profile or reduction kinetics.

See [docs/equilibrium-method.md](docs/equilibrium-method.md) and [docs/defect-model.md](docs/defect-model.md) for the full formulations.

## Run locally

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt

python3 scripts/export_equations.py
python3 scripts/build_ti_solver.py
python3 scripts/build_portal.py
python3 -m pytest tests/ -q
python3 -m http.server 8000 --directory docs
```

The Python package is in `solidgas/`. Browser sources are in `web/`, generated site files in `docs/`, build and oracle commands in `scripts/`, and committed inputs and reference outputs in `data/`.

## Verification

Production results are checked against an independent 80-digit thermodynamic oracle and a 50-digit statistical-mechanical oracle. The test suite also enforces Python-to-JavaScript parity, elemental balance, active-set KKT conditions, deterministic Monte Carlo trajectories, embedded-page freshness, and self-contained browser builds.

Numerical methods and measured error bounds are documented in [docs/verification.md](docs/verification.md).

## Limitations

- Equilibrium calculations contain no kinetic or transport information.
- The rutile TiO2-x solution phase is not implemented.
- CH4 calculations exclude graphite and should not be quoted.
- CO2-refill kinetic parameters are placeholders.
- The unreconstructed rutile (110) surface Hamiltonian is extrapolative beyond reconstruction onset.

Data scope and planned extensions are tracked in [docs/thermodynamic-data.md](docs/thermodynamic-data.md) and [docs/roadmap.md](docs/roadmap.md).
