# Solid-gas thermosolver

## What this repository answers

This repository calculates whether a reacting gas can reduce rutile TiO2 at equilibrium. It also distributes a measured oxygen-vacancy inventory among surface, first-subsurface, and bulk rutile sites.

[Open the browser solver](https://robin-yk.github.io/Solid-gas-thermosolver/)

## Main results

An equilibrated CO2/H2 = 1:1 feed does not reduce rutile from 500 to 1500 °C. The nearest reduced phase remains 70.40 to 27.30 kJ mol-O⁻¹ above stability over this range.

| Gas charge | Equilibrium result |
| --- | --- |
| CO2/H2 = 1:1 | Rutile remains stable from 500 to 1500 °C |
| Pure H2 | Reaches the TiO2/Ti10O19 boundary |
| Sealed N2 | Produces only trace reduction at the highest temperatures |

Detailed tables and stability margins are in [docs/results.md](docs/results.md).

## Three workspaces

### Gas-solid equilibrium

Active-set Gibbs minimization determines the stable gas composition and Ti-O phase assemblage at fixed temperature, pressure, volume, and elemental inventory. Condensed-phase stability is evaluated from KKT reduced costs. The model is a closed, finite-inventory equilibrium calculation.

The workspace at `ti_solver.html#thermo` draws that stability twice: one bar per excluded phase at the condition on the panel, and the smallest of those bars swept across temperature with the winning assemblage as a band underneath, which is where the Magnéli ladder becomes visible. Both use the reduced cost per formula unit — the basis the KKT test is stated on, and the only one that stays positive once the host, feed and temperature are all free to move.

### Defect statistical mechanics

Given a measured oxygen-deficiency inventory, a canonical lattice model computes its conditional-equilibrium partition among rutile (110) surface, first-subsurface, and bulk sites, using particle geometry, site energetics, and site exclusion. The (1x2) surface reconstruction and the shear-plane solubility ceiling enter as competing phases (lever rule at pinned chemical potential), so nothing is clipped and the phase boundaries are certified numbers. A coarse-grained KMC layer adds reoxidation transport. It reports thermodynamic tendency, not a diffusion profile or reduction kinetics.

### Steady-state redox

A reductant that needs lattice oxygen and an oxidant that needs vacancies, written against one common vacancy fraction: one falls and one rises, so they cross once and the crossing is the steady state. Both barriers move with the same descriptor in opposite directions, so sweeping it gives a volcano whose peak has a closed form, and a family of solutions is available before any material is fitted. Feeding the layer-resolved partition in as the average-to-surface mapping leaves the steady rate untouched, moves the steady degree of reduction by three orders of magnitude, and removes the steady state entirely on the reducible half of the axis.

The model runs in the browser at `ti_solver.html#redox`, alongside a two-reservoir cycle model that takes measured per-cycle integrals and returns the recovered fraction, how much of each residual locked, and the prediction for the experiment that separates the two: two samples at the same total V_O reached by different paths.

See [docs/equilibrium-method.md](docs/equilibrium-method.md), [docs/defect-model.md](docs/defect-model.md) and [docs/redox-steady-state.md](docs/redox-steady-state.md) for the full formulations.

## Run locally

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt

python3 scripts/export_equations.py
python3 scripts/build_ti_solver.py
python3 scripts/build_portal.py
python3 scripts/export_plates.py       # freeze the manuscript figure data
python3 scripts/build_plates.py        # docs/figures.html proof sheet
python3 -m pytest tests/ -q
python3 -m http.server 8000 --directory docs
```

The Python package is in `solidgas/`. Browser sources are in `web/`, generated site files in `docs/`, build and oracle commands in `scripts/`, and committed inputs and reference outputs in `data/`.

Manuscript figures are drawn to print size from frozen solver output: `scripts/export_plates.py` extracts every value the plates carry into `data/figure_data.json` and stamps it with the commit, and `scripts/build_plates.py` renders the proof sheet at [docs/figures.html](docs/figures.html). Captions are written with slots rather than numbers, so a caption cannot drift from the model, and the test suite fails while the extract is stale.

## Verification

Production results are checked against an independent 80-digit thermodynamic oracle and a 50-digit statistical-mechanical oracle. The test suite also enforces Python-to-JavaScript parity, elemental balance, active-set KKT conditions, deterministic Monte Carlo trajectories, embedded-page freshness, and self-contained browser builds.

Numerical methods and measured error bounds are documented in [docs/verification.md](docs/verification.md).

## Limitations

- Equilibrium calculations contain no kinetic or transport information.
- The rutile TiO2-x solution phase is not implemented.
- CH4 calculations exclude graphite and should not be quoted.
- CO2-refill kinetic parameters are placeholders.
- The unreconstructed rutile (110) surface Hamiltonian is extrapolative beyond reconstruction onset.
- The subsurface class thickness is a declared choice (one oxygen layer by default, four O per surface cell); the shell/bulk split moves with it and the workspace exposes the ladder.

Data scope and planned extensions are tracked in [docs/thermodynamic-data.md](docs/thermodynamic-data.md) and [docs/roadmap.md](docs/roadmap.md).
