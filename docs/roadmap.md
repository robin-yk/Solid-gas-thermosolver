# Roadmap

## Thermodynamic model

1. Add a rutile TiO2-x solution phase with assessed vacancy thermodynamics.
2. Add graphite before releasing any CH4 equilibrium result.
3. Extend the gas and condensed registries only with source-traceable data and oracle coverage.

## Defect model

1. Replace literature placeholder energetics with project DFT values through `data/rutile_dft.json`.
2. Add reconstructed rutile (110) energetics rather than extrapolating the unreconstructed surface Hamiltonian.
3. Add explicit Ti3+ polaron variables where required by the data.

## Kinetics and transport

1. Replace the CO2-refill placeholders with fitted class-specific parameters.
2. Add radial coarse cells and Arrhenius exchange rates in detailed balance with the equilibrium model.
3. Verify that the long-time kinetic limit reproduces the canonical equilibrium distribution.

These additions must preserve the existing separation between model data, production solvers, independent oracles, browser ports, and generated pages.
