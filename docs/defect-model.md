# Defect statistical-mechanical model

## Question and ensemble

The defect workspace starts from a measured total oxygen-vacancy inventory. It asks how that fixed inventory partitions among rutile (110) bridging-oxygen sites, the first subsurface oxygen layer, and bulk oxygen sites.

This is a canonical equilibrium calculation. It does not predict how the inventory was created and does not interpret the particle graphic as a diffusion profile.

## Particle geometry

Particle diameter or BET area determines the real number of sites in each class. For a 0.9 µm spherical rutile particle, the outer bridging-oxygen layer is about 0.05% of all oxygen sites. The representative Monte Carlo slab contains a much larger surface fraction, so slab class ratios are never used as particle-scale capacities.

## Site thermodynamics

Surface vacancies interact through short-range pair terms along and across bridging-oxygen rows. Canonical swap Metropolis sampling generates the surface isotherm, and Widom insertion supplies the vacancy chemical potential. First-subsurface and bulk vacancies use exact site-exclusion isotherms in the current version.

The particle-scale calculation finds the common vacancy chemical potential that satisfies

```text
N_surface theta_surface(mu_V)
+ N_subsurface theta_subsurface(mu_V)
+ N_bulk theta_bulk(mu_V)
= N_vacancy,total
```

The site energies, interaction terms, geometry constants, and model defaults are stored in `data/rutile_dft.json`.

## Surface reconstruction

The reported 17 to 20% bridging-site coverage is an onset interval for rutile (110) reconstruction, not a hard vacancy capacity. The solver does not clip surface coverage at 20%. At or above this interval it marks the result as a conditional extrapolation of the unreconstructed surface Hamiltonian because reconstructed-phase energetics are not implemented.

## Particle view

The spherical cross-section uses concentric site classes. Surface positions come from the sampled surface configuration, the first-subsurface shell reflects its calculated class occupancy, and the core fill represents the bulk-class occupancy. It is a compact equilibrium occupancy map, not a random placement of vacancies through a diffusion field.

## Thermodynamic bridge

The oxygen-environment card re-solves the gas and solid charge at the defect temperature. Rutile alone is classified as stable, rutile plus a reduced phase as a reduction boundary, and an assemblage without rutile as outside the defect model's validity.

## Placeholder layers

The CO2-accessibility calculation applies class-specific first-order refill probabilities. Its kinetic parameters are placeholders. The dielectric card evaluates the stored empirical fit over its calibration range and labels extrapolation outside that range.

## Verification

`scripts/oracle_statmech.py` provides an independent 50-digit reference for analytic matching, exact finite-lattice canonical ensembles, and enumerated interacting lattices. The JavaScript engine follows the Python production engine move for move under the same seed.
