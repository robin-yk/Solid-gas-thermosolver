# Defect statistical-mechanical model

## Question and ensemble

The defect workspace starts from a measured total oxygen-vacancy inventory. It asks how that fixed inventory partitions among rutile (110) bridging-oxygen sites, the first subsurface oxygen layer, and bulk oxygen sites.

This is a canonical equilibrium calculation. It does not predict how the inventory was created and does not interpret the particle graphic as a diffusion profile.

## Particle geometry

Particle diameter or BET area determines the real number of sites in each class. For a 0.9 µm spherical rutile particle, the outer bridging-oxygen layer is about 0.05% of all oxygen sites. The representative Monte Carlo slab contains a much larger surface fraction, so slab class ratios are never used as particle-scale capacities.

## Site thermodynamics

The layer partition is analytic and uses the polaron-relaxed formation energies alone: every class follows its exact site-exclusion isotherm, and the particle-scale calculation finds the common vacancy chemical potential that satisfies

```text
N_surface theta_surface^phase(mu_V)
+ N_subsurface theta_subsurface(mu_V)
+ N_bulk theta_bulk(mu_V)
+ N_extended
= N_vacancy,total
```

Surface lateral ordering is a deliberately separate calculation: effective pair interactions along and across the bridging-oxygen rows, constrained to reproduce the vacancy-correlation statistics reported by Birschitzky et al., are sampled by canonical swap Metropolis with Widom insertion. They inform the arrangement exhibits only and never enter the layer partition — they parameterize arrangement statistics, not layer energetics.

The site energies, phase parameters, interaction terms, geometry constants, and model defaults are stored in `data/rutile_dft.json`.

## Phase construction

Two departures from dilute point-defect behavior are treated as thermodynamic phases, not as validity cutoffs or caps.

**Surface reconstruction.** The reported 17 to 20% coverage interval is the onset of the rutile (110)-(1x2) reconstruction. The solver treats the reconstruction as a competing surface phase: at theta_t = 0.17 the (1x1) vacancy lattice gas undergoes a first-order transition to the added-row Ti2O3 line phase (Onishi-Iwasawa stoichiometry: one added Ti2O3 unit per (1x2) cell is one oxygen vacancy per two bridging sites, an areal deficiency of 0.5 bridging-monolayers), handled by the lever rule at pinned mu_V. The high-coverage branch of the ideal (1x1) gas is excluded from the phase set on physical grounds: the repulsive ordering energies penalize it and the observed strongly reduced surface is reconstructed, not bridging-stripped. Declared idealizations: the line phase has zero width, cross-linked (1x2) variants are not distinguished, and no surface phase beyond theta_eff = 0.5 is modeled.

**Shear-plane ceiling.** Beyond the approximate rutile point-defect solubility (x_sol ~ 0.008 in TiO2-x, an estimate carried as a sensitivity parameter), mu_V pins at the estimated rutile/CS-phase coexistence and the excess inventory precipitates as extended defects by the same lever rule, reported as its own output class.

At 600 C and the default 0.9 um geometry the certified boundary ladder is: (1x2) onset at 2.31, reconstruction complete at 6.78, CS ceiling at ~160 umol-O/g. The flagship 95 umol-O/g inventory therefore sits in the reconstructed regime: surface (1x2) phase 6.78 (7.1%), first-subsurface layer 51.8 (54.5%, at 95% layer occupancy — a heavily reduced shell, not isolated point defects), bulk 36.5 (38.4%), just below the CS ceiling. The measured BET area (1.60 m2/g) and the diameter model agree on the ladder to about 2%.

## Particle view

The spherical cross-section uses concentric site classes. Surface positions come from the sampled surface configuration, the first-subsurface shell reflects its calculated class occupancy, and the core fill represents the bulk-class occupancy. It is a compact equilibrium occupancy map, not a random placement of vacancies through a diffusion field.

## Thermodynamic bridge

The oxygen-environment card re-solves the gas and solid charge at the defect temperature. Rutile alone is classified as stable, rutile plus a reduced phase as a reduction boundary, and an assemblage without rutile as outside the defect model's validity.

## Kinetic layers

The CO2-accessibility calculation applies class-specific first-order refill probabilities. Its kinetic parameters are placeholders, and every derived recoverable fraction is labeled illustrative until fitted to reoxidation measurements. The dielectric card evaluates the stored empirical fit over its calibration range and labels extrapolation outside that range.

A coarse-grained kinetic Monte Carlo transport layer (`solidgas/cgmc.py`, `web/cgmc.js`) follows Katsoulakis & Vlachos, J. Chem. Phys. 119, 9412 (2003): radial coarse cells over the particle, exchange rates in detailed balance with the (1x1) site-class Hamiltonian (the closed system relaxes to the exact product-binomial equilibrium), the dilute limit reducing to spherical diffusion, and the CO2 surface refill as the one irreversible channel. With the neutral-vacancy bulk barrier (0.65 eV, Iddir et al.) transport is not rate-limiting at 600 C, so a measured partial recovery reads back as an effective migration barrier, which the panel can fit. Reconstruction kinetics are not modeled.

## Verification

`scripts/oracle_statmech.py` provides an independent 50-digit reference for the phase-aware analytic matching (including the closed-form boundary ladder), exact finite-lattice canonical ensembles, enumerated interacting lattices, and a high-precision matrix exponential of the dilute transport master equation. The JavaScript engines follow the Python production engines move for move under the same seed.
