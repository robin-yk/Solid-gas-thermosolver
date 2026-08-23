# Defect statistical-mechanical model

## Question and ensemble

The defect workspace starts from a measured total oxygen-vacancy inventory. It asks how that fixed inventory partitions among rutile (110) bridging-oxygen sites, the first subsurface oxygen layer, and bulk oxygen sites.

This is a canonical equilibrium calculation. It does not predict how the inventory was created and does not interpret the particle graphic as a diffusion profile.

## Particle geometry

Particle diameter or BET area determines the real number of sites in each class. For a 0.9 µm spherical rutile particle, the outer bridging-oxygen layer is about 0.05% of all oxygen sites. The representative Monte Carlo slab contains a much larger surface fraction, so slab class ratios are never used as particle-scale capacities.

The rutile (110) 1x1 cell (c x a*sqrt(2) = 19.22 A^2) carries one bridging oxygen. Every oxygen slab of thickness d110 = a/sqrt(2) = 3.25 A below the bridging plane carries four oxygens per cell, which reproduces the bulk oxygen density exactly. The subsurface class is therefore not one atomic row: one declared layer is four bridging rows of capacity, 54.23 against 13.56 umol-O/g at 0.9 um, and the near-surface shell holds 1 + 4*n_layers oxygens per cell.

## Shell thickness

`subsurface_layers` (default 1) is a declared model choice, not a fitted quantity. The default follows the DFT set, which resolves a single first-subsurface formation energy and takes deeper layers as bulk. Because that energy (0.59 eV) lies far below the bulk value (1.31 eV), the shell/bulk split depends on where the line is drawn. At 600 C, 95 umol-O/g and 0.9 um:

| n_layers | shell depth | subsurface | theta_ss | bulk | mu_V (eV) | shell share |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 0.41 nm | 51.76 | 0.954 | 36.47 | 0.8188 | 61.6% |
| 2 | 0.73 nm | 82.65 | 0.762 | 5.57 | 0.6776 | 94.1% |
| 3 | 1.06 nm | 86.26 | 0.530 | 1.96 | 0.5991 | 97.9% |
| 4 | 1.38 nm | 87.06 | 0.401 | 1.16 | 0.5599 | 98.8% |

The surface class is fixed at the reconstructed line-phase deficiency (6.78 umol-O/g) throughout, so the ladder moves inventory between the subsurface shell and the bulk only. The reconstruction boundaries barely move (2.31 to 2.32 and 6.78 to 6.80 umol-O/g), but the CS ceiling does: it is a chemical-potential condition on the bulk, so a thicker cheap shell holds more inventory before the bulk reaches it (160 / 213 / 266 / 319 umol-O/g for one to four layers). Real formation energies relax smoothly towards the bulk value over the first few layers, so the physical answer lies between the first two rows: a reduced region of the order of 0.5 to 1 nm at 50 to 95% local depletion, with a correspondingly smaller bulk remainder than the default row reports. What the choice does not put at stake is the conclusion that the near-surface shell carries most of the inventory - that grows from 62% to 99% along the ladder. The full ladder is certified by the oracle (`shell_ladder` in `data/reference_statmech.json`) and is reproduced by the shell-thickness control in the workspace.

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
