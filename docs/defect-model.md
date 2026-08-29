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

## Layer resolution and the crowding term

Two degrees of freedom were added to the site thermodynamics above. Both ship switched off, so the partition the workspace reports is the one already certified; `resolved_layers` and `coverage_penalty` in `data/rutile_dft.json` turn them on.

**Layer resolution.** The near-surface shell can be split into up to three classes, each one d110 oxygen slab. Their energies come from `layer_profile.span_fraction`, which is a *shape* rather than a table: `eps_i = eps_subsurface + span_i * (eps_bulk - eps_subsurface)`, so the first layer is the preset's own subsurface energy exactly and the profile carries over to a different functional without being re-tabulated. At one resolved layer the class is the declared shell and the legacy capacity is reproduced exactly.

What this buys is a profile the one-class model cannot represent. At 600 C, 95 umol-O/g and 0.9 um the single class reports 95% occupancy for the whole shell; resolved into three it reports 92% / 28% / 3%, and mu_V rises from 0.599 to 0.779 eV because the one-class model had been charging all three slabs the cheap first-subsurface energy. The bulk remainder moves accordingly. The conclusion that the shell carries most of the inventory is unaffected; what changes is the depth over which it is spread.

**The crowding term.** Each class may carry an effective mean-field penalty, so its free energy per site is `eps_i theta_i + (omega_i/2) theta_i^2 + kT s(theta_i)` and its chemical potential is

```text
mu = eps_i + omega_i theta_i + kT ln(theta_i/(1-theta_i))
```

which is implicit in theta and has no closed form. It is solved in the logit, where the sigmoid bounds the root between `(mu-eps-omega)/kT` and `(mu-eps)/kT` so no bracket has to be searched for, and `d(mu)/d(theta) = omega + kT/(theta(1-theta))` is strictly positive for `omega >= 0` so the root is unique. At `omega = 0` the bracket collapses to a point and the site-exclusion isotherm is returned directly rather than bisected to, which is what makes the legacy numbers reproduce bit for bit rather than to a tolerance.

`omega_i` is an *effective* quantity, `omega_i ~ sum_j z_ij J_ij`, and not a vacancy pair energy: the surface ordering pair table must not be summed into it, because coordination and pair counting would be double counted. It has no literature value and is carried as a sampling range, not a measurement. Only repulsive values are admitted — `omega < -4kT` makes `mu(theta)` non-monotone and splits a class into two phases, a construction this branch does not carry, so a negative value is refused rather than allowed to produce a silently wrong root.

**What was deliberately not touched.** `omega_surface` is pinned at zero in the dataset, so the (1x2) phase ladder cannot move and any change in a deeper class has one possible cause. The surface Monte Carlo is untouched for the same reason. The bulk interaction *does* enter the CS coexistence, `mu_cs = eps_bulk + omega_bulk theta_sol + kT ln(theta_sol/(1-theta_sol))`, because that boundary is a condition on the bulk class: leaving it out would put the ceiling at the non-interacting potential while the bulk sitting at it was interacting, and the two would disagree about what `x_sol` means.

Certified by `theta_of_mu` and `layer_ladder` in `data/reference_statmech.json` — the mpmath oracle re-derives the implicit root with a different algorithm and the whole construction at every resolution, across all four phase branches — and by three-way parity against the browser engine in `tests/test_statmech_port.py`.

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

## What this model does not carry

Three limits that matter for reading any number above, stated here rather than left to be inferred.

**One transport channel.** The CGMC moves oxygen vacancies and nothing else. Titanium interstitials are a real and, on reduced rutile, sometimes dominant mass-transport route (Henderson, Surf. Sci. 419, 174 (1999)): under reducing conditions Ti moves out through the lattice and regrows surface layers, which is a different process from a vacancy hopping inward. Nothing in this workspace represents it, so a fitted effective barrier absorbs whatever share of the real transport runs through that channel.

**The dielectric correlation has no depth.** The stored power law relates the measured loss to the *total* vacancy inventory. It was fitted against a single number per sample and therefore cannot distinguish a surface-concentrated inventory from a uniform one. The partition on this page says those two states are very different; the correlation cannot see the difference, so it must not be read as a probe of where the vacancies are.

**Polarons are not variables.** The site energies are polaron-relaxed values taken from DFT, which is not the same as tracking the Ti(3+) polarons themselves. The vacancy is most stable in the surface layer while the electrons it donates prefer subsurface Ti (Reticcioli et al., Phys. Rev. B 98, 045306 (2018)), and that separation is baked into the numbers rather than represented. In particular the fitted `E_m,eff` from the CGMC panel is a vacancy migration barrier absorbing everything not modelled; it is not a polaron hopping rate and should not be compared with one.

## Cycling, and what one scalar cannot carry

The partition above has a single scalar state, the total inventory, and it re-equilibrates that total on every call. That is the right model for one measurement at one condition and the wrong one for a sequence of them: it cannot tell a vacancy that has just been made from one that already refused to come out.

The reported two-cycle experiment shows the difference. Against the inventory present at the start of each reoxidation, recovery falls from 0.860 to 0.694. A memoryless partition cannot produce that: with a smaller total it re-equilibrates into a distribution with a *larger* accessible share, so it can only predict the same fraction or a higher one. The sign of the change is the evidence, and the fix is a second reservoir.

`solidgas/cycling.py` carries a pair — accessible `A` and locked `L` — with one parameter joining them. Reduction adds to `A`; reoxidation removes `f*A` and locks a fraction `lock` of what is left. At `lock = 0` nothing ever locks and the model is the memoryless one exactly, so the memory is a declared addition rather than a new assumption.

Fitted to the two reported cycles it gives `f = 0.860` and `lock` pinned at its bound of 1: not only did none of the 13.3 umol-O/g residual return, complete locking is the least recovery the construction can predict and it still over-shoots the second cycle by 0.84 umol-O/g, 1.6%. So the fresh inventory itself recovered slightly less well the second time. That also settles a question of definition: normalising recovery against the *newly created* oxygen rather than the total present is the correct reading here, and it is correct because the residual is inert, not by assumption.

`lock` is phenomenological. It is equally consistent with vacancies migrating beyond reach, with aggregation into an ordered defect, and with a local relaxation this model does not resolve; cycle integrals alone cannot choose. An unlocking rate is the obvious third parameter and is left out rather than fitted to data that cannot see it.

**The experiment that measures it.** Two samples brought to the same total inventory by different paths — one long reduction against several reduce/reoxidise cycles — then reoxidised and compared. A memoryless inventory predicts a gap of exactly zero at any `f`, because total V_O is its whole state. At `f = 0.860` and complete locking the model predicts the cycled sample recovers 7.6 umol-O/g less out of 95, a 9.4% gap, and the gap is monotone in `lock`, so the measurement reads the parameter directly. `two_sample_prediction()` produces the number for any pair of values.

**What still cannot be tested this way.** The dielectric correlation is fitted over roughly 95-5000 umol-O/g and the residuals of interest are 13 and 23, four to seven times below its calibration floor, so a residual-memory test against it is an extrapolation until low-inventory points exist. And because that correlation reads total V_O with no depth resolution, it cannot separate "deep vacancies remain" from "vacancies refilled but the structure did not recover" - both give the same loss at the same inventory.

## Steady-state redox

`solidgas/redox.py` is a separate, smaller model that uses this one as a closure. It writes a reduction and an oxidation half-reaction against one common vacancy fraction and takes the steady state as their crossing. See `docs/redox-steady-state.md`; the connection is that the rate laws need the vacancy fraction on the *reacting face*, and the partition above is what supplies it.

## Verification

`scripts/oracle_statmech.py` provides an independent 50-digit reference for the phase-aware analytic matching (including the closed-form boundary ladder), exact finite-lattice canonical ensembles, enumerated interacting lattices, and a high-precision matrix exponential of the dilute transport master equation. The JavaScript engines follow the Python production engines move for move under the same seed.
