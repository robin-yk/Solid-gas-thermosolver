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
| 1 | 0.41 nm | 3.59 | 0.066 | 84.63 | 0.2372 | 10.9% |
| 2 | 0.73 nm | 6.96 | 0.064 | 81.26 | 0.2284 | 14.5% |
| 3 | 1.06 nm | 10.13 | 0.062 | 78.09 | 0.2198 | 17.8% |
| 4 | 1.38 nm | 13.11 | 0.060 | 75.12 | 0.2115 | 20.9% |

The surface class is fixed at the reconstructed line-phase deficiency (6.78 umol-O/g) throughout, so the ladder moves inventory between the subsurface shell and the bulk only. The reconstruction boundaries move a little (32.75 to 36.81 and 37.23 to 41.28 umol-O/g), and the CS ceiling moves with them (110.8 / 114.6 / 118.5 / 122.4 umol-O/g for one to four layers), because it is a chemical-potential condition on the bulk and a thicker cheap shell holds more inventory before the bulk reaches it.

The reading of this table reversed when the compensating polarons went into the free energy. It used to say that a single shell held 62% of the inventory and four shells held 99%; it now says 11% and 21%. The surface-to-bulk enrichment is the cube root of what the neutral-vacancy form gave - exp(dE/3kT) = 333 rather than exp(dE/kT) = 3.7e7 - and against 1846 times as many bulk sites that is not enough to win. **Most of the oxygen is in the bulk, at every loading.** The conclusion the old table supported, that the near-surface shell carries most of the inventory, was an artifact of a free energy that did not pay for the vacancy's electrons.

What survives is that the shell is enriched: at 0.9 um the first oxygen layer runs twenty times the bulk plateau, and the outer nanometre holds twelve percent of the oxygen in seven tenths of a percent of the volume. That is the depth-resolved engine's statement rather than this one's, and it is drawn in [particle-model.md](particle-model.md).

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

## What was removed, and what replaced it

Two degrees of freedom used to sit here and both are gone.

**Layer resolution.** The near-surface shell could be split into up to three classes, with their energies taken from a `span_fraction` between the subsurface and bulk values. That fraction was a declared shape, labelled *assumed* in the dataset, and it is superseded: `solidgas/particle/` derives the depth profile

```
E(z) = E_bulk - dE_seg exp(-z/xi),   xi = -d110 / ln[(E_b - E_ss)/(E_b - E_s)]
```

from the same three DFT energies, which fixes both parameters and leaves no shape to assume. It lands within a few hundredths of the fractions that were being declared, which is a check on both.

**The crowding term.** Each class could carry an effective mean-field penalty `omega_i`. It had no literature value, and its own provenance note said the magnitudes were chosen to bring the first subsurface layer out of near-saturation. That near-saturation was the symptom of a missing term, not of crowding: the free energy was written on the oxygen sublattice alone, so it had no cap and put the first subsurface layer at theta = 0.954. With the cation sublattice counted the same layer sits at 0.065 and there is nothing for `omega` to fix.

The evidence that both were standing in for the missing entropy rather than adding physics: with the compensation in place, resolving the shell into three layers moves the bulk inventory by 2%, where before it moved it by 30%.

Both blocks are removed from `data/rutile_dft.json`, which carries a `retired` note in their place, and both controls are gone from the workspace. What sits there now is the depth-resolved card described in [particle-model.md](particle-model.md).

## What this engine is still for

The three-class partition is the coarse reading, and it is the one the site-resolved Monte Carlo, the CO2 accessibility model, the dielectric correlation and the coarse-grained KMC are all written against. Its free energy is not its own: `solidgas/statmech.py` re-exports `solidgas/particle/freeenergy.py`, so there is one compensated relation in the repository and the coarse partition cannot drift from the depth-resolved one. The two agree to about 1.5% at the shipped condition, the difference being that the depth-resolved engine gives the layers below the first their own energies instead of the bulk value.

The shell thickness is still a declared choice, and still matters, because it moves *sites*: every slab in the shell is charged the first-subsurface energy, so making the shell thicker makes it cheaper than it should be. That is the reason more than one slab is coarse, and the reason the depth-resolved engine exists.

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

The card is in the browser, on the redox workspace (`ti_solver.html#redox`, *Cycling*): the per-cycle integrals are editable, `f` and `lock` follow, and the two-sample gap is quoted for any total and any number of build cycles. `web/cycling.js` is a line-for-line port and is held to equality with the Python module by `tests/test_redox_port.py`, including a five-cycle round trip that generates a schedule at (0.7, 0.35) and fits it back.

**What still cannot be tested this way.** The dielectric correlation is fitted over roughly 95-5000 umol-O/g and the residuals of interest are 13 and 23, four to seven times below its calibration floor, so a residual-memory test against it is an extrapolation until low-inventory points exist. And because that correlation reads total V_O with no depth resolution, it cannot separate "deep vacancies remain" from "vacancies refilled but the structure did not recover" - both give the same loss at the same inventory.

## Steady-state redox

`solidgas/redox.py` is a separate, smaller model that uses this one as a closure. It writes a reduction and an oxidation half-reaction against one common vacancy fraction and takes the steady state as their crossing. See `docs/redox-steady-state.md`; the connection is that the rate laws need the vacancy fraction on the *reacting face*, and the partition above is what supplies it.

## Verification

`scripts/oracle_statmech.py` provides an independent 50-digit reference for the phase-aware analytic matching (including the closed-form boundary ladder), exact finite-lattice canonical ensembles, enumerated interacting lattices, and a high-precision matrix exponential of the dilute transport master equation. The JavaScript engines follow the Python production engines move for move under the same seed.
