# Steady-state redox

## Question

An oxygen carrier sits in a gas that both reduces and oxidises it. A reductant pulls lattice oxygen out and leaves a vacancy; an oxidant fills a vacancy and puts oxygen back. Where does the solid settle, and how fast does it turn over there?

This is a **Mars–van Krevelen** mechanism, and the model is the Mars–van Krevelen site balance at steady state: nothing adsorbs, there are no surface intermediates, and the two steps are written straight onto the lattice-oxygen sites. The model is one nonlinear equation. Nothing in it is new — it is 1970s reaction engineering — and the point is not the method but that the question has not been carried through for this class of material with the solid described honestly.

## What theta counts

Vacancies on the **two-coordinate bridging oxygen rows of rutile (110)**, as a fraction of those sites.

That is not a free choice; it is the site the rest of the model is already referenced to. The descriptor is `E_vac` for (110) bridging O — 4.39 eV, sX hybrid, Li, Guo & Robertson, *J. Phys. Chem. C* **119** (2015) — and the face ceiling `theta_max = 0.5` is the (1x2) added-row reconstruction, which is a bridging-row phase. Both numbers are bridging-O numbers, so `theta` has to be a bridging-O coverage or the two do not belong in the same equation.

Three-coordinate in-plane oxygen is **not** counted. It is bound more strongly, its removal energy is not the descriptor's, and averaging the two into one `theta` would produce a number describing neither site. The cost of keeping them apart is that the reacting inventory here is the bridging rows alone: a mechanism that consumes in-plane oxygen is outside this model, and so is any face other than (110).

The particle-average `theta` is a different kind of quantity. It is stoichiometric bookkeeping, `x/2` in `MO(2-x)`, and carries no claim about which defect holds the reduction — in rutile at low `x` the majority point defect is the titanium interstitial, not the oxygen vacancy, and past the solubility the solid stops making point defects at all and shears instead. The face `theta` is a coverage on a named site; the particle `theta` is a composition. Keeping them separate is the whole content of the enrichment factor below.

## The two half-reactions

Per reacting oxygen site,

```text
r_red(theta) = A_r p_r^m_r (1 - theta)^n_r exp(-E_r/kT)     falls with theta
r_ox (theta) = A_o p_o^m_o      theta ^n_o exp(-E_o/kT)     rises with theta
```

Reduction needs oxygen to remove, so it slows as the surface empties. Oxidation needs vacancies to fill, so it speeds up. One falls and one rises over `theta` in `[0, 1]`, so they cross exactly once and the crossing is the steady state. The extent of reduction enters the rate as a multiplicative factor, not as a shifted baseline.

`steady_state()` bisects for the crossing. It bisects in the logit rather than in `theta`, because `1 - theta` cannot be resolved by subtraction once `theta` approaches one and the rate residual degrades to about 1e-5 at the ends of the descriptor axis; in the logit both tails keep full relative precision and the residual stays at 1e-12 everywhere. At unit site orders the answer is available in closed form,

```text
theta* = k_red / (k_red + k_ox)
```

which is kept public as `theta_first_order()` so the shipped numeric path can be checked against it rather than replaced by it.

## Why it is a volcano and not a ranking

Both barriers move with the same material descriptor — the oxygen vacancy formation energy — and they move in opposite directions. An oxide that gives up oxygen easily is quick to reduce and slow to reoxidise. With linear (Bronsted-Evans-Polanyi) coupling,

```text
E_red = E_red0 + alpha * dev
E_ox  = E_ox0  - beta  * dev
```

the steady rate is the harmonic mean of two opposed exponentials, and a harmonic mean of opposed exponentials peaks. Maximising it gives `beta k_red = alpha k_ox`, hence

```text
theta at the peak = alpha / (alpha + beta)
dev at the peak   = [ (E_ox0 - E_red0) - kT ln(alpha/beta)
                      + kT ln(C_r / C_o) ] / (alpha + beta)
```

The peak *coverage* carries no temperature, no pressure and no prefactor: only the two slopes set it. That is what makes a family of solutions available before any material is fitted — sweep the descriptor and every row is a hypothetical oxide placed by one number.

## Parameters

Everything lives in `data/redox.json` with a label. Nothing is fitted.

The descriptor is an **offset** from a reference, not an absolute energy, because only differences enter the two barriers. The reference is quoted for orientation as rutile (110) bridging O at 4.39 eV (Li, Guo & Robertson, J. Phys. Chem. C 119 (2015)) and is strongly functional dependent: the same authors report 3.71 eV with GGA for the same site. That 0.68 eV spread is the systematic uncertainty on placing a material along the axis, and it is comparable to the width of the volcano — so a computed formation energy alone does not rank two materials that sit within about 0.7 eV of each other.

## The axis with nothing assumed in it

Everything above is written in terms of a material descriptor and two linear free-energy slopes, and each of those is a parameter that has to be assumed. None of them is needed to state the steady state. All the material dependence collapses into one dimensionless number,

```text
K = k_reduction / k_oxidation        theta* = K / (1 + K)
```

`from_ratio()` takes K and nothing else — no temperature, no pressure, no prefactor, no descriptor, no slopes — and returns the same crossing the descriptor route finds when handed the same ratio. The Bronsted-Evans-Polanyi layer is best read as an *overlay* on this axis: it is what turns K into a material ranking, and it is where the assumptions live.

**The design window.** A carrier can be run while two things hold: the solid stays one phase, and the reacting face stays below the coverage its own phase pins it at. Both are closed forms in K, and both depend on the surface enrichment `E = theta_face / theta_particle`:

```text
second phase      K > E*theta_sol / (1 - E*theta_sol)
no steady state   K > theta_face,max / (1 - theta_face,max)
```

The uniform assumption is `E = 1`. For rutile at 600 C that gives a usable window of `K < 0.0040` — reduction must be 250 times slower than oxidation or the solid shears. The depth-weighted enrichment is `E = 269`, and `theta_sol * E = 1.08 > 1`, so the bulk solubility stops binding at all: the face saturates first and the window becomes `K < 1`, oxidation merely having to keep up with reduction. **The uniform assumption discards a factor of 249 of the design space**, and it does so silently, because nothing in it announces that a surface concentration was being used as a bulk composition.

That conclusion survived a change in `E` of nearly an order of magnitude and it survived it narrowly. `E` was 1844 while the vacancy free energy was written on the oxygen sublattice alone; counting the cation sublattice makes the dilute isotherm go as `exp(mu/3kT)` instead of `exp(mu/kT)`, so the enrichment is the cube root of what it was. The margin on `theta_sol * E > 1` went from 7.4 to 1.08. It is still on the right side of one, but it is no longer comfortable, and it is now sensitive to the solubility estimate `x_sol = 0.008`, which the dataset labels approximate.

The failure mode changes with it. Below `E = theta_face,max / theta_sol = 125` the solid shears before the face saturates; above it the face saturates first and the bulk never approaches the solubility. Which of the two binds decides what a carrier should be engineered against, and the uniform reading gets that wrong as well as the number.

One consequence is worth stating rather than leaving implicit: with the face pinned at 0.5, an *oxidation-limited* steady state — the solid sitting reduced, theta above 0.75 — is unreachable for this carrier at any rate constants, because that state lies past the boundary where no steady state exists at all.

The two figures are `docs/figures/redox_crossing.svg` and `docs/figures/redox_phase.svg`, drawn from rows the engine computes; the gates check that every boundary on the page is the closed form the engine returns.

## Domain

`theta` is a point-defect concentration only while the reduced solid stays one phase. Past the declared solubility the vacancies order into crystallographic shear planes and mass action in `theta` stops meaning anything. `steady_state()` reports that as a flag and a warning; it never clips.

This is where the uniform model runs into itself. With symmetric slopes the optimum sits at `theta = 1/2`. Half the oxygen sublattice vacant is a coverage a surface reaches routinely — the rutile (1x2) added-row phase is exactly that — and a bulk composition no oxide holds as a solid solution. In the uniform model those are the same number, so the model puts its own optimum outside the domain of every carrier. That is a defect of forcing one number to be both quantities, not a finding about oxides.

## What the depth-weighted correction actually does

`surface_correction()` takes a mapping from the particle-average vacancy fraction to the fraction on the reacting face — the partition in `solidgas/statmech.py` supplies one for rutile, and the declared shell thickness changes that mapping without any of the rate laws below being aware of it — and runs the same crossing twice. Three separate things come out, and only one of them is the one usually expected.

**The steady rate does not move.** Any strictly increasing mapping is a reparameterisation of the same crossing in surface coverage, and that coverage is fixed by the rate constants alone. The turnover is identical to ten digits. The uniform model gets the rate right.

**The steady state moves by two orders of magnitude.** What changes is the inventory that holds the face at that coverage. At 600 C the rutile surface runs about 269x above the particle average in the dilute limit — well below its geometric ceiling of `N_O_total / N_s = 1847`, because the compensated isotherm only enriches the surface by `exp(dE/3kT) = 333` and there are 1846 times as many bulk sites to dilute it into. The carrier therefore sits at roughly a two-hundred-and-seventieth of the reduction extent the uniform model assigns it. Degree of reduction is the measured quantity — thermogravimetry, conductivity, diffraction all see the particle, not the face — so this is the half that matters for comparison with experiment.

**On half the descriptor axis the steady state stops existing.** The rutile mapping saturates: above a small inventory the surface is pinned at the (1x2) line-phase deficiency, so neither half-rate depends any more on how reduced the particle is and the imbalance keeps its sign. Every oxide easier to reduce than the balance point then has no crossing at all — it reduces until something outside this pair of rate laws stops it. The boundary coincides with the Sabatier descriptor to machine precision, because the surface ceiling and the peak coverage are both 1/2 here. The correction removes the root rather than moving it, which is not something a perturbative reading of the uniform model anticipates.

The same separation resolves the objection above: once the face may differ from the average, the carrier runs its surface near the (1x2) deficiency while the particle average stays inside the declared point-defect solubility, and the flag clears.

## How much of this depends on the parts that are assumed

The defect model has one degree of freedom left that is declared rather than measured: how thick the subsurface class is taken to be. It matters at the inventories the experiments sit at — going from one oxygen slab to four moves the bulk remainder by more than a tenth and `mu_V` by tens of meV at 95 umol-O/g — because every slab in the shell is charged the first-subsurface energy.

It changes nothing above. The redox steady state does not sit at 95 umol-O/g; it sits at a particle average of order 1e-5, and in that dilute regime the deeper classes are Boltzmann-suppressed, so the enrichment is set by the energies and the site counts rather than by how the shell is sliced: it moves from 269 to 240 across the whole range of declared thicknesses, and none of the three claims turns on that.

Two other knobs used to be here — resolving the shell into separate layers with a declared depth profile, and a mean-field crowding penalty — and both have been removed. The profile shape is now derived from the same DFT triple rather than declared (see [particle-model.md](particle-model.md)), and the penalty was standing in for the cation-sublattice entropy that is now in the free energy itself. The evidence: with the compensation in place, resolving the shell into three layers moved the bulk inventory by 2%, where before it moved it by 30%.

The runaway boundary is untouched for a different reason: it is set by the surface's own ceiling, which the phase ladder fixes and no deeper class can move.

So the three claims above rest on the surface enrichment ratio and the (1x2) line-phase deficiency, both of which are geometry and a reported phase boundary, and not on the parts of the defect model that are least constrained. Checked across all six combinations of resolution and interaction, together with a companion test that the knobs do move the measurement regime — otherwise the robustness claim would be vacuous.

## Verification

`tests/test_redox.py`, 21 gates, none of them a frozen number. The bisection is checked against the first-order closed form across temperature, pressure and descriptor; the swept volcano peak against the analytic Sabatier descriptor and the `alpha/(alpha+beta)` coverage; uniqueness against the strict monotonicity of the two curves; the rate invariance against mappings that have nothing to do with rutile, because it is algebra rather than chemistry. Leaving the point-defect domain, extrapolating the free-energy relation to a negative barrier, running away instead of settling, and finding a root that does not locate the state are each checked to be reported rather than hidden.

## What this model does not carry

It is a single well-mixed particle at one composition with two lumped half-reactions and no adsorbed intermediates, no site heterogeneity beyond the one mapping, and no transport.

Holding `E` fixed while the carrier turns over is a statement about timescales, not only about geometry: the face and the interior are taken to be in equilibrium partition at every instant, so oxygen exchange between shell and bulk is fast compared with both half-reactions. If it is not — if the interior cannot resupply the face on the turnover timescale — this is the wrong model and the problem becomes transport into the particle rather than a boundary condition on it. The partition is the defect workspace's equilibrium answer, so the assumption is inherited from there rather than made a second time here.

The reductant and the oxidant are deliberately unnamed. The rate parameters in `data/redox.json` are identical placeholders for both half-reactions (`log10_A = 8`, `E0 = 1.0 eV`, both labelled *assumed*), so attaching a gas pair would assert chemistry the numbers do not carry. What the model does require of them is stated: unit order in gas and in site, for both halves. A pair whose order is not one — methane is the obvious case — breaks the claim that temperature, pressure and prefactor cancel out of the crossing, because that cancellation is what unit orders buy. The reduction and oxidation channels are assumed simultaneous; a cycled carrier alternating between two gas feeds is a different problem with the same crossing hiding inside it. The linear free-energy slopes are assumed, not measured, and they are the parameters the volcano's width and position are most sensitive to.

## In the browser

`ti_solver.html#redox` is the workspace. Four inputs — the rate-constant ratio K, the surface enrichment E, the point-defect solubility and the face ceiling — and the two figures redraw live: the crossing, and the usable region over (K, E) with a ring on the current setting. The verdict line says which of the three things is true at that point: a usable steady state, a crossing outside the single-phase solid, or no steady state because the face has saturated. Opening the page at K = 0.5 and pressing *depth-weighted, E = 269* is the whole argument in one click — the same K goes from unrunnable to runnable, and the turnover does not move.

The browser engine is `web/redox.js`, and it is deliberately not a port. It carries the closed forms and nothing else; the bisection, the descriptor axis and the linear free-energy overlay stay in Python. `tests/test_redox_port.py` runs the real Python bisection over 455 (K, E) points — with the enrichment entered as the average-to-surface mapping the solver takes — and holds the page's arithmetic against it. The occupancies agree to 1e-12; the rate agrees to the solver's own convergence residual, which is the honest bound because the rate is the mean of the two numbers whose difference the bisection drove to zero.

The two figures are drawn by the same `web/figures_redox.js` that writes `docs/figures/redox_crossing.svg` and `redox_phase.svg`, from the same payload, so the interactive plot and the manuscript file cannot disagree. The one thing the page adds is the "here" ring; the static export never sets it, and a gate checks the exported SVG does not contain it.
