# Steady-state redox

## Question

An oxygen carrier sits in a gas that both reduces and oxidises it. A reductant pulls lattice oxygen out and leaves a vacancy; an oxidant fills a vacancy and puts oxygen back. Where does the solid settle, and how fast does it turn over there?

The model is one nonlinear equation. Nothing in it is new — it is 1970s reaction engineering — and the point is not the method but that the question has not been carried through for this class of material with the solid described honestly.

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

## Domain

`theta` is a point-defect concentration only while the reduced solid stays one phase. Past the declared solubility the vacancies order into crystallographic shear planes and mass action in `theta` stops meaning anything. `steady_state()` reports that as a flag and a warning; it never clips.

This is where the uniform model runs into itself. With symmetric slopes the optimum sits at `theta = 1/2`. Half the oxygen sublattice vacant is a coverage a surface reaches routinely — the rutile (1x2) added-row phase is exactly that — and a bulk composition no oxide holds as a solid solution. In the uniform model those are the same number, so the model puts its own optimum outside the domain of every carrier. That is a defect of forcing one number to be both quantities, not a finding about oxides.

## What the layer-resolved correction actually does

`surface_correction()` takes a mapping from the particle-average vacancy fraction to the fraction on the reacting face — the partition in `solidgas/statmech.py` supplies one for rutile — and runs the same crossing twice. Three separate things come out, and only one of them is the one usually expected.

**The steady rate does not move.** Any strictly increasing mapping is a reparameterisation of the same crossing in surface coverage, and that coverage is fixed by the rate constants alone. The turnover is identical to ten digits. The uniform model gets the rate right.

**The steady state moves by three orders of magnitude.** What changes is the inventory that holds the face at that coverage. At 600 C the rutile surface runs about 1844x above the particle average in the dilute limit, close to its geometric ceiling of `N_O_total / N_s = 1847`, so the carrier sits at roughly a two-thousandth of the reduction extent the uniform model assigns it. Degree of reduction is the measured quantity — thermogravimetry, conductivity, diffraction all see the particle, not the face — so this is the half that matters for comparison with experiment.

**On half the descriptor axis the steady state stops existing.** The rutile mapping saturates: above a small inventory the surface is pinned at the (1x2) line-phase deficiency, so neither half-rate depends any more on how reduced the particle is and the imbalance keeps its sign. Every oxide easier to reduce than the balance point then has no crossing at all — it reduces until something outside this pair of rate laws stops it. The boundary coincides with the Sabatier descriptor to machine precision, because the surface ceiling and the peak coverage are both 1/2 here. The correction removes the root rather than moving it, which is not something a perturbative reading of the uniform model anticipates.

The same separation resolves the objection above: once the face may differ from the average, the carrier runs its surface near the (1x2) deficiency while the particle average stays inside the declared point-defect solubility, and the flag clears.

## Verification

`tests/test_redox.py`, 21 gates, none of them a frozen number. The bisection is checked against the first-order closed form across temperature, pressure and descriptor; the swept volcano peak against the analytic Sabatier descriptor and the `alpha/(alpha+beta)` coverage; uniqueness against the strict monotonicity of the two curves; the rate invariance against mappings that have nothing to do with rutile, because it is algebra rather than chemistry. Leaving the point-defect domain, extrapolating the free-energy relation to a negative barrier, running away instead of settling, and finding a root that does not locate the state are each checked to be reported rather than hidden.

## What this model does not carry

It is a single well-mixed particle at one composition with two lumped half-reactions and no adsorbed intermediates, no site heterogeneity beyond the one mapping, and no transport. The reduction and oxidation channels are assumed simultaneous; a cycled carrier alternating between two gas feeds is a different problem with the same crossing hiding inside it. The linear free-energy slopes are assumed, not measured, and they are the parameters the volcano's width and position are most sensitive to.
