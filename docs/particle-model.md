# The particle model

What it answers: a titration says the sample lost so many micromoles of
oxygen per gram.  Where is that oxygen missing from?

This is the second engine in the repository to answer that question.  The
first split the particle into three site classes — bridging row, one
subsurface slab, everything else — and solved a three-term lattice-gas
partition.  This one keeps the measurement and the material and replaces
the physics underneath, following the construction used for reduced
ceria nanoparticles.  The reasons are in **What changed and why** at the
end; they are not cosmetic.

## The construction

**One chemical potential.**  Equilibrium means every oxygen site in the
particle shares one vacancy chemical potential μ.  μ is not an input:
it is whatever value makes the sites together hold the measured
inventory.  Everything else follows from it.

**Vacancies pay for their electrons.**  Pulling a neutral oxygen out of
rutile leaves two electrons, and in rutile they localise as two Ti³⁺
small polarons rather than delocalising.  A vacancy is therefore one
event on the oxygen sublattice and two on the cation sublattice, and all
three carry configurational entropy.  Rutile has one Ti per two O, so
θ(Ti³⁺) = 4 θ(V), and the chemical potential is

    μ(θ) = E(z) + kT [ ln(θ/(1−θ)) + 2 ln(4θ/(1−4θ)) ].

**The cation sublattice runs out at θ = 1/4.**  That composition is
TiO₁.₅, which is Ti₂O₃.  Nothing in the code imposes it — the second
logarithm diverges there because there is no fourth Ti³⁺ to make.  Halve
the gap to the quarter and μ rises by 2kT ln 2, without bound.

**The dilute isotherm is p(O₂)^(−1/6).**  The two logarithms add to
3 kT ln θ when θ is small, which is the textbook exponent for a doubly
ionised vacancy compensated by electrons.  A free energy written on the
oxygen sublattice alone gives p(O₂)^(−1/2), the exponent for a *neutral*
vacancy, and an enrichment between two sites that is the cube of the
right one.

**Depth, not classes.**  The formation energy relaxes towards the bulk
value over a decay length:

    E(z) = E_bulk − ΔE_seg · exp(−z/ξ).

Two parameters, and the DFT set has three energies with the bulk as the
asymptote, so both are determined and no shape is assumed:

| preset | ΔE_seg | ξ | in d₁₁₀ layers |
|---|---|---|---|
| sX (Li, Guo & Robertson 2015) | 1.31 eV | 0.543 nm | 1.67 |
| GGA (same reference) | 1.08 eV | 0.725 nm | 2.23 |

The ceria study assumed 0.7 nm on the grounds that the near-surface
relaxation depth is two or three atomic layers.  Ours is derived and
lands in the same place, which is a check on both.

**Layers, not a smear.**  The interior is a stack of d₁₁₀ oxygen layers,
layer *k* centred at depth *k*·d₁₁₀, and the engine sums it exactly:
explicitly near the surface, in closed form once the layers stop
differing from the bulk.  Smearing the stack into a continuum — which
is what the ceria formulation does — costs 0.5 % in the interior mean
here, because the energy changes by a large factor from one layer to the
next.  Both forms are implemented and a gate holds them apart.

**The bridging row is a quarter of a layer.**  A full d₁₁₀ slab of
rutile (110) holds four oxygens per 1×1 cell; the bridging row holds
one.  It is carried as an areal site density, not folded into the
interior, because it is the cheapest site in the particle and giving it
a bulk density would quadruple it.

## What rutile adds that ceria does not have

Two things cut the profile off before the cation sublattice fills, and
both are real crystallography rather than model choices.

**The (110) face reconstructs.**  Past a bridging-row coverage near 0.17
the 1×1 vacancy lattice gas gives way to the added-row Ti₂O₃ (1×2)
phase.  That phase is a line compound, so the transition is first order
and μ pins while the surface converts.  At 0.9 μm and 600 °C the window
is 33.2 → 37.7 μmol-O/g; above it the row is fully reconstructed and
stops responding.

**The bulk throws shear planes.**  Rutile does not dissolve point
defects all the way to Ti₂O₃.  Past x ≈ 0.008 crystallographic shear
planes form, so μ pins again and further inventory leaves the
point-defect count as extended defects.  At 600 °C the point-defect
ceiling for a 0.9 μm particle is **112.4 μmol-O/g**.

That last number is worth reading twice.  The shipped condition is
95 μmol-O/g — 85 % of the way to the ceiling.  The particle cannot hold
much more oxygen vacancy as point defects whatever the gas does.

## What the model says at the shipped condition

0.9 μm, 600 °C, 95 μmol-O/g, sX energies:

| | |
|---|---|
| μ | 0.2336 eV |
| bridging row | fully reconstructed, θ = 0.5, holding 6.8 μmol-O/g |
| first subsurface layer | θ = 0.0654 |
| bulk plateau | θ = 0.00334 |
| interior | 88.2 μmol-O/g |
| within 1 nm of the surface | 0.7 % of the volume, 12.4 % of the oxygen |

The surface saturates long before the bulk does — it costs 1.31 eV less
— so by the time a titration registers anything it is already
reconstructed and is a small fixed reservoir, not the main one.  The
near-surface enrichment is real but it is a factor of twenty over one
nanometre, not a shell that holds the inventory.

## The assumptions, and what each is worth

**Is the equilibrium reachable?**  With the literature bulk migration
barrier of 0.65 eV, a 0.9 μm particle homogenises 99 % of the way in
0.54 ms at 600 °C — six orders below a 300 s exposure.  The margin is
the number to quote: the barrier would have to reach **1.645 eV**, about
an eV above the literature value, before transport became the limit.
The group's own recovery measurement implies a higher effective barrier
than the literature one, so the margin is not infinite.

**Is local electroneutrality allowed?**  The model assumes each vacancy
drags its two polarons with it pointwise.  The ceria study tested this
against a full Poisson–Boltzmann treatment and found it held over more
than 95 % of the radius — and said explicitly that the collapse was a
consequence of how concentrated that system was.  This one is 13× more
dilute, which is the direction space charge grows.  Computing it rather
than inheriting it:

    Debye length at the bulk plateau   0.57 nm
    segregation depth ξ                0.54 nm

They are the same length.  Screening and chemistry cannot be separated
here, and over ε_r from 50 to 170 the ratio only moves from 0.74 to 1.37.
So the outer few nanometres of this profile carry an error the model
does not quantify, and the ceria result says that error runs towards
*more* surface enrichment, not less.

**Does the polaron ratio matter?**  At the bridging row the local cation
reservoir is one Ti³⁺ site per vacancy, not four; ratio four is the
statement that the polarons drain into the subsurface, which is what DFT
reports for reduced rutile (110).  In the dilute limit the ratio enters
only as an additive constant, so it shifts μ by 2kT ln 4 and leaves the
distribution alone.  Where it is not a gauge is near the cap: the first
subsurface layer sits close enough to feel it, and the enrichment moves
by 15 % between ratio 1 and ratio 4.

**Facets.**  The DFT energies are for (110) and the particle is a
sphere.  The ceria study makes the same idealisation with (111)
energies.  A faceted particle would weight the surface term by facet
area and give each facet its own ΔE_seg; nothing here does that.

**θ counts oxygen sites.**  TiO₂₋ₓ has x = 2θ.  The ceria report calls
this out as the easiest factor of two in the subject to lose, and it is.

## What changed and why

The three-class engine that came before this one wrote the free energy
on the oxygen sublattice alone.  Two things followed, and both were
wrong in ways that mattered:

- its dilute isotherm was p(O₂)^(−1/2), a neutral vacancy, and its
  surface-to-bulk enrichment was the cube of the compensated one;
- it had no cap, and at the shipped condition it put the first
  subsurface layer at θ = 0.954 — four times past Ti₂O₃, locally past
  TiO, on a branch that had no business existing.  It is now 0.065.

The depth profile changed for a smaller reason: the old one placed
intermediate layers at a declared fraction of the subsurface-to-bulk
interval, and the dataset labelled that fraction *assumed*.  The
exponential reproduces it to within a few hundredths without being told
to, and has no free shape parameter.

Two things were dropped rather than ported.  The mean-field coverage
penalty ω had no literature value and its own provenance note said its
magnitudes were chosen to bring the first subsurface layer out of
near-saturation — which is exactly the symptom the missing cation
entropy was causing, so it was standing in for the term that is now
there.  And the sampled-surface-isotherm extrapolation, which continued
the 1×1 branch linearly in logit past the sampled range, is not needed:
the surface phase construction says what happens there instead.

## Verification

Three independent implementations have to agree before a number is
quoted.  For this engine:

- `solidgas/particle/` — the engine, double precision, analytic
  one-sided brackets on two branches;
- `scripts/oracle_particle.py` — fifty decimal places, no engine import,
  one logistic substitution and plain bisection, summing every layer of
  the particle to the centre;
- `tests/test_particle.py` — the gates, including full oracle parity.

The Crank sphere series independently reproduces the Fourier number the
ceria study solves for at 20 % uptake (τ = 0.003912), which is a check
on the transport half against a source outside this repository.

## Reading

- Crank, *The Mathematics of Diffusion*, 2nd ed. — the sphere solution.
- Skorodumova et al., *Phys. Rev. Lett.* **89**, 166601 — vacancy
  formation, polaron localisation, and the sesquioxide limit.
- Maier, *Physical Chemistry of Ionic Materials* — the space-charge
  framework; Guo & Waser for the grain-boundary case.
- Li, Guo & Robertson, *J. Phys. Chem. C* **119** (2015), DOI
  10.1021/acs.jpcc.5b02430 — the rutile (110) formation energies.
- Reticcioli et al., *Phys. Rev. X* **7**, 031053 and *Phys. Rev. B*
  **98**, 045306 — polaron localisation and the (1×2) under reduction.
- Onishi & Iwasawa, *Surf. Sci.* **313** (1994) L783 — the added-row
  structure.
- Iddir et al., *J. Phys. Chem. C* (2010), DOI 10.1021/jp107986a — the
  bulk migration barrier.
