# Vacancy inventory against surface capacity

## The question a titration can answer on its own

Reducing rutile under H2 and integrating the evolved H2O gives one number per
sample: the total oxygen removed, in µmol-O per gram. It does not say where in
the particle that oxygen came from, and for catalysis the difference matters —
surface oxygen is taken by the gas directly, oxygen below the surface has to
reach it first.

There is no measurement in this work that resolves depth. What the measurement
*can* support is a comparison against how much oxygen the surface could hold,
which is lattice geometry.

## The capacity

Rutile (110) exposes two kinds of oxygen: three-coordinate oxygen in the plane,
and a row of two-coordinate bridging oxygen standing above it. Only the bridging
row is counted here. It is the oxygen the gas takes, and counting the in-plane
oxygen as well would inflate the capacity in the direction that flatters the
conclusion.

The (110) 1×1 cell is `c · a√2` = 19.22 Å² and carries one bridging oxygen, so
the areal density is 5.20 × 10¹⁴ cm⁻². Multiplying by the exposed area gives the
capacity per gram. The area comes from a measured BET value when there is one,
and otherwise from the smooth-sphere estimate 6/(ρd), which is a lower bound on
the real area of a rough particle and therefore the conservative choice.

| | 0.9 µm | 1.6 µm |
| --- | ---: | ---: |
| Exposed area, m² g⁻¹ | 1.569 | 0.883 |
| Whole bridging row, µmol-O g⁻¹ | 13.56 | 7.63 |
| (1×2) reconstructed surface, µmol-O g⁻¹ | 6.78 | 3.81 |

The second row is the added-row Ti2O3 stoichiometry: one added unit per (1×2)
cell is one vacancy per two bridging sites, so half a bridging monolayer. That
0.5 is the only structural literature value in the calculation. Everything else
is arithmetic on the lattice constants.

## The bound

For a measured removal *V*,

    f_surface ≤ C / V,        f_below ≥ 1 − C / V

with *C* the capacity above. It is an inequality, not an estimate, and both
capacities are worth quoting: the reconstructed surface is what is physically
reachable, and the whole bridging row is the weakest claim that needs no
structural input at all.

At 95 µmol-O g⁻¹ on a 0.9 µm particle the measured removal is 7.0× the whole
bridging row. At most 7.1 % of it can be surface oxygen; at least 92.9 % has to
be below the surface. Across the reduction series the bound falls from 22.6 %
at 30 µmol-O g⁻¹ to 3.1 % at 220.

## What is deliberately not done

The below-surface oxygen is not divided among subsurface point defects, bulk
point defects and extended defects. Doing so needs vacancy formation energies as
a function of depth, and the resulting split would be a model output presented
where a measurement is expected. Literature DFT gives a surface / first
subsurface / bulk triple of 0.0 / 0.59 / 1.31 eV, which is a strong depth
dependence; it is recorded in `data/literature_dft.json` for context and is not
read by any code in this repository. A gate enforces that.

For the same reason there is no surface Monte Carlo, no radial profile, no
transport clock and no shear-plane solubility ceiling here. Each of those was
implemented and removed: they reproduced the bound above and added assumptions
the measurement cannot check.
