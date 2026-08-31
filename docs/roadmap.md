# Roadmap

## Thermodynamic model

1. Add a rutile TiO2-x solution phase with assessed vacancy thermodynamics.
2. Add graphite before releasing any CH4 equilibrium result.
3. Extend the gas and condensed registries only with source-traceable data and oracle coverage.

## Defect model

1. Replace literature placeholder energetics with project DFT values through `data/rutile_dft.json`.
2. Add reconstructed rutile (110) energetics rather than extrapolating the unreconstructed surface Hamiltonian.
3. Put the defect workspace, the JavaScript mirror and the figure plates on `solidgas/particle/` so one engine answers the partition question. Until then two engines coexist; they agree to about 1.5 percent at the shipped condition, the difference being layer resolution.
4. Reconcile the coarse-grained KMC layer with the compensated free energy. Its exchange rates satisfy detailed balance against a product-binomial measure on the oxygen sublattice alone, which is the neutral-vacancy equilibrium; the thermodynamic model is no longer that, so the long-time kinetic limit and the equilibrium partition currently disagree.
5. Weight facets. The DFT energies are (110) and the particle is a sphere.

## Where this stands against the ceria treatment

The particle engine follows a four-model study of vacancy distribution in
a pre-reduced 100 nm ceria nanoparticle, and the mapping to rutile is
exact for the part that matters: fluorite CeO2 and rutile TiO2 both have
one cation per two oxygens, so a vacancy converts four times its own
site fraction of the cation sublattice in both, and both cap at a
quarter — Ce2O3 there, Ti2O3 here.

| ceria model | here |
| --- | --- |
| 1. Transient diffusion, Dirichlet surface | `transport.py`; the Crank series is shared and reproduces its Fourier number for 20 percent uptake to six figures |
| 2. Closed-system relaxation, zero flux | not implemented; the homogenisation time it measures is reported instead, and at 0.54 ms against a 300 s exposure the relaxation has no room to be interesting |
| 3. Free-energy segregation, local electroneutrality | `profile.py`, with the decay length derived from DFT rather than assumed, layers summed rather than smeared, and two rutile-specific cutoffs the ceria system has no analogue for |
| 4. Maier space charge, Poisson-Boltzmann | not implemented; `spacecharge.py` computes the screening length that decides whether it is needed, and finds it comparable to the segregation depth |

Model 4 is the open item, and it is more pressing here than it was
there. The ceria study found local electroneutrality held over 95
percent of the radius and said explicitly that the collapse followed
from how concentrated that system was. This one runs 13 times more
dilute, and the Debye length at the bulk plateau (0.57 nm) comes out
the same size as the segregation depth (0.54 nm). A full self-consistent
solve is what would turn that from a caveat into a number.

Two smaller items from the same source: convert between reduction
conventions carefully — theta counts oxygen sites and TiO(2-x) has
x = 2 theta — and replace the illustrative segregation gap with values
fitted to depth-resolved measurement. The predicted enrichment shell is
in principle distinguishable from a flat profile by XPS depth
profiling, which is the natural validation route for both systems.

## Kinetics and transport

1. Replace the CO2-refill placeholders with fitted class-specific parameters.
2. Add radial coarse cells and Arrhenius exchange rates in detailed balance with the equilibrium model.
3. Verify that the long-time kinetic limit reproduces the canonical equilibrium distribution.

These additions must preserve the existing separation between model data, production solvers, independent oracles, browser ports, and generated pages.
