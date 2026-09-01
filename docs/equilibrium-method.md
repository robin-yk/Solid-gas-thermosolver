# Gas-solid equilibrium method

## Scope

The production solver calculates closed, finite-inventory equilibrium for a CO2/H2/CO/H2O/N2/O2 gas in contact with titanium oxides. Temperature, pressure, gas volume, oxide mass, initial oxide, and gas composition are fixed inputs. The output is the stable gas composition and condensed-phase assemblage.

Hydrocarbons, graphite, nitrides, and a rutile TiO2-x solution phase are outside the production registry.

## Active-set formulation

`solidgas/activeset.py` solves the gas in logarithmic element-potential coordinates and enumerates feasible active sets of condensed phases. Every single phase and every phase pair is tested. Three-phase sets are used only to detect degeneracy because a binary Ti-O condensed system admits two generic line compounds at fixed temperature and pressure.

For condensed phase `j`, the reduced cost is

```text
r_j = mu0_j - t_j lambda_Ti - o_j lambda_O
```

It is reported per formula unit and per mole of oxygen removed. An inactive phase is not a numerical optimization variable. Its amount is exactly zero, and its reduced cost measures the distance to stability.

The feasible candidate with clean KKT conditions is the global minimum of the convex problem. Candidates are classified by feasibility and KKT conditions rather than tiny differences in total Gibbs energy.

## Gas and solid chemical potentials

Gas species use ideal-mixture chemical potentials at the specified total pressure. Gas standard states come from NIST-JANAF Shomate coefficients. Titanium oxide line compounds use the Waldner and Eriksson assessment.

Element balances are enforced for C, H, O, N, and Ti. Solid amounts are recovered from these balances after the element potentials and gas composition have been solved.

## Oxygen-potential dual

`solidgas/potential.py` calculates critical oxygen potentials for phase boundaries. Its margin is algebraically the active-set reduced cost per mole of oxygen. It is useful as a fast dual representation, but it is not an independent verification of the Gibbs calculation.

## The reduction line, without the minimisation

There is a second way to ask "does the host survive", and it does not need
the minimisation at all. Titanium does not enter the gas below about
2000 K, so the solid and the gas trade nothing but oxygen, and the whole
question is whether the gas sits above or below one number: the highest
critical potential among the phases more reduced than the host. That
number is a difference of two Waldner free energies — a plain function of
temperature. The gas side is one monotone scalar root, and inverting it
back to a feed composition is closed form.

`activeset.reduction_boundary` returns the phase that would come first and
the H2/CO2 feed that sits exactly on the line;
`activeset.gas_oxygen_potential` returns what a feed is worth before the
solid has touched it, and `activeset.equivalent_co2_fraction` puts that on
the same axis. The cost is a handful of exponentials against sixty Newton
solves for a bisection, which is why the workspace can draw the line at
181 temperatures for free while the margin curve needs a button.

Two things were easy to get wrong here and are gated rather than argued.

**Free O2 is carried.** With hydrogen in the feed it is worth nothing —
1.6e-5 kJ per mol O at 1650 °C, and the first draft left it out on that
basis. Take the hydrogen away and it is worth 5.8 kJ: a CO2/CO mixture at
1650 °C runs 0.2 mol% O2, and pure CO2 at 400 °C is off by 57 kJ without
it, which is enough to put a feed on the wrong side of its own line. It
costs one more exponential, and it brings the total pressure in with it,
because a partial pressure is not a property of a mixture on its own.
`tests/test_activeset_boundary.py` measures both halves: the pressure
moves a hydrogen-free answer by tens of kilojoules over a decade of P and
an H2/CO2 answer by under a millijoule, so the boundary curve itself can
be drawn without asking.

**The marker is the feed, not the equilibrium.** Once a closed charge has
equilibrated against a reducible solid, the gas is buffered onto this very
line and no longer says how reducing it started out. So the figure places
the gas equilibrated on its own — water-gas shift among its own members,
no solid — which is also the only version of the question that has a
temperature-independent answer worth drawing.

Against the long way round the shortcut lands within 2e-4 relative of a
bisection of the full minimisation over feed composition, and within 1e-10
kJ per mol O of that minimisation's own oxygen multiplier wherever the
host still stands. The one family held to a looser number is the one where
the *solver* is the coarser of the two: a gas with no reduced member in it
reaches equilibrium with CO at 3e-12 of the mixture, and a Newton solve on
a scaled residual of 1e-13 is working at its own floor, while the closed
form has no trace species to resolve.

## Outputs

The browser workspace reports:

- active condensed phases
- gas composition and CO2 conversion
- Ti3+ fraction and oxygen deficit
- oxygen transferred between the solid and gas
- KKT reduced cost for every inactive phase
- element-balance residuals
- comparison with the high-precision oracle

The complete equation sheet is generated by `scripts/export_equations.py` and inlined into `docs/ti_solver.html`.

## In the browser

`ti_solver.html#thermo` is the workspace. Four figures on the Results tab draw what the KKT tab reports in numbers.

The first is one bar per excluded phase: the reduced cost it would have to give up to appear at the condition on the panel. The shortest bar is the phase that would come first, and it is drawn in a second colour because it is the only one that matters for whether the answer is safe. It redraws on every run and costs no extra solve — the numbers are already in the result on screen.

The second sweeps temperature at the same charge and plots the smallest of those bars at each point, with the winning assemblage as a band along the bottom: shaded where two phases coexist, pale where one does. The curve touches zero at each edge of that band, so the sawtooth and the band are the same crossings seen two ways. In the Magnéli window the band walks the ladder down — TiO2+Ti10O19, Ti10O19, Ti10O19+Ti8O15, Ti8O15, and on — which is the headline behaviour of this solver and had never been drawn from live numbers. It is behind a button because it is 121 solves, about half a second; once drawn it follows every later run. The band is sampled, so a tie line narrower than one step — about 11 °C at 121 points — goes unshaded even though the curve still comes down to zero there; `tests/test_thermo_figures.py` refines every place two single-phase fields appear to meet and finds the tie line the grid stepped over.

Both are drawn per formula unit, not per mole of oxygen. Per mol O divides by the oxygen deficiency relative to the host, `def_j = rho_host*t_j - o_j`, and that denominator changes sign as soon as a candidate is less deficient than the host — over the conditions this page admits, 41% of reduced costs come back negative for that reason alone, and it is exactly zero for a candidate that differs from the host by no oxygen at all. The KKT test is stated per formula, where `r_j > 0` means the phase is absent, and `tests/test_thermo_figures.py` pins both halves of that: per formula never goes negative at a reported solution, and per mol O often does. Plate e2 keeps the per-mol-O axis, which is the comparable quantity at the one frozen condition it reports.

The third is the reduction line in the units of a mass-flow controller rather than kilojoules: how much CO2 the feed has to carry before the host stops giving up oxygen, against temperature, with the panel's own feed marked as a square. It answers the same question as the margin curve — the header on one says "54.9 kJ from the nearest phase", the header on the other says "7918 times more CO2 than it would take" — but nobody sets a kilojoule on a rig. Every gas can be placed on that axis: a mixture of H2 and CO2 lands at the fraction it was mixed at, exactly, because the axis is the same oxygen balance run backwards; anything else lands at the H2/CO2 mixture that would be equally reducing, and the payload marks which of the two it is showing. Two feeds have no place on it at all and say so rather than guessing — a gas with no carbon and no hydrogen has no couple, and a gas short of the oxygen to make every carbon a CO is not a C-H-O mixture in this registry. The curve is closed form, so it redraws with every run and needs no sweep.

The fourth is what the gas itself did: CO2 conversion against temperature, out of the same sweep the margin curve runs, with the runs where the assemblage is no longer the host alone shaded. Conversion is counted against the feed, so oxygen that left the gas is counted whether it went into CO or into the solid, and the shaded run is where part of the answer is the solid being eaten.

Both curves are sampled and the workspace is not, so both are handed the exact value at the panel's temperature rather than reading the nearest sample — eleven degrees is worth a fifth of a percent of conversion and eight percent of the boundary, and a header that far from the panel would be reporting a condition nobody set.

The payload builders in `web/figures_thermo.js` only reshape a solver result; the gate writes the same reshaping out longhand in Python and holds the two to equality, so nothing can grow its own arithmetic between the engine and the drawing. The sweep itself goes through `window.ThermoBridge.solveAt`, which re-runs the panel's own charge with one argument changed, and the boundary marker goes through `Solver.feedOnBoundaryAxis`, so the page owns no physics of its own — including the reasons a feed can have no answer, which are gated on both sides.
