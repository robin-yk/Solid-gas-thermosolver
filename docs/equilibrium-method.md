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

`activeset.critical_potential` gives the oxygen potential at which one line
compound gives way to another. It is algebraically the reduced cost per mole of
oxygen, so it is a second reading of the same optimality condition and not an
independent verification of the Gibbs calculation. A standalone module that
carried this route separately was removed: it answered a question the active set
already answers, and its gas registry disagreed with the solver's below 800 °C.

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

## In the browser

`docs/index.html` displays the result; it does not compute it. Three figures on
the equilibrium workspace draw the record the Python package wrote: the
equilibrium CO2 conversion against temperature, the margin to the nearest
reduced phase, and the CO2 fraction at which rutile stops surviving. Changing
the temperature selects a different row of that record. There is no solver in
the browser.
