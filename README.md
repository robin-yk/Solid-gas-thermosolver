# solid-gas-thermosolver

Equilibrium solver for a gas phase in contact with a reducible oxide, built to
answer one question: **can a reverse water-gas-shift feed reduce rutile TiO2 to
a Magneli phase?**

The answer is no across 500-1500 C, but the margin is not uniform, and the
place where it is thinnest is not yet fully tested. See *Where this is soft*.

## One route reports, the other checks

Every number in `results.json` is a **full-composition Gibbs minimisation**:
the total Gibbs energy minimised over the whole species set, gases ideal and
mixing, condensed phases pure at unit activity, subject to elemental balance.
Nothing is read off a phase-boundary table, and no point is filled in by the
faster route.

`solidgas/gibbs.py` does that in reaction-extent coordinates. The change of
variables is exact - same objective, same feasible set - but the constant part
of G is removed analytically rather than carried through the arithmetic, and
that is what makes the inert case computable at all. Sealed N2 pulls rutile off
stoichiometry by a mole fraction near 1e-38; in species space that is destroyed
before the optimiser sees it, because `n_TiO2 = n_Ti - 10*xi` with `n_Ti` at
1e-3 rounds to `n_Ti` exactly. Minimising the reaction Gibbs energy directly
keeps the extent as a primary variable at its own magnitude.

`solidgas/potential.py` is the **oxygen-potential route** and is kept as the
cross-check. Titanium does not enter the gas below about 2000 K, so the solid
exchanges nothing but oxygen with it; that lets each condensed phase's critical
potential be computed as a plain function of temperature, with no optimiser in
it. It is roughly 200x faster, and it is still not what gets reported.

The reason is not that it is wrong. It is that the two agree nearly everywhere,
which is exactly what makes a value from the wrong one impossible to spot by
reading it. So the routes are separated in code, every row of `results.json`
carries a `method` field, and `tests/test_results.py` refuses any row that does
not say `gibbs_min`. A rule stated in prose survives one sitting; a rule with a
test behind it survives the next person.

The cross-check runs on every row and is recorded there. It is a real test, not
a comparison of a number against itself: where the minimisation leaves rutile
sitting with a reduced phase, coexistence fixes the gas oxygen potential
exactly, and the boundary value comes from the phase data alone. The minimised
gas has to land on it.

```bash
python3 run_all.py      # 24 points -> results.json, every row stamped gibbs_min
python3 build_site.py   # embeds that same file in index.html
python3 -m pytest tests/ -q
```

## Results

100 mg TiO2, 25 mL, 1 atm. Every number below is `gibbs_min`, straight out of
`results.json`.

**Ti(3+) as a percentage of total titanium**

| feed | 500 C | 700 C | 900 C | 1100 C | 1300 C | 1500 C |
|---|---|---|---|---|---|---|
| CO2/H2 = 1:1 | 1.0e-17 | 1.0e-17 | 4.6e-16 | 6.7e-18 | 5.0e-18 | 1.1e-42 |
| pure H2 | 0.0005 | 0.0168 | 0.1702 | 0.8694 | 2.82 | 6.37 |
| pure CH4 | 0.1678 | 5.59 | 33.33 | 36.59 | 33.33 | 33.33 |
| sealed N2 | 1.5e-36 | 1.7e-26 | 7.2e-20 | 3.9e-15 | 1.4e-11 | 8.7e-09 |

**Phases that appear**

| feed | reduced phases |
|---|---|
| CO2/H2 = 1:1 | none - rutile throughout |
| pure H2 | Ti10O19 |
| pure CH4 | Ti10O19, Ti6O11, Ti5O9 |
| sealed N2 | Ti10O19 |

RWGS never reduces rutile anywhere in the range. Pure H2 sits on the
TiO2/Ti10O19 buffer at every temperature, converting more of the charge as it
gets hotter but never getting past the first step. Pure CH4 is the only feed
that drives the solid off rutile entirely.

Sealed inert gas does nothing, and the minimisation now says so in its own
right rather than by a capacity argument: the oxide gives up oxygen until the
gas reaches the boundary potential, and that leaves between 1e-36 and 1e-8
percent of the titanium reduced across 500-1500 C.

**The cross-check.** Where the minimisation leaves two phases together the gas
must sit exactly on their boundary, a value computed from the phase data with
no optimiser in it. Across the twelve buffered rows the worst disagreement is
**4e-5 kJ per mol O**. The inert case agrees with the oxygen-potential route on
pO2 to five significant figures at all six temperatures.

**Where the two routes part company.** At 1100 C under CH4 the minimisation
puts Ti5O9 + TiO2 below Ti6O11 by 0.6% in reaction Gibbs energy, and the
oxygen-potential route's ladder walk cannot represent that: it assumes the
series is monotonic in n, and here it is not - the intermediate members lie
above the convex hull. The Gibbs value is the reported one. The gap is well
inside the assessment's own scatter, so the *phase identity* at the CH4 hot end
is not resolved by this data; the Ti(3+) level, 33 to 37%, is robust either way.

### How far RWGS sits from reduction

Margin in kJ per mole of oxygen removed, above the phase that forms most
easily out of untouched rutile. Positive means nothing reduces.

| T (C) | 500 | 700 | 900 | 1100 | 1300 | 1500 |
|---|---|---|---|---|---|---|
| **margin, from the minimised gas** | **+78.44** | **+62.92** | **+54.74** | **+46.17** | **+37.03** | **+27.30** |

These are measured against Ti10O19, the shallowest member in `ACTIVE_TI_PHASES`
and so the easiest phase to form out of untouched rutile. Against Ti20O39, the
highest-n member the assessment carries, the 1500 C figure falls to +18.3 - and
that spread is the line-compound model reaching its limit, not a numerical
problem.

Two corrections cancel in part, and both matter:

- Bringing in the shallow Magneli members costs about 14 kJ/mol-O at 1500 C.
  They form more easily than Ti4O7, so any margin measured against Ti4O7
  alone is an overestimate.
- Moving the condensed phases from NIST to Waldner *gains* about 18 kJ/mol-O.
  The two assessments put rutile 6.0 kJ/mol apart and Ti4O7 12.0 kJ/mol apart,
  and the difference does not cancel in the reduction reaction.

Taking the NIST margin and applying a Waldner-derived shift would give 8.7 at
1500 C, which is not a number this repository produces: it mixes the two
bases in exactly the way that the offsets above make invalid. On a single
consistent basis the answer is +18.3.

**RWGS does not reduce rutile anywhere in 500-1500 C**, and the margin at the
hot end is now measured rather than assumed.

## Layout

```
solidgas/
  shomate.py      NIST-JANAF coefficients; H, S, Cp, mu0, Kp
  waldner.py      Waldner & Eriksson Ti-O assessment: the Magneli series
  equilibrium.py  species, element matrix, G_total, constrained solve
  gibbs.py        the reported route: Gibbs minimisation in extent coordinates
  potential.py    oxygen-potential route, kept as the cross-check only
run_all.py              the 24 reported points -> results.json
run_Ti4O7_RWGS.py       CO2/H2 sweep, writes Ti4O7_RWGS.png + results_rwgs.json
run_Ti4O7_reducing.py   pure-H2 and pure-CH4 sweeps, writes results_reducing.json
solver.js               the same thermodynamics for the browser
export_thermo.py        dumps the coefficients as thermo_data.json
site_template.html      page shell; build_site.py inlines the two into index.html
tests/                  data checks, solver checks, export-drift checks
```

```bash
pip install -r requirements.txt
python3 -m pytest tests/ -q
```

## Thermodynamic data

Two sources, kept strictly apart. **Gases** are NIST-JANAF 4th ed. (Chase,
1998). **Titanium-bearing condensed phases** are Waldner & Eriksson, Calphad 23
(1999) 189-218, Appendix pp. 216-218, which carries the whole Magneli series
where NIST stops at Ti4O7.

They are not interchangeable: the assessments put rutile 6.0 kJ/mol apart,
Ti4O7 12.0 kJ/mol apart and Ti2O3 5.3 kJ/mol apart, all comparable to the
margins being measured. `solid_mu0(name, T, source=...)` makes the choice
explicit and a test pins the offsets so the two can never quietly merge.

Heat capacity, which is measured rather than assessed, agrees between them to
under 1% across the working range for every phase both carry - that is the
check that the Waldner integration is right.

The NIST entries below are still used for the `equilibrium.py` route and
verified against the published tables.

| phase | O/Ti | dHf298 (kJ/mol) | fit range (K) |
|---|---|---|---|
| TiO2 rutile | 2.000 | -938.722 | 298-2000 |
| Ti4O7 | 1.750 | -3404.525 | 298-1950 |
| Ti3O5 (beta) | 1.667 | -2446.192 | 298-2050 |
| Ti2O3 (solid) | 1.500 | -1520.884 | 470-2115 |

Notes on the two added this round:

- **Ti3O5** has two polymorphs in NIST. Beta overtakes alpha at 450 K, far
  below the 773 K floor here, so beta is the one in `SPECIES`; alpha is kept as
  `Ti3O5_alpha` for the polymorph test only.
- **Ti2O3** is listed twice: a solid (dHf = -1520.884) and a liquid/glass
  (dHf = -1418.46). Only the solid belongs in a solid-gas equilibrium. The
  liquid value is 102 kJ/mol less stable and would keep Ti2O3 from ever
  appearing. The solid's 470-2115 K fit covers the whole working range; Cp, S
  and H each jump at the 470 K join across the metal-insulator transition, but
  G is continuous there to 0.0006 kJ/mol, which is what equilibrium uses.

## Where this is soft

- **The series does not converge in n, and that is physics rather than a
  gap.** Each higher member forms more easily than the last, all the way to
  Ti20O39. That is not a hint that another line compound is missing: the
  n -> infinity limit is oxygen-deficient rutile TiO2-x, which no set of line
  compounds can represent. Waldner models rutile as a solution,
  (Ti4+)(O2-,Va2-)2 with two interaction parameters, and that solution is not
  implemented here. Ti20O39 should be read as a bound on how far the margin
  can move, not as the phase that would actually appear.
- **Ti7O13 cannot be placed.** Its formation enthalpy is misprinted in the
  paper as -3415327 J/mol, which is Ti4O7's value to four figures; the series
  is linear in n to under 600 J and gives -6251260 instead. That repair is
  good enough to keep Ti7O13 out of the way but not to seat it exactly - the
  residual is the same size as the gaps between neighbouring steps, and it
  shows up as a small reversal in the stepwise ladder. It does not touch any
  margin quoted here, all of which are measured straight from rutile.
- **No carbon deposition.** Graphite is not a species, so the CH4 results
  assume all carbon leaves as CO. The equilibrium gas comes out at H2:CO = 2:1
  exactly, which is what you get when carbon has nowhere else to go - that
  ratio may be an artefact of the species list rather than a result. The CH4
  numbers should not be quoted until graphite is in.
- **Condensed phases are pure and immiscible.** Real Magneli formation runs
  through a homologous series with mutual solubility, and at the reduced end
  through oxygen-deficient rutile TiO2-x, which no line-compound set can
  represent. The n -> infinity limit of the series is not another phase, it is
  rutile non-stoichiometry.
- **Equilibrium only.** Nothing here says how fast any of it happens.

## The pages

**<https://robin-yk.github.io/Solid-gas-thermosolver/>** is a portal over two tools.

**`oxide_tool.html` &mdash; oxide reduction under a reacting gas.** The Gibbs minimisation, in the
browser. Ti&ndash;O and Ce&ndash;O under a CO2/H2 feed; the phase assemblage is an output. Carries
its own verification tab (the same condition solved by the oxygen-potential route, with the
difference shown), the coefficients it is actually running on, the element balance, the solver's
own state, the full derivation, and the solver source printed from the running page. Self-contained:
no external requests.

**`tiox.html` &mdash; the oxygen-potential dashboard.** The fast route, swept continuously across
temperature with the phase boundaries drawn. Its reported values are read from `results.json`, not
recomputed there, and a test compares the embedded copy against the file so the page and the table
cannot drift apart.

```bash
python3 run_all.py         # results.json   - the 24 reported points
python3 export_thermo.py   # thermo_data.json
python3 export_oxide.py    # oxide_data.json
python3 oxide_reference.py # oxide_reference.json - what the JS port must reproduce
python3 build_site.py      # tiox.html
python3 build_oxide.py     # oxide_tool.html
python3 build_portal.py    # index.html
```

### Porting the engine to the browser

`oxide.js` is `solidgas/gibbs.py` carried across: reaction-extent coordinates, the constant part of
G removed analytically, the analytic gradient. Two things could not be copied and are worth naming.

SciPy's SLSQP does not exist in a browser. The change of variables has already removed every
equality constraint, so what replaces it is a projected gradient with a Barzilai-Borwein step and
backtracking, over the log of the extents; `G_rel` returns a large number on an infeasible
composition, so a line search that only accepts a decrease cannot step outside the feasible set.
`tests/test_oxide_port.py` holds the two implementations to five significant figures - the worst
disagreement on any species actually present is 3.5e-6, across both oxide systems and six
temperatures each.

The WGS page's `shomate(c, T)` reads seven coefficients and lands on the absolute potential because
NIST's F already carries the formation term; this package keeps eight and adds dHf explicitly.
Rather than generalise that function - the gas path is the part already checked against an Excel
workbook - `export_oxide.py` folds the constant in as `F' = F - H + dHf`, and the evaluator is
copied unchanged.

### A defect, stated

An absent phase should have exactly zero moles. In logarithmic extent coordinates it cannot: the
logarithm of zero does not exist, so absent phases come back at 1e-26 to 1e-17 percent instead of at
nothing. Telling that residue apart from a phase genuinely present in a trace amount then needs a
second argument - under sealed N2 a share of 1e-36 percent is the entire answer, and under RWGS the
same size is nothing - which is a thing that has already been got wrong once here.

The fix is to enumerate which phases are present, hold the absent ones at exactly zero, and compare
the total Gibbs energy across those cases; the degrees of freedom are small enough that the
enumeration is cheap. Until it lands, both pages state a detection limit and report below it as
below it, rather than printing a number that looks like a measurement.
