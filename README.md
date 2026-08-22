# solid-gas-thermosolver

Equilibrium solver for a gas phase in contact with a reducible oxide, built to
answer one question: **can this gas reduce this titanium oxide?**

For a reverse water-gas-shift feed the answer is no across 500-1500 C, now
with the margin reported as what it mathematically is: the KKT reduced cost
of the absent phase.

## The production route: active-set Gibbs minimisation

`solidgas/activeset.py` is the solver the released page runs, and it exists
to enforce one discipline: **phase stability is decided by KKT reduced costs,
never by the size of a mole number.** Gas species live in logarithmic
element-potential coordinates; condensed phases are enumerated as active sets
(every single and every pair - the phase rule allows two generic line
compounds in an M-O system at fixed T and P, and triples are probed only as
degeneracy detectors). A phase left out of a candidate is not a variable at
all, so its mole number in the result is exactly 0.0, and every inactive
phase carries r_j = mu0_j - t_j*lambda_Ti - o_j*lambda_O, reported per
formula unit and per mole of oxygen. The program is convex, so the candidate
that is feasible with clean KKT *is* the global minimum - candidates are
never ranked by total-G values double precision cannot resolve.

Verification is layered, and the layers are genuinely independent:

- `oracle_tio.py` - an 80-digit mpmath oracle sharing only the coefficient
  tables with production, solving a different formulation (element-potential
  root solve per candidate). Its output is
  `reference_results_high_precision.json`; `tests/test_activeset.py` holds
  production to it (worst deviation observed: 2e-13 on gas fractions,
  8e-12 kJ on reduced costs, including solid amounts down to 1e-41 mol).
- `activeset.js` - the browser port, held to production by
  `tests/test_activeset_port.py` (worst observed 6e-15; trace amounts agree
  bit for bit in log10).
- NIST-vs-Waldner data comparison and the Ellingham direction of the
  boundary pO2, inside the reference file - data against data, no solver.

**What the old "cross-check" actually was.** The oxygen-potential route
(`solidgas/potential.py`) is the *dual* of the same minimisation: its
"margin" is algebraically the KKT reduced cost per mole of oxygen. The
earlier README presented it as an independent check; it is not one, and the
repository no longer claims it is. It survives as the dual identity - useful,
fast, and now correctly named.

## Legacy route

Every number in `results.json` is a **full-composition Gibbs minimisation**
by the legacy reaction-extent SLSQP engine (`solidgas/gibbs.py`): total Gibbs
energy minimised over the whole species set, gases ideal and mixing,
condensed phases pure at unit activity, subject to elemental balance.
Nothing is read off a phase-boundary table, and no point is filled in by the
faster route.

The legacy engine is kept, isolated, for regression - production imports
nothing from it, and a test enforces that. Its known defect (absent phases
at a log-coordinate floor instead of at zero) is what the active-set rewrite
removed; see *A defect, resolved* below. Its gas registry also includes CH4
(without graphite), so its low-temperature carbon-feed rows are conditional
in a way the new page states and excludes.

`solidgas/gibbs.py` works in reaction-extent coordinates. The change of
variables is exact - same objective, same feasible set - but the constant part
of G is removed analytically rather than carried through the arithmetic. That
made the inert case *representable*; the active-set solver is what finally
made it *clean*, because the trace there is a genuine unknown, not a residue.

`solidgas/potential.py` is the oxygen-potential route: each condensed phase's
critical potential as a plain function of temperature, no optimiser. Roughly
200x faster, and - see above - the dual of the reported route, not an
independent check of it. Every row of `results.json` carries a `method`
field, and `tests/test_results.py` refuses any row that does not say
`gibbs_min`.

```bash
# the production pipeline
python3 oracle_tio.py         # 38 points at 80 digits -> reference_results_high_precision.json
python3 oracle_statmech.py    # defect stat-mech reference -> reference_statmech.json
python3 export_activeset.py   # activeset_data.json, straight from the package
python3 export_equations.py   # equations_method.html, typeset from equations_registry.py
python3 build_ti_solver.py    # ti_solver.html, self-contained (both workspaces)
python3 -m pytest tests/ -q   # release gate: no skips, no xfails

# legacy pipeline (kept for regression)
python3 run_all.py            # 24 points -> results.json
python3 build_site.py         # tiox.html
```

## Results

100 mg TiO2, 25 mL, 1 atm. Active-set values from
`reference_results_high_precision.json` (production agrees to 2e-13); the CH4
row is legacy-only, from `results.json`, because the v0.1 gas registry
excludes CH4 rather than carry it without graphite.

**Ti(3+) as a percentage of total titanium**

| feed | 500 C | 700 C | 900 C | 1100 C | 1300 C | 1500 C |
|---|---|---|---|---|---|---|
| CO2/H2 = 1:1 | **0, exact** | **0, exact** | **0, exact** | **0, exact** | **0, exact** | **0, exact** |
| pure H2 | 0.0005 | 0.0168 | 0.1702 | 0.8694 | 2.82 | 6.37 |
| pure CH4 (legacy engine) | 0.1678 | 5.59 | 33.33 | 36.59 | 33.33 | 33.33 |
| sealed N2 | 1.5e-36 | 1.7e-26 | 7.2e-20 | 3.9e-15 | 1.4e-11 | 8.7e-09 |

The RWGS row used to print 1e-17-ish numbers here: those were the legacy
optimiser's log-floor residue, and the KKT reduced costs (+78 to +27 kJ/mol O
across the range) prove the true value is exactly zero. The sealed-N2 row is
the opposite case - the same magnitudes are *real*, the two-phase buffer
under an oxygen-free gas, and the active-set solver resolves them as
first-class unknowns (m(Ti10O19) = 9.66e-42 mol at 500 C, bit-identical
between Python and the browser port in log10).

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

## Defect statistical mechanics (second workspace)

The page carries a second workspace, `Defect stat mech`: given the total
oxygen-vacancy inventory measured in experiment (H2O integration), it computes
where those vacancies sit inside rutile - (110) bridging-O surface, first
subsurface O layer, or bulk. Canonical, not grand-canonical: a 30-minute H2
reduction is not guaranteed to be at equilibrium with the gas, so the total is
taken as a fixed input rather than predicted.

The surface Monte Carlo is independent of particle size. Surface vacancies
interact through a short-range pair table along and across the bridging-O
rows), so the surface isotherm theta_s(mu_V) is sampled by canonical swap
Metropolis on a representative slab with Widom insertion supplying mu_V. The
particle-scale split then equalises mu_V over the *real* geometric site
counts: at 0.9 um the surface layer is 13.6 umol-O/g of capacity (0.05% of
all O sites - the slab's 10% would misplace the inventory by two orders of
magnitude). Subsurface and bulk are non-interacting in this version and use
the exact site-exclusion isotherm. Because the isotherm depends only on
temperature and energetics, one run serves every inventory and particle size;
inventory and geometry changes re-solve instantly in the browser.

At the flagship condition (95 umol-O/g, 600 C, sX energies), the calculated
surface occupation crosses the reported 17-20% reconstruction-onset interval.
The engine does not clip the coverage. It reports the distribution as a
conditional extrapolation of the unreconstructed (1x1) surface Hamiltonian,
because reconstructed-phase energetics have not been implemented. Additional
warnings mark subsurface occupation beyond the dilute regime (a heavily
reduced shell, not isolated point defects) and x approaching the shear-plane
range.

A particle cross-section panel shows the three equilibrium site classes as
concentric regions of a spherical particle. The surface row comes from the
sampled Monte Carlo configuration. The subsurface row uses its calculated
class occupancy. The core uses a continuous fill keyed to the bulk-class
occupancy. A depth plot reports the same three class averages. A toggle
recolours vacancies by CO2 refill probability, and the figure exports as a
standalone SVG.

The two workspaces are bridged. An "Oxygen environment" card in the defect
rail takes the equilibrium workspace's gas, pressure, volume and solid,
then re-solves the charge with the gas-solid Gibbs engine at the defect
temperature, and tiers the verdict for the defect model
(solidgas/statmech.py::phase_validity, ported to the browser and
parity-tested on real equilibria): rutile alone is "stable" with the
smallest reduced-phase KKT cost per mole of oxygen as the margin; rutile
coexisting with a reduced phase is "on the reduction boundary"; an
assemblage without rutile is "invalid". Thermodynamics then predicts a
Magneli phase, and a fixed measured inventory in rutile is metastable. The
results panel reports that condition in red.

A CO2-accessibility layer sits on top of the distribution: first-order
refill per class, k_c = A_c p^m exp(-Ea_c/kT), P_c = 1 - exp(-k_c t), with
the recoverable fraction f_rec as the headline output. The JSON and page
identify the kinetics as placeholders. The prefactors are effective values
because deeper classes are transport-limited. A test preserves that
disclosure. The workspace's main figure sweeps total
inventory: surface follows the sampled isotherm, subsurface and bulk take
the balance, and the accessible curve falls away from the total; a slot for
measured CO2-recovery points overlays experiment on the same axes. The
dielectric card carries the correlation measured in this work: a linear
regression of log10(eps'') on log10(VO) - i.e. the power law
eps'' = 10^-3.07 * VO^1.224 (OriginLab, no weighting, N = 4, Pearson r
0.966, R^2 0.934, fitted over ~95-5000 umol-O/g; below that range the card
flags the value as an extrapolation). Equations (S1)-(S12) are typeset
into the Model panel from equations_registry_statmech.py. The planned
kinetic transport layer will follow the group's coarse-grained kinetic
Monte Carlo framework (Katsoulakis & Vlachos, J. Chem. Phys. 119, 9412
(2003), DOI 10.1063/1.1616513; Katsoulakis, Majda & Vlachos, PNAS 100,
782 (2003)): radial coarse cells over the particle, Arrhenius exchange
rates in detailed balance with the same Gibbs measure - so its t -> inf
limit must reproduce this tab's equilibrium split - and q^2 time
coarse-graining to reach 300 s of real exposure.

Every parameter is data, not code: `rutile_dft.json` holds the layer energies
(Li, Guo & Robertson, J. Phys. Chem. C 119 (2015), sX and GGA presets - the
surface < subsurface < bulk ordering survives the functional), the surface
ordering constraints the pair table is calibrated to (Birschitzky et al., npj
Comput. Mater. 10 (2024): spacings of 1-2 sites suppressed, modal 4-5,
experimental average reconstruction onset ~17%, DFT transition ~20%, local
vacancy-rich coverage up to ~20%, cutoff ~10 A), and diffusion
barriers stored for a later kinetic version (Iddir et al. 2010; Wang et al.
2010). When project DFT numbers arrive they replace the JSON; nothing is
re-fitted in code. Explicit Ti3+ polaron variables, CO2-accessibility
kinetics and the grand-canonical bridge to the equilibrium workspace's
oxygen potential are the next stages.

Verification mirrors the thermodynamic side. `oracle_statmech.py` (mpmath,
50 digits, no shared solver code) certifies the analytic matching grid, exact
finite-lattice canonical ensembles via generating polynomials, and exact
enumerations of interacting micro-lattices - including the chemical potential
that Widom insertion estimates. The JS engine is a line-for-line port pinned
move-for-move: same-seed swap trajectories are compared as integers, class
occupancies and histograms bitwise, mu to 1e-9. The shipped default
configuration reproduces the oracle flagship end-to-end in node.

```bash
python3 oracle_statmech.py     # reference_statmech.json
python3 -m pytest tests/test_statmech.py tests/test_statmech_port.py -q
```

## Layout

```
solidgas/
  shomate.py      NIST-JANAF coefficients; H, S, Cp, mu0, Kp
  waldner.py      Waldner & Eriksson Ti-O assessment: the Magneli series
  activeset.py    PRODUCTION: active-set Gibbs minimisation, KKT decisions
  statmech.py     PRODUCTION: canonical defect lattice model, isotherm + matching
  equilibrium.py  legacy: species, element matrix, G_total, constrained solve
  gibbs.py        legacy: extent-coordinate minimisation (regression only)
  potential.py    oxygen-potential route = the dual/KKT identity
oracle_tio.py           80-digit mpmath oracle -> reference_results_high_precision.json
oracle_statmech.py      50-digit stat-mech oracle -> reference_statmech.json
rutile_dft.json         defect dataset: literature energetics, swapped for project DFT
export_activeset.py     data layer for the browser -> activeset_data.json
activeset.js            the production solver, ported for the browser
statmech.js             the defect engine, ported for the browser (+ _worker, _ui)
ti_solver_template.html + ti_solver_page.js + build_ti_solver.py -> ti_solver.html
run_all.py              legacy: the 24 results.json points
solver.js / page.js     legacy browser engine (tiox.html)
export_thermo.py        legacy data dump (thermo_data.json)
tests/                  release gate: oracle parity, KKT invariants, JS parity,
                        data checks, export-drift checks - no skips, no xfails
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

**<https://robin-yk.github.io/Solid-gas-thermosolver/>** opens the active-set solver
directly (`index.html` is a copy of `ti_solver.html`, made by `build_portal.py`).
The legacy pages below remain reachable from the solver's **Legacy** tab, labelled
with why they are retired; `portal.html` stays in the tree for history but is no
longer the root.

**`oxide_tool.html`: oxide reduction under a reacting gas.** The Gibbs minimisation, in the
browser. Ti&ndash;O and Ce&ndash;O under a CO2/H2 feed; the phase assemblage is an output. Carries
its own verification tab (the same condition solved by the oxygen-potential route, with the
difference shown), the coefficients it is actually running on, the element balance, the solver's
own state, the full derivation, and the solver source printed from the running page. Self-contained:
no external requests.

**`tiox.html`: the oxygen-potential dashboard.** The fast route, swept continuously across
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

### A defect, resolved

An absent phase should have exactly zero moles, and in the legacy logarithmic
extent coordinates it could not: absent phases came back at 1e-26 to 1e-17
percent, and telling that residue apart from a genuine trace needed a second
argument that was already misread once.

The active-set solver removed the defect the way the note above proposed,
with one correction learned on the way: candidates are *not* compared by
total Gibbs energy, because a 1e-34 kJ difference does not survive double
precision - the winner is the candidate whose KKT conditions are clean, which
for this convex program is equivalent and numerically robust. Absent phases
are exactly zero because they are not variables; the sealed-N2 trace at
1e-41 mol is a first-class unknown recovered linearly from the balances; and
`ti_solver.html` shows the reduced cost of every absent phase instead of a
detection-limit apology. The legacy pages (`oxide_tool.html`, `tiox.html`)
retain the old engine and its stated limits; `oxide_tool.html` is no longer
linked from the portal.
