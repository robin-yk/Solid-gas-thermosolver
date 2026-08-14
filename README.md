# solid-gas-thermosolver

Equilibrium solver for a gas phase in contact with a reducible oxide, built to
answer one question: **can a reverse water-gas-shift feed reduce rutile TiO2 to
a Magneli phase?**

The answer is no across 500-1500 C, but the margin is not uniform, and the
place where it is thinnest is not yet fully tested. See *Where this is soft*.

## Two routes to the same answer

**Full Gibbs minimisation** (`solidgas/equilibrium.py`). Total Gibbs energy is
minimised over the mole vector, gases ideal and mixing, condensed phases pure
at unit activity, subject to elemental balance on C, H, O and Ti. Multiple
starts per temperature, because each solid sits in its own basin.

**Oxygen potential** (`solidgas/potential.py`). Titanium does not enter the gas
below about 2000 K, so the solid exchanges nothing with the gas but oxygen.
That splits the problem: equilibrate the C-H-O gas alone, read off its oxygen
potential, and compare against each condensed phase's critical potential, which
is a plain function of temperature with no optimisation in it at all.

The second route is **238x faster** (1.35 s versus 321 s for an 80-point sweep)
and **more accurate** - it holds the RWGS quotient to 4e-5 against Kp, where the
full nine-species minimisation manages 1e-3. Both are kept: the fast route is
the one to use, and `test_potential_route_matches_the_full_minimisation`
pins them to each other so the fast one cannot drift.

The speed is not the interesting part. Line compounds are present or absent,
not continuous unknowns, so making them variables forces the optimiser to walk
an absent phase's mole number down to its bound through a log singularity that
never needed to exist. Removing the solids as unknowns removes the reason for
multistart.

## Results

100 mg TiO2, 25 mL, 1 atm, 500-1500 C.

| feed | Ti(3+) at 900 C | Ti(3+) at 1500 C | deepest phase |
|---|---|---|---|
| CO2/H2 = 1:1 | 0.000% | 0.000% | none - rutile throughout |
| pure H2 | 0.17% | 6.4% | Ti10O19 |
| pure CH4 | 33.3% | 33.0% | Ti5O9 / Ti6O11 |
| sealed N2 | 0.000% | 0.000% | none |

Both solvers were run on the pure-H2 case and agree to the printed digit.
An apparent disagreement at first turned out to be only that one phase list
carried Ti20O39 and the other did not, which is worth more than the agreement:
truncating the series at Ti10O19 gives 6.4% at 1500 C and including Ti20O39
gives 9.8%, and that spread is the line-compound model reaching its limit
rather than a numerical problem.

Sealed inert gas does nothing, and not because the temperature is too low: the
oxide has nowhere to put the oxygen. At 1500 C the boundary sits at
pO2 = 1.6e-10 atm, and 25 mL at that pressure holds 2.7e-14 mol of O2 while
converting even 1% of the charge would release 6.3e-7 mol.

Under RWGS the gas equilibrates to ordinary shift behaviour, 25% to 66% CO2
conversion, and the solid never enters. Kp from standard potentials and Q from
the converged mole fractions agree at every temperature.

**Ti3O5 and Ti2O3 never appear.** They were added to test whether the phase set
was limiting the answer. It was not: under CH4 the solid parks on the
TiO2/Ti4O7 two-phase buffer, and the gas there is 12 to 74 times too oxidising
to carry the next step. The 92% Ti4O7 reached at 1500 C is a buffer, not
saturation.

### How far RWGS sits from reduction

Margin in kJ per mole of oxygen removed, above the phase that forms most
easily out of untouched rutile. Positive means nothing reduces.

| T (C) | 500 | 700 | 900 | 1100 | 1300 | 1500 |
|---|---|---|---|---|---|---|
| vs Ti4O7 only (NIST) | +69.2 | +47.1 | +45.9 | +35.4 | +24.4 | +22.7 |
| vs Ti4O7 only (Waldner) | +85.0 | +71.0 | +64.3 | +57.2 | +49.5 | +41.2 |
| **vs the full series (Waldner)** | **+74.7** | **+58.2** | **+48.9** | **+39.3** | **+29.1** | **+18.3** |

The last row is the answer. It is measured against Ti20O39, the highest-n
member assessed and the easiest phase to form; against the more conventional
Ti10O19 the 1500 C figure is +27.3.

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
  potential.py    oxygen-potential route: gas-only solve + critical potentials
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

## The interactive page

**<https://robin-yk.github.io/Solid-gas-thermosolver/>**

Pick an atmosphere and conditions; the page equilibrates the C-H-O gas, reads
off its oxygen potential, and compares it against every member of the series.
Nothing is a stored answer - the whole calculation runs in the browser, which
is only possible because the oxygen-potential route has no optimiser in it.

```bash
python3 export_thermo.py   # thermo_data.json, straight from the package
python3 build_site.py      # inlines solver.js + page.js + the data into index.html
```

**Series colour is two different jobs, and merging them fails by measurement.**
Phase boundaries are ordered by O/Ti, so swapping them changes the meaning:
they take an ordinal one-hue ramp stepped off `--accent`. Gas species are
nominal - CO2 is not "more" than H2 - so identity needs separate hues. Five
steps of one hue score CVD ΔE 5.4 and a normal-vision ΔE 5.6 against a floor
of 15; the separate-hue set clears every check in both modes. Both palettes
were run through the validator rather than eyeballed.

`--danger` and `--safe` are a reserved status channel and never carry a
series. Reduction is the outcome under test, so it wears `--danger`; rutile
surviving is the safe reading.

The coefficients are exported rather than retyped, and `tests/test_export.py`
fails if `thermo_data.json` or `index.html` goes stale against the package -
a drifted page would keep working and quietly answer with old numbers.

`solver.js` reimplements the thermodynamics in the browser and agrees with the
Python package to the printed digit on every case. It is in one respect
better: the gas equilibrium is solved by nested bisection rather than SLSQP,
which holds the RWGS quotient to machine precision against Kp where SLSQP
manages 4e-5.

The CH4 case carries a visible caveat in the page itself, since graphite is
still missing.
