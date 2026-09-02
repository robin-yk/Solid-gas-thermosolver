# Solid–gas thermosolver

Two calculations used in the accompanying manuscript: gas–solid equilibrium for
the C–H–O–N / Ti–O system, and the isolated surface-vacancy population that sets
the CO formation rate across the reduction series.

The interactive version is at
[robin-yk.github.io/Solid-gas-thermosolver](https://robin-yk.github.io/Solid-gas-thermosolver/).

## What it computes

### Gas–solid equilibrium

Given a gas charge and a temperature, Gibbs-energy minimisation returns the
stable Ti–O phase assemblage, the equilibrium gas composition, the CO2
conversion, and the margin to the nearest reduced phase. Condensed phases are
carried as an active set, so a phase that is not stable has exactly zero moles
rather than a numerical floor, and its distance from stability is reported as a
KKT reduced cost instead. The feed at which rutile stops surviving comes out in
closed form.

Gas data are NIST-JANAF; the titanium oxides through Ti2O3 are the Waldner and
Eriksson CALPHAD assessment.

### Vacancy population model

Titration gives the total oxygen removed per sample; it rises eightfold across
the reduction series while the initial CO rate peaks after 600 °C. The smallest
model that produces that shape counts only isolated vacancies in a surface
active region of capacity nₛ, generated uniformly during reduction and lost by
a thermally activated second-order association. Writing q = n_iso + n_assoc the
loss term cancels from the sum, so q and n_below are closed form and one scalar
equation is integrated; the mass balance is an identity. Three parameters
(nₛ, E_loss, ν_eff) are fitted to five rates, and the residual definition is
reported with every fit because the rates span 86-fold.

## Reproducing the manuscript numbers

```bash
python3 scripts/reproduce_paper.py    # -> paper_outputs/*.csv and site_data.json
python3 scripts/build_site.py         # -> docs/index.html
python3 -m pytest tests/ -q
```

`paper_outputs/` is committed, and a gate fails if regenerating it changes a
byte.

## Verification

Three implementations must agree before a number is quoted: the Python package,
an mpmath oracle at 80 significant digits in `scripts/oracle_tio.py` that shares
only the raw coefficient tables and never imports the solver, and the browser
mirror. The page solves rather than displaying a precomputed grid, because a
grid cannot answer a composition or a parameter set the user actually has — and every
browser engine is paired with a parity gate that holds it to its Python
original: `web/activeset.js` to 4 ulp (the difference between glibc's `exp` and
V8's), `web/population.js` to a relative 1e-6, the difference between two
integrators of the same equation.

The values the manuscript quotes are not recomputed in the page. They ride along
as committed JSON, so what a reader types and what the paper says cannot drift.

`scripts/check_page.py` loads the built page in a real browser and holds every
number it displays against a fresh calculation.

## Layout

```
solidgas/       activeset · vacancy_population · species · shomate · waldner
analysis/       dielectric_correlation — an empirical regression, not thermodynamics
scripts/        reproduce_paper · build_site · oracle_tio · check_page
data/           thermodynamic tables, the reduction series, the 80-digit reference
paper_outputs/  every computed number in the manuscript, as committed CSV
web/            the browser mirrors, a drawing kit, figure modules, one page
docs/           the built page and the method notes
```

See [docs/equilibrium-method.md](docs/equilibrium-method.md) and
[docs/verification.md](docs/verification.md).
