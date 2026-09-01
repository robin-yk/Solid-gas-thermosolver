# Solid–gas thermosolver

Two calculations used in the accompanying manuscript: gas–solid equilibrium for
the C–H–O–N / Ti–O system, and a geometric comparison between the measured
oxygen-removal inventory and the oxygen capacity of the rutile surface. It does
not assign the below-surface inventory among subsurface point defects, bulk
defects, or extended defects.

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

### Vacancy inventory against surface capacity

A titration measures one number, the total oxygen removed. It does not say where
that oxygen came from. From the rutile lattice constants and the exposed area,
this returns how much oxygen the (110) bridging row could hold, and therefore a
**bound**: the largest share of the measured inventory that could be surface
oxygen, and the smallest share that must lie below it.

For a 0.9 µm particle the bridging row holds 13.6 µmol-O g⁻¹ in total and
6.8 µmol-O g⁻¹ once reconstructed. A measured removal of 95 µmol-O g⁻¹ is
therefore at most 7.1 % surface oxygen and at least 92.9 % below it.

No defect formation energy enters this calculation. Literature DFT does indicate
a strong depth dependence of vacancy stability; those values are tabulated in
`data/literature_dft.json` and are deliberately not used to assign a depth
distribution, because no measurement in this work resolves depth.

## Reproducing the manuscript numbers

```bash
python3 scripts/reproduce_paper.py    # -> paper_outputs/*.csv and site_data.json
python3 scripts/build_site.py         # -> docs/index.html
python3 -m pytest tests/ -q
```

`paper_outputs/` is committed, and a gate fails if regenerating it changes a
byte.

## Verification

Two independent implementations must agree before a number is quoted: the Python
package, and an mpmath oracle at 80 significant digits in `scripts/oracle_tio.py`
that shares only the raw coefficient tables and never imports the solver. The
browser does not compute anything — it displays the JSON the package wrote — so
there is no third implementation to keep in step.

`scripts/check_page.py` loads the built page in a real browser and holds every
number it displays against a fresh calculation.

## Layout

```
solidgas/       activeset · vacancy_inventory · species · shomate · waldner
analysis/       dielectric_correlation — an empirical regression, not thermodynamics
scripts/        reproduce_paper · build_site · oracle_tio · check_page
data/           thermodynamic tables, rutile geometry, tabulated literature DFT
paper_outputs/  every computed number in the manuscript, as committed CSV
web/            the display layer: a drawing kit, two figure modules, one page
docs/           the built page and the method notes
```

See [docs/equilibrium-method.md](docs/equilibrium-method.md),
[docs/surface-capacity.md](docs/surface-capacity.md) and
[docs/verification.md](docs/verification.md).
