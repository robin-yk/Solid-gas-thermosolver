# Legacy solvers

## Reaction-extent Gibbs solver

`solidgas/gibbs.py` minimizes total Gibbs energy in reaction-extent coordinates with SciPy SLSQP. It remains for regression and for the older Ti-O and Ce-O browser tool. Production active-set code does not import it.

The legacy logarithmic coordinates cannot represent an absent phase as exact zero. Inactive phases can return small floor residues that look like trace reduction. The active-set rewrite removed this ambiguity by excluding inactive phases from the variable set and reporting their KKT reduced costs.

## Oxygen-potential route

`solidgas/potential.py` evaluates the phase ladder directly from oxygen chemical potential. It is the dual form of the active-set stability test, not an independent numerical oracle. It remains useful for fast sweeps and visual phase-boundary plots.

## Legacy browser pages

- `docs/oxide_tool.html`: first-generation browser Gibbs solver
- `docs/tiox.html`: oxygen-potential dashboard with committed regression results

Their sources are in `web/`. Their build commands and reference generators are in `scripts/`.

## Browser port

The old browser Gibbs engine replaces SciPy SLSQP with projected gradient steps, Barzilai-Borwein scaling, and backtracking in transformed extent coordinates. Gas Shomate constants are folded during export so the JavaScript evaluator reproduces the Python standard-state convention. `tests/test_oxide_port.py` compares the port with Python across both oxide systems and six temperatures.

## Regression data

`data/results.json` contains 24 legacy Gibbs-minimization points and records `method: gibbs_min` on every row. `tests/test_results.py` ensures the table and generated dashboard embed the same committed data.
