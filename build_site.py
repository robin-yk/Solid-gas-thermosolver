"""Build index.html for GitHub Pages from the solver's JSON output.

No backend and no build toolchain: the numbers are computed here, embedded in
the page as JSON, and drawn client-side as inline SVG. GitHub Pages serves the
single file as-is.

    python3 run_Ti4O7_RWGS.py       # writes results_rwgs.json
    python3 run_Ti4O7_reducing.py   # writes results_reducing.json
    python3 build_site.py           # writes index.html
"""

import json
import os

import numpy as np

from solidgas import (ACTIVE_TI_PHASES, FORMULAS, SHOMATE, VALID_RANGE,
                      R_KJ, mu0)

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    path = os.path.join(HERE, name)
    if not os.path.exists(path):
        raise SystemExit(f"missing {name} - run the solver scripts first")
    with open(path) as fh:
        return json.load(fh)


def reduction_ladder(temps_C):
    """Stepwise reduction thermodynamics, computed straight from mu0."""
    steps = [('4 TiO2 + H2 -> Ti4O7 + H2O', {'TiO2': -4, 'Ti4O7': 1}),
             ('3 Ti4O7 + H2 -> 4 Ti3O5 + H2O', {'Ti4O7': -3, 'Ti3O5': 4}),
             ('2 Ti3O5 + H2 -> 3 Ti2O3 + H2O', {'Ti3O5': -2, 'Ti2O3': 3})]
    out = []
    for label, rx in steps:
        dG, thr = [], []
        for T_C in temps_C:
            T = T_C + 273.15
            g = mu0('H2O', T) - mu0('H2', T) + sum(v * mu0(k, T) for k, v in rx.items())
            dG.append(round(g, 3))
            thr.append(float(np.exp(-g / (R_KJ * T))))
        out.append({'label': label, 'T_C': list(temps_C), 'dG': dG, 'threshold': thr})
    return out


def disproportionation(temps_C):
    """Ti4O7 + Ti2O3 <-> 2 Ti3O5. Negative dG means Ti3O5 holds together."""
    vals = []
    for T_C in temps_C:
        T = T_C + 273.15
        vals.append(round(2 * mu0('Ti3O5', T) - mu0('Ti4O7', T) - mu0('Ti2O3', T), 3))
    return {'T_C': list(temps_C), 'dG': vals}


def phase_table():
    rows = []
    for p in ACTIVE_TI_PHASES:
        f = FORMULAS[p]
        rows.append({'name': p, 'O_per_Ti': round(f['O'] / f['Ti'], 4),
                     'dHf': SHOMATE[p][2], 'range_K': list(VALID_RANGE[p])})
    return rows


def main():
    rwgs = load('results_rwgs.json')
    red = load('results_reducing.json')
    grid = list(range(500, 1501, 50))
    payload = {
        'conditions': rwgs['conditions'],
        'rwgs': rwgs['rows'],
        'reducing': red['feeds'],
        'ladder': reduction_ladder(grid),
        'dispro': disproportionation(grid),
        'phases': phase_table(),
    }
    with open(os.path.join(HERE, 'site_template.html')) as fh:
        template = fh.read()
    html = template.replace('__DATA__', json.dumps(payload, separators=(',', ':')))
    with open(os.path.join(HERE, 'index.html'), 'w') as fh:
        fh.write(html)
    print(f"index.html written ({len(html) / 1024:.0f} kB)")


if __name__ == '__main__':
    main()
