"""Render the two figures the redox model is for, at print size.

The physics stays in solidgas/redox.py: this computes the rows there and
writes them out, then hands them to web/figures_redox.js, which draws them
through the shared kit and knows nothing about the model. A figure cannot
disagree with the engine because it never evaluates the engine's formulae.
"""

import json
import math
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from solidgas import redox as R                       # noqa: E402

OUT = os.path.join(ROOT, 'docs', 'figures')
DATA = os.path.join(ROOT, 'data', 'redox_figure_data.json')

K_SHOWN = 0.5           # a carrier well inside the balanced band
E_RUTILE = 1844.0       # the 0.9 um layer-resolved enrichment


def build(p):
    st = p['structure']
    sol = st['x_solubility'] / st['oxygen_per_formula']
    th_max = st.get('theta_surface_max', 0.5)

    curves = R.crossing_curves(p, n=201)
    k = R.rate_constants(p)
    scale = k['k_ox_s1']
    ratio = R.from_ratio(p, K_SHOWN, window_n=0)
    # the drawn curves are in units of k_ox at the shown ratio, so the
    # crossing on the page is the crossing the solver reports
    rows = []
    for i in range(201):
        th = i / 200.0
        rows.append({'theta': th,
                     'r_red': K_SHOWN * (1.0 - th),
                     'r_ox': th})

    # The second-phase boundary: K at which the average composition first
    # leaves the point-defect range. It rises with enrichment and meets the
    # face-ceiling line at E = theta_surface_max / theta_sol, above which it
    # bounds nothing - so it is emitted only that far.
    cross_e = th_max / sol
    band = []
    for j in range(241):
        e = 10.0 ** (4.0 * j / 240.0)
        if e > cross_e:
            break
        t = e * sol
        band.append({'E': e, 'K': t / (1.0 - t)})
    band.append({'E': cross_e,
                 'K': th_max / (1.0 - th_max)})

    return {
        'generated_by': 'export_redox_figures.py',
        'crossing': {
            'K': K_SHOWN,
            'theta_star': ratio['theta'],
            'rate_over_k_ox': ratio['rate_s1'],
            'regime': ratio['regime'],
            'rows': rows,
            # the exported crossing is the uniform reading, E = 1, so the
            # particle's single-phase limit lands on the face axis at
            # theta_sol itself - which is the figure's own comment on the
            # assumption it is drawn under
            'theta_sol': sol,
            'enrichment': 1.0,
            'theta_face_at_solubility': min(1.0, sol),
            'theta_surface_max': th_max,
            'descriptor_k_ox_s1': scale},
        'phase': {
            'theta_sol': sol, 'theta_surface_max': th_max,
            'K_lo': 1e-4, 'K_hi': 1e4, 'E_lo': 1.0, 'E_hi': 1e4,
            'K_no_steady_state': R.usable_window(p, 1.0)[
                'K_no_steady_state'],
            'crossover_E': th_max / sol,
            'second_phase': band,
            'points': [
                {'E': 1.0, 'K': R.usable_window(p, 1.0)['K_max'],
                 'label': 'uniform'},
                {'E': E_RUTILE, 'K': R.usable_window(p, E_RUTILE)['K_max'],
                 'label': 'layer-resolved'}]},
        'sweep': [{'K': r['K'], 'theta': r['theta'],
                   'rate': r['rate_over_k_ox']}
                  for r in R.ratio_sweep(p, n=161)]}


def main():
    p = R.load_params()
    doc = build(p)
    with open(DATA, 'w') as fh:
        json.dump(doc, fh, indent=1)
    os.makedirs(OUT, exist_ok=True)
    r = subprocess.run(
        ['node', os.path.join(ROOT, 'scripts', 'render_redox_figures.js'),
         ROOT, DATA, OUT],
        capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit('render failed:\n' + r.stderr[-2000:])
    print(r.stdout.strip())
    print('redox_figure_data.json written')


if __name__ == '__main__':
    main()
