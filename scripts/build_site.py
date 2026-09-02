"""Assemble the site: one self-contained page, no external requests.

    python3 scripts/export_activeset.py   # the coefficient tables
    python3 scripts/oracle_tio.py         # the 80-digit reference
    python3 scripts/reproduce_paper.py    # the committed manuscript values
    python3 scripts/build_site.py         # -> docs/index.html

Both workspaces solve in the page. The equilibrium half runs the active-set
solver mirrored in web/activeset.js, and the surface half runs the dozen
lines of lattice arithmetic mirrored in web/vacancy.js; parity gates hold
each against its Python original. The committed manuscript values ride
along as site_data.json so the figures and the paper cannot drift.
"""

import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WEB = os.path.join(ROOT, 'web')
DATA = os.path.join(ROOT, 'data')
OUT = os.path.join(ROOT, 'docs', 'index.html')
SITE = os.path.join(ROOT, 'paper_outputs', 'site_data.json')

CASES_FEEDS = {
    'rwgs_1_1': {'CO2': 1, 'H2': 1},
    'h2_rich': {'CO2': 1, 'H2': 3},
    'product_mix': {'CO2': 3, 'H2': 3, 'CO': 2, 'H2O': 2},
    'pure_h2': {'H2': 1},
    'pure_n2': {'N2': 1},
    'n2_10ppm_o2': {'N2': 999990, 'O2': 10},
}

PARTS = [
    ('/*CSS*/', os.path.join(WEB, 'site.css')),
    ('/*FIGKIT*/', os.path.join(WEB, 'figkit.js')),
    ('/*ENGINE*/', os.path.join(WEB, 'activeset.js')),
    ('/*VACANCY*/', os.path.join(WEB, 'vacancy.js')),
    ('/*THERMO*/', os.path.join(WEB, 'figures_thermo.js')),
    ('/*SURFACE*/', os.path.join(WEB, 'figures_surface.js')),
    ('/*PAGE*/', os.path.join(WEB, 'ti_solver_page.js')),
    ('/*THERMOUI*/', os.path.join(WEB, 'thermo_ui.js')),
    ('/*UI*/', os.path.join(WEB, 'site_ui.js')),
]


def read(path):
    with open(path) as fh:
        return fh.read()


def slim_reference():
    """Winner-level fields only: what the Validation tab diffs against."""
    doc = json.loads(read(os.path.join(DATA,
                                       'reference_results_high_precision.json')))
    rows = []
    for r in doc['rows']:
        if r['case'] not in CASES_FEEDS:
            continue
        rows.append({
            'case': r['case'], 'T_C': r['T_C'],
            'active': r['active_condensed_phases'],
            'gas_fractions': r['gas_fractions'],
            'r_per_mol_O': {s: rc['per_mol_O_kJ'] for s, rc
                            in r['inactive_phase_reduced_costs'].items()},
            'log10_traces': r['log10_trace_amounts'],
            'reduced_pct': r['reduced_pct'],
        })
    return {'precision': doc['precision'], 'method': doc['method'],
            'cases_feeds': CASES_FEEDS, 'rows': rows,
            'h2_h2o_boundary': doc['boundary_validation']['h2_h2o_boundary']}


def build():
    if not os.path.exists(SITE):
        raise SystemExit('run scripts/reproduce_paper.py first')
    html = read(os.path.join(WEB, 'template.html'))
    html = html.replace('/*ASDATA*/',
                        read(os.path.join(DATA, 'activeset_data.json')).strip())
    html = html.replace('/*REFDATA*/',
                        json.dumps(slim_reference(), separators=(',', ':')))
    html = html.replace('/*DATA*/',
                        json.dumps(json.loads(read(SITE)),
                                   separators=(',', ':')))
    for token, path in PARTS:
        html = html.replace(token, read(path))
    for token, _ in PARTS + [('/*DATA*/', None), ('/*ASDATA*/', None),
                             ('/*REFDATA*/', None)]:
        if token in html:
            raise SystemExit('unsubstituted token left in the page: ' + token)
    if '<script src=' in html or 'href="http' in html.replace(
            'href="https://github.com', ''):
        raise SystemExit('the page must not reach outside itself')
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, 'w') as fh:
        fh.write(html)
    return OUT, len(html)


if __name__ == '__main__':
    path, size = build()
    print('%s written (%d kB)' % (os.path.relpath(path, ROOT), size // 1024))
