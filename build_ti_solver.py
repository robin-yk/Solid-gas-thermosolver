"""Assemble the self-contained Ti-O active-set solver page.

Everything ends up inside one file: the coefficient tables, a slim copy of the
80-digit reference (winner-level values for the Validation tab), the engine,
and the UI. No external requests - a hard requirement, checked below.

    python3 export_activeset.py   # activeset_data.json from the package
    python3 oracle_tio.py         # reference_results_high_precision.json
    python3 build_ti_solver.py    # ti_solver.html
"""

import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))

CASES_FEEDS = {
    'rwgs_1_1': {'CO2': 1, 'H2': 1},
    'h2_rich': {'CO2': 1, 'H2': 3},
    'product_mix': {'CO2': 3, 'H2': 3, 'CO': 2, 'H2O': 2},
    'pure_h2': {'H2': 1},
    'pure_n2': {'N2': 1},
    'n2_10ppm_o2': {'N2': 999990, 'O2': 10},
}


def slim_reference():
    """Winner-level fields only: what the Validation tab diffs against."""
    with open(os.path.join(HERE, 'reference_results_high_precision.json')) as fh:
        doc = json.load(fh)
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


def main():
    with open(os.path.join(HERE, 'ti_solver_template.html')) as fh:
        html = fh.read()
    with open(os.path.join(HERE, 'activeset_data.json')) as fh:
        data = fh.read().strip()
    with open(os.path.join(HERE, 'activeset.js')) as fh:
        engine = fh.read()
    with open(os.path.join(HERE, 'ti_solver_page.js')) as fh:
        page = fh.read()
    with open(os.path.join(HERE, 'equations_method.html')) as fh:
        eqs = fh.read()
    with open(os.path.join(HERE, 'rutile_dft.json')) as fh:
        sm_data = fh.read().strip()
    with open(os.path.join(HERE, 'statmech.js')) as fh:
        sm_engine = fh.read()
    with open(os.path.join(HERE, 'statmech_worker.js')) as fh:
        sm_worker = fh.read()
    with open(os.path.join(HERE, 'statmech_ui.js')) as fh:
        sm_ui = fh.read()
    ref = json.dumps(slim_reference(), separators=(',', ':'))

    for token, body in (('__DATA__', data), ('__REF__', ref),
                        ('__ENGINE__', engine), ('__PAGE__', page),
                        ('__EQS__', eqs), ('__SM_DATA__', sm_data),
                        ('__SM_ENGINE__', sm_engine),
                        ('__SM_WORKER__', sm_worker), ('__SM_UI__', sm_ui)):
        assert token in html, f'{token} missing from the template'
        html = html.replace(token, body)

    out = os.path.join(HERE, 'ti_solver.html')
    with open(out, 'w') as fh:
        fh.write(html)
    # XML namespace URIs in the inlined SVGs are identifiers, not requests.
    checked = html.replace('http://www.w3.org/2000/svg', '') \
                  .replace('http://www.w3.org/1999/xlink', '')
    if 'http://' in checked:
        raise SystemExit('external reference found: http://')
    n_ext = checked.count('https://') - checked.count('https://github.com')
    if n_ext:
        raise SystemExit('external request found beyond the repo link')
    print(f'ti_solver.html written ({os.path.getsize(out) / 1024:.0f} kB)')


if __name__ == '__main__':
    main()
