"""Assemble the self-contained Ti-O active-set solver page.

Everything ends up inside one file: the coefficient tables, a slim copy of the
80-digit reference (winner-level values for the Validation tab), the engine,
and the UI. No external requests - a hard requirement, checked below.

    python3 scripts/export_activeset.py
    python3 scripts/oracle_tio.py
    python3 scripts/build_ti_solver.py
"""

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
WEB = ROOT / 'web'
DATA = ROOT / 'data'
DOCS = ROOT / 'docs'

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
    with (DATA / 'reference_results_high_precision.json').open() as fh:
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
    with (WEB / 'ti_solver_template.html').open() as fh:
        html = fh.read()
    with (DATA / 'activeset_data.json').open() as fh:
        data = fh.read().strip()
    with (WEB / 'activeset.js').open() as fh:
        engine = fh.read()
    with (WEB / 'ti_solver_page.js').open() as fh:
        page = fh.read()
    with (WEB / 'generated' / 'equations_method.html').open() as fh:
        eqs = fh.read()
    with (WEB / 'generated' / 'equations_statmech.html').open() as fh:
        sm_eqs = fh.read()
    with (DATA / 'rutile_dft.json').open() as fh:
        sm_data = fh.read().strip()
    with (WEB / 'statmech.js').open() as fh:
        sm_engine = fh.read()
    with (WEB / 'cgmc.js').open() as fh:
        cg_engine = fh.read()
    with (WEB / 'statmech_worker.js').open() as fh:
        sm_worker = fh.read()
    with (WEB / 'figkit.js').open() as fh:
        figkit = fh.read()
    with (WEB / 'figures_defect.js').open() as fh:
        sm_figs = fh.read()
    with (WEB / 'statmech_ui.js').open() as fh:
        sm_ui = fh.read()
    with (DATA / 'redox.json').open() as fh:
        rx_data = fh.read().strip()
    with (WEB / 'redox.js').open() as fh:
        rx_engine = fh.read()
    with (WEB / 'cycling.js').open() as fh:
        cy_engine = fh.read()
    with (WEB / 'figures_redox.js').open() as fh:
        rx_figs = fh.read()
    with (WEB / 'redox_ui.js').open() as fh:
        rx_ui = fh.read()
    with (WEB / 'figures_thermo.js').open() as fh:
        th_figs = fh.read()
    with (WEB / 'thermo_ui.js').open() as fh:
        th_ui = fh.read()
    with (WEB / 'particle.js').open() as fh:
        px_engine = fh.read()
    with (WEB / 'figures_particle.js').open() as fh:
        px_figs = fh.read()
    with (WEB / 'particle_ui.js').open() as fh:
        px_ui = fh.read()
    ref = json.dumps(slim_reference(), separators=(',', ':'))

    for token, body in (('__DATA__', data), ('__REF__', ref),
                        ('__ENGINE__', engine), ('__PAGE__', page),
                        ('__EQS__', eqs), ('__SM_EQS__', sm_eqs),
                        ('__SM_DATA__', sm_data),
                        ('__SM_ENGINE__', sm_engine),
                        ('__CG_ENGINE__', cg_engine),
                        ('__SM_WORKER__', sm_worker),
                        ('__FIGKIT__', figkit), ('__SM_FIGS__', sm_figs),
                        ('__SM_UI__', sm_ui),
                        ('__RX_DATA__', rx_data),
                        ('__RX_ENGINE__', rx_engine),
                        ('__CY_ENGINE__', cy_engine),
                        ('__RX_FIGS__', rx_figs), ('__RX_UI__', rx_ui),
                        ('__TH_FIGS__', th_figs), ('__TH_UI__', th_ui),
                        ('__PX_ENGINE__', px_engine),
                        ('__PX_FIGS__', px_figs), ('__PX_UI__', px_ui)):
        assert token in html, f'{token} missing from the template'
        html = html.replace(token, body)

    out = DOCS / 'ti_solver.html'
    with out.open('w') as fh:
        fh.write(html)
    # XML namespace URIs in the inlined SVGs are identifiers, not requests.
    checked = html.replace('http://www.w3.org/2000/svg', '') \
                  .replace('http://www.w3.org/1999/xlink', '')
    if 'http://' in checked:
        raise SystemExit('external reference found: http://')
    n_ext = checked.count('https://') - checked.count('https://github.com')
    if n_ext:
        raise SystemExit('external request found beyond the repo link')
    print(f'{out.relative_to(ROOT)} written ({out.stat().st_size / 1024:.0f} kB)')


if __name__ == '__main__':
    main()
