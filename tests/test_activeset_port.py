"""The browser engine against the Python production solver, point by point.

activeset.js is a port, not a reimplementation, and the two run the same
algorithm in the same precision - so the tolerance here is tight. This file
also carries, on the new engine, the two properties the old port could only
xfail: excluded phases exactly zero, and the full assemblage agreeing, not
just the dominant phase. They pass now because absent phases are not
variables at all.

These are release gates: a missing interpreter or artefact fails the test
rather than skipping it, because a gate that can be skipped is not a gate.
"""

import json
import math
import pathlib
import shutil
import subprocess

import pytest

from solidgas import activeset as A

ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA = ROOT / 'data'
WEB = ROOT / 'web'
DOCS = ROOT / 'docs'
HARNESS = ROOT / 'tests' / 'js' / 'activeset_harness.js'

CASES = {
    'rwgs_1_1': {'CO2': 1, 'H2': 1},
    'h2_rich': {'CO2': 1, 'H2': 3},
    'product_mix': {'CO2': 3, 'H2': 3, 'CO': 2, 'H2O': 2},
    'pure_h2': {'H2': 1},
    'pure_n2': {'N2': 1},
    'n2_10ppm_o2': {'N2': 999990, 'O2': 10},
}


@pytest.fixture(scope='module')
def pair():
    assert shutil.which('node') is not None, \
        'node is required to run the release gate'
    for f in ('activeset_data.json', 'reference_results_high_precision.json'):
        assert (DATA / f).exists(), f'{f} missing - run the exporters first'
    r = subprocess.run(['node', str(HARNESS), str(ROOT)],
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, f'harness failed:\n{r.stderr[-2000:]}'
    js = {}
    for row in json.loads(r.stdout):
        js[(row['case'], row.get('T_C_label'))] = row
    py = {}
    doc = json.loads((DATA / 'reference_results_high_precision.json').read_text())
    for case, feed in CASES.items():
        for t in doc['temperatures_C']:
            py[(case, t)] = A.solve(feed, t + 273.15)
    ratio = doc['boundary_validation']['h2_h2o_boundary']
    ratio = float(ratio['H2O_over_H2_critical_ratio'])
    py[('h2_h2o_below_boundary', 900)] = A.solve(
        {'H2': 1, 'H2O': ratio * (1 - 1e-3)}, 1173.15)
    py[('h2_h2o_above_boundary', 900)] = A.solve(
        {'H2': 1, 'H2O': ratio * (1 + 1e-3)}, 1173.15)
    return py, js


def test_every_point_ran_on_both_sides(pair):
    py, js = pair
    for k in py:
        assert k in js, k


def test_method_is_gibbs_min_on_both_sides(pair):
    py, js = pair
    for k in py:
        assert py[k]['method'] == 'gibbs_min', k
        assert js[k]['method'] == 'gibbs_min', k


def test_active_sets_agree(pair):
    py, js = pair
    for k in py:
        assert js[k]['active_condensed_phases'] == \
            py[k]['active_condensed_phases'], k


def test_full_assemblage_agrees_not_just_the_dominant_phase(pair):
    """The property the old port had to xfail. Exact zeros make it exact."""
    py, js = pair
    for k in py:
        pa = {s for s, v in py[k]['solid_moles'].items() if v > 0}
        ja = {s for s, v in js[k]['solid_moles'].items() if v > 0}
        assert pa == ja, (k, sorted(pa), sorted(ja))


def test_absent_phases_are_exactly_zero_in_the_port(pair):
    """The second xfail, retired the same way."""
    py, js = pair
    for k in py:
        for s, v in js[k]['solid_moles'].items():
            if s not in js[k]['active_condensed_phases']:
                assert v == 0.0, (k, s, v)


def test_gas_fractions_agree(pair):
    py, js = pair
    for k in py:
        for g, v in py[k]['gas_fractions'].items():
            if v < 1e-25:
                continue
            d = abs(js[k]['gas_fractions'][g] - v) / v
            assert d < 1e-9, (k, g, d)


def test_trace_solids_agree_in_log10(pair):
    py, js = pair
    for k in py:
        for s, v in py[k]['log10_trace_amounts'].items():
            assert abs(js[k]['log10_trace_amounts'][s] - v) < 1e-9, (k, s)


def test_reduced_costs_agree(pair):
    py, js = pair
    for k in py:
        for s, rc in py[k]['inactive_phase_reduced_costs'].items():
            v = rc['per_mol_O_kJ']
            jv = js[k]['inactive_phase_reduced_costs'][s]['per_mol_O_kJ']
            if v is None:
                assert jv is None, (k, s)
            elif math.isinf(v):
                assert jv is None, (k, s)   # JSON carries -inf as null
            else:
                assert abs(jv - v) < 1e-9, (k, s, jv, v)


def test_gibbs_energies_agree(pair):
    py, js = pair
    for k in py:
        p = py[k]
        d = abs(js[k]['G_total_kJ'] - p['G_total_kJ']) / abs(p['G_total_kJ'])
        assert d < 1e-12, (k, d)
        assert js[k]['common_reference_G_kJ'] == \
            pytest.approx(p['common_reference_G_kJ'], rel=1e-12), k
        if not p['G_rel_below_double_resolution']:
            assert js[k]['G_rel_kJ'] == \
                pytest.approx(p['G_rel_kJ'], rel=1e-9), k


def test_balance_closes_in_the_port(pair):
    py, js = pair
    for k in py:
        for e, v in js[k]['balance_by_element'].items():
            assert abs(v) < 1e-15, (k, e, v)


def test_boundary_flip_reproduces_in_the_port(pair):
    py, js = pair
    assert js[('h2_h2o_below_boundary', 900)]['active_condensed_phases'] == \
        ['TiO2', 'Ti10O19']
    assert js[('h2_h2o_above_boundary', 900)]['active_condensed_phases'] == \
        ['TiO2']


def test_permutation_invariance_in_the_port(pair):
    _, js = pair
    assert js[('permutation_probe', 900)]['identical'] is True


def test_exported_data_is_current():
    """The JSON the browser runs on must match the live package tables."""
    from scripts import export_activeset
    live = export_activeset.build()
    on_disk = json.loads((DATA / 'activeset_data.json').read_text())
    assert live == on_disk, 'run python3 export_activeset.py'


def test_page_is_built_and_fresh():
    """index.html inlines the data, both engines and the UI verbatim."""
    page = DOCS / 'index.html'
    assert page.exists(), 'run python3 scripts/build_site.py'
    html = page.read_text()
    assert (DATA / 'activeset_data.json').read_text().strip() in html, \
        'page is stale against activeset_data.json - run build_site.py'
    sources = (WEB / 'activeset.js', WEB / 'population.js',
               WEB / 'ti_solver_page.js', WEB / 'site_ui.js')
    for src in sources:
        assert src.read_text() in html, \
            f'page is stale against {src.name} - run build_site.py'
    for token in ('/*ASDATA*/', '/*REFDATA*/', '/*DATA*/', '/*ENGINE*/',
                  '/*PAGE*/', '/*POPENGINE*/'):
        assert token not in html
    assert 'gibbs_min' in html
    for phrase in ('CH', 'graphite', 'nitrides', 'not a flow reactor'):
        assert phrase in html, f'scope statement lost: {phrase}'


def test_page_has_no_external_requests():
    html = (DOCS / 'index.html').read_text()
    for bad in ('src="http', "src='http", 'href="http://', '@import',
                'fetch(', 'XMLHttpRequest'):
        assert bad not in html, f'page reaches outside: {bad}'
