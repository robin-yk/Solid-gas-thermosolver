"""web/population.js against solidgas/vacancy_population.py.

Two different integrators of one equation: Python hands the scalar ODE to
Radau, the browser uses Strang splitting with both sub-flows in closed
form. The budget is a relative tolerance rather than the ulp, because
pretending two schemes agree bit for bit would be a fiction.
"""

import json
import pathlib
import shutil
import subprocess

import pytest

from solidgas import vacancy_population as VP

ROOT = pathlib.Path(__file__).resolve().parents[1]
HARNESS = ROOT / 'tests' / 'js' / 'population_harness.js'
TOL = 1e-6

CASES = [
    ('reported', VP.REPORTED['ns'], VP.REPORTED['E_loss'],
     10.0 ** VP.REPORTED['log10_nu']),
    ('raw_optimum', 15.19, 1.6365, 10.0 ** 4.1834),
    ('nu_1e13', 14.92, 3.760, 1e13),
    ('capacity_at_the_top_bound', 67.8, 2.18, 5.5e6),
    ('cold_surface', 15.3, 4.0, 1e13),
]


@pytest.fixture(scope='module')
def js():
    assert shutil.which('node') is not None, 'node is required for this gate'
    spec = [{'label': c[0], 'ns': c[1], 'E_loss': c[2], 'nu_eff': c[3]}
            for c in CASES]
    r = subprocess.run(['node', str(HARNESS), str(ROOT), json.dumps(spec)],
                       capture_output=True, text=True, timeout=180)
    assert r.returncode == 0, r.stderr[-2000:]
    return json.loads(r.stdout)


@pytest.fixture(scope='module')
def series():
    return VP.load_series()


def test_the_two_mirrors_report_the_same_reported_parameters(js):
    assert js['reported']['ns'] == pytest.approx(VP.REPORTED['ns'], rel=1e-12)
    assert js['reported']['E_loss'] == pytest.approx(VP.REPORTED['E_loss'], rel=1e-12)
    assert js['reported']['log10_nu'] == pytest.approx(VP.REPORTED['log10_nu'],
                                                       rel=1e-12)


def test_every_population_matches_across_the_parameter_grid(js, series):
    worst = 0.0
    for case in js['cases']:
        py = {q['sample']: q for q in
              VP.partition(series, case['ns'], case['E_loss'], case['nu_eff'])}
        for row in case['rows']:
            q = py[row['sample']]
            for a, b in (('n_iso', 'n_iso_umol_g'),
                         ('n_assoc', 'n_assoc_umol_g'),
                         ('n_below', 'n_below_umol_g')):
                want = q[b]
                # a near-cancelling population (n_assoc is a difference of
                # two much larger numbers) has no meaningful relative
                # error, so the scale is the sample's own inventory
                scale = max(abs(want), 1e-3 * q['nV_total_umol_g'])
                rel = abs(row[a] - want) / scale
                worst = max(worst, rel)
                assert rel < TOL, (case['label'], row['sample'], a, rel)
    assert worst < TOL, worst


def test_the_closure_is_exact_on_both_sides(js):
    for case in js['cases']:
        for row in case['rows']:
            assert abs(row['closure']) <= VP.CLOSURE_TOL, (case['label'], row)


def test_the_fitted_rates_match(js, series):
    for case in js['cases']:
        py = {q['sample']: q['r_CO_fitted'] for q in
              VP.partition(series, case['ns'], case['E_loss'], case['nu_eff'])}
        for row in case['rows']:
            assert row['r_CO_fitted'] == pytest.approx(py[row['sample']],
                                                       rel=TOL), case['label']


def test_all_three_residual_definitions_match(js, series):
    for case in js['cases']:
        rows = VP.partition(series, case['ns'], case['E_loss'], case['nu_eff'])
        for kind in VP.RESIDUALS:
            want = sum(v * v for v in VP.residuals(rows, kind))
            assert case['objective'][kind] == pytest.approx(want, rel=1e-5), \
                (case['label'], kind)


def test_the_browser_fit_finds_the_reported_parameters(js):
    """The page can refit, and it lands where Python lands."""
    f = js['fit']
    assert f['ns'] == pytest.approx(VP.REPORTED['ns'], rel=3e-3)
    assert f['E_loss'] == pytest.approx(VP.REPORTED['E_loss'], rel=3e-3)
    assert f['log10_nu'] == pytest.approx(VP.REPORTED['log10_nu'], rel=3e-3)


def test_the_mirror_refuses_what_the_module_refuses(js):
    assert js['rejects_zero_inventory'] is True
    assert js['rejects_zero_capacity'] is True
