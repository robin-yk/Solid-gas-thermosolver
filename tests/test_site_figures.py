"""The figures on the page, built the way the page builds them.

The browser holds no solver, so there is no parity to check any more. What
is left worth gating is that every figure draws, that its payload carries
the value the Python package produced rather than the nearest grid point,
and that the labels say what the figure is.
"""

import json
import pathlib
import shutil
import subprocess

import pytest

from solidgas import activeset as A
from solidgas import vacancy_inventory as V

ROOT = pathlib.Path(__file__).resolve().parents[1]
HARNESS = ROOT / 'tests' / 'js' / 'site_figures_harness.js'
FEED = {'CO2': 10, 'H2': 10, 'N2': 80}


def _run(T_C):
    assert shutil.which('node') is not None, 'node is required for this gate'
    r = subprocess.run(['node', str(HARNESS), str(ROOT), str(T_C)],
                       capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stderr[-2000:]
    return json.loads(r.stdout)


@pytest.fixture(scope='module')
def js():
    return _run(600)


def test_every_figure_draws(js):
    for name, n in js['bytes'].items():
        assert n > 2000, '%s came out %d bytes - it did not draw' % (name, n)


def test_every_figure_is_the_same_square(js):
    """One geometry for the whole page, in typographic points."""
    assert set(js['viewBoxes'].values()) == {'0 0 252 252'}


def test_the_conversion_payload_is_the_package_number(js):
    r = A.solve(FEED, 873.15)
    assert js['payload']['conversion']['current']['conv_pct'] == pytest.approx(
        r['conversion_CO2_pct'], rel=1e-9)


def test_the_boundary_payload_is_the_closed_form(js):
    b = A.reduction_boundary(873.15)
    assert js['payload']['boundary']['at']['pct'] == pytest.approx(
        b['y_CO2'] * 100.0, rel=1e-9)
    assert js['payload']['boundary']['at']['phase'] == b['phase']


def test_the_margin_payload_is_the_reduced_cost(js):
    r = A.solve(FEED, 873.15)
    costs = {k: v['per_formula_kJ']
             for k, v in r['inactive_phase_reduced_costs'].items()}
    near = min(costs, key=costs.get)
    row = [q for q in js['payload']['margin']['rows']
           if abs(q['T_C'] - 600) < 1e-6][0]
    assert row['margin_kJ'] == pytest.approx(costs[near], rel=1e-9)
    assert row['nearest'] == near


def test_the_capacity_payload_is_the_geometry(js):
    g = V.load()
    cap = js['payload']['capacity']
    assert cap['ceiling'] == pytest.approx(
        V.reconstruction_deficiency(g, d_um=0.9), rel=1e-12)
    assert cap['full_row'] == pytest.approx(
        V.bridging_oxygen_capacity(g, d_um=0.9), rel=1e-12)


def test_every_bar_splits_at_the_surface_ceiling(js):
    cap = js['payload']['capacity']
    for b in cap['bars']:
        assert b['surface_max'] == pytest.approx(cap['ceiling'], rel=1e-12)
        assert b['surface_max'] + b['below_min'] == pytest.approx(
            b['measured'], rel=1e-12)


def test_the_bound_payload_is_the_bound(js):
    g = V.load()
    at = js['payload']['bound']['at']
    want = V.surface_fraction_bound(g, 95.0, d_um=0.9)
    assert at['pct'] == pytest.approx(want['surface_fraction_max'] * 100, rel=1e-9)
    assert at['pct_full'] == pytest.approx(
        want['surface_fraction_max_full_row'] * 100, rel=1e-9)


def test_the_bound_figure_says_it_is_a_bound(js):
    """A reader must not be able to read the curve as an estimate."""
    text = ' '.join(js['texts']['bound'])
    assert 'at most' in text
    assert 'largest possible surface share' in text


def test_the_capacity_figure_names_what_the_bands_are(js):
    text = ' '.join(js['texts']['capacity'])
    assert 'could be surface (at most)' in text
    assert 'must be below the surface' in text
    assert 'whole bridging row' in text


@pytest.mark.parametrize('T_C', [400, 800, 1200, 1500])
def test_the_figures_track_the_temperature_they_are_asked_for(T_C):
    out = _run(T_C)
    r = A.solve(FEED, T_C + 273.15)
    assert out['payload']['conversion']['current']['conv_pct'] == pytest.approx(
        r['conversion_CO2_pct'], rel=1e-9)
    for name, n in out['bytes'].items():
        assert n > 2000, name
