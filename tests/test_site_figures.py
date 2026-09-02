"""The figures on the page, built the way the page builds them.

The browser holds no solver, so there is no parity to check any more. What
is left worth gating is that every figure draws, that its payload carries
the value the Python package produced rather than the nearest grid point,
and that the labels say what the figure is.
"""

import json
import math
import pathlib
import shutil
import subprocess

import pytest

from solidgas import activeset as A
from solidgas import vacancy_population as VP

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
    assert set(js['viewBoxes'].values()) == {'0 0 360 360'}


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


@pytest.mark.parametrize('T_C', [400, 800, 1200, 1500])
def test_the_figures_track_the_temperature_they_are_asked_for(T_C):
    out = _run(T_C)
    r = A.solve(FEED, T_C + 273.15)
    assert out['payload']['conversion']['current']['conv_pct'] == pytest.approx(
        r['conversion_CO2_pct'], rel=1e-9)
    for name, n in out['bytes'].items():
        assert n > 2000, name


def test_the_population_payload_is_the_package_partition(js):
    """The scatter must plot what the model produced, not a rounded copy."""
    rows = VP.partition(VP.load_series(), VP.REPORTED['ns'],
                        VP.REPORTED['E_loss'], 10.0 ** VP.REPORTED['log10_nu'])
    want = {q['sample']: q for q in rows}
    got = js['payload']['population']['points']
    assert [q['sample'] for q in got] == [q['sample'] for q in rows]
    for q in got:
        w = want[q['sample']]
        assert q['total'] == pytest.approx(w['nV_total_umol_g'], rel=1e-12)
        for a, b in (('iso', 'n_iso_umol_g'), ('assoc', 'n_assoc_umol_g'),
                     ('below', 'n_below_umol_g')):
            assert q[a] == pytest.approx(w[b], rel=1e-6, abs=1e-9)


def test_the_population_figure_plots_every_assignment(js):
    """Three shapes, one per population, and nothing silently dropped.

    A log axis cannot show a non-positive value, so a figure that quietly
    skipped one would understate the partition. Every assignment in the
    payload is positive and every one of them is on the plot; the extra
    mark of each shape is its legend key.
    """
    pts = js['payload']['population']['points']
    for q in pts:
        assert min(q['iso'], q['assoc'], q['below']) > 0
    m = js['marks']['population']
    assert m == {'circle': len(pts) + 1, 'square': len(pts) + 1,
                 'triangle': len(pts) + 1}


def test_the_population_axes_bracket_the_data(js):
    """Decade-bracketed axes, so no computed point can fall outside them.

    The page recomputes this figure from whatever parameters the reader
    sets, so the range cannot be hand-written: it is derived from the data
    every time. The gate is that the labelled decades are exactly the ones
    the data spans - no missing edge, no empty padding.
    """
    def span(lo, hi):
        return [10 ** e for e in range(math.floor(math.log10(lo)),
                                       math.ceil(math.log10(hi)) + 1)]

    pts = js['payload']['population']['points']
    xs = [q['total'] for q in pts]
    ys = [v for q in pts for v in (q['iso'], q['assoc'], q['below'])]
    want_x = span(min(xs), max(xs))
    want_y = span(min(ys), max(ys))
    assert want_x[0] <= min(xs) and want_x[-1] >= max(xs)
    assert want_y[0] <= min(ys) and want_y[-1] >= max(ys)

    labels = [t for t in js['texts']['population'] if t.startswith('10')]
    assert len(labels) == len(want_x) + len(want_y)
    for v in set(want_x) | set(want_y):
        e = round(math.log10(v))
        sup = str(e).replace('-', '\u2212').translate(
            str.maketrans('0123456789\u2212', '\u2070\u00b9\u00b2\u00b3'
                          '\u2074\u2075\u2076\u2077\u2078\u2079\u207b'))
        assert '10' + sup in labels, 'decade 1e%d has no tick label' % e


def test_the_population_figure_names_the_three_populations(js):
    text = ' '.join(js['texts']['population'])
    assert 'Isolated accessible' in text
    assert 'Associated surface-region' in text
    assert 'Below active region' in text
    assert 'Total Vₒ' in text and 'Assigned Vₒ' in text
