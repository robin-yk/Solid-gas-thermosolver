"""The two redox figures, and the promise that they cannot lie.

The physics is computed in solidgas/redox.py and written to
data/redox_figure_data.json; web/figures_redox.js only draws. So the gates
here are of two kinds: the house style the kit is supposed to enforce, and
the agreement between what the exporter wrote and what the engine says -
if a boundary on the page ever disagreed with usable_window(), one of the
two would have to have grown its own copy of the algebra.

Nothing here may skip: a figure that will not build is a failure.
"""

import json
import math
import pathlib
import re
import shutil
import subprocess
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

from solidgas import redox as R                        # noqa: E402
import figure_registry as REG                          # noqa: E402

SQUARE = 'width="3.5in" height="3.5in" viewBox="0 0 252 252"'
ATTR = re.compile(r'([\w-]+)="([^"]*)"')
DATA = ROOT / 'data' / 'redox_figure_data.json'
OUT = ROOT / 'docs' / 'figures'


@pytest.fixture(scope='module')
def built():
    assert shutil.which('node') is not None, 'node is required'
    r = subprocess.run([sys.executable,
                        str(ROOT / 'scripts' / 'export_redox_figures.py')],
                       capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, r.stderr[-2000:]
    doc = json.loads(DATA.read_text())
    svg = {k: (OUT / ('redox_%s.svg' % k)).read_text()
           for k in ('crossing', 'phase')}
    return doc, svg


def _attrs(s):
    return dict(ATTR.findall(s))


def _boxes(svg):
    out = []
    for a in (_attrs(m.group(1))
              for m in re.finditer(r'<rect ([^>]*fill="none"[^>]*)/>', svg)):
        if a.get('stroke') != '#000000':
            continue
        if abs(float(a.get('stroke-width', 0)) - 0.9) > 1e-9:
            continue
        x, y = float(a['x']), float(a['y'])
        out.append((x, y, x + float(a['width']), y + float(a['height'])))
    return out


# ------------------------------------------------------------ house style

def test_both_are_exactly_three_and_a_half_inches_square(built):
    _, svg = built
    for fid, s in svg.items():
        assert SQUARE in s, fid


def test_each_has_one_closed_square_axis_box(built):
    _, svg = built
    for fid, s in svg.items():
        boxes = _boxes(s)
        assert len(boxes) == 1, (fid, len(boxes))
        x0, y0, x1, y1 = boxes[0]
        assert abs((x1 - x0) - (y1 - y0)) < 1e-9, fid


def test_ticks_point_inward_and_labels_stay_bottom_left(built):
    _, svg = built
    tol = 0.01
    for fid, s in svg.items():
        b = _boxes(s)[0]
        for a in (_attrs(m.group(1))
                  for m in re.finditer(r'<line ([^>]*)/>', s)):
            if abs(float(a.get('stroke-width', 0)) - 0.9) > 1e-9:
                continue
            xs = [float(a['x1']), float(a['x2'])]
            ys = [float(a['y1']), float(a['y2'])]
            assert (b[0] - tol <= min(xs) and max(xs) <= b[2] + tol
                    and b[1] - tol <= min(ys) and max(ys) <= b[3] + tol), \
                (fid, xs, ys)
        for m in re.finditer(r'<text ([^>]*)>([^<]*)</text>', s):
            a, body = _attrs(m.group(1)), m.group(2)
            if 'rotate' in (a.get('transform') or ''):
                continue
            if not re.fullmatch(r'[-−0-9.,⁰-₟ ]+', body or 'x'):
                continue
            assert not (float(a['y']) < b[1] - 1), (fid, body)
            assert not (float(a['x']) > b[2] + 1), (fid, body)


def test_type_and_palette_are_the_declared_ones(built):
    _, svg = built
    hues = {h for _, h, _ in REG.ROLES}
    neutral = {'#FFFFFF', '#000000'}

    def rgb(c):
        n = int(c[1:], 16)
        return ((n >> 16) & 255, (n >> 8) & 255, n & 255)

    def tint_of(c, hue):
        cr, hr = rgb(c), rgb(hue)
        best = max(zip(cr, hr), key=lambda ab: 255 - ab[1])
        if best[1] == 255:
            return all(a == 255 for a in cr)
        f = (255 - best[0]) / (255 - best[1])
        if not 0.0 < f <= 1.0:
            return False
        return all(abs(255 - (255 - b) * f - a) <= 2.0
                   for a, b in zip(cr, hr))

    for fid, s in svg.items():
        sizes = {float(v) for v in re.findall(r'font-size="([\d.]+)"', s)}
        assert sizes <= {8.0, 9.0}, (fid, sorted(sizes))
        assert min(sizes) >= REG.SPEC['type_floor_pt'], fid
        for c in set(re.findall(r'#[0-9A-Fa-f]{6}', s)):
            if c.upper() in neutral or c in hues:
                continue
            assert any(tint_of(c, h) for h in hues), (fid, c)


def test_no_grid_and_no_title_inside_either_panel(built):
    _, svg = built
    # Bold marks a region name on the phase map and nothing else; a
    # figure title inside the panel would look the same, so the count is
    # gated rather than the mere absence.
    bold = re.findall(r'<text[^>]*font-weight="bold"[^>]*>([^<]*)</text>',
                      svg['phase'])
    assert bold == ['usable'], bold
    assert 'font-weight="bold"' not in svg['crossing']
    for fid, s in svg.items():
        b = _boxes(s)[0]
        for a in (_attrs(m.group(1))
                  for m in re.finditer(r'<line ([^>]*)/>', s)):
            if abs(float(a.get('stroke-width', 0)) - 0.9) > 1e-9:
                continue
            span = abs(float(a['x2']) - float(a['x1']))
            assert span < (b[2] - b[0]) - 1, (fid, 'grid line')


# --------------------------------------- the drawing agrees with the engine

def test_the_crossing_drawn_is_the_crossing_solved(built):
    """The marked point is where the two drawn curves actually meet."""
    doc, _ = built
    d = doc['crossing']
    p = R.load_params()
    assert d['theta_star'] == pytest.approx(
        R.theta_of_ratio(d['K']), abs=1e-15)
    # the two columns cross once, and there
    rows = d['rows']
    sign = [r['r_red'] - r['r_ox'] for r in rows]
    flips = sum(1 for a, b in zip(sign, sign[1:]) if a > 0 >= b)
    assert flips == 1, flips
    near = min(rows, key=lambda r: abs(r['theta'] - d['theta_star']))
    assert near['r_red'] == pytest.approx(near['r_ox'], abs=5e-3)


def test_the_phase_boundaries_drawn_are_the_closed_forms(built):
    """Every point of the drawn curve is the engine's own K_max there."""
    doc, _ = built
    d = doc['phase']
    p = R.load_params()
    assert d['K_no_steady_state'] == pytest.approx(1.0, rel=1e-12)
    assert d['crossover_E'] == pytest.approx(
        d['theta_surface_max'] / d['theta_sol'], rel=1e-12)
    for q in d['second_phase']:
        w = R.usable_window(p, q['E'])
        assert q['K'] == pytest.approx(w['K_second_phase'], rel=1e-12), q
    # the curve stops exactly where it stops bounding anything
    last = d['second_phase'][-1]
    assert last['E'] == pytest.approx(d['crossover_E'], rel=1e-12)
    assert last['K'] == pytest.approx(d['K_no_steady_state'], rel=1e-12)


def test_the_two_marked_points_are_the_two_readings_of_the_same_particle(
        built):
    """Uniform and layer-resolved, and the window between them."""
    doc, _ = built
    pts = {q['label']: q for q in doc['phase']['points']}
    assert set(pts) == {'uniform', 'layer-resolved'}
    assert pts['uniform']['E'] == 1.0
    assert pts['uniform']['K'] == pytest.approx(0.0040161, rel=1e-4)
    assert pts['layer-resolved']['E'] == pytest.approx(1844.0)
    assert pts['layer-resolved']['K'] == pytest.approx(1.0, rel=1e-12)
    assert pts['layer-resolved']['K'] / pts['uniform']['K'] > 200.0


def test_the_sweep_written_out_is_the_engine_sweep(built):
    doc, _ = built
    p = R.load_params()
    for row in doc['sweep'][::20]:
        assert row['theta'] == pytest.approx(
            R.theta_of_ratio(row['K']), abs=1e-12), row
