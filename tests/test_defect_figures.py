"""The defect workspace figures, and the export contract they carry.

Every quantitative plot on the workspace page is a publication figure, not
a screen chart: exactly 3.5 in square, a closed axis box with inward
ticks, no grid, no title inside the panel, and a 600 dpi raster of exactly
2100 px on a side. What the page previews is the file that downloads, so
the gate measures the same SVG string the browser mounts.

Nothing here may skip: a figure that will not build is a failure.
"""

import json
import pathlib
import re
import shutil
import subprocess
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
HARNESS = ROOT / 'tests' / 'js' / 'defect_figures_harness.js'

sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

import figure_registry as REG                 # noqa: E402

SQUARE = 'width="3.5in" height="3.5in" viewBox="0 0 252 252"'
WIDE = 'width="7in" height="3.5in" viewBox="0 0 504 252"'
ATTR = re.compile(r'([\w-]+)="([^"]*)"')

# four square plots and one landscape schematic
PLOTS = ('distribution', 'accessible', 'recovery', 'radial')
SCHEMATIC = 'particle'


@pytest.fixture(scope='module')
def figs():
    assert shutil.which('node') is not None, 'node is required'
    r = subprocess.run(['node', str(HARNESS), str(ROOT)],
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, f'harness failed:\n{r.stderr[-2000:]}'
    out = json.loads(r.stdout)
    assert not out['errors'], out['errors']
    assert len(out['ids']) == 5, out['ids']
    return out


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


def plots(figs):
    return [(fid, figs['svg'][fid]) for fid in PLOTS]


def test_every_plot_is_exactly_three_and_a_half_inches_square(figs):
    for fid, svg in plots(figs):
        assert SQUARE in svg, fid
        assert figs['px'][fid] == [2100, 2100], fid
    assert figs['px_wide'] == [4200, 2100]


def test_the_cross_section_is_seven_by_three_and_a_half_inches(figs):
    """The one landscape figure: the disc and its magnification have to
    sit side by side for the magnification to read as one."""
    svg = figs['svg'][SCHEMATIC]
    assert WIDE in svg
    assert figs['px'][SCHEMATIC] == [4200, 2100]


def test_the_cross_section_keeps_the_house_axis_style(figs):
    """Three panels, two of them framed, every tick turned inward."""
    svg = figs['svg'][SCHEMATIC]
    boxes = _boxes(svg)
    assert len(boxes) == 2, len(boxes)      # the window and the depth plot
    letters = re.findall(r'<text[^>]*font-weight="bold"[^>]*>([^<]*)</text>',
                         svg)
    assert letters == ['a', 'b', 'c'], letters
    tol = 0.01
    for a in (_attrs(m.group(1))
              for m in re.finditer(r'<line ([^>]*)/>', svg)):
        if abs(float(a.get('stroke-width', 0)) - 0.9) > 1e-9:
            continue
        xs = [float(a['x1']), float(a['x2'])]
        ys = [float(a['y1']), float(a['y2'])]
        assert any(b[0] - tol <= min(xs) and max(xs) <= b[2] + tol
                   and b[1] - tol <= min(ys) and max(ys) <= b[3] + tol
                   for b in boxes), (xs, ys)


def test_the_cross_section_draws_one_dot_per_oxygen_site(figs):
    """The magnified window must not read as one atomic row emptying.

    Each (1x1) column carries one bridging O and four O per d110 layer, so
    the subsurface class is drawn as 2*n_layers sub-rows of two sites and
    the horizontal and vertical scales are equal."""
    js = (ROOT / 'web' / 'figures_defect.js').read_text()
    assert 'var pxNm = sitePx / cell;' in js                 # isotropic
    assert 'var hPx = t1 * pxNm / 2;' in js                  # O sub-rows
    assert 'var rows = 2 * nLay;' in js
    assert 'one dot = one vacancy on one O site' in js
    assert "'shell ' + p.shell_nm.toFixed(2)" in js          # engine value
    assert 'O per cell' in js
    # the bulk class is a tint, never placed dots: the partition carries
    # no depth structure inside it
    assert 'p.theta.bulk' in js
    assert 'vacDot(f, siteX' in js
    assert "vacDot(f, siteX2(jj, s), subY(1 + rr), C.bulk" not in js


def test_every_plot_has_one_closed_axis_box(figs):
    for fid, svg in plots(figs):
        boxes = _boxes(svg)
        assert len(boxes) == 1, (fid, len(boxes))
        x0, y0, x1, y1 = boxes[0]
        assert abs((x1 - x0) - (y1 - y0)) < 1e-9, \
            (fid, 'the plotting panel must be square')


def test_ticks_point_into_the_panel(figs):
    """A tick drawn outward would fall outside the only box there is."""
    tol = 0.01
    for fid, svg in plots(figs):
        b = _boxes(svg)[0]
        for a in (_attrs(m.group(1))
                  for m in re.finditer(r'<line ([^>]*)/>', svg)):
            if abs(float(a.get('stroke-width', 0)) - 0.9) > 1e-9:
                continue
            xs = [float(a['x1']), float(a['x2'])]
            ys = [float(a['y1']), float(a['y2'])]
            assert (b[0] - tol <= min(xs) and max(xs) <= b[2] + tol
                    and b[1] - tol <= min(ys) and max(ys) <= b[3] + tol), \
                (fid, xs, ys)


def test_tick_labels_appear_on_the_bottom_and_left_only(figs):
    """No numeral sits above the top spine or right of the right one."""
    for fid, svg in plots(figs):
        b = _boxes(svg)[0]
        for m in re.finditer(r'<text ([^>]*)>([^<]*)</text>', svg):
            a, body = _attrs(m.group(1)), m.group(2)
            if 'rotate' in (a.get('transform') or ''):
                continue
            if not re.fullmatch(r'[-−0-9.,⁰-₟ ]+', body or 'x'):
                continue
            x, y = float(a['x']), float(a['y'])
            assert not (y < b[1] - 1), (fid, 'label above the box', body)
            assert not (x > b[2] + 1), (fid, 'label right of the box', body)


def test_no_grid_and_no_title_inside_the_panel(figs):
    for fid, svg in plots(figs):
        b = _boxes(svg)[0]
        assert 'font-weight="bold"' not in svg, (fid, 'in-panel title')
        # a grid rule would span the panel at the spine weight and be
        # neither a tick nor a boundary; boundaries are drawn at the guide
        # weight, so any full-width spine-weight rule is a grid line
        for a in (_attrs(m.group(1))
                  for m in re.finditer(r'<line ([^>]*)/>', svg)):
            if abs(float(a.get('stroke-width', 0)) - 0.9) > 1e-9:
                continue
            span = abs(float(a['x2']) - float(a['x1']))
            assert span < (b[2] - b[0]) - 1, (fid, 'grid line')


def test_type_is_the_declared_figure_scale(figs):
    """Tick labels at 8 pt, axis labels at 9 pt, and nothing else."""
    for fid, svg in figs['svg'].items():
        sizes = {float(v) for v in re.findall(r'font-size="([\d.]+)"', svg)}
        assert sizes, fid
        assert sizes <= {8.0, 9.0}, (fid, sorted(sizes))
        assert min(sizes) >= REG.SPEC['type_floor_pt'], fid


def test_figures_use_the_role_palette_only(figs):
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

    for fid, svg in figs['svg'].items():
        for c in set(re.findall(r'#[0-9A-Fa-f]{6}', svg)):
            if c.upper() in neutral or c in hues:
                continue
            assert any(tint_of(c, h) for h in hues), (fid, c)


def test_the_raster_carries_its_resolution(figs):
    """A canvas PNG has no pHYs, so the kit splices one in after IHDR.

    Without it the file opens at 2100 px of nothing rather than at 3.5 in.
    """
    p = figs['phys']
    assert p['tag'] == 'pHYs'
    assert p['grew_by'] == 21             # 4 length + 4 tag + 9 body + 4 CRC
    assert p['px_per_m'] == 23622         # 600 dpi
    assert p['unit'] == 1                 # metres
    assert p['ihdr_intact'] == 'IHDR'


def test_the_page_mounts_the_figures_it_builds():
    """The template holds a container for each figure and one menu each.

    Counted per workspace, not over the file: a figure added to one
    workspace must not be able to satisfy the count owed by another."""
    html = (ROOT / 'web' / 'ti_solver_template.html').read_text()
    for fid in ('figDistribution', 'figAccessible', 'figRecovery',
                'figParticle',
                'figRadial'):
        assert f'id="{fid}"' in html, fid
    defect = html.split('<div id="ws-defect"')[1].split('<div id="ws-redox"')[0]
    redox = html.split('<div id="ws-redox"')[1]
    # the defect workspace: four square figures and one 7-inch cross-section
    assert defect.count('class="pubfig"') == 4
    assert defect.count('class="pubfig wide"') == 1
    assert defect.count('data-fmt="svg"') == 5
    assert defect.count('data-fmt="png"') == 5
    # the redox workspace: the crossing and the phase map
    for fid in ('figCrossing', 'figPhase'):
        assert f'id="{fid}"' in redox, fid
    assert redox.count('class="pubfig"') == 2
    assert redox.count('data-fmt="svg"') == 2
    assert redox.count('data-fmt="png"') == 2
    built = (ROOT / 'docs' / 'ti_solver.html').read_text()
    for src in ('figkit.js', 'figures_defect.js'):
        body = (ROOT / 'web' / src).read_text()
        assert body in built, f'page is stale against {src}'
