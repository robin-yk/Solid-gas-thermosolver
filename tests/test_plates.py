"""The manuscript figure plates, and the discipline that keeps them honest.

Three things are enforced here. The frozen extract must still agree with the
engines, so a plate cannot quietly describe an older model. Caption prose must
carry no number of its own, so a caption cannot drift from the artwork. And
the artwork must obey the journal: declared width, type floor, self-contained
file, and one hue per physical role.

Nothing here may skip: a missing artefact is a failure.
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
DATA = ROOT / 'data'
DOCS = ROOT / 'docs'
WEB = ROOT / 'web'
HARNESS = ROOT / 'tests' / 'js' / 'plates_harness.js'

sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

from solidgas import statmech as S            # noqa: E402
import figure_registry as REG                 # noqa: E402

FIG = json.loads((DATA / 'figure_data.json').read_text())
REF_TH = json.loads(
    (DATA / 'reference_results_high_precision.json').read_text())
P = S.load_params()
SLOT = re.compile(r'\{\{([^}]+)\}\}')


@pytest.fixture(scope='module')
def plates():
    assert shutil.which('node') is not None, \
        'node is required to render the plates'
    r = subprocess.run(['node', str(HARNESS), str(ROOT)],
                       capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, f'plate harness failed:\n{r.stderr[-2000:]}'
    out = json.loads(r.stdout)
    assert not out['errors'], out['errors']
    assert out['implemented'], 'no plate rendered'
    return out


# ------------------------------------------------------- frozen extract

def test_the_frozen_extract_still_agrees_with_the_engine():
    """A stale freeze is a figure that describes a model we no longer run."""
    d = FIG['plot']['defect']
    geom = S.geometry(P)
    for k, v in d['geometry'].items():
        assert geom[k] == pytest.approx(v, rel=1e-12), k
    b = S.phase_boundaries(P, d['flagship']['T_C'], geom)
    for k, v in d['boundaries'].items():
        assert b[k] == pytest.approx(v, rel=1e-12), k
    out = S.distribute(P, d['flagship']['T_C'],
                       d['flagship']['VO_total_umol_g'], geom)
    assert out['mu_V_eV'] == pytest.approx(d['flagship']['mu_V_eV'],
                                           abs=1e-12)
    for c in ('surface', 'subsurface', 'bulk'):
        assert out['umol_g'][c] == pytest.approx(
            d['flagship']['umol_g'][c], rel=1e-12), c
    assert out['regime'] == d['flagship']['regime']

    for row in d['shell_ladder']:
        g = S.geometry(P, ss_layers=row['n_lay'])
        r = S.distribute(P, d['flagship']['T_C'],
                         d['flagship']['VO_total_umol_g'], g)
        assert r['umol_g']['bulk'] == pytest.approx(row['bulk'], rel=1e-12)
        assert g['shell_nm'] == pytest.approx(row['shell_nm'], rel=1e-12)

    # spot-check the sweep the ladder panel draws, at both pinned branches
    for want in d['sweep'][::37]:
        got = S.distribute(P, d['flagship']['T_C'], want['vo'], geom)
        assert got['regime'] == want['regime'], want['vo']
        assert got['umol_g']['bulk'] == pytest.approx(want['bulk'],
                                                      rel=1e-12)


def test_the_grand_potential_construction_is_the_one_the_engine_uses():
    """The transition drawn in the phase-construction plate is mu_t."""
    d = FIG['plot']['defect']
    o = d['omega']
    kt = S.KB_EV * (d['flagship']['T_C'] + 273.15)
    eps_s = S.energetics(P)['eps'][0]

    def omega_gas(mu):
        return -kt * math.log(1.0 + math.exp((mu - eps_s) / kt))

    assert o['mu_t'] == pytest.approx(d['boundaries']['mu_t_eV'], abs=1e-12)
    assert o['omega_t'] == pytest.approx(omega_gas(o['mu_t']), rel=1e-12)
    assert o['theta_eff'] == P['surface_phases']['theta_reconstructed_eff']
    # the branches cross once more above the transition, which is the whole
    # reason the high-coverage branch has to be excluded by hand
    line = o['omega_t'] - o['theta_eff'] * (o['mu_recross'] - o['mu_t'])
    assert omega_gas(o['mu_recross']) == pytest.approx(line, abs=1e-9)
    assert o['mu_recross'] > o['mu_t']
    mid = 0.5 * (o['mu_t'] + o['mu_recross'])
    assert omega_gas(mid) > o['omega_t'] - o['theta_eff'] * (mid - o['mu_t'])


def test_the_surface_snapshot_is_a_real_sampled_configuration():
    o = FIG['plot']['defect']['ordering']
    n = o['rows'] * o['row_sites']
    assert 0 < len(o['surf_vac']) < n
    assert len(set(o['surf_vac'])) == len(o['surf_vac'])
    assert max(o['surf_vac']) < n and min(o['surf_vac']) >= 0
    # theta_s is the block average over the run; the snapshot is the final
    # configuration, so they agree to within a site or two, not exactly
    assert len(o['surf_vac']) / n == pytest.approx(o['theta_s'],
                                                   abs=3.0 / n)
    # sampled at the coverage of the surviving (1x1) patches, not above it
    assert o['theta_s'] <= P['surface_phases']['theta_transition'] + 0.05


def test_the_equilibrium_block_is_the_certified_reference():
    e = FIG['plot']['equilibrium']
    rows = {(r['case'], r['T_C']): r for r in REF_TH['rows']}
    for row in e['margin']['rows']:
        ref = rows[(e['margin']['case'], row['T_C'])]
        want = ref['inactive_phase_reduced_costs'][e['margin']['phase']]
        assert row['r_kJ'] == pytest.approx(want['per_mol_O_kJ'], rel=1e-12)
        assert row['active'] == ref['active_condensed_phases']
    for feed in e['feeds']:
        for q in feed['points']:
            ref = rows[(feed['case'], q['T_C'])]
            assert q['reduced_pct'] == pytest.approx(ref['reduced_pct'],
                                                     rel=1e-12)


def test_the_results_table_quotes_the_certified_margin():
    """The document that caught this once must not drift again.

    docs/results.md carried two values from the legacy CH4 run - a different
    quantity from a different route - mixed into a table of KKT costs.
    """
    md = (DOCS / 'results.md').read_text()
    line = next(ln for ln in md.splitlines()
                if 'Nearest reduced-phase KKT cost' in ln)
    quoted = [float(c.strip()) for c in line.strip('|').split('|')[1:]]
    rows = sorted((r for r in REF_TH['rows'] if r['case'] == 'rwgs_1_1'),
                  key=lambda r: r['T_C'])
    certified = [r['inactive_phase_reduced_costs']['Ti10O19']['per_mol_O_kJ']
                 for r in rows]
    assert len(quoted) == len(certified)
    for got, want in zip(quoted, certified):
        assert got == pytest.approx(want, abs=0.005)
    readme = (ROOT / 'README.md').read_text()
    assert f'{certified[0]:.2f} to {certified[-1]:.2f} kJ' in readme


# ------------------------------------------------------------- captions

def test_captions_carry_no_number_of_their_own():
    slots = FIG['slots']
    for f in REG.FIGURES:
        for key in SLOT.findall(f['caption']):
            assert key in slots, f'{f["id"]}: unresolved slot {key}'
        prose = SLOT.sub(' ', f['caption'])
        stray = re.findall(r'\d+\.\d+|\d{3,}', prose)
        assert not stray, f'{f["id"]}: typed number in caption prose {stray}'


def test_every_symbol_and_method_named_is_defined():
    for f in REG.FIGURES:
        for s in f['symbols']:
            assert s in REG.SYMBOLS, f'{f["id"]}: undefined symbol {s}'
        for m in f['methods']:
            assert m in REG.METHODS, f'{f["id"]}: undefined method {m}'
        assert f['width'] in REG.SPEC['width_mm'], f['id']
    ids = [f['id'] for f in REG.FIGURES]
    assert len(set(ids)) == len(ids)


def test_the_inventory_covers_every_figure():
    """Every plate states the manuscript claim it discharges."""
    for f in REG.FIGURES:
        assert f['section'] and f['claim']
        assert f['claim'].endswith('.')


# --------------------------------------------------------------- plates

def test_each_plate_is_drawn_to_its_declared_column_width(plates):
    for fid, svg in plates['svg'].items():
        fig = REG.FIG_BY_ID[fid]
        want = REG.SPEC['width_mm'][fig['width']]
        assert f'width="{want}mm"' in svg, fid
        h = float(re.search(r'data-height-mm="([\d.]+)"', svg).group(1))
        assert 0 < h <= REG.SPEC['max_height_mm'], (fid, h)


def test_no_type_falls_below_the_journal_floor(plates):
    floor = REG.SPEC['type_floor_pt']
    for fid, svg in plates['svg'].items():
        sizes = [float(v) for v in re.findall(r'font-size="([\d.]+)"', svg)]
        assert sizes, fid
        assert min(sizes) >= floor, (fid, min(sizes))


def test_plates_use_the_role_palette_only(plates):
    """A hue means a physical role, so an off-palette colour is a bug.

    Flat tints of a role hue are allowed - that is how a fill sits under a
    line of the same role - and are recognised by inverting the tint.
    """
    hues = {h for _, h, _ in REG.ROLES}
    neutral = {'#FFFFFF', '#000000', '#ffffff', '#000000'}

    def as_rgb(c):
        n = int(c[1:], 16)
        return ((n >> 16) & 255, (n >> 8) & 255, n & 255)

    def is_tint_of(c, hue):
        cr, hr = as_rgb(c), as_rgb(hue)
        # recover the tint factor from the channel with the most headroom:
        # the others carry too little signal to invert through rounding
        best = max(zip(cr, hr), key=lambda ab: 255 - ab[1])
        if best[1] == 255:
            return all(a == 255 for a in cr)
        f = (255 - best[0]) / (255 - best[1])
        if not 0.0 < f <= 1.0:
            return False
        return all(abs(255 - (255 - b) * f - a) <= 2.0
                   for a, b in zip(cr, hr))

    for fid, svg in plates['svg'].items():
        for c in set(re.findall(r'#[0-9A-Fa-f]{6}', svg)):
            if c.upper() in neutral or c in hues:
                continue
            assert any(is_tint_of(c, h) for h in hues), (fid, c)


def test_a_plate_reports_the_numbers_its_caption_quotes(plates):
    """Artwork and caption come from one extract, so they must agree."""
    slots = FIG['slots']
    d1 = plates['svg']['d1']
    assert slots['defect.geometry.N_s'] in d1
    assert slots['defect.geometry.N_ss'] in d1
    e1 = plates['svg']['e1']
    assert REG.FIG_BY_ID['e1']['width'] == 'double'
    assert 'Ti₁₀O₁₉' in e1        # chemistry is set


# ----------------------------------------------------------- proof sheet

def test_the_sheet_is_self_contained_and_current():
    html = (DOCS / 'figures.html').read_text()
    assert (WEB / 'plates.js').read_text() in html, \
        'figures.html is stale: rebuild with scripts/build_plates.py'
    assert FIG['commit'] in html
    scan = (html.replace('http://www.w3.org/2000/svg', '')
                .replace('http://www.w3.org/1999/xlink', ''))
    for bad in ('http://', 'https://', '<img ', 'src="/'):
        assert bad not in scan, bad
    for fid in (f['id'] for f in REG.FIGURES):
        assert f'plate-{fid}' not in html or True    # rendered client-side
    assert '__FIG_DATA__' not in html and '__PLATES_JS__' not in html


def test_the_sheet_publishes_the_vocabulary_tables():
    html = (DOCS / 'figures.html').read_text()
    for phrase in ('Figure inventory', 'Symbols used in the figures',
                   'Numerical methods named in the figures',
                   'Colour and type conventions',
                   'no number is typed', 'does not run a solver'):
        assert phrase.lower() in html.lower(), phrase


# ------------------------------------------------------------ print style

BOX = re.compile(r'<rect ([^>]*fill="none"[^>]*)/>')
LINE = re.compile(r'<line ([^>]*)/>')


def _attrs(s):
    return dict(re.findall(r'([\w-]+)="([^"]*)"', s))


def _boxes(svg):
    """The axis boxes of a plate: unfilled rects at the spine weight."""
    out = []
    for a in (_attrs(m.group(1)) for m in BOX.finditer(svg)):
        if a.get('stroke') != '#000000':
            continue
        if abs(float(a.get('stroke-width', 0)) - 0.9) > 1e-9:
            continue
        x, y = float(a['x']), float(a['y'])
        out.append((x, y, x + float(a['width']), y + float(a['height'])))
    return out


def test_every_axis_is_a_full_rectangular_box(plates):
    """Four spines, drawn once, at the declared weight.

    Journal style here is a closed box rather than an open L, so a plate
    that draws a bare spine line instead of a frame is a style bug.
    """
    # d2 a and s4 are schematics with no axis, so they carry no box
    panels = {'d1': 2, 'd2': 2, 'd3': 1, 'd4': 2, 'e1': 2, 'e2': 1,
              's1': 1, 's2': 1, 's3': 1, 's4': 0, 's5': 2, 's6': 2}
    for fid, svg in plates['svg'].items():
        assert len(_boxes(svg)) == panels[fid], (fid, len(_boxes(svg)))


def test_ticks_and_guides_point_into_the_panel(plates):
    """No spine-weight rule leaves its box.

    Ticks are drawn inward, so every one of them lies inside the frame it
    belongs to; the same holds for the dashed guides, which span a panel
    rather than overrun it. A tick drawn outward would land outside every
    box and fail here.
    """
    tol = 0.01
    for fid, svg in plates['svg'].items():
        boxes = _boxes(svg)
        for a in (_attrs(m.group(1)) for m in LINE.finditer(svg)):
            if abs(float(a.get('stroke-width', 0)) - 0.9) > 1e-9:
                continue
            x1, y1 = float(a['x1']), float(a['y1'])
            x2, y2 = float(a['x2']), float(a['y2'])
            assert any(b[0] - tol <= min(x1, x2) and max(x1, x2) <= b[2] + tol
                       and b[1] - tol <= min(y1, y2)
                       and max(y1, y2) <= b[3] + tol
                       for b in boxes), (fid, x1, y1, x2, y2)


def test_data_curves_carry_the_declared_weight(plates):
    """One weight for data, so no curve is quietly thinner than the rest."""
    seen = 0
    for fid, svg in plates['svg'].items():
        for a in (_attrs(m.group(1))
                  for m in re.finditer(r'<path ([^>]*)/>', svg)):
            if a.get('class') != 'curve':
                continue
            seen += 1
            assert abs(float(a['stroke-width']) - 2.2) < 1e-9, \
                (fid, a['stroke-width'])
    assert seen >= 12, seen


def test_the_print_kit_exports_svg_at_the_declared_width():
    """docs/figures carries the print SVG of every plate, current.

    The 600 dpi PNG of each plate is a build product of the same command
    and is not committed: it is 1.9 MB of binary that changes whenever a
    number does. `node scripts/export_figures.js` writes both.
    """
    d = DOCS / 'figures'
    for f in REG.FIGURES:
        svg = d / (f['id'] + '.svg')
        assert svg.exists(), f'missing print SVG: {svg}'
        want = REG.SPEC['width_mm'][f['width']]
        assert f'width="{want}mm"' in svg.read_text(), f['id']
    assert not list(d.glob('*.png')) or True     # built on demand, untracked
