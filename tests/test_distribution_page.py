"""The vacancy-distribution workspace reads paper_outputs/tof_range.json
and draws it. The gates: that file is what the one-model run wrote, it is
current, the page carries exactly it, and the rendered panels were made
from its current values."""

import csv
import json
import pathlib
import subprocess
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
TOF = ROOT / 'paper_outputs' / 'tof_range.json'
PILOT = ROOT / 'pilot' / 'titania-super-multiscale' / 'outputs'
PAGE = ROOT / 'docs' / 'index.html'


@pytest.fixture(scope='module')
def doc():
    return json.loads(TOF.read_text())


def test_regenerating_the_json_is_byte_identical():
    before = TOF.read_bytes()
    r = subprocess.run([sys.executable, 'scripts/export_tof_range.py'],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert TOF.read_bytes() == before, 'run scripts/export_tof_range.py and commit'


def test_the_ranges_are_the_run_summary_as_written(doc):
    with (PILOT / 'sample_tof_range.csv').open() as fh:
        rows = {r['sample']: r for r in csv.DictReader(fh)}
    for s in doc['samples']:
        assert s['summary'] == rows[s['sample']]


def test_every_case_is_the_run_row_to_four_figures(doc):
    grid = doc['axes']
    eps = {f"{e['value']:g} ({e['key']})": i for i, e in enumerate(grid['eps'])}
    by = {s['sample']: s for s in doc['samples']}
    with (PILOT / 'model_cases.csv').open() as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 6 * 720 * 5
    for r in rows[::97]:
        i = grid['energy_map'].index(r['energy_map'])
        i = i * len(grid['cutoff_nm']) + grid['cutoff_nm'].index(float(r['cutoff_nm']))
        i = i * len(grid['dG_eV']) + grid['dG_eV'].index(float(r['dG_eV']))
        i = i * len(grid['f110']) + grid['f110'].index(float(r['f110']))
        i = i * len(grid['eps']) + eps[r['eps']]
        c, j = by[r['sample']]['cases'], doc['reactive'].index(r['reactive'])
        assert c['sites'][i][j] == pytest.approx(float(r['N_react_umol_g']), rel=5e-4)
        assert c['theta'][i] == pytest.approx(float(r['theta_bri']), rel=5e-4)
        assert c['below'][i][j] == (r['site_class'] != 'OK')
        tof = by[r['sample']]['rate_co_umol_g_s'] / c['sites'][i][j]
        assert tof == pytest.approx(float(r['TOF_s_1']), rel=1e-3)


def test_the_page_carries_the_json_unchanged(doc):
    assert json.dumps(doc, separators=(',', ':')) in PAGE.read_text(), \
        'the page is stale - run scripts/build_site.py'


RENDER = ROOT / 'web' / 'render'


def test_the_surface_renders_use_the_current_site_fractions(doc):
    """slab_sites.json records the fractions each slab image was rendered
    at; they must be the starting-parameter fractions in tof_range.json."""
    meta = json.loads((RENDER / 'slab_sites.json').read_text())
    st, ax = meta['start'], doc['axes']
    i = ax['energy_map'].index(st['energy_map'])
    i = i * len(ax['cutoff_nm']) + ax['cutoff_nm'].index(st['cutoff_nm'])
    i = i * len(ax['dG_eV']) + ax['dG_eV'].index(st['dG_eV'])
    i = i * len(ax['f110']) + ax['f110'].index(st['f110'])
    i = i * len(ax['eps']) + [e['key'] for e in ax['eps']].index(st['eps'])
    assert meta['point'] == i
    cap = doc['capacity_umol_g']['%g' % st['f110']]
    for s in doc['samples']:
        pools = dict(zip(doc['pools'], s['cases']['pools'][i]))
        want = dict(BRI=s['cases']['theta'][i], IPL=pools['basal_L1'] / cap['basal_L1'],
                    SBR=pools['L1_subbridging'] / cap['L1_subbridging'],
                    L24=pools['subsurface_L2_4'] / cap['subsurface_L2_4'])
        assert meta['fractions'][s['sample']] == pytest.approx(want, rel=1e-12), s['sample']
        assert (RENDER / ('slab_%s.webp' % s['sample'])).exists()


def test_the_slab_keeps_each_vacancy_as_the_inventory_rises(doc):
    """Fixed site ranks: a sample with larger fractions keeps every vacancy
    of a smaller one, so the vacant counts never fall along the series."""
    meta = json.loads((RENDER / 'slab_sites.json').read_text())
    order = sorted(doc['samples'], key=lambda s: meta['fractions'][s['sample']]['BRI'])
    counts = [meta['vacant'][s['sample']]['BRI'] for s in order]
    assert counts == sorted(counts)
    assert all(0 <= meta['vacant'][k]['BRI'] <= meta['sites']['BRI'] for k in meta['vacant'])


def test_the_page_carries_the_renders():
    html = PAGE.read_text()
    assert 'RENDER:' not in html and '/*RENDERDATA*/' not in html
    for name in ('particle_default.webp', 'slab_R600.webp', 'particle_mask.png'):
        assert (RENDER / name).exists(), name
    assert html.count('data:image/webp;base64,') >= 10
