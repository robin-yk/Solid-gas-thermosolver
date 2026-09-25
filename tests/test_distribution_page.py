"""The vacancy-distribution workspace reads paper_outputs/tof_range.json
and draws it. The gates: that file is what the one-model run wrote, it is
current, the page carries exactly it, and the 3D bundle is built from the
sources committed next to it."""

import csv
import hashlib
import json
import pathlib
import subprocess
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
TOF = ROOT / 'paper_outputs' / 'tof_range.json'
PILOT = ROOT / 'pilot' / 'titania-super-multiscale' / 'outputs'
PAGE = ROOT / 'docs' / 'index.html'
SRC = ROOT / 'web' / '3d'


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


def test_the_3d_bundle_is_built_from_the_committed_sources():
    """Same hash as web/3d/build.mjs: file name, a NUL, the bytes."""
    names = ['common.js', 'main.js', 'particle.js', 'slab.js', 'assets/studio.jpg', 'package-lock.json']
    assert f"'{names[0]}'" in (SRC / 'build.mjs').read_text()
    h = hashlib.sha256()
    for n in names:
        h.update(n.encode() + b'\0' + (SRC / n).read_bytes())
    head = (ROOT / 'web' / 'vacancy3d.js').read_text().split('\n', 1)[0]
    assert h.hexdigest() in head, 'run npm run build in web/3d and commit web/vacancy3d.js'


def test_the_3d_figures_draw_counts_and_fractions_only():
    """The 3D code gets pool sizes and site fractions from the page and
    must not carry a physical number of its own beyond the rutile cell."""
    src = (SRC / 'slab.js').read_text() + (SRC / 'particle.js').read_text()
    assert 'const A = 4.594, C = 2.959, U = 0.305' in src
    for word in ('TOF', 'rate_co', 'inventory'):
        assert word not in src
