"""The browser page must not drift from the Python package.

index.html recomputes everything client-side from thermo_data.json. If that
file goes stale, the page keeps working and quietly answers with old
coefficients, which is the failure mode worth a test.
"""

import json
import os
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, 'data')
WEB = os.path.join(ROOT, 'web')
DOCS = os.path.join(ROOT, 'docs')
SCRIPTS = os.path.join(ROOT, 'scripts')


@pytest.fixture(scope='module')
def exported():
    with open(os.path.join(DATA, 'thermo_data.json')) as fh:
        return json.load(fh)


def test_exported_file_is_current(exported):
    """Regenerating must reproduce the committed file byte for byte."""
    out = subprocess.run([sys.executable, os.path.join(SCRIPTS, 'export_thermo.py')], cwd=ROOT,
                         capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    with open(os.path.join(DATA, 'thermo_data.json')) as fh:
        assert json.load(fh) == exported, 'run export_thermo.py and commit the result'


def test_every_coefficient_matches_the_package(exported):
    from solidgas import shomate as S
    from solidgas import waldner as W

    for name, e in exported['shomate'].items():
        assert e['lo'] == list(S.SHOMATE[name][0])
        assert e['hi'] == list(S.SHOMATE[name][1])
        assert e['dHf'] == S.SHOMATE[name][2]
        assert e['brk'] == S.BREAKPOINT.get(name, S.DEFAULT_BREAKPOINT)

    for name, e in exported['waldner'].items():
        assert e['dHf'] == pytest.approx(W.dHf(name))
        assert e['S298'] == W.WALDNER[name]['S298']
        assert len(e['ranges']) == len(W.WALDNER[name]['ranges'])
        for got, want in zip(e['ranges'], W.WALDNER[name]['ranges']):
            assert (got['lo'], got['hi'], got['form'], got['dHt']) == \
                   (want[0], want[1], want[2], want[4])
            assert got['c'] == list(want[3])


def test_index_html_embeds_the_current_data_and_solver(exported):
    """index.html is built by inlining both, so it must contain them verbatim."""
    path = os.path.join(DOCS, 'tiox.html')
    if not os.path.exists(path):
        pytest.skip('index.html not built')
    html = open(path).read()
    with open(os.path.join(DATA, 'thermo_data.json')) as fh:
        assert fh.read() in html, 'index.html is stale - run build_site.py'
    for js in ('solver.js', 'page.js'):
        with open(os.path.join(WEB, js)) as fh:
            assert fh.read() in html, f'index.html is stale against {js} - run build_site.py'
    for token in ('__THERMO__', '__SOLVER__', '__PAGE__'):
        assert token not in html


def test_page_has_no_external_requests():
    """GitHub Pages serves it flat; nothing may be fetched from elsewhere."""
    path = os.path.join(DOCS, 'tiox.html')
    if not os.path.exists(path):
        pytest.skip('index.html not built')
    html = open(path).read()
    for bad in ('src="http', "src='http", 'href="http://', '@import',
                'fetch(', 'XMLHttpRequest'):
        assert bad not in html, f'page reaches outside: {bad}'
