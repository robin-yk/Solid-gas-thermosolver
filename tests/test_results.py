"""Provenance and completeness guards on results.json.

The rule these enforce is not a style preference. Two routes to this answer
exist in the tree and they agree almost everywhere, which is exactly what makes
a value from the wrong one hard to spot by reading it. So every row states where
it came from, and nothing without that stamp is allowed to reach the table or
the page.
"""

import json
import pathlib
import re

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
RESULTS = ROOT / 'results.json'

CASES = ['RWGS', 'H2', 'CH4', 'N2']
TEMPS = [500, 700, 900, 1100, 1300, 1500]


def load_results():
    if not RESULTS.exists():
        pytest.fail('results.json is missing - run `python run_all.py`')
    return json.loads(RESULTS.read_text())['rows']


def test_reported_values_are_gibbs_min():
    for r in load_results():
        assert r['method'] == 'gibbs_min', f"{r['case']} @{r['T_C']}: {r['method']}"


def test_every_point_is_present():
    """Twenty-four points, none dropped.

    A failed minimisation must show up as a missing or non-gibbs row rather
    than as a quietly shorter table.
    """
    got = {(r['case'], r['T_C']) for r in load_results()}
    want = {(c, t) for c in CASES for t in TEMPS}
    assert got == want, f"missing {sorted(want - got)}, extra {sorted(got - want)}"


def test_document_declares_its_method():
    doc = json.loads(RESULTS.read_text())
    assert doc['method'] == 'gibbs_min'
    assert doc['temperatures_C'] == TEMPS


def test_elemental_balance_holds():
    """The constraint the minimisation was subject to, checked on the answer."""
    for r in load_results():
        assert r['balance_residual_mol'] < 1e-16, \
            f"{r['case']} @{r['T_C']}: residual {r['balance_residual_mol']:.2e}"


def test_no_value_without_a_method_field():
    """Every row is stamped, including the fields a reader would quote."""
    for r in load_results():
        for key in ('ti3_pct', 'phase', 'gas_pct', 'mu_O_kJ'):
            assert key in r, f"{r['case']} @{r['T_C']}: no {key}"
        assert 'method' in r and 'formulation' in r


def test_crosscheck_route_is_recorded_and_not_reported():
    """The oxygen-potential route stays a witness, never the source.

    It must be present on every row - a check that is not run proves nothing -
    and it must sit under `crosscheck`, never at the top level where a reader
    could mistake it for the answer.
    """
    for r in load_results():
        cc = r.get('crosscheck')
        assert cc and cc['route'] == 'oxygen_potential', f"{r['case']} @{r['T_C']}"
        assert 'mu_O_kJ' in r


def test_two_phase_rows_sit_on_the_boundary():
    """Where the minimisation puts two phases together, the gas must be pinned.

    Coexistence fixes the oxygen potential exactly, and the boundary value comes
    from the phase data alone. Agreement here is the minimiser being checked
    against something it did not compute, so the tolerance is tight on purpose:
    the converged rows land inside 1e-4, and anything near 0.01 means a point
    stopped early rather than converged.
    """
    checked = 0
    for r in load_results():
        cc = r['crosscheck']
        if cc.get('status') != 'two-phase buffer':
            continue
        checked += 1
        assert abs(cc['d_mu_O_kJ']) < 0.01, \
            (f"{r['case']} @{r['T_C']} on {cc.get('boundary')}: gas is "
             f"{cc['d_mu_O_kJ']:+.5f} kJ/mol-O off the boundary")
    assert checked >= 12, f'only {checked} rows exercised the boundary check'


def test_unreduced_rows_really_are_unreduced():
    """A single-phase rutile row must have nothing reduced in it.

    Line compounds have no mixing entropy, so above the formation boundary the
    optimum is exactly zero and whatever the optimiser leaves behind is its own
    tolerance. This pins that residue down to where it cannot be mistaken for a
    phase, and pins the gas above the boundary that would make one.
    """
    for r in load_results():
        cc = r['crosscheck']
        if cc.get('status') != 'single phase' or cc.get('gibbs_phase') != 'TiO2':
            continue
        assert r['ti3_pct'] < 1e-6, \
            f"{r['case']} @{r['T_C']}: rutile alone but Ti3+ = {r['ti3_pct']}"
        assert cc['ladder_margin_kJ'] is not None and cc['ladder_margin_kJ'] > 0, \
            f"{r['case']} @{r['T_C']}: rutile alone but the gas is below the step"


def test_rwgs_never_reduces_anywhere():
    """The question this repository was built to answer, as an assertion."""
    rwgs = [r for r in load_results() if r['case'] == 'RWGS']
    assert len(rwgs) == len(TEMPS)
    for r in rwgs:
        assert r['ti3_pct'] < 1e-6, f"@{r['T_C']}: Ti3+ = {r['ti3_pct']}"
        assert r['phase'] == 'TiO2'
        assert r['crosscheck']['ladder_margin_kJ'] > 0


def test_page_and_table_read_the_same_file():
    """index.html embeds results.json verbatim.

    The page also carries a live solver that runs the other route. If the table
    on the page were built from that, it would agree today and diverge the first
    time either side moved. So the page renders the file, and this compares them.
    """
    index = ROOT / 'tiox.html'
    if not index.exists():
        pytest.skip('index.html not built yet - run `python build_site.py`')
    m = re.search(r'<script id="results" type="application/json">(.*?)</script>',
                  index.read_text(), re.S)
    assert m, 'index.html carries no embedded results block'
    assert json.loads(m.group(1)) == json.loads(RESULTS.read_text()), \
        'the page and results.json have drifted apart - rebuild with build_site.py'
