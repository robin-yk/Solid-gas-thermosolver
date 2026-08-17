"""The browser engine must reproduce the Python answer.

`oxide.js` is a port, not a reimplementation, and the one piece that could not
be carried across is the optimiser: SciPy's SLSQP has no browser equivalent, so
the JavaScript runs a projected gradient instead. That substitution is the whole
risk in the port, and this is what holds it down.

Comparisons are made only on species that are actually present. A relative test
on a mole fraction of 1e-11 measures where the optimiser stopped, not whether
the port is faithful, and demanding agreement there would make this test fail
for a reason that has nothing to do with what it is checking.
"""

import json
import os
import pathlib
import shutil
import subprocess

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
HARNESS = ROOT / 'tests' / 'js' / 'port_harness.js'
REFERENCE = ROOT / 'oxide_reference.json'
TEMPS = [500, 700, 900, 1100, 1300, 1500]
PRESENT_PCT = 1e-4          # a gas species below this is not being compared
TOL = 1e-5                  # five significant figures


@pytest.fixture(scope='module')
def pair():
    if shutil.which('node') is None:
        pytest.skip('node is not available')
    for f in (REFERENCE, ROOT / 'oxide_data.json', ROOT / 'oxide.js'):
        if not f.exists():
            pytest.skip(f'{f.name} missing - run export_oxide.py and oxide_reference.py')
    r = subprocess.run(['node', str(HARNESS), str(ROOT)],
                       capture_output=True, text=True, timeout=900)
    if r.returncode != 0:
        pytest.fail(f'harness failed:\n{r.stderr[-2000:]}')
    js = {(x['system'], x['T_C']): x for x in json.loads(r.stdout)}
    py = {(x['system'], x['T_C']): x for x in json.loads(REFERENCE.read_text())['rows']}
    return py, js


def rel(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-12)


def test_every_reference_point_was_run(pair):
    py, js = pair
    assert set(js) == set(py)
    assert {t for _, t in py} == set(TEMPS)
    assert {s for s, _ in py} == {'Ti-O', 'Ce-O'}


def test_method_is_gibbs_min_on_both_sides(pair):
    py, js = pair
    for k in py:
        assert py[k]['method'] == 'gibbs_min', k
        assert js[k]['method'] == 'gibbs_min', k
        assert js[k]['formulation'] == py[k]['formulation'], k


def test_gas_composition_agrees(pair):
    py, js = pair
    for k, p in py.items():
        for sp, v in p['gas_pct'].items():
            if v <= PRESENT_PCT:
                continue
            d = rel(js[k]['gas_pct'][sp], v)
            assert d < TOL, f'{k} {sp}: py={v:.8g} js={js[k]["gas_pct"][sp]:.8g} rel={d:.2e}'


def test_conversion_and_quotient_agree(pair):
    py, js = pair
    for k, p in py.items():
        for field in ('conversion_CO2_pct', 'Q_rwgs'):
            d = rel(js[k][field], p[field])
            assert d < TOL, f'{k} {field}: py={p[field]:.8g} js={js[k][field]:.8g} rel={d:.2e}'


def test_reduction_extent_agrees(pair):
    """Offset by one so a shared zero compares as equal rather than as 0/0."""
    py, js = pair
    for k, p in py.items():
        d = rel(js[k]['reduced_pct'] + 1, p['reduced_pct'] + 1)
        assert d < TOL, f'{k} reduced: py={p["reduced_pct"]:.4g} js={js[k]["reduced_pct"]:.4g}'


def test_dominant_phase_agrees(pair):
    """Which phase actually holds the metal. This is the answer to the question."""
    py, js = pair
    for k, p in py.items():
        a = max(p['phase_split_pct'], key=p['phase_split_pct'].get)
        b = max(js[k]['phase_split_pct'], key=js[k]['phase_split_pct'].get)
        assert a == b, f'{k}: py={a} js={b}'


@pytest.mark.xfail(reason='the residue under the dominant phase is arbitrary in both '
                          'implementations; the active-set rewrite is what makes this pass',
                   strict=False)
def test_phase_assemblage_agrees(pair):
    """The full set of phases present, not just the dominant one.

    This cannot pass while absent phases sit at a floor instead of at zero, and
    the reason is worth keeping rather than deleting the test. Across the six
    ceria points the residue in Ce2O3 runs 0 to 1.8e-6 percent on the JavaScript
    side and 0 to 9.2e-29 on the Python side, uncorrelated - two optimisers
    stopping in different places on a variable whose true value is exactly zero.

    So no threshold separates residue from a phase genuinely present in trace
    amounts, because the residue is not small in a predictable way. Under sealed
    N2 a share of 1e-36 percent is the whole answer; here 1.8e-6 percent is
    nothing. Only holding absent phases at zero makes this question answerable.
    """
    py, js = pair
    for k, p in py.items():
        a = {q for q, v in p['phase_split_pct'].items() if v > 1e-6}
        b = {q for q, v in js[k]['phase_split_pct'].items() if v > 1e-6}
        assert a == b, f'{k}: py={sorted(a)} js={sorted(b)}'


def test_total_gibbs_energy_agrees(pair):
    py, js = pair
    for k, p in py.items():
        d = rel(js[k]['G_rel_kJ'], p['G_rel_kJ'])
        assert d < 1e-4, f'{k} G_rel: py={p["G_rel_kJ"]:.8g} js={js[k]["G_rel_kJ"]:.8g} rel={d:.2e}'


def test_elemental_balance_holds_in_the_port(pair):
    """Per element, not just the largest - the constraint is on each of them."""
    py, js = pair
    for k, x in js.items():
        for e, resid in x['balance_by_element'].items():
            assert abs(resid) < 1e-16, f'{k} element {e}: residual {resid:.2e}'
        assert set(x['elements']) == ({'C', 'H', 'O', 'Ti'} if k[0] == 'Ti-O'
                                      else {'C', 'H', 'O', 'Ce'})


@pytest.mark.xfail(reason='absent phases sit at the log-coordinate floor, not at zero; '
                          'the active-set rewrite is what makes this pass',
                   strict=False)
def test_absent_phases_are_exactly_zero(pair):
    """The defect, written down as a failing test rather than as a comment.

    A line compound has no mixing entropy, so above its formation boundary its
    optimum is exactly zero. Logarithmic coordinates cannot represent that, and
    the residue left behind has already been misread once as a phase.
    """
    py, js = pair
    for k, x in js.items():
        for phase, v in x['phase_split_pct'].items():
            assert v == 0.0 or v > 1e-6, f'{k} {phase}: {v:.3e} is neither zero nor a phase'
