"""The two live figures on the equilibrium workspace.

Three separate claims are gated here.

The first is about the quantity the figures are drawn on. e2 plots the
KKT reduced cost per mol of oxygen, which is the comparable number at the
single frozen condition the plate reports. Per mol O divides by the
oxygen deficiency relative to the host, and that denominator changes sign
as soon as a candidate is less deficient than the host - so on a page
whose host, feed and temperature all move, it goes negative and a bar
chart of it would look like a KKT violation that is not one. The live
figures use the per-formula cost instead, and the first test is the
reason: over the conditions this workspace admits, per formula is never
negative and per mol O often is.

The second is parity. web/figures_thermo.js carries two payload builders
that only reshape a solver result; the Python mirror below is written out
longhand and the two are held to equality, so a page that quietly grew
its own arithmetic between the engine and the drawing fails here.

The third is that the page actually mounts them, refreshes them, and did
not disturb the frozen manuscript plates while doing it.

Nothing may skip: a module that will not load is a failure.
"""

import hashlib
import json
import pathlib
import re
import shutil
import subprocess
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from solidgas import activeset as A                    # noqa: E402

HARNESS = ROOT / 'tests' / 'js' / 'thermo_figures_harness.js'
PAGE = ROOT / 'docs' / 'ti_solver.html'
TEMPLATE = ROOT / 'web' / 'ti_solver_template.html'


@pytest.fixture(scope='module')
def js():
    assert shutil.which('node') is not None, 'node is required'
    r = subprocess.run(['node', str(HARNESS), str(ROOT)],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-2000:]
    return json.loads(r.stdout)


@pytest.fixture(scope='module')
def html():
    return PAGE.read_text()


CASES = {
    'rwgs': dict(feed={'CO2': 1, 'H2': 1}, mass_g=0.100, T_C=900,
                 initial_solid='TiO2'),
    'magneli': dict(feed={'H2': 1}, mass_g=0.002, T_C=1200,
                    initial_solid='TiO2'),
    'reox': dict(feed={'CO2': 1}, mass_g=0.100, T_C=600,
                 initial_solid='Ti4O7'),
    'h2co2': dict(feed={'H2': 99, 'CO2': 1}, mass_g=0.002, T_C=1000,
                  initial_solid='TiO2'),
    'n2': dict(feed={'N2': 1}, mass_g=0.100, T_C=1400,
               initial_solid='Ti3O5'),
}


def solve_case(c, T_C):
    return A.solve(c['feed'], T_C + 273.15, P_atm=1.0, V_cm3=25.0,
                   mass_g=c['mass_g'], initial_solid=c['initial_solid'])


def finite(v):
    return v is not None and v == v and abs(v) != float('inf')


# --------------------------------------------------- the drawn quantity

def test_per_formula_is_the_only_reduced_cost_a_live_figure_can_draw():
    """Per formula is positive at every reported solution; per mol O is not.

    Both statements have to hold for the figures to be drawn the way they
    are. If per formula ever went negative the bar chart would be lying
    about the KKT test; if per mol O never went negative there would be
    no reason to depart from the plate."""
    feeds = [{'CO2': 1, 'H2': 1}, {'H2': 1}, {'CO2': 1}, {'N2': 1},
             {'H2': 99, 'CO2': 1}, {'CO2': 3, 'H2': 3, 'CO': 2, 'H2O': 2}]
    hosts = ['TiO2', 'Ti10O19', 'Ti5O9', 'Ti3O5']
    seen = neg_formula = neg_per_O = 0
    worst_formula = float('inf')
    for feed in feeds:
        for host in hosts:
            for T_C in range(300, 1651, 225):
                r = solve_case(dict(feed=feed, mass_g=0.002,
                                    initial_solid=host), T_C)
                for rc in r['inactive_phase_reduced_costs'].values():
                    pf, po = rc['per_formula_kJ'], rc['per_mol_O_kJ']
                    if not finite(pf):
                        continue
                    seen += 1
                    worst_formula = min(worst_formula, pf)
                    neg_formula += pf < -1e-9
                    neg_per_O += finite(po) and po < -1e-9
    assert seen > 1000, seen
    assert neg_formula == 0, (neg_formula, worst_formula)
    assert worst_formula > 0.0
    # and the normalisation the plate uses does flip, which is the point
    assert neg_per_O > 0.2 * seen, (neg_per_O, seen)


def test_the_basis_string_still_names_the_denominator_that_flips():
    r = solve_case(CASES['rwgs'], 900)
    assert 'def_j = rho_host*t_j - o_j' in r['reduced_cost_basis']


# ------------------------------------------------------------- parity

def python_optimality(R):
    """The Python mirror of ThermoFigures.optimalityData."""
    tol = R['solver_tolerance']['degeneracy_kJ_per_mol_O']
    costs, unbounded = [], []
    for ph, rc in R['inactive_phase_reduced_costs'].items():
        pf = rc['per_formula_kJ']
        if not finite(pf):
            unbounded.append(ph)
            continue
        po = rc['per_mol_O_kJ']
        costs.append({'phase': ph, 'r_kJ': pf, 'r_per_mol_O_kJ': po,
                      'degenerate': finite(po) and abs(po) < tol})
    costs.sort(key=lambda q: q['r_kJ'])
    cand = R['candidates']
    return {
        'T_C': R['T_C'],
        'winner': list(R['active_condensed_phases']),
        'costs': costs,
        'unbounded': unbounded,
        'candidates': {
            'total': len(cand),
            'feasible': sum(1 for c in cand if c['feasible']),
            'optimal': sum(1 for c in cand if c['feasible'] and c['kkt_ok']),
        },
    }


def python_margin(rows, current_T_C):
    """The Python mirror of ThermoFigures.marginData."""
    pts = []
    for r in rows:
        best, who = None, None
        for ph, rc in r['inactive_phase_reduced_costs'].items():
            pf = rc['per_formula_kJ']
            if not finite(pf):
                continue
            if best is None or pf < best:
                best, who = pf, ph
        pts.append({'T_C': r['T_C'], 'margin_kJ': best, 'nearest': who,
                    'winner': list(r['active_condensed_phases'])})
    segs = []
    last = None
    for q in pts:
        key = '+'.join(q['winner'])
        if key != last:
            if segs:
                segs[-1]['to_C'] = q['T_C']
            segs.append({'from_C': q['T_C'], 'to_C': q['T_C'],
                         'phases': list(q['winner'])})
            last = key
        else:
            segs[-1]['to_C'] = q['T_C']
    cur = min(pts, key=lambda q: abs(q['T_C'] - current_T_C)) if pts else None
    return {
        'T_lo': pts[0]['T_C'], 'T_hi': pts[-1]['T_C'],
        'rows': pts, 'segments': segs,
        'current': cur and {'T_C': current_T_C, 'margin_kJ': cur['margin_kJ'],
                            'winner': list(cur['winner']),
                            'nearest': cur['nearest']},
    }


def close(a, b, tol=1e-9):
    scale = max(abs(a), abs(b), 1e-300)
    return abs(a - b) / scale <= tol


def test_the_bar_payload_is_what_python_would_build(js):
    for name, c in CASES.items():
        want = python_optimality(solve_case(c, c['T_C']))
        got = js['cases'][name]['optimality']
        assert got['winner'] == want['winner'], name
        assert got['unbounded'] == want['unbounded'], name
        assert got['candidates'] == want['candidates'], name
        assert close(got['T_C'], want['T_C']), name
        assert len(got['costs']) == len(want['costs']), name
        for g, w in zip(got['costs'], want['costs']):
            assert g['phase'] == w['phase'], (name, g, w)
            assert g['degenerate'] == w['degenerate'], (name, g, w)
            assert close(g['r_kJ'], w['r_kJ']), (name, g, w)
            # null on both sides is the agreement, not a missing number:
            # it is how each engine says the deficiency is exactly zero
            gp, wp = g['r_per_mol_O_kJ'], w['r_per_mol_O_kJ']
            assert (gp is None) == (wp is None), (name, g, w)
            assert gp is None or close(gp, wp), (name, g, w)


def test_the_bar_is_sorted_so_the_phase_that_appears_first_is_first(js):
    for name, case in js['cases'].items():
        costs = case['optimality']['costs']
        assert costs, name
        assert costs == sorted(costs, key=lambda q: q['r_kJ']), name
        assert all(q['r_kJ'] > 0.0 for q in costs), name


def test_the_sweep_payload_is_what_python_would_build(js):
    n, lo, hi = js['sweep_n'], js['T_lo'], js['T_hi']
    for name, c in CASES.items():
        rows = [solve_case(c, lo + (hi - lo) * i / (n - 1)) for i in range(n)]
        want = python_margin(rows, c['T_C'])
        got = js['cases'][name]['margin']
        assert close(got['T_lo'], want['T_lo']) and close(got['T_hi'],
                                                          want['T_hi']), name
        assert len(got['rows']) == n, name
        for g, w in zip(got['rows'], want['rows']):
            assert g['winner'] == w['winner'], (name, g, w)
            assert g['nearest'] == w['nearest'], (name, g, w)
            assert close(g['margin_kJ'], w['margin_kJ']), (name, g, w)
        assert len(got['segments']) == len(want['segments']), name
        for g, w in zip(got['segments'], want['segments']):
            assert g['phases'] == w['phases'], (name, g, w)
            assert close(g['from_C'], w['from_C']), (name, g, w)
            assert close(g['to_C'], w['to_C']), (name, g, w)
        assert got['current'] == pytest.approx(want['current'], rel=1e-9) \
            or (got['current']['winner'] == want['current']['winner']
                and got['current']['nearest'] == want['current']['nearest']
                and close(got['current']['margin_kJ'],
                          want['current']['margin_kJ'])), name


def test_a_segment_edge_is_where_the_assemblage_actually_changes(js):
    """The ladder strip and the curve's zeros are the same crossings.

    Every segment boundary must separate two different winners, the
    segments must tile the swept range without a gap, and the margin must
    come down near zero at each edge - that is what makes the sawtooth an
    explanation of the strip rather than noise beside it."""
    for name, case in js['cases'].items():
        m = case['margin']
        segs, rows = m['segments'], m['rows']
        assert segs[0]['from_C'] == pytest.approx(m['T_lo'])
        assert segs[-1]['to_C'] == pytest.approx(m['T_hi'])
        for a, b in zip(segs, segs[1:]):
            assert a['to_C'] == pytest.approx(b['from_C']), name
            assert a['phases'] != b['phases'], (name, a, b)
        if len(segs) == 1:
            continue
        step = (m['T_hi'] - m['T_lo']) / (len(rows) - 1)
        floor = max(q['margin_kJ'] for q in rows)
        for edge in [s['from_C'] for s in segs[1:]]:
            near = [q['margin_kJ'] for q in rows
                    if abs(q['T_C'] - edge) <= 1.01 * step]
            assert min(near) < 0.25 * floor, (name, edge, min(near), floor)


def test_two_fields_only_ever_meet_across_a_tie_line_the_grid_missed(js):
    """A one-phase field can only be left through a two-phase tie line.

    The strip is shaded on that distinction and nothing else, so a place
    where two single-phase runs meet directly is either a violation of
    the phase rule or a tie line narrower than one grid step. It is
    always the second, and this refines every such place to prove it -
    which is also the honest statement of the figure's resolution: at the
    page's 121 points a tie line under about 11 °C wide goes unshaded,
    though the curve still comes down to zero there."""
    checked = 0
    for name, case in js['cases'].items():
        c = CASES[name]
        m = case['margin']
        step = (m['T_hi'] - m['T_lo']) / (len(m['rows']) - 1)
        for a, b in zip(m['segments'], m['segments'][1:]):
            if len(a['phases']) > 1 or len(b['phases']) > 1:
                continue
            assert a['phases'] != b['phases'], (name, a, b)
            found = None
            for i in range(41):
                T = a['to_C'] - step + 2 * step * i / 40
                w = solve_case(c, T)['active_condensed_phases']
                if len(w) > 1:
                    found = (T, w)
                    break
            assert found, (name, a, b, 'no tie line between two fields')
            assert set(found[1]) == set(a['phases']) | set(b['phases']), \
                (name, a, b, found)
            checked += 1
    assert checked, 'no field-to-field meeting in the grid to refine'


# ---------------------------------------------------------------- page

def test_the_page_carries_both_modules_and_neither_token_survives(html):
    for token in ('__TH_FIGS__', '__TH_UI__'):
        assert token in TEMPLATE.read_text(), token
        assert token not in html, token
    for src in ('figures_thermo.js', 'thermo_ui.js'):
        body = (ROOT / 'web' / src).read_text()
        needle = body.split('\n')[-4].strip()
        assert needle and needle in html, src
    assert 'ThermoFigures' in html
    assert 'optimalityData' in html and 'marginData' in html


def test_the_equilibrium_workspace_mounts_two_figures(html):
    thermo = html.split('<div id="ws-thermo"')[1].split('<div id="ws-defect"')[0]
    assert thermo.count('class="pubfig"') == 2, 'two figures, no more'
    assert thermo.count('data-fmt="svg"') == 2
    assert thermo.count('data-fmt="png"') == 2
    for host in ('figOptimality', 'figMargin'):
        assert 'id="%s"' % host in thermo, host


def test_the_sweep_is_behind_a_button_over_the_range_the_panel_admits(html):
    ui = (ROOT / 'web' / 'thermo_ui.js').read_text()
    assert 'id="sweepBtn"' in html and 'id="sweepNote"' in html
    assert "btn.addEventListener('click', sweep)" in ui
    # the range comes off the input, not out of a constant in the figure
    assert "el.getAttribute('min')" in ui and "el.getAttribute('max')" in ui
    thermo = html.split('<div id="ws-thermo"')[1].split('<div id="ws-defect"')[0]
    tin = re.search(r'<input[^>]*id="T"[^>]*>', thermo).group(0)
    assert 'min="300"' in tin and 'max="1650"' in tin, tin


def test_a_solve_redraws_the_bar_and_the_curve_that_was_already_drawn(html):
    page = (ROOT / 'web' / 'ti_solver_page.js').read_text()
    ui = (ROOT / 'web' / 'thermo_ui.js').read_text()
    assert 'onSolve: function (fn) { LISTENERS.push(fn); }' in page
    assert 'last: function () { return LAST; }' in page
    assert 'LISTENERS.forEach' in page
    assert 'TB.onSolve(function (R) {' in ui
    assert 'if (swept) {' in ui
    assert 'drawOptimality(TB.last());' in ui


def test_the_figures_disclose_which_reduced_cost_they_draw(html):
    thermo = html.split('<div id="ws-thermo"')[1].split('<div id="ws-defect"')[0]
    flat = re.sub(r'\s+', ' ', thermo)
    for phrase in (
        'per formula unit',
        'positive at every reported solution',
        'changes sign whenever a candidate is less deficient than the host',
        'shaded where two phases coexist, pale where one does',
    ):
        assert phrase in flat, phrase


def test_a_phase_with_no_oxygen_difference_is_not_reported_as_degenerate(html):
    """Math.abs(null) is 0, and that is not a phase boundary.

    A candidate whose oxygen deficiency relative to the host is exactly
    zero has no denominator to normalise by, so per mol O comes back
    null. Both engines guard that before the degeneracy tolerance; the
    KKT table's reading and this figure's payload have to as well, or a
    phase kJ per formula away from appearing is reported as sitting on
    the boundary. The condition is not exotic - it is what the shipped
    deep-Magneli preset produces."""
    r = solve_case(CASES['magneli'], 1200)
    rc = r['inactive_phase_reduced_costs']['TiO2']
    assert rc['per_mol_O_kJ'] is None
    assert rc['oxygen_removed_per_extent'] == 0.0
    assert rc['per_formula_kJ'] > 1.0
    assert 'TiO2' not in r['degenerate_phases']

    page = (ROOT / 'web' / 'ti_solver_page.js').read_text()
    guard = page.index("var normalisable = po !== null && isFinite(po);")
    used = page.index("Math.abs(po) < R.solver_tolerance", guard)
    assert page.index("!normalisable ?", guard) < used, \
        'the reading tests the tolerance before it tests the guard'
    assert 'var normalisable = po !== null && isFinite(po);' in html

    fig = (ROOT / 'web' / 'figures_thermo.js').read_text()
    assert 'degenerate: po !== null && isFinite(po)' in fig


def test_the_frozen_manuscript_plates_are_untouched():
    """The live figures are new drawings, not edits to the plate kit.

    Every plate is re-rendered from the shared kit and held against the
    file already on disk. If a change here moved e1 or e2 - or anything
    else the kit draws - it shows up as a diff in a manuscript figure
    rather than in a proof."""
    assert shutil.which('node') is not None, 'node is required'
    r = subprocess.run(['node', str(ROOT / 'tests' / 'js' / 'plates_harness.js'),
                        str(ROOT)], capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, r.stderr[-2000:]
    out = json.loads(r.stdout)
    assert not out['errors'], out['errors']
    frozen = ROOT / 'docs' / 'figures'
    moved = []
    for fid, svg in out['svg'].items():
        onfile = frozen / (fid + '.svg')
        assert onfile.exists(), fid
        # the exporter ends the file with a newline; the kit does not
        if onfile.read_text().rstrip('\n') != svg.rstrip('\n'):
            moved.append(fid)
    assert not moved, moved
    assert {'e1', 'e2'} <= set(out['svg']), 'the equilibrium plates are gated'
