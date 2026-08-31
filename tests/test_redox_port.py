"""The browser redox and cycling engines against the Python ones.

Two different claims are gated here, and they are not the same claim.

web/cycling.js is a line-for-line port of solidgas/cycling.py, so it is
held to equality: the same schedule, the same fit, the same discriminating
prediction, digit for digit.

web/redox.js is NOT a port of the Python solver. It carries the closed
forms - theta* = K/(1+K), the two window limits - and leaves the
bisection, the descriptor axis and the Bronsted-Evans-Polanyi overlay in
Python where they belong. That is only safe if the closed form and the
solver agree, so the gates below run the real bisection, with the surface
enrichment entered as the average -> surface mapping the solver takes, and
compare it against what the page would draw. A page that quietly grew its
own algebra would fail here.

Nothing may skip: an engine that will not load is a failure.
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

from solidgas import cycling as CY                     # noqa: E402
from solidgas import redox as R                        # noqa: E402

HARNESS = ROOT / 'tests' / 'js' / 'redox_harness.js'


@pytest.fixture(scope='module')
def js():
    assert shutil.which('node') is not None, 'node is required'
    r = subprocess.run(['node', str(HARNESS), str(ROOT)],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-2000:]
    return json.loads(r.stdout)


@pytest.fixture(scope='module')
def p():
    return R.load_params()


def rel(a, b):
    scale = max(abs(a), abs(b), 1e-300)
    return abs(a - b) / scale


# ------------------------------------------------------------------ redox

def test_the_page_solves_the_same_crossing_the_bisection_finds(js, p):
    """Closed form on the page, real root-finding in Python, same answer.

    The enrichment enters the solver the way the module documents it: as a
    callable from the particle average to the reacting face. So this is
    not the closed form compared against itself - it is 455 bisections
    against 455 evaluations of K/(1+K) divided by E."""
    worst_t, worst_s, worst_r, worst_res = 0.0, 0.0, 0.0, 0.0
    for row in js['steady']:
        E = row['E']
        mapping = None if E == 1.0 else (lambda t, e=E: min(1.0, e * t))
        s = R.from_ratio(p, row['K'], theta_of_average=mapping, window_n=0)
        worst_t = max(worst_t, rel(s['theta'], row['theta']))
        worst_s = max(worst_s, rel(s['theta_surface'], row['theta_surface']))
        worst_r = max(worst_r, rel(s['rate_s1'], row['rate_over_k_ox']))
        worst_res = max(worst_res, s['residual'])
        assert s['regime'] == row['regime'], row
        # E = 125 and K = 1 is the exact corner of the diagram: the
        # crossover enrichment and the face ceiling both land on
        # theta = theta_sol, so the two engines are comparing a strict
        # inequality against a tie and either answer is defensible. The
        # occupancy itself is still held to 1e-12 above.
        if abs(s['theta'] - s['theta_solubility']) > 1e-12:
            assert s['beyond_solubility'] == row['beyond_solubility'], row
    assert worst_t < 1e-12, worst_t
    assert worst_s < 1e-12, worst_s
    # The rate is the mean of the two half-rates, i.e. of two numbers whose
    # DIFFERENCE is what the bisection drove to zero, so it inherits the
    # solver's own convergence rather than the occupancy's. Holding it to a
    # round number would be inventing a tolerance; the honest bound is the
    # residual the solver itself reports, worst over the same grid, which
    # at the K = E = 1e4 corner is about 9e-12.
    assert worst_r <= worst_res, (worst_r, worst_res)
    assert worst_res < 1e-10, worst_res


def test_the_closed_form_is_the_one_the_python_module_publishes(js):
    for row in js['steady']:
        assert rel(R.theta_of_ratio(row['K']), row['theta_surface']) < 1e-15


def test_the_design_window_matches_at_every_enrichment(js, p):
    for w in js['windows']:
        ref = R.usable_window(p, w['E'])
        assert w['binding'] == ref['binding'], w
        for key in ('K_max', 'K_second_phase', 'K_no_steady_state'):
            a, b = ref[key], w[key]
            # JSON has no infinity; the harness sends it as a sentinel
            if isinstance(b, str) or math.isinf(a):
                assert math.isinf(a) and b == 'inf', (key, w)
            else:
                assert rel(a, b) < 1e-13, (key, w)


def test_the_binding_limit_swaps_where_the_algebra_says_it_does(js, p):
    """E = theta_max/theta_sol is where widening the phase field stops
    buying anything, and it has to be the same number on both sides."""
    st = p['structure']
    sol = st['x_solubility'] / st['oxygen_per_formula']
    th_max = st.get('theta_surface_max', 0.5)
    assert rel(js['crossover_E'], th_max / sol) < 1e-15
    by_e = {w['E']: w['binding'] for w in js['windows']}
    assert by_e[124] == 'second_phase'
    assert by_e[126] == 'no_steady_state'


def test_the_uniform_assumption_is_the_249_fold_narrowing(js):
    by_e = {w['E']: w for w in js['windows']}
    uniform, resolved = by_e[1]['K_max'], by_e[1844]['K_max']
    assert uniform == pytest.approx(0.0040160642570281124, rel=1e-12)
    assert resolved == pytest.approx(1.0, rel=1e-12)
    assert resolved / uniform == pytest.approx(249.0, abs=0.5)


def test_no_steady_state_is_the_ceiling_and_says_so(js):
    """The flag is a statement about the face, not about the root: the
    crossing above K = 1 exists mathematically and the surface cannot go
    there to realise it."""
    for row in js['steady']:
        assert row['no_steady_state'] == (row['K'] > 1.0), row
        if row['no_steady_state']:
            assert not row['usable']
            assert row['theta_surface'] > 0.5


def test_the_page_draws_the_rows_the_exporter_froze(js):
    """Same figure payload from the browser engine and from Python, so the
    interactive plot and the file in docs/figures cannot diverge."""
    frozen = json.loads(
        (ROOT / 'data' / 'redox_figure_data.json').read_text())
    a, b = frozen['crossing'], js['crossing']
    assert a['K'] == b['K']
    assert rel(a['theta_star'], b['theta_star']) < 1e-15
    assert len(a['rows']) == len(b['rows'])
    for ra, rb in zip(a['rows'], b['rows']):
        for key in ('theta', 'r_red', 'r_ox'):
            assert rel(ra[key], rb[key]) < 1e-15, key

    a, b = frozen['phase'], js['phase']
    for key in ('theta_sol', 'theta_surface_max', 'K_no_steady_state',
                'crossover_E'):
        assert rel(a[key], b[key]) < 1e-13, key
    assert len(a['second_phase']) == len(b['second_phase'])
    for ra, rb in zip(a['second_phase'], b['second_phase']):
        assert rel(ra['E'], rb['E']) < 1e-13
        assert rel(ra['K'], rb['K']) < 1e-13


# ---------------------------------------------------------------- cycling

def test_the_cycle_fit_is_the_same_fit(js):
    ref = CY.fit([{'created': 94.8, 'recovered': 81.5},
                  {'created': 61.0, 'recovered': 51.6}])
    got = js['cycling']['measured']
    assert rel(ref['f_rec'], got['f_rec']) < 1e-15
    assert rel(ref['lock'], got['lock']) < 1e-15
    assert rel(ref['rms_umol_g'], got['rms_umol_g']) < 1e-13
    assert ref['lock_at_bound'] == got['lock_at_bound']
    assert len(ref['notes']) == len(got['notes'])
    for a, b in zip(ref['residual_umol_g'], got['residual_umol_g']):
        assert abs(a - b) < 1e-12


def test_the_schedule_runs_the_same_way(js):
    ref = CY.run([100, 80, 60, 50, 40], 0.7, 0.35)
    got = js['cycling']['synthetic_rows']
    assert len(ref) == len(got)
    for a, b in zip(ref, got):
        for key in ('recovered_umol_g', 'residual_umol_g',
                    'accessible_umol_g', 'locked_umol_g',
                    'pre_accessible_umol_g', 'pre_locked_umol_g'):
            assert rel(a[key], b[key]) < 1e-14, key


def test_the_browser_fit_inverts_the_browser_run(js):
    """Five cycles generated at (0.7, 0.35) and fitted back in the browser
    engine alone: if the two disagreed, the parity above would only prove
    that both are wrong in the same way."""
    back = js['cycling']['synthetic_fit']
    assert back['f_rec'] == pytest.approx(0.7, rel=1e-12)
    assert back['lock'] == pytest.approx(0.35, rel=1e-12)
    assert back['rms_umol_g'] < 1e-12


def test_the_discriminating_experiment_predicts_the_same_gap(js):
    for key, split in (('two_sample', 2), ('two_sample_4', 4)):
        got = js['cycling'][key]
        ref = CY.two_sample_prediction(got['total_umol_g'], got['f_rec'],
                                       got['lock'], split=split)
        assert ref['cycles_to_build'] == got['cycles_to_build']
        assert rel(ref['recovery_gap_umol_g'],
                   got['recovery_gap_umol_g']) < 1e-13
        assert rel(ref['recovery_ratio'], got['recovery_ratio']) < 1e-13
        assert rel(ref['cycled']['locked'], got['cycled']['locked']) < 1e-13


def test_switching_the_memory_off_reproduces_the_memoryless_model(js):
    """lock = 0 has to be the one-scalar model exactly, in the browser as
    in Python, or the memory is not a switch but a new assumption."""
    r = subprocess.run(
        ['node', '-e',
         'const C=require(process.argv[1]);'
         'const rows=C.run([100,80,60],0.6,0);'
         'const t=C.twoSamplePrediction(150,0.6,0,3);'
         'process.stdout.write(JSON.stringify({rows:rows,two:t}));',
         str(ROOT / 'web' / 'cycling.js')],
        capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-2000:]
    got = json.loads(r.stdout)
    ref = CY.run([100, 80, 60], 0.6, 0.0)
    for a, b in zip(ref, got['rows']):
        assert rel(a['recovered_umol_g'], b['recovered_umol_g']) < 1e-14
        assert b['locked_umol_g'] == 0.0
    assert got['two']['recovery_gap_umol_g'] == pytest.approx(0.0, abs=1e-12)


# --------------------------------------------------------------- the page

DOCS = ROOT / 'docs'
WEB = ROOT / 'web'


def test_the_page_inlines_the_redox_stack():
    """A workspace that fetches anything is not self-contained."""
    html = (DOCS / 'ti_solver.html').read_text()
    assert (ROOT / 'data' / 'redox.json').read_text().strip() in html, \
        'page is stale against redox.json - run build_ti_solver.py'
    for src in (WEB / 'redox.js', WEB / 'cycling.js',
                WEB / 'figures_redox.js', WEB / 'redox_ui.js'):
        assert src.read_text() in html, \
            f'page is stale against {src.name} - run build_ti_solver.py'
    for token in ('__RX_DATA__', '__RX_ENGINE__', '__CY_ENGINE__',
                  '__RX_FIGS__', '__RX_UI__'):
        assert token not in html


def test_the_workspace_is_reachable_and_says_what_it_assumes():
    template = (WEB / 'ti_solver_template.html').read_text()
    portal = (WEB / 'portal.html').read_text()
    index = (DOCS / 'index.html').read_text()

    # a workspace nothing links to is not shipped
    assert 'data-ws="ws-redox"' in template
    assert 'id="ws-redox"' in template
    assert 'id="head-redox"' in template
    assert "h.indexOf('redox')" in template, 'no hash route to the workspace'
    assert 'ti_solver.html#redox' in portal, 'the portal does not link it'
    assert portal == index, 'run build_portal.py'
    assert '03 / 03' in portal and '01 / 03' in portal

    # the controls that make it a workspace rather than a picture
    for ident in ('id="rxK"', 'id="rxKs"', 'id="rxE"', 'id="rxSol"',
                  'id="rxMax"', 'id="figCrossing"', 'id="figPhase"',
                  'id="cyIn"', 'id="cyOut"', 'id="cyTotal"',
                  'id="cySplit"'):
        assert ident in template, f'{ident} missing'

    # And the assumptions, in the page's own words. Two haystacks: the
    # built HTML for the prose, and the UI source with its adjacent string
    # literals rejoined, because a sentence the page assembles at runtime
    # is still a sentence the reader gets.
    ui = (WEB / 'redox_ui.js').read_text()
    hay = ' '.join((DOCS / 'ti_solver.html').read_text().split()).lower() \
        + ' ' + ' '.join(re.sub(r"'\s*\+\s*'", '', ui).split()).lower()
    for phrase in ('unit site orders',
                   'k = k(reduction) / k(oxidation)',
                   'only the <i>ratio</i>',
                   'sets to 1 without saying so',
                   'the turnover at a given k is <i>unchanged</i>',
                   'nothing here is methodologically new',
                   'no transport limitation',
                   'carried as a sensitivity parameter',
                   'lock = 0 is the memoryless model exactly',
                   'two cycles constrain two parameters and no more',
                   'a memoryless inventory predicts a gap of',
                   'cycle integrals alone do not separate them'):
        assert ' '.join(phrase.split()) in hay, \
            f'disclosure missing: {phrase!r}'


def test_the_interactive_phase_map_marks_the_current_point_and_the_file_does_not():
    """The 'here' ring is the one thing the page adds to the figure, and
    the exported SVG must not grow it - the manuscript figure is frozen."""
    figs = (WEB / 'figures_redox.js').read_text()
    assert 'd.current' in figs
    ui = (WEB / 'redox_ui.js').read_text()
    assert 'current: { K: q.K, E: q.E, usable: ok }' in ui
    exported = (ROOT / 'docs' / 'figures' / 'redox_phase.svg').read_text()
    assert '>here<' not in exported
    frozen = json.loads(
        (ROOT / 'data' / 'redox_figure_data.json').read_text())
    assert 'current' not in frozen['phase']


# ------------------------------------------------- what the model is

def _prose():
    """Page, docs and engine as one haystack, with the page's JS string
    literals rejoined so a phrase split across a concatenation still
    counts as written."""
    page = (ROOT / 'docs' / 'ti_solver.html').read_text()
    ui = re.sub(r"'\s*\+\s*'", '', (WEB / 'redox_ui.js').read_text())
    doc = (ROOT / 'docs' / 'redox-steady-state.md').read_text()
    eng = (ROOT / 'solidgas' / 'redox.py').read_text()
    data = (ROOT / 'data' / 'redox.json').read_text()
    return {'page': ' '.join((page + ui).split()),
            'doc': ' '.join(doc.split()),
            'engine': ' '.join(eng.split()),
            'data': ' '.join(data.split())}


def test_the_mechanism_is_named():
    """Mars-van Krevelen is the one phrase that tells a reaction engineer
    what this is in a single word, and it was missing from every file."""
    p = _prose()
    assert 'Mars&ndash;van Krevelen' in p['page'] or 'Mars-van Krevelen' in p['page']
    assert 'Mars–van Krevelen' in p['doc']
    assert 'Mars-van Krevelen' in p['engine']
    assert 'Mars-van Krevelen' in p['data']
    for where in ('page', 'doc', 'engine'):
        assert 'no surface intermediates' in p[where], where


def test_theta_says_which_face_and_which_oxygen():
    """A vacancy fraction with no site attached is not a coverage.

    The site is not a free choice: the descriptor is referenced to (110)
    bridging O and theta_max = 0.5 is the (1x2) added-row phase, which is
    a bridging-row reconstruction. Both numbers are already in the
    dataset, so theta has to be that site or they do not belong in the
    same equation - and the in-plane oxygen has to be excluded out loud,
    because averaging two removal energies describes neither site."""
    p = _prose()
    params = R.load_params()
    assert 'bridging O' in params['descriptor']['reference_material']
    assert '(110)' in params['descriptor']['reference_material']
    for where in ('page', 'doc', 'engine', 'data'):
        assert 'bridging oxygen rows of rutile (110)' in p[where] \
            or 'BRIDGING OXYGEN ROWS of\nrutile (110)'.replace('\n', ' ') \
            in p[where], where
    for where in ('page', 'doc', 'engine'):
        assert 'in-plane oxygen is' in p[where].replace('<i>', '') \
            .replace('</i>', '').replace('**', ''), where
    # and the particle number is named as a composition, not a defect count
    for where in ('page', 'doc', 'engine'):
        assert 'titanium interstitial' in p[where], where


def test_a_fixed_enrichment_is_declared_as_a_timescale_assumption():
    """Holding E constant says shell and bulk stay in equilibrium
    partition while the carrier turns over. If that fails the problem is
    transport into the particle, not a boundary condition on it, and the
    page has to say so rather than leave it as geometry."""
    p = _prose()
    for where in ('page', 'doc'):
        assert 'equilibrium partition' in p[where], where
        assert 'fast compared with both half-reactions' in p[where], where
        assert 'transport into the particle' in p[where], where


def test_the_gases_stay_unnamed_and_the_page_says_why():
    """A deliberate omission, pinned so it is not quietly filled in.

    Both half-reactions carry identical placeholder rate parameters, so a
    named pair would assert chemistry the numbers do not hold. What the
    model does require of the pair - unit order in gas and in site - is
    stated instead."""
    params = R.load_params()
    red = params['half_reactions']['reduction']
    ox = params['half_reactions']['oxidation']
    for key in ('log10_A_s1', 'E0_eV', 'order_gas', 'order_site'):
        assert red[key] == ox[key], key
    assert red['label'] == ox['label'] == 'assumed'
    assert red['stoichiometry'].startswith('R(g)')
    assert ox['stoichiometry'].startswith('OX(g)')
    p = _prose()
    for where in ('page', 'doc', 'engine', 'data'):
        assert 'unnamed' in p[where], where
    for where in ('page', 'doc', 'engine'):
        assert 'unit order in gas and in site' in p[where], where


def test_the_crossing_figure_says_where_its_axis_stops_being_the_material(js):
    """theta on that axis is the face, so the particle's single-phase
    limit lands at theta_sol * E and moves with the enrichment: at E = 1
    it is at 0.004 and the crossing at theta* = 1/3 is inside the second
    phase, at E = 1844 it leaves the axis entirely. The face ceiling does
    not move. Both are the phase map's own limits, in the phase map's own
    colours, seen along one line."""
    params = R.load_params()
    st = params['structure']
    sol = st['x_solubility'] / st['oxygen_per_formula']
    th_max = st.get('theta_surface_max', 0.5)

    for key, E in (('crossing', 1.0), ('crossing_resolved', 1844.0)):
        c = js[key]
        assert rel(c['theta_sol'], sol) < 1e-15, key
        assert rel(c['theta_surface_max'], th_max) < 1e-15, key
        assert c['enrichment'] == E
        assert rel(c['theta_face_at_solubility'], min(1.0, sol * E)) < 1e-15

    assert js['crossing']['theta_face_at_solubility'] < \
        js['crossing']['theta_star'], 'the uniform crossing must be past it'
    assert js['crossing_resolved']['theta_face_at_solubility'] == 1.0, \
        'the layer-resolved reading takes the limit off the axis'

    figs = (WEB / 'figures_redox.js').read_text()
    assert 'd.theta_face_at_solubility' in figs
    assert "bandLabel(xSol, 'second phase', C.bulk)" in figs
    assert "bandLabel(xMax, 'no steady state', C.subsurface)" in figs
    # the same two tints the phase map uses for the same two statements
    assert figs.count('tint(C.bulk, 0.13)') == 2
    assert figs.count('tint(C.subsurface, 0.13)') == 2

    exported = (ROOT / 'docs' / 'figures' / 'redox_crossing.svg').read_text()
    assert '>second phase<' in exported
    assert '>no steady state<' in exported
    frozen = json.loads(
        (ROOT / 'data' / 'redox_figure_data.json').read_text())
    assert frozen['crossing']['enrichment'] == 1.0, \
        'the exported crossing is the uniform reading'
    assert rel(frozen['crossing']['theta_face_at_solubility'], sol) < 1e-15

    hay = ' '.join((ROOT / 'docs' / 'ti_solver.html').read_text().split())
    assert 'a different compound, not a reduced rutile' in hay
