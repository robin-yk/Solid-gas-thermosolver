"""Gates for the steady-state redox model.

The model is one nonlinear equation, so almost every gate here is the
shipped numeric path checked against a closed form that is derivable by
hand. Where no closed form exists the gate is a structural property -
monotonicity, uniqueness, invariance under reparameterisation - rather than
a frozen number, so the test says why the answer is right and not merely
what it was last time."""

import json
import math
import os
import pathlib

import pytest

from solidgas import redox as R
from solidgas import statmech as S

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, 'data')

DEVS = (-1.5, -0.8, -0.3, -0.05, 0.0, 0.05, 0.3, 0.8, 1.5)
TEMPS = (400.0, 600.0, 900.0)
PRESSURES = ((0.02, 0.5), (0.1, 0.1), (0.8, 0.05))


@pytest.fixture(scope='module')
def p():
    return R.load_params()


@pytest.fixture(scope='module')
def rutile_surface():
    """theta_surface(theta_average) for a 0.9 um rutile particle at 600 C.

    The layer-resolved partition read as the mapping the redox model needs:
    give it how reduced the particle is overall, it returns how reduced the
    reacting face is."""
    d = json.load(open(os.path.join(DATA, 'rutile_dft.json')))
    geom = S.geometry(d, d_um=0.90, ss_layers=1)
    n_o = geom['N_O_total']

    def theta_surface(t_bar):
        vo = max(0.0, t_bar) * n_o
        return S.distribute(d, 600.0, vo, geom)['theta']['surface']
    theta_surface.geom = geom
    return theta_surface


# ------------------------------------------------------- the closed forms

def test_the_bisection_reproduces_the_first_order_crossing(p):
    """theta* = k_red / (k_red + k_ox) at unit site orders, everywhere.

    The solver never takes this shortcut - it bisects on every call - so
    agreement here is the numeric path being checked, not restated."""
    assert p['half_reactions']['reduction']['order_site'] == 1.0
    assert p['half_reactions']['oxidation']['order_site'] == 1.0
    worst = 0.0
    for t_c in TEMPS:
        for p_red, p_ox in PRESSURES:
            for dev in DEVS:
                s = R.steady_state(p, t_c=t_c, dev=dev,
                                   p_red=p_red, p_ox=p_ox, window_n=0)
                want = R.theta_first_order(s['k_red_s1'], s['k_ox_s1'])
                worst = max(worst, abs(s['theta'] - want))
    assert worst < 1e-12, 'bisection drifts from the closed form: %g' % worst


def test_the_two_half_rates_are_equal_at_the_crossing(p):
    for t_c in TEMPS:
        for dev in DEVS:
            s = R.steady_state(p, t_c=t_c, dev=dev, window_n=0)
            assert s['interior_crossing']
            assert s['residual'] < 1e-12, (dev, t_c, s['residual'])


def test_the_volcano_peaks_where_the_closed_form_says(p):
    """The swept maximum sits at the analytic Sabatier descriptor.

    dev* = [ (E_ox0 - E_red0) - kT ln(a/b) + kT ln(C_r/C_o) ] / (a + b)."""
    for t_c in TEMPS:
        for p_red, p_ox in PRESSURES:
            opt = R.sabatier_optimum(p, t_c=t_c, p_red=p_red, p_ox=p_ox)
            rows = R.descriptor_sweep(p, t_c=t_c, p_red=p_red, p_ox=p_ox,
                                      n=4001, dev_lo=opt['dev_eV'] - 1.0,
                                      dev_hi=opt['dev_eV'] + 1.0, window_n=0)
            peak = max(rows, key=lambda r: r['rate_s1'])
            step = 2.0 / 4000.0
            assert abs(peak['dev_eV'] - opt['dev_eV']) <= step, \
                (t_c, p_red, p_ox, peak['dev_eV'], opt['dev_eV'])


def test_the_peak_coverage_is_set_by_the_two_slopes_alone(p):
    """theta at the volcano peak is alpha / (alpha + beta) - no T, no p, no A.

    This is the reason the model can be swept before any material is fitted:
    the optimum coverage is a property of the two linear free-energy slopes,
    and only the descriptor at which it occurs moves with the conditions."""
    a = p['half_reactions']['reduction']['bep_slope']
    b = p['half_reactions']['oxidation']['bep_slope']
    want = a / (a + b)
    for t_c in TEMPS:
        for p_red, p_ox in PRESSURES:
            opt = R.sabatier_optimum(p, t_c=t_c, p_red=p_red, p_ox=p_ox)
            assert opt['theta'] == pytest.approx(want, abs=1e-15)
            s = R.steady_state(p, t_c=t_c, dev=opt['dev_eV'],
                               p_red=p_red, p_ox=p_ox, window_n=0)
            assert s['theta'] == pytest.approx(want, abs=1e-12)


def test_a_skewed_slope_pair_moves_the_optimum_off_half(p):
    """alpha != beta must move the peak coverage, or the gate above is vacuous."""
    q = json.loads(json.dumps(p))
    q['half_reactions']['reduction']['bep_slope'] = 0.8
    q['half_reactions']['oxidation']['bep_slope'] = 0.2
    opt = R.sabatier_optimum(q)
    assert opt['theta'] == pytest.approx(0.8, abs=1e-15)
    rows = R.descriptor_sweep(q, n=4001, dev_lo=opt['dev_eV'] - 1.0,
                              dev_hi=opt['dev_eV'] + 1.0, window_n=0)
    peak = max(rows, key=lambda r: r['rate_s1'])
    assert peak['theta'] == pytest.approx(0.8, abs=1e-3)


# ------------------------------------------------------------- uniqueness

def test_one_curve_falls_while_the_other_rises(p):
    """Why the crossing is unique: strict monotonicity in opposite senses."""
    for dev in DEVS:
        rows = R.crossing_curves(p, dev=dev, n=401)
        for a, b in zip(rows, rows[1:]):
            assert b['r_red_s1'] <= a['r_red_s1']
            assert b['r_ox_s1'] >= a['r_ox_s1']
        assert rows[0]['r_red_s1'] > rows[-1]['r_red_s1']
        assert rows[-1]['r_ox_s1'] > rows[0]['r_ox_s1']


def test_the_crossing_curves_meet_at_the_solved_root(p):
    """The plot and the solver are looking at the same two functions."""
    for dev in (-0.4, 0.0, 0.4):
        s = R.steady_state(p, dev=dev, window_n=0)
        k = R.rate_constants(p, dev=dev)
        r_red, r_ox = R.half_rates(p, s['theta'], k)
        assert r_red == pytest.approx(r_ox, rel=1e-12)
        assert r_red == pytest.approx(s['rate_s1'], rel=1e-12)


def test_higher_site_orders_still_cross_once(p):
    q = json.loads(json.dumps(p))
    q['half_reactions']['reduction']['order_site'] = 2.0
    q['half_reactions']['oxidation']['order_site'] = 2.0
    for dev in DEVS:
        s = R.steady_state(q, dev=dev, window_n=0)
        assert 0.0 < s['theta'] < 1.0
        assert s['residual'] < 1e-12
    assert R.sabatier_optimum(q) is None, \
        'the closed form must decline to answer away from unit orders'


def test_a_zero_site_order_reports_a_runaway_rather_than_a_root(p):
    """With nothing to slow the reduction there is no steady state, and the
    model must say so instead of returning the bracket end as an answer."""
    q = json.loads(json.dumps(p))
    q['half_reactions']['reduction']['order_site'] = 0.0
    q['half_reactions']['oxidation']['bep_slope'] = 0.0
    q['half_reactions']['oxidation']['E0_eV'] = 3.0
    s = R.steady_state(q, dev=0.0)
    assert s['theta'] == 1.0
    assert s['interior_crossing'] is False
    assert s['runaway'] == 'reduction'
    assert not s['localised']
    assert any(w['kind'] == 'runaway_reduction' for w in s['warnings'])


# ------------------------------------------------------------ the domain

def test_the_solubility_flag_never_clips_the_answer(p):
    """Leaving the point-defect domain is reported, not enforced."""
    for sol in (0.001, 0.004, 0.4):
        s = R.steady_state(p, dev=0.0, theta_sol=sol, window_n=0)
        assert s['theta'] == pytest.approx(0.5, abs=1e-12), \
            'theta must not depend on the declared solubility'
        assert s['beyond_solubility'] is (0.5 > sol)
        flagged = any(w['kind'] == 'beyond_solubility' for w in s['warnings'])
        assert flagged is (0.5 > sol)


def test_the_uniform_optimum_demands_a_bulk_no_oxide_sustains(p):
    """Why the uniform coupling fails, stated without reference to a material.

    With symmetric slopes the volcano peaks at theta = 1/2, and in the
    uniform model theta is simultaneously the surface coverage and the bulk
    composition. Half the oxygen sublattice vacant is a coverage a surface
    reaches routinely - the rutile (1x2) row phase is exactly that - and a
    bulk composition no oxide holds as a solid solution. So the uniform
    model puts its own optimum outside the domain of every carrier, which
    is a defect of forcing one number to be both quantities rather than a
    finding about oxides."""
    opt = R.sabatier_optimum(p)
    peak = R.steady_state(p, dev=opt['dev_eV'], window_n=0)
    assert peak['theta'] == pytest.approx(0.5, abs=1e-12)
    generous = 0.125          # a fluorite-width delta, far past rutile
    for sol in list(p['structure']['x_solubility_range']) + [generous]:
        assert R.steady_state(p, dev=opt['dev_eV'], theta_sol=sol,
                              window_n=0)['beyond_solubility'], sol

    rows = [r for r in R.descriptor_sweep(p, n=3001, window_n=0)
            if not r['beyond_solubility']]
    assert rows, 'some of the axis must stay inside the domain'
    best_ok = max(rows, key=lambda r: r['rate_s1'])
    assert best_ok['rate_s1'] < peak['rate_s1'], \
        'the reachable best must sit below the unreachable peak'


def test_the_layer_resolved_model_puts_the_optimum_back_in_reach(
        p, rutile_surface):
    """And why separating the two quantities resolves it.

    Once the surface is allowed to differ from the average, the carrier can
    run its face at the optimum coverage while the bulk stays a dilute
    solid solution. Just past the balance point the rutile surface sits near
    the (1x2) deficiency while the particle average is three orders of
    magnitude lower - inside the declared point-defect solubility, with no
    flag raised. The uniform model cannot represent that state at all."""
    opt = R.sabatier_optimum(p)
    s = R.steady_state(p, dev=opt['dev_eV'] + 0.05,
                       theta_of_average=rutile_surface)
    sol = p['structure']['x_solubility'] / p['structure']['oxygen_per_formula']
    assert s['theta_surface'] > 40.0 * sol, 'the face must do the work'
    assert s['theta'] < sol, 'while the bulk stays a solid solution'
    assert not s['beyond_solubility']
    uniform = R.steady_state(p, dev=opt['dev_eV'] + 0.05, window_n=0)
    assert uniform['beyond_solubility'], 'the same point is out of domain uniform'
    assert uniform['theta'] == pytest.approx(s['theta_surface'], rel=1e-9), \
        'the two models agree on the face and differ on the particle'


def test_the_bep_relation_flags_its_own_extrapolation(p):
    """A negative barrier is reported, never quietly floored."""
    k = R.rate_constants(p, dev=-4.0)
    assert k['E_red_eV'] < 0.0
    assert any(w['kind'] == 'bep_out_of_range' for w in k['warnings'])
    assert k['k_red_s1'] > 0.0


# ---------------------------------------------- the layer-resolved correction

def test_the_identity_mapping_is_the_uniform_model_exactly(p):
    c = R.surface_correction(p, lambda t: t, dev=0.4)
    assert c['rate_ratio'] == pytest.approx(1.0, abs=1e-15)
    assert c['theta_ratio'] == pytest.approx(1.0, abs=1e-15)
    assert c['corrected']['theta'] == pytest.approx(
        c['uniform']['theta'], abs=1e-15)


def test_a_monotone_surface_mapping_leaves_the_steady_rate_untouched(p):
    """The invariance that says which half of the correction matters.

    Any strictly increasing average -> surface mapping reparameterises the
    same crossing in surface coverage, and that coverage is fixed by the
    rate constants alone. So the turnover cannot move; only the inventory
    holding the surface there can. Checked against a mapping that has
    nothing to do with rutile, because this is algebra, not chemistry."""
    for gamma in (0.5, 2.0, 7.0):
        c = R.surface_correction(p, lambda t, g=gamma: t ** g, dev=0.4)
        assert c['rate_ratio'] == pytest.approx(1.0, rel=1e-10), gamma
        assert c['corrected']['theta_surface'] == pytest.approx(
            c['uniform']['theta'], rel=1e-10), gamma


def test_the_rutile_surface_runs_far_above_the_particle_average(p, rutile_surface):
    """The enrichment factor, and its geometric ceiling.

    In the dilute limit every class sits on its own isotherm at one mu, so
    the average is the surface diluted by the Boltzmann-weighted site
    count, and the factor is bounded by N_O_total / N_s. It would reach
    that ceiling only if the deeper classes froze out entirely.

    They do not, and that is what changed. A compensated vacancy has
    three configurational logarithms, so the surface-to-bulk ratio is
    exp(dE/3kT) = 333 rather than exp(dE/kT) = 3.7e7 - the cube root -
    and against 1846 times as many bulk sites the deep classes carry most
    of the inventory even in the dilute limit. The factor is 270, not the
    1844 the neutral-vacancy form gave, and it sits well below the
    ceiling instead of just under it."""
    geom = rutile_surface.geom
    ceiling = geom['N_O_total'] / geom['N_s']
    assert ceiling == pytest.approx(1847.0, abs=1.0)
    for t_bar in (1e-7, 1e-6, 1e-5):
        factor = rutile_surface(t_bar) / t_bar
        assert 250.0 < factor < 280.0, (t_bar, factor)
        assert factor < 0.2 * ceiling, (t_bar, factor)


def test_the_correction_moves_the_state_by_two_orders_not_the_rate(
        p, rutile_surface):
    """Both halves of the headline, on the real depth-weighted mapping.

    The rate is untouched exactly - that is the point, and it is exact
    because the correction rescales theta and the two half-reactions
    rescale with it. The state moves by the enrichment, which is 270 in
    the dilute limit and falls below it as the surface starts to
    saturate, so the check is a range with a direction rather than one
    number."""
    moved_by_dev = []
    for dev in (0.1, 0.3, 0.6, 1.0):
        c = R.surface_correction(p, rutile_surface, dev=dev)
        assert c['rate_ratio'] == pytest.approx(1.0, rel=1e-9), dev
        moved = c['uniform']['theta'] / c['corrected']['theta']
        assert 100.0 < moved < 275.0, (dev, moved)
        moved_by_dev.append(moved)
    # the strongly reducing end is where the surface is furthest up its
    # isotherm and can supply the least extra enrichment; as the
    # descriptor climbs the state thins out and the factor relaxes onto
    # its dilute limit from below
    assert all(y >= x * 0.999 for x, y in
               zip(moved_by_dev, moved_by_dev[1:])), moved_by_dev
    assert moved_by_dev[0] < 0.7 * moved_by_dev[-1]
    assert moved_by_dev[-1] == pytest.approx(270.0, abs=3.0)


def test_a_pinned_surface_removes_the_steady_state_on_half_the_axis(
        p, rutile_surface):
    """The structural change the uniform model cannot see.

    The rutile surface saturates at the (1x2) line-phase deficiency, so
    above a small inventory neither half-rate depends on how reduced the
    particle is. Every oxide easier to reduce than the balance point then
    has no crossing at all: it reduces until something outside this pair of
    rate laws stops it. The boundary is the Sabatier descriptor itself,
    because the surface ceiling and the peak coverage are both 1/2 here.
    The uniform model reports an ordinary crossing at every one of these
    points."""
    opt = R.sabatier_optimum(p)
    assert R.steady_state(
        p, dev=0.0, theta_of_average=rutile_surface)['surface_buffered'], \
        'the rutile mapping must be recognised as buffered'

    for dev in (opt['dev_eV'] - 0.5, opt['dev_eV'] - 0.05):
        s = R.steady_state(p, dev=dev, theta_of_average=rutile_surface)
        assert s['runaway'] == 'reduction', dev
        assert not s['localised'], dev
        assert any(w['kind'] == 'surface_buffered' for w in s['warnings'])
        assert R.steady_state(p, dev=dev, window_n=0)['localised'], \
            'the uniform model must see nothing wrong here'

    for dev in (opt['dev_eV'] + 0.05, opt['dev_eV'] + 0.5):
        s = R.steady_state(p, dev=dev, theta_of_average=rutile_surface)
        assert s['runaway'] is None, dev
        assert s['localised'], dev


def test_the_degenerate_point_is_reported_rather_than_resolved(
        p, rutile_surface):
    """Exactly at the balance point every buffered inventory is a steady
    state, and the solver must say so instead of returning the near edge as
    though it were the answer."""
    s = R.steady_state(p, dev=R.sabatier_optimum(p)['dev_eV'],
                       theta_of_average=rutile_surface)
    assert s['interior_crossing'], 'a root is found'
    assert not s['localised'], 'but it does not locate the state'
    assert any(w['kind'] == 'root_not_localised' for w in s['warnings'])


# ---------------------------------------------------------------- dataset

def test_every_parameter_block_declares_where_it_came_from(p):
    for name in ('reduction', 'oxidation'):
        hr = p['half_reactions'][name]
        assert hr['label'] in ('assumed', 'literature', 'measured', 'fitted')
        assert hr['note']
        for key in ('order_gas', 'order_site', 'log10_A_s1', 'E0_eV',
                    'bep_slope'):
            assert isinstance(hr[key], float), (name, key)
    d = p['descriptor']
    assert 'DOI' in d['reference_provenance']
    assert d['offset_range_eV'][0] < 0.0 < d['offset_range_eV'][1]
    st = p['structure']
    assert st['x_solubility_range'][0] <= st['x_solubility'] \
        <= st['x_solubility_range'][1]
    assert 'shear' in st['note'].lower()


def test_the_descriptor_reference_matches_the_rutile_dataset(p):
    """The absolute anchor is the same number the defect engine cites."""
    d = json.load(open(os.path.join(DATA, 'rutile_dft.json')))
    prov = d['vacancy_energetics']['provenance']
    assert 'sX 4.39' in prov
    assert p['descriptor']['reference_eV'] == 4.39
    assert 'Li, Guo & Robertson' in p['descriptor']['reference_provenance']
    # the functional spread is the uncertainty on placing a material
    assert '0.68 eV' in p['descriptor']['note']
    assert abs(4.39 - 3.71 - 0.68) < 1e-12


# --------------------------------------------------------------- the writeup

def test_the_documented_claims_are_the_ones_the_code_makes():
    """The three results in docs/redox-steady-state.md, re-derived here.

    A document that drifts from the engine is worse than no document, so
    the numbers it quotes are computed rather than transcribed."""
    doc = (pathlib.Path(ROOT) / 'docs' / 'redox-steady-state.md').read_text()
    assert '269x above the particle average' in doc
    assert 'N_O_total / N_s = 1847' in doc
    assert 'The steady rate does not move' in doc
    assert 'the steady state stops existing' in doc
    assert '0.68 eV spread' in doc

    # the enrichment the document quotes, recomputed
    d = json.load(open(os.path.join(DATA, 'rutile_dft.json')))
    g = S.geometry(d, d_um=0.90, ss_layers=1)
    n_o = g['N_O_total']
    e = S.distribute(d, 600.0, 1e-6 * n_o, g)['theta']['surface'] / 1e-6
    assert ('%.0fx above the particle average' % e) in doc
    # and the margin the window conclusion now rests on
    th_sol = d['saturation']['x_max_shear'] / 2.0
    assert th_sol * e > 1.0
    assert ('%.2f' % (th_sol * e)) in doc
    # and the model it describes is the one that ships
    p = R.load_params()
    a = p['half_reactions']['reduction']['bep_slope']
    b = p['half_reactions']['oxidation']['bep_slope']
    assert R.sabatier_optimum(p)['theta'] == pytest.approx(a / (a + b))


def test_the_defect_model_states_what_it_does_not_carry():
    """Three limits that a reader would otherwise have to infer."""
    doc = (pathlib.Path(ROOT) / 'docs' / 'defect-model.md').read_text()
    for claim in ('Henderson, Surf. Sci. 419, 174 (1999)',
                  'Titanium interstitials',
                  'cannot distinguish a surface-concentrated inventory',
                  'not a polaron hopping rate',
                  'docs/redox-steady-state.md'):
        assert claim in doc, claim


# ----------------------------------- robustness to the shell thickness

@pytest.mark.parametrize('ss_layers', [1, 2, 4])
def test_the_redox_conclusions_survive_the_shell_thickness(p, ss_layers):
    """The least-constrained knob left in the defect model changes nothing
    here.

    The declared thickness of the subsurface class acts below the bridging
    plane and matters at the inventories the measurements sit at. The
    redox steady state does not sit there: it sits at theta_bar of order
    1e-5, deep in the dilute regime, where the deeper classes are
    Boltzmann-suppressed. So the enrichment is set by the energies and the
    site counts, not by how the shell is sliced, and the three claims that
    rest on it do not move.

    The enrichment itself moved once, and by a lot: 1844 to 270. That is
    the compensating polarons entering the free energy, which turns the
    dilute exponent from exp(dE/kT) into exp(dE/3kT). The claims below are
    the ones that survive it, and they survive because they are about
    where the crossing sits and not about how large the factor is."""
    d = json.load(open(os.path.join(DATA, 'rutile_dft.json')))
    geom = S.geometry(d, d_um=0.90, ss_layers=ss_layers)
    n_o = geom['N_O_total']

    def theta_surface(t_bar):
        return S.distribute(d, 600.0, max(0.0, t_bar) * n_o,
                            geom)['theta']['surface']

    # the enrichment, and its geometric ceiling
    ceiling = n_o / geom['N_s']
    factor = theta_surface(1e-6) / 1e-6
    # the shell thickness moves the factor by a tenth at most: it trades
    # subsurface sites against bulk ones, and both are far below the row
    assert 230.0 < factor < 275.0, factor
    assert factor < ceiling

    # the rate is still untouched, and the state still moves by the factor
    c = R.surface_correction(p, theta_surface, dev=0.3)
    assert c['rate_ratio'] == pytest.approx(1.0, rel=1e-9)
    moved = c['uniform']['theta'] / c['corrected']['theta']
    assert moved == pytest.approx(factor, rel=0.05), moved

    # and the runaway boundary is still the Sabatier descriptor
    opt = R.sabatier_optimum(p)['dev_eV']
    assert R.steady_state(p, dev=opt - 0.05,
                          theta_of_average=theta_surface)['runaway'] \
        == 'reduction'
    assert R.steady_state(p, dev=opt + 0.05,
                          theta_of_average=theta_surface)['runaway'] is None


def test_the_shell_thickness_does_move_the_measurement_regime(p):
    """The sensitivity check above would be vacuous if the knob were inert.

    At the inventory the experiments sit at it moves the split a lot; the
    point of the test above is that the redox steady state is not at that
    inventory.

    What no longer moves it is resolving the shell into separate layers
    with a declared depth profile, and that mechanism has been removed:
    with the compensation in the free energy it changed the bulk by two
    percent, where before it changed it by thirty. The thickness of the
    single shell still matters, because it moves sites rather than
    energies."""
    d = json.load(open(os.path.join(DATA, 'rutile_dft.json')))
    thin = S.distribute(d, 600.0, 95.0,
                        S.geometry(d, d_um=0.90, ss_layers=1))
    thick = S.distribute(d, 600.0, 95.0,
                         S.geometry(d, d_um=0.90, ss_layers=4))
    assert thin['umol_g']['bulk'] / thick['umol_g']['bulk'] > 1.1
    assert abs(thick['mu_V_eV'] - thin['mu_V_eV']) > 0.02


# --------------------------------------------- the assumption-free K axis

def test_the_ratio_entry_point_needs_nothing_but_the_ratio(p):
    """theta* = K/(1+K), to the last digit, with no BEP anywhere in it.

    This is the whole material dependence of the steady state expressed as
    one dimensionless number. Everything else the module offers is a way of
    predicting K, and each of those steps costs an assumed parameter."""
    for kk in (1e-6, 1e-3, 0.004, 0.1, 0.5, 1.0, 2.0, 10.0, 1e4):
        s = R.from_ratio(p, kk)
        assert s['theta'] == pytest.approx(R.theta_of_ratio(kk), abs=1e-15), kk
        assert s['K'] == kk
        assert s['rate_unit'] == 'k_ox'
        assert s['residual'] < 1e-12, kk
        # rate in units of k_ox is theta* itself at unit orders
        assert s['rate_s1'] == pytest.approx(s['theta'], rel=1e-12), kk
    assert R.from_ratio(p, 0.0)['theta'] == 0.0
    with pytest.raises(ValueError):
        R.from_ratio(p, -1.0)


def test_the_ratio_and_descriptor_axes_are_the_same_model(p):
    """Two parameterisations, one crossing.

    Whatever the descriptor route computes for the two rate constants, the
    ratio route reproduces its state exactly when handed their ratio - so
    the BEP layer is an overlay on this axis rather than a different model
    underneath."""
    for dev in (-0.5, -0.1, 0.0, 0.3, 0.9):
        d = R.steady_state(p, dev=dev, window_n=0)
        r = R.from_ratio(p, d['k_red_s1'] / d['k_ox_s1'], window_n=0)
        assert r['theta'] == pytest.approx(d['theta'], abs=1e-12), dev
        assert r['regime'] == d['regime'], dev


def test_the_sweep_is_monotone_and_spans_both_limits(p):
    rows = R.ratio_sweep(p, n=81)
    th = [r['theta'] for r in rows]
    assert all(b > a for a, b in zip(th, th[1:]))
    assert th[0] < 1e-3 and th[-1] > 0.999
    assert {r['regime'] for r in rows} == {
        'reduction_limited', 'balanced', 'oxidation_limited'}


# ------------------------------------------------------ the phase diagram

def test_the_usable_window_opens_with_surface_enrichment(p):
    """What the uniform assumption costs, as a design criterion.

    A carrier is usable while it stays one phase and its face stays below
    the coverage its own phase pins it at. Assuming the particle uniformly
    reduced sets the enrichment to 1 and demands a 250:1 rate asymmetry;
    the layer-resolved enrichment of a 0.9 um rutile particle demands only
    that oxidation keep up with reduction."""
    uni = R.usable_window(p, 1.0)
    lay = R.usable_window(p, 1844.0)
    assert uni['K_max'] == pytest.approx(0.0040161, rel=1e-4)
    assert lay['K_max'] == pytest.approx(1.0, rel=1e-12)
    assert lay['K_max'] / uni['K_max'] > 200.0
    # and the failure mode is not the same failure mode
    assert uni['binding'] == 'second_phase'
    assert lay['binding'] == 'no_steady_state'


def test_the_binding_limit_switches_at_a_computable_enrichment(p):
    """The two limits cross where E*theta_sol = theta_surface_max.

    Below it the solid shears before the face saturates; above it the face
    saturates first and the bulk never gets near the solubility. Which one
    binds decides what a carrier should be engineered against."""
    st = p['structure']
    sol = st['x_solubility'] / st['oxygen_per_formula']
    cross = 0.5 / sol
    assert cross == pytest.approx(125.0, rel=1e-12)
    assert R.usable_window(p, cross * 0.9)['binding'] == 'second_phase'
    assert R.usable_window(p, cross * 1.1)['binding'] == 'no_steady_state'
    assert R.usable_window(p, cross)['K_max'] == pytest.approx(
        R.usable_window(p, cross)['K_no_steady_state'], rel=1e-12)


def test_the_phase_map_agrees_with_the_closed_form_boundaries(p):
    """Every cell of the swept map is on the side the closed form says."""
    m = R.phase_map(p, nk=41, ne=31)
    assert m['K_no_steady_state'] == pytest.approx(1.0, rel=1e-12)
    seen = set()
    lgk0, lgk1 = math.log10(m['K_lo']), math.log10(m['K_hi'])
    for row in m['rows']:
        e = row['E']
        kmax = R.usable_window(p, e)['K_max']
        for i, cell in enumerate(row['cells']):
            kk = 10.0 ** (lgk0 + (lgk1 - lgk0) * i / float(m['nk'] - 1))
            seen.add(cell)
            usable = cell not in ('no_steady_state', 'second_phase')
            assert usable == (kk <= kmax + 1e-12), (e, kk, cell)
    assert seen == {'no_steady_state', 'second_phase',
                    'reduction_limited', 'balanced'}


def test_an_oxidation_limited_steady_state_is_out_of_reach_here(p):
    """A consequence of the face ceiling worth stating rather than hiding.

    Oxidation-limited means the solid sits reduced, theta above 0.75. The
    rutile face is pinned at 0.5 by its own phase, so that state is on the
    far side of the boundary where no steady state exists at all: this
    carrier cannot be run oxidation-limited, whatever the rate constants.
    Checked against its own cause by lifting the ceiling and watching the
    region appear."""
    m = R.phase_map(p, nk=41, ne=31)
    assert 'oxidation_limited' not in {c for row in m['rows']
                                       for c in row['cells']}
    assert m['theta_surface_max'] < p['regime_thresholds']['theta_hi']
    lifted = R.phase_map(p, nk=41, ne=31, theta_surface_max=0.9)
    assert 'oxidation_limited' in {c for row in lifted['rows']
                                   for c in row['cells']}


def test_the_no_steady_state_region_does_not_depend_on_enrichment(p):
    """It is set by the face's own ceiling, so it is a vertical boundary.

    Worth gating because it is the one region the uniform model cannot
    produce at all, and a diagram that let it drift with E would be
    telling a different story than the algebra."""
    m = R.phase_map(p, nk=41, ne=31)
    cols = [[row['cells'][i] == 'no_steady_state' for row in m['rows']]
            for i in range(m['nk'])]
    for col in cols:
        assert len(set(col)) == 1, 'the boundary bends with enrichment'
