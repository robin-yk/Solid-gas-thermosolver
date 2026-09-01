"""The reduction line: closed form, oracle, browser mirror, and figures.

The question "how reducing does the gas have to be before the solid gives
way" already had an answer - bisect `activeset.solve` over feed
composition, sixty solves a temperature.  `activeset.reduction_boundary`
answers it in two exponentials instead, and the first test here is the
one that matters: the shortcut has to land where the long way lands.

The shortcut leaves nothing out.  Free O2 was tempting to drop - with
hydrogen in the feed it is worth 1.6e-5 kJ per mol O even at 1650 C - and
the second test is here because dropping it would have been wrong: take
the hydrogen away and it is worth 5.8 kJ, so a CO2/CO feed would have been
placed on the wrong side of its own line.

These are release gates: a missing interpreter or artefact fails the test
rather than skipping it.
"""

import json
import math
import pathlib
import struct

import pytest

from solidgas import activeset as A

ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA = ROOT / 'data'
REF = DATA / 'reference_results_high_precision.json'

PANEL_LO, PANEL_HI = 300.0, 1650.0          # the range the workspace admits

FEEDS = ({'CO2': 0.5, 'H2': 0.5}, {'CO2': 0.9, 'H2': 0.1},
         {'CO2': 0.1, 'H2': 0.9}, {'CO2': 0.99, 'CO': 0.01},
         {'H2': 0.98, 'H2O': 0.02}, {'CO2': 0.3, 'H2': 0.3, 'N2': 0.4},
         {'CO2': 1.0}, {'CO2': 0.9, 'CO': 0.1})


def _num(v):
    return float(v) if isinstance(v, str) else v


def _bisect_with_the_solver(T_C, host='TiO2', iters=70):
    """The long way round: solve the whole equilibrium, over and over."""
    T_K = T_C + 273.15

    def intact(y):
        r = A.solve({'CO2': y, 'H2': 1.0 - y}, T_K, initial_solid=host)
        return r['active_condensed_phases'] == [host]

    lo, hi = 1e-12, 0.5
    assert not intact(lo) and intact(hi), (T_C, 'bracket does not straddle')
    for _ in range(iters):
        mid = math.sqrt(lo * hi)            # geometric: the answer spans
        if intact(mid):                     # seven decades over this range
            hi = mid
        else:
            lo = mid
    return math.sqrt(lo * hi)


# ------------------------------------------------------ the shortcut holds


@pytest.mark.parametrize('T_C', [400.0, 600.0, 900.0, 1200.0, 1500.0])
def test_the_closed_form_lands_where_the_solver_bisects_to(T_C):
    want = _bisect_with_the_solver(T_C)
    got = A.reduction_boundary(T_C + 273.15)['y_CO2']
    rel = abs(got - want) / want
    assert rel < 2e-4, (
        '%.0f C: closed form %.6g, solver bisection %.6g, %.2e apart'
        % (T_C, got, want, rel))


def test_the_closed_form_is_the_solvers_own_answer_where_nothing_reduces():
    """Where the host still stands, the multiplier the solver puts on
    oxygen is the gas potential itself, so the two are directly comparable
    and the closed form has nowhere to hide.

    One family has to be held to a looser number, and it is the solver
    that is loose there, not the closed form. A gas with no reduced member
    in it - pure CO2 at 400 C - reaches equilibrium with CO at 3e-12 of
    the mixture, and a Newton solve on a scaled residual of 1e-13 is
    working at its own floor. The closed form has no trace species to
    resolve: it solves for the potential and reads the composition off it.
    """
    tight = loose = 0.0
    n_tight = n_loose = 0
    for feed in FEEDS:
        for T_C in (400.0, 900.0, 1500.0, 1650.0):
            r = A.solve(feed, T_C + 273.15)
            if r['active_condensed_phases'] != ['TiO2']:
                continue
            mu = A.gas_oxygen_potential(feed, T_C + 273.15)
            gap = abs(mu - r['lambda_kJ_per_mol']['O'])
            reduced = min(r['gas_fractions'][g] for g in ('CO', 'H2')
                          if r['gas_fractions'][g] > 0.0)
            if reduced < 1e-10:
                loose, n_loose = max(loose, gap), n_loose + 1
            else:
                tight, n_tight = max(tight, gap), n_tight + 1
    assert n_tight and n_loose, 'both branches must be exercised'
    assert tight < 1e-10, 'closed form is %.2e kJ/mol-O off' % tight
    assert loose < 1e-3, 'even the trace-limited cases are %.2e off' % loose


def test_free_O2_is_worth_nothing_in_the_rwgs_feed_and_kilojoules_without_it():
    """Why the O2 term is carried at all, measured rather than argued.

    The comparison is against the same closed form with the O2 term
    switched off, which is what this code looked like before the term
    went in. In the feed this whole project is about - CO2 and H2, one to
    one - it is worth a hundredth of a joule even at 1650 C, and leaving
    it out really would have been harmless. Strip the hydrogen out and
    the reason it is here shows up: a CO2/CO mixture runs 0.2 mol% O2 at
    1650 C and the two answers part by kilojoules, which is enough to put
    the feed on the wrong side of its own line.
    """
    def without_o2(feed, T_K):
        saved = A._o2_per_unit
        A._o2_per_unit = lambda mu, T, P: 0.0
        try:
            return A.gas_oxygen_potential(feed, T_K)
        finally:
            A._o2_per_unit = saved

    def worth(feed, T_C):
        T_K = T_C + 273.15
        return abs(A.gas_oxygen_potential(feed, T_K) - without_o2(feed, T_K))

    assert worth({'CO2': 0.5, 'H2': 0.5}, 1650.0) < 1e-4
    assert worth({'CO2': 0.99, 'CO': 0.01}, 1650.0) > 1.0
    assert worth({'CO2': 1.0}, 400.0) > 10.0
    # The feed's oxygen count is fixed, so giving it another place to put
    # oxygen can only lower the potential it settles at - never raise it.
    for feed in FEEDS:
        for T_C in (400.0, 1650.0):
            T_K = T_C + 273.15
            a, b = A.gas_oxygen_potential(feed, T_K), without_o2(feed, T_K)
            if a is None or b is None:
                continue
            assert a <= b + 1e-12, (feed, T_C, a, b)


def test_the_total_pressure_reaches_exactly_the_feeds_it_should():
    """A partial pressure is not a property of a mixture on its own.

    Wherever O2 carries the potential the answer moves with the total
    pressure, and wherever a couple carries it the answer does not. Both
    halves matter: the first is why P is an argument at all, the second is
    why the boundary curve can be drawn without asking about it.
    """
    T_K = 1923.15
    dry = [A.gas_oxygen_potential({'CO2': 1.0}, T_K, P)
           for P in (0.1, 1.0, 10.0)]
    assert dry[2] - dry[0] > 20.0, dry
    wet = [A.gas_oxygen_potential({'CO2': 0.5, 'H2': 0.5}, T_K, P)
           for P in (0.1, 1.0, 10.0)]
    assert abs(wet[2] - wet[0]) < 1e-3, wet
    line = [A.reduction_boundary(T_K, P_atm=P)['y_CO2']
            for P in (0.1, 1.0, 10.0)]
    assert abs(line[2] - line[0]) / line[0] < 1e-6, line


# ---------------------------------------------------------------- the axis


@pytest.mark.parametrize('y', [0.5, 0.1, 0.01, 1e-4, 1e-6, 0.9, 0.99])
def test_the_axis_inverts_exactly_for_a_two_component_feed(y):
    """An H2/CO2 mixture must land back on the fraction it was mixed at."""
    for T_C in (400.0, 900.0, 1500.0):
        T_K = T_C + 273.15
        mu = A.gas_oxygen_potential({'CO2': y, 'H2': 1.0 - y}, T_K)
        back = A.equivalent_co2_fraction(mu, T_K)
        assert abs(back - y) / y < 1e-12, (T_C, y, back)


def test_a_feed_that_sets_no_potential_says_so_instead_of_guessing():
    """Two ways to have no answer, and each is named rather than swallowed.

    A gas with no carbon and no hydrogen has no couple to set a potential
    with - the oxygen in it is the whole story, and where it sits is a
    property of the vessel. A gas short of the oxygen to make every carbon
    a CO is not a C-H-O mixture in this registry at all, it is that
    mixture plus solid carbon.
    """
    for feed in ({'N2': 1.0}, {'O2': 1.0}, {'N2': 1.0, 'O2': 1e-6}):
        pt = A.feed_on_boundary_axis(feed, 873.15)
        assert pt['y_CO2'] is None and 'couple' in pt['why'], feed
    for feed in ({'CO': 1.0}, {'H2': 1.0}, {'CO': 1.0, 'H2': 1.0}):
        pt = A.feed_on_boundary_axis(feed, 873.15)
        assert pt['y_CO2'] is None and 'carbon' in pt['why'], feed
    for feed in ({'CO2': 1.0}, {'H2O': 1.0}, {'CO2': 1.0, 'H2O': 1.0},
                 {'CO': 1.0, 'CO2': 1e-9}, {'H2': 1.0, 'H2O': 1e-9}):
        assert A.feed_on_boundary_axis(feed, 873.15)['y_CO2'] is not None, feed


def test_the_axis_marks_which_feeds_it_was_built_for():
    """H2 and CO2 land at the fraction they were mixed at; the rest land
    at the mixture that would be equally reducing, and are told apart so
    the page can say which it is showing."""
    for feed in ({'CO2': 0.5, 'H2': 0.5}, {'CO2': 1.0},
                 {'CO2': 0.3, 'H2': 0.3, 'N2': 0.4}):
        assert A.feed_on_boundary_axis(feed, 873.15)['literal'] is True, feed
    for feed in ({'CO2': 0.4, 'CO': 0.1, 'H2': 0.5}, {'H2': 0.98, 'H2O': 0.02}):
        assert A.feed_on_boundary_axis(feed, 873.15)['literal'] is False, feed


# --------------------------------------------------------------- the shape


def test_the_line_rises_with_temperature_over_the_whole_panel_range():
    prev = None
    T_C = PANEL_LO
    while T_C <= PANEL_HI + 1e-9:
        y = A.reduction_boundary(T_C + 273.15)['y_CO2']
        assert 0.0 < y < 1.0, (T_C, y)
        if prev is not None:
            assert y > prev, 'the line dips at %.0f C' % T_C
        prev = y
        T_C += 5.0


def test_the_first_phase_out_is_the_one_the_kkt_solve_also_names():
    """Two routes to the same phase.

    Per mole of oxygen, an excluded phase's reduced cost is the gap
    between the gas potential and that phase's critical potential, so
    ranking phases by reduced cost and ranking them by critical potential
    have to give the same winner. Per formula unit they need not, because
    the deficiency each is divided by differs, and that column is not the
    one to compare against.
    """
    for T_C in (400.0, 700.0, 1000.0, 1300.0, 1650.0):
        r = A.solve({'CO2': 0.5, 'H2': 0.5}, T_C + 273.15)
        rc = r['inactive_phase_reduced_costs']
        nearest = min(rc, key=lambda s: rc[s]['per_mol_O_kJ'])
        assert A.reduction_boundary(T_C + 273.15)['phase'] == nearest, T_C


def test_a_more_reduced_host_needs_a_more_reducing_gas():
    """Each rung of the ladder is harder to get off than the one above."""
    prev = None
    for host in ('TiO2', 'Ti4O7', 'Ti3O5'):
        b = A.reduction_boundary(873.15, host=host)
        assert b is not None and b['host'] == host
        if prev is not None:
            assert b['y_CO2'] < prev, host
        prev = b['y_CO2']
    assert A.reduction_boundary(873.15, host='Ti2O3') is None, \
        'nothing in the registry is more reduced than Ti2O3'


def test_the_line_and_the_solver_agree_on_which_side_a_feed_is_on():
    """The figure's whole claim, checked against the full equilibrium.

    The marker sits above the line exactly when the solver leaves the host
    standing - including a pair of feeds set either side of the line at
    900 C, five percent apart, where nothing but the real answer separates
    them.
    """
    cases = [({'CO2': 0.5, 'H2': 0.5}, 600.0),
             ({'CO2': 1e-4, 'H2': 1.0}, 900.0),
             ({'CO2': 0.02, 'H2': 0.98}, 1400.0),
             ({'CO2': 0.0040, 'H2': 0.9960}, 900.0),
             ({'CO2': 0.0042, 'H2': 0.9958}, 900.0),
             ({'CO': 1.0, 'CO2': 1.0}, 700.0)]
    for feed, T_C in cases:
        T_K = T_C + 273.15
        mu = A.gas_oxygen_potential(feed, T_K)
        y = A.equivalent_co2_fraction(mu, T_K)
        line = A.reduction_boundary(T_K)['y_CO2']
        intact = A.solve(feed, T_K)['active_condensed_phases'] == ['TiO2']
        assert (y > line) == intact, (feed, T_C, y, line, intact)


# --------------------------------------------------------------- the oracle


def test_the_oracle_agrees_on_the_line():
    doc = json.loads(REF.read_text())
    blk = doc['co2_h2_reduction_line']
    assert blk['host'] == 'TiO2'
    for key, row in blk['boundary_by_T'].items():
        T_K = float(key) + 273.15
        got = A.reduction_boundary(T_K, P_atm=blk['P_atm'])
        assert got['phase'] == row['phase'], key
        assert abs(got['mu_crit_kJ_per_mol_O']
                   - _num(row['mu_crit_kJ_per_mol_O'])) < 1e-11, key
        want = _num(row['y_CO2_boundary'])
        assert abs(got['y_CO2'] - want) / want < 1e-12, key


def test_the_oracle_agrees_on_what_a_feed_is_worth():
    doc = json.loads(REF.read_text())
    blk = doc['co2_h2_reduction_line']
    P = blk['P_atm']
    for row in blk['feeds']:
        T_K = row['T_C'] + 273.15
        mu = A.gas_oxygen_potential(row['feed'], T_K, P)
        assert mu is not None, row['feed']
        assert abs(mu - _num(row['mu_O_kJ_per_mol'])) < 1e-10, row['feed']
        want = _num(row['y_CO2_equivalent'])
        got = A.equivalent_co2_fraction(mu, T_K, P)
        assert abs(got - want) / want < 1e-12, row['feed']
