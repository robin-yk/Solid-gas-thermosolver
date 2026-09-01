"""Checks on the oxygen-potential route.

These are written as invariants rather than as fixed numbers: adding a phase
to the series must not break them, and they should fail only when something is
actually wrong.
"""

import numpy as np
import pytest

from solidgas import (ACTIVE_TI_PHASES, FORMULAS, Kp, R_KJ, mu0, moles_of_gas,
                      mu_O_critical, mu_O_from_gas, phase_ladder,
                      stable_phase, gas_equilibrium, survey, ti_phases)
from solidgas import activeset as A
from solidgas.potential import _ratio

RWGS = {'CO2': -1, 'H2': -1, 'CO': 1, 'H2O': 1}


def rwgs_gas_basis(T):
    n0 = moles_of_gas(1.0, 25.0, T) / 2
    return np.array([n0, 2 * n0, 2 * n0]), n0


# ------------------------------------------------------------------- the gas

@pytest.mark.parametrize('T_C', [500, 800, 1100, 1500])
def test_gas_only_equilibrium_satisfies_kp(T_C):
    T = T_C + 273.15
    b, _ = rwgs_gas_basis(T)
    n = gas_equilibrium(b, T)
    tot = sum(n.values())
    y = {k: v / tot for k, v in n.items()}
    Q = (y['CO'] * y['H2O']) / (y['CO2'] * y['H2'])
    assert Q == pytest.approx(Kp(RWGS, T), rel=1e-4)


@pytest.mark.parametrize('T_C', [500, 800, 1100, 1500])
def test_both_redox_couples_report_the_same_oxygen_potential(T_C):
    """H2/H2O and CO/CO2 must agree once the gas is equilibrated."""
    T = T_C + 273.15
    r = survey(*rwgs_gas_basis(T)[:1], T)
    assert r['couple_spread'] < 0.05, r['couple_spread']


def test_oxygen_potential_falls_as_temperature_rises():
    """A fixed RWGS charge gets more reducing with temperature."""
    mus = []
    for T_C in np.linspace(500, 1500, 20):
        T = T_C + 273.15
        mus.append(survey(rwgs_gas_basis(T)[0], T)['mu_O'])
    assert all(b < a for a, b in zip(mus, mus[1:])), mus


# ------------------------------------------------------------------ the ladder

@pytest.mark.parametrize('source', ['waldner', 'nist'])
def test_ladder_is_ordered_by_oxygen_content(source):
    for T_C in (500, 1000, 1500):
        ladder = phase_ladder(T_C + 273.15, source=source)
        ratios = [_ratio(h) for h, _l, _c in ladder]
        assert ratios == sorted(ratios, reverse=True)


@pytest.mark.parametrize('source', ['waldner', 'nist'])
def test_each_step_needs_a_lower_potential_than_the_one_before(source):
    """Reducing further must always cost more, at every temperature.

    Ti7O13 is left out: its formation enthalpy is misprinted in the paper and
    the value reconstructed from the rest of the series carries about 400 J of
    uncertainty, against step gaps of one to two kJ. That is not enough to
    place it, and a reversal there says nothing about the assessment. See
    test_ti7o13_cannot_be_placed_from_the_reconstruction.
    """
    phases = [p for p in ti_phases(source) if p != 'Ti7O13']
    for T_C in (700, 1000, 1500):
        crits = [c for _h, _l, c in phase_ladder(T_C + 273.15, phases=phases,
                                                 source=source)]
        assert all(b < a for a, b in zip(crits, crits[1:])), (T_C, crits)


def test_ti7o13_cannot_be_placed_from_the_reconstruction():
    """Documents the one thing the repaired value is not good enough for.

    It does not touch the headline: the phase that forms first is the highest-n
    member, and every margin quoted is measured against that.
    """
    full = [c for _h, _l, c in phase_ladder(973.15, source='waldner')]
    assert not all(b < a for a, b in zip(full, full[1:]))
    without = [c for _h, _l, c in phase_ladder(
        973.15, phases=[p for p in ti_phases('waldner') if p != 'Ti7O13'],
        source='waldner')]
    assert all(b < a for a, b in zip(without, without[1:]))


def test_critical_potential_is_symmetric_in_its_definition():
    """mu_O_critical must not depend on how the reaction was scaled."""
    T = 1273.15
    direct = mu_O_critical('Ti4O7', 'TiO2', T, source='nist')
    # 4 TiO2 -> Ti4O7 + O, written out by hand
    by_hand = 4 * mu0('TiO2', T) - mu0('Ti4O7', T)
    assert direct == pytest.approx(by_hand, rel=1e-12)


def test_rejects_a_step_that_is_not_a_reduction():
    with pytest.raises(ValueError):
        mu_O_critical('TiO2', 'Ti4O7', 1273.15)


# --------------------------------------------------- agreement with the solver

@pytest.mark.parametrize('T_C', [600, 900, 1200, 1500])
def test_potential_route_matches_the_full_minimisation(T_C):
    """The decomposition must reproduce what the production solver returns.

    Both routes are run on the same four gases.  The potential route's own
    registry also carries CH4, which the active-set registry does not, and
    below about 800 C that difference is worth several mole percent - so the
    comparison names its species instead of inheriting them.  Once they are
    solving the same problem the two agree to eight decimals, not four.
    """
    T = T_C + 273.15
    b_gas, _n0 = rwgs_gas_basis(T)
    four = ['CO2', 'H2', 'CO', 'H2O']

    n = gas_equilibrium(b_gas, T, species=four)
    tot = sum(n.values())
    y = {k: v / tot for k, v in n.items()}
    full = A.solve({'CO2': 1, 'H2': 1}, T)

    for s in four:
        assert y[s] == pytest.approx(full['gas_fractions'][s], abs=1e-7), s

    # and both agree that nothing reduces
    mu_O, _spread = mu_O_from_gas(y, T)
    assert stable_phase(mu_O, T, ACTIVE_TI_PHASES)[0] == 'TiO2'
    assert full['active_condensed_phases'] == ['TiO2']
    assert full['reduced_pct'] < 1e-3


def test_rwgs_leaves_rutile_across_the_whole_sweep():
    """The headline result, stated as an invariant on the margin."""
    for T_C in np.linspace(500, 1500, 25):
        r = survey(rwgs_gas_basis(T_C + 273.15)[0], T_C + 273.15)
        assert r['phase'] == 'TiO2'
        assert r['margin'] > 0


def test_margin_shrinks_with_temperature():
    """Where the conclusion is weakest is the top of the range, not the bottom."""
    lo = survey(rwgs_gas_basis(773.15)[0], 773.15)['margin']
    hi = survey(rwgs_gas_basis(1773.15)[0], 1773.15)['margin']
    assert hi < lo
