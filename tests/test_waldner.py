"""Checks on the Waldner & Eriksson dataset and the way it is integrated."""

import numpy as np
import pytest

from solidgas import waldner as W
from solidgas import shomate as N
from solidgas import activeset as A

MAGNELI = ['Ti4O7', 'Ti5O9', 'Ti6O11', 'Ti7O13', 'Ti8O15', 'Ti9O17', 'Ti10O19']
REGISTRY = [m for m in MAGNELI if m in A.STOICH]
WORKING = (773.15, 1273.15, 1773.15)


# ------------------------------------------------------- the table as printed

def test_magneli_k0_steps_by_one_rutile_unit():
    """Each member adds a TiO2, so k0 must step by rutile's k0."""
    rutile_k0 = W.WALDNER['TiO2']['ranges'][0][3][0]
    k0 = [W.WALDNER[m]['ranges'][0][3][0] for m in MAGNELI]
    for a, b in zip(k0, k0[1:]):
        assert b - a == pytest.approx(rutile_k0, abs=1e-4)


def test_magneli_shear_term_is_independent_of_n():
    """k1 is the crystallographic shear contribution and does not scale."""
    k1 = {W.WALDNER[m]['ranges'][0][3][1] for m in MAGNELI}
    assert len(k1) == 1


def test_ti7o13_printed_enthalpy_is_a_misprint():
    """It is Ti4O7's value to four figures and breaks a linear series."""
    fit, slope, intercept = W._fit_ti7o13()
    assert abs(W.TI7O13_PRINTED - W.WALDNER['Ti4O7']['dHf298']) < 1500
    for m in MAGNELI:
        if W.WALDNER[m]['dHf298'] is None:
            continue
        n = W._stoich(m)[0]
        assert abs(W.WALDNER[m]['dHf298'] - (slope * n + intercept)) < 600, m
    assert abs(W.TI7O13_PRINTED - fit) > 2.5e6
    assert W.dHf('Ti7O13') == pytest.approx(fit)


def test_formation_enthalpy_is_monotonic_in_n():
    hs = [W.dHf(m) for m in MAGNELI]
    assert all(b < a for a, b in zip(hs, hs[1:]))


# ----------------------------------------------------- agreement with NIST

@pytest.mark.parametrize('name', ['TiO2', 'Ti4O7', 'Ti3O5', 'Ti2O3'])
def test_heat_capacity_agrees_with_nist_in_the_working_range(name):
    """Different assessments, but Cp is measured and both must land on it."""
    for T in WORKING:
        w, n = W.heat_capacity(name, T), N.heat_capacity(name, T)
        assert abs(w - n) / n < 0.01, f'{name} at {T}: {w:.2f} vs {n:.2f}'


@pytest.mark.parametrize('name,offset', [
    ('TiO2', -6.03), ('Ti4O7', -12.00), ('Ti2O3', -5.29),
])
def test_known_offsets_against_nist_are_stable(name, offset):
    """The sources differ, and by how much is worth pinning: it is the same
    size as the margins, which is why they are never mixed."""
    got = W.dHf(name) / 1000.0 - N.SHOMATE[name][2]
    assert got == pytest.approx(offset, abs=0.05)


# ------------------------------------------------------- piecewise integration

@pytest.mark.parametrize('name,T_join', [('Ti2O3', 470.0), ('Ti3O5', 450.0)])
def test_gibbs_energy_is_continuous_across_cp_ranges(name, T_join):
    """G must not jump, even though H and S each do."""
    eps = 1e-6
    lo, hi = W.mu0(name, T_join - eps), W.mu0(name, T_join + eps)
    assert lo == pytest.approx(hi, abs=1e-5)


@pytest.mark.parametrize('name,T_join,dH', [('Ti2O3', 470.0, 1138.0),
                                            ('Ti3O5', 450.0, 11757.0)])
def test_enthalpy_and_entropy_step_by_the_transition(name, T_join, dH):
    """The second enthalpy on the row is a latent heat, not a formation term.

    Enforcing continuity in H instead leaves G continuous too, so the mistake
    hides - but it displaces mu by dH*(1 - T/T_trans), which is -35 kJ/mol for
    Ti3O5 by 1500 C and moves the phase out of position entirely.
    """
    eps = 1e-6
    assert W.enthalpy(name, T_join + eps) - W.enthalpy(name, T_join - eps) \
        == pytest.approx(dH, abs=0.1)
    assert W.entropy(name, T_join + eps) - W.entropy(name, T_join - eps) \
        == pytest.approx(dH / T_join, abs=1e-3)


def test_ti3o5_lands_near_nist_once_the_transition_is_carried():
    """Two assessments should differ by a few kJ, not by tens."""
    for T in WORKING:
        assert abs(W.mu0('Ti3O5', T) - N.mu0('Ti3O5', T)) < 10.0


@pytest.mark.parametrize('name', ['TiO2'] + MAGNELI)
def test_entropy_rises_and_heat_capacity_stays_positive(name):
    S = [W.entropy(name, T) for T in np.linspace(773.15, 1773.15, 20)]
    assert all(b > a for a, b in zip(S, S[1:]))
    assert all(W.heat_capacity(name, T) > 0 for T in np.linspace(773.15, 1773.15, 20))


# ------------------------------------------------------------- what it changes

def test_shallow_magneli_phases_form_more_easily_than_ti4o7():
    """The whole reason for bringing this dataset in."""
    for T in WORKING:
        crits = [A.critical_potential(m, 'TiO2', T) for m in REGISTRY]
        assert all(b > a for a, b in zip(crits, crits[1:])), crits


def test_rwgs_survives_the_full_series():
    """The headline result, measured against the easiest phase in the set."""
    for T_C in np.linspace(500, 1500, 21):
        r = A.solve({'CO2': 1, 'H2': 1}, T_C + 273.15)
        assert r['active_condensed_phases'] == ['TiO2'], (T_C, r['active_condensed_phases'])
        costs = {k: v['per_formula_kJ']
                 for k, v in r['inactive_phase_reduced_costs'].items()}
        assert min(costs.values()) > 0, (T_C, costs)


def test_the_easiest_phase_is_the_highest_n_member():
    """Among the phases carried, the shallowest reduction is the cheapest."""
    for T in WORKING:
        r = A.solve({'CO2': 1, 'H2': 1}, T)
        costs = {k: v['per_formula_kJ']
                 for k, v in r['inactive_phase_reduced_costs'].items()}
        assert min(costs, key=costs.get) == 'Ti10O19', (T, costs)


def test_sources_disagree_enough_to_matter():
    """If these ever came out equal, mixing them would be harmless and the
    single-source rule could be dropped. They do not."""
    T = 1773.15
    w = W.mu0('Ti4O7', T) - 4 * W.mu0('TiO2', T)
    n = N.mu0('Ti4O7', T) - 4 * N.mu0('TiO2', T)
    assert abs(w - n) > 10.0


def test_the_phase_registry_is_ordered_by_oxygen_content():
    r = [A.STOICH[p][1] / A.STOICH[p][0] for p in A.SOLIDS]
    assert r == sorted(r, reverse=True)
