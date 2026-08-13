"""Checks on the thermodynamic data and on the equilibrium solver.

Run with:  python3 -m pytest tests/ -q
"""

import numpy as np
import pytest

from solidgas import (SHOMATE, R_KJ, mu0, enthalpy, entropy, heat_capacity, Kp,
                      SPECIES, ELEM, SOLIDS, gas_idx, solve, residual,
                      ti3_percent, moles_of_gas, G_total)

RWGS = {'CO2': -1, 'H2': -1, 'CO': 1, 'H2O': 1}


# ---------------------------------------------------------------- Shomate data

@pytest.mark.parametrize('name', list(SHOMATE))
def test_enthalpy_at_298_returns_formation_enthalpy(name):
    """The (H - H298) part of the Shomate form must vanish at 298.15 K.

    The condensed-phase fits carry more residual than the gases: their
    (H - H298) term crosses zero at 298.34 K (TiO2) and 298.77 K (Ti4O7)
    rather than exactly 298.15, leaving 0.01 and 0.13 kJ/mol. That is fit
    slop, not a data error, and it is four orders of magnitude below the
    formation enthalpies themselves.
    """
    if name == 'H2O':
        pytest.skip('H2O entry uses the 500-1700 K fit, not valid at 298 K')
    dHf = SHOMATE[name][2]
    tol = 0.2 if name in SOLIDS else 0.01
    assert enthalpy(name, 298.15) == pytest.approx(dHf, abs=tol)


@pytest.mark.parametrize('name,S298', [
    ('CO2', 213.79), ('CO', 197.66), ('H2', 130.68), ('CH4', 186.25),
])
def test_gas_entropy_at_298_matches_janaf(name, S298):
    assert entropy(name, 298.15) == pytest.approx(S298, abs=0.2)


@pytest.mark.parametrize('name,Cp298', [
    ('CO2', 37.13), ('CO', 29.14), ('H2', 28.84), ('CH4', 35.69), ('TiO2', 55.06),
])
def test_heat_capacity_at_298_matches_janaf(name, Cp298):
    """Cp is the cleanest check on the solid fits: rutile lands within 0.2."""
    assert heat_capacity(name, 298.15) == pytest.approx(Cp298, abs=0.2)


@pytest.mark.parametrize('name,natoms', [('TiO2', 3), ('Ti4O7', 11)])
def test_solid_heat_capacity_approaches_dulong_petit(name, natoms):
    """Above the Debye temperature a solid's Cp sits near 3R per atom."""
    limit = 3 * 8.314 * natoms
    for T in (773.15, 1073.15, 1773.15):
        Cp = heat_capacity(name, T)
        assert limit * 0.93 < Cp < limit * 1.20, f'{name} Cp({T}) = {Cp:.1f}'


def test_properties_are_monotonic_across_the_working_range():
    """Entropy rises with temperature for every species in 773-1773 K."""
    for name in SHOMATE:
        S = [entropy(name, T) for T in np.linspace(773.15, 1773.15, 25)]
        assert all(b > a for a, b in zip(S, S[1:])), name
        assert all(heat_capacity(name, T) > 0 for T in np.linspace(773.15, 1773.15, 25))


def test_mu0_is_h_minus_ts():
    for name in SHOMATE:
        for T in (600.0, 1200.0, 1700.0):
            assert mu0(name, T) == pytest.approx(
                enthalpy(name, T) - T * entropy(name, T) / 1000.0, rel=1e-12)


def test_shomate_branch_switch_is_continuous():
    """Both coefficient sets must agree at the 1000 K switch point."""
    for name in SHOMATE:
        lo, hi, dHf = SHOMATE[name]
        vals = []
        for c in (lo, hi):
            A, B, C, D, E, F, G, H = c
            t = 1.0
            Ht = dHf + (A * t + B * t**2 / 2 + C * t**3 / 3 + D * t**4 / 4 - E / t + F - H)
            St = A * np.log(t) + B * t + C * t**2 / 2 + D * t**3 / 3 - E / (2 * t**2) + G
            vals.append(Ht - 1000.0 * St / 1000.0)
        assert abs(vals[1] - vals[0]) < 0.05, f'{name} jumps {vals[1] - vals[0]:.4f} kJ/mol'


def test_rwgs_is_endothermic_and_crosses_unity():
    """RWGS is uphill at low T and turns over near 1100 K."""
    assert Kp(RWGS, 800.0) < 1.0
    assert Kp(RWGS, 1400.0) > 1.0
    lo, hi = 900.0, 1300.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if Kp(RWGS, mid) < 1.0:
            lo = mid
        else:
            hi = mid
    assert 1050.0 < mid < 1150.0


def test_rwgs_enthalpy_near_literature():
    """dH of RWGS at 298 K is about +41 kJ/mol, from formation enthalpies."""
    dH = sum(nu * SHOMATE[n][2] for n, nu in RWGS.items())
    assert dH == pytest.approx(41.2, abs=0.5)


def test_magneli_reduction_is_uphill_at_low_temperature():
    """4 TiO2 + H2 -> Ti4O7 + H2O must be strongly unfavourable at 773 K.

    This is what keeps rutile intact under the RWGS feed.
    """
    for T, expect_positive in [(773.15, True), (1773.15, True)]:
        dG = (mu0('Ti4O7', T) + mu0('H2O', T)
              - 4 * mu0('TiO2', T) - mu0('H2', T))
        assert bool(dG > 0) == expect_positive, f'dG({T}) = {dG:.1f}'
    dG_lo = mu0('Ti4O7', 773.15) + mu0('H2O', 773.15) - 4 * mu0('TiO2', 773.15) - mu0('H2', 773.15)
    dG_hi = mu0('Ti4O7', 1773.15) + mu0('H2O', 1773.15) - 4 * mu0('TiO2', 1773.15) - mu0('H2', 1773.15)
    assert dG_hi < dG_lo, 'reduction must get easier as temperature rises' 


# ------------------------------------------------------------- element balance

def test_element_matrix_matches_formulas():
    expected = {
        'CO2': (1, 0, 2, 0), 'H2': (0, 2, 0, 0), 'CO': (1, 0, 1, 0),
        'H2O': (0, 2, 1, 0), 'CH4': (1, 4, 0, 0),
        'TiO2': (0, 0, 2, 1), 'Ti4O7': (0, 0, 7, 4),
    }
    for i, s in enumerate(SPECIES):
        assert tuple(ELEM[i].astype(int)) == expected[s]


def test_ti3_percent_endpoints():
    n = np.zeros(len(SPECIES))
    n_Ti = 1.0
    n[SPECIES.index('Ti4O7')] = 0.0
    assert ti3_percent(n, n_Ti) == 0.0
    # all Ti as Ti4O7 -> 0.25 mol Ti4O7 per mol Ti -> half the Ti is Ti(3+)
    n[SPECIES.index('Ti4O7')] = 0.25
    assert ti3_percent(n, n_Ti) == pytest.approx(50.0)


# ------------------------------------------------------------------- the solver

@pytest.fixture(scope='module')
def rwgs_solution():
    T = 1073.15
    n0_TiO2 = 0.100 / 79.87
    n0 = moles_of_gas(1.0, 25.0, T) / 2
    b = np.array([n0, 2 * n0, 2 * n0 + 2 * n0_TiO2, n0_TiO2])
    inits = [np.array([n0 * a, n0 * b_, n0 * c, n0 * d, n0 * 1e-6,
                       n0_TiO2 * 0.99, n0_TiO2 * 0.0025])
             for a, b_, c, d in [(0.5, 0.5, 0.01, 0.01), (0.1, 0.1, 0.4, 0.4)]]
    n = solve(b, T, inits)
    assert n is not None
    return T, n, b, n0_TiO2


def test_solution_closes_element_balance(rwgs_solution):
    _T, n, b, _ = rwgs_solution
    assert residual(n, b) < 1e-10


def test_solution_reproduces_rwgs_equilibrium_constant(rwgs_solution):
    """The strongest check: the converged gas must satisfy Q == Kp."""
    T, n, _b, _ = rwgs_solution
    ng = sum(n[i] for i in gas_idx)
    y = {s: n[SPECIES.index(s)] / ng for s in ('CO2', 'H2', 'CO', 'H2O')}
    Q = (y['CO'] * y['H2O']) / (y['CO2'] * y['H2'])
    assert Q == pytest.approx(Kp(RWGS, T), rel=1e-3)


def test_rwgs_feed_does_not_reduce_titania(rwgs_solution):
    """A 1:1 CO2/H2 feed leaves rutile untouched."""
    _T, n, _b, n_Ti = rwgs_solution
    assert ti3_percent(n, n_Ti) < 1e-6


def test_solution_is_a_minimum(rwgs_solution):
    """Perturbing along a balance-preserving direction must not lower G."""
    T, n, b, _ = rwgs_solution
    G0 = G_total(n, T)
    # CO2 + H2 -> CO + H2O conserves every element
    d = np.zeros(len(SPECIES))
    for s, nu in RWGS.items():
        d[SPECIES.index(s)] = nu
    for eps in (1e-7, -1e-7):
        trial = n + eps * d
        if (trial > 0).all():
            assert residual(trial, b) < 1e-10
            assert G_total(trial, T) >= G0 - 1e-12


def test_methane_opens_the_magneli_branch():
    """Pure CH4 at 900 C reduces a large fraction of the Ti, unlike RWGS."""
    T = 1173.15
    n0_TiO2 = 0.100 / 79.87
    n0g = moles_of_gas(1.0, 25.0, T)
    b = np.array([n0g, 4 * n0g, 2 * n0_TiO2, n0_TiO2])
    inits = []
    for fr in [0.001, 0.05, 0.2, 0.5]:
        nT4 = n0_TiO2 * fr / 4 * 0.999
        inits.append(np.array([n0g * fr * 0.3, n0g * fr, n0g * fr, n0g * fr * 0.1,
                               n0g * (1 - fr), n0_TiO2 - 4 * nT4, nT4]))
    n = solve(b, T, inits)
    assert n is not None
    assert residual(n, b) < 1e-10
    assert ti3_percent(n, n0_TiO2) > 20.0
