"""Checks on the thermodynamic data and on the equilibrium solver.

Run with:  python3 -m pytest tests/ -q
"""

import numpy as np
import pytest

from solidgas import (SHOMATE, VALID_RANGE, R_KJ, mu0, enthalpy, entropy,
                      heat_capacity, Kp, SPECIES, ELEM, ELEMENTS, SOLIDS,
                      FORMULAS, TI_PHASES, ACTIVE_TI_PHASES, gas_idx, solve,
                      residual, ti3_percent, phase_split, mean_valence,
                      phase_seeds, moles_of_gas, G_total)

RWGS = {'CO2': -1, 'H2': -1, 'CO': 1, 'H2O': 1}


# ---------------------------------------------------------------- Shomate data

# Parametrised over exactly the entries the property applies to, instead of
# skipping at run time: H2O and Ti2O3 have fits that start above 298.15 K, and
# the ceria pair is pinned as known-bad in test_ceria_fit_is_not_self_consistent.
_H298_NAMES = [n for n in SHOMATE
               if VALID_RANGE[n][0] <= 298.15 and n not in ('CeO2', 'Ce2O3')]


@pytest.mark.parametrize('name', _H298_NAMES)
def test_enthalpy_at_298_returns_formation_enthalpy(name):
    """The (H - H298) part of the Shomate form must vanish at 298.15 K.

    The condensed-phase fits carry more residual than the gases: their
    (H - H298) term crosses zero at 298.34 K (TiO2) and 298.77 K (Ti4O7)
    rather than exactly 298.15, leaving 0.01 and 0.13 kJ/mol. That is fit
    slop, not a data error, and it is four orders of magnitude below the
    formation enthalpies themselves.
    """
    dHf = SHOMATE[name][2]
    condensed = 'Ti' in FORMULAS.get(name, {}) or name.startswith('Ti')
    tol = 0.2 if condensed else 0.01
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
    """Entropy rises with temperature for every species, inside its own fit range."""
    for name in SHOMATE:
        lo, hi = VALID_RANGE[name]
        grid = np.linspace(max(773.15, lo), min(1773.15, hi), 25)
        S = [entropy(name, T) for T in grid]
        assert all(b > a for a, b in zip(S, S[1:])), name
        assert all(heat_capacity(name, T) > 0 for T in grid), name



def test_ceria_fit_is_not_self_consistent():
    """The ceria coefficients do not close on their own formation enthalpy.

    Written down rather than corrected. These come from the notebook this work
    started in, and replacing them with a better-sourced set would break parity
    with the calculation everything else here is anchored to - so they are used
    as they arrive and the defect is pinned here, where a change to it fails a
    test instead of passing silently.

    The size matters. Reducing CeO2 to Ce2O3 removes one oxygen per two cerium,
    so the gap propagates as 2*(-20.5) - (-29.9) = -11.1 kJ into the reduction
    enthalpy, which is about 11 kJ per mole of oxygen on the phase boundary
    itself. Real NIST gas entries close to 0.003 kJ/mol.
    """
    for name, gap in (('CeO2', -20.513), ('Ce2O3', -29.929)):
        assert enthalpy(name, 298.15) - SHOMATE[name][2] == pytest.approx(gap, abs=0.01), name
    for name in ('CO2', 'H2O', 'CO', 'CH4'):
        assert abs(enthalpy(name, 298.15) - SHOMATE[name][2]) < 0.01, name


def test_ceria_reduction_boundary_carries_that_error():
    """How far the ceria phase boundary moves because of the gap above."""
    T = 1173.15
    exact = 2 * mu0('CeO2', T) - mu0('Ce2O3', T)
    assert exact == pytest.approx(-402.669, abs=0.5)


# --------------------------------------------------------------- the TiOx set

def test_ti2o3_entry_is_the_solid_not_the_liquid():
    """NIST lists both; the liquid (dHf = -1418.46) must not be used here."""
    assert SHOMATE['Ti2O3'][2] == pytest.approx(-1520.884, abs=0.01)
    assert entropy('Ti2O3', 298.15) != pytest.approx(127.12, abs=1.0)


def test_ti3o5_beta_is_the_stable_polymorph_in_range():
    """alpha wins below about 450 K, beta above - the whole sweep is above."""
    assert mu0('Ti3O5_alpha', 350.0) < mu0('Ti3O5', 350.0)
    for T in (773.15, 1173.15, 1773.15):
        assert mu0('Ti3O5', T) < mu0('Ti3O5_alpha', T)


def test_magneli_phases_carry_exactly_two_ti3_each():
    """Charge balance on Ti_n O_(2n-1) gives x = 2 for every n."""
    for name, (n_ti, n_ti3) in TI_PHASES.items():
        f = FORMULAS[name]
        assert f['Ti'] == n_ti
        charge_on_ti = 2 * f['O']
        assert 3 * n_ti3 + 4 * (n_ti - n_ti3) == charge_on_ti


def test_reduction_ladder_gets_harder_at_high_temperature():
    """Each successive step must cost more than the last where it matters."""
    steps = [{'TiO2': -4, 'Ti4O7': 1}, {'Ti4O7': -3, 'Ti3O5': 4}, {'Ti3O5': -2, 'Ti2O3': 3}]
    for T in (1173.15, 1573.15, 1773.15):
        dG = [mu0('H2O', T) - mu0('H2', T) + sum(v * mu0(k, T) for k, v in rx.items())
              for rx in steps]
        assert dG[0] < dG[1] < dG[2], f'T={T}: {dG}'
        assert all(g > 0 for g in dG), 'every step is uphill against pure H2'


def test_ti3o5_is_stable_against_disproportionation_above_700C():
    """Ti4O7 + Ti2O3 <-> 2 Ti3O5; negative dG means Ti3O5 is a real phase."""
    for T_C in (700, 900, 1100, 1300, 1500):
        T = T_C + 273.15
        dG = 2 * mu0('Ti3O5', T) - mu0('Ti4O7', T) - mu0('Ti2O3', T)
        assert dG < 0, f'{T_C} C: dG = {dG:.2f}'
    # and it does come apart at the cold end
    T = 500 + 273.15
    assert 2 * mu0('Ti3O5', T) - mu0('Ti4O7', T) - mu0('Ti2O3', T) > 0


def test_mu0_is_h_minus_ts():
    for name in SHOMATE:
        for T in (600.0, 1200.0, 1700.0):
            assert mu0(name, T) == pytest.approx(
                enthalpy(name, T) - T * entropy(name, T) / 1000.0, rel=1e-12)


def test_shomate_branch_switch_is_continuous():
    """Both coefficient sets must agree where that species actually switches.

    NIST does not use one breakpoint for everything - O2 changes at 700 K and
    N2 at 500 K - so checking every entry at 1000 K tests the wrong join.
    """
    from solidgas.shomate import BREAKPOINT, DEFAULT_BREAKPOINT
    for name in SHOMATE:
        lo, hi, dHf = SHOMATE[name]
        T_sw = BREAKPOINT.get(name, DEFAULT_BREAKPOINT)
        vals = []
        for c in (lo, hi):
            A, B, C, D, E, F, G, H = c
            t = T_sw / 1000.0
            Ht = dHf + (A * t + B * t**2 / 2 + C * t**3 / 3 + D * t**4 / 4 - E / t + F - H)
            St = A * np.log(t) + B * t + C * t**2 / 2 + D * t**3 / 3 - E / (2 * t**2) + G
            vals.append(Ht - T_sw * St / 1000.0)
        assert abs(vals[1] - vals[0]) < 0.05, \
            f'{name} jumps {vals[1] - vals[0]:.4f} kJ/mol at {T_sw:.0f} K'


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
        'TiO2': (0, 0, 2, 1), 'Ti20O39': (0, 0, 39, 20),
        'Ti10O19': (0, 0, 19, 10), 'Ti9O17': (0, 0, 17, 9),
        'Ti8O15': (0, 0, 15, 8), 'Ti7O13': (0, 0, 13, 7),
        'Ti6O11': (0, 0, 11, 6), 'Ti5O9': (0, 0, 9, 5),
        'Ti4O7': (0, 0, 7, 4), 'Ti3O5': (0, 0, 5, 3), 'Ti2O3': (0, 0, 3, 2),
    }
    for i, s in enumerate(SPECIES):
        assert tuple(ELEM[i].astype(int)) == expected[s]


def test_every_active_phase_has_data_and_a_seed():
    """Every active phase needs stoichiometry, a seed, and at least one source.

    Not every one is in NIST - the shallow Magneli members exist only in the
    Waldner assessment, which is why they were brought in.
    """
    from solidgas import waldner
    for p in ACTIVE_TI_PHASES:
        assert p in FORMULAS and p in TI_PHASES
        assert p in SHOMATE or p in waldner.WALDNER, p
    assert len(phase_seeds(1.0)) == 1 + 2 * (len(ACTIVE_TI_PHASES) - 1)


def test_phase_seeds_conserve_titanium():
    for sv in phase_seeds(1.0):
        total = sum(v * TI_PHASES[p][0] for v, p in zip(sv, ACTIVE_TI_PHASES))
        assert total == pytest.approx(1.0, abs=1e-9)


def test_ti3_percent_endpoints():
    n = np.zeros(len(SPECIES))
    assert ti3_percent(n, 1.0) == 0.0
    # all Ti as Ti4O7 -> 0.25 mol Ti4O7 per mol Ti -> half the Ti is Ti(3+)
    n[SPECIES.index('Ti4O7')] = 0.25
    assert ti3_percent(n, 1.0) == pytest.approx(50.0)
    assert mean_valence(n, 1.0) == pytest.approx(3.5)


@pytest.mark.parametrize('phase,expect_ti3,expect_valence', [
    ('TiO2', 0.0, 4.0), ('Ti4O7', 50.0, 3.5),
    ('Ti3O5', 200 / 3, 4 - 2 / 3), ('Ti2O3', 100.0, 3.0),
])
def test_ti3_percent_for_each_pure_phase(phase, expect_ti3, expect_valence):
    """One mole of Ti entirely as one phase."""
    n = np.zeros(len(SPECIES))
    n[SPECIES.index(phase)] = 1.0 / TI_PHASES[phase][0]
    assert ti3_percent(n, 1.0) == pytest.approx(expect_ti3)
    assert mean_valence(n, 1.0) == pytest.approx(expect_valence)


def test_phase_split_sums_to_all_titanium():
    n = np.zeros(len(SPECIES))
    n[SPECIES.index('TiO2')] = 0.4
    n[SPECIES.index('Ti4O7')] = 0.1     # 0.4 mol Ti
    n[SPECIES.index('Ti3O5')] = 0.2 / 3  # 0.2 mol Ti
    assert sum(phase_split(n, 1.0).values()) == pytest.approx(100.0)


# ------------------------------------------------------------------- the solver

@pytest.fixture(scope='module')
def rwgs_solution():
    T = 1073.15
    n0_TiO2 = 0.100 / 79.87
    n0 = moles_of_gas(1.0, 25.0, T) / 2
    b = np.array([n0, 2 * n0, 2 * n0 + 2 * n0_TiO2, n0_TiO2])
    inits = [np.array([n0 * a, n0 * b_, n0 * c, n0 * d, n0 * 1e-6] + list(sv))
             for a, b_, c, d in [(0.5, 0.5, 0.01, 0.01), (0.1, 0.1, 0.4, 0.4)]
             for sv in phase_seeds(n0_TiO2)]
    n = solve(b, T, inits)
    assert n is not None
    return T, n, b, n0_TiO2


def test_solution_closes_element_balance(rwgs_solution):
    """Absolute closure, against a titanium basis of 1.25e-3 mol."""
    _T, n, b, _ = rwgs_solution
    assert residual(n, b) < 1e-9


def test_solution_reproduces_rwgs_equilibrium_constant(rwgs_solution):
    """The strongest check: the converged gas must satisfy Q == Kp."""
    T, n, _b, _ = rwgs_solution
    ng = sum(n[i] for i in gas_idx)
    y = {s: n[SPECIES.index(s)] / ng for s in ('CO2', 'H2', 'CO', 'H2O')}
    Q = (y['CO'] * y['H2O']) / (y['CO2'] * y['H2'])
    # SLSQP settles to about a tenth of a percent on the quotient with nine
    # species; the seven-species problem used to reach 1e-4. That is solver
    # convergence, not physics - the mole fractions themselves are stable.
    assert Q == pytest.approx(Kp(RWGS, T), rel=5e-3)


def test_rwgs_feed_does_not_reduce_titania(rwgs_solution):
    """A 1:1 CO2/H2 feed leaves rutile untouched."""
    _T, n, _b, n_Ti = rwgs_solution
    # the residue is the 1e-30 lower bound working its way up, not reduction
    assert ti3_percent(n, n_Ti) < 1e-3


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
            assert residual(trial, b) < 1e-9
            assert G_total(trial, T) >= G0 - 1e-12


def test_methane_opens_the_magneli_branch():
    """Pure CH4 at 900 C reduces a large fraction of the Ti, unlike RWGS."""
    T = 1173.15
    n0_TiO2 = 0.100 / 79.87
    n0g = moles_of_gas(1.0, 25.0, T)
    b = np.array([n0g, 4 * n0g, 2 * n0_TiO2, n0_TiO2])
    inits = [np.array([n0g * gf * 0.3, n0g * gf, n0g * gf, n0g * gf * 0.1,
                       n0g * (1 - gf)] + list(sv))
             for gf in (0.05, 0.4, 0.9) for sv in phase_seeds(n0_TiO2)]
    n = solve(b, T, inits)
    assert n is not None
    assert residual(n, b) < 1e-10
    assert ti3_percent(n, n0_TiO2) > 20.0
    # it stops at Ti4O7: the deeper phases never turn up
    split = phase_split(n, n0_TiO2)
    assert split['Ti3O5'] < 0.01 and split['Ti2O3'] < 0.01
