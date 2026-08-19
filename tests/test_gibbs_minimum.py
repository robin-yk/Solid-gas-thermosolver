"""The production answer IS the Gibbs minimum - measured, not asserted.

Three legs, each attacking the claim a different way:

  1. Rank every feasible candidate by its actual Gibbs energy at oracle
     precision: the KKT-clean candidate must be the lowest. This is the
     direct form of the claim, and it also documents *why* the production
     decision is made by KKT signs instead - under sealed N2 the runner-up
     sits ~1e-24 of |G| above the winner, a gap double precision cannot rank.
  2. Perturb the production solution along feasible reaction directions,
     both signs, several magnitudes: G must rise every time. A point where
     G rises in every feasible direction is a minimum by definition.
  3. Re-solve the mass/volume ladder sweep against the 80-digit oracle.
     These points make the deep single-phase windows (Ti10O19, Ti5O9,
     Ti4O7, Ti3O5) and a host-free pair (Ti8O15+Ti6O11) win - branches of
     the enumeration that the six canonical feeds never exercise.

The oracle module is imported by the tests only; the oracle itself still
imports nothing from production, and test_activeset.py enforces the reverse
direction with an AST scan.
"""

import math

import pytest
from mpmath import mp, mpf, log as mlog

import oracle_tio as O
from solidgas import activeset as A

TEMPS_SWEEP_MG = [100, 20, 8.5, 8, 5, 2, 0.8, 0.3, 0.1]
TEMPS_SWEEP_ML = [100, 300, 1000, 3000, 10000]


def _oracle_point(feed, T_C, V=25.0, mass=0.100):
    T = mpf(str(T_C)) + mpf('273.15')
    pt = O.build_point(feed, T, V=mpf(repr(float(V))),
                       mass=mpf(repr(float(mass))))
    mu = O.mu_table(T)
    passing, rows = [], []
    for cand in O.candidates():
        r = O.solve_candidate(cand, pt, mu)
        rows.append(r)
        if r.get('feasible') and r.get('kkt_ok'):
            passing.append(r)
    assert passing, 'oracle found no KKT-clean candidate'
    passing.sort(key=lambda r: len(r['active']))
    return passing[0], rows, pt, mu, T


# ------------------------------------------------- 1. G-ranking, 80 digits


@pytest.mark.parametrize('feed,T_C', [
    ({'N2': 1}, 900),
    ({'H2': 1}, 1200),
    ({'CO2': 1, 'H2': 1}, 900),
])
def test_kkt_winner_is_the_lowest_G_feasible_candidate(feed, T_C):
    win, rows, pt, mu, T = _oracle_point(feed, T_C)
    ranked = []
    for r in rows:
        if not r.get('feasible'):
            continue
        m = {s: (v if v is not None else mpf(0)) for s, v in r['m'].items()}
        g = O.g_total(r['n_gas'], m, mu, T, pt['P'])
        ranked.append((g, r['kkt_ok'], '+'.join(r['active'])))
    ranked.sort(key=lambda x: x[0])
    assert ranked[0][1], f'lowest-G candidate {ranked[0][2]} is not KKT-clean'
    assert ranked[0][2] == '+'.join(win['active'])
    for g, kkt, name in ranked[1:]:
        assert g > ranked[0][0], f'{name} ties the winner without degeneracy'


def test_the_pure_n2_gap_is_why_kkt_decides_not_G():
    """The runner-up under sealed N2 sits ~1e-24 of |G| above the winner:
    real at 80 digits, unrankable in double. The production rule 'decide by
    KKT signs, never by comparing G' exists because of this number."""
    win, rows, pt, mu, T = _oracle_point({'N2': 1}, 900)
    gs = []
    for r in rows:
        if not r.get('feasible'):
            continue
        m = {s: (v if v is not None else mpf(0)) for s, v in r['m'].items()}
        gs.append(O.g_total(r['n_gas'], m, mu, T, pt['P']))
    gs.sort()
    gap = gs[1] - gs[0]
    assert gap > 0
    assert gap / abs(gs[0]) < mpf('1e-16'), \
        'the gap became double-rankable; revisit the decision-rule note'


# ------------------------------------ 2. G rises in every feasible direction


def test_G_rises_in_every_feasible_direction_around_the_answer():
    T = mpf('1473.15')
    P = mpf(1)
    r = A.solve({'H2': 1}, 1473.15)
    mu = O.mu_table(T)
    n0 = {s: mpf(repr(v)) for s, v in r['species_moles'].items() if v > 0}
    n0['TiO2'] = mpf(repr(r['solid_moles']['TiO2']))
    n0['Ti10O19'] = mpf(repr(r['solid_moles']['Ti10O19']))

    def G(n):
        Ng = sum(n[s] for s in ('H2', 'H2O', 'O2') if n.get(s, 0) > 0)
        g = mpf(0)
        for s, v in n.items():
            if v < 0:
                return None
            if v == 0:
                continue
            g += v * (mu[s] + O.R_KJ * T * (mlog(v / Ng * P)
                                            if s in ('H2', 'H2O', 'O2')
                                            else 0))
        return g

    def shifted(d, eps):
        n = dict(n0)
        for s, c in d.items():
            n[s] = n.get(s, mpf(0)) + eps * c
        return n

    G0 = G(n0)
    m_eq, o_eq = n0['Ti10O19'], n0['O2']
    directions = [
        # 10 TiO2 + H2 -> Ti10O19 + H2O : the reduction extent
        ({'TiO2': -10, 'Ti10O19': 1, 'H2': -1, 'H2O': 1},
         [m_eq * f for f in (mpf('-0.5'), mpf('-0.01'), mpf('0.01'),
                             mpf('0.5'), mpf('3'))]),
        # H2O -> H2 + 1/2 O2 : oxygen redistribution inside the gas
        ({'H2': 1, 'H2O': -1, 'O2': mpf(1) / 2},
         [o_eq * f for f in (mpf('-0.9'), mpf('0.5'), mpf('50'))]),
        # 10 TiO2 -> Ti10O19 + 1/2 O2 : reduction without the couple
        ({'TiO2': -10, 'Ti10O19': 1, 'O2': mpf(1) / 2},
         [o_eq * f for f in (mpf('-0.9'), mpf('0.9'), mpf('20'))]),
    ]
    checked = 0
    for d, steps in directions:
        for eps in steps:
            g = G(shifted(d, eps))
            assert g is not None, 'perturbation left the feasible set'
            assert g > G0, f'G fell along {d} at eps={mp.nstr(eps, 3)}'
            checked += 1
    assert checked >= 11


# --------------------------------------- 3. the ladder sweep vs the oracle


@pytest.fixture(scope='module')
def sweep_pairs():
    out = []
    for mg in TEMPS_SWEEP_MG:
        p = A.solve({'H2': 1}, 1473.15, mass_g=mg / 1000)
        w, rows, pt, mu, T = _oracle_point({'H2': 1}, 1200, mass=mg / 1000)
        out.append((f'{mg} mg / 25 mL', p, w, pt))
    for mL in TEMPS_SWEEP_ML:
        p = A.solve({'H2': 1}, 1473.15, V_cm3=mL)
        w, rows, pt, mu, T = _oracle_point({'H2': 1}, 1200, V=mL)
        out.append((f'100 mg / {mL} mL', p, w, pt))
    return out


def test_sweep_active_sets_match_the_oracle(sweep_pairs):
    for label, p, w, pt in sweep_pairs:
        assert p['active_condensed_phases'] == w['active'], label


def test_sweep_reduction_matches_the_oracle(sweep_pairs):
    for label, p, w, pt in sweep_pairs:
        m = {s: (v if v is not None else mpf(0)) for s, v in w['m'].items()}
        red_o = float(100 * sum(O.TI3[s] * v for s, v in m.items())
                      / pt['n_Ti'])
        assert p['reduced_pct'] == pytest.approx(red_o, rel=1e-9), label


def test_sweep_gas_and_reduced_costs_match_the_oracle(sweep_pairs):
    for label, p, w, pt in sweep_pairs:
        tot = sum(w['n_gas'].values())
        for g, v in p['gas_fractions'].items():
            ov = float(w['n_gas'][g] / tot)
            if ov > 1e-25:
                assert v == pytest.approx(ov, rel=1e-9), (label, g)
        for s, rc in p['inactive_phase_reduced_costs'].items():
            v = rc['per_mol_O_kJ']
            if v is None or math.isinf(v):
                continue
            assert abs(v - float(w['r_costs'][s]['per_mol_O'])) < 1e-8, \
                (label, s)


def test_sweep_exercises_the_branches_the_canonical_points_never_win(sweep_pairs):
    """Deep single-phase windows and a pair without the host: the half of
    the enumeration the six reference feeds leave dormant."""
    winners = {tuple(p['active_condensed_phases']) for _, p, _, _ in sweep_pairs}
    for needed in (('Ti10O19',), ('Ti8O15',), ('Ti5O9',), ('Ti4O7',),
                   ('Ti3O5',), ('Ti8O15', 'Ti6O11'), ('TiO2', 'Ti10O19')):
        assert needed in winners, f'sweep no longer reaches {needed}'
    for label, p, w, pt in sweep_pairs:
        assert p['method'] == 'gibbs_min', label
        for s, v in p['solid_moles'].items():
            if s not in p['active_condensed_phases']:
                assert v == 0.0, (label, s)
