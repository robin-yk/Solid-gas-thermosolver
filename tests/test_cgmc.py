"""The coarse-grained KMC transport layer against its anchors.

Four properties pin the implementation, mirroring the paper it follows
(Katsoulakis & Vlachos, J. Chem. Phys. 119, 9412 (2003)):
detailed balance - the closed system must land on the equilibrium
workspace's split; the dilute limit must match a 50-digit matrix
exponential of the linear master equation; observables must be invariant
under the level of coarse-graining; and mass must be conserved to
machine precision at every step.
"""

import json
import pathlib

import pytest

from solidgas import cgmc as C
from solidgas import statmech as S

ROOT = pathlib.Path(__file__).resolve().parent.parent
REF = json.loads((ROOT / 'data' / 'reference_statmech.json').read_text())
P = S.load_params()

NO_FILL = {'A_eff_per_s_atm': {'surface': 0.0}}
ALL_BULK = {'surface': 0.0, 'subsurface': 0.0, 'bulk': 1.0}


def _linear_sys(eps_ctrl=0.02):
    s = REF['cgmc_linear']['spec']
    kt = S.KB_EV * (s['T_C'] + 273.15)
    return {'m': 3, 'q': [float(x) for x in s['q']], 'eps_c': s['eps'],
            'lam': s['lam'], 'beta': 1.0 / kt, 'k_fill': s['k_fill'],
            'eta0': s['eta0'], 'V0': sum(s['eta0']), 'R_nm': 3.0,
            'a_nm': 1.0, 'centers': [2.5, 1.5, 0.5],
            'T_reox_C': s['T_C'], 'E_m_eV': 0.0, 'eps_ctrl': eps_ctrl,
            'column_counts': 1e5}


def test_dilute_limit_matches_the_matrix_exponential_oracle():
    """The linear master equation, expm'd at 50 digits by the oracle."""
    for row in REF['cgmc_linear']['rows']:
        r = C.run(_linear_sys(), row['t_s'], n_out=8)
        for a, b in zip(r['final_eta'], row['eta']):
            assert a == pytest.approx(b, rel=6e-3, abs=1e-6)


def test_closed_system_lands_on_the_equilibrium_split():
    """Detailed balance in practice: start with every vacancy in the
    core, shut the gas channel, and the network must relax to the exact
    product-binomial equilibrium of its own Hamiltonian - every cell on
    the common-mu site-exclusion isotherm over the cell capacities.

    (This is the (1x1) three-class equilibrium; the particle workspace's
    phase-aware partition adds the reconstruction and CS phases ON TOP of
    this Hamiltonian, so it is not the comparison point here.)"""
    sys0 = C.build(P, kin=NO_FILL, dist_fractions=ALL_BULK)
    r = C.run(sys0, 300.0)
    kt = 1.0 / sys0['beta']
    v0 = sum(sys0['eta0'])

    def total(mu):
        return sum(q * S.fd(mu, e, kt)
                   for q, e in zip(sys0['q'], sys0['eps_c']))
    lo, hi = -3.0, 3.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if total(mid) < v0:
            lo = mid
        else:
            hi = mid
    mu = 0.5 * (lo + hi)
    for got, q, e in zip(r['final_eta'], sys0['q'], sys0['eps_c']):
        assert got / q == pytest.approx(S.fd(mu, e, kt), abs=5e-3)
    assert abs(r['mass_check']) < 1e-12


def test_mass_is_conserved_to_machine_precision():
    for cg in ({}, {'E_m_bulk_eV': 1.9}, {'E_m_bulk_eV': 2.1}):
        r = C.run(C.build(P, cg=cg), 300.0)
        assert abs(r['mass_check']) < 1e-12
        assert r['completed']


def test_neutral_barrier_transport_is_not_rate_limiting():
    """With Iddir's neutral-vacancy barrier the particle drains long
    before the 300 s exposure - the finding that turns a measured partial
    recovery into an effective-barrier readout."""
    r = C.run(C.build(P), 300.0)
    assert r['recovered_fraction'] > 0.99
    t50 = next(t for t, v in zip(r['t_s'], r['recovery']) if v > 0.5)
    assert t50 < 1e-3


def test_recovery_falls_monotonically_with_the_barrier():
    recs = [C.run(C.build(P, cg={'E_m_bulk_eV': em}),
                  300.0)['recovered_fraction']
            for em in (0.65, 1.7, 1.85, 1.96, 2.1)]
    assert all(a >= b - 1e-9 for a, b in zip(recs, recs[1:]))
    assert recs[0] > 0.99
    assert recs[-1] < 0.40


def test_observables_are_invariant_under_the_coarse_graining_level():
    """The paper's collapse: changing the number of coarse cells must not
    change the physics."""
    a = C.run(C.build(P, cg={'E_m_bulk_eV': 1.9}),
              300.0)['recovered_fraction']
    b = C.run(C.build(P, cg={'E_m_bulk_eV': 1.9, 'n_cells': 56}),
              300.0)['recovered_fraction']
    c = C.run(C.build(P, cg={'E_m_bulk_eV': 1.9, 'n_cells': 16}),
              300.0)['recovered_fraction']
    assert b == pytest.approx(a, abs=0.03)
    assert c == pytest.approx(a, abs=0.05)


def test_step_control_is_converged():
    a = C.run(C.build(P, cg={'E_m_bulk_eV': 1.9}),
              300.0)['recovered_fraction']
    b = C.run(C.build(P, cg={'E_m_bulk_eV': 1.9, 'eps_ctrl': 0.02}),
              300.0)['recovered_fraction']
    assert b == pytest.approx(a, abs=0.01)


def test_effective_barrier_inversion_round_trips():
    em = C.effective_barrier(P, 0.40, 300.0)
    assert 1.5 < em < 2.2
    r = C.run(C.build(P, cg={'E_m_bulk_eV': em}), 300.0)
    assert r['recovered_fraction'] == pytest.approx(0.40, abs=0.02)


def test_stochastic_mode_agrees_with_the_mean_where_the_leap_is_valid():
    """tau-leap against the deterministic mean, on both sides of validity.

    Where the leap is valid the two agree cell for cell. Where it is not,
    they do not, and the way to show the difference is sampling and not a
    bug is to close it by adding columns.

    Which configurations fall on which side moved when the compensating
    polarons went into the free energy: the subsurface cell now starts
    with 3.6 umol-O/g instead of 52, so a slow barrier that used to move
    enough material to leap over now does not move a single discrete
    unit in the time given."""
    def cfg(em, cols):
        return dict(kin=NO_FILL,
                    cg={'E_m_bulk_eV': em, 'n_cells': 8,
                        'stochastic_column_counts': cols})

    # valid: fast enough that every cell exchanges many units
    kw = cfg(1.0, 50000)
    det = C.run(C.build(P, **kw), 60.0)
    st = C.run(C.build(P, **kw), 60.0, mode='stoch', seed=902)
    assert st['completed']
    assert abs(st['mass_check']) < 1e-12
    td, ts = sum(det['final_eta']), sum(st['final_eta'])
    for i, (a, b) in enumerate(zip(det['final_eta'], st['final_eta'])):
        assert b / ts == pytest.approx(a / td, abs=0.01), i
    st2 = C.run(C.build(P, **kw), 60.0, mode='stoch', seed=902)
    assert st2['final_eta'] == st['final_eta']

    # not valid: a slow barrier, where the mean creeps and the leap
    # cannot. The gap has to close as the columns are refined, which is
    # what says it is discreteness rather than a broken rate.
    gaps = []
    for cols in (50000, 5000000, 50000000):
        k = cfg(1.9, cols)
        d = C.run(C.build(P, **k), 60.0)
        r = C.run(C.build(P, **k), 60.0, mode='stoch', seed=902)
        gaps.append(abs(r['final_eta'][1] / sum(r['final_eta'])
                        - d['final_eta'][1] / sum(d['final_eta'])))
    assert gaps[0] > 0.05, gaps
    assert all(y < x for x, y in zip(gaps, gaps[1:])), gaps
    assert gaps[-1] < 0.5 * gaps[0], gaps


def test_the_kinetic_layer_relaxes_to_a_measure_the_thermodynamics_no_longer_uses():
    """An open item, pinned so it cannot be forgotten.

    The exchange rates satisfy detailed balance against a product-binomial
    measure on the oxygen sublattice with H = sum eps_k eta_k. That is the
    equilibrium of a NEUTRAL vacancy. The thermodynamic model is no longer
    that - it counts the cation sublattice too - so the long-time limit of
    this layer and the equilibrium partition disagree, and the size of the
    disagreement is the same cube that separates the two isotherms.

    This test does not assert they agree. It asserts they do not, so that
    the day someone reconciles them the gate fails and gets updated
    rather than the inconsistency being rediscovered. The size of it is
    visible in the test above: run the closed system to its long-time
    limit and the subsurface cell settles near 0.54 of the inventory,
    where the equilibrium partition puts 0.038."""
    kt = S.KB_EV * 873.15
    eps = S.energetics(P)['eps']
    mu = S.distribute(P, 600.0, 95.0, S.geometry(P))['mu_V_eV']
    neutral = S.fd(mu, eps[2], kt)
    compensated = S.theta_of_mu(mu, eps[2], kt)
    assert neutral != pytest.approx(compensated, rel=0.1)
    note = P['kinetics_cgmc']['note']
    assert 'detailed balance' in note

def test_geometry_capacities_are_consistent_with_the_per_gram_model():
    sys0 = C.build(P)
    g = S.geometry(P)
    tot = sum(sys0['q'])
    assert sys0['q'][0] / tot == pytest.approx(
        g['N_s'] / g['N_O_total'], rel=0.02)
    assert sys0['q'][1] / tot == pytest.approx(
        g['N_ss'] / g['N_O_total'], rel=0.02)
    assert sum(sys0['eta0']) == pytest.approx(sys0['V0'], rel=1e-12)


def test_provenance_is_pinned():
    note = P['kinetics_cgmc']['note']
    prov = P['kinetics_cgmc']['provenance']
    assert 'detailed balance' in note
    assert 'Katsoulakis' in prov and '10.1063/1.1616513' in prov
    assert 'Iddir' in prov
