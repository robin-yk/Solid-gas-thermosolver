"""The three-class construction, and what has to stay true about it.

This engine splits the measured inventory over a bridging row, one
subsurface slab and the bulk, at a common vacancy chemical potential,
with the (1x2) reconstruction and the shear ceiling as competing phases.
Its free energy is not its own: it is imported from
solidgas/particle/freeenergy.py, so there is one compensated relation in
the repository and the coarse partition cannot drift from the
depth-resolved one.

This file replaces tests/test_statmech_layers.py, whose subject was
resolving the shell into more classes with an assumed depth profile and
a mean-field crowding penalty.  Both are gone.  The profile shape is
derived rather than assumed in solidgas/particle/, and the penalty was
standing in for the cation-sublattice entropy that is now in the free
energy - with the compensation in place, resolving the shell into three
layers moved the bulk inventory by two percent rather than the thirty it
used to.  The gates below are the ones from that file whose subject
survived, brought over with their expectations recomputed.
"""

import json
import math
import os

import pytest

from solidgas import statmech as S
from solidgas.particle import freeenergy as F

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF = json.load(open(os.path.join(ROOT, 'data', 'reference_statmech.json')))

T_C = 600.0
VO = 95.0
D_UM = 0.90


@pytest.fixture(scope='module')
def p():
    return S.load_params()


def geom(p, **kw):
    kw.setdefault('d_um', D_UM)
    kw.setdefault('ss_layers', 1)
    return S.geometry(p, **kw)


# ------------------------------------------------- one free energy, shared

def test_the_engine_does_not_own_a_second_free_energy(p):
    """statmech re-exports particle's relation rather than repeating it.

    Two implementations of the same equation is how a repository ends up
    with two answers, and the redox model reads this one while the
    workspace card reads the other."""
    assert S.mu_of_theta is F.mu_of_theta
    assert S.theta_of_mu is F.theta_of_mu
    assert S.THETA_CAP == F.THETA_CAP
    src = open(os.path.join(ROOT, 'solidgas', 'statmech.py')).read()
    assert 'def mu_of_theta' not in src
    assert 'def theta_of_mu' not in src


def test_the_retired_machinery_is_gone_and_says_why(p):
    """The layer resolution and the crowding penalty are not switched
    off, they are removed - along with the dataset blocks that fed them,
    which leave a note behind rather than a hole."""
    src = open(os.path.join(ROOT, 'solidgas', 'statmech.py')).read()
    for gone in ('omega', 'resolved_layers', 'eps_layers'):
        assert gone not in src, 'statmech.py still carries ' + gone
    ve = p['vacancy_energetics']
    assert 'layer_profile' not in ve
    assert 'coverage_penalty' not in ve
    assert 'retired' in ve and 'exp(-z/xi)' in ve['retired']
    assert 'resolved_layers' not in p['defaults']


# ------------------------------------------------------- the pinned answer

def test_the_flagship_partition_is_pinned(p):
    """The shipped condition, to the last digit the oracle prints.

    These numbers moved once, when the compensating polarons went into
    the free energy. They should not move again without a reason written
    down next to them."""
    r = S.distribute(p, T_C, VO, geom(p))
    assert round(r['mu_V_eV'], 9) == 0.237153143
    assert round(r['umol_g']['surface'], 4) == 6.7785
    assert round(r['umol_g']['subsurface'], 4) == 3.5894
    assert round(r['umol_g']['bulk'], 4) == 84.6321
    assert r['regime'] == 'reconstructed'
    assert r['extended_defects_umol_g'] == 0.0


def test_every_frozen_row_still_comes_out_of_the_engine(p):
    """The oracle rows, through the production path."""
    f = REF['flagship']
    r = S.distribute(p, f['T_C'], f['VO_total_umol_g'], geom(p),
                     preset=f['preset'])
    assert r['mu_V_eV'] == pytest.approx(f['mu_V_eV'], abs=1e-11)
    for k in ('surface', 'subsurface', 'bulk'):
        assert r['theta'][k] == pytest.approx(f['theta'][k], rel=1e-11)
        assert r['umol_g'][k] == pytest.approx(f['umol_g'][k], rel=1e-11)
    assert r['extended_defects_umol_g'] == pytest.approx(
        f['extended_defects_umol_g'], abs=1e-9)


def test_the_inventory_residual_is_numerical_zero(p):
    """Every branch of the construction, closed."""
    g = geom(p)
    worst = 0.0
    for vo in (0.5, 2.31, 4.0, 6.79, 30.0, 34.0, 95.0, 155.0, 400.0):
        r = S.distribute(p, T_C, vo, g)
        tot = (r['umol_g']['surface'] + r['umol_g']['subsurface']
               + r['umol_g']['bulk'] + r['extended_defects_umol_g'])
        worst = max(worst, abs(tot - vo) / vo)
    assert worst < 1e-10, worst


def test_the_solver_inverts_its_own_chemical_potential(p):
    kt = S.KB_EV * (T_C + 273.15)
    worst, checked = 0.0, 0
    for eps in (0.0, 0.59, 1.05, 1.31):
        for mu in (-0.4, 0.0, 0.5, 0.85, 1.2, 2.0):
            th = S.theta_of_mu(mu, eps, kt)
            # both 1 - theta and cap - theta are subtractions, so the
            # window is where the check measures the solver rather than
            # its own arithmetic
            if not 1e-4 < th < S.THETA_CAP - 1e-6:
                continue
            worst = max(worst, abs(S.mu_of_theta(th, eps, kt) - mu))
            checked += 1
    assert checked > 8, checked
    assert worst < 1e-11, worst


def test_the_implicit_root_matches_the_oracle(p):
    """Two algorithms on the same equation, not one checked against itself.

    The oracle inverts by bisecting a logistic coordinate at fifty digits;
    the engine uses analytic one-sided brackets on two branches."""
    rows = REF['theta_of_mu']
    assert len(rows) == 45, len(rows)
    worst = 0.0
    for r in rows:
        kt = S.KB_EV * (r['T_C'] + 273.15)
        got = S.theta_of_mu(r['mu_eV'], r['eps_eV'], kt)
        worst = max(worst, abs(got - r['theta'])
                    / max(r['theta'], 1e-300))
    assert worst < 1e-13, worst


# ------------------------------------------------- the two pinned windows

def test_mu_is_continuous_across_the_cs_boundary(p):
    """Approaching the ceiling from below must land on the pinned value."""
    g = geom(p)
    b = S.phase_boundaries(p, T_C, g)
    cs = b['VO_cs_umol_g']
    below = S.distribute(p, T_C, cs * (1 - 1e-9), g)
    at = S.distribute(p, T_C, cs, g)
    above = S.distribute(p, T_C, cs * 1.5, g)
    assert below['mu_V_eV'] == pytest.approx(b['mu_cs_eV'], abs=1e-9)
    assert at['mu_V_eV'] == pytest.approx(b['mu_cs_eV'], abs=1e-9)
    assert above['mu_V_eV'] == pytest.approx(b['mu_cs_eV'], abs=1e-12)
    assert above['extended_defects_umol_g'] > 0.0
    assert below['extended_defects_umol_g'] == 0.0


def test_the_cs_boundary_is_the_bulk_at_the_declared_solubility(p):
    """mu_cs is a condition on the bulk class and nothing else."""
    g = geom(p)
    kt = S.KB_EV * (T_C + 273.15)
    th_sol = p['saturation']['x_max_shear'] / 2.0
    en = S.energetics(p)
    b = S.phase_boundaries(p, T_C, g)
    assert b['mu_cs_eV'] == pytest.approx(
        S.mu_of_theta(th_sol, en['eps'][2], kt), abs=1e-15)
    assert S.theta_of_mu(b['mu_cs_eV'], en['eps'][2], kt) == pytest.approx(
        th_sol, abs=1e-12)


def test_mu_is_pinned_across_the_surface_coexistence(p):
    """A first-order transition holds mu fixed while the lever arm moves."""
    g = geom(p)
    b = S.phase_boundaries(p, T_C, g)
    lo, hi = b['VO_onset_umol_g'], b['VO_recon_complete_umol_g']
    assert hi > lo
    phis = []
    for f in (0.001, 0.25, 0.5, 0.75, 0.999):
        r = S.distribute(p, T_C, lo + f * (hi - lo), g)
        assert r['regime'] == 'surface_coexistence', f
        assert r['mu_V_eV'] == pytest.approx(b['mu_t_eV'], abs=1e-12)
        phis.append(r['surface_phase']['phi_reconstructed'])
    assert all(y > x for x, y in zip(phis, phis[1:])), phis
    assert phis[0] < 0.01 and phis[-1] > 0.99


def test_the_shell_thickness_is_still_a_declared_choice(p):
    """Thickening the shell moves inventory out of the bulk, and the
    engine gives every slab in it the first-subsurface energy - which is
    exactly why more than one slab is coarse, and why the depth-resolved
    engine exists."""
    thin = S.distribute(p, T_C, VO, geom(p, ss_layers=1))
    thick = S.distribute(p, T_C, VO, geom(p, ss_layers=4))
    assert thick['umol_g']['subsurface'] > thin['umol_g']['subsurface']
    assert thick['umol_g']['bulk'] < thin['umol_g']['bulk']
    assert thick['theta']['subsurface'] < thin['theta']['subsurface']
