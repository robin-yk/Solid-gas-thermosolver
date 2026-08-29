"""The layer-resolved free energy, and the eight things it must not break.

Stage 1 replaces the site-exclusion isotherm with an implicit occupancy,

    mu = eps_i + omega_i theta_i + kT ln(theta_i/(1-theta_i))

and lets the near-surface shell resolve into more than one class. Both are
off by default in data/rutile_dft.json, so the shipped partition is the one
that was already certified; these gates hold the reduction exact and check
that the new degrees of freedom do what the physics says they should.

The surface phase ladder and the surface Monte Carlo are deliberately
untouched at this stage, so any change in a deeper class is attributable to
the layer interaction alone and to nothing else."""

import json
import math
import os

import pytest

from solidgas import statmech as S

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


# ------------------------------------------------------ 1. the legacy triple

def test_the_flagship_partition_is_unchanged_to_the_last_digit(p):
    """6.7785 / 51.7557 / 36.4658 umol-O/g, and the same mu that made them.

    The default dataset ships one resolved layer and no coverage penalty,
    so the new code path has to walk through the old answer exactly - not
    to a tolerance, since a tolerance would hide a solver that is merely
    close."""
    r = S.distribute(p, T_C, VO, geom(p))
    assert round(r['umol_g']['surface'], 4) == 6.7785
    assert round(r['umol_g']['subsurface'], 4) == 51.7557
    assert round(r['umol_g']['bulk'], 4) == 36.4658
    ref = REF['flagship']
    assert r['mu_V_eV'] == pytest.approx(float(ref['mu_V_eV']), abs=1e-12)
    for k in ('surface', 'subsurface', 'bulk'):
        assert r['umol_g'][k] == pytest.approx(
            float(ref['umol_g'][k]), abs=1e-12), k


# ------------------------------------------- 2. the frozen references hold

def test_every_frozen_row_still_comes_out_of_the_new_path(p):
    """The oracle rows the legacy engine was certified against.

    Checked here as well as in test_statmech.py because those tests could
    in principle be satisfied by a code path this change no longer takes;
    this one calls the same public entry points the workspace calls."""
    for row in REF['analytic_grid']:
        g = geom(p)
        r = S.distribute(p, float(row['T_C']),
                         float(row['VO_total_umol_g']), g,
                         preset=row['preset'])
        assert r['mu_V_eV'] == pytest.approx(float(row['mu_V_eV']), abs=1e-11)
        assert r['regime'] == row['regime']
        for k in ('surface', 'subsurface', 'bulk'):
            assert r['umol_g'][k] == pytest.approx(
                float(row['umol_g'][k]), abs=1e-9), (row['T_C'], k)
    for row in REF['shell_ladder']:
        g = geom(p, ss_layers=int(row['ss_layers']))
        r = S.distribute(p, float(row['T_C']),
                         float(row['VO_total_umol_g']), g,
                         preset=row['preset'])
        assert r['mu_V_eV'] == pytest.approx(float(row['mu_V_eV']), abs=1e-11)
        for k in ('surface', 'subsurface', 'bulk'):
            assert r['umol_g'][k] == pytest.approx(
                float(row['umol_g'][k]), abs=1e-9), (row['ss_layers'], k)


# ------------------------------------------------- 3. the inventory closes

def test_the_inventory_residual_is_numerical_zero(p):
    """Every branch of the construction, at every resolution, closed."""
    worst = 0.0
    for n in (1, 2, 3):
        for cp in ('off', 'on'):
            g = geom(p, resolved_layers=n)
            for vo in (0.5, 2.31, 4.0, 6.79, 30.0, 95.0, 155.0, 400.0):
                r = S.distribute(p, T_C, vo, g, omega=cp)
                tot = (r['umol_g']['surface'] + r['umol_g']['subsurface']
                       + r['umol_g']['bulk']
                       + r['extended_defects_umol_g'])
                worst = max(worst, abs(tot - vo) / vo)
    assert worst < 1e-10, worst


# --------------------------------------- 4. crowding pushes occupancy down

def test_a_larger_penalty_empties_the_class_it_is_applied_to(p):
    """mu = eps + omega theta + kT logit(theta) is increasing in omega at
    fixed mu, so theta must fall. Checked on the solver, not the algebra."""
    kt = S.KB_EV * (T_C + 273.15)
    last = 1.0
    for w in (0.0, 0.1, 0.3, 0.6, 1.0):
        th = S.theta_of_mu(0.80, 0.59, w, kt)
        assert th < last, w
        last = th
    # and through the whole partition: the shell gives inventory back
    g = geom(p)
    off = S.distribute(p, T_C, VO, g, omega='off')
    on = S.distribute(p, T_C, VO, g, omega='on')
    assert on['theta']['subsurface'] < off['theta']['subsurface']
    assert on['umol_g']['bulk'] > off['umol_g']['bulk']
    assert on['mu_V_eV'] > off['mu_V_eV'], \
        'a crowded shell has to be paid for with a higher potential'


# ------------------------------------------------- 5. capacity is conserved

def test_resolving_the_shell_moves_no_sites(p):
    """Splitting a class into layers redistributes inventory, never sites.

    The site counts are geometry; only where the inventory sits is model."""
    base = geom(p, ss_layers=3, resolved_layers=1)
    split = geom(p, resolved_layers=3)
    assert split['N_s'] == base['N_s']
    assert sum(split['N_layers']) == pytest.approx(base['N_ss'], rel=1e-15)
    assert split['N_b'] == pytest.approx(base['N_b'], rel=1e-15)
    for g in (base, split):
        assert (g['N_s'] + g['N_ss'] + g['N_b']) == pytest.approx(
            g['N_O_total'], rel=1e-15)
    # and the reported class total is the sum of the resolved layers
    r = S.distribute(p, T_C, VO, split, omega='on')
    assert r['umol_g']['subsurface'] == pytest.approx(
        sum(l['umol_g'] for l in r['layers']), rel=1e-14)


# ------------------------------- 6. omega = 0 is the site-exclusion isotherm

def test_no_penalty_is_the_analytic_isotherm(p):
    """Bit for bit at omega exactly zero, and converging to it as omega -> 0.

    The zero case returns fd() rather than bisecting to it, so the legacy
    numbers cannot move. The limit case proves the bisection agrees with
    the closed form rather than merely being short-circuited past it: the
    root is bracketed inside a window of width omega/kT in the logit, so
    the two answers cannot differ by more than that, and the gap has to
    shrink with omega rather than sit at a floor."""
    for t_c in (300.0, 600.0, 1000.0):
        kt = S.KB_EV * (t_c + 273.15)
        for eps in (0.0, 0.59, 1.31):
            for mu in (-0.5, 0.0, 0.4, 0.82, 1.5):
                exact = S.fd(mu, eps, kt)
                assert S.theta_of_mu(mu, eps, 0.0, kt) == exact
                prev = None
                for w in (1e-6, 1e-9, 1e-12):
                    gap = abs(S.theta_of_mu(mu, eps, w, kt) - exact)
                    assert gap <= w / kt, (t_c, eps, mu, w, gap)
                    if prev is not None:
                        assert gap <= prev
                    prev = gap


def test_the_solver_inverts_its_own_chemical_potential(p):
    """The strongest statement available without a second implementation:
    put theta back into mu = eps + omega theta + kT logit(theta)."""
    kt = S.KB_EV * (T_C + 273.15)
    worst = 0.0
    checked = 0
    for w in (0.15, 0.3, 0.6, 1.0):
        for eps in (0.0, 0.59, 1.05, 1.31):
            for mu in (-0.4, 0.0, 0.5, 0.85, 1.2, 2.0):
                th = S.theta_of_mu(mu, eps, w, kt)
                # 1 - theta is a subtraction, so the round trip runs out of
                # digits before the solver does; the window is where the
                # check measures the solver and not its own arithmetic
                if not 1e-4 < th < 1.0 - 1e-4:
                    continue
                back = eps + w * th + kt * math.log(th / (1.0 - th))
                worst = max(worst, abs(back - mu))
                checked += 1
    assert checked > 20, checked
    assert worst < 1e-11, worst


# ----------------------------------------- 7. mu is continuous at CS onset

def test_mu_is_continuous_across_the_cs_boundary(p):
    """Approaching the ceiling from below must land on the pinned value.

    The construction pins mu at the rutile / CS coexistence above the
    ceiling; if the branch below it did not converge to the same number
    there would be a jump in a quantity that has no reason to jump."""
    for n in (1, 3):
        for cp in ('off', 'on'):
            g = geom(p, resolved_layers=n)
            b = S.phase_boundaries(p, T_C, g, omega=cp)
            cs = b['VO_cs_umol_g']
            below = S.distribute(p, T_C, cs * (1 - 1e-9), g, omega=cp)
            at = S.distribute(p, T_C, cs, g, omega=cp)
            above = S.distribute(p, T_C, cs * 1.5, g, omega=cp)
            assert below['mu_V_eV'] == pytest.approx(
                b['mu_cs_eV'], abs=1e-9), (n, cp)
            assert at['mu_V_eV'] == pytest.approx(b['mu_cs_eV'], abs=1e-9)
            assert above['mu_V_eV'] == pytest.approx(b['mu_cs_eV'], abs=1e-12)
            assert above['extended_defects_umol_g'] > 0.0
            assert below['extended_defects_umol_g'] == 0.0


def test_the_cs_boundary_carries_the_bulk_interaction(p):
    """mu_cs is a condition on the bulk class, so omega_bulk belongs in it.

    Leaving it out would put the ceiling at the non-interacting potential
    while the bulk sitting at it was interacting - the two would disagree
    about what x_sol means."""
    g = geom(p)
    kt = S.KB_EV * (T_C + 273.15)
    th_sol = p['saturation']['x_max_shear'] / 2.0
    en = S.energetics(p, coverage_penalty='on')
    b_on = S.phase_boundaries(p, T_C, g, omega='on')
    b_off = S.phase_boundaries(p, T_C, g, omega='off')
    want = (en['eps'][2] + en['omega_bulk'] * th_sol
            + kt * math.log(th_sol / (1.0 - th_sol)))
    assert b_on['mu_cs_eV'] == pytest.approx(want, abs=1e-15)
    assert b_on['mu_cs_eV'] > b_off['mu_cs_eV']
    # the bulk at the ceiling is exactly at the declared solubility
    assert S.theta_of_mu(b_on['mu_cs_eV'], en['eps'][2], en['omega_bulk'],
                         kt) == pytest.approx(th_sol, abs=1e-12)


# ------------------------------- 8. mu is pinned through the surface riser

def test_mu_is_pinned_across_the_surface_coexistence(p):
    """A first-order transition holds mu fixed while the lever arm moves.

    Between the (1x2) onset and its completion the potential must not move
    at all, and the reconstructed area fraction must sweep 0 to 1."""
    for n in (1, 3):
        for cp in ('off', 'on'):
            g = geom(p, resolved_layers=n)
            b = S.phase_boundaries(p, T_C, g, omega=cp)
            lo, hi = b['VO_onset_umol_g'], b['VO_recon_complete_umol_g']
            assert hi > lo
            phis = []
            for f in (0.001, 0.25, 0.5, 0.75, 0.999):
                r = S.distribute(p, T_C, lo + f * (hi - lo), g, omega=cp)
                assert r['regime'] == 'surface_coexistence', (n, cp, f)
                assert r['mu_V_eV'] == pytest.approx(b['mu_t_eV'], abs=1e-12)
                phis.append(r['surface_phase']['phi_reconstructed'])
            assert all(b > a for a, b in zip(phis, phis[1:])), phis
            assert phis[0] < 0.01 and phis[-1] > 0.99


# ------------------------------------------------------ the declared state

def test_the_dataset_ships_the_new_machinery_switched_off(p):
    """Nothing above changes the workspace until someone asks it to."""
    d = p['defaults']
    assert d['resolved_layers'] == 1
    assert d['coverage_penalty'] == 'off'
    en = S.energetics(p)
    assert en['omega_layers'] == [0.0]
    assert en['omega_surface'] == 0.0 and en['omega_bulk'] == 0.0
    assert en['eps_layers'] == [en['eps'][1]]


def test_the_surface_penalty_is_held_at_zero_by_the_dataset(p):
    """Stage 1's design constraint, written into the data rather than the
    code: with omega_surface pinned at zero the phase ladder cannot move,
    so a change in a deeper layer has one possible cause."""
    om = p['vacancy_energetics']['coverage_penalty']['omega_eV']
    assert om['surface'] == 0.0
    g = geom(p)
    on = S.phase_boundaries(p, T_C, g, omega='on')
    off = S.phase_boundaries(p, T_C, g, omega='off')
    assert on['mu_t_eV'] == off['mu_t_eV']
    assert on['theta_transition'] == off['theta_transition']


def test_a_negative_penalty_is_refused(p):
    """omega < -4kT makes mu(theta) non-monotone and splits a class into
    two phases. That construction is not carried here, so the door is shut
    rather than left to produce a silently wrong root."""
    q = json.loads(json.dumps(p))
    q['vacancy_energetics']['coverage_penalty']['omega_eV']['layer1'] = -0.2
    with pytest.raises(ValueError, match='repulsive'):
        S.energetics(q, coverage_penalty='on')


def test_the_layer_profile_is_a_shape_not_a_table(p):
    """The stored span carries over to a different functional.

    eps_layers is the preset's own subsurface energy plus a fraction of its
    own subsurface-to-bulk gap, so the first layer is always the preset
    value exactly and the GGA column is not silently given sX numbers."""
    for preset in ('sx', 'gga'):
        en = S.energetics(p, preset=preset, resolved_layers=3)
        eps = en['eps']
        assert en['eps_layers'][0] == eps[1]
        assert eps[1] < en['eps_layers'][1] < en['eps_layers'][2] < eps[2]
    sx = S.energetics(p, preset='sx', resolved_layers=3)['eps_layers']
    illus = p['vacancy_energetics']['layer_profile']['eps_eV_sx_illustrative']
    for got, stored in zip(sx, illus):
        assert got == pytest.approx(stored, abs=5e-3)


def test_resolving_the_shell_finds_the_depth_structure_it_hid(p):
    """What the feature is for, stated as an outcome rather than a number.

    One class at one energy has to report a single occupancy for the whole
    shell. Resolved, the first layer stays crowded while the third is
    nearly empty - a profile the one-class model cannot represent, and the
    reason the bulk remainder moves."""
    one = S.distribute(p, T_C, VO, geom(p, ss_layers=3, resolved_layers=1))
    three = S.distribute(p, T_C, VO, geom(p, resolved_layers=3))
    assert len(three['layers']) == 3
    th = [l['theta'] for l in three['layers']]
    assert th[0] > th[1] > th[2], th
    assert th[0] > one['theta']['subsurface'] > th[2]
    # and it costs more: the one-class model charged all three slabs the
    # first subsurface energy, so it was holding the inventory too cheaply
    assert three['mu_V_eV'] > one['mu_V_eV']
    assert three['umol_g']['bulk'] > one['umol_g']['bulk']


# --------------------------------------------------- the 50-digit oracle

def test_the_implicit_root_matches_the_oracle(p):
    """Two algorithms on the same equation, not one checked against itself.

    scripts/oracle_statmech.py hands the residual to mpmath's root finder
    at 50 digits; the engine bisects it in double precision. Agreement to
    1e-13 is the two methods landing on the same number."""
    rows = REF['theta_of_mu']
    assert len(rows) == 180, len(rows)
    worst = 0.0
    for r in rows:
        kt = S.KB_EV * (r['T_C'] + 273.15)
        got = S.theta_of_mu(r['mu_eV'], r['eps_eV'], r['omega_eV'], kt)
        worst = max(worst, abs(got - r['theta']))
    assert worst < 1e-13, worst


def test_the_layered_ladder_matches_the_oracle(p):
    """Every resolution and interaction, through the whole construction.

    The oracle re-derives the phase boundaries, the branch selection and
    the layer occupancies independently, so this covers the parts of the
    construction the closed forms above cannot reach."""
    rows = REF['layer_ladder']
    assert len(rows) == 30, len(rows)
    for r in rows:
        g = geom(p, resolved_layers=r['resolved_layers'])
        d = S.distribute(p, r['T_C'], r['VO_total_umol_g'], g,
                         preset=r['preset'], omega=r['coverage_penalty'])
        tag = (r['resolved_layers'], r['coverage_penalty'],
               r['VO_total_umol_g'])
        assert d['regime'] == r['regime'], tag
        assert d['mu_V_eV'] == pytest.approx(r['mu_V_eV'], abs=1e-11), tag
        for i, cap in enumerate(r['N_layers']):
            assert g['N_layers'][i] == pytest.approx(cap, rel=1e-13), tag
        for i, th in enumerate(r['layer_theta']):
            assert d['layers'][i]['theta'] == pytest.approx(
                th, abs=1e-12), (tag, i)
            assert d['layers'][i]['umol_g'] == pytest.approx(
                r['layer_umol_g'][i], abs=1e-9), (tag, i)
        for k in ('surface', 'subsurface', 'bulk'):
            assert d['umol_g'][k] == pytest.approx(
                r['umol_g'][k], abs=1e-9), (tag, k)
        assert d['extended_defects_umol_g'] == pytest.approx(
            r['extended_defects_umol_g'], abs=1e-9), tag
        for k, v in r['boundaries'].items():
            assert d['phase_boundaries'][k] == pytest.approx(
                v, abs=1e-9), (tag, k)


def test_the_oracle_covers_all_four_phase_branches(p):
    """A ladder that never left the dilute branch would certify nothing."""
    seen = {r['regime'] for r in REF['layer_ladder']}
    assert seen == {'dilute', 'surface_coexistence', 'reconstructed',
                    'cs_coexistence'}, seen
