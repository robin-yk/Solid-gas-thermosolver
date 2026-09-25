"""Gates for the single model: inputs, site counting, Monte Carlo,
electrostatics, reconstruction bookkeeping and reproducible outputs."""
import csv
import math
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))

import run                                                # noqa: E402
from tofrange import aggregates as ag                     # noqa: E402
from tofrange import model as md                          # noqa: E402
from tofrange import particle as pt                       # noqa: E402
from tofrange.engine import COULOMB, Model                # noqa: E402

T = md.T_EQ
KT = pt.KB * T


def registry():
    with open(HERE / 'specification' / 'parameter_registry.csv') as f:
        return {r['id']: r['value'] for r in csv.DictReader(f)}


def floats(v):
    return tuple(float(x) for x in v.split(';'))


def cases():
    with open(HERE / 'outputs' / 'model_cases.csv') as f:
        return list(csv.DictReader(f))


def small(**kw):
    """The model on a coarse bulk grid (tests only)."""
    args = dict(cutoff=0.34, dG=0.0, f110=0.75, eps_r=64.0, n_bulk=60)
    args.update(kw)
    return md.particle(args.pop('energy_map', 'PAB'), args.pop('d_nm', 900.0), **args)


# ------------------------------------------------------------------ inputs
def test_constants_match_the_registry():
    r = registry()
    assert (pt.A_NM, pt.C_NM, pt.U) == (float(r['rutile_a_nm']), float(r['rutile_c_nm']), float(r['rutile_u']))
    assert (pt.RHO, pt.M) == (float(r['rutile_density_g_cm3']), float(r['rutile_molar_mass_g_mol']))
    assert (pt.NA, pt.KB) == (float(r['avogadro']), float(r['boltzmann_eV_K']))
    assert md.T_EQ == float(r['equilibrium_T_K'])
    assert md.ZHA == (float(r['zha2017_A_eV_A']), float(r['zha2017_B_eV_A2']), float(r['zha2017_C_eV_A6']))
    assert md.S0_PENALTY == float(r['s0_polaron_penalty_eV'])
    assert md.EPS == {'a_axis': float(r['permittivity_873K_a_axis']), 'c_axis': float(r['permittivity_873K_c_axis'])}
    assert md.MATSUNAGA == float(r['matsunaga_basal_minus_bridging_eV'])
    assert md.CUTOFFS == floats(r['mc_cutoff_nm'])
    assert md.DG == floats(r['recon_dG_scan_eV'])
    assert md.F110 == floats(r['facet_110_fraction'])
    lo, hi, n = floats(r['size_mixture_nm'])
    assert md.DIAMETERS == pytest.approx(tuple(np.linspace(lo, hi, int(n))))
    assert ag.BOX == tuple(int(x) for x in floats(r['mc_box_cells']))
    assert (run.MIN_SITES, run.MIN_FRACTION) == (float(r['denominator_absolute_threshold_umol_g']),
                                                 float(r['denominator_capacity_fraction_threshold']))
    assert run.STRONG_REDUCTION_C == float(r['strong_reduction_T_C'])


def test_every_o_and_ti_site_is_counted_once():
    assert pt.bridging_capacity(900) == pytest.approx(13.558, abs=5e-4)      # Note 2a: ~13.6
    for name in md.MAPS:
        m, lay = small(energy_map=name)
        o = ti = 0.0
        for i, t in enumerate(m.tags):
            if t == 'Ti':
                ti += m.C[i]
            elif i == lay.idx['cell']:
                o += 2 * m.C[i]; ti += 4 * m.C[i]
            elif t == 'bulk':
                o += m.C[i] * 4 * int(np.prod(ag.BOX))
            else:
                o += m.C[i]
        assert o == pytest.approx(pt.O_TOTAL, rel=1e-12) and ti == pytest.approx(pt.TI_TOTAL, rel=1e-12)


# ------------------------------------------------------ aggregates (MC)
def test_monte_carlo_is_exact_for_one_and_two_vacancies():
    beta = 1 / KT
    lat = ag.lattice((4, 4, 6), 0.6)
    n, J, st, idx = lat['nO'], lat['J'], lat['start'], lat['idx']
    q2 = sum(math.exp(-beta * J[p]) for i in range(n) for p in range(st[i], st[i + 1]) if idx[p] > i)
    q2 += n * (n - 1) / 2 - len(idx) / 2
    lnq = ag.ln_q(T, 0.6, (4, 4, 6), moves_per_n=2000)
    assert lnq[1] == pytest.approx(math.log(n), abs=1e-9)
    assert lnq[2] == pytest.approx(math.log(q2), abs=1e-9)


def test_monte_carlo_is_reproducible_across_seeds():
    # The last vacancy (N = 175, a perfect Ti matching) is not always placed.
    a, b = ag.ln_q(T, 1.0, seed=1), ag.ln_q(T, 1.0, seed=2)
    assert min(len(a), len(b)) >= 175
    assert -KT * a[170] / 170 == pytest.approx(-KT * b[170] / 170, abs=0.05)


# ---------------------------------------------------------- electrostatics
def test_capacitance_matrix_is_the_exact_inverse():
    rng = np.random.default_rng(1)
    r = np.sort(rng.uniform(1, 450, 12))[::-1]
    doms = [dict(C=1.0, E=[0, 0], v=[0, 1], region=k, shell=k) for k in range(12)]
    m = Model(doms, electrostatics=dict(radii=dict(enumerate(r)), eps=64, mass_g=1.62e-12))
    H = m.es['K'] / np.maximum.outer(m.es['r'], m.es['r'])
    q = rng.normal(size=12)
    assert m.potentials(q) == pytest.approx(H @ q, rel=1e-13)
    P = np.diag(m.es['main']) + np.diag(m.es['off'], 1) + np.diag(m.es['off'], -1)
    exact = KT * np.linalg.inv(H)
    assert np.abs(KT / m.es['K'] * P - exact).max() <= 1e-10 * np.abs(exact).max()


def test_split_domain_gradient_and_hessian_match_finite_differences():
    """A domain with its vacancy charge on one shell and its electrons on the
    next: the analytic dual derivatives equal central differences."""
    doms = [dict(C=1.0, E=[0.0, -0.3, 0.1, 0.2], v=[0, 1, 1, 0], e=[0, 0, 2, 1], g=[1, 2, 1, 2],
                 region=0, shell=0, shell2=1),
            dict(C=2.0, E=[0.0, 0.1], v=[0, 0], e=[0, 1], region=0, shell=1),
            dict(C=3.0, E=[0.0, 0.05], v=[0, 1], region=1, shell=2),
            dict(C=1.5, E=[0.0, 0.0], v=[0, 0], e=[0, 1], region=1, shell=2)]
    m = Model(doms, electrostatics=dict(radii={0: 10.0, 1: 9.8, 2: 5.0}, eps=64, mass_g=1e-15))
    beta = 1 / KT
    y, u = np.array([0.3, -0.2]), np.array([0.1, -0.05, 0.02])
    args = (beta, 1.0, 2.0, np.zeros(3))
    psi, _, g, Hyy, B, diag, off = m._global_eval(y, u, *args)
    H = np.zeros((5, 5)); H[:2, :2] = Hyy; H[2:, :2] = B; H[:2, 2:] = B.T
    H[2:, 2:] = np.diag(diag) + np.diag(off, 1) + np.diag(off, -1)
    z0, h = np.r_[y, u], 1e-5
    f = lambda z: m._global_eval(z[:2], z[2:], *args, need_hess=False)[0]        # noqa: E731
    fg = lambda z: m._global_eval(z[:2], z[2:], *args)[2]                         # noqa: E731
    for k in range(5):
        e = np.zeros(5); e[k] = h
        assert (f(z0 + e) - f(z0 - e)) / (2 * h) == pytest.approx(g[k], abs=1e-7)
        assert (fg(z0 + e) - fg(z0 - e)) / (2 * h) == pytest.approx(H[:, k], abs=1e-6)


def test_model_meets_inventory_neutrality_and_poisson():
    for name in ('PAB', 'LI_L2'):
        m, lay = small(energy_map=name)
        for N in (13.45, 1343.0):
            s = m.solve(N, T)
            assert abs(s.inventory_residual) < 1e-9 * N and abs(s.charge_residual) < 1e-9 * N
            assert s.poisson_residual < 1e-7 and abs(s.shell_charge.sum()) < 1e-8
            assert sum(md.populations(s, lay).values()) == pytest.approx(N, rel=1e-10)


def test_model_reaches_the_electron_reservoir_limit_as_eps_grows():
    """eps -> infinity removes the field: every state then follows the Gibbs
    form with one vacancy and one electron chemical potential."""
    m, lay = small(eps_r=1e14)          # the deviation falls as 1/eps (2.6e-6 at 1e12)
    for N in (13.45, 1343.0):
        s = m.solve(N, T)
        L = m.logg - (m.E - s.mu_eV * m.v - s.mu_e * m.e) / KT
        p = np.exp(L - L.max(axis=1)[:, None]); p /= p.sum(axis=1)[:, None]
        assert np.abs(s.x / m.C[:, None] - p).max() < 1e-7


def test_screening_length_is_debye_huckel():
    """A small surface perturbation decays over the Debye-Huckel length
    computed from the bulk site statistics."""
    R, eps = 450.0, 64.0
    edges = np.concatenate([[0.0], np.geomspace(1e-3, R, 2500)])
    doms, radii = [], {}
    for j in range(len(edges) - 1):
        C = pt.O_TOTAL * pt.shell_fraction(R, edges[j], edges[j + 1])
        doms.append(dict(C=C, E=[0, -0.05 if j == 0 else 0.0], v=[0, 1], region=j, shell=j))
        doms.append(dict(C=C / 2, E=[0, 0], v=[0, 0], e=[0, 1], region=j, shell=j))
        radii[j] = R - 0.5 * (edges[j] + edges[j + 1])
    m = Model(doms, electrostatics=dict(radii=radii, eps=eps, mass_g=pt.particle_mass(900)))
    s = m.solve(94.0, T)
    th = s.occupancy(4000)
    y = float(s.x[4001] @ m.e[4001]) / m.C[4001]
    n_o = 2 * pt.RHO / pt.M * pt.NA * 1e-21
    lam = 1 / math.sqrt(4 * math.pi * COULOMB * (n_o * 4 * th * (1 - th) + n_o / 2 * y * (1 - y)) / (eps * KT))
    z = R - m.es['r']
    sel = (z > 1.0) & (z < 3.0)
    fit = -1 / np.polyfit(z[sel], np.log(np.abs(s.phi[sel] - s.phi[-1])), 1)[0]
    assert fit == pytest.approx(lam, rel=3e-3)


# -------------------------------------------------------- reconstruction
def test_reconstruction_keeps_its_own_ti_and_vanishes_at_high_energy():
    m, lay = small(dG=5.0)
    assert md.surface(m.solve(94.0, T), lay)[2] < 1e-20
    m, lay = small(dG=-3.0)
    i, o = lay.idx['cell'], lay.cell_o
    free = m.e[i][:len(o)] - np.where(o == 3, 2, 0)
    assert free.max() == 4 and free[o == 3].max() == 2 and free.min() == 0
    assert not any(t == 'Ti' and m.region[k] == m.region[i] for k, t in enumerate(m.tags))
    s = m.solve(94.0, T)
    assert md.surface(s, lay)[2] > 0.5
    assert float(s.x[i] @ m.e[i]) <= 4 * m.C[i]


# ----------------------------------------------------------------- outputs
def test_q_is_the_ratio_to_r600_at_the_same_point():
    rows = cases()
    ref = {run.key(r): float(r['TOF_s_1']) for r in rows if r['sample'] == 'R600' and r['TOF_s_1'] != 'inf'}
    n = 0
    for r in rows:
        if r['Q_vs_R600']:
            assert float(r['Q_vs_R600']) == pytest.approx(float(r['TOF_s_1']) / ref[run.key(r)], rel=1e-5)
            n += 1
    assert n > 0


def test_pools_add_up_to_the_measured_inventory():
    inv = {s['sample']: float(s['inventory_umol_g']) for s in run.samples()}
    for r in cases():
        total = sum(float(r[f'N_{q}_umol_g']) for q in md.POOLS)
        assert total == pytest.approx(inv[r['sample']], rel=1e-5)


def test_saved_rows_regenerate():
    """Two parameter points solved again give the saved rows exactly."""
    saved = {(r['sample'],) + run.key(r): r for r in cases()}
    S = run.samples()
    for p in [('PAB', 0.34, 0.0, 1.0, 'a_axis'), ('LI_L2', 0.45, -0.4, 0.5, 'c_axis')]:
        for r in run.point((p, S)):
            old = saved[(r['sample'],) + run.key(r)]
            assert {k: old[k] for k in r} == r
