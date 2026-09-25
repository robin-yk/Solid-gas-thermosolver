"""Gates for the fixed-inventory TOF range: thermodynamics, electrostatics,
pairs, transport, classification and byte-stable outputs."""
import csv
import filecmp
import itertools
import math
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))

from tofrange import build as bd                          # noqa: E402
from tofrange import kinetics as kin                      # noqa: E402
from tofrange import particle as pt                       # noqa: E402
from tofrange import scenarios as sc                      # noqa: E402
from tofrange.engine import COULOMB, Model                # noqa: E402

T = 873.15
KT = pt.KB * T


def registry():
    with open(HERE / 'specification' / 'parameter_registry.csv') as f:
        return {r['id']: r['value'] for r in csv.DictReader(f)}


def floats(v):
    return tuple(float(x) for x in v.split(';'))


# ------------------------------------------------------------------ inputs
def test_constants_match_the_registry():
    r = registry()
    assert (pt.A_NM, pt.C_NM, pt.U) == (float(r['rutile_a_nm']), float(r['rutile_c_nm']), float(r['rutile_u']))
    assert (pt.RHO, pt.M) == (float(r['rutile_density_g_cm3']), float(r['rutile_molar_mass_g_mol']))
    assert (pt.NA, pt.KB, kin.H_PLANCK) == (float(r['avogadro']), float(r['boltzmann_eV_K']), float(r['planck_eV_s']))
    assert sc.T_EQ == float(r['equilibrium_T_K'])
    assert bd.A_LI == float(r['li_continuum_A_eV'])
    assert (0.25, 1.0) == floats(r['li_continuum_xi_sensitivity_nm'])
    assert bd.ZHA == (float(r['zha2017_A_eV_A']), float(r['zha2017_B_eV_A2']), float(r['zha2017_C_eV_A6']))
    assert bd.PAIR_RMAX == float(r['zha2017_range_nm'])
    assert bd.S0_PENALTY == float(r['s0_polaron_penalty_eV'])
    assert sc.EPS == {'a_axis': float(r['permittivity_873K_a_axis']), 'c_axis': float(r['permittivity_873K_c_axis'])}
    assert sc.THETA_C == float(r['reconstruction_onset_ML'])
    assert sc.STRONG_REDUCTION_C == float(r['strong_reduction_T_C'])
    assert sc.PAIR_BOUNDARY == float(r['pair_dilute_boundary'])
    assert (sc.MIN_SITES, sc.MIN_FRACTION) == (float(r['denominator_absolute_threshold_umol_g']),
                                               float(r['denominator_capacity_fraction_threshold']))
    import run
    assert run.T_OBS == floats(r['observation_time_grid_s'])
    diam = {s[2]['d_nm'] for s in sc.model_specs()}
    assert diam == {float(r['diameter_baseline_nm'])} | set(floats(r['diameter_sensitivity_nm']))
    lo, hi, n = floats(r['size_mixture_nm'])
    assert sc.SIZE_MIX == pytest.approx(tuple(np.linspace(lo, hi, int(n))))
    assert sc.RECON_F == floats(r['recon_fixed_fraction'])
    assert sc.RECON_DG == floats(r['recon_dG_scan_eV'])
    assert bd.RECON_ROW_VACANCIES == float(r['recon_row_vacancies_per_1x1'])
    from tofrange import aggregates as ag
    assert ag.BOX == tuple(int(x) for x in floats(r['mc_box_cells']))
    assert sc.MC_CUTOFFS == floats(r['mc_cutoff_nm'])
    assert sc.FACET_110 == floats(r['facet_110_fraction'])
    assert sc.MATSUNAGA == float(r['matsunaga_basal_minus_bridging_eV'])


def test_transport_edges_come_from_the_table():
    # Atomic O layer 5 is the in-plane O of trilayer 2; Jug keeps the lower
    # (downhill) barrier of each pair it reports in both directions.
    assert ((1, 'SBR'), (2, 'IPL'), 0.71) in kin.SURFACE_SETS['WU']
    assert dict(((a, b), B) for a, b, B in kin.SURFACE_SETS['JUG'])[((1, 'BRI'), (1, 'SBR'))] == pytest.approx(0.5804, abs=1e-4)
    assert kin.BULK_BARRIERS == {'IDDIR_2007': 1.10, 'UBERUAGA_2011': 1.5}


def test_capacities_partition_every_site_once():
    assert pt.bridging_capacity(900) == pytest.approx(13.558, abs=5e-4)      # Note 2a: ~13.6
    for name in sc.MAPS:
        m, lay = bd.build(name, 'LOCAL', 900.0)
        ti = sum(m.C[i] for i, t in enumerate(m.tags) if t == 'Ti')
        o = sum(r['C'] for r in lay.o)
        assert o == pytest.approx(pt.O_TOTAL, rel=1e-12) and ti == pytest.approx(pt.TI_TOTAL, rel=1e-12)


# ------------------------------------------------------------ equilibrium
def test_neutral_solution_matches_an_independent_bisection():
    m, lay = bd.build('PAB', 'NEUTRAL', 900.0)
    E, C = m.E[:, 1], m.C
    lo, hi = -5.0, 5.0
    for _ in range(200):
        mu = 0.5 * (lo + hi)
        lo, hi = (lo, mu) if np.sum(C / (1 + np.exp((E - mu) / KT))) > 94.0 else (mu, hi)
    sol = m.solve(94.0, T)
    assert sol.mu_eV == pytest.approx(mu, abs=1e-12)
    assert sol.x[:, 1] == pytest.approx(C / (1 + np.exp((E - mu) / KT)), rel=1e-9)


@pytest.mark.parametrize('name', ['PAB', 'HAM', 'LI_SBR1', 'LI_L2'])
def test_local_solution_satisfies_the_stationarity_condition(name):
    """eps + kT ln theta/(1-theta) + 2 kT ln y/(1-y) = mu at every O site."""
    m, lay = bd.build(name, 'LOCAL', 900.0)
    sol = m.solve(193.1, T)
    for o in lay.o:
        th = sol.occupancy(o['idx'])
        i_ti, c_ti = lay.ti[o['region']]
        y = float(sol.x[i_ti] @ m.e[i_ti]) / c_ti
        lhs = o['eps'] + KT * math.log(th / (1 - th)) + 2 * KT * math.log(y / (1 - y))
        assert lhs == pytest.approx(sol.mu_eV, abs=1e-9)
    assert abs(sol.inventory_residual) < 1e-9 and sol.charge_residual < 1e-8


def test_note_2b_is_reproduced():
    """Manuscript SI Note 2b: x(R), R600 interior and top-2-nm share, xi scan."""
    m, lay = bd.build('LI_CONT', 'LOCAL', 900.0)
    x = [bd.surface_coverage(m.solve(n, T), lay) for n in (13.45, 36.85, 94.0, 193.1, 781.2)]
    assert [round(100 * v, 2) for v in x[:3]] == [10.13, 17.59, 22.49]
    assert [round(100 * v, 1) for v in x[3:]] == [24.1, 24.9]

    def top2(xi):
        mm, ll = bd.build('LI_CONT', 'LOCAL', 900.0, xi=xi)
        s = mm.solve(94.0, T)
        return sum(s.occupancy(o['idx']) * o['C'] for k, o in enumerate(ll.o) if ll.edges[k + 1] <= 2.0) / 94.0
    sol = m.solve(94.0, T)
    assert round(sol.occupancy(lay.o[-1]['idx']), 5) == 0.00338
    assert (round(top2(0.25), 3), round(top2(0.5), 3), round(top2(1.0), 3)) == (0.064, 0.111, 0.193)
    ratio = [x[4] * (1 - x[4]) ** k / (x[2] * (1 - x[2]) ** k) for k in (2, 4, 8)]
    assert (round(min(ratio), 2), round(max(ratio), 2)) == (0.86, 1.04)


def _chain(n, periodic):
    """All 2^n patterns of an n-site bridging row with additive energy -1."""
    E, v, iso = [], [], []
    for s in itertools.product((0, 1), repeat=n):
        E.append(-1.0 * sum(s)); v.append(sum(s))
        nb = lambda i: [(i - 1) % n, (i + 1) % n] if periodic else [j for j in (i - 1, i + 1) if 0 <= j < n]  # noqa: E731
        iso.append(sum(s[i] and not any(s[j] for j in nb(i)) for i in range(n)))
    return dict(C=1.0, E=E, v=v), np.array(iso, float)


def test_isolated_count_is_exact_for_independent_sites():
    """A periodic 12-site row (4096 states) gives theta(1-theta)^2 isolated
    vacancies per site; the coupled_v1 4-site open chain gives
    (2-theta)/(2(1-theta)) times more."""
    for n, periodic in ((12, True), (4, False)):
        dom, iso = _chain(n, periodic)
        sol = Model([dom]).solve(0.3 * n, T)
        th = 0.3
        expect = th * (1 - th) ** 2 * (1 if periodic else (2 - th) / (2 * (1 - th)))
        assert float(sol.x[0] @ iso) / n == pytest.approx(expect, rel=1e-10)


def test_coupled_v1_neutral_r600_is_reproduced():
    """coupled_v1 inputs/project.json, Model.equilibrium(94.0, 873.15)."""
    cb = pt.bridging_capacity(900)
    dom, _ = _chain(4, False)
    dom = dict(dom, C=cb / 4, E=[1.31 * e for e in dom['E']], tag='surface')
    m = Model([dom, dict(C=cb, E=[0, -0.72], v=[0, 1], tag='sub'),
               dict(C=pt.O_TOTAL - 2 * cb, E=[0, 0], v=[0, 1], tag='bulk')])
    sol = m.solve(94.0, T)
    got = [sol.vacancies(t) for t in ('surface', 'sub', 'bulk')]
    assert got == pytest.approx([13.557923761144504, 13.215478361519827, 67.2265978773728], rel=1e-9)


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


def test_global_converges_to_poisson_and_neutrality():
    for name in ('PAB', 'LI_CONT'):
        m, lay = bd.build(name, 'GLOBAL', 900.0, eps_r=64.0)
        for N in (13.45, 781.2):
            s = m.solve(N, T)
            assert abs(s.inventory_residual) < 1e-9 * N and abs(s.charge_residual) < 1e-9 * N
            assert s.poisson_residual < 1e-8 and abs(s.shell_charge.sum()) < 1e-8


def test_global_tends_to_local_as_eps_vanishes():
    mg, _ = bd.build('LI_CONT', 'GLOBAL', 900.0, eps_r=1e-5, s0=0.0)
    ml, _ = bd.build('LI_CONT', 'LOCAL', 900.0)
    for N in (13.45, 94.0):
        assert np.abs(mg.solve(N, T).x - ml.solve(N, T).x).max() < 1e-5


def test_global_screening_length_is_debye_huckel():
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


# ------------------------------------------------------------------ pairs
@pytest.mark.parametrize('closure', ['NEUTRAL', 'LOCAL'])
def test_pairs_only_lower_the_surface_population(closure):
    for name in sc.MAPS:
        m0, l0 = bd.build(name, closure, 900.0)
        m1, l1 = bd.build(name, closure, 900.0, pairs=True)
        for N in (13.45, 193.1):
            assert bd.surface_coverage(m1.solve(N, T), l1) <= bd.surface_coverage(m0.solve(N, T), l0) + 1e-12


def test_pair_weights_are_the_mayer_second_virial():
    ps = bd.pair_states(0.0, T)
    for (r, n), E, g in zip(pt.oxygen_shells(1.0), ps['E'][1:], ps['g'][1:]):
        assert g * math.exp(-E / KT) == pytest.approx(0.5 * n * (math.exp(-bd.zha2017(r) / KT) - 1), rel=1e-12)


def test_frozen_pairs_at_the_same_temperature_reproduce_equilibrium():
    spec = dict(name='PAB', closure='LOCAL', d_nm=900.0)
    m1, l1 = bd.build(**spec, pairs=True)
    s1 = m1.solve(193.1, T)
    v_fix, q = 0.0, {}
    for i in l1.pair_idx:
        n = float(s1.x[i][1:].sum())
        v_fix += 2 * n
        r = m1.region_ids[m1.region[i]]
        q[r] = q.get(r, 0.0) + 4 * n
    m0, l0 = bd.build(**spec)
    s0 = m0.solve(193.1, T, fixed=dict(v=v_fix, q_region=q))
    for a, b in zip(l0.o, l1.o):
        assert s0.occupancy(a['idx']) == pytest.approx(s1.occupancy(b['idx']), rel=1e-8)


# -------------------------------------------------------------- transport
@pytest.mark.parametrize('name,closure', [('PAB', 'LOCAL'), ('LI_CONT', 'LOCAL'), ('HAM', 'NEUTRAL')])
def test_relaxation_ends_at_the_equilibrium_solver(name, closure):
    kw = dict(n_cont=600, z0_cont=1e-2) if name == 'LI_CONT' else {}
    m, lay = bd.build(name, closure, 900.0, **kw)
    r = kin.relax(m, lay, 781.2, 1273.15, T, 1.5, t_obs=(1.0, 1e5))
    assert np.abs(r['theta_end'] - r['theta_eq']).max() < 1e-9
    assert r['mass_drift'] < 1e-9


@pytest.mark.parametrize('closure', ['NEUTRAL', 'LOCAL'])
def test_transport_reproduces_sphere_diffusion(closure):
    """Flat energies: the slowest radial mode decays at D k1^2 (tan k1 R = k1 R).
    LOCAL adds the ambipolar factor 1 + 2(1 - theta)/(1 - y)."""
    B, R = 1.5, 450.0
    lam, zc = pt.bulk_hop()
    D = zc * lam ** 2 / 6 * KT / kin.H_PLANCK * math.exp(-B / KT)
    k1 = 4.493409457909064 / R
    m, lay = bd.build('LI_CONT', closure, 900.0, n_cont=800, z0_cont=1e-1, amp=0.0)
    o = sorted(lay.o, key=lambda r: r['z'])
    z = np.array([r['z'] for r in o]); C = np.array([r['C'] for r in o])
    th = 94.0 / pt.O_TOTAL * (1 + 0.02 * np.cos(np.pi * z / R))
    th *= 94.0 / (C @ th)
    factor = 1.0
    if closure == 'LOCAL':
        tb = 94.0 / pt.O_TOTAL
        factor = 1 + 2 * (1 - tb) / (1 - 4 * tb)
    tau = 1 / (factor * D * k1 ** 2)
    r = kin.relax(m, lay, 94.0, T, T, B, t_obs=(3 * tau, 6 * tau), theta0=th, rtol=1e-10)
    dev = lambda k: math.sqrt(C @ (r['theta'][:, k] - r['theta_eq']) ** 2)      # noqa: E731
    i1, i2 = np.searchsorted(r['t'], 3 * tau), np.searchsorted(r['t'], 6 * tau)
    rate = math.log(dev(i1) / dev(i2)) / (r['t'][i2] - r['t'][i1])
    assert rate == pytest.approx(factor * D * k1 ** 2, rel=2e-3)


def test_free_energy_falls_along_the_trajectory():
    m, lay = bd.build('PAB', 'LOCAL', 900.0)
    r = kin.relax(m, lay, 781.2, 1273.15, T, 1.5, t_obs=(1.0, 10.0))
    o = sorted(lay.o, key=lambda q: q['z'])
    eps = np.array([q['eps'] for q in o]); C = np.array([q['C'] for q in o])
    reg = np.array([q['region'] for q in o])
    h = lambda u: u * np.log(u) + (1 - u) * np.log1p(-u)                        # noqa: E731
    F = []
    for k in range(r['theta'].shape[1]):
        th = r['theta'][:, k]
        f = C @ (th * eps) + KT * (C @ h(th))
        for g, (i_ti, c_ti) in lay.ti.items():
            y = 2 * (C[reg == g] @ th[reg == g]) / c_ti
            f += KT * c_ti * h(y)
        F.append(f)
    assert np.all(np.diff(F) <= 1e-9 * abs(F[0]))


# ---------------------------------------------------------- classification
def cases():
    with open(HERE / 'outputs' / 'cases.csv') as f:
        return list(csv.DictReader(f))


def test_classification_and_bounds():
    rows = cases()
    by = lambda **k: [r for r in rows if all(r[a] == b for a, b in k.items())]           # noqa: E731
    core = by(sample='R600', family='PRIMARY')
    assert len(core) == 40 and all('RECON_CAP' in r['flags'] for r in core)
    # At the cap, the bridging definition gives the SI 2a denominator on this sphere.
    bri = by(sample='R600', family='PRIMARY', reactive='BRI')
    assert [float(r['N_react_umol_g']) for r in bri] == pytest.approx([0.17 * pt.bridging_capacity(900)] * len(bri), rel=1e-5)
    a600 = by(sample='A600', family='PRIMARY', energy_map='LI_CONT xi 0.5', closure='LOCAL', reactive='BRI')[0]
    assert a600['flags'] == '' and float(a600['theta_bri']) == pytest.approx(0.1013, abs=1e-4)
    assert all(r['cls'] in ('BOUNDARY', 'INVALID') for r in by(sample='R1000'))
    assert not by(sample='R1100')
    # Explicit reconstruction is never capped.
    assert not any('RECON_CAP' in r['flags'] for r in rows if r['family'].startswith('recon_'))


# ------------------------------------------------------ aggregates (MC)
def test_monte_carlo_is_exact_for_one_and_two_vacancies():
    from tofrange import aggregates as ag
    beta = 1 / KT
    lat = ag.lattice((4, 4, 6), 0.6)
    n, J, st, idx = lat['nO'], lat['J'], lat['start'], lat['idx']
    q2 = sum(math.exp(-beta * J[p]) for i in range(n) for p in range(st[i], st[i + 1]) if idx[p] > i)
    q2 += n * (n - 1) / 2 - len(idx) / 2
    lnq = ag.ln_q(T, 0.6, (4, 4, 6), moves_per_n=2000)
    assert lnq[1] == pytest.approx(math.log(n), abs=1e-9)
    assert lnq[2] == pytest.approx(math.log(q2), abs=1e-9)


def test_monte_carlo_is_reproducible_across_seeds():
    from tofrange import aggregates as ag
    # The last vacancy (N = 175, a perfect Ti matching) is not always placed.
    a, b = ag.ln_q(T, 1.0, seed=1), ag.ln_q(T, 1.0, seed=2)
    assert min(len(a), len(b)) >= 175
    assert -KT * a[170] / 170 == pytest.approx(-KT * b[170] / 170, abs=0.05)


def test_aggregates_empty_the_surface():
    """Pairwise-additive clusters bind each vacancy by eV, far below any
    surface site energy, so the bridging row empties."""
    m, lay = bd.build('PAB', 'LOCAL', 900.0, aggregates=1.0)
    s = m.solve(94.0, T)
    assert bd.surface_state(s, lay)[0] < 1e-3
    assert abs(s.inventory_residual) < 1e-8


# -------------------------------------------------------- reconstruction
def test_fixed_reconstruction_bookkeeping():
    m, lay = bd.build('PAB', 'LOCAL', 900.0, recon=('fixed', 0.5))
    s = m.solve(94.0, T, fixed=lay.fixed)
    free = sum(float(s.x[i] @ m.v[i]) for i in range(len(m.C)))
    assert free + lay.fixed['v'] == pytest.approx(94.0, rel=1e-12)
    assert lay.fixed['v'] == pytest.approx(0.5 * 0.5 * pt.bridging_capacity(900))
    assert bd.surface_state(s, lay)[1] == pytest.approx(0.5 * pt.bridging_capacity(900))


def test_reconstruction_state_is_monotone_in_its_energy():
    f = []
    for dG in (-0.4, -0.2, 0.0, 0.2, 0.4, 5.0):
        m, lay = bd.build('PAB', 'LOCAL', 900.0, recon=('state', dG))
        f.append(bd.surface_state(m.solve(94.0, T), lay)[2])
    assert all(a >= b for a, b in zip(f, f[1:])) and f[-1] < 1e-20
    m0, l0 = bd.build('PAB', 'LOCAL', 900.0)
    m5, l5 = bd.build('PAB', 'LOCAL', 900.0, recon=('state', 5.0))
    assert bd.surface_state(m5.solve(94.0, T), l5)[0] == pytest.approx(
        bd.surface_state(m0.solve(94.0, T), l0)[0], rel=1e-9)


def test_size_mixture_is_the_mass_average():
    rows = cases()
    mix = [r for r in rows if r['family'] == 'size_mix' and r['sample'] == 'R600'
           and r['energy_map'] == 'HAM' and r['closure'] == 'LOCAL' and r['reactive'] == 'BRI'][0]
    n = []
    for d in sc.SIZE_MIX:
        m, lay = bd.build('HAM', 'LOCAL', d)
        n.append(sc.site_counts(m.solve(94.0, T), lay, 'size_mix')[0]['BRI'][0])
    assert float(mix['N_react_umol_g']) == pytest.approx(np.mean(n), rel=1e-5)


def test_facet_share_scales_the_surface_only():
    m1, l1 = bd.build('PAB', 'LOCAL', 900.0)
    mf, lf = bd.build('PAB', 'LOCAL', 900.0, f110=1.0)
    assert np.array_equal(m1.C, mf.C)
    m5, l5 = bd.build('PAB', 'LOCAL', 900.0, f110=0.5)
    assert l5.c_bri == pytest.approx(0.5 * l1.c_bri)
    assert sum(o['C'] for o in l5.o) == pytest.approx(pt.O_TOTAL, rel=1e-12)


def test_q_is_the_ratio_to_r600_in_the_same_scenario():
    rows = cases()
    key = lambda r: (r['family'], r['variant'], r['diameter_nm'], r['energy_map'], r['closure'], r['reactive'])  # noqa: E731
    ref = {key(r): float(r['TOF_s_1']) for r in rows if r['sample'] == 'R600' and r['cls'] != 'INVALID'
           and r['variant'] != 'R600 94.6'}
    for r in rows:
        if r['Q_vs_R600']:
            assert float(r['Q_vs_R600']) == pytest.approx(float(r['TOF_s_1']) / ref[key(r)], rel=1e-5)


def test_outputs_regenerate_byte_for_byte(tmp_path):
    work = tmp_path / 'p'
    shutil.copytree(HERE, work, ignore=shutil.ignore_patterns('outputs', '__pycache__', '.pytest_cache'))
    (work / 'outputs').mkdir()
    subprocess.run([sys.executable, 'run.py'], cwd=work, check=True)
    for name in ('cases.csv', 'sample_tof_range.csv', 'sensitivity_effects.csv', 'transport_gate.csv'):
        assert filecmp.cmp(work / 'outputs' / name, HERE / 'outputs' / name, shallow=False), name
