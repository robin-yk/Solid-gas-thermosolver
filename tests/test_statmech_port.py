"""The browser stat-mech engine against the Python reference, move by move.

statmech.js is a port, not a reimplementation: with the same seed the swap
trajectory is identical, so final vacancy positions, class occupancies and
histograms are compared for equality, not tolerance. The only library call
whose last ulp can differ between the runtimes is exp/log, which is compared
against (never accumulated) in the Metropolis step and enters results only
through the Widom average - hence the 1e-9 tolerance on mu alone.

Release gates: a missing interpreter or artefact fails, never skips.
"""

import json
import re
import pathlib
import shutil
import subprocess

import pytest

from solidgas import statmech as S

ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA = ROOT / 'data'
WEB = ROOT / 'web'
DOCS = ROOT / 'docs'
HARNESS = ROOT / 'tests' / 'js' / 'statmech_harness.js'
REF = json.loads((DATA / 'reference_statmech.json').read_text())
P = S.load_params()

NO_PAIRS = {'pair_eV': {'along_row': [], 'cross_row': []}}
PAIRS = {'pair_eV': P['surface_ordering']['pair_eV']}
ROW_PAIRS = {'pair_eV': {
    'along_row': P['surface_ordering']['pair_eV']['along_row'],
    'cross_row': []}}

MC_CASES = [
    {'name': 'ideal_slab', 'rows': 4, 'row_sites': 12, 'ss_layers': 1,
     'bulk_layers': 2, 'ordering': NO_PAIRS,
     'eps': REF['ideal_ensemble']['eps_eV'], 'T_C': 600, 'n_vac': 14,
     'sweeps_eq': 400, 'sweeps_sample': 1200, 'widom_every': 5,
     'hist_every': 20, 'blocks': 12, 'seed': 902},
    {'name': 'interacting_ring', 'rows': 1, 'row_sites': 10, 'ss_layers': 0,
     'bulk_layers': 0, 'ordering': ROW_PAIRS, 'eps': [0.0, 0.0, 0.0],
     'T_C': 600, 'n_vac': 3, 'sweeps_eq': 1000, 'sweeps_sample': 8000,
     'widom_every': 5, 'hist_every': 2, 'blocks': 16, 'seed': 903},
    {'name': 'interacting_slab', 'rows': 4, 'row_sites': 25, 'ss_layers': 1,
     'bulk_layers': 2, 'ordering': PAIRS,
     'eps': REF['ideal_ensemble']['eps_eV'], 'T_C': 600, 'n_vac': 20,
     'sweeps_eq': 500, 'sweeps_sample': 1500, 'widom_every': 10,
     'hist_every': 10, 'blocks': 15, 'seed': 904},
]

MC_SMALL = {'rows': 4, 'row_sites': 12, 'ss_layers': 1, 'bulk_layers': 2,
            'fillings': [2, 4, 8, 14, 22, 34], 'equil_sweeps': 600,
            'warm_equil_sweeps': 300, 'sample_sweeps': 1600,
            'widom_every': 5, 'hist_every': 20, 'blocks': 16, 'seed': 12345}

PIPELINES = [
    {'name': 'small_default', 'T_C': 600.0, 'VO_total': 95.0,
     'zero_pairs': False, 'mc': MC_SMALL},
    {'name': 'small_j0', 'T_C': 600.0, 'VO_total': 95.0,
     'zero_pairs': True, 'mc': MC_SMALL},
    {'name': 'shipped_default', 'T_C': 600.0, 'VO_total': 95.0,
     'zero_pairs': False, 'mc': None},
]


def _py_mc(c):
    lat = S.build_lattice(c['rows'], c['row_sites'], c['ss_layers'],
                          c['bulk_layers'], c['ordering'])
    rng = S.make_rng(c['seed'])
    occ = bytearray(lat['n'])
    vac = []
    S.adjust_filling(lat, occ, vac, c['n_vac'], rng)
    kt = S.KB_EV * (c['T_C'] + 273.15)
    return S.mc_run(lat, c['eps'], kt, occ, vac, c['sweeps_eq'],
                    c['sweeps_sample'], rng, widom_every=c['widom_every'],
                    hist_every=c['hist_every'], blocks=c['blocks'])


SWEEP_VO = [10.0, 60.0, 95.0, 150.0]

_HP_REF = json.loads(
    (DATA / 'reference_results_high_precision.json').read_text())
_RATIO = float(_HP_REF['boundary_validation']['h2_h2o_boundary']
               ['H2O_over_H2_critical_ratio'])

VALIDITY_CASES = [
    {'name': 'stable_rwgs', 'feed': {'CO2': 1, 'H2': 1},
     'T_K': 873.15, 'mass_g': 0.1},
    {'name': 'invalid_hot_h2', 'feed': {'H2': 1},
     'T_K': 1473.15, 'mass_g': 0.002},
    {'name': 'boundary_h2_h2o', 'feed': {'H2': 1, 'H2O': _RATIO * (1 - 1e-3)},
     'T_K': 1173.15, 'mass_g': 0.1},
]


CGMC_BARRIERS = [0.65, 1.9, 2.1]
CGMC_NO_FILL = {'A_eff_per_s_atm': {'surface': 0.0}}
CGMC_ALL_BULK = {'surface': 0.0, 'subsurface': 0.0, 'bulk': 1.0}
CGMC_STOCH_CG = {'E_m_bulk_eV': 1.9, 'n_cells': 8,
                 'stochastic_column_counts': 50000}
LADDER_VO = [1.0, 3.5, 9.0, 95.0, 200.0]
SHELL_LAYERS = [1, 2, 3, 4]


def _cgmc_linear_sys():
    s = REF['cgmc_linear']['spec']
    kt = S.KB_EV * (s['T_C'] + 273.15)
    return {'m': 3, 'q': [float(x) for x in s['q']], 'eps_c': s['eps'],
            'lam': s['lam'], 'beta': 1.0 / kt, 'k_fill': s['k_fill'],
            'eta0': s['eta0'], 'V0': sum(s['eta0']), 'R_nm': 3.0,
            'a_nm': 1.0, 'centers': [2.5, 1.5, 0.5],
            'T_reox_C': s['T_C'], 'E_m_eV': 0.0, 'eps_ctrl': 0.02,
            'column_counts': 1e5}


@pytest.fixture(scope='module')
def pair(tmp_path_factory):
    assert shutil.which('node') is not None, \
        'node is required to run the release gate'
    spec = tmp_path_factory.mktemp('spec') / 'statmech_spec.json'
    spec.write_text(json.dumps({'mc_cases': MC_CASES,
                                'pipelines': PIPELINES,
                                'sweep_vo': SWEEP_VO,
                                'ladder_vo': LADDER_VO,
                                'shell_layers': SHELL_LAYERS,
                                'validity_cases': VALIDITY_CASES,
                                'cgmc': {
                                    'sys': _cgmc_linear_sys(),
                                    'times': [r['t_s'] for r in
                                              REF['cgmc_linear']['rows']],
                                    'barriers': CGMC_BARRIERS,
                                    'no_fill': CGMC_NO_FILL,
                                    'all_bulk': CGMC_ALL_BULK,
                                    'stoch_cg': CGMC_STOCH_CG}}))
    r = subprocess.run(['node', str(HARNESS), str(ROOT), str(spec)],
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, f'harness failed:\n{r.stderr[-2000:]}'
    js = json.loads(r.stdout)
    py = {'mc': {c['name']: _py_mc(c) for c in MC_CASES},
          'pipelines': {}}
    for c in PIPELINES:
        if c['mc'] is None:
            continue                       # shipped size: JS only, see below
        out = S.run_full(P, t_c=c['T_C'], vo_total=c['VO_total'],
                         zero_pairs=c['zero_pairs'], mc=c['mc'])
        eps = S.energetics(P)['eps']
        out['sweep'] = S.sweep(P, c['T_C'], out['geometry'], SWEEP_VO,
                               eps=eps)
        py['pipelines'][c['name']] = out
    return py, js


def test_rng_streams_are_identical(pair):
    _, js = pair
    rng = S.make_rng(12345)
    assert js['rng'] == [rng() for _ in range(8)]
    assert js['rng'] == REF['rng_sfc32_seed12345_first8']


def test_swap_trajectories_end_in_the_same_state(pair):
    """Final vacancy positions equal as integers: the whole accept/reject
    history matched, not just the averages."""
    py, js = pair
    for row in js['mc_cases']:
        assert row['vac_sorted'] == py['mc'][row['name']]['vac_sorted'], \
            row['name']


def test_occupancies_and_histograms_are_bitwise_equal(pair):
    py, js = pair
    for row in js['mc_cases']:
        ref = py['mc'][row['name']]
        for k in ('theta_s', 'theta_ss', 'theta_b', 'E_mean_eV'):
            assert row[k] == ref[k], (row['name'], k)
        assert row['hist'] == ref['hist'], row['name']
        for k in ('err_s', 'err_ss', 'err_b'):
            assert row[k] == pytest.approx(ref[k], abs=1e-15), (row['name'], k)
        assert row['E_drift_eV'] < 1e-9


def test_widom_mu_agrees_to_library_ulp(pair):
    py, js = pair
    for row in js['mc_cases']:
        assert row['mu_eV'] == pytest.approx(
            py['mc'][row['name']]['mu_eV'], abs=1e-9), row['name']


def test_pipelines_agree_end_to_end(pair):
    py, js = pair
    for row in js['pipelines']:
        if row['name'] not in py['pipelines']:
            continue
        ref = py['pipelines'][row['name']]
        assert row['mu_V_eV'] == pytest.approx(ref['mu_V_eV'], abs=1e-12)
        for k in ('surface', 'subsurface', 'bulk', 'extended'):
            assert row['fractions'][k] == pytest.approx(
                ref['fractions'][k], rel=1e-12, abs=1e-15), k
        assert row['regime'] == ref['regime']
        assert row['surface_phase']['phase'] == \
            ref['surface_phase']['phase']
        assert row['surface_phase']['phi_reconstructed'] == pytest.approx(
            ref['surface_phase']['phi_reconstructed'], abs=1e-12)
        assert row['extended_defects_umol_g'] == pytest.approx(
            ref['extended_defects_umol_g'], rel=1e-12, abs=1e-12)
        for k, v in ref['phase_boundaries'].items():
            assert row['phase_boundaries'][k] == pytest.approx(
                v, rel=1e-12, abs=1e-15), k
        assert row['surface_reconstruction_regime'] == \
            ref['surface_reconstruction_regime']
        assert row['warn_kinds'] == [w['kind'] for w in ref['warnings']]
        for a, b in zip(row['isotherm'], ref['isotherm']):
            assert a['filling'] == b['filling']
            for k in ('theta_s', 'theta_ss', 'theta_b', 'E_mean_eV'):
                assert a[k] == b[k], k
            assert a['hist'] == b['hist']
            assert a['mu_eV'] == pytest.approx(b['mu_eV'], abs=1e-9)
            # the sampled surface configuration the particle view draws
            assert a['surf_vac'] == b['surf_vac']
            assert a['rows'] == b['rows']
            assert a['row_sites'] == b['row_sites']


def test_shipped_configuration_reproduces_the_oracle_flagship(pair):
    """The shipped default (95 umol/g at 600 C), exactly as the page runs
    it: the reconstructed regime, with the 50-digit oracle certifying the
    partition, the phase, and the whole boundary ladder."""
    _, js = pair
    ref = REF['flagship']
    fl = js['analytic_flagship']
    assert fl['mu_V_eV'] == pytest.approx(ref['mu_V_eV'], abs=1e-9)
    for k in ('surface', 'subsurface', 'bulk', 'extended'):
        assert fl['fractions'][k] == pytest.approx(
            ref['fractions'][k], rel=1e-9, abs=1e-12), k
    for k in ('surface', 'subsurface', 'bulk'):
        assert fl['umol_g'][k] == pytest.approx(
            ref['umol_g'][k], rel=1e-9), k
    assert fl['regime'] == ref['regime'] == 'reconstructed'
    assert fl['surface_phase']['phase'] == '1x2'
    for k, v in ref['boundaries'].items():
        assert fl['phase_boundaries'][k] == pytest.approx(
            v, rel=1e-9, abs=1e-12), k
    assert fl['warn_kinds'] == ref['warn_kinds']
    # the full shipped pipeline (with MC attached) reaches the same state
    row = next(r for r in js['pipelines'] if r['name'] == 'shipped_default')
    assert row['mu_V_eV'] == pytest.approx(ref['mu_V_eV'], abs=1e-9)
    assert row['regime'] == 'reconstructed'
    assert sum(row['umol_g'].values()) == pytest.approx(95.0, rel=1e-12)
    # surface ordering is a separate exhibit, sampled at the coverage of
    # the surviving (1x1) patches (theta_t), not at the line phase
    sp = row['spacing']
    assert sp is not None
    assert sp['at_theta_s'] <= P['surface_phases']['theta_transition'] + 0.05
    assert sp['modal_gap_sites'] in (3, 4, 5)
    assert sp['P_gap_le_2'] < 0.15


def test_ladder_parity_walks_the_same_regimes(pair):
    """The JS engine walks the four regimes with the same numbers."""
    _, js = pair
    g = S.geometry(P)
    for got in js['ladder']:
        d = S.distribute(P, 600.0, got['vo'], g)
        assert got['regime'] == d['regime'], got['vo']
        assert got['mu_V_eV'] == pytest.approx(d['mu_V_eV'], abs=1e-12)
        for k in ('surface', 'subsurface', 'bulk'):
            assert got['umol_g'][k] == pytest.approx(
                d['umol_g'][k], rel=1e-12, abs=1e-12), (got['vo'], k)
        assert got['extended'] == pytest.approx(
            d['extended_defects_umol_g'], rel=1e-12, abs=1e-12)
        assert got['phi'] == pytest.approx(
            d['surface_phase']['phi_reconstructed'], abs=1e-12)
        assert got['warn_kinds'] == [w['kind'] for w in d['warnings']]
    assert [r['regime'] for r in js['ladder']] == [
        'dilute', 'surface_coexistence', 'reconstructed', 'reconstructed',
        'cs_coexistence']


def test_shell_thickness_parity_and_oracle_ladder(pair):
    """The shell control is a real parameter in both engines.

    The declared subsurface thickness moves the shell/bulk split, so the
    JS engine has to reproduce the Python numbers and both have to land on
    the oracle ladder - otherwise the browser control would drift from the
    certified table in docs/defect-model.md."""
    _, js = pair
    assert [r['ss_layers'] for r in js['shell']] == SHELL_LAYERS
    ref_by_n = {r['ss_layers']: r for r in REF['shell_ladder']}
    for got in js['shell']:
        n = got['ss_layers']
        g = S.geometry(P, ss_layers=n)
        for k in ('N_s', 'N_ss', 'N_b', 'N_O_total', 'area_m2_g',
                  'layer_nm', 'shell_nm'):
            assert got['geometry'][k] == pytest.approx(g[k], rel=1e-12), k
        d = S.distribute(P, 600.0, 95.0, g)
        assert got['regime'] == d['regime'] == ref_by_n[n]['regime']
        assert got['mu_V_eV'] == pytest.approx(d['mu_V_eV'], abs=1e-12)
        assert d['mu_V_eV'] == pytest.approx(
            ref_by_n[n]['mu_V_eV'], abs=1e-9)
        for k in ('surface', 'subsurface', 'bulk'):
            assert got['umol_g'][k] == pytest.approx(
                d['umol_g'][k], rel=1e-12, abs=1e-12), (n, k)
            assert d['umol_g'][k] == pytest.approx(
                ref_by_n[n]['umol_g'][k], rel=1e-9, abs=1e-9), (n, k)
        assert got['warn_kinds'] == [w['kind'] for w in d['warnings']]
    # the browser default is the shipped default
    assert P['defaults']['subsurface_layers'] == 1


def test_bet_sensitivity_case_matches_the_oracle():
    """The measured-BET geometry (1.60 m2/g) against the diameter model:
    the boundaries move by ~2%, the verdict does not."""
    gd = S.geometry(P)
    gb = S.geometry(P, bet_m2_g=REF['geometry_bet_1p60']['area_m2_g'])
    for k in ('N_s', 'N_ss', 'N_b', 'N_O_total'):
        assert gb[k] == pytest.approx(REF['geometry_bet_1p60'][k], rel=1e-9)
    ref = REF['bet_case']
    d = S.distribute(P, ref['T_C'], ref['VO_total_umol_g'], gb)
    assert d['mu_V_eV'] == pytest.approx(ref['mu_V_eV'], abs=1e-9)
    assert d['regime'] == ref['regime'] == 'reconstructed'
    vd = S.distribute(P, ref['T_C'], ref['VO_total_umol_g'], gd)
    for k in ('VO_onset_umol_g', 'VO_recon_complete_umol_g'):
        shift = abs(d['phase_boundaries'][k] - vd['phase_boundaries'][k])
        assert shift / vd['phase_boundaries'][k] < 0.05, k


def test_accessibility_and_sweep_parity(pair):
    """The CO2 layer and the sweep rows, port vs reference."""
    py, js = pair
    for row in js['pipelines']:
        if row['name'] not in py['pipelines']:
            continue
        ref = py['pipelines'][row['name']]
        a, b = row['accessibility'], ref['accessibility']
        for c in ('surface', 'subsurface', 'bulk'):
            assert a['P'][c] == pytest.approx(b['P'][c], rel=1e-12,
                                              abs=1e-15), c
        assert a['accessible_umol_g'] == pytest.approx(
            b['accessible_umol_g'], rel=1e-12)
        assert a['f_recoverable'] == pytest.approx(
            b['f_recoverable'], rel=1e-12)
        for ra, rb in zip(row['sweep'], ref['sweep']):
            assert ra['VO_total_umol_g'] == rb['VO_total_umol_g']
            assert ra['regime'] == rb['regime']
            for k in ('surface', 'subsurface', 'bulk', 'extended',
                      'accessible', 'f_recoverable'):
                assert ra[k] == pytest.approx(rb[k], rel=1e-12,
                                              abs=1e-15), k


def test_shipped_accessibility_matches_the_oracle(pair):
    """The placeholder kinetics on the phase-aware flagship split: the
    300 s oracle row certifies the JS chain (illustrative until the
    recovery experiments land, and the page says so)."""
    _, js = pair
    ref = next(r for r in REF['co2_flagship'] if r['exposure_s'] == 300)
    a = js['analytic_flagship']['access']
    for c in ('surface', 'subsurface', 'bulk'):
        assert a['P'][c] == pytest.approx(ref['P'][c], rel=1e-9,
                                          abs=1e-12), c
    assert a['accessible_umol_g'] == pytest.approx(
        ref['accessible_umol_g'], rel=1e-9)
    assert a['f_recoverable'] == pytest.approx(
        ref['f_recoverable'], rel=1e-9)
    # the 95 row of the shipped sweep is the flagship distribution
    row = next(r for r in js['pipelines'] if r['name'] == 'shipped_default')
    r95 = next(r for r in row['sweep'] if r['VO_total_umol_g'] == 95.0)
    for k in ('surface', 'subsurface', 'bulk'):
        assert r95[k] == pytest.approx(REF['flagship']['umol_g'][k],
                                       rel=1e-9), k
    assert r95['regime'] == 'reconstructed'


def test_phase_validity_parity_on_real_equilibria(pair):
    """The bridge's verdict, computed by the JS thermo engine + JS tiering
    against the Python pair - the exact chain the oxygen-environment card
    runs in the browser."""
    from solidgas import activeset as A
    _, js = pair
    is_red = {s: A.STOICH[s][1] / A.STOICH[s][0] < 2.0 - 1e-12
              for s in A.SOLIDS}
    by_name = {r['name']: r for r in js['validity']}
    for c in VALIDITY_CASES:
        r = A.solve(c['feed'], c['T_K'], mass_g=c['mass_g'])
        v = S.phase_validity(r['active_condensed_phases'],
                             r['inactive_phase_reduced_costs'], is_red)
        jrow = by_name[c['name']]
        assert jrow['tier'] == v['tier'], c['name']
        assert jrow['active'] == r['active_condensed_phases'], c['name']
        assert jrow['nearest_phase'] == v['nearest_phase'], c['name']
        if v['margin_per_mol_O_kJ'] is None:
            assert jrow['margin_per_mol_O_kJ'] is None
        else:
            assert jrow['margin_per_mol_O_kJ'] == pytest.approx(
                v['margin_per_mol_O_kJ'], abs=1e-9), c['name']
        assert jrow['mu_O_kJ_mol'] == pytest.approx(
            r['lambda_kJ_per_mol']['O'], abs=1e-9), c['name']
    assert by_name['stable_rwgs']['tier'] == 'stable'
    assert by_name['invalid_hot_h2']['tier'] == 'invalid'
    assert by_name['boundary_h2_h2o']['tier'] == 'boundary'


def test_cgmc_port_matches_the_reference_engine(pair):
    """Deterministic transport runs, move for move: same recovery, same
    final occupancies, same number of controller steps."""
    from solidgas import cgmc as C
    _, js = pair
    for em, row in zip(CGMC_BARRIERS, js['cgmc']['ladder']):
        r = C.run(C.build(P, cg={'E_m_bulk_eV': em}), 300.0)
        assert row['recovery'] == pytest.approx(
            r['recovered_fraction'], rel=1e-12, abs=1e-15), em
        assert row['steps'] == r['steps'], em
        assert abs(row['mass']) < 1e-12
    sysc = C.build(P, kin=CGMC_NO_FILL, dist_fractions=CGMC_ALL_BULK)
    rc = C.run(sysc, 300.0)
    assert js['cgmc']['closed']['steps'] == rc['steps']
    for a, b in zip(js['cgmc']['closed']['final_eta'], rc['final_eta']):
        assert a == pytest.approx(b, rel=1e-12, abs=1e-9)


def test_cgmc_linear_port_matches_the_expm_oracle(pair):
    _, js = pair
    for got, row in zip(js['cgmc']['linear'], REF['cgmc_linear']['rows']):
        for a, b in zip(got, row['eta']):
            assert a == pytest.approx(b, rel=6e-3, abs=1e-6)


def test_cgmc_stochastic_same_seed_parity(pair):
    """Same seed, same event counts: a single flipped Poisson draw would
    shift a cell by a whole count (~1e-4 relative), far above the float
    noise this tolerance admits."""
    from solidgas import cgmc as C
    _, js = pair
    sys0 = C.build(P, kin=CGMC_NO_FILL, cg=CGMC_STOCH_CG)
    r = C.run(sys0, 60.0, mode='stoch', seed=902)
    for a, b in zip(js['cgmc']['stoch'], r['final_eta']):
        assert a == pytest.approx(b, rel=1e-9, abs=1e-9)


def test_dielectric_parity(pair):
    _, js = pair
    for v, row in zip([95.0, 10.0, 1000.0, 0.0], js['dielectric']):
        ref = S.dielectric(P, v)
        if ref is None:
            assert row is None
            continue
        assert row['eps2'] == pytest.approx(ref['eps2'], rel=1e-12)
        assert row['extrapolated'] == ref['extrapolated']


def test_the_layer_resolved_partition_matches_the_python_engine(pair):
    """Three-way parity for Stage 1: both engines, every branch.

    Resolution and interaction are both off in the shipped dataset, so
    these combinations would otherwise never be exercised across the port
    boundary. mu is compared to 1e-11 rather than exactly because the
    interacting root is bisected against a residual containing exp(), and
    the two runtimes may differ in its last ulp; the occupancies that
    follow from it are compared just as tightly."""
    _, js = pair
    rows = js['layers']
    assert len(rows) == 30, len(rows)
    for r in rows:
        g = S.geometry(P, resolved_layers=r['n'])
        d = S.distribute(P, 600.0, r['vo'], g, omega=r['omega'])
        tag = (r['n'], r['omega'], r['vo'])
        assert d['regime'] == r['regime'], tag
        assert d['mu_V_eV'] == pytest.approx(r['mu_V_eV'], abs=1e-11), tag
        assert g['shell_nm'] == pytest.approx(r['shell_nm'], rel=1e-15), tag
        for a, b in zip(g['N_layers'], r['N_layers']):
            assert a == pytest.approx(b, rel=1e-15), tag
        for k in ('surface', 'subsurface', 'bulk'):
            assert d['umol_g'][k] == pytest.approx(
                r['umol_g'][k], abs=1e-9), (tag, k)
            assert d['theta'][k] == pytest.approx(
                r['theta'][k], abs=1e-11), (tag, k)
        assert d['extended_defects_umol_g'] == pytest.approx(
            r['extended'], abs=1e-9), tag
        got = d['layers']
        assert len(got) == len(r['layer_theta']) == r['n'], tag
        for i, lay in enumerate(got):
            assert lay['theta'] == pytest.approx(
                r['layer_theta'][i], abs=1e-11), (tag, i)
            assert lay['umol_g'] == pytest.approx(
                r['layer_umol'][i], abs=1e-9), (tag, i)
            assert lay['eps_eV'] == pytest.approx(
                r['layer_eps'][i], rel=1e-15), (tag, i)
            assert lay['omega_eV'] == r['layer_omega'][i], (tag, i)
        for k in ('mu_t_eV', 'mu_cs_eV', 'VO_onset_umol_g',
                  'VO_recon_complete_umol_g', 'VO_cs_umol_g'):
            assert d['phase_boundaries'][k] == pytest.approx(
                r['boundaries'][k], abs=1e-9), (tag, k)


def test_the_implicit_occupancy_matches_across_the_port(pair):
    """The solver on its own, including the omega = 0 short circuit.

    At zero both engines return the site-exclusion isotherm rather than a
    bisection, so what separates them is one call to exp() and nothing
    else - a last-ulp difference, which is the whole reason this file
    tolerates 1e-9 on mu elsewhere. Above zero the bisection carries its
    own error on top, and the looser bound is the honest one."""
    _, js = pair
    kt = S.KB_EV * 873.15
    rows = js['theta_of_mu']
    assert len(rows) == 60, len(rows)
    for r in rows:
        got = S.theta_of_mu(r['mu'], r['eps'], r['w'], kt)
        if r['w'] == 0.0:
            assert got == pytest.approx(r['theta'], rel=1e-15), r
        else:
            assert got == pytest.approx(r['theta'], abs=1e-12), r


# ------------------------------------------------------------ page gates

def test_page_inlines_the_statmech_stack():
    """ti_solver.html carries the dataset, engine, worker and UI verbatim."""
    page = DOCS / 'ti_solver.html'
    assert page.exists(), 'run python3 scripts/build_ti_solver.py'
    html = page.read_text()
    assert (DATA / 'rutile_dft.json').read_text().strip() in html, \
        'page is stale against rutile_dft.json - run build_ti_solver.py'
    sources = (WEB / 'statmech.js', WEB / 'statmech_worker.js',
               WEB / 'statmech_ui.js',
               WEB / 'generated' / 'equations_statmech.html')
    for src in sources:
        assert src.read_text() in html, \
            f'page is stale against {src.name} - run build_ti_solver.py'
    for token in ('__SM_DATA__', '__SM_ENGINE__', '__CG_ENGINE__',
                  '__SM_WORKER__', '__SM_UI__', '__SM_EQS__'):
        assert token not in html
    assert (WEB / 'cgmc.js').read_text() in html, \
        'page is stale against cgmc.js - run build_ti_solver.py'


def test_page_carries_the_defect_workspace_and_disclosures():
    html = (DOCS / 'ti_solver.html').read_text()
    assert 'Defect stat mech' in html
    for phrase in ('Li, Guo &amp; Robertson', 'Birschitzky',
                   'not reduction kinetics', 'conditional-equilibrium',
                   'rutile_dft.json', 'Widom insertion',
                   'Placeholder kinetics', 'Recoverable fraction',
                   '-accessible inventory', '(S1)', '(S15)',
                   'Oxygen environment', 'ThermoBridge', 'phaseValidity',
                   'Particle cross-section', 'one vacancy on one oxygen site',
                   'Monte Carlo arrangement', 'Katsoulakis',
                   'OriginLab, no weighting',
                   # the phase construction and its honesty clauses
                   'added-row', 'lever rule', 'Onishi',
                   'line compound', 'extended defects',
                   'effective pair interactions',
                   'illustrative until fitted',
                   'high-coverage branch',
                   # the declared shell thickness and its sensitivity
                   'Subsurface shell thickness', 'id="smSS"',
                   'declared model choice', 'four O per cell',
                   '51.8 / 82.6 / 86.3 / 87.1'):
        assert phrase.lower() in html.lower(), f'disclosure lost: {phrase}'
    # the dataset is swappable, and the page says so
    assert 'calibrated to those constraints' in html
    assert 'will replace them directly' in html


def test_particle_view_uses_spherical_shells_without_random_bulk_dots():
    ui = (WEB / 'statmech_ui.js').read_text()
    template = (WEB / 'ti_solver_template.html').read_text()
    assert 'bulk occupation' in ui
    assert 'tint = bulk occupation' in ui
    assert "'bulk', acc.P.bulk" not in ui
    assert 'Spherical particle cross-section' in template
    assert 'core fill tracks the' in template


# The magnified slab window moved into web/figures_defect.js with the rest
# of the drawing, and its gate went with it:
# tests/test_defect_figures.py::test_the_cross_section_draws_one_dot_per_oxygen_site
# holds the isotropic scale, the 2*n_layers sub-rows and the one-dot-per-
# oxygen-site claim against the file that now owns them.


def test_portal_fits_the_desktop_viewport_and_keeps_mobile_scroll():
    portal = (WEB / 'portal.html').read_text()
    index = (DOCS / 'index.html').read_text()
    assert portal == index
    assert 'height:calc(100svh - 58px)' in portal
    assert 'grid-template-rows:auto minmax(0,1fr) auto' in portal
    assert '@media (max-width:720px){body{overflow:auto}' in portal
