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
import pathlib
import shutil
import subprocess

import pytest

from solidgas import statmech as S

ROOT = pathlib.Path(__file__).resolve().parent.parent
HARNESS = ROOT / 'tests' / 'js' / 'statmech_harness.js'
REF = json.loads((ROOT / 'reference_statmech.json').read_text())
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
    (ROOT / 'reference_results_high_precision.json').read_text())
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


@pytest.fixture(scope='module')
def pair(tmp_path_factory):
    assert shutil.which('node') is not None, \
        'node is required to run the release gate'
    spec = tmp_path_factory.mktemp('spec') / 'statmech_spec.json'
    spec.write_text(json.dumps({'mc_cases': MC_CASES,
                                'pipelines': PIPELINES,
                                'sweep_vo': SWEEP_VO,
                                'validity_cases': VALIDITY_CASES}))
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
        kt = S.KB_EV * (c['T_C'] + 273.15)
        eps = S.energetics(P)['eps']
        fn = S.theta_s_interp(out['isotherm'], eps[0], kt)
        out['sweep'] = S.sweep(P, c['T_C'], out['geometry'], SWEEP_VO,
                               theta_s_fn=fn, eps=eps)
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
        for k in ('surface', 'subsurface', 'bulk'):
            assert row['fractions'][k] == pytest.approx(
                ref['fractions'][k], rel=1e-12, abs=1e-15), k
        assert row['surface_cap_bound'] == ref['surface_cap_bound']
        assert row['warn_kinds'] == [w['kind'] for w in ref['warnings']]
        for a, b in zip(row['isotherm'], ref['isotherm']):
            assert a['filling'] == b['filling']
            for k in ('theta_s', 'theta_ss', 'theta_b', 'E_mean_eV'):
                assert a[k] == b[k], k
            assert a['hist'] == b['hist']
            assert a['mu_eV'] == pytest.approx(b['mu_eV'], abs=1e-9)


def test_shipped_configuration_reproduces_the_oracle_flagship(pair):
    """The default lattice and sweep counts, exactly as the page runs them.
    At the flagship condition the reconstruction cap binds, so the answer
    is analytic and the 50-digit oracle certifies the full JS pipeline."""
    _, js = pair
    row = next(r for r in js['pipelines'] if r['name'] == 'shipped_default')
    ref = REF['flagship']
    assert row['mu_V_eV'] == pytest.approx(ref['mu_V_eV'], abs=1e-9)
    for k in ('surface', 'subsurface', 'bulk'):
        assert row['fractions'][k] == pytest.approx(
            ref['fractions'][k], rel=1e-9), k
        assert row['umol_g'][k] == pytest.approx(
            ref['umol_g'][k], rel=1e-9), k
    assert row['surface_cap_bound'] is True
    assert row['warn_kinds'] == ref['warn_kinds']
    sp = row['spacing']
    assert sp is not None
    assert sp['modal_gap_sites'] in (3, 4, 5)
    assert sp['P_gap_le_2'] < 0.15


def test_geometry_port_matches_oracle(pair):
    """Independent of MC: the JS geometry numbers against the oracle."""
    _, js = pair
    row = next(r for r in js['pipelines'] if r['name'] == 'shipped_default')
    # fractions already checked; anchor the derived quantities too
    assert row['theta']['subsurface'] == pytest.approx(
        REF['flagship']['theta']['subsurface'], rel=1e-9)
    assert row['theta']['bulk'] == pytest.approx(
        REF['flagship']['theta']['bulk'], rel=1e-9)


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
            for k in ('surface', 'subsurface', 'bulk', 'accessible',
                      'f_recoverable'):
                assert ra[k] == pytest.approx(rb[k], rel=1e-12,
                                              abs=1e-15), k


def test_shipped_accessibility_matches_the_oracle(pair):
    """The default kinetics on the default distribution: the 300 s oracle
    row certifies the JS chain (the flagship split is cap-bound analytic)."""
    _, js = pair
    row = next(r for r in js['pipelines'] if r['name'] == 'shipped_default')
    ref = next(r for r in REF['co2_flagship'] if r['exposure_s'] == 300)
    a = row['accessibility']
    for c in ('surface', 'subsurface', 'bulk'):
        assert a['P'][c] == pytest.approx(ref['P'][c], rel=1e-9,
                                          abs=1e-12), c
    assert a['accessible_umol_g'] == pytest.approx(
        ref['accessible_umol_g'], rel=1e-9)
    assert a['f_recoverable'] == pytest.approx(
        ref['f_recoverable'], rel=1e-9)
    # the 95 row of the shipped sweep is the flagship distribution
    r95 = next(r for r in row['sweep'] if r['VO_total_umol_g'] == 95.0)
    for k in ('surface', 'subsurface', 'bulk'):
        assert r95[k] == pytest.approx(REF['flagship']['umol_g'][k],
                                       rel=1e-9), k


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


# ------------------------------------------------------------ page gates

def test_page_inlines_the_statmech_stack():
    """ti_solver.html carries the dataset, engine, worker and UI verbatim."""
    page = ROOT / 'ti_solver.html'
    assert page.exists(), 'run python3 build_ti_solver.py'
    html = page.read_text()
    assert (ROOT / 'rutile_dft.json').read_text().strip() in html, \
        'page is stale against rutile_dft.json - run build_ti_solver.py'
    for src in ('statmech.js', 'statmech_worker.js', 'statmech_ui.js',
                'equations_statmech.html'):
        assert (ROOT / src).read_text() in html, \
            f'page is stale against {src} - run build_ti_solver.py'
    for token in ('__SM_DATA__', '__SM_ENGINE__', '__SM_WORKER__',
                  '__SM_UI__', '__SM_EQS__'):
        assert token not in html


def test_page_carries_the_defect_workspace_and_disclosures():
    html = (ROOT / 'ti_solver.html').read_text()
    assert 'Defect stat mech' in html
    for phrase in ('Li, Guo &amp; Robertson', 'Birschitzky',
                   'reconstruction threshold', 'not reduction kinetics',
                   'rutile_dft.json', 'Widom insertion',
                   'Placeholder kinetics', 'Recoverable fraction',
                   '-accessible inventory', '(S1)', '(S12)',
                   'Oxygen environment', 'ThermoBridge', 'phaseValidity'):
        assert phrase in html, f'disclosure lost: {phrase}'
    # the dataset is swappable, and the page says so
    assert 'replaced' in html and 're-fitted' in html
