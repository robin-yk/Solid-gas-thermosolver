"""The defect stat-mech engine against its independent reference.

Three certified layers, mirroring the thermodynamic side: an mpmath oracle
(analytic matching + exact finite-lattice ensembles), the Python production
engine, and - in test_statmech_port.py - the JS port. Nothing here may skip:
a missing artefact is a failure.
"""

import json
import math
import pathlib

import pytest

from solidgas import statmech as S

ROOT = pathlib.Path(__file__).resolve().parent.parent
REF = json.loads((ROOT / 'reference_statmech.json').read_text())
P = S.load_params()


# ------------------------------------------------------------ freshness

def test_oracle_reproduces_committed_reference():
    """Re-derive the flagship row in-process: catches a stale JSON."""
    import oracle_statmech as O
    geom = O.geometry(d_um=P['defaults']['particle_diameter_um'])
    row = O.match(P['defaults']['T_C'], P['defaults']['VO_total_umol_g'],
                  geom, 'sx')
    ref = REF['flagship']
    assert row['mu_V_eV'] == pytest.approx(ref['mu_V_eV'], abs=1e-12)
    for k in ('surface', 'subsurface', 'bulk'):
        assert row['fractions'][k] == pytest.approx(
            ref['fractions'][k], rel=1e-12)
    assert row['warn_kinds'] == ref['warn_kinds']
    assert row['regime'] == 'unresolved'
    assert row['VO_valid_max_umol_g'] == pytest.approx(
        ref['VO_valid_max_umol_g'], rel=1e-12)


def test_rng_matches_the_specification():
    rng = S.make_rng(12345)
    got = [rng() for _ in range(8)]
    assert got == REF['rng_sfc32_seed12345_first8']


# ------------------------------------------------------------ geometry

def test_geometry_matches_oracle():
    g = S.geometry(P)
    for k, v in REF['geometry_0p9um'].items():
        assert g[k] == pytest.approx(v, rel=1e-10), k
    gb = S.geometry(P, bet_m2_g=1.60)
    for k, v in REF['geometry_bet_1p60'].items():
        assert gb[k] == pytest.approx(v, rel=1e-10), k


def test_hand_anchors():
    """Numbers checked by hand: they pin the unit system."""
    g = S.geometry(P)
    assert g['area_m2_g'] == pytest.approx(1.569, abs=0.002)
    assert g['N_s'] == pytest.approx(13.56, abs=0.02)       # umol-O/g
    assert g['N_O_total'] == pytest.approx(25042.0, abs=1.0)
    r = S.distribute(P, 600.0, 95.0, g)
    assert r['x_TiO2mx'] == pytest.approx(0.00759, abs=1e-5)
    assert r['Ti3_frac'] == pytest.approx(0.0152, abs=5e-5)


# ------------------------------------------------- analytic matching

def _assert_row(row, ref):
    assert row['mu_V_eV'] == pytest.approx(ref['mu_V_eV'], abs=1e-9)
    for k in ('surface', 'subsurface', 'bulk'):
        assert row['fractions'][k] == pytest.approx(
            ref['fractions'][k], rel=1e-9, abs=1e-12), k
        assert row['umol_g'][k] == pytest.approx(
            ref['umol_g'][k], rel=1e-9, abs=1e-12), k
    assert row['validity']['regime'] == ref['regime']
    assert row['validity']['binding'] == ref['binding']
    assert row['validity']['VO_valid_max_umol_g'] == pytest.approx(
        ref['VO_valid_max_umol_g'], rel=1e-6)
    assert [w['kind'] for w in row['warnings']
            if w['kind'] != 'not_converged'] == ref['warn_kinds']


def test_analytic_grid_matches_oracle():
    g = S.geometry(P)
    for ref in REF['analytic_grid']:
        row = S.distribute(P, ref['T_C'], ref['VO_total_umol_g'], g,
                           preset=ref['preset'])
        _assert_row(row, ref)


def test_flagship_both_presets():
    g = S.geometry(P)
    _assert_row(S.distribute(P, 600.0, 95.0, g, preset='sx'),
                REF['flagship'])
    _assert_row(S.distribute(P, 600.0, 95.0, g, preset='gga'),
                REF['flagship_gga'])


def test_bet_mode_matches_oracle():
    g = S.geometry(P, bet_m2_g=1.60)
    _assert_row(S.distribute(P, 600.0, 95.0, g), REF['bet_case'])


def test_presets_keep_the_literature_ordering():
    for key in P['vacancy_energetics']['presets']:
        e = S.energetics(P, key)['eps']
        assert e[0] == 0.0 and 0.0 < e[1] < e[2], key


# ------------------------------------------- exact finite-lattice MC

def _mc_case(rows, row_sites, ss_layers, bulk_layers, ordering, eps, t_c,
             n_vac, sweeps_eq, sweeps_sample, seed=902, **kw):
    lat = S.build_lattice(rows, row_sites, ss_layers, bulk_layers, ordering)
    rng = S.make_rng(seed)
    occ = bytearray(lat['n'])
    vac = []
    S.adjust_filling(lat, occ, vac, n_vac, rng)
    kt = S.KB_EV * (t_c + 273.15)
    return S.mc_run(lat, eps, kt, occ, vac, sweeps_eq, sweeps_sample, rng,
                    **kw)


NO_PAIRS = {'pair_eV': {'along_row': [], 'cross_row': []}}


def test_ideal_ensemble_occupancies_and_widom_mu_are_exact():
    """48/48/96 ideal lattice vs the generating-polynomial ensemble."""
    ref = REF['ideal_ensemble']
    res = _mc_case(4, 12, 1, 2, NO_PAIRS, ref['eps_eV'], ref['T_C'],
                   ref['n_vac'], 1500, 8000, widom_every=5)
    for k, name in ((0, 'theta_s'), (1, 'theta_ss'), (2, 'theta_b')):
        exact = ref['theta_exact'][k]
        err = res['err_' + name.split('_')[1]]
        assert abs(res[name] - exact) < max(6 * err, 2e-4), (name, exact)
    assert abs(res['mu_eV'] - ref['mu_exact_eV']) < 0.006
    assert res['E_drift_eV'] < 1e-9


def test_interacting_ring_matches_exact_enumeration():
    ref = REF['ring_case']
    ordering = {'pair_eV': {'along_row': ref['pair_eV'], 'cross_row': []}}
    res = _mc_case(1, ref['row_sites'], 0, 0, ordering, [0.0, 0.0, 0.0],
                   ref['T_C'], ref['n_vac'], 3000, 40000,
                   widom_every=5, hist_every=2)
    assert abs(res['E_mean_eV'] - ref['E_mean_eV']) < 0.004
    tot = sum(res['hist'][1:])
    for gap in range(1, ref['row_sites']):
        got = res['hist'][gap] / tot
        assert abs(got - ref['gap_probs'][gap]) < 0.01, gap
    assert abs(res['mu_eV'] - ref['mu_exact_eV']) < 0.01
    assert res['E_drift_eV'] < 1e-9


def test_interacting_three_class_micro_lattice_is_exact():
    ref = REF['micro3_case']
    ordering = {'pair_eV': {'along_row': ref['pair_eV'], 'cross_row': []}}
    res = _mc_case(1, ref['row_sites'], 1, 1, ordering, ref['eps_eV'],
                   ref['T_C'], ref['n_vac'], 3000, 60000, widom_every=5)
    for k, name in ((0, 'theta_s'), (1, 'theta_ss'), (2, 'theta_b')):
        err = res['err_' + name.split('_')[1]]
        assert abs(res[name] - ref['theta_exact'][k]) < max(6 * err, 2e-3), name
    assert abs(res['E_mean_eV'] - ref['E_mean_eV']) < 0.006
    assert abs(res['mu_eV'] - ref['mu_exact_eV']) < 0.01


def test_hard_first_neighbour_repulsion_empties_the_gap():
    ordering = {'pair_eV': {'along_row': [5.0], 'cross_row': []}}
    res = _mc_case(1, 12, 0, 0, ordering, [0.0, 0.0, 0.0], 600.0, 4,
                   2000, 20000, hist_every=2)
    tot = sum(res['hist'][1:])
    assert res['hist'][1] / tot < 1e-3


# ------------------------------------------------------ full pipeline

MC_SMALL = {'rows': 4, 'row_sites': 12, 'ss_layers': 1, 'bulk_layers': 2,
            'fillings': [2, 4, 8, 14, 22, 34], 'equil_sweeps': 600,
            'warm_equil_sweeps': 300, 'sample_sweeps': 1600,
            'widom_every': 5, 'hist_every': 20, 'blocks': 16, 'seed': 12345}


def test_pipeline_with_zero_pairs_reproduces_the_analytic_solution():
    """Inside the resolved regime and the sampled mu-range, the J = 0
    Monte Carlo isotherm is the ideal one, so the full pipeline lands on
    the analytic partition to sampling noise."""
    out = S.run_full(P, vo_total=1.8, zero_pairs=True, mc=MC_SMALL)
    ref = S.distribute(P, 600.0, 1.8, S.geometry(P))
    assert out['validity']['regime'] == 'valid'
    for k in ('surface', 'subsurface', 'bulk'):
        assert out['fractions'][k] == pytest.approx(
            ref['fractions'][k], abs=0.02), k


def test_pipeline_is_deterministic_for_a_seed():
    a = S.run_full(P, mc=MC_SMALL)
    b = S.run_full(P, mc=MC_SMALL)
    assert a['isotherm'] == b['isotherm']
    assert a['mu_V_eV'] == b['mu_V_eV']
    assert a['fractions'] == b['fractions']


def test_low_loading_lives_on_the_surface():
    """Below the validity boundary the equilibrium expectation is
    surface-first: the sX energies put subsurface 0.59 eV up."""
    g = S.geometry(P)
    r = S.distribute(P, 600.0, 1.0, g)
    assert r['fractions']['surface'] > 0.995
    assert r['fractions']['bulk'] < 1e-3
    assert r['validity']['regime'] == 'valid'


def test_validity_boundary_brackets_the_regimes():
    """The resolved regime ends where the unconstrained surface coverage
    crosses the reconstruction boundary - nothing is pinned beyond it,
    the split is simply declared unresolved."""
    g = S.geometry(P)
    vmax = S.validity_boundary(P, 600.0, g)
    assert vmax == pytest.approx(REF['flagship']['VO_valid_max_umol_g'],
                                 rel=1e-6)
    assert vmax == pytest.approx(2.309, abs=0.01)
    below = S.distribute(P, 600.0, vmax * 0.98, g)
    above = S.distribute(P, 600.0, vmax * 1.02, g)
    assert below['validity']['regime'] == 'valid'
    assert above['validity']['regime'] == 'unresolved'
    assert above['validity']['binding'] == 'surface'
    kinds = [w['kind'] for w in above['warnings']]
    assert 'surface_reconstruction_regime' in kinds
    # hotter surface holds relatively less before the boundary at fixed
    # coverage threshold, but the deeper classes activate: the boundary
    # inventory grows with temperature
    assert S.validity_boundary(P, 1000.0, g) > vmax


def test_unresolved_extrapolation_is_labelled_not_pinned():
    """At the flagship inventory the model must say 'unresolved', carry
    the breakdown warnings verbatim, and leave the unconstrained numbers
    as an extrapolation - no pinning, no reassignment."""
    g = S.geometry(P)
    r = S.distribute(P, 600.0, 95.0, g)
    assert r['validity']['regime'] == 'unresolved'
    assert r['validity']['binding'] == 'surface'
    texts = ' '.join(w['text'] for w in r['warnings'])
    assert 'cannot be assigned quantitatively' in texts
    assert 'dilute-defect approximation fails' in texts
    # unconstrained: the surface is NOT pinned at any threshold
    assert r['theta']['surface'] > 0.9
    assert r['matched_umol_g'] == pytest.approx(95.0, rel=1e-9)


def test_high_temperature_limit_is_proportional_to_site_counts():
    g = S.geometry(P)
    r = S.distribute(P, 1.0e6, 2000.0, g)       # kT ~ 86 eV >> every gap
    for k, n in (('surface', 'N_s'), ('subsurface', 'N_ss'), ('bulk', 'N_b')):
        assert r['fractions'][k] == pytest.approx(
            g[n] / g['N_O_total'], rel=0.03), k


def test_deeper_bulk_level_moves_inventory_out_of_the_bulk():
    g = S.geometry(P)
    lo = S.distribute(P, 600.0, 95.0, g, eps=[0.0, 0.59, 1.10])
    hi = S.distribute(P, 600.0, 95.0, g, eps=[0.0, 0.59, 1.60])
    assert hi['fractions']['bulk'] < lo['fractions']['bulk']


def test_warning_fires_beyond_the_shear_limit():
    g = S.geometry(P)
    r = S.distribute(P, 600.0, 120.0, g)
    assert 'x_max' in [w['kind'] for w in r['warnings']]


def test_isotherm_theta_s_is_monotone_in_filling():
    pts = S.isotherm_scan(P, 600.0, mc=MC_SMALL)
    for a, b in zip(pts, pts[1:]):
        assert b['theta_s'] >= a['theta_s'] - 3 * (a['err_s'] + b['err_s'])


def test_isotherm_points_carry_a_valid_surface_snapshot():
    pts = S.isotherm_scan(P, 600.0, mc=MC_SMALL)
    n_surf = MC_SMALL['rows'] * MC_SMALL['row_sites']
    for q in pts:
        assert q['rows'] == MC_SMALL['rows']
        assert q['row_sites'] == MC_SMALL['row_sites']
        sv = q['surf_vac']
        assert sv == sorted(set(sv))
        assert all(0 <= v < n_surf for v in sv)
        # the closing configuration sits near the run's own average
        assert abs(len(sv) / n_surf - q['theta_s']) < 0.12


def test_spacing_statistics_sit_in_the_literature_window():
    """Near critical coverage: modal gap 4-5 bridging sites, close pairs
    suppressed - the constraints the pair table was calibrated to."""
    mc = {'rows': 4, 'row_sites': 25, 'ss_layers': 1, 'bulk_layers': 2,
          'fillings': [12, 17, 22], 'equil_sweeps': 800,
          'warm_equil_sweeps': 400, 'sample_sweeps': 2000,
          'widom_every': 10, 'hist_every': 4, 'blocks': 16, 'seed': 12345}
    pts = S.isotherm_scan(P, 600.0, mc=mc)
    q = min(pts, key=lambda r: abs(r['theta_s'] - 0.17))
    hist = q['hist']
    tot = sum(hist[1:])
    modal = max(range(1, len(hist)), key=lambda gp: hist[gp])
    assert modal in (4, 5)
    assert (hist[1] + hist[2]) / tot < 0.08


# --------------------------------------------------------- CO2 layer

def test_accessibility_matches_oracle_at_flagship():
    for ref in REF['co2_flagship']:
        a = S.accessibility(P, REF['flagship']['umol_g'],
                            kin={'exposure_s': ref['exposure_s'],
                                 'T_reox_C': ref['T_reox_C'],
                                 'p_CO2_atm': ref['p_CO2_atm']})
        for c in ('surface', 'subsurface', 'bulk'):
            assert a['P'][c] == pytest.approx(ref['P'][c], rel=1e-12,
                                              abs=1e-15), c
        assert a['accessible_umol_g'] == pytest.approx(
            ref['accessible_umol_g'], rel=1e-12)
        assert a['f_recoverable'] == pytest.approx(
            ref['f_recoverable'], rel=1e-12)


def test_accessibility_limits():
    u = REF['flagship']['umol_g']
    z = S.accessibility(P, u, kin={'exposure_s': 0})
    assert all(v == 0.0 for v in z['P'].values())
    assert z['accessible_umol_g'] == 0.0
    full = S.accessibility(P, u, kin={'exposure_s': 1e30})
    assert full['P']['surface'] == 1.0
    assert full['P']['subsurface'] == 1.0
    off = S.accessibility(P, u, kin={'A_eff_per_s_atm': {
        'surface': 0.0, 'subsurface': 0.0, 'bulk': 0.0}})
    assert off['accessible_umol_g'] == 0.0


def test_accessibility_rises_with_reoxidation_temperature():
    u = REF['flagship']['umol_g']
    lo = S.accessibility(P, u, kin={'T_reox_C': 500})
    hi = S.accessibility(P, u, kin={'T_reox_C': 700})
    assert hi['P']['subsurface'] > lo['P']['subsurface']
    assert hi['f_recoverable'] > lo['f_recoverable']


def test_sweep_agrees_with_pointwise_solutions():
    g = S.geometry(P)
    rows = S.sweep(P, 600.0, g, [10.0, 95.0, 150.0])
    for row in rows:
        d = S.distribute(P, 600.0, row['VO_total_umol_g'], g)
        a = S.accessibility(P, d['umol_g'])
        assert row['surface'] == d['umol_g']['surface']
        assert row['subsurface'] == d['umol_g']['subsurface']
        assert row['bulk'] == d['umol_g']['bulk']
        assert row['accessible'] == a['accessible_umol_g']


def test_sweep_carries_the_validity_split():
    """The sweep is the paper's panel: solid within the resolved regime,
    declared unresolved beyond the boundary - and inside the resolved
    regime the inventory is essentially all surface."""
    g = S.geometry(P)
    vmax = S.validity_boundary(P, 600.0, g)
    rows = S.sweep(P, 600.0, g,
                   [0.5 * vmax, 0.9 * vmax, 20.0, 95.0, 150.0])
    assert rows[0]['valid'] and rows[1]['valid']
    assert not rows[2]['valid'] and not rows[4]['valid']
    assert rows[0]['surface'] / (0.5 * vmax) > 0.99
    assert rows[0]['f_recoverable'] > 0.95
    acc = [r['accessible'] for r in rows]
    assert all(b >= a - 1e-9 for a, b in zip(acc, acc[1:]))


def test_dielectric_power_law_of_this_work():
    """The measured loss-factor fit: a log10-log10 regression, so the
    correlation is eps'' = 10^a * VO^b, and the provenance stays intact."""
    import copy
    c = P['dielectric_correlation']
    d = S.dielectric(P, 95.0)
    assert d['eps2'] == pytest.approx(
        10.0 ** (c['log10_intercept']
                 + c['exponent'] * math.log10(95.0)), rel=1e-12)
    assert d['eps2'] == pytest.approx(0.2238, abs=2e-3)
    assert d['extrapolated'] is False
    assert 'tan_delta' not in d                 # eps_prime not measured yet
    assert S.dielectric(P, 10.0)['extrapolated'] is True
    assert S.dielectric(P, 0.0) is None
    q = copy.deepcopy(P)
    q['dielectric_correlation']['eps_prime'] = 2.0
    assert S.dielectric(q, 95.0)['tan_delta'] == pytest.approx(
        d['eps2'] / 2.0)
    note = c['note']
    for token in ('OriginLab', '1.22355', '-3.07028', '0.93388'):
        assert token in note, token
    assert 'this work' in c['provenance']


def test_kinetics_are_declared_placeholders():
    """The JSON must keep saying these are not literature numbers."""
    note = P['co2_kinetics']['note']
    assert 'Placeholder' in note or 'placeholder' in note
    assert 'EFFECTIVE' in note or 'effective' in note


# ------------------------------------------------- thermodynamic bridge

def _is_reduced():
    from solidgas import activeset as A
    return {s: A.STOICH[s][1] / A.STOICH[s][0] < 2.0 - 1e-12
            for s in A.SOLIDS}


def test_phase_validity_stable_under_rwgs():
    from solidgas import activeset as A
    r = A.solve({'CO2': 1, 'H2': 1}, 873.15)
    v = S.phase_validity(r['active_condensed_phases'],
                         r['inactive_phase_reduced_costs'], _is_reduced())
    assert v['tier'] == 'stable'
    assert v['nearest_phase'] == 'Ti10O19'
    assert v['margin_per_mol_O_kJ'] == pytest.approx(66.6399, abs=0.01)
    assert r['lambda_kJ_per_mol']['O'] < 0


def test_phase_validity_invalid_under_hot_dry_hydrogen():
    from solidgas import activeset as A
    r = A.solve({'H2': 1}, 1473.15, mass_g=0.002)
    v = S.phase_validity(r['active_condensed_phases'],
                         r['inactive_phase_reduced_costs'], _is_reduced())
    assert v['tier'] == 'invalid'
    assert 'TiO2' not in r['active_condensed_phases']
    assert v['reduced_active'] == r['active_condensed_phases']


def test_phase_validity_boundary_at_the_h2_h2o_flip():
    import json as _json
    from solidgas import activeset as A
    ref = _json.loads(
        (pathlib.Path(__file__).resolve().parent.parent
         / 'reference_results_high_precision.json').read_text())
    ratio = float(ref['boundary_validation']['h2_h2o_boundary']
                  ['H2O_over_H2_critical_ratio'])
    r = A.solve({'H2': 1, 'H2O': ratio * (1 - 1e-3)}, 1173.15)
    v = S.phase_validity(r['active_condensed_phases'],
                         r['inactive_phase_reduced_costs'], _is_reduced())
    assert v['tier'] == 'boundary'
    assert v['margin_per_mol_O_kJ'] == 0.0


def test_phase_validity_skips_divergent_costs():
    """-inf sentinels (dry-gas divergence) must not become the margin."""
    v = S.phase_validity(
        ['TiO2'],
        {'Ti10O19': {'per_mol_O_kJ': 12.5},
         'Ti4O7': {'per_mol_O_kJ': -math.inf},
         'Ti2O3': {'per_mol_O_kJ': None}},
        {'TiO2': False, 'Ti10O19': True, 'Ti4O7': True, 'Ti2O3': True})
    assert v['tier'] == 'stable'
    assert v['nearest_phase'] == 'Ti10O19'
    assert v['margin_per_mol_O_kJ'] == 12.5


def test_mode_note_discloses_scope():
    out = S.run_full(P, mc=MC_SMALL)
    assert out['method'] == 'canonical_statmech'
    assert out['validity']['regime'] == 'unresolved'
    kinds = [w['kind'] for w in out['warnings']]
    assert 'surface_reconstruction_regime' in kinds
