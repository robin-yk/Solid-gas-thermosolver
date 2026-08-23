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
DATA = ROOT / 'data'
REF = json.loads((DATA / 'reference_statmech.json').read_text())
P = S.load_params()


# ------------------------------------------------------------ freshness

def test_oracle_reproduces_committed_reference():
    """Re-derive the flagship row in-process: catches a stale JSON."""
    from scripts import oracle_statmech as O
    geom = O.geometry(d_um=P['defaults']['particle_diameter_um'])
    row = O.match(P['defaults']['T_C'], P['defaults']['VO_total_umol_g'],
                  geom, 'sx')
    ref = REF['flagship']
    assert row['mu_V_eV'] == pytest.approx(ref['mu_V_eV'], abs=1e-12)
    for k in ('surface', 'subsurface', 'bulk'):
        assert row['fractions'][k] == pytest.approx(
            ref['fractions'][k], rel=1e-12)
    assert row['warn_kinds'] == ref['warn_kinds']


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


def test_the_shell_is_four_oxygens_per_cell_per_declared_layer():
    """The subsurface class is a slab, not an atomic row.

    One d110 = a/sqrt(2) oxygen slab carries four O per (1x1) cell against
    the bridging row's one, and that ratio is what makes a 51.8 umol-O/g
    subsurface population fit inside a 0.33 nm layer at all."""
    g1 = S.geometry(P, ss_layers=1)
    assert g1['N_ss'] == pytest.approx(4.0 * g1['N_s'], rel=1e-12)
    d110_nm = P['lattice_constants_A']['a'] / math.sqrt(2.0) / 10.0
    assert g1['layer_nm'] == pytest.approx(d110_nm, rel=1e-12)
    # the shell holds 1 + 4n O per cell, i.e. (1+4n)/4 slabs of material
    for n in (1, 2, 3, 4):
        g = S.geometry(P, ss_layers=n)
        assert g['N_ss'] == pytest.approx(4.0 * n * g['N_s'], rel=1e-12)
        assert g['shell_nm'] == pytest.approx(
            (1.0 + 4.0 * n) / 4.0 * d110_nm, rel=1e-12)
        assert g['N_s'] + g['N_ss'] + g['N_b'] == pytest.approx(
            g['N_O_total'], rel=1e-12)
    # the areal layer density is the bulk oxygen density over one slab
    sig_b, sig_l = S.site_densities(P)
    a = P['lattice_constants_A']['a']
    c = P['lattice_constants_A']['c']
    assert sig_l / (d110_nm * 1e-7) == pytest.approx(
        4.0 / (a * a * c * 1e-24), rel=1e-12)
    assert sig_l == pytest.approx(4.0 * sig_b, rel=1e-12)


def test_the_shell_thickness_ladder_matches_the_oracle():
    """Declared shell thickness moves the shell/bulk split, so certify it.

    The surface sits on the reconstructed line phase throughout, so the
    ladder only trades subsurface against bulk; the shell share rises
    monotonically and mu_V falls as the cheap class is made thicker."""
    assert len(REF['shell_ladder']) == 4
    shares, mus = [], []
    for ref in REF['shell_ladder']:
        g = S.geometry(P, ss_layers=ref['ss_layers'])
        for k, v in ref['geometry'].items():
            assert g[k] == pytest.approx(v, rel=1e-10), k
        row = S.distribute(P, ref['T_C'], ref['VO_total_umol_g'], g,
                           preset=ref['preset'])
        _assert_row(row, ref)
        assert row['umol_g']['surface'] == pytest.approx(
            0.5 * g['N_s'], rel=1e-12)
        shares.append((row['umol_g']['surface'] + row['umol_g']['subsurface'])
                      / row['VO_total_umol_g'])
        mus.append(row['mu_V_eV'])
    assert shares == sorted(shares)
    assert mus == sorted(mus, reverse=True)
    assert shares[0] == pytest.approx(0.616, abs=0.005)
    assert shares[1] == pytest.approx(0.941, abs=0.005)


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
    assert row['regime'] == ref['regime']
    assert row['surface_phase']['phase'] == ref['surface_phase']['phase']
    assert row['surface_phase']['phi_reconstructed'] == pytest.approx(
        ref['surface_phase']['phi_reconstructed'], abs=1e-9)
    assert row['extended_defects_umol_g'] == pytest.approx(
        ref['extended_defects_umol_g'], rel=1e-9, abs=1e-9)
    for k, v in ref['boundaries'].items():
        assert row['phase_boundaries'][k] == pytest.approx(
            v, rel=1e-9, abs=1e-12), k
    assert row['surface_reconstruction_regime'] == \
        ref['surface_reconstruction_regime']
    assert [w['kind'] for w in row['warnings']
            if w['kind'] != 'unresolved'] == ref['warn_kinds']


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
    out = S.run_full(P, zero_pairs=True, mc=MC_SMALL)
    ref = REF['flagship']
    assert out['mu_V_eV'] == pytest.approx(ref['mu_V_eV'], abs=1e-9)
    for k in ('surface', 'subsurface', 'bulk'):
        assert out['fractions'][k] == pytest.approx(
            ref['fractions'][k], rel=1e-9), k


def test_pipeline_is_deterministic_for_a_seed():
    a = S.run_full(P, mc=MC_SMALL)
    b = S.run_full(P, mc=MC_SMALL)
    assert a['isotherm'] == b['isotherm']
    assert a['mu_V_eV'] == b['mu_V_eV']
    assert a['fractions'] == b['fractions']


def test_low_loading_lives_on_the_surface():
    """Below the surface capacity the equilibrium expectation is
    surface-first: the sX energies put subsurface 0.59 eV up."""
    g = S.geometry(P)
    r = S.distribute(P, 600.0, 1.0, g)
    assert r['fractions']['surface'] > 0.995
    assert r['fractions']['bulk'] < 1e-3
    assert not r['surface_reconstruction_regime']


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


def test_the_phase_ladder_and_its_lever_rules():
    """The four regimes in order, mass exact in each, mu pinned on the
    coexistence branches, and continuity across every boundary."""
    g = S.geometry(P)
    b = S.phase_boundaries(P, 600.0, g)
    cases = [(1.0, 'dilute'), (3.5, 'surface_coexistence'),
             (9.0, 'reconstructed'), (95.0, 'reconstructed'),
             (200.0, 'cs_coexistence')]
    for vo, regime in cases:
        d = S.distribute(P, 600.0, vo, g)
        assert d['regime'] == regime, vo
        assert d['matched_umol_g'] == pytest.approx(vo, rel=1e-9), vo
    co = S.distribute(P, 600.0, 3.5, g)
    assert co['mu_V_eV'] == b['mu_t_eV']                 # riser: mu pinned
    assert 0.0 < co['surface_phase']['phi_reconstructed'] < 1.0
    cs = S.distribute(P, 600.0, 200.0, g)
    assert cs['mu_V_eV'] == b['mu_cs_eV']                # ceiling: mu pinned
    assert cs['extended_defects_umol_g'] == pytest.approx(
        200.0 - b['VO_cs_umol_g'], rel=1e-6)
    assert 'cs_precipitation' in [w['kind'] for w in cs['warnings']]
    # continuity at each boundary: the two sides agree to the solve tol
    for edge in ('VO_onset_umol_g', 'VO_recon_complete_umol_g',
                 'VO_cs_umol_g'):
        lo = S.distribute(P, 600.0, b[edge] * (1 - 1e-9), g)
        hi = S.distribute(P, 600.0, b[edge] * (1 + 1e-9), g)
        for k in ('surface', 'subsurface', 'bulk'):
            assert hi['umol_g'][k] == pytest.approx(
                lo['umol_g'][k], abs=5e-6), (edge, k)


def test_boundary_ladder_matches_the_hand_calculation():
    """mu_t = eps_s + kT ln(th/(1-th)) and friends, checked as numbers:
    2.31 / 6.78 / 160 umol-O/g at 600 C for the shipped dataset."""
    g = S.geometry(P)
    b = S.phase_boundaries(P, 600.0, g)
    assert b['mu_t_eV'] == pytest.approx(-0.11931, abs=2e-5)
    assert b['mu_cs_eV'] == pytest.approx(0.89485, abs=2e-5)
    assert b['VO_onset_umol_g'] == pytest.approx(2.3092, abs=2e-4)
    assert b['VO_recon_complete_umol_g'] == pytest.approx(6.7830, abs=2e-4)
    assert b['VO_cs_umol_g'] == pytest.approx(159.98, abs=0.05)
    hot = S.phase_boundaries(P, 800.0, g)
    assert hot['VO_onset_umol_g'] > b['VO_onset_umol_g']


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


def test_sweep_walks_the_phase_ladder():
    """Along the inventory axis: surface rises through the coexistence
    riser, holds at the line-phase deficiency, extended defects appear
    only beyond the CS ceiling, and the recoverable fraction falls as the
    inventory is forced deep."""
    g = S.geometry(P)
    rows = S.sweep(P, 600.0, g, [2.0, 4.0, 20.0, 60.0, 95.0, 150.0, 200.0])
    b = S.phase_boundaries(P, 600.0, g)
    line_umol = b['theta_reconstructed_eff'] * g['N_s']
    surface = [r['surface'] for r in rows]
    assert all(y >= a - 1e-12 for a, y in zip(surface, surface[1:]))
    assert surface[0] < line_umol                     # dilute point
    assert rows[1]['regime'] == 'surface_coexistence'
    for r in rows[2:]:
        assert r['surface'] == pytest.approx(line_umol, rel=1e-12)
    assert [r['regime'] for r in rows][-1] == 'cs_coexistence'
    assert rows[-1]['extended'] > 0.0
    assert all(r['extended'] == 0.0 for r in rows[:-1])
    fr = [r['f_recoverable'] for r in rows]
    assert fr[0] > 0.95
    assert fr[-1] < fr[-2] < fr[2]
    acc = [r['accessible'] for r in rows[:-1]]
    assert all(y >= a for a, y in zip(acc, acc[1:]))


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
        (DATA / 'reference_results_high_precision.json').read_text())
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
    assert out['surface_reconstruction_regime'] is True
    assert out['regime'] == 'reconstructed'
    kinds = [w['kind'] for w in out['warnings']]
    assert 'surface_phase' in kinds and 'subsurface_dense' in kinds
    phase_note = next(w['text'] for w in out['warnings']
                      if w['kind'] == 'surface_phase')
    assert 'added-row' in phase_note
    assert 'line compound' in phase_note
    shell = next(w['text'] for w in out['warnings']
                 if w['kind'] == 'subsurface_dense')
    assert 'dilute-defect approximation fails' in shell
    assert 'heavily reduced' in shell
