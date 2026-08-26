"""Freeze every number the manuscript plates draw into data/figure_data.json.

The proof sheet does not run a solver. It renders one frozen data object,
extracted here from the production engines and the certified references and
stamped with the commit it came from, so a figure and the manuscript text
cannot drift apart. Regenerating after an engine change means re-running this
script; tests/test_plates.py fails while the freeze is stale.

Two products:

  slots  flat, preformatted strings keyed by the dotted path a caption uses,
         so caption prose never carries a typed number.
  plot   the arrays and records the plate kit draws.

    python3 scripts/export_plates.py
"""

import json
import math
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / 'data'
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from solidgas import statmech as S            # noqa: E402
from solidgas import cgmc as CG               # noqa: E402
import figure_registry as REG                 # noqa: E402

P = S.load_params()
REF_SM = json.loads((DATA / 'reference_statmech.json').read_text())
REF_TH = json.loads(
    (DATA / 'reference_results_high_precision.json').read_text())

FLAG_T = P['defaults']['T_C']
FLAG_VO = P['defaults']['VO_total_umol_g']
CGMC_TARGET_PCT = 40.0
CGMC_TARGET_S = 300.0
MARGIN_CASE = 'rwgs_1_1'
MARGIN_PHASE = 'Ti10O19'


def commit():
    try:
        out = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'],
                             cwd=ROOT, capture_output=True, text=True,
                             timeout=20)
        return out.stdout.strip() or 'unknown'
    except Exception:                                    # pragma: no cover
        return 'unknown'


def f(v, n=2):
    return f'{v:.{n}f}'


# ------------------------------------------------------------- defect

def defect_block():
    geom = S.geometry(P)
    b = S.phase_boundaries(P, FLAG_T, geom)
    out = S.distribute(P, FLAG_T, FLAG_VO, geom)
    eps = S.energetics(P)['eps']
    kt = S.KB_EV * (FLAG_T + 273.15)
    th_t = P['surface_phases']['theta_transition']
    th_r = P['surface_phases']['theta_reconstructed_eff']

    shell_u = out['umol_g']['surface'] + out['umol_g']['subsurface']
    shell_cap = geom['N_s'] + geom['N_ss']

    # inventory sweep for the ladder panel: dense enough to show the riser
    hi = max(220.0, b['VO_cs_umol_g'] * 1.32)
    grid = [0.02 + (hi - 0.02) * (i / 239.0) ** 1.6 for i in range(240)]
    for extra in (b['VO_onset_umol_g'], b['VO_recon_complete_umol_g'],
                  b['VO_cs_umol_g'], FLAG_VO):
        grid += [extra * 0.9995, extra, extra * 1.0005]
    grid = sorted(v for v in grid if 0 < v <= hi)
    sweep = []
    for vo in grid:
        d = S.distribute(P, FLAG_T, vo, geom)
        sweep.append({'vo': vo,
                      'surface': d['umol_g']['surface'],
                      'subsurface': d['umol_g']['subsurface'],
                      'bulk': d['umol_g']['bulk'],
                      'extended': d['extended_defects_umol_g'],
                      'regime': d['regime']})

    # grand potential per bridging site: ideal (1x1) branch against the
    # added-row line phase, and the re-crossing that forces the exclusion
    def omega_gas(mu):
        return -kt * math.log(1.0 + math.exp((mu - eps[0]) / kt))

    mu_t = b['mu_t_eV']
    w_t = omega_gas(mu_t)

    def omega_line(mu):
        return w_t - th_r * (mu - mu_t)

    lo, hi_mu = mu_t + 1e-6, 3.0
    for _ in range(200):                       # first re-crossing above mu_t
        mid = 0.5 * (lo + hi_mu)
        if omega_gas(mid) < omega_line(mid):
            hi_mu = mid
        else:
            lo = mid
    mu_recross = 0.5 * (lo + hi_mu)
    mu_axis = [-0.55 + (0.85 + 0.55) * i / 180.0 for i in range(181)]
    omega_curve = [{'mu': m, 'gas': omega_gas(m), 'line': omega_line(m)}
                   for m in mu_axis]

    # shell-thickness ladder, the certified table
    ladder = []
    for n_lay in (1, 2, 3, 4):
        g = S.geometry(P, ss_layers=n_lay)
        d = S.distribute(P, FLAG_T, FLAG_VO, g)
        sh = d['umol_g']['surface'] + d['umol_g']['subsurface']
        ladder.append({'n_lay': n_lay, 'shell_nm': g['shell_nm'],
                       'surface': d['umol_g']['surface'],
                       'subsurface': d['umol_g']['subsurface'],
                       'bulk': d['umol_g']['bulk'],
                       'theta_ss': d['theta']['subsurface'],
                       'mu_V_eV': d['mu_V_eV'],
                       'shell_pct': 100.0 * sh / FLAG_VO})

    # measured BET geometry against the equivalent sphere
    bet_area = 1.60
    g_bet = S.geometry(P, bet_m2_g=bet_area)
    b_bet = S.phase_boundaries(P, FLAG_T, g_bet)
    keys = ('VO_onset_umol_g', 'VO_recon_complete_umol_g', 'VO_cs_umol_g')
    shift = max(abs(b_bet[k] - b[k]) / b[k] for k in keys) * 100.0

    # transport: neutral barrier and the barrier the measurement implies
    e_m0 = P['kinetics_cgmc']['E_m_bulk_eV']
    prof_t = [CGMC_TARGET_S * x for x in (1e-4, 1e-2, 0.2, 1.0)]
    sys0 = CG.build(P, vo_total=FLAG_VO, dist_fractions=out['fractions'])
    run0 = CG.run(sys0, CGMC_TARGET_S, profile_times=prof_t)
    e_fit = CG.effective_barrier(P, CGMC_TARGET_PCT / 100.0, CGMC_TARGET_S,
                                 vo_total=FLAG_VO,
                                 dist_fractions=out['fractions'])
    sys1 = CG.build(P, vo_total=FLAG_VO, dist_fractions=out['fractions'],
                    cg={'E_m_bulk_eV': e_fit})
    run1 = CG.run(sys1, CGMC_TARGET_S, profile_times=prof_t)

    # sampled surface arrangement at the coverage of the surviving (1x1)
    # patches: the ordering exhibit, kept separate from the partition
    th_target = min(out['theta']['surface'], th_t)
    pts = S.isotherm_scan(P, FLAG_T)
    best = min(pts, key=lambda q: abs(q['theta_s'] - th_target))
    ordering = {
        'theta_target': th_target,
        'rows': best['rows'], 'row_sites': best['row_sites'],
        'surf_vac': best['surf_vac'], 'theta_s': best['theta_s'],
        'hist': best['hist'], 'mu_eV': best['mu_eV'],
        'spacing': S.spacing_summary(pts, th_target),
        'pair_eV': P['surface_ordering']['pair_eV']['along_row'],
    }

    def cg_rows(r):
        return [{'t_s': t, 'rec': v}
                for t, v in zip(r['t_s'], r['recovery']) if t > 0]

    plot = {
        'geometry': geom,
        'boundaries': b,
        'flagship': {
            'T_C': FLAG_T, 'VO_total_umol_g': FLAG_VO,
            'mu_V_eV': out['mu_V_eV'], 'regime': out['regime'],
            'umol_g': out['umol_g'], 'theta': out['theta'],
            'fractions': out['fractions'],
            'extended': out['extended_defects_umol_g'],
            'shell_umol_g': shell_u,
            'shell_pct': 100.0 * shell_u / FLAG_VO,
            'shell_depletion': shell_u / shell_cap,
            'surface_phase': out['surface_phase'],
        },
        'sweep': sweep,
        'omega': {'curve': omega_curve, 'mu_t': mu_t, 'omega_t': w_t,
                  'theta_eff': th_r, 'theta_t': th_t,
                  'mu_recross': mu_recross, 'mu_cs': b['mu_cs_eV']},
        'shell_ladder': ladder,
        'ordering': ordering,
        'bet': {'area_m2_g': bet_area, 'boundaries': b_bet,
                'max_shift_pct': shift},
        'cgmc': {
            'T_C': P['co2_kinetics']['defaults']['T_reox_C'],
            'E_m_eV': e_m0, 'E_m_fit_eV': e_fit,
            'target_pct': CGMC_TARGET_PCT, 'target_s': CGMC_TARGET_S,
            'neutral': {'rows': cg_rows(run0),
                        'profiles': run0['profiles'],
                        'recovered': run0['recovered_fraction']},
            'fitted': {'rows': cg_rows(run1),
                       'profiles': run1['profiles'],
                       'recovered': run1['recovered_fraction']},
            'depth_nm': run0['depth_nm'], 'R_nm': sys0['R_nm'],
        },
    }

    slots = {
        'defect.geometry.d_um': f(P['defaults']['particle_diameter_um'], 1),
        'defect.geometry.area_m2_g': f(geom['area_m2_g']),
        'defect.geometry.N_s': f(geom['N_s']),
        'defect.geometry.N_ss': f(geom['N_ss']),
        'defect.geometry.layer_nm': f(geom['layer_nm'], 3),
        'defect.geometry.surface_pct_of_O':
            f(100.0 * geom['N_s'] / geom['N_O_total'], 3),
        'defect.flagship.T_C': f'{FLAG_T:g}',
        'defect.flagship.VO_total_umol_g': f'{FLAG_VO:g}',
        'defect.flagship.mu_V_eV': f(out['mu_V_eV'], 3),
        'defect.flagship.umol.surface': f(out['umol_g']['surface']),
        'defect.flagship.umol.subsurface': f(out['umol_g']['subsurface'], 1),
        'defect.flagship.umol.bulk': f(out['umol_g']['bulk'], 1),
        'defect.flagship.shell_pct': f(100.0 * shell_u / FLAG_VO, 1),
        'defect.flagship.shell_depletion_pct':
            f(100.0 * shell_u / shell_cap, 0),
        'defect.boundaries.mu_t_eV': f(mu_t, 3),
        'defect.boundaries.VO_onset_umol_g': f(b['VO_onset_umol_g']),
        'defect.boundaries.VO_recon_complete_umol_g':
            f(b['VO_recon_complete_umol_g']),
        'defect.boundaries.VO_cs_umol_g': f(b['VO_cs_umol_g'], 0),
        'defect.omega.theta_eff': f(th_r, 1),
        'defect.omega.mu_recross_eV': f(mu_recross, 2),
        'defect.shell_ladder_text.subsurface':
            ' / '.join(f(r['subsurface'], 1) for r in ladder),
        'defect.shell_ladder_text.bulk':
            ' / '.join(f(r['bulk'], 1) for r in ladder),
        'defect.shell_ladder_text.first_pct': f(ladder[0]['shell_pct'], 0),
        'defect.shell_ladder_text.last_pct': f(ladder[-1]['shell_pct'], 0),
        'defect.bet.area_m2_g': f(bet_area),
        'defect.bet.max_shift_pct': f(shift, 1),
        'defect.cgmc.T_C': f'{P["co2_kinetics"]["defaults"]["T_reox_C"]:g}',
        'defect.cgmc.E_m_eV': f(e_m0),
        'defect.cgmc.E_m_fit_eV': f(e_fit),
        'defect.cgmc.target_pct': f'{CGMC_TARGET_PCT:g}',
        'defect.cgmc.target_s': f'{CGMC_TARGET_S:g}',
    }
    return plot, slots


# -------------------------------------------------------- equilibrium

def equilibrium_block():
    rows = REF_TH['rows']
    cond = REF_TH['conditions']

    def case_rows(case):
        return sorted((r for r in rows if r['case'] == case),
                      key=lambda r: r['T_C'])

    margin = []
    for r in case_rows(MARGIN_CASE):
        rc = r['inactive_phase_reduced_costs'].get(MARGIN_PHASE)
        margin.append({'T_C': r['T_C'],
                       'r_kJ': rc['per_mol_O_kJ'] if rc else None,
                       'active': r['active_condensed_phases']})

    feeds = []
    for case, label in (('rwgs_1_1', 'CO2/H2 = 1:1'),
                        ('pure_h2', 'pure H2'),
                        ('pure_n2', 'sealed N2')):
        pts = [{'T_C': r['T_C'], 'reduced_pct': r['reduced_pct'],
                'active': r['active_condensed_phases']}
               for r in case_rows(case)]
        feeds.append({'case': case, 'label': label, 'points': pts})

    h2_last = case_rows('pure_h2')[-1]
    h2_phase = next((s for s in h2_last['active_condensed_phases']
                     if s != 'TiO2'), 'none')

    # candidate ladder for the method plate: every enumerated active set at
    # one temperature, with the winner marked
    method_T = 900
    row = next(r for r in case_rows(MARGIN_CASE) if r['T_C'] == method_T)
    costs = sorted(
        ((s, rc['per_mol_O_kJ'])
         for s, rc in row['inactive_phase_reduced_costs'].items()
         if rc['per_mol_O_kJ'] is not None),
        key=lambda kv: kv[1])
    n_solids = len(row['solid_moles'])

    plot = {
        'margin': {'case': MARGIN_CASE, 'phase': MARGIN_PHASE,
                   'label': REF_TH['cases'][MARGIN_CASE], 'rows': margin},
        'feeds': feeds,
        'method': {'T_C': method_T, 'n_solids': n_solids,
                   'active': row['active_condensed_phases'],
                   'costs': [{'phase': s, 'r_kJ': v} for s, v in costs]},
        'conditions': cond,
    }

    first, last = margin[0], margin[-1]
    slots = {
        'equilibrium.margin.phase': MARGIN_PHASE,
        'equilibrium.margin.feed_label': 'CO2/H2 = 1:1',
        'equilibrium.margin.first_kJ': f(first['r_kJ']),
        'equilibrium.margin.last_kJ': f(last['r_kJ']),
        'equilibrium.margin.first_T_C': f'{first["T_C"]:g}',
        'equilibrium.margin.last_T_C': f'{last["T_C"]:g}',
        'equilibrium.feeds.h2_phase': f'TiO2/{h2_phase}',
        'equilibrium.conditions.mass_mg':
            f'{cond.get("mass_g", 0.1) * 1000:g}',
        'equilibrium.conditions.volume_mL': f'{cond.get("V_cm3", 25.0):g}',
        'equilibrium.conditions.pressure_atm':
            f'{cond.get("P_atm", 1.0):g}',
        'equilibrium.method.n_solids': f'{n_solids:d}',
        'equilibrium.method.T_C': f'{method_T:g}',
    }
    return plot, slots


# -------------------------------------------------------- verification

def verification_block():
    """Deviations measured here, live, against the committed references."""
    geom = S.geometry(P)
    layers = []

    worst = 0.0
    for ref in REF_SM['analytic_grid'] + [REF_SM['flagship']]:
        row = S.distribute(P, ref['T_C'], ref['VO_total_umol_g'], geom,
                           preset=ref['preset'])
        worst = max(worst, abs(row['mu_V_eV'] - ref['mu_V_eV']))
    layers.append({'layer': 'partition vs 50-digit oracle',
                   'measured': worst, 'tolerance': 1e-9,
                   'unit': 'eV in mu_V'})

    worst_sh = 0.0
    for ref in REF_SM['shell_ladder']:
        g = S.geometry(P, ss_layers=ref['ss_layers'])
        row = S.distribute(P, ref['T_C'], ref['VO_total_umol_g'], g,
                           preset=ref['preset'])
        for k in ('surface', 'subsurface', 'bulk'):
            worst_sh = max(worst_sh, abs(row['umol_g'][k]
                                         - ref['umol_g'][k]))
    layers.append({'layer': 'shell ladder vs oracle', 'measured': worst_sh,
                   'tolerance': 1e-9, 'unit': 'umol-O/g'})

    ideal = REF_SM['ideal_ensemble']
    kt = S.KB_EV * (ideal['T_C'] + 273.15)
    lo, hi = -40.0, 40.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        n = sum(m * S.fd(mid, e, kt)
                for m, e in zip(ideal['m_by_class'], ideal['eps_eV']))
        if n < ideal['n_vac']:
            lo = mid
        else:
            hi = mid
    layers.append({'layer': 'finite-lattice mu vs exact ensemble',
                   'measured': abs(0.5 * (lo + hi) - ideal['mu_exact_eV']),
                   'tolerance': 5e-3, 'unit': 'eV (finite-size)'})

    worst_cg = 0.0
    spec = REF_SM['cgmc_linear']['spec']
    sysl = {'m': 3, 'q': [float(x) for x in spec['q']], 'eps_c': spec['eps'],
            'lam': spec['lam'], 'beta': 1.0 / (S.KB_EV * (spec['T_C']
                                                          + 273.15)),
            'k_fill': spec['k_fill'], 'eta0': list(spec['eta0']),
            'V0': sum(spec['eta0']), 'R_nm': 3.0, 'a_nm': 1.0,
            'centers': [2.5, 1.5, 0.5], 'T_reox_C': spec['T_C'],
            'E_m_eV': 0.0, 'eps_ctrl': 0.02, 'column_counts': 1e5}
    for ref in REF_SM['cgmc_linear']['rows']:
        r = CG.run(dict(sysl), ref['t_s'], n_out=8)
        for got, want in zip(r['final_eta'], ref['eta']):
            worst_cg = max(worst_cg, abs(got - want) / max(abs(want), 1e-30))
    layers.append({'layer': 'coarse transport vs matrix exponential',
                   'measured': worst_cg, 'tolerance': 1e-2,
                   'unit': 'relative'})

    n_tests = collected_tests()
    thermo_dps = int(re.search(r'(\d+)\s+certified',
                               REF_TH['precision']).group(1))
    plot = {'layers': layers, 'n_tests': n_tests,
            'statmech_dps': REF_SM['precision_dps'],
            'thermo_dps': thermo_dps}
    slots = {
        'verification.statmech_dps': f'{REF_SM["precision_dps"]:g}',
        'verification.thermo_dps': str(thermo_dps),
        'verification.n_tests': str(n_tests),
        'verification.trace_floor_log10': '-41',
    }
    return plot, slots


def collected_tests():
    """How many tests the suite actually has - never a typed number."""
    out = subprocess.run([sys.executable, '-m', 'pytest', '--collect-only',
                          '-q', 'tests/'], cwd=ROOT, capture_output=True,
                         text=True, timeout=900)
    for line in reversed(out.stdout.strip().splitlines()):
        line = line.strip()
        if line.endswith('tests collected') or ' tests collected' in line:
            return int(line.split()[0])
        if line.endswith('test collected'):
            return int(line.split()[0])
    raise RuntimeError('could not read the collected test count:\n'
                       + out.stdout[-2000:])


def main():
    d_plot, d_slots = defect_block()
    e_plot, e_slots = equilibrium_block()
    v_plot, v_slots = verification_block()
    doc = {
        'generated_by': 'scripts/export_plates.py',
        'commit': commit(),
        'spec': REG.SPEC,
        'roles': [{'key': k, 'hue': h, 'where': w} for k, h, w in REG.ROLES],
        'plot': {'defect': d_plot, 'equilibrium': e_plot,
                 'verification': v_plot},
        'slots': dict(**d_slots, **e_slots, **v_slots),
    }
    with (DATA / 'figure_data.json').open('w') as fh:
        json.dump(doc, fh, indent=1)
    print(f'figure_data.json written at {doc["commit"]}: '
          f'{len(doc["slots"])} slots, '
          f'{len(d_plot["sweep"])} sweep points, '
          f'E_m,eff = {d_plot["cgmc"]["E_m_fit_eV"]:.2f} eV')


if __name__ == '__main__':
    main()
