"""Every computed number in the manuscript, in one command.

    python3 scripts/reproduce_paper.py

Writes paper_outputs/*.csv and paper_outputs/site_data.json. The CSVs are
what the manuscript and its supplement quote; the JSON is what the web
page displays, so the page cannot show a number the package did not
produce.

Two calculations run here and nothing else:

  gas-solid equilibrium over the C-H-O-N/Ti-O system, for the reaction
  feed the paper uses, across the temperature range it covers;

  the isolated-vacancy population fitted to the reduction series: the
  effective partition of each sample's oxygen removal and the CO rate
  it predicts.

The dielectric regression is carried through verbatim from the fit
published with the paper; it is a measurement analysis and no model here
feeds it or reads it.
"""

import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from solidgas import activeset as A
from solidgas import vacancy_population as VP
from analysis import dielectric_correlation as D

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, 'paper_outputs')

# The reaction feed of the manuscript: H2/CO2/He = 10/10/80 sccm.
FEED = {'CO2': 10, 'H2': 10, 'N2': 80}
FEED_LABEL = 'CO2/H2/He = 10/10/80'
T_SWEEP = list(range(400, 1501, 25))
T_QUOTED = [500, 600, 700, 800, 900, 1000, 1100, 1300, 1500]



def _rows_equilibrium():
    rows = []
    for T_C in T_SWEEP:
        r = A.solve(FEED, T_C + 273.15)
        costs = {k: v['per_mol_O_kJ']
                 for k, v in r['inactive_phase_reduced_costs'].items()}
        first = min(costs, key=costs.get) if costs else ''
        rows.append({
            'T_C': float(T_C),
            'X_CO2_pct': r['conversion_CO2_pct'],
            'y_CO2': r['gas_fractions']['CO2'],
            'y_H2': r['gas_fractions']['H2'],
            'y_CO': r['gas_fractions']['CO'],
            'y_H2O': r['gas_fractions']['H2O'],
            'stable_phases': '+'.join(r['active_condensed_phases']),
            'nearest_reduced_phase': first,
            'margin_kJ_per_mol_O': costs[first] if first else '',
        })
    return rows


def _rows_phase_stability():
    rows = []
    for T_C in T_SWEEP:
        b = A.reduction_boundary(T_C + 273.15)
        rows.append({
            'T_C': T_C,
            'critical_y_CO2_mol_pct': b['y_CO2'] * 100.0,
            'mu_O_crit_kJ_per_mol': b['mu_crit_kJ_per_mol_O'],
            'first_phase': b.get('phase', ''),
        })
    return rows


def _rows_dielectric():
    s = D.stored()
    return [{
        'form': s['form'],
        'log10_intercept': s['log10_intercept'],
        'log10_intercept_err': s['log10_intercept_err'],
        'exponent': s['exponent'],
        'exponent_err': s['exponent_err'],
        'r_squared': s['fit_stats']['r_squared'],
        'n_points': s['fit_stats']['n_points'],
        'fit_range_umol_g': '-'.join(str(v) for v in s['fit_range_umol_g']),
        'source': 'fit published with the manuscript; not recomputed here',
    }]


def _write(name, rows):
    path = os.path.join(OUT, name)
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    return path, len(rows)


def main():
    os.makedirs(OUT, exist_ok=True)
    pop_series = VP.load_series()
    R = VP.REPORTED
    pop_rows = VP.partition(pop_series, R['ns'], R['E_loss'],
                            10.0 ** R['log10_nu'])
    pop_csv = [{
        'sample': q['sample'], 'treatment': q['treatment'],
        'T_red_C': q['T_red_C'], 't_red_s': q['t_red_s'],
        'nV_total_umol_g': q['nV_total_umol_g'],
        'r_CO_measured': q['r_CO_measured'], 'r_CO_fitted': q['r_CO_fitted'],
        'n_iso_umol_g': q['n_iso_umol_g'], 'n_assoc_umol_g': q['n_assoc_umol_g'],
        'n_below_umol_g': q['n_below_umol_g'],
        'below_fraction': q['below_fraction'],
    } for q in pop_rows]

    eq = _rows_equilibrium()
    ps = _rows_phase_stability()
    di = _rows_dielectric()

    for name, rows in (('equilibrium_conversion.csv', eq),
                       ('phase_stability.csv', ps),
                       ('vacancy_population.csv', pop_csv),
                       ('dielectric_correlation.csv', di)):
        path, n = _write(name, rows)
        print('%-34s %4d rows' % (os.path.basename(path), n))

    # The figure code reads solver records, so the page is handed the
    # records themselves rather than a reshaped copy of them - one less
    # place for the page and the package to disagree.
    fig_rows, boundary = [], []
    for T_C in T_SWEEP:
        r = A.solve(FEED, T_C + 273.15)
        fig_rows.append({
            # T_K - 273.15 is not exact in binary; the page keys off this
            # value, so it carries the temperature that was asked for.
            'T_C': float(T_C),
            'conversion_CO2_pct': r['conversion_CO2_pct'],
            'active_condensed_phases': r['active_condensed_phases'],
            'gas_fractions': r['gas_fractions'],
            'solver_tolerance': r['solver_tolerance'],
            'inactive_phase_reduced_costs': {
                k: {'per_formula_kJ': float(v['per_formula_kJ']),
                    'per_mol_O_kJ': (None if v['per_mol_O_kJ'] is None
                                     else float(v['per_mol_O_kJ'])),
                    'oxygen_removed_per_extent': v['oxygen_removed_per_extent']}
                for k, v in r['inactive_phase_reduced_costs'].items()},
        })
        b = A.reduction_boundary(T_C + 273.15)
        boundary.append({'T_C': T_C, 'y_CO2': b['y_CO2'], 'phase': b.get('phase')})

    # The supplement does not record which residual it minimised and the
    # rates span 86-fold, so every one is reported.
    pop_fit = {k: VP.fit(pop_series, kind=k, iters=250) for k in VP.RESIDUALS}

    site = {
        'generated_by': 'scripts/reproduce_paper.py',
        'feed': {'label': FEED_LABEL, 'sccm': FEED},
        'equilibrium': {
            'rows': fig_rows,
            'boundary': boundary,
            'T_C': [r['T_C'] for r in eq],
            'X_CO2_pct': [r['X_CO2_pct'] for r in eq],
            'margin_kJ_per_mol_O': [r['margin_kJ_per_mol_O'] for r in eq],
            'stable_phases': [r['stable_phases'] for r in eq],
            'nearest_reduced_phase': [r['nearest_reduced_phase'] for r in eq],
            'critical_y_CO2_mol_pct': [r['critical_y_CO2_mol_pct'] for r in ps],
            'quoted': {str(T): A.solve(FEED, T + 273.15)['conversion_CO2_pct']
                       for T in T_QUOTED},
        },
        'population': {
            'series': VP.load_series(),
            'ns_bounds': list(VP.NS_BOUNDS),
            'reported': {'ns_umol_g': VP.REPORTED['ns'],
                         'E_loss_eV': VP.REPORTED['E_loss'],
                         'log10_nu_eff': VP.REPORTED['log10_nu'],
                         'nu_eff_s1': 10.0 ** VP.REPORTED['log10_nu']},
            'residual_kinds': list(VP.RESIDUALS),
            'rows': pop_rows,
            'fit': pop_fit,
        },
        'dielectric': D.stored(),
        'scope': ('This repository computes gas-solid equilibrium and an '
                  'isolated-vacancy population fitted to the reduction series. '
                  'The population model infers n_iso from the measured rates; '
                  'it does not predict them independently.'),
    }
    path = os.path.join(OUT, 'site_data.json')
    with open(path, 'w') as fh:
        json.dump(site, fh, indent=1, sort_keys=True)
        fh.write('\n')
    print('%-34s %4d kB' % ('site_data.json', os.path.getsize(path) // 1024))


if __name__ == '__main__':
    main()
