"""Python reference for the browser solver: RWGS over Ti-O and over Ce-O.

The page runs the same calculation in JavaScript. This writes the answer the
port has to reproduce, and `tests/test_oxide_port.py` compares them to five
significant figures.

Both systems are full-composition Gibbs minimisation in reaction-extent
coordinates - the same engine, only the condensed phases differ.
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / 'data'
sys.path.insert(0, str(ROOT))

import numpy as np

from solidgas.gibbs import GibbsSystem
from solidgas.equilibrium import ACTIVE_TI_PHASES, moles_of_gas
from solidgas.shomate import R_KJ, mu0

P_TOT = 1.0
V_CM3 = 25.0
MASS_G = 0.100
TEMPS_C = [500, 700, 900, 1100, 1300, 1500]

GASES = ['CO2', 'H2', 'CO', 'H2O', 'CH4', 'O2']

SYSTEMS = {
    'Ti-O': dict(metal='Ti', host='TiO2', mw=79.87, solids=ACTIVE_TI_PHASES,
                 source='waldner',
                 note='Waldner & Eriksson, Calphad 23 (1999) 189-218'),
    'Ce-O': dict(metal='Ce', host='CeO2', mw=172.12, solids=['CeO2', 'Ce2O3'],
                 source='nist',
                 note='NIST-JANAF, one coefficient set across the range'),
}


def build(name):
    cfg = SYSTEMS[name]
    return GibbsSystem(GASES, cfg['solids'],
                       ['CO2', 'H2', 'H2O', cfg['host']], source=cfg['source']), cfg


def seeds_for(S, ng, n_metal):
    """Extent starts: a ladder of gas conversion crossed with each phase basin."""
    free = S.free
    ref = np.array([ng if S.is_gas[i] else n_metal for i in free])
    gas_pos = [j for j, i in enumerate(free) if S.is_gas[i]]
    sol_pos = [j for j, i in enumerate(free) if not S.is_gas[i]]
    out = []
    for g in (1e-6, 0.05, 0.2, 0.4, 0.7):
        base = np.full(len(free), 1e-48)
        for j in gas_pos:
            base[j] = ref[j] * g
        out.append(base.copy())
        for j in sol_pos:
            for m in (1e-1, 1e-4, 1e-8, 1e-14, 1e-22, 1e-32, 1e-42):
                v = base.copy()
                v[j] = ref[j] * m
                out.append(v)
    return out


def kp_rwgs(T):
    return float(np.exp(-(mu0('CO', T) + mu0('H2O', T)
                          - mu0('CO2', T) - mu0('H2', T)) / (R_KJ * T)))


def run(name):
    S, cfg = build(name)
    n_metal = MASS_G / cfg['mw']
    rows = []
    for T_C in TEMPS_C:
        T = T_C + 273.15
        ng = moles_of_gas(P_TOT, V_CM3, T)
        n0 = ng / 2
        b = S.b_from(**{'CO2': n0, 'H2': n0, cfg['host']: n_metal})
        r = S.minimise(b, T, seeds_for(S, ng, n_metal), P_TOT)
        n = r['n']
        y = S.gas_fractions(n)
        Q = ((y['CO'] * y['H2O']) / (y['CO2'] * y['H2'])
             if y['CO2'] > 0 and y['H2'] > 0 else float('nan'))
        rows.append({
            'system': name, 'T_C': T_C, 'T_K': round(T, 2),
            'method': 'gibbs_min', 'formulation': r['formulation'],
            'reduced_pct': S.reduced_percent(n, n_metal),
            'phase_split_pct': {p: v for p, v in S.phase_split(n, n_metal).items()
                                if v > 1e-30},
            'gas_pct': {k: v * 100.0 for k, v in y.items()},
            'conversion_CO2_pct': (n0 - n[S.species.index('CO2')]) / n0 * 100.0,
            'Kp_rwgs': kp_rwgs(T), 'Q_rwgs': Q,
            'G_rel_kJ': r['G_rel_kJ'],
            'balance_residual_mol': r['residual'],
            'n_gas_mol': ng, 'n_metal_mol': n_metal,
        })
    return rows


def main():
    doc = {'conditions': {'P_atm': P_TOT, 'V_cm3': V_CM3, 'mass_g': MASS_G},
           'method': 'gibbs_min', 'temperatures_C': TEMPS_C,
           'systems': {k: {kk: vv for kk, vv in v.items() if kk != 'solids'}
                       for k, v in SYSTEMS.items()},
           'rows': []}
    for name in SYSTEMS:
        rows = run(name)
        doc['rows'] += rows
        print(f"--- {name} ---")
        for r in rows:
            ph = ' + '.join(f'{p} {v:.3g}%' for p, v in r['phase_split_pct'].items())
            print(f"  {r['T_C']:5d}C  red={r['reduced_pct']:10.6f}%  conv={r['conversion_CO2_pct']:6.2f}%  "
                  f"Kp={r['Kp_rwgs']:7.4f} Q={r['Q_rwgs']:7.4f}  resid={r['balance_residual_mol']:.1e}  {ph}")
    with (DATA / 'oxide_reference.json').open('w') as fh:
        json.dump(doc, fh, indent=1)
    print(f"\n{len(doc['rows'])} rows -> oxide_reference.json")
    return 0


if __name__ == '__main__':
    sys.exit(main())
