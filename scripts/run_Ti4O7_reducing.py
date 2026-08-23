"""Strongly reducing feeds over the TiO2-Ti4O7-Ti3O5-Ti2O3 series, 500-1500 C.

Pure H2 and pure CH4, the two reducing extremes. CH4 opens the Magneli branch
roughly 500 C below H2 because carbon removes oxygen as CO rather than H2O.

Neither feed gets past Ti4O7: the solid settles on the TiO2/Ti4O7 two-phase
buffer and the gas is pinned there, an order of magnitude too oxidising to
carry the next step. See README.
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / 'data'
sys.path.insert(0, str(ROOT))

import numpy as np

from solidgas import (ACTIVE_TI_PHASES, gas_idx, mean_valence, moles_of_gas,
                      phase_seeds, phase_split, solve, ti3_percent)

P_tot = 1.0       # atm
V = 25.0          # cm3
m = 0.100         # g TiO2
MW = 79.87        # g/mol TiO2
n0_TiO2 = m / MW


def seeds(n0g, feed, previous=None, cold=True):
    """Starting points for one temperature.

    `previous` is the neighbouring temperature's answer, carried forward as a
    warm start; along a smooth branch it lands almost on the solution. The full
    cold sweep over every solid basin is what catches a branch change, so it is
    run periodically rather than everywhere - it costs about seven times as
    much, and the answer only moves where a phase boundary is crossed.
    """
    def gas_at(gf):
        if feed == 'H2':
            return [1e-12, 2 * n0g * (1 - gf), 1e-12, 2 * n0g * gf, 1e-12]
        return [n0g * gf * 0.3, n0g * gf, n0g * gf, n0g * gf * 0.1, n0g * (1 - gf)]

    out = [] if previous is None else [previous.copy()]
    solids = phase_seeds(n0_TiO2)
    if cold:
        for sv in solids:
            for gf in (0.05, 0.4, 0.9):
                out.append(np.array(gas_at(gf) + list(sv)))
    else:
        for gf in (0.05, 0.4, 0.9):
            out.append(np.array(gas_at(gf) + list(solids[0])))
        for sv in solids[1:]:           # still probe each reduced basin once
            out.append(np.array(gas_at(0.4) + list(sv)))
    return out


def run(feed):
    out, previous = [], None
    for step, T_C in enumerate(np.linspace(500, 1500, 40)):
        T = T_C + 273.15
        n0g = moles_of_gas(P_tot, V, T)
        if feed == 'H2':
            #             C       H            O          Ti
            b = np.array([0.0, 2 * n0g, 2 * n0_TiO2, n0_TiO2])
        else:
            b = np.array([n0g, 4 * n0g, 2 * n0_TiO2, n0_TiO2])

        # full cold restart every fifth point, warm start in between
        n = solve(b, T, seeds(n0g, feed, previous, cold=(step % 5 == 0)), P_tot)
        if n is None:
            continue
        previous = n
        ng = sum(n[i] for i in gas_idx)
        rec = dict(T_C=T_C, Ti3=ti3_percent(n, n0_TiO2),
                   valence=mean_valence(n, n0_TiO2),
                   y_H2=n[1] / ng * 100, y_H2O=n[3] / ng * 100,
                   y_CO=n[2] / ng * 100, y_CH4=n[4] / ng * 100)
        rec.update(phase_split(n, n0_TiO2))
        out.append(rec)
    return out


all_results = {}
for feed in ['H2', 'CH4']:
    res = run(feed)
    all_results[feed] = res
    print(f"\n===== pure {feed} feed | {' / '.join(ACTIVE_TI_PHASES)} =====")
    head = (f"{'T(C)':>6} {'Ti3+%':>8} {'valence':>8} {'H2%':>7} {'H2O%':>7} "
            f"{'CO%':>7} {'CH4%':>7}") + ''.join(f'{p:>8}' for p in ACTIVE_TI_PHASES)
    print(head)
    for r in res:
        if r['T_C'] % 100 < 26:
            print(f"{r['T_C']:6.0f} {r['Ti3']:8.4f} {r['valence']:8.3f} {r['y_H2']:7.2f} "
                  f"{r['y_H2O']:7.2f} {r['y_CO']:7.2f} {r['y_CH4']:7.2f}"
                  + ''.join(f"{r[p]:8.2f}" for p in ACTIVE_TI_PHASES))
    best = max(res, key=lambda r: r['Ti3'])
    print(f"  -> Ti3+ max = {best['Ti3']:.4f}% at {best['T_C']:.0f} C")
    seen = [p for p in ACTIVE_TI_PHASES if any(r[p] > 0.01 for r in res)]
    print(f"  -> phases that actually appear: {', '.join(seen)}")

with (DATA / 'results_reducing.json').open('w') as fh:
    json.dump({'conditions': {'P_atm': P_tot, 'V_cm3': V, 'm_TiO2_g': m,
                              'n0_TiO2_mol': n0_TiO2, 'phases': ACTIVE_TI_PHASES},
               'feeds': all_results}, fh, indent=1)
print("\ndata saved -> data/results_reducing.json")
