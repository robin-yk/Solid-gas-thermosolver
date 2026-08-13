"""Strongly reducing feeds over TiO2/Ti4O7: pure H2 and pure CH4, 500-1500 C.

The RWGS run leaves TiO2 untouched, so this sweeps the two reducing extremes
to find what it actually takes to open the Magneli branch. CH4 does it well
below H2 because carbon removes oxygen as CO rather than H2O.
"""

import numpy as np

from solidgas import gas_idx, moles_of_gas, solve, ti3_percent

P_tot = 1.0       # atm
V = 25.0          # cm3
m = 0.100         # g TiO2
MW = 79.87        # g/mol TiO2
n0_TiO2 = m / MW


def run(feed):
    out = []
    for T_C in np.linspace(500, 1500, 80):
        T = T_C + 273.15
        n0g = moles_of_gas(P_tot, V, T)
        if feed == 'H2':
            #             C       H            O          Ti
            b = np.array([0.0, 2 * n0g, 2 * n0_TiO2, n0_TiO2])
        else:
            b = np.array([n0g, 4 * n0g, 2 * n0_TiO2, n0_TiO2])

        # Seed across a spread of reduction extents so the solver sees both the
        # untouched-TiO2 basin and the Ti4O7 basin.
        inits = []
        for fr in [0.001, 0.05, 0.2, 0.5]:
            nT4 = n0_TiO2 * fr / 4 * 0.999
            nTi = n0_TiO2 - 4 * nT4
            if feed == 'H2':
                inits.append(np.array([1e-12, 2 * n0g * (1 - fr), 1e-12,
                                       2 * n0g * fr, 1e-12, nTi, nT4]))
            else:
                inits.append(np.array([n0g * fr * 0.3, n0g * fr, n0g * fr,
                                       n0g * fr * 0.1, n0g * (1 - fr), nTi, nT4]))

        n = solve(b, T, inits, P_tot)
        if n is None:
            continue
        ng = sum(n[i] for i in gas_idx)
        out.append(dict(T_C=T_C, Ti3=ti3_percent(n, n0_TiO2), nT4=n[6],
                        y_H2=n[1] / ng * 100, y_H2O=n[3] / ng * 100,
                        y_CO=n[2] / ng * 100, y_CH4=n[4] / ng * 100))
    return out


for feed in ['H2', 'CH4']:
    res = run(feed)
    print(f"\n===== pure {feed} feed | TiO2/Ti4O7 =====")
    print(f"{'T(C)':>6} {'Ti3+%':>8} {'H2%':>7} {'H2O%':>7} {'CO%':>7} {'CH4%':>7}")
    for r in res:
        if r['T_C'] % 100 < 13:
            print(f"{r['T_C']:6.0f} {r['Ti3']:8.4f} {r['y_H2']:7.2f} {r['y_H2O']:7.2f} "
                  f"{r['y_CO']:7.2f} {r['y_CH4']:7.2f}")
    print(f"  -> Ti3+ max = {max(r['Ti3'] for r in res):.4f}% "
          f"at {max(res, key=lambda r: r['Ti3'])['T_C']:.0f} C")
