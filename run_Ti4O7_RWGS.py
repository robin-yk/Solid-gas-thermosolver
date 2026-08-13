"""RWGS over TiO2/Ti4O7: CO2/H2 = 1:1 across 500-1500 C.

Asks whether a stoichiometric RWGS feed can pull rutile TiO2 down to the
Magneli phase Ti4O7. It cannot -- the CO/CO2 and H2O/H2 ratios the gas
settles at are nowhere near reducing enough. See README.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from solidgas import (ELEM, gas_idx, mu0, R_KJ, moles_of_gas, solve,
                      ti3_percent, gas_fractions)

matplotlib.rcParams.update({'figure.dpi': 120, 'font.size': 11})

P_tot = 1.0       # atm
V = 25.0          # cm3
m = 0.100         # g TiO2
MW = 79.87        # g/mol TiO2
n0_TiO2 = m / MW

T_C_arr = np.linspace(500, 1500, 80)
results = []
for T_C in T_C_arr:
    T = T_C + 273.15
    n0g = moles_of_gas(P_tot, V, T)
    n0 = n0g / 2                      # equimolar CO2 and H2
    #           C     H          O                Ti
    b = np.array([n0, 2 * n0, 2 * n0 + 2 * n0_TiO2, n0_TiO2])

    inits = [np.array([n0 * fCO2, n0 * fH2, n0 * fCO, n0 * fH2O,
                       n0 * 1e-6, n0_TiO2 * 0.99, n0_TiO2 * 0.0025])
             for fCO2, fH2, fCO, fH2O in [(0.5, 0.5, 0.01, 0.01),
                                          (0.1, 0.1, 0.4, 0.4),
                                          (0.01, 0.01, 0.48, 0.48)]]
    n = solve(b, T, inits, P_tot)
    if n is None:
        continue

    y = gas_fractions(n)[:5]
    conv = (n0 - n[0]) / n0 * 100 if n0 > 0 else 0
    Ti3 = ti3_percent(n, n0_TiO2)
    Kp_rw = np.exp(-(mu0('CO', T) + mu0('H2O', T)
                     - mu0('CO2', T) - mu0('H2', T)) / (R_KJ * T))
    Q_rw = (y[2] * y[3]) / (y[0] * y[1]) if y[0] > 0 and y[1] > 0 else float('nan')
    results.append(dict(T_C=T_C, y_CO2=y[0], y_H2=y[1], y_CO=y[2], y_H2O=y[3], y_CH4=y[4],
                        Ti3=Ti3, conv=conv, Kp_rw=Kp_rw, Q_rw=Q_rw,
                        n_Ti4O7=n[6], n_TiO2=n[5]))

# -- table --
print(f"n0_TiO2 = {n0_TiO2:.4e} mol (100 mg / 79.87)")
print(f"{'T(C)':>6} {'CO2%':>6} {'H2%':>6} {'CO%':>6} {'H2O%':>6} {'CH4%':>7} "
      f"{'conv%':>6} {'Ti3+%':>7} {'Kp':>6} {'Q':>6}")
for r in results:
    if r['T_C'] % 100 < 13:
        print(f"{r['T_C']:6.0f} {r['y_CO2']:6.2f} {r['y_H2']:6.2f} {r['y_CO']:6.2f} {r['y_H2O']:6.2f} "
              f"{r['y_CH4']:7.3f} {r['conv']:6.1f} {r['Ti3']:7.3f} {r['Kp_rw']:6.3f} {r['Q_rw']:6.3f}")

print(f"\nTi3+ max over range = {max(r['Ti3'] for r in results):.4f}%")

# -- figure --
T_p = np.array([r['T_C'] for r in results])
fig, axes = plt.subplots(2, 1, figsize=(8, 8))
fig.suptitle('CO2/H2=1:1 + TiO2/Ti4O7 (RWGS) - Gibbs Minimization\n'
             '100 mg TiO2 | 25 mL | 1 atm | NIST-JANAF', fontsize=11)
ax2 = axes[0]
for key, col, lab in [('y_CO2', 'purple', 'CO2'), ('y_H2', 'steelblue', 'H2'),
                      ('y_CO', 'tomato', 'CO'), ('y_H2O', 'orange', 'H2O'),
                      ('y_CH4', 'brown', 'CH4')]:
    ax2.plot(T_p, [r[key] for r in results], lw=2, label=lab, color=col)
ax2.set_ylabel('Gas mol%'); ax2.set_xlabel('T (C)')
ax2.legend(fontsize=9, ncol=2); ax2.grid(True, alpha=0.3)
ax3 = axes[1]
Ti3p = np.array([r['Ti3'] for r in results])
ax3.plot(T_p, 100 - Ti3p, 'darkgreen', lw=2.5, label='Ti4+ %')
ax3.plot(T_p, Ti3p, 'darkorange', lw=2.5, label='Ti3+ % (via Ti4O7)')
ax3.set_ylim([-2, 102]); ax3.legend(fontsize=9)
ax3.set_ylabel('Ti oxidation state (%)'); ax3.set_xlabel('T (C)')
ax3.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('Ti4O7_RWGS.png', dpi=120)
print("\nfigure saved -> Ti4O7_RWGS.png")
