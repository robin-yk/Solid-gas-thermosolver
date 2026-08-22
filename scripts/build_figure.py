"""Verification figure for the Ti-O active-set solver -> verification_figure.png.

Panels A-C are computed live by the same production/oracle code the tests
run; panel D's ledger values are the worst deviations measured this session
(the gates themselves are enforced in tests/, not here).
"""
import os
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(SCRIPT_DIR))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from mpmath import mp, mpf, log as mlog
import oracle_tio as O
from solidgas import activeset as A

# palette / chrome (dataviz reference instance, light mode)
BLUE, ORANGE = '#2a78d6', '#eb6834'
SEQ100 = '#cde2fb'
SURF, INK, INK2, MUTED = '#fcfcfb', '#0b0b0b', '#52514e', '#898781'
GRID, AXIS = '#e1e0d9', '#c3c2b7'

plt.rcParams.update({
    'figure.facecolor': SURF, 'axes.facecolor': SURF,
    'font.family': 'DejaVu Sans', 'font.size': 9.5,
    'axes.edgecolor': AXIS, 'axes.linewidth': 1.0,
    'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.8,
    'xtick.color': MUTED, 'ytick.color': MUTED,
    'axes.labelcolor': INK2, 'text.color': INK,
})


def strip(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_axisbelow(True)


fig, axes = plt.subplots(2, 2, figsize=(13.2, 9.6), dpi=160)
fig.subplots_adjust(left=0.065, right=0.975, top=0.865, bottom=0.145,
                    hspace=0.62, wspace=0.26)
(axA, axB), (axC, axD) = axes

# ---------------------------------------------------------------- panel A
T = mpf('1473.15')
r0 = A.solve({'H2': 1}, 1473.15)
mu = O.mu_table(T)
n0 = {s: mpf(repr(v)) for s, v in r0['species_moles'].items() if v > 0}
n0['TiO2'] = mpf(repr(r0['solid_moles']['TiO2']))
n0['Ti10O19'] = mpf(repr(r0['solid_moles']['Ti10O19']))


def G(n):
    Ng = sum(n[s] for s in ('H2', 'H2O', 'O2') if n.get(s, 0) > 0)
    g = mpf(0)
    for s, v in n.items():
        if v <= 0:
            continue
        g += v * (mu[s] + O.R_KJ * T * (mlog(v / Ng) if s in ('H2', 'H2O', 'O2') else 0))
    return g


G0 = G(n0)
m_eq = n0['Ti10O19']
d_red = {'TiO2': -10, 'Ti10O19': 1, 'H2': -1, 'H2O': 1}
ds = np.logspace(-3.2, np.log10(0.85), 16)
branches = []
for sign in (-1, 1):
    xs, ys = [], []
    for d in ds:
        eps = m_eq * mpf(repr(float(sign * d)))
        n = dict(n0)
        for s, c in d_red.items():
            n[s] = n.get(s, mpf(0)) + eps * c
        dg = G(n) - G0
        xs.append(d)
        ys.append(float(dg) * 1000.0)      # kJ -> J
    branches.append((xs, ys))

axA.loglog(branches[0][0], branches[0][1], color=BLUE, lw=2, zorder=3)
axA.loglog(branches[1][0], branches[1][1], color=ORANGE, lw=2, zorder=3)
guide = [(1e-3, None), (0.5, None)]
gx = np.array([1.3e-3, 0.6])
gy = branches[1][1][6] / branches[1][0][6] ** 2 * gx ** 2
axA.loglog(gx, gy * 2.6, color=MUTED, lw=1.2, ls=(0, (4, 3)), zorder=2)
axA.text(6e-3, 3.2e-4, 'slope 2:\nquadratic bowl', color=MUTED, fontsize=8.5)
axA.text(1.35e-3, 2.4e-8, 'less reduced', color=BLUE, fontsize=9)
axA.text(1.35e-3, 3.2e-9, 'more reduced', color=ORANGE, fontsize=9)
axA.set_xlabel('|m(Ti$_{10}$O$_{19}$) / m$_{eq}$ − 1|   (displacement from the answer)')
axA.set_ylabel('G − G(answer)   (J)')
axA.set_title('A · The answer sits at the bottom of the bowl',
              loc='left', fontsize=11, fontweight='bold', color=INK, pad=26)
axA.text(0, 1.028, 'pure H$_2$, 1200 °C — G re-evaluated at 80 digits along the reduction extent, both directions',
         transform=axA.transAxes, fontsize=8.3, color=INK2)
strip(axA)

# ---------------------------------------------------------------- panel B
dense_x, dense_y = [], []
for mg in np.logspace(np.log10(0.04), np.log10(400), 80):
    r = A.solve({'H2': 1}, 1473.15, mass_g=mg / 1000)
    dense_x.append(r['n_gas_charge_mol'] / r['n_Ti_mol'])
    dense_y.append(r['reduced_pct'])
axB.plot(dense_x, dense_y, color=MUTED, lw=1.4, zorder=2)

mass_pts = [(mg, A.solve({'H2': 1}, 1473.15, mass_g=mg / 1000))
            for mg in (100, 20, 8.5, 8, 5, 2, 0.8, 0.3, 0.1)]
vol_pts = [(mL, A.solve({'H2': 1}, 1473.15, V_cm3=mL))
           for mL in (25, 100, 300, 1000, 3000, 10000)]
axB.scatter([r['n_gas_charge_mol'] / r['n_Ti_mol'] for _, r in mass_pts],
            [r['reduced_pct'] for _, r in mass_pts],
            s=68, color=BLUE, zorder=4, edgecolors=SURF, linewidths=1.2)
axB.scatter([r['n_gas_charge_mol'] / r['n_Ti_mol'] for _, r in vol_pts],
            [r['reduced_pct'] for _, r in vol_pts],
            s=74, color=ORANGE, marker='D', zorder=5,
            edgecolors=SURF, linewidths=1.2)
axB.set_xscale('log')
for x, y, lab in ((2.9, 21.3, 'Ti$_{10}$O$_{19}$'), (13.5, 27.6, 'Ti$_8$O$_{15}$'),
                  (26, 41.3, 'Ti$_5$O$_9$'), (66, 53.0, 'Ti$_4$O$_7$'),
                  (180, 62.5, 'Ti$_3$O$_5$')):
    axB.text(x, y, lab, fontsize=8.5, color=INK2, ha='center')
axB.text(0.34, 6.2, 'TiO$_2$ + Ti$_{10}$O$_{19}$\nbuffer', fontsize=8.5,
         color=INK2)
axB.set_xlabel('gas charge / titanium   (mol H$_2$ per mol Ti, log)')
axB.set_ylabel('Ti$^{3+}$  (% of total Ti)')
axB.set_title('B · One variable walks the whole ladder — the oxygen budget',
              loc='left', fontsize=11, fontweight='bold', color=INK, pad=26)
axB.text(0, 1.028, 'pure H$_2$, 1200 °C fixed — both sweeps collapse onto one curve; the 14 points are oracle-verified',
         transform=axB.transAxes, fontsize=8.3, color=INK2)
axB.legend(handles=[
    Line2D([], [], color=BLUE, marker='o', ls='none', ms=8.5,
           markeredgecolor=SURF, label='mass sweep (25 mL fixed)'),
    Line2D([], [], color=ORANGE, marker='D', ls='none', ms=8,
           markeredgecolor=SURF, label='volume sweep (100 mg fixed)'),
], loc='upper left', frameon=False, fontsize=8.5, labelcolor=INK2)
strip(axB)

# ---------------------------------------------------------------- panel C
Ts = list(range(500, 1501, 100))
crit = []
for T_C in Ts:
    lo, hi = 0.0, 0.5
    for _ in range(46):
        mid = 0.5 * (lo + hi)
        rr = A.solve({'H2': 1 - mid, 'CO2': mid}, T_C + 273.15)
        if len(rr['active_condensed_phases']) == 2:
            lo = mid
        else:
            hi = mid
    crit.append(0.5 * (lo + hi) * 100)
axC.semilogy(Ts, crit, color=BLUE, lw=2, zorder=3)
axC.fill_between(Ts, 1e-5, crit, color=SEQ100, alpha=0.55, zorder=1)
axC.set_ylim(1e-5, 60)
axC.set_xlim(500, 1560)
axC.text(640, 3e-4, 'reduces:\nTiO$_2$ + Ti$_{10}$O$_{19}$ buffer',
         fontsize=9, color=INK2)
axC.text(660, 4.5, 'no reduction — Ti$^{3+}$ = 0, exact\n(positive KKT reduced cost)',
         fontsize=9, color=INK2)
for T_on, y_on, lab, dx, dy in ((778, 0.1, 'H$_2$:CO$_2$ = 999:1', 12, -0.055),
                                (992, 1.0, '99:1', 18, -0.62),
                                (1329, 10.0, '9:1', 20, -6.0)):
    axC.scatter([T_on], [y_on], s=64, color=BLUE, zorder=4,
                edgecolors=SURF, linewidths=1.2)
    axC.text(T_on + dx, y_on + dy, lab, fontsize=8.5, color=INK2)
axC.annotate('3:1 (25% CO$_2$) → ~1574 °C, extrapolated', xy=(1500, 20.1),
             xytext=(560, 26), fontsize=8.5, color=INK2,
             arrowprops=dict(arrowstyle='-', color=MUTED, lw=1.0))
axC.set_xlabel('temperature  (°C)')
axC.set_ylabel('CO$_2$ in the H$_2$ feed  (mol %, log)')
axC.set_title('C · Where H$_2$-rich feeds cross into reduction',
              loc='left', fontsize=11, fontweight='bold', color=INK, pad=26)
axC.text(0, 1.028, '100 mg TiO$_2$, 25 mL, 1 atm — the boundary composition at each temperature, bisected on the active set',
         transform=axC.transAxes, fontsize=8.3, color=INK2)
strip(axC)

# ---------------------------------------------------------------- panel D
rows = [
    ('gas fractions, Python ↔ oracle', 3.5e-13, 1e-9),
    ('reduced costs r$_j$, Python ↔ oracle (kJ/mol O)', 5.1e-11, 1e-8),
    ('trace amounts log$_{10}$n, Python ↔ oracle', 1.7e-13, 1e-6),
    ('gas fractions, Python ↔ JS', 6.1e-15, 1e-9),
    ('reduced costs r$_j$, Python ↔ JS (kJ/mol O)', 1.8e-12, 1e-9),
    ('G$_{total}$ relative, Python ↔ JS', 4.5e-12, 1e-9),
]
ypos = np.arange(len(rows))[::-1]
for y, (lab, v, gate) in zip(ypos, rows):
    axD.plot([1e-16, v], [y, y], color=GRID, lw=1.2, zorder=1)
    axD.scatter([v], [y], s=72, color=BLUE, zorder=3,
                edgecolors=SURF, linewidths=1.2)
    axD.scatter([gate], [y], s=110, color=MUTED, marker='|', zorder=2)
    axD.text(v * 2.4, y, f'{v:.1e}', fontsize=8.5, color=INK2, va='center')
axD.set_xscale('log')
axD.set_xlim(1e-16, 3e-5)
axD.set_yticks(ypos)
axD.set_yticklabels([r[0] for r in rows], fontsize=8.8, color=INK2)
axD.set_ylim(-0.7, len(rows) - 0.3)
axD.grid(axis='y', visible=False)
axD.text(1.1e-9, len(rows) - 0.62, 'release gate ‖', fontsize=8.5,
         color=MUTED)
axD.set_xlabel('worst deviation across all 52 verified points  (log)')
axD.set_title('D · Verification ledger — worst case vs the release gate',
              loc='left', fontsize=11, fontweight='bold', color=INK, pad=26)
axD.text(0, 1.028, '38 reference + 14 sweep points; ticks ‖ mark the gate each row must beat (enforced in tests/)',
         transform=axD.transAxes, fontsize=8.3, color=INK2)
strip(axD)

fig.suptitle('Ti–O active-set solver, verified: a true Gibbs minimum, an oxygen-budget ladder',
             x=0.065, y=0.978, ha='left', fontsize=13, fontweight='bold', color=INK)
fig.text(0.065, 0.938, 'solidgas/activeset.py · gases NIST-JANAF, solids Waldner & Eriksson (1999) · '
                       'every point method = gibbs_min, decided by KKT reduced costs',
         fontsize=9, color=INK2)
fig.text(0.065, 0.055,
         'A: every feasible displacement raises G — the two branches collapse on the slope-2 line, i.e. a clean quadratic minimum with zero gradient.   '
         'B: mass and volume sweeps land on one curve —\nreduction depth is set by mol gas per mol Ti, not by temperature; treads are single-phase windows, ramps are two-phase buffers.   '
         'C: the same boundary seen from the composition axis.\n'
         'D: solver-to-solver and solver-to-oracle agreement, with the gates the test suite enforces.',
         fontsize=8.6, color=INK2, va='top')

out = ROOT / 'docs' / 'assets' / 'verification_figure.png'
fig.savefig(out, facecolor=SURF)
print('written', out)
