"""Typeset equation sheet for the Ti-O active-set Gibbs solver.

matplotlib mathtext with the Computer Modern fontset: real LaTeX-look
typesetting, no external TeX needed. Output: PNG + PDF.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ASSETS = ROOT / 'docs' / 'assets'

BLUE = '#2a78d6'
SEQ100 = '#cde2fb'
SURF, INK, INK2, MUTED = '#fcfcfb', '#0b0b0b', '#52514e', '#898781'
GRID = '#e1e0d9'

plt.rcParams.update({
    'mathtext.fontset': 'cm',
    'font.family': 'DejaVu Sans',
    'figure.facecolor': SURF,
})

from equations_registry import R

# ---------------------------------------------------------------- layout
H = {'title': 0.55, 'sub': 0.34, 'head': 0.62, 'eq': 0.55, 'boxeq': 0.60,
     'note': 0.34}
total = sum(H.get(r[0], r[1] if r[0] == 'gap' else 0.4) for r in R) + 0.7
W = 11.4
fig = plt.figure(figsize=(W, total), dpi=150)
fig.patch.set_facecolor(SURF)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(0, total)
ax.axis('off')

y = total - 0.45
LM, NUMX = 0.75, W - 0.55
for row in R:
    kind = row[0]
    if kind == 'gap':
        y -= row[1]
        continue
    h = H[kind]
    if kind == 'title':
        ax.text(0.55, y, row[1], fontsize=19, fontweight='bold', color=INK,
                va='center')
    elif kind == 'sub':
        ax.text(0.55, y, row[1], fontsize=10.5, color=INK2, va='center')
    elif kind == 'head':
        y -= 0.16
        ax.text(0.55, y, row[1], fontsize=12.5, fontweight='bold',
                color=BLUE, va='center')
        ax.plot([0.55, W - 0.55], [y - 0.21, y - 0.21], color=GRID, lw=1.1)
        y -= 0.10
    elif kind in ('eq', 'boxeq'):
        if kind == 'boxeq':
            ax.add_patch(plt.Rectangle((0.58, y - h / 2 + 0.02), W - 1.16, h - 0.06,
                                       facecolor=SEQ100, alpha=0.38,
                                       edgecolor='none', zorder=0))
            ax.plot([0.58, 0.58], [y - h / 2 + 0.02, y + h / 2 - 0.04],
                    color=BLUE, lw=3, zorder=1)
        ax.text(LM, y, row[1], fontsize=13, color=INK, va='center')
        ax.text(NUMX, y, f'({row[2]})', fontsize=10, color=MUTED,
                va='center', ha='right')
    elif kind == 'note':
        ax.text(LM, y, row[1], fontsize=9.8, color=INK2, va='center')
    y -= h

ASSETS.mkdir(parents=True, exist_ok=True)
fig.savefig(ASSETS / 'equations.png', facecolor=SURF)
fig.savefig(ASSETS / 'equations.pdf', facecolor=SURF)
print('written docs/assets/equations.png and .pdf, height', round(total, 1), 'in')
