"""Draw web/schemes/vacancy-distribution.svg, the workspace scheme.

    python3 scripts/scheme_vacancy_distribution.py     (matplotlib 3.9.4)

Same plate as the other two schemes: 5 x 4 in, Helvetica, role hues, the
first three texts being the plate letter, title and subtitle that
build_site.py strips. Conceptual geometry: the surface band is drawn far
thicker than its 1.3 nm, and the atoms are not to scale."""
import math
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                      # noqa: E402
from matplotlib.patches import Circle, Wedge, Rectangle, FancyArrowPatch  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, 'web', 'schemes', 'vacancy-distribution.svg')
SURF, SUB, BULK, GREY, INK = '#0072b2', '#d55e00', '#6a51a3', '#777777', '#222222'
plt.rcParams.update({'svg.fonttype': 'none', 'font.family': 'Helvetica',
                     'svg.hashsalt': 'vacancy-distribution'})

fig = plt.figure(figsize=(5, 4))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 360)
ax.set_ylim(288, 0)
ax.axis('off')
ax.text(20.4, 33.8, 'b', fontsize=15, fontweight='bold', color=INK)
ax.text(180, 33.8, 'Vacancy distribution', fontsize=15, ha='center', color=INK)
ax.text(180, 67.0, 'One inventory, every pool', fontsize=12, ha='center', color=INK)

# particle: bulk interior, a (110) band drawn thick, bulk-like arcs left bare
cx, cy, R = 118, 176, 64
ax.add_patch(Circle((cx, cy), R, fc='#f6f6f6', ec='#9b9b9b', lw=1.0))
for a0, a1 in ((-60, 70), (110, 200), (235, 285)):
    ax.add_patch(Wedge((cx, cy), R, a0, a1, width=9, fc='#dcebf5', ec='none'))
ax.add_patch(Circle((cx, cy), R, fc='none', ec='#9b9b9b', lw=1.0))
pts = [(-30, 12), (-8, -30), (18, 20), (26, -12), (-24, -10), (4, 36), (-40, 28), (38, 8), (6, 2), (-10, 18)]
for dx, dy in pts:
    ax.add_patch(Circle((cx + dx, cy + dy), 2.6, fc=BULK, ec='none'))
for ang in (-40, 10, 50, 130, 170, 250, 270):
    t = math.radians(ang)
    ax.add_patch(Circle((cx + (R - 4) * math.cos(t), cy - (R - 4) * math.sin(t)), 3.0, fc='none', ec=SURF, lw=1.4))
for ang in (-20, 150):
    t = math.radians(ang)
    ax.add_patch(Circle((cx + (R - 12) * math.cos(t), cy - (R - 12) * math.sin(t)), 2.6, fc='none', ec=SUB, lw=1.4))

# zoom: the (110) band as four trilayers, one bridging vacancy
zx, zy = 222, 118
ax.add_patch(Rectangle((zx, zy), 118, 104, fc='white', ec='#c6c6c6', lw=1.0))
ax.add_patch(FancyArrowPatch((cx + R * 0.5, cy - R * 0.87 + 4), (zx, zy + 40),
                             arrowstyle='-', color='#c6c6c6', lw=1.0, linestyle=(0, (3, 2))))
for k in range(4):
    y = zy + 22 + 22 * k
    for i in range(6):
        x = zx + 14 + 18 * i
        ax.add_patch(Circle((x, y), 3.2, fc='#9b9b9b', ec='none'))
        if k == 0 and i % 2 == 0:
            if i == 2:
                ax.add_patch(Circle((x, y - 9), 4.2, fc='none', ec=SURF, lw=1.6))
            else:
                ax.add_patch(Circle((x, y - 9), 4.2, fc='white', ec='#9b9b9b', lw=1.0))
ax.add_patch(Circle((zx + 14 + 18 * 3, zy + 22 + 22 + 9), 3.6, fc='none', ec=SUB, lw=1.6))

ax.text(zx + 50, zy - 8, 'Reactive site', fontsize=10, ha='center', color=SURF)
ax.text(zx + 59, zy + 118, 'Subsurface', fontsize=10, ha='center', color=SUB)
ax.text(cx, cy + R + 16, 'Bulk', fontsize=10, ha='center', color=BULK)
ax.text(cx - R - 6, cy - R + 6, 'Surface', fontsize=10, ha='right', color=SURF)
ax.text(zx + 59, zy + 132, 'Band not to scale', fontsize=10, ha='center', color=GREY)
fig.savefig(OUT, metadata={'Date': None})
print('wrote', os.path.relpath(OUT, ROOT))
