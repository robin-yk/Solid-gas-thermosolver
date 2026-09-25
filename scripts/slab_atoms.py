"""Atoms of the rutile (110) surface drawn in the vacancy-distribution
workspace: four O-Ti2O2-O trilayers of 8 x 12 cells, generated from
scripts/render/rutile.cif by its symmetry operations. Unrelaxed bulk
positions.

Every O site carries a fixed random rank and is vacant when its rank is
below the calculated site fraction of its class, so a sample with larger
fractions keeps every vacancy of a smaller one. Each O also lists its two
nearest Ti, where the page draws Ti3+ (the model resolves Ti3+ by layer,
not by site).

Writes web/render/slab_atoms.json. `vacant` is the reference the page's
mirror (web/slab_canvas.js) is tested against.

    python3 scripts/slab_atoms.py"""
import json
import math
import os
import random
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CIF = os.path.join(ROOT, 'scripts', 'render', 'rutile.cif')
OUT = os.path.join(ROOT, 'web', 'render', 'slab_atoms.json')
NX, NY, NL = 8, 12, 4
SEED = 2959
CLASSES = ('BRI', 'IPL', 'SBR', 'L24')


def read_cif(path):
    t = open(path).read()
    a = float(re.search(r'_cell_length_a\s+(\S+)', t).group(1))
    c = float(re.search(r'_cell_length_c\s+(\S+)', t).group(1))
    ops = re.findall(r'^\s*([-\d/xyz+]+,[-\d/xyz+]+,[-\d/xyz+]+)\s*$', t, re.M)
    sites = re.findall(r'^\s*\w+\s+(Ti|O)\s+(\S+)\s+(\S+)\s+(\S+)\s*$', t, re.M)
    frac = set()
    for el, x, y, z in sites:
        env = dict(x=float(x), y=float(y), z=float(z))
        for op in ops:
            v = [eval(e, {}, env) % 1.0 for e in op.split(',')]      # noqa: S307 (CIF symmetry strings)
            frac.add((el,) + tuple(round(w, 5) % 1.0 for w in v))
    return a, c, sorted(frac)


def atoms():
    """[(element, site class, trilayer, x, y, z in angstrom)], z up, surface at 0."""
    a, c, basis = read_cif(CIF)
    assert len(basis) == 6, basis
    px, d = a * math.sqrt(2), a / math.sqrt(2)
    out, seen = [], set()
    for i in range(-40, 41):
        for j in range(-40, 41):
            for k in range(NY):
                for el, x, y, z in basis:
                    X = ((x + i) - (y + j)) * d
                    Z = ((x + i) + (y + j)) * d
                    Y = (z + k) * c
                    home = round(Z / d)
                    if home > 0 or home < 1 - NL or X < -1e-6 or X >= NX * px - 1e-6:
                        continue
                    key = (el, round(X, 2), round(Y, 2), round(Z, 2))
                    if key in seen:
                        continue
                    seen.add(key)
                    dz = Z - home * d
                    site = el if el == 'Ti' else ('L24' if home < 0 else 'BRI' if dz > 0.5
                                                  else 'SBR' if dz < -0.5 else 'IPL')
                    out.append((el, site, 1 - home, X - NX * px / 2, Y - NY * c / 2, Z))
    return out


def vacant(doc, fractions):
    """Indices of vacant O for one set of site fractions {class: fraction}."""
    return [i for i, (s, r) in enumerate(zip(doc['site'], doc['rank']))
            if s in fractions and r < fractions[s]]


def build():
    at = atoms()
    rnd = random.Random(SEED)
    rank = [round(rnd.random(), 6) for _ in at]
    ti = [k for k, a in enumerate(at) if a[0] == 'Ti']
    near = []
    for a in at:
        if a[0] != 'O':
            near.append([])
            continue
        d = sorted(ti, key=lambda k: (at[k][3] - a[3]) ** 2 + (at[k][4] - a[4]) ** 2 + (at[k][5] - a[5]) ** 2)
        near.append(d[:2])
    return dict(cells=[NX, NY], trilayers=NL, seed=SEED, classes=list(CLASSES),
                el=[a[0] for a in at], site=[a[1] for a in at], layer=[a[2] for a in at],
                xyz=[[round(a[3], 3), round(a[4], 3), round(a[5], 3)] for a in at],
                rank=[r if a[0] == 'O' else None for a, r in zip(at, rank)],
                ti3=near)


if __name__ == '__main__':
    doc = build()
    with open(OUT, 'w') as fh:
        json.dump(doc, fh, separators=(',', ':'))
        fh.write('\n')
    print(OUT, len(doc['el']), 'atoms', os.path.getsize(OUT) // 1024, 'kB')
