"""Surface panel of the vacancy-distribution workspace: rutile (110), four
O-Ti2O2-O trilayers, atoms generated from scripts/render/rutile.cif by
its symmetry operations. Unrelaxed bulk positions.

One image per sample at the page's starting parameters, read from
paper_outputs/tof_range.json. Every O site carries a fixed random rank and
is vacant when its rank is below the calculated site fraction of its
class, so a higher-inventory sample keeps the vacancies of a lower one and
adds more. Ti3+ is drawn on the two Ti nearest each vacancy; the model
resolves Ti3+ by layer, not by site.

Writes web/render/slab_<sample>.webp and web/render/slab_sites.json.

    <python with the bpy module> scripts/render/slab_render.py"""
import json
import math
import os
import random
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bpy  # noqa: E402
from mathutils import Vector  # noqa: E402

from blender_common import area, camera, catcher, linear, mat, render_to, reset  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(ROOT, 'web', 'render')
NX, NY, NL = 8, 12, 4
START = dict(energy_map='PAB', cutoff_nm=0.34, dG_eV=0.0, f110=0.75, eps='a_axis')
R_O, R_TI = 1.15, 0.72


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


def slab_atoms():
    a, c, basis = read_cif(os.path.join(HERE, 'rutile.cif'))
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
                    out.append(dict(el=el, site=site, layer=1 - home,
                                    p=Vector((X - NX * px / 2, Y - NY * c / 2, Z))))
    return out


def fractions():
    doc = json.load(open(os.path.join(ROOT, 'paper_outputs', 'tof_range.json')))
    ax = doc['axes']
    i = ax['energy_map'].index(START['energy_map'])
    i = i * len(ax['cutoff_nm']) + ax['cutoff_nm'].index(START['cutoff_nm'])
    i = i * len(ax['dG_eV']) + ax['dG_eV'].index(START['dG_eV'])
    i = i * len(ax['f110']) + ax['f110'].index(START['f110'])
    i = i * len(ax['eps']) + [e['key'] for e in ax['eps']].index(START['eps'])
    cap = doc['capacity_umol_g']['%g' % START['f110']]
    out = {}
    for s in doc['samples']:
        pools = dict(zip(doc['pools'], s['cases']['pools'][i]))
        out[s['sample']] = dict(BRI=s['cases']['theta'][i],
                                IPL=pools['basal_L1'] / cap['basal_L1'],
                                SBR=pools['L1_subbridging'] / cap['L1_subbridging'],
                                L24=pools['subsurface_L2_4'] / cap['subsurface_L2_4'])
    return i, out


def render_sample(name, atoms, frac, rank):
    reset(res=1000, samples=128)
    vac = [a for a, r in zip(atoms, rank) if a['el'] == 'O' and r < frac[a['site']]]
    ti3 = []
    for v in vac:
        near = sorted((a for a in atoms if a['el'] == 'Ti'), key=lambda a: (a['p'] - v['p']).length)[:2]
        ti3 += [id(a) for a in near]
    mO = mat('O', linear('#e9e4dc'), rough=0.38, coat=0.35, sss=0.05)
    mT = mat('Ti', linear('#6f7c87'), rough=0.28, coat=0.3)
    mT3 = mat('Ti3', linear('#CC79A7'), rough=0.3, coat=0.4)
    ring = {'surface': mat('ring_s', linear('#0072B2'), rough=0.2, emit=(linear('#2f9bea'), 8.0)),
            'subsurface': mat('ring_b', linear('#D55E00'), rough=0.2, emit=(linear('#f08a3c'), 8.0))}
    bpy.ops.mesh.primitive_uv_sphere_add(segments=48, ring_count=24, radius=1)
    sph = bpy.context.object.data
    bpy.ops.object.shade_smooth()
    bpy.data.objects.remove(bpy.context.object)
    meshes = {}
    for key, m in (('O', mO), ('Ti', mT), ('Ti3', mT3)):
        meshes[key] = sph.copy()
        meshes[key].materials.append(m)
    gone = set(id(v) for v in vac)
    col = bpy.context.scene.collection
    for a in atoms:
        if id(a) in gone:
            continue
        key = 'O' if a['el'] == 'O' else ('Ti3' if id(a) in ti3 else 'Ti')
        ob = bpy.data.objects.new('atom', meshes[key])
        col.objects.link(ob)
        r = R_O if a['el'] == 'O' else R_TI
        ob.location = a['p']
        ob.scale = (r, r, r)
    for v in vac:
        surface = v['layer'] == 1 and v['site'] in ('BRI', 'IPL')
        bpy.ops.mesh.primitive_torus_add(major_radius=1.5 if surface else 1.1,
                                         minor_radius=0.12, location=v['p'] + Vector((0, 0, 0.1)))
        bpy.context.object.data.materials.append(ring['surface' if surface else 'subsurface'])
    catcher(-11.5)
    area((-20, -18, 34), (0, 0, -3), 22, 8000, (1, 0.97, 0.93))
    area((30, -20, 6), (0, 0, -3), 30, 2500, (0.9, 0.95, 1))
    area((5, 38, 16), (0, 0, -3), 18, 5000, (0.92, 0.96, 1))
    foc = Vector((0, 0, 1.3))
    cam = camera((foc.x - 18, foc.y - 56, foc.z + 34), (foc.x, foc.y + 2, foc.z - 2.5), lens=62)
    cam.data.dof.use_dof = True
    cam.data.dof.focus_distance = (cam.location - foc).length
    cam.data.dof.aperture_fstop = 0.6
    render_to(os.path.join(OUT, 'slab_%s.webp' % name))
    return {k: sum(1 for v in vac if v['site'] == k) for k in ('BRI', 'IPL', 'SBR', 'L24')}


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    atoms = slab_atoms()
    rnd = random.Random(2959)
    rank = [rnd.random() for _ in atoms]
    point, frac = fractions()
    only = sys.argv[1:]
    meta = dict(point=point, start=START, cells=[NX, NY], trilayers=NL, fractions=frac,
                sites={k: sum(1 for a in atoms if a['site'] == k) for k in ('BRI', 'IPL', 'SBR', 'L24')},
                vacant={})
    for name, f in frac.items():
        if only and name not in only:
            continue
        meta['vacant'][name] = render_sample(name, atoms, f, rank)
    if not only:
        with open(os.path.join(OUT, 'slab_sites.json'), 'w') as fh:
            json.dump(meta, fh, indent=1, sort_keys=True)
            fh.write('\n')
    print(json.dumps(meta['vacant']))
