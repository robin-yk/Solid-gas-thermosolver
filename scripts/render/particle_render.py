"""Particle panel of the vacancy-distribution workspace: rounded, irregular
rutile grains clumped as in the TEM images (manuscript Fig. S10d), one
grain cut open to show its surface band, subsurface band and bulk.

Writes web/render/particle_<state>.webp for the states default, surface,
subsurface and bulk (the named layer highlighted), and
web/render/particle_mask.png, whose red, green and blue pixels mark the
surface, subsurface and bulk sections for hover lookup.

    <python with the bpy module> scripts/render/particle_render.py

Display geometry: band thicknesses are not to scale."""
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bpy  # noqa: E402
import bmesh  # noqa: E402
from mathutils import Vector, noise  # noqa: E402

from blender_common import area, camera, catcher, flat, linear, mat, render_to, reset  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, 'web', 'render')
ROLE = {'surface': '#0072B2', 'subsurface': '#D55E00', 'bulk': '#6A51A3'}
LAYERS = (('surface', 1.0, 0.0), ('subsurface', 0.9, 0.003), ('bulk', 0.8, 0.006))
GRAINS = ((11, 0.8, (-1.45, 0.6, -0.25)), (17, 0.7, (0.2, 1.55, -0.3)),
          (23, 0.62, (-0.9, 1.6, 0.35)), (29, 0.55, (1.35, 0.9, -0.45)),
          (31, 0.5, (-1.7, -0.6, -0.5)))


def blob(seed, scale):
    """Noisy icosphere with a few blunt planar cuts, bevelled and smoothed."""
    bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=5, radius=1)
    ob = bpy.context.object
    off = Vector((seed * 7.1, seed * 3.3, seed * 5.7))
    for v in ob.data.vertices:
        v.co *= 1 + noise.noise(v.co * 1.3 + off) * 0.16 + noise.noise(v.co * 3.1 + off) * 0.04
    bm = bmesh.new()
    bm.from_mesh(ob.data)
    rnd = random.Random(seed)
    for _ in range(6):
        d = Vector((rnd.uniform(-1, 1), rnd.uniform(-1, 1), rnd.uniform(-1, 1))).normalized()
        cut = rnd.uniform(0.78, 0.92)
        bmesh.ops.bisect_plane(bm, geom=bm.verts[:] + bm.edges[:] + bm.faces[:], dist=1e-4,
                               plane_co=d * cut, plane_no=d, clear_outer=True)
        edges = [e for e in bm.edges if e.is_boundary]
        if edges:
            bmesh.ops.holes_fill(bm, edges=edges, sides=0)
    bm.to_mesh(ob.data)
    bm.free()
    ob.scale = (scale, scale * 0.85, scale * 0.9)
    bpy.ops.object.transform_apply(scale=True)
    for kind, kw in (('BEVEL', dict(width=0.03 * scale, segments=4, limit_method='ANGLE')),
                     ('SUBSURF', dict(levels=1, render_levels=1))):
        mod = ob.modifiers.new(kind, kind)
        for k, v in kw.items():
            setattr(mod, k, v)
        bpy.ops.object.modifier_apply(modifier=kind)
    bpy.ops.object.shade_smooth()
    return ob


def wedge_cut(ob, cutmat, offset):
    bpy.ops.mesh.primitive_cube_add(size=1)
    k = bpy.context.object
    k.scale = (6, 6, 6)
    k.location = (3 + offset, -3 - offset, 3 + offset)
    k.data.materials.append(cutmat)
    mod = ob.modifiers.new('cut', 'BOOLEAN')
    mod.operation = 'DIFFERENCE'
    mod.solver = 'EXACT'
    mod.material_mode = 'TRANSFER'
    mod.object = k
    bpy.context.view_layer.objects.active = ob
    bpy.ops.object.modifier_apply(modifier='cut')
    bpy.data.objects.remove(k)


def scene(state, mask=False):
    reset(res=500 if mask else 1000, samples=1 if mask else 128)
    sc = bpy.context.scene
    if mask:
        sc.cycles.use_denoising = False
        sc.render.filter_size = 0.01
        sc.view_settings.view_transform = 'Standard'
        sc.view_settings.look = 'None'
        skin = flat('skin', (0, 0, 0))
        cut = {k: flat(k, rgb) for k, rgb in
               (('surface', (1, 0, 0)), ('subsurface', (0, 1, 0)), ('bulk', (0, 0, 1)))}
    else:
        skin = mat('skin', linear('#e9e4dc'), rough=0.45, coat=0.3, sss=0.06)
        cut = {}
        for name, hexcol in ROLE.items():
            if state == 'default':
                cut[name] = mat(name, linear(hexcol, 0.45), rough=0.6, coat=0.1)
            elif state == name:
                cut[name] = mat(name, linear(hexcol, 0.08), rough=0.5, coat=0.1,
                                emit=(linear(hexcol, 0.2), 0.35))
            else:
                cut[name] = mat(name, linear(hexcol, 0.82), rough=0.7, coat=0.1)
    base = blob(3, 1.0)
    me = base.data.copy()
    bpy.data.objects.remove(base)
    for name, f, off in LAYERS:
        ob = bpy.data.objects.new(name, me.copy())
        bpy.context.scene.collection.objects.link(ob)
        ob.scale = (f, f, f)
        bpy.ops.object.select_all(action='DESELECT')
        ob.select_set(True)
        bpy.context.view_layer.objects.active = ob
        bpy.ops.object.transform_apply(scale=True)
        ob.data.materials.clear()
        ob.data.materials.append(skin)
        wedge_cut(ob, cut[name], off)
    for seed, s, loc in GRAINS:
        g = blob(seed, s)
        g.location = loc
        g.data.materials.append(skin)
    if not mask:
        catcher(-1.05)
    camera((4.6, -6.2, 3.9), (-0.25, 0.4, -0.1), lens=58)
    area((-3.5, -4.5, 6.5), (0, 0, 0), 4.5, 420, (1, 0.96, 0.9))
    area((6, -1, 1.5), (0, 0, 0), 6, 110, (0.88, 0.94, 1))
    area((-1, 6, 3.5), (0, 0, 0), 4, 380, (0.9, 0.95, 1))


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    for state in ('default', 'surface', 'subsurface', 'bulk'):
        scene(state)
        render_to(os.path.join(OUT, 'particle_%s.webp' % state))
    scene('default', mask=True)
    render_to(os.path.join(OUT, 'particle_mask.png'))
