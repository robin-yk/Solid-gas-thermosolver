"""Scene, material, light and camera helpers shared by the two renders of
the vacancy-distribution workspace. Needs the bpy module (Blender 5)."""
import os

import bpy
from mathutils import Vector
from PIL import Image


def reset(res, samples):
    bpy.ops.wm.read_factory_settings(use_empty=True)
    sc = bpy.context.scene
    sc.render.engine = 'CYCLES'
    sc.cycles.device = 'CPU'
    sc.cycles.samples = samples
    sc.cycles.use_denoising = True
    sc.cycles.seed = 1
    sc.render.resolution_x = sc.render.resolution_y = res
    sc.render.film_transparent = True
    sc.view_settings.view_transform = 'AgX'
    sc.view_settings.look = 'AgX - Medium High Contrast'
    w = bpy.data.worlds.new('w')
    w.use_nodes = True
    bg = w.node_tree.nodes['Background']
    bg.inputs[0].default_value = (0.92, 0.94, 0.96, 1)
    bg.inputs[1].default_value = 0.5
    sc.world = w
    return sc


def linear(hexcol, mix_white=0.0):
    """sRGB hex -> linear RGB, optionally mixed toward white first."""
    h = hexcol.lstrip('#')
    c = [int(h[i:i + 2], 16) / 255 for i in (0, 2, 4)]
    c = [x + (1 - x) * mix_white for x in c]
    return tuple(((x + 0.055) / 1.055) ** 2.4 if x > 0.04045 else x / 12.92 for x in c)


def mat(name, rgb, rough=0.45, coat=0.2, emit=None, alpha=1.0, sss=0.0):
    m = bpy.data.materials.new(name)
    m.use_nodes = True
    p = m.node_tree.nodes['Principled BSDF']
    p.inputs['Base Color'].default_value = (*rgb, 1)
    p.inputs['Roughness'].default_value = rough
    p.inputs['Coat Weight'].default_value = coat
    p.inputs['Coat Roughness'].default_value = 0.25
    if sss:
        p.inputs['Subsurface Weight'].default_value = sss
        p.inputs['Subsurface Radius'].default_value = (0.3, 0.3, 0.3)
        p.inputs['Subsurface Scale'].default_value = 0.05
    if emit:
        p.inputs['Emission Color'].default_value = (*emit[0], 1)
        p.inputs['Emission Strength'].default_value = emit[1]
    if alpha < 1:
        p.inputs['Alpha'].default_value = alpha
    return m


def flat(name, rgb):
    """Emission-only material for the hover-region mask."""
    m = bpy.data.materials.new(name)
    m.use_nodes = True
    nt = m.node_tree
    nt.nodes.clear()
    e = nt.nodes.new('ShaderNodeEmission')
    e.inputs['Color'].default_value = (*rgb, 1)
    e.inputs['Strength'].default_value = 1.0
    o = nt.nodes.new('ShaderNodeOutputMaterial')
    nt.links.new(e.outputs[0], o.inputs[0])
    return m


def area(loc, target, size, power, rgb=(1, 1, 1)):
    bpy.ops.object.light_add(type='AREA', location=loc)
    L = bpy.context.object
    L.data.size = size
    L.data.energy = power
    L.data.color = rgb
    L.rotation_euler = (Vector(target) - Vector(loc)).to_track_quat('-Z', 'Y').to_euler()
    return L


def camera(loc, target, lens):
    bpy.ops.object.camera_add(location=loc)
    c = bpy.context.object
    c.rotation_euler = (Vector(target) - Vector(loc)).to_track_quat('-Z', 'Y').to_euler()
    c.data.lens = lens
    bpy.context.scene.camera = c
    return c


def catcher(z):
    bpy.ops.mesh.primitive_plane_add(size=60, location=(0, 0, z))
    g = bpy.context.object
    g.is_shadow_catcher = True
    return g


def render_to(path, quality=86):
    """Render, flatten onto white, write WebP (or PNG for masks)."""
    tmp = path + '.tmp.png'
    bpy.context.scene.render.filepath = tmp
    bpy.ops.render.render(write_still=True)
    im = Image.open(tmp).convert('RGBA')
    os.remove(tmp)
    if path.endswith('.png'):
        im.save(path, optimize=True)
        return
    bg = Image.new('RGBA', im.size, (255, 255, 255, 255))
    bg.alpha_composite(im)
    bg.convert('RGB').save(path, quality=quality, method=6)
