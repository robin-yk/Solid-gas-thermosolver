"""Render the studio environment the 3D view reflects: a white cyclorama
with two soft boxes and a rim strip, as a 512 x 256 equirectangular image.

    <python with the bpy module> blender/studio_env.py

Lighting only; no data reaches this image."""
import math
import os
import bpy

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'assets', 'studio.png')

bpy.ops.wm.read_factory_settings(use_empty=True)
sc = bpy.context.scene
sc.render.engine = 'CYCLES'
sc.cycles.device = 'CPU'
sc.cycles.samples = 96
sc.render.resolution_x, sc.render.resolution_y = 512, 256
sc.render.image_settings.file_format = 'PNG'
sc.view_settings.view_transform = 'Standard'

world = bpy.data.worlds.new('w')
world.use_nodes = True
bg = world.node_tree.nodes['Background']
bg.inputs['Color'].default_value = (0.82, 0.85, 0.88, 1)
bg.inputs['Strength'].default_value = 0.55
sc.world = world


def softbox(name, loc, rot, size, strength, colour=(1, 1, 1)):
    bpy.ops.mesh.primitive_plane_add(size=1, location=loc, rotation=rot)
    ob = bpy.context.object
    ob.name = name
    ob.scale = (size[0], size[1], 1)
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()
    em = nt.nodes.new('ShaderNodeEmission')
    em.inputs['Color'].default_value = (*colour, 1)
    em.inputs['Strength'].default_value = strength
    out = nt.nodes.new('ShaderNodeOutputMaterial')
    nt.links.new(em.outputs[0], out.inputs[0])
    ob.data.materials.append(mat)


# floor and back sweep, slightly darker than the sky so objects read grounded
bpy.ops.mesh.primitive_plane_add(size=60, location=(0, 0, -3))
floor = bpy.context.object
fm = bpy.data.materials.new('floor')
fm.use_nodes = True
fm.node_tree.nodes['Principled BSDF'].inputs['Base Color'].default_value = (0.55, 0.58, 0.62, 1)
floor.data.materials.append(fm)

softbox('key', (-4, -3, 4), (math.radians(50), 0, math.radians(-50)), (5, 3.5), 9.0)
softbox('fill', (5, -2, 1.5), (math.radians(80), 0, math.radians(65)), (4, 5), 3.0, (0.93, 0.97, 1.0))
softbox('rim', (0, 5, 3), (math.radians(-60), 0, 0), (9, 0.8), 14.0)
softbox('top', (0, 0, 7), (0, 0, 0), (3, 3), 5.0)

bpy.ops.object.camera_add(location=(0, 0, 0), rotation=(math.radians(90), 0, 0))
cam = bpy.context.object
cam.data.type = 'PANO'
cam.data.panorama_type = 'EQUIRECTANGULAR'
sc.camera = cam
sc.render.filepath = OUT
bpy.ops.render.render(write_still=True)
# Blender's panorama puts the zenith at the bottom row; three.js expects it on top.
from PIL import Image
Image.open(OUT).transpose(Image.FLIP_TOP_BOTTOM).convert('RGB').save(OUT.replace('.png', '.jpg'), quality=88)
os.remove(OUT)
print('wrote', os.path.abspath(OUT.replace('.png', '.jpg')))
