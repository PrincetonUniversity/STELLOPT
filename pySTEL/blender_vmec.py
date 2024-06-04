import sys,os
import bpy
import numpy as np
sys.path.append("/Users/lazerson/Sims_Work/STELLOPT/pySTEL/")
os.environ["STELLOPT_PATH"] = "/Users/lazerson/Sims_Work/STELLOPT/"
wout_file = '/Users/lazerson/Sims_Work/GAUSS/GIGA_v120/wout_GIGA_v120.nc'

from libstell.vmec import VMEC


# Define the RGB color (values between 0 and 1)
rgb_color = (0.1, 0.8, 0.2)  # Example: a shade of green

vmec_wout = VMEC()
vmec_wout.read_wout(wout_file)
nu = 360
nv = 360
theta = np.ndarray((nu,1))
phi   = np.ndarray((nv,1))
for j in range(nu):  theta[j]=2.0*np.pi*j/(nu)
for j in range(nv):    phi[j]=2.0*np.pi*j/(nv)
r = vmec_wout.cfunct(theta,phi,vmec_wout.rmnc,vmec_wout.xm,vmec_wout.xn)
z = vmec_wout.sfunct(theta,phi,vmec_wout.zmns,vmec_wout.xm,vmec_wout.xn)
[vertices,faces] = vmec_wout.blenderSurface(r,z,phi,surface=vmec_wout.ns-1)

# Create mesh data
mesh = bpy.data.meshes.new(name="VMECMesh")
mesh.from_pydata(vertices, [], faces)
mesh.update()

# Create object
vmec_obj = bpy.data.objects.new("VMECObject", mesh)

# Create a new material
material = bpy.data.materials.new(name="PlasmaMaterial")
material.use_nodes = True

# Clear default nodes
nodes = material.node_tree.nodes
for node in nodes:
    nodes.remove(node)

# Create new nodes
output_node = nodes.new(type='ShaderNodeOutputMaterial')
emission_node = nodes.new(type='ShaderNodeEmission')
transparent_node = nodes.new(type='ShaderNodeBsdfTransparent')
mix_shader_node = nodes.new(type='ShaderNodeMixShader')
noise_texture_node = nodes.new(type='ShaderNodeTexNoise')
color_ramp_node = nodes.new(type='ShaderNodeValToRGB')

# Arrange nodes
noise_texture_node.location = (-600, 0)
color_ramp_node.location = (-400, 0)
emission_node.location = (-200, 0)
transparent_node.location = (-200, -200)
mix_shader_node.location = (0, 0)
output_node.location = (200, 0)

# Link nodes
links = material.node_tree.links
links.new(noise_texture_node.outputs['Fac'], color_ramp_node.inputs['Fac'])
links.new(color_ramp_node.outputs['Color'], emission_node.inputs['Color'])
links.new(emission_node.outputs['Emission'], mix_shader_node.inputs[2])
links.new(transparent_node.outputs['BSDF'], mix_shader_node.inputs[1])
links.new(mix_shader_node.outputs['Shader'], output_node.inputs['Surface'])

# Adjust noise texture
noise_texture_node.inputs['Scale'].default_value = 20.0
noise_texture_node.inputs['Detail'].default_value = 16.0
noise_texture_node.inputs['Distortion'].default_value = 1.0

# Adjust color ramp
color_ramp = color_ramp_node.color_ramp
color_ramp.interpolation = 'EASE'
color_ramp.elements[0].position = 0.0
color_ramp.elements[0].color = (0.0, 0.1, 0.0, 1.0)  # Dark green
color_ramp.elements.new(0.3)
color_ramp.elements[1].position = 0.3
color_ramp.elements[1].color = (0.0, 0.5, 0.0, 1.0)  # Medium green
color_ramp.elements.new(0.6)
color_ramp.elements[2].position = 0.6
color_ramp.elements[2].color = (0.0, 1.0, 0.0, 1.0)  # Bright green
color_ramp.elements.new(1.0)
color_ramp.elements[3].position = 1.0
color_ramp.elements[3].color = (0.5, 1.0, 0.5, 1.0)  # Light green


# Ensure emission strength is sufficient
emission_node.inputs['Strength'].default_value = 1.0

# Adjust mix shader to balance between emission and transparency
mix_shader_node.inputs['Fac'].default_value = 0.5  # Adjust this value to change transparency

# Ensure the object has an active material slot
if vmec_obj.data.materials:
    # Assign to the first material slot
    vmec_obj.data.materials[0] = material
else:
    # Create a new material slot and assign the material
    vmec_obj.data.materials.append(material)


# Link object to the scene
scene = bpy.context.scene
scene.collection.objects.link(vmec_obj)

# Set render engine to Eevee and enable Bloom
bpy.context.scene.render.engine = 'BLENDER_EEVEE'
bpy.context.scene.eevee.use_bloom = True
bpy.context.scene.eevee.bloom_intensity = 0.2  # Adjust bloom intensity as needed

# Set object as active and select it
bpy.context.view_layer.objects.active = vmec_obj
vmec_obj.select_set(True)

#####