import sys,os
import bpy
import numpy as np
sys.path.append("/Users/lazerson/Sims_Work/STELLOPT/pySTEL/")
os.environ["STELLOPT_PATH"] = "/Users/lazerson/Sims_Work/STELLOPT/"
from libstell.coils import COILSET

coil_file = '/Users/lazerson/Sims_Work/GAUSS/GIGA_v120_FOCUS/GIGA_v120_c100.coils'

# Generate Coil
coils = COILSET()
coils.read_coils_file(coil_file)
[vertices_coil,faces_coil]=coils.blenderCoil(dist=0.5)

# Create mesh        
mesh_coil = bpy.data.meshes.new(name="CoilMesh")
mesh_coil.from_pydata(vertices_coil, [], faces_coil)
mesh_coil.update()

# Create object
coil_obj = bpy.data.objects.new("CoilObject", mesh_coil)

# Create a new material
material = bpy.data.materials.new(name="CopperMaterial")
material.use_nodes = True

# Clear default nodes
nodes = material.node_tree.nodes
nodes.clear()

# Create Principled BSDF shader node
principled_node = nodes.new(type='ShaderNodeBsdfPrincipled')

# Set Principled BSDF parameters for copper
principled_node.inputs['Base Color'].default_value = (0.8, 0.4, 0.2, 1.0)  # Reddish-brown color
principled_node.inputs['Metallic'].default_value = 1.0  # Fully metallic
principled_node.inputs['Roughness'].default_value = 0.3  # Moderate roughness
principled_node.inputs['Specular'].default_value = 0.5  # Specular intensity
principled_node.inputs['Anisotropic'].default_value = 0.2  # Anisotropic reflection

# Create material output node
output_node = nodes.new(type='ShaderNodeOutputMaterial')

# Link nodes
links = material.node_tree.links
links.new(principled_node.outputs['BSDF'], output_node.inputs['Surface'])

if coil_obj.data.materials:
    # Assign to the first material slot
    coil_obj.data.materials[0] = material
else:
    # Create a new material slot and assign the material
    coil_obj.data.materials.append(material)


# Link object to the scene
scene = bpy.context.scene
scene.collection.objects.link(coil_obj)


# Set object as active and select it
bpy.context.view_layer.objects.active = coil_obj
coil_obj.select_set(True)

