#!/usr/bin/env python3

import sys
import gmsh

if __name__=="__main__":
	from argparse import ArgumentParser
	parser = ArgumentParser(description='''A tool for generating meshes.''')
	parser.add_argument("-i","--input", dest="in_file",
		help="Input file for conversion.", default = None)
	parser.add_argument("-o","--output", dest="out_file",
		help="Output file for conversion.", default = 'output.stl')
	parser.add_argument("--min", dest="min_val", type = float, 
		help="Min value for meshes.", default = 0.0)
	parser.add_argument("--max", dest="max_val", type = float,
		help="Max value for meshes.", default = 1E22)
	parser.add_argument("--3d", dest="lthreed", action='store_true',
		help="Create 3D mesh", default = False)
	args = parser.parse_args()

	gmsh.initialize()
	gmsh.open(args.in_file)

	entities = gmsh.model.getEntities()
	num_entities = len(entities)
	print(f"Number of entities in the model: {num_entities}")

	# Iterate through entities and print their properties
	for entity in entities:
		entity_type = gmsh.model.getType(entity[0],entity[1])
		print(f"Entity Type:{entity_type}")

	gmsh.option.setNumber("Mesh.CharacteristicLengthMin", args.min_val)
	gmsh.option.setNumber("Mesh.CharacteristicLengthMax", args.max_val)

	gmsh.model.geo.synchronize()

	out_file = args.out_file

	if args.lthreed:
		out_file.replace('.stl','.msh')
		gmsh.model.mesh.generate(3)
	else:
		gmsh.model.mesh.generate(2)

	gmsh.write(out_file)

	gmsh.finalize()

	sys.exit(0)

