#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Main routine
if __name__=="__main__":
	import sys
	from argparse import ArgumentParser
	from libstell.mumat import MUMAT
	from libstell.plot3D import PLOT3D
	import numpy as np
	import gmsh
	parser = ArgumentParser(description= 
		'''Provides a tool creating and plotting mumaterail data.''')
	parser.add_argument("-m", "--mumat", dest="mumat_file",
		help="Mumaterial file name", default = None)
	parser.add_argument("--plot", dest="lplot", action='store_true',
		help="Plot the mumat data in 3D.", default = False)
	parser.add_argument("--plot_state", dest="lplot_state", action='store_true',
		help="Plot the mumat state functions.", default = False)
	parser.add_argument("--step", dest="step_file",
		help="STEP CAD model file name", default = None)
	parser.add_argument("--min", dest="min_dist", type = float, 
		help="Min value for meshes. [m]", default = 0.05)
	parser.add_argument("--max", dest="max_dist", type = float,
		help="Max value for meshes. [m]", default = 0.10)
	parser.add_argument("--mm2m", dest="lmm2m", action='store_true',
		help="Convert STEP file from mm to m.", default = False)
	parser.add_argument("--save", dest="lsave", action='store_true',
		help="Save the plots with ext names.", default = False)
	args = parser.parse_args()
	mumat_data = MUMAT()
	if args.mumat_file:
		out_name = args.mumat_file.replace('.dat','')
		mumat_data.read_mumat(args.mumat_file)
	if args.step_file:
		out_name = args.step_file.replace('.step','').replace('.stp','')
		out_file = args.step_file.replace('.step','.dat').replace('.stp','.dat')
		gmsh.initialize()
		gmsh.open(args.step_file)
		scale = 1.0
		if args.lmm2m: scale = 1.0E-3
		gmsh.option.setNumber("Mesh.CharacteristicLengthMin", args.min_dist/scale)
		gmsh.option.setNumber("Mesh.CharacteristicLengthMax", args.max_dist/scale)
		gmsh.model.geo.synchronize()
		gmsh.model.mesh.generate(3)
		mumat_data.load_gmsh(gmsh)
		mumat_data.vertex = mumat_data.vertex*scale
		state_arr = np.array([[232.07003497538471, 0.04872389791183283],
							[406.831659753871, 0.10232018561484922],
							[649.4138669272695, 0.20220417633410695],
							[855.4771124381402, 0.29965197215777256],
							[1046.669744246662, 0.3970997679814387],
							[1243.295150667669, 0.496983758700696],
							[1462.3816681297542, 0.5968677494199535],
							[1686.4991091042011, 0.7016241299303941],
							[1916.491039470617, 0.7941995359628773],
							[2177.7977464857613, 0.8965197215777262],
							[2523.9859109091462, 0.9964037122969841],
							[3118.629589691876, 1.0962877030162412],
							[3988.5853433458387, 1.1912993039443156],
							[4832.08588970385, 1.2960556844547564],
							[5739.7996458882735, 1.3983758700696056]]) # Grade 91 (maybe)
		H = state_arr[:,0] 
		M = state_arr[:,1] 
		mumat_data.add_state(2,H=H,M=M)
		tetdex = list(range(0,mumat_data.ntet))
		mumat_data.set_state(1,tetdex)
		mumat_data.write_mumat(out_file)
	if args.lplot:
		plt3d = PLOT3D()
		mumat_data.plot_mesh(plot3D=plt3d)
		plt3d.render()
		if (args.lsave): plt3d.saveImage(f'mumat_mesh3D_{out_name}.png')
	if args.lplot_state:
		px = 1/pyplot.rcParams['figure.dpi']
		fig=pyplot.figure(figsize=(1024*px,768*px))
		ax=fig.add_subplot(111)
		mumat_data.plot_state()
		
	sys.exit(0)