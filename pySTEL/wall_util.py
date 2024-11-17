#!/usr/bin/env python3
# -*- coding: utf-8 -*-


if __name__=="__main__":
	import sys
	from argparse import ArgumentParser
	from libstell.wall import WALL
	from libstell.plot3D import PLOT3D
	import numpy as np
	from datetime import datetime
	from stl import mesh
	parser = ArgumentParser(description= 
		'''Provides tool for plotting wall files and working with them.''')
	parser.add_argument("-w", "--wall", dest="wall_file",
		help="Wall file for input", default = None)
	parser.add_argument("-p", "--plot", dest="lplot", action='store_true',
		help="Plot the wall file.", default = False)
	parser.add_argument("--stl", dest="lstl", action='store_true',
		help="Generate STL of the wall.", default = False)
	parser.add_argument("--clean", dest="lclean", action='store_true',
		help="Generated a cleaned version of the wall file.", default = False)
	parser.add_argument("--refine", dest="refine", nargs=2, metavar=('dlmin', 'dlmax'),
		help="Refine the mesh using dlmin, dlmax", default = [None,None], type=float)
	parser.add_argument("--kisslinger", dest="lkisslinger", action='store_true',
		help="Write the kisslinger version of the wall.", default = False)
	parser.add_argument("--nphi", dest="nphi",
		help="Number of toroidal points to use (when applicable)", default = 180, type=int)
	parser.add_argument("--nfp", dest="nfp",
		help="Number of field periods to use (when applicable)", default = 1, type=int)
	args = parser.parse_args()
	wall = WALL()
	if args.wall_file: 
		wall.read_wall(args.wall_file)
		temp = args.wall_file.split('/')
		fileout = temp[-1]
		# Any wall modification comes here
		if args.lclean: 
			wall.wallClean()
			fileout=fileout.replace('.dat','_clean.dat')
		if args.refine[0]: 
			wall.refineWall(args.refine[0],args.refine[1])
			fileout=fileout.replace('.dat','_refine.dat')
		print(fileout)
		# Plots go here
		if args.lplot: wall.plot_wall_3D()
		# Outputting to other grids goes here.
		if args.refine[0]: wall.write_wall(fileout.replace('.dat','_new.dat'))
		if args.lstl: wall.write_wall_stl(fileout.replace('.dat','.stl'))
		if args.lkisslinger: wall.write_wall_kisslinger(fileout.replace('.dat','.kis'),args.nphi,args.nfp)
	sys.exit(0)

