#!/usr/bin/env python3
# -*- coding: utf-8 -*-


if __name__=="__main__":
	import sys
	from argparse import ArgumentParser
	from libstell.nescoil import NESCOIL
	from libstell.bnorm import BNORM
	import matplotlib.pyplot as pyplot
	import numpy as np
	from datetime import datetime
	parser = ArgumentParser(description= 
		'''Provides class for accessing nescoil files''')
	parser.add_argument("--output", dest="nescout_file",
		help="NESCOIL output file", default = None)
	parser.add_argument("-c", "--cut_coils", dest="lcut_coils", action='store_true',
		help="Cut modular coils from the potential", default = False)
	parser.add_argument("--ncoil", dest="ncoil",
		help="Number of coils per field period", default = 5, type=int)
	parser.add_argument("-ch", "--cut_helical_coils", dest="lcut_helical_coils", action='store_true',
		help="Cut helical coils from the potential", default = False)
	parser.add_argument("-p", "--plot", dest="lplot", action='store_true',
		help="Make 2D plots.", default = False)
	parser.add_argument("-p3d", "--plot_3d", dest="lplot_3d", action='store_true',
		help="Make 3D plots.", default = False)
	args = parser.parse_args()
	nescout = NESCOIL()
	if args.nescout_file: 
		nescout.read_nescout(args.nescout_file)
		if args.lcut_coils:
			coil = nescout.cutcoils(args.ncoil,lplot=args.lplot,npts=256)
			coil_txt = args.nescout_file.split('.',1)
			coil.rescalecoils(256)
			coil.write_coils_file(f'coils.{coil_txt[1]}')
			if args.lplot_3d: coil.plotcoilsHalfFP()
		elif args.lcut_helical_coils:
			coil = nescout.cutcoils_helical(2,lplot=args.lplot)
			coil_txt = args.nescout_file.split('.',1)
			coil.rescalecoils(128)
			coil.write_coils_file(f'coils.{coil_txt[1]}')
			if args.lplot_3d: coil.plotcoils()
		else:
			if args.lplot: 
				px = 1/pyplot.rcParams['figure.dpi']
				fig=pyplot.figure(figsize=(1024*px,768*px))
				ax=fig.add_subplot(121)
				pyplot.subplots_adjust(hspace=0.4,wspace=0.3)
				nescout.plotpotential(ax=ax)
				ax=fig.add_subplot(122)
				nescout.plottotalpotential(ax=ax)
				pyplot.show()
			if args.lplot_3d: nescout.plotsurfaces()

	sys.exit(0)
