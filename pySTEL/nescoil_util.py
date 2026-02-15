#!/usr/bin/env python3
# -*- coding: utf-8 -*-


if __name__=="__main__":
	import sys
	from argparse import ArgumentParser
	from libstell.nescoil import NESCOIL
	from libstell.bnorm import BNORM
	from libstell.plot3D import PLOT3D
	import matplotlib.pyplot as pyplot
	import numpy as np
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
	parser.add_argument("--HUD", dest="lhud", action='store_true',
		help="Add camera HUD to 3D plots.", default = False)
	parser.add_argument("--save", dest="lsave", action='store_true',
		help="Save the plots with ext names.", default = False)
	args = parser.parse_args()
	nescout = NESCOIL()
	if args.nescout_file: 
		nescout.read_nescout(args.nescout_file)
		ext_txt = args.nescout_file.split('.',1)[1]
		if args.lcut_coils:
			coil = nescout.cutcoils(args.ncoil,lplot=args.lplot,npts=256)
			coil_txt = args.nescout_file.split('.',1)
			coil.rescalecoils(256)
			coil.write_coils_file(f'coils.{ext_txt}')
			if args.lplot_3d: 
				plt3d = PLOT3D()
				coil.plotcoilsHalfFP(plot3D=plt3d)
				if args.lhud: plt3d.addCameraHUD()
				plt3d.setCamera(pos=[-2.543,-19.896,-2.512],focus=[8.068,-6.088,-1.409],camup=[0,0,1])
				plt3d.render()
				if (args.lsave): 
					plt3d.saveImage(f'nescoil_coils_{ext_txt}.png')
		elif args.lcut_helical_coils:
			coil = nescout.cutcoils_helical(2,lplot=args.lplot)
			coil_txt = args.nescout_file.split('.',1)
			coil.rescalecoils(128)
			coil.write_coils_file(f'coils.{ext_txt}')
			if args.lplot_3d: 
				plt3d = PLOT3D()
				coil.plotcoils(plot3D=plt3d)
				plt3d.render()
				if (args.lsave): plt3d.saveImage(f'nescoil_coils_{ext_txt}.png')
		else:
			if args.lplot: 
				px = 1/pyplot.rcParams['figure.dpi']
				fig=pyplot.figure(figsize=(1024*px,768*px))
				ax=fig.add_subplot(121)
				pyplot.subplots_adjust(hspace=0.4,wspace=0.3)
				nescout.plotpotential(ax=ax,cmap='Greens')
				ax=fig.add_subplot(122)
				nescout.plottotalpotential(ax=ax,cmap='Greens')
				pyplot.show()
				if (args.lsave): fig.savefig(f'nescoil_potential_{ext_txt}.png', dpi=fig.dpi)
			if args.lplot_3d: 
				plt3d = PLOT3D()
				nescout.plotsurfaces(plot3D=plt3d)
				if args.lhud: plt3d.addCameraHUD()
				plt3d.render()
				if (args.lsave): plt3d.saveImage(f'nescoil_surfaces_{ext_txt}.png')

	sys.exit(0)
