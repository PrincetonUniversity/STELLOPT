#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Main routine
if __name__=="__main__":
	import sys
	from argparse import ArgumentParser
	import numpy as np
	import matplotlib.pyplot as pyplot
	from libstell.vmec import VMEC
	from libstell.focus import FOCUS
	from libstell.coils import COILSET
	parser = ArgumentParser(description= 
		'''Provides class for plotting FOCUS simulation results.''')
	parser.add_argument("-f", "--focus", dest="focus_ext",
		help="FOCUS file extension", default = None)
	parser.add_argument("-v", "--vmec", dest="vmec_ext",
		help="VMEC file extension", default = None)
	parser.add_argument("-p", "--plot", dest="lplot", action='store_true',
		help="Plot the FOCUS data.", default = False)
	parser.add_argument("--plot3d", dest="lplot3d", action='store_true',
		help="Plot the plams surface and coil in 3D", default = False)
	parser.add_argument("--plotcoildist", dest="lplotcoildist", action='store_true',
		help="Plot the coil-plasma distance.", default = False)
	focus_data = FOCUS()
	coil_data=COILSET()
	args = parser.parse_args()
	if args.focus_ext:
		try:
			focus_data.read_focus_HDF5(args.focus_ext)
		except:
			print(f'Could not file FOCUS HDF5 file: {args.focus_ext}')
			sys.exit(-1)
		if args.lplot:
			px = 1/pyplot.rcParams['figure.dpi']
			fig,((ax1,ax2),(ax3,ax4)) = pyplot.subplots(2,2,figsize=(1024*px,768*px))
			pyplot.subplots_adjust(hspace=0.3,wspace=0.3)
			focus_data.plotConvergence(ax1)
			focus_data.plotBNormal(ax2)
			focus_data.plotPoincare(ax3)
			focus_data.plotIota(ax4)
			pyplot.show()
		if args.lplot3d:
			px = 1/pyplot.rcParams['figure.dpi']
			fig=pyplot.figure(figsize=(1024*px,768*px))
			ax1=fig.add_subplot(111,projection='3d')
			focus_data.plotBN3D(ax1)
			try:
				coil_data.read_coils_file(args.focus_ext+'.coils')
				coil_data.plotcoilsHalfFP(ax1)
			except:
				i=1
			pyplot.show()
		if args.lplotcoildist:
			px = 1/pyplot.rcParams['figure.dpi']
			fig=pyplot.figure(figsize=(1024*px,768*px))
			ax1=fig.add_subplot(111,projection='3d')
			coil_data.read_coils_file(args.focus_ext+'.coils')
			coil_data.coilSurfDist(focus_data.xsurf.flatten(),\
					focus_data.ysurf.flatten(),\
					focus_data.zsurf.flatten())
			coil_data.plotcoilsDist(ax=ax1)
			ax1.set_axis_off()
			pyplot.show()
