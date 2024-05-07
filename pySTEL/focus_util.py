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
			ax4.remove()
			ax4=fig.add_subplot(2,2,4,projection='3d')
			focus_data.plotBN3D(ax4)
			#try:
			coil_data.read_coils_file(args.focus_ext+'.coils')
			coil_data.plotcoilsHalfFP(ax4)
			#except:
			#	i=1
			pyplot.show()