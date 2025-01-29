#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Main routine
if __name__=="__main__":
	import sys
	from argparse import ArgumentParser
	import numpy as np
	import matplotlib.pyplot as pyplot
	from libstell.bootsj import BOOTSJ
	parser = ArgumentParser(description= 
		'''serves as a simple tool for assessing bootsj output files.''')
	parser.add_argument("-b", "--bootsj", dest="bootsj_ext",
		help="BOOTSJ file extension", default = None)
	parser.add_argument("-p", "--plot", dest="lplot", action='store_true',
		help="Plot the BOOTSJ output.", default = False)
	bootsj_data = BOOTSJ()
	args = parser.parse_args()
	if args.bootsj_ext:
		try:
			bootsj_data.read_answers_plot(rf'answers_plot.{args.bootsj_ext}')
		except:
			print(f'Could not file answers_plot file: {args.booz_ext}')
			sys.exit(-1)
		if args.lplot:
			px = 1/pyplot.rcParams['figure.dpi']
			fig,(ax1,ax2) = pyplot.subplots(1,2,figsize=(1024*px,768*px))
			pyplot.subplots_adjust(hspace=0.1,wspace=0.3)
			ax1.plot(bootsj_data.rhoar,bootsj_data.tempe1,'b',label=rf'$T_e$')
			ax1.plot(bootsj_data.rhoar,bootsj_data.tempi1,'r',label=rf'$T_i$')
			ax1.plot(bootsj_data.rhoar,bootsj_data.dense*10.0,'b--',label=rf'$n_e$')
			ax1.plot(bootsj_data.rhoar,bootsj_data.densi*10.0,'r--',label=rf'$n_i$')
			ax1.set_xlabel('Norm. Toroidal Flux (s)')
			ax1.set_ylabel(rf'T [keV]; n $x10^{19}$ [$m^{-3}$]')
			ax2.plot(bootsj_data.rhoar,bootsj_data.dibs)
			ax2.set_xlabel('Norm. Toroidal Flux (s)')
			ax2.set_ylabel(rf'dI/ds [A]')
			ax2.text(0.02,0.05,rf'I = {bootsj_data.Itotal*1.0E-6:4.3f} [MA]', horizontalalignment='left',\
				verticalalignment='center', transform=ax2.transAxes)
			pyplot.show()
	sys.exit(0)