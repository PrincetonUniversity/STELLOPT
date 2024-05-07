#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Main routine
if __name__=="__main__":
	import sys
	from argparse import ArgumentParser
	import numpy as np
	import matplotlib.pyplot as pyplot
	from libstell.boozer import BOOZER
	parser = ArgumentParser(description= 
		'''Provides class for accessing boozer data also servers as a
		   simple tool for assessing boozer boozmn files.''')
	parser.add_argument("-b", "--boozer", dest="booz_ext",
		help="BOOZER boozmn file extension", default = None)
	parser.add_argument("-p", "--plot", dest="lplot", action='store_true',
		help="Plot the boozer file.", default = False)
	parser.add_argument("-s", "--surf", dest="sval",
		help="Equilibrium surface to plot.", default = None)
	parser.add_argument("--bootsj", dest="lbootsj", action='store_true',
		help="Output a bootsj file.", default = False)
	boozmn_data = BOOZER()
	args = parser.parse_args()
	if args.booz_ext:
		try:
			boozmn_data.read_boozer(args.booz_ext)
		except:
			print(f'Could not file boozmn file: {args.booz_ext}')
			sys.exit(-1)
		if args.lplot:
			if not args.sval:
				sval = int(0.25*boozmn_data.ns_b)
			else:
				sval = int(args.sval)-1
			px = 1/pyplot.rcParams['figure.dpi']
			fig,(ax1,ax2) = pyplot.subplots(1,2,figsize=(1024*px,768*px))
			pyplot.subplots_adjust(hspace=0.1,wspace=0.3)
			boozmn_data.plotBmnSpectrum(sval,ax=ax1)
			boozmn_data.plotBsurf(sval,ax=ax2)
			pyplot.show()
		if args.lbootsj:
			filename = 'in_bootsj.'+args.booz_ext
			f = open(filename,'w')
			f.write(f'{args.booz_ext}\n')
			for i in range(1,boozmn_data.ns_b):
				f.write(f' {i+1}')
			f.write('\n')
			f.close()
			# Setup bootin dict
			bootin_dict={'mbuse':int(max(np.squeeze(boozmn_data.ixm_b))),\
				'nbuse':int(max(np.squeeze(boozmn_data.ixn_b))/boozmn_data.nfp_b),\
				'teti':1.0,'zeff1':1.0,'dens0':1.0,'tempres':1.0,\
				'damp_bs':0.001,\
				'ate':np.array([1.0,-1.0,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]),\
				'ati':np.array([1.0,-1.0,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]),\
				}
			boozmn_data.libStell.set_bootin(bootin_dict)
			boozmn_data.libStell.write_bootin('input.'+args.booz_ext)
	sys.exit(0)