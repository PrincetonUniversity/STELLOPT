#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Main routine
if __name__=="__main__":
	import sys
	from argparse import ArgumentParser
	import matplotlib.pyplot as pyplot
	from libstell.vmec import VMEC
	from libstell.fieldlines import FIELDLINES
	from libstell.plot3D import PLOT3D
	import numpy as np
	parser = ArgumentParser(description= 
		'''Provides a tool for accessing fieldlines data.''')
	parser.add_argument("-f", "--fieldlines", dest="fieldlines_ext",
		help="FIELDLINES file extension", default = None)
	parser.add_argument("-p", "--plot", dest="lplot", action='store_true',
		help="Plot the fieldlines file.", default = False)
	parser.add_argument("--plot3d", dest="k3d",
		help="Plot a fieldline in 3D.", default = None, type=int)
	parser.add_argument("-v", "--vmec", dest="vmec_ext", 
		help="Add VMEC equilbrium to plot", default = None)
	parser.add_argument("--plotheat", dest="heatfactor",
		help="Plot the heatflux scaled to a total power in [W]", default = None, type=float)
	args = parser.parse_args()
	field_data = FIELDLINES()
	if args.fieldlines_ext:
		field_data.read_fieldlines('fieldlines_'+args.fieldlines_ext+'.h5')
		if args.lplot:
			px = 1/pyplot.rcParams['figure.dpi']
			fig,(ax1,ax2,ax3) = pyplot.subplots(1,3,sharey=True,figsize=(1024*px,512*px))
			pyplot.subplots_adjust(hspace=0.1,wspace=0.15)
			k = 0
			field_data.plot_poincare(k,6,ax=ax1)
			k = field_data.phiaxis[int(field_data.nphi/4)]
			field_data.plot_poincare(k,6,ax=ax2)
			k = field_data.phiaxis[int(field_data.nphi/2)]
			field_data.plot_poincare(k,6,ax=ax3)
			if args.vmec_ext:
				vmec_wout = VMEC()
				vmec_wout.read_wout(args.vmec_ext)
				theta = np.ndarray((360,1))
				zeta  = np.ndarray((3,1))
				for j in range(360): theta[j]=2.0*np.pi*j/359.0
				for j in range(3):   zeta[j]=     np.pi*j/2.0
				r = vmec_wout.cfunct(theta,zeta,vmec_wout.rmnc,vmec_wout.xm,vmec_wout.xn/vmec_wout.nfp)
				z = vmec_wout.sfunct(theta,zeta,vmec_wout.zmns,vmec_wout.xm,vmec_wout.xn/vmec_wout.nfp)
				j = vmec_wout.ns-1
				ax1.plot(r[j,:,0],z[j,:,0],'r')
				ax2.plot(r[j,:,1],z[j,:,1],'r')
				ax3.plot(r[j,:,2],z[j,:,2],'r')
			pyplot.show()
		if args.k3d:
			plt3d = PLOT3D()
			field_data.plot_cloud(args.k3d,plot3D=plt3d,pointsize=0.1)
			plt3d.render()
		if args.heatfactor:
			plt3d = PLOT3D()
			fact = args.heatfactor/field_data.nlines
			field_data.plot_heatflux(factor=args.heatfactor/field_data.nlines,colormap='hot',plot3D=plt3d)
			plt3d.colorbar(title=rf'Q [W/$m^2$]')
			plt3d.render()
	sys.exit(0)