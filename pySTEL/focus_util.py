#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Main routine
if __name__=="__main__":
	import sys
	from argparse import ArgumentParser
	import numpy as np
	import matplotlib.pyplot as pyplot
	import vtk
	from libstell.vmec import VMEC
	from libstell.focus import FOCUS
	from libstell.bnorm import BNORM
	from libstell.coils import COILSET
	from libstell.plot3D import PLOT3D
	parser = ArgumentParser(description= 
		'''Provides class for plotting FOCUS simulation results.''')
	parser.add_argument("-f", "--focus", dest="focus_ext",
		help="FOCUS file extension", default = None)
	parser.add_argument("-v", "--vmec", dest="vmec_ext",
		help="VMEC file extension", default = None)
	parser.add_argument("--bnorm", dest="bnorm_ext",
		help="BNORM file extension", default = None)
	parser.add_argument("-p", "--plot", dest="lplot", action='store_true',
		help="Plot the FOCUS data.", default = False)
	parser.add_argument("--plot3d", dest="lplot3d", action='store_true',
		help="Plot the plams surface and coil in 3D", default = False)
	parser.add_argument("--plotcoildist", dest="lplotcoildist", action='store_true',
		help="Plot the coil-plasma distance.", default = False)
	parser.add_argument("--gensurf", dest="lgensurf", action='store_true',
		help="Generate plasma.boundary VMEC.", default = False)
	parser.add_argument("--limiter_dist", dest="lim_dist",
		help="Generate limiter.boundary at offset distance of lim_dist.", default = None)
	focus_data = FOCUS()
	coil_data=COILSET()
	args = parser.parse_args()
	if args.focus_ext:
		try:
			focus_data.read_focus_HDF5(args.focus_ext)
		except:
			print(f'Could not file FOCUS HDF5 file: {args.focus_ext}')
			sys.exit(-1)
		if args.lgensurf and args.vmec_ext:
			wout = VMEC()
			wout.read_wout(args.vmec_ext)
			rmns = None;	zmnc = None
			xm_b = None;	xn_b = None
			bmnc = None;	bmns = None
			if args.bnorm_ext:
				bnorm = BNORM()
				bnorm.read_bnorm(bnorm_filename)
				curpol = wout.getCurrentPoloidal()
				xn_b = -bnorm.xn
				xm_b =  bnorm.xm
				bnmnc = curpol*bnorm.bnmnc[0,:]
				bnmns = curpol*bnorm.bnmns[0,:]
			k = wout.ns-1
			xm = wout.xm[:,0]
			xn = -wout.xn[:,0]/wout.nfp
			rmnc = wout.rmnc[k,:]
			zmns = wout.zmns[k,:]
			if wout.iasym == 1:
				rmns = wout.rmns[k,:]
				zmnc = wout.zmnc[k,:]
			focus_data.write_focus_plasma(wout.nfp,xm,xn,rmnc,zmns,rmns=rmns,zmnc=zmnc,xm_b=xm_b,xn_b=xn_b,bmnc=bnmnc,bmns=bnmns,filename='plasma.boundary')
			if args.lim_dist:
				[rmnc,zmns,rmns,zmnc]=wout.fitSurface(dist=-dist)
				focus_data.write_focus_plasma(wout.nfp,xm,xn,rmnc,zmns,rmns=rmns,zmnc=zmnc,filename='limiter.boundary')
				# Make a plot
				px = 1/pyplot.rcParams['figure.dpi']
				fig=pyplot.figure(figsize=(1024*px,768*px))
				ax=fig.add_subplot(111)
				theta = np.ndarray((360,1))
				zeta  = np.ndarray((3,1))
				for j in range(360): theta[j]=2.0*np.pi*j/359.0
				for j in range(3):   zeta[j]=     np.pi*j/2.0
				r = wout.cfunct(theta,zeta,wout.rmnc,wout.xm,wout.xn/wout.nfp)
				z = wout.sfunct(theta,zeta,wout.zmns,wout.xm,wout.xn/wout.nfp)
				r2 = wout.cfunct(theta,zeta,np.broadcast_to(rmnc,(1,wout.mnmax)),wout.xm,wout.xn/wout.nfp)
				z2 = wout.sfunct(theta,zeta,np.broadcast_to(zmns,(1,wout.mnmax)),wout.xm,wout.xn/wout.nfp)
				ax.plot(r[1,1,0],z[1,1,0],'+r')
				ax.plot(r[1,1,1],z[1,1,1],'+g')
				ax.plot(r[1,1,2],z[1,1,2],'+b')
				j = wout.ns-1
				ax.plot(r[j,:,0],z[j,:,0],'r',label='Plasma')
				ax.plot(r[j,:,1],z[j,:,1],'g')
				ax.plot(r[j,:,2],z[j,:,2],'b')
				j = 0
				ax.plot(r2[j,:,0],z2[j,:,0],'--r',label='Limiter')
				ax.plot(r2[j,:,1],z2[j,:,1],'--g')
				ax.plot(r2[j,:,2],z2[j,:,2],'--b')
				ax.set_xlabel(r"R [m]$")
				ax.set_ylabel(r"Z [m]")
				ax.set_title("FOCUS Limiter Surface")
				ax.set_aspect('equal', adjustable='box')
				ax.legend()
				pyplot.show()
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
			plt3d = PLOT3D()
			focus_data.plotBN3D(plt3d)
			focus_data.plotLimiter3D(plt3d)
			try:
				coil_data.read_coils_file(args.focus_ext+'.coils')
				coil_data.plotcoilsHalfFP(plt3d)
			except:
				i=1
			plt3d.render()
		if args.lplotcoildist:
			coil_data.read_coils_file(args.focus_ext+'.coils')
			coil_data.coilSurfDist(focus_data.xsurf.flatten(),\
					focus_data.ysurf.flatten(),\
					focus_data.zsurf.flatten())
			coil_data.plotcoilplasmaDist()
