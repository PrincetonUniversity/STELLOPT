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
	from libstell.bnorm import BNORM
	#from libstell.plot3D import PLOT3D
	parser = ArgumentParser(description= 
		'''Utility for generating FOCUS boundary and limiter information''')
	parser.add_argument("-f", "--focus", dest="focus_ext",
		help="FOCUS file extension", default = None)
	parser.add_argument("-v", "--vmec", dest="vmec_ext",
		help="VMEC file extension", default = None)
	parser.add_argument("--bnorm", dest="bnorm_ext",
		help="BNORM file extension", default = None)
	parser.add_argument("--dist", dest="dist",
		help="Distance in m to extrapolate.", default = 0.0, type=float)
	parser.add_argument("--fitsurf", dest="lfitsurf", action='store_true',
		help="Fit a surface (default is to extrapolate)", default = False)
	parser.add_argument("--width", dest="width",
		help="Island width in [m].", default = 0.1, type=float)
	parser.add_argument("-m","--misland", dest="misland",
		help="Island poloidal mode.", default = 4, type=int)
	parser.add_argument("-n","--nisland", dest="nisland",
		help="Island toroidal mode.", default = 4, type=int)
	parser.add_argument("--plot", dest="lplot", action='store_true',
		help="Make a 2D plot.", default = False)
	focus = FOCUS()
	wout  = VMEC()
	bnorm = BNORM()
	args = parser.parse_args()
	# Read VMEC
	wout = VMEC()
	wout.read_wout(args.vmec_ext)
	curpol = wout.getCurrentPoloidal()
	# Bnorm
	if (args.bnorm_ext):
		print(' Using BNORM Field on plasma boundary.')
		print(f'    VMEC Normalization I(poloidal): {curpol:10.3e}')
		bnorm.read_bnorm(args.bnorm_ext)
		xm_b = np.squeeze(bnorm.xm)
		xn_b = np.squeeze(bnorm.xn)
		bmnc = np.squeeze(bnorm.bnmnc*curpol)
		bmns = np.squeeze(bnorm.bnmns*curpol)
	else:
		xm_b = np.array([0])
		xn_b = np.array([0])
		bmnc = np.array([0])
		bmns = np.array([0])
	# Write plasma boundary
	focus.write_focus_plasma(wout.nfp,np.squeeze(wout.xm),-np.squeeze(wout.xn)/wout.nfp, 
		np.squeeze(wout.rmnc[wout.ns-1,:]),np.squeeze(wout.zmns[wout.ns-1,:]),
		xm_b = xm_b, xn_b = xn_b, bmnc = bmnc, bmns = bmns)
	if (abs(args.dist) > 0):
		# Extrapolate or
		if args.lfitsurf:
			print(' Fitting surface')
			[rmnc,zmns,rmns,zmnc]=wout.fitSurface(dist=-args.dist)
		else:
			print(' Extrapolating surface')
			[rmnc,zmns,rmns,zmnc]=wout.extrapSurface(dist=args.dist)
		# Calculte Bmn from island width w=sqrt(R0*Bmn/(m*B0*iota))
		Bnorm = np.squeeze(args.width**2 * args.misland * wout.b0 * wout.iotaf[-1] / wout.rmajor)
		xm_b = np.array([args.misland])
		xn_b = np.array([args.nisland])/wout.nfp
		bmnc = np.array([0])
		bmns = np.array([Bnorm])
		print(f'    Surface Distance: {args.dist:5.3f} [m]')
		print(f'        Island Width: {args.width:5.3f} [m]')
		print(f'          Island n/m: {args.nisland:d}/{args.misland:d}')
		print(f'        Island Bnorm: {Bnorm:5.3f} [T]')
		print(f'             VMEC B0: {wout.b0:5.3f} [T]')
		print(f'     VMEC iota(edge): {np.squeeze(wout.iotaf[-1]):5.3f}')
		print(f'         VMEC Rmajor: {wout.rmajor:5.3f} [m]')
		focus.write_focus_plasma(wout.nfp,np.squeeze(wout.xm),-np.squeeze(wout.xn)/wout.nfp,
			np.squeeze(rmnc),np.squeeze(zmns),
			xm_b = xm_b, xn_b = xn_b, bmnc = bmnc, bmns = bmns, filename='limiter.boundary')
	# plots
	if args.lplot:
		px = 1/pyplot.rcParams['figure.dpi']
		fig=pyplot.figure(figsize=(1024*px,768*px))
		ax=fig.add_subplot(111)
		theta = np.ndarray((360,1))
		zeta  = np.ndarray((3,1))
		for j in range(360): theta[j]=2.0*np.pi*j/359.0
		for j in range(3):   zeta[j]=     np.pi*j/2.0
		r = wout.cfunct(theta,zeta,wout.rmnc,wout.xm,wout.xn/wout.nfp)
		z = wout.sfunct(theta,zeta,wout.zmns,wout.xm,wout.xn/wout.nfp)
		ax.plot(r[1,1,0],z[1,1,0],'+r')
		ax.plot(r[1,1,1],z[1,1,1],'+g')
		ax.plot(r[1,1,2],z[1,1,2],'+b')
		j = wout.ns-1
		ax.plot(r[j,:,0],z[j,:,0],'r')
		ax.plot(r[j,:,1],z[j,:,1],'g')
		ax.plot(r[j,:,2],z[j,:,2],'b')
		ax.set_xlabel('R [m]')
		ax.set_ylabel('Z [m]')
		if abs(args.dist) > 0:
			r2 = wout.cfunct(theta,zeta,np.broadcast_to(rmnc,(1,wout.mnmax)),wout.xm,wout.xn/wout.nfp)
			z2 = wout.sfunct(theta,zeta,np.broadcast_to(zmns,(1,wout.mnmax)),wout.xm,wout.xn/wout.nfp)
			j = 0
			ax.plot(r2[j,:,0],z2[j,:,0],'--r')
			ax.plot(r2[j,:,1],z2[j,:,1],'--g')
			ax.plot(r2[j,:,2],z2[j,:,2],'--b')
			ax.set_title(f'Extrapolated Surface ({abs(args.dist)/wout.aminor:4.2f} Aminor)')
		pyplot.show()
