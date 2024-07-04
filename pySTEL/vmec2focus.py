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
	parser.add_argument("--fitsurf", dest="lfitsurf",
		help="Fit a surface (default is to extrapolate)", default = False)
	parser.add_argument("--width", dest="width",
		help="Island width in [m].", default = 0.1, type=float)
	parser.add_argument("-m","--misland", dest="misland",
		help="Island poloidal mode.", default = 4, type=int)
	parser.add_argument("-n","--nisland", dest="nisland",
		help="Island toroidal mode.", default = 4, type=int)
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
	focus.write_focus_plasma(wout.nfp,-np.squeeze(wout.xn)/wout.nfp,np.squeeze(wout.xm), 
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
		xn_b = np.array([args.nisland])
		bmnc = np.array([0])
		bmns = np.array([Bnorm])
		print(f'    Surface Distance: {args.dist:5.3f} [m]')
		print(f'        Island Width: {args.width:5.3f} [m]')
		print(f'          Island n/m: {args.nisland:d}/{args.misland:d}')
		print(f'        Island Bnorm: {Bnorm:5.3f} [T]')
		print(f'             VMEC B0: {wout.b0:5.3f} [T]')
		print(f'     VMEC iota(edge): {np.squeeze(wout.iotaf[-1]):5.3f}')
		print(f'         VMEC Rmajor: {wout.rmajor:5.3f} [m]')
		focus.write_focus_plasma(wout.nfp,-np.squeeze(wout.xn)/wout.nfp,np.squeeze(wout.xm), 
			np.squeeze(rmnc),np.squeeze(zmns),
			xm_b = xm_b, xn_b = xn_b, bmnc = bmnc, bmns = bmns, filename='limiter.boundary')
