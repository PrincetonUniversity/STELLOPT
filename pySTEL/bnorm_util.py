#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Main routine
if __name__=="__main__":
	import sys
	from argparse import ArgumentParser
	import numpy as np
	import matplotlib.pyplot as pyplot
	from libstell.boozer import BOOZER
	from libstell.vmec import VMEC
	from libstell.libstell import FourierRep
	parser = ArgumentParser(description= 
		'''Utility for computing the normal magnetic field of a plasma.
		Similar in function to the BNORM code.''')
	parser.add_argument("-b", "--boozer", dest="booz_ext",
		help="BOOZER boozmn file extension", default = None)
	parser.add_argument("-v", "--vmec", dest="vmec_ext",
		help="VMEC wout file extension", default = None)
	parser.add_argument("--ntheta", dest="ntheta",
		help="Number of poloidal gridpoints to use.", default = 256)
	parser.add_argument("--nzeta", dest="nzeta",
		help="Number of toroidal gridpoints to use.", default = 91)
	args = parser.parse_args()
	fourier = FourierRep()
	pi2 = 2.0*np.pi
	# First process input
	if args.vmec_ext:
		data = VMEC()
		data.read_wout(args.vmec_ext)
		nfp = data.nfp
		lasym = data.lasym
		curpol = data.getCurrentPoloidal()
		print(curpol)
		xm       = data.xm
		xn       = data.xn/nfp
		xm_nyq   = data.xm_nyq
		xn_nyq   = data.xn_nyq/nfp
		rmnc     = data.rmnc
		zmns     = data.zmns
		bsubumnc = pi2*data.bsubumnc/curpol
		bsubvmnc = pi2*data.bsubvmnc/(curpol*nfp)
		pmns     = np.zeros_like(rmnc)
		if lasym:
			rmns     = data.rmns
			zmnc     = data.zmnc
			pmnc     = np.zeros_like(rmns)
			bsubumns = data.bsubumns
			bsubvmns = data.bsubvmns
	elif args.booz_ext:
		data = BOOZER()
		data.read_boozer(args.booz_ext)
		nfp = data.nfp_b
		lasym = data.lasym_b
		curpol = data.getCurrentPoloidal()
		xm       = data.ixm_b
		xn       = data.ixn_b/nfp
		xm_nyq   = data.ixm_b
		xn_nyq   = data.ixn_b/nfp
		rmnc     = data.rmnc_b
		zmns     = data.zmns_b
		pmns     = data.pmns_b
		bsubumnc = np.zeros_like(rmnc)
		bsubvmnc = np.zeros_like(rmnc)
		bsubumnc[:,mn00] = data.getCurrentToroidal()/pi2
		bsubvmnc[:,mn00] = data.getCurrentPoloidal()/pi2
		if lasym:
			rmns     = data.rmns_b
			zmnc     = data.zmnc_b
			pmnc     = data.pmnc_b
			bsubumns = np.zeros_like(rmnc)
			bsubvmns = np.zeros_like(rmnc)
	if (lasym):
		print('Not fully Implemented')
		sys.exit(0)
	# Realspace arrays
	ntheta = args.ntheta
	nzeta  = args.nzeta
	theta = np.deg2rad(np.linspace([0],[360],ntheta+1))
	zeta = np.deg2rad(np.linspace([0],[360],nzeta+1))
	theta = theta[:-1]
	zeta = zeta[:-1]
	# Derivatives (note both vmec and boozer classes in mu+nv)
	rumns = np.zeros_like(rmnc)
	rvmns = np.zeros_like(rmnc)
	zumnc = np.zeros_like(rmnc)
	zvmnc = np.zeros_like(rmnc)
	pumnc = np.zeros_like(rmnc)
	pvmnc = np.zeros_like(rmnc)
	for mn in range(rmnc.shape[1]):
		rumns[:,mn] = -xm[mn]*rmnc[:,mn]
		rvmns[:,mn] = -xn[mn]*rmnc[:,mn]
		zumnc[:,mn] =  xm[mn]*zmns[:,mn]
		zvmnc[:,mn] =  xn[mn]*zmns[:,mn]
		pumnc[:,mn] =  xm[mn]*pmns[:,mn]
		pvmnc[:,mn] =  xn[mn]*pmns[:,mn]
	# Fourier transform
	r  = fourier.cfunct(theta,zeta,rmnc,xm,xn)
	z  = fourier.sfunct(theta,zeta,zmns,xm,xn)
	p  = fourier.sfunct(theta,zeta,pmns,xm,xn)
	ru = fourier.sfunct(theta,zeta,rumns,xm,xn)
	rv = fourier.sfunct(theta,zeta,rvmns,xm,xn)
	zu = fourier.cfunct(theta,zeta,zumnc,xm,xn)
	zv = fourier.cfunct(theta,zeta,zvmnc,xm,xn)
	pu = fourier.cfunct(theta,zeta,pumnc,xm,xn)
	pv = fourier.cfunct(theta,zeta,pvmnc,xm,xn)
	bu = fourier.cfunct(theta,zeta,bsubumnc,xm_nyq,xn_nyq)
	bv = fourier.cfunct(theta,zeta,bsubvmnc,xm_nyq,xn_nyq)
	# Integral helper
	indu = np.zeros((ntheta,nzeta))
	indv = np.zeros((ntheta,nzeta))
	for u in range(ntheta): indu[u,:] = u
	for v in range(nzeta):  indv[:,v] = v
	# Get the toroidal angle array
	zeta0 = np.zeros_like(r)
	for v in range(r.shape[2]):
		zeta0[:,:,v] = p[:,:,v] + zeta[v]
	# We want just the edge values
	r = r[-1,:,:].flatten(); z = z[-1,:,:].flatten(); p = p[-1,:,:].flatten()
	ru = ru[-1,:,:].flatten(); zu = zu[-1,:,:].flatten(); pu = pu[-1,:,:].flatten()
	rv = rv[-1,:,:].flatten(); zv = zv[-1,:,:].flatten(); pv = pv[-1,:,:].flatten()
	bu = bu[-1,:,:].flatten(); bv = bv[-1,:,:].flatten()
	indu = indu.flatten(); indv = indv.flatten()
	zeta0 = zeta0[-1,:,:].flatten()
	# Get x y
	x = r * np.cos(zeta0)
	y = r * np.sin(zeta0)
	xu = ru * np.cos(zeta0) - r * np.sin(zeta0) * pu
	yu = ru * np.sin(zeta0) + r * np.cos(zeta0) * pu
	xv = rv * np.cos(zeta0) - pi2 * y / nfp - r * np.sin(zeta0) * pv
	yv = rv * np.sin(zeta0) + pi2 * x / nfp + r * np.cos(zeta0) * pv
	# Normals
	snx = yu * zv - zu * yv
	sny = zu * xv - xu * zv
	snz = xu * yv - yu * xv
	sqf = np.sqrt(snx * snx + sny * sny + snz * snz)
	#print(2*np.pi*sum(sqf)/(ntheta*nzeta))
	guu = xu * xu + yu * yu + zu * zu
	guv = xu * xv + yu * yv + zu * zv
	gvv = xv * xv + yv * yv + zv * zv
	dju = bu * guv + bv * guu
	djv = bu * gvv + bv * guv
	djx = bu * xv - bv * xu
	djy = bu * yv - bv * yu
	djz = bu * zv - bv * zu
	# Correct guv and gvv
	guv = guv * nfp
	gvv = gvv * nfp * nfp
	# Compute vector potential (non-sigular part)
	nuv = ntheta*nzeta
	ax  = np.zeros_like(x)
	ay  = np.zeros_like(x)
	az  = np.zeros_like(x)
	for i in range(1,nuv-1):
		# Up to i
		i1 = 0
		i2 = i-1
		dx = x[i] - x[i1:i2]
		dy = y[i] - y[i1:i2]
		dz = z[i] - z[i1:i2]
		sq = 1.0/np.sqrt(dx * dx + dy * dy + dz * dz)
		th1 = np.pi * (indu[i1:i2] - indu[i])/ntheta
		zt1 = np.pi * (indv[i1:i2] - indv[i])/(nzeta*nfp)
		tu = np.tan(th1)/np.pi
		tv = np.tan(zt1)/np.pi
		sqs = 1.0 / np.sqrt( guu[i]*tu*tu + 2.0 * guv[i] * tu * tv + gvv[i] * tv * tv)
		sqsum = sum(sqs)
		ax = ax + sum(djx[i1:i2]*sq) - djx[i] * sqsum
		ay = ay + sum(djy[i1:i2]*sq) - djy[i] * sqsum
		az = az + sum(djz[i1:i2]*sq) - djz[i] * sqsum
		# From i
		i1 = i+1
		i2 = ntheta*nzeta
		dx = x[i] - x[i1:i2]
		dy = y[i] - y[i1:i2]
		dz = z[i] - z[i1:i2]
		sq = 1.0/np.sqrt(dx * dx + dy * dy + dz * dz)
		th1 = np.pi * (indu[i1:i2] - indu[i])/ntheta
		zt1 = np.pi * (indv[i1:i2] - indv[i])/(nzeta*nfp)
		tu = np.tan(th1)/np.pi
		tv = np.tan(zt1)/np.pi
		sqs = 1.0 / np.sqrt( guu[i]*tu*tu + 2.0 * guv[i] * tu * tv + gvv[i] * tv * tv)
		sqsum = sum(sqs)
		ax = ax + sum(djx[i1:i2]*sq) - djx[i] * sqsum
		ay = ay + sum(djy[i1:i2]*sq) - djy[i] * sqsum
		az = az + sum(djz[i1:i2]*sq) - djz[i] * sqsum
	# Now do the regular integral
	Integral = np.zeros_like(ax)
	for i in range(nuv):
		sqp = np.sqrt( guu[i] + 2.0 * guv[i] + gvv [i] )
		sqm = np.sqrt( guu[i] - 2.0 * guv[i] + gvv [i] )
		sqa = np.sqrt( guu[i] )
		sqc = np.sqrt( gvv[i] )
		top = np.log( ( sqc * sqp + gvv[i] + guv[i] ) / ( sqa * sqp - guu[i] - guv[i] ) ) / sqp
		tom = np.log( ( sqc * sqm + gvv[i] - guv[i] ) / ( sqa * sqm - guu[i] + guv[i] ) ) / sqm
		Integral[i] = top + tom
	# Now compute AU and AV
	au = np.zeros_like(ax)
	av = np.zeros_like(ax)
	for i in range(nuv):
		dintu = ( ax[i] * xu[i] + ay[i] * yu[i] + az[i] * zu[i] ) / nuv
		dintv = ( ax[i] * xv[i] + ay[i] * yv[i] + az[i] * zv[i] ) / nuv
		au[i] = ( dintu + dju[i] * Integral[i] * nfp ) / (pi2*2.0)
		av[i] = ( dintv + djv[i] * Integral[i] * nfp ) / (pi2*2.0)
	# Fouier transform AU
	mpol = int(max(xm))
	ntor = int(max(abs(xn)))
	aumnc = np.zeros((mpol+1,2*ntor+1))
	avmnc = np.zeros((mpol+1,2*ntor+1))
	faz   = np.ones((mpol+1)); faz[0] = 2.0; faz = faz/nuv
	for m in range(mpol+1):
		for n in range(-ntor,ntor+1):
			ndex = n + ntor
			for i in range(nuv):
				temp = indu[i] / ntheta + indv[i] / nzeta
				aumnc[m,ndex] = aumnc[m,ndex] + au[i] * np.cos(pi2*temp)*faz[m]
				avmnc[m,ndex] = avmnc[m,ndex] + av[i] * np.cos(pi2*temp)*faz[m]
	print(aumnc)
	print(avmnc)
	# Compute B-normal
	bn = np.zeros((nuv))
	for m in range(mpol+1):
		for n in range(-ntor,ntor+1):
			ndex = n + ntor
			for i in range(nuv):
				temp = indu[i] / ntheta + indv[i] / nzeta
				bn[i] = bn[i] + pi2 * ( m * avmnc[m,ndex] - n * aumnc[m,ndex] ) * np.sin(pi2*temp)
	# Fouier transform BN
	bnmns = np.zeros((mpol+1,2*ntor+1))
	for m in range(mpol+1):
		for n in range(-ntor,ntor+1):
			ndex = n + ntor
			for i in range(nuv):
				temp = indu[i] / ntheta + indv[i] / nzeta
				bnmns[m,ndex] = bnmns[m,ndex] + bn[i] * np.sin(pi2*temp) * faz[m] / sqf[i]
			print(m,n,bnmns[m,ndex])
	sys.exit(0)