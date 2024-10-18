#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for reading and handling VMEC
equilibrium data.
"""

# Libraries
from libstell.libstell import LIBSTELL, FourierRep

# Constants

# VMEC Class
class VMEC(FourierRep):
	"""Class for working with VMEC equilibria

	"""
	def __init__(self):
		super().__init__()
		self.nfp = None
		self.libStell = LIBSTELL()

	def read_wout(self,filename):
		"""Reads a VMEC wout_file

		This routine reads and initilizes the VMEC class
		with variable information from a VMEC wout file.

		Parameters
		----------
		file : str
			Path to wout file.
		"""
		import numpy as np
		wout_dict = self.libStell.read_wout(filename)
		for key in wout_dict:
			setattr(self, key, wout_dict[key])
		# (mu-nv) -> (mu+nv)
		self.xn = -self.xn
		self.xn_nyq = -self.xn_nyq
		# Half to full grid
		self.buco = self.h2f(self.buco)
		self.bvco = self.h2f(self.bvco)
		self.vp = self.h2f(self.vp)
		self.overr = self.h2f(self.overr)
		self.specw = self.h2f(self.specw)
		for mn in range(self.mnmax):
			self.lmns[:,mn] = self.h2f(self.lmns[:,mn])
		for mn in range(self.mnmax_nyq):
			self.bmnc[:,mn] = self.h2f(self.bmnc[:,mn])
			self.gmnc[:,mn] = self.h2f(self.gmnc[:,mn])
			self.bsupumnc[:,mn] = self.h2f(self.bsupumnc[:,mn])
			self.bsupvmnc[:,mn] = self.h2f(self.bsupvmnc[:,mn])
			self.bsubsmns[:,mn] = self.h2f(self.bsubsmns[:,mn])
			self.bsubumnc[:,mn] = self.h2f(self.bsubumnc[:,mn])
			self.bsubvmnc[:,mn] = self.h2f(self.bsubvmnc[:,mn])
		if self.iasym==1:
			for mn in range(self.mnmax):
				self.lmnc[:,mn] = self.h2f(self.lmnc[:,mn])
			for mn in range(self.mnmax_nyq):
				self.bmns[:,mn] = self.h2f(self.bmns[:,mn])
				self.gmns[:,mn] = self.h2f(self.gmns[:,mn])
				self.bsupumns[:,mn] = self.h2f(self.bsupumns[:,mn])
				self.bsupvmns[:,mn] = self.h2f(self.bsupvmns[:,mn])
				self.bsubsmnc[:,mn] = self.h2f(self.bsubsmnc[:,mn])
				self.bsubumns[:,mn] = self.h2f(self.bsubumns[:,mn])
				self.bsubvmns[:,mn] = self.h2f(self.bsubvmns[:,mn])
		# Calc Eplasma
		self.eplasma = 1.5*4*np.pi*np.pi*sum( self.vp * self.presf ) / self.ns
		# Get mn00
		self.mn00 = None
		for mn in range(self.mnmax):
			if self.xm[mn]==0 and self.xn[mn]==0:
				self.mn00 = mn


	def h2f(self,var_half):
		"""Half to full grid

		This routine takes a 1D field and interpolates it from the half
		to the full grid. For an ns sized array we assumes that the
		first index [0]=0 and is just a placeholder.

		Parameters
		----------
		var_half : list
			Variable on half grid
		Returns
		----------
		var_full : list
			Variable on full grid
		"""
		temp = var_half.copy()
		temp[0] = 1.5 * temp[1] - 0.5 * temp[2]
		temp[1:-1] = 0.5 * (temp[1:-1] + temp[2:])
		temp[-1] = 2.0 * temp[-1] - 1.0 * temp[-2]
		return temp

	def calc_jll(self, theta, phi ):
		"""Compute parallel current density

		This routine computes the paralle current density 

		Parameters
		----------
		theta : ndarray
			Polidal angle grid [rad]
		phi : ndarray
			Toroidal angle grid [rad]
		Returns
		----------
		jll : ndarray
			Parallel current density [A/m^2]
		"""
		b  = self.cfunct(theta,phi,self.bmnc,self.xm_nyq,self.xn_nyq)
		g  = self.cfunct(theta,phi,self.gmnc,self.xm_nyq,self.xn_nyq)
		bu = self.cfunct(theta,phi,self.bsubumnc,self.xm_nyq,self.xn_nyq)
		bv = self.cfunct(theta,phi,self.bsubvmnc,self.xm_nyq,self.xn_nyq)
		ju = self.cfunct(theta,phi,self.currumnc,self.xm_nyq,self.xn_nyq)
		jv = self.cfunct(theta,phi,self.currvmnc,self.xm_nyq,self.xn_nyq)

		if (self.iasym==1):
			b  = b + self.sfunct(theta,phi,self.bmns,self.xm_nyq,self.xn_nyq)
			g  = g + self.sfunct(theta,phi,self.gmns,self.xm_nyq,self.xn_nyq)
			bu = bu + self.sfunct(theta,phi,self.bsubumns,self.xm_nyq,self.xn_nyq)
			bv = bv + self.sfunct(theta,phi,self.bsubvmns,self.xm_nyq,self.xn_nyq)
			ju = ju + self.sfunct(theta,phi,self.currumns,self.xm_nyq,self.xn_nyq)
			jv = jv + self.sfunct(theta,phi,self.currvmns,self.xm_nyq,self.xn_nyq)
		jll = (bu*ju+bv*jv)/(g*b)
		return jll

	def calc_grad_rhosq(self):
		"""Compute <|grad(rho)|^2> 
		This routine flux surface average of |grad(rho)|^2 

		Parameters
		----------
		theta : ndarray
			Polidal angle grid [rad]
		phi : ndarray
			Toroidal angle grid [rad]
		Returns
		----------
		avgrho2 : ndarray
			Flux surface average of |grad(rho)|^2
		"""
		import numpy as np
		nu = 64
		nv = 128
		theta = np.linspace(0,2*np.pi,nu)
		phi   = np.linspace(0,2*np.pi,nv)
		# Create derivatives
		xm2d  = np.broadcast_to(self.xm,(self.ns,self.mnmax))
		xn2d  = np.broadcast_to(self.xn,(self.ns,self.mnmax))
		rumns = - xm2d * self.rmnc
		rvmns = - xn2d * self.rmnc
		zumnc =   xm2d * self.zmns
		zvmnc =   xn2d * self.zmns
		if self.iasym==1:
			rumnc =   xm2d * self.rmns
			rvmnc =   xn2d * self.rmns
			zumns = - xm2d * self.zmnc
			zvmns = - xn2d * self.zmnc
		r = self.cfunct(theta,zeta,self.rmnc,self.xm,self.xn)
		g = self.cfunct(theta,zeta,self.gmnc,self.xm_nyq,self.xn_nyq)
		ru = self.sfunct(theta,zeta,rumns,self.xm,self.xn)
		rv = self.sfunct(theta,zeta,rvmns,self.xm,self.xn)
		zu = self.cfunct(theta,zeta,zumnc,self.xm,self.xn)
		zv = self.cfunct(theta,zeta,zvmnc,self.xm,self.xn)
		if self.iasym==1:
			r  = r  + self.sfunct(theta,zeta,self.rmns,self.xm,self.xn)
			g  = g  + self.sfunct(theta,zeta,self.gmns,self.xm_nyq,self.xn_nyq)
			ru = ru + self.cfunct(theta,zeta,rumnc,self.xm,self.xn)
			rv = rv + self.cfunct(theta,zeta,rvmnc,self.xm,self.xn)
			zu = zu + self.sfunct(theta,zeta,zumns,self.xm,self.xn)
			zv = zv + self.sfunct(theta,zeta,zvmns,self.xm,self.xn)
		# Calc metrics
		gsr = - zu * r
		gsp = zu * rv - ru * zv
		gsz = ru * r
		gs  = ( gsr * gsr + gsp * gsp + gsz * gsz ) / ( g * g )
		for i in range(self.ns):
			gs[i,:,:] = 0.25 * gs[i,:,:] / ( rho[i] * rho[i] )
		vp = np.sum(g,axis=(1,2))
		avgrho2 = np.sum(gs * g,axis=(1,2)) / vp
		avgrho2[1] = 2.0 * avgrho2[1] - avgrho2[2]
		return avgrho2

	def calc_susceptance(self):
		"""Compute susceptance matrix elements 
		This routine calculates the susceptance matrix elements
		S11, S12, S21, S22.

		Returns
		----------
		S11 : ndarray
			Susceptance matrix elements S11
		S12 : ndarray
			Susceptance matrix elements S11
		S21 : ndarray
			Susceptance matrix elements S11
		S22 : ndarray
			Susceptance matrix elements S11
		"""
		import numpy as np
		nu = 64
		nv = 128
		theta = np.linspace(0,2*np.pi,nu)
		phi   = np.linspace(0,2*np.pi,nv)
		# Create derivatives
		xm2d  = np.broadcast_to(self.xm,(self.ns,self.mnmax))
		xn2d  = np.broadcast_to(self.xn,(self.ns,self.mnmax))
		rumns = - xm2d * self.rmnc
		rvmns = - xn2d * self.rmnc
		zumnc =   xm2d * self.zmns
		zvmnc =   xn2d * self.zmns
		lumnc =   xm2d * self.lmns
		lvmnc =   xn2d * self.lmns
		if self.iasym==1:
			rumnc =   xm2d * self.rmns
			rvmnc =   xn2d * self.rmns
			zumns = - xm2d * self.zmnc
			zvmns = - xn2d * self.zmnc
			lumns = - xm2d * self.lmnc
			lvmns = - xn2d * self.lmnc
		r = self.cfunct(theta,zeta,self.rmnc,self.xm,self.xn)
		g = self.cfunct(theta,zeta,self.gmnc,self.xm_nyq,self.xn_nyq)
		ru = self.sfunct(theta,zeta,rumns,self.xm,self.xn)
		rv = self.sfunct(theta,zeta,rvmns,self.xm,self.xn)
		zu = self.cfunct(theta,zeta,zumnc,self.xm,self.xn)
		zv = self.cfunct(theta,zeta,zvmnc,self.xm,self.xn)
		lu = self.cfunct(theta,zeta,lumnc,self.xm,self.xn)
		lv = self.cfunct(theta,zeta,lvmnc,self.xm,self.xn)
		if self.iasym==1:
			r  = r  + self.sfunct(theta,zeta,self.rmns,self.xm,self.xn)
			g  = g  + self.sfunct(theta,zeta,self.gmns,self.xm_nyq,self.xn_nyq)
			ru = ru + self.cfunct(theta,zeta,rumnc,self.xm,self.xn)
			rv = rv + self.cfunct(theta,zeta,rvmnc,self.xm,self.xn)
			zu = zu + self.sfunct(theta,zeta,zumns,self.xm,self.xn)
			zv = zv + self.sfunct(theta,zeta,zvmns,self.xm,self.xn)
			lu = lu + self.sfunct(theta,zeta,lumns,self.xm,self.xn)
			lv = lv + self.sfunct(theta,zeta,lvmns,self.xm,self.xn)
		# Calc suscpetance matrices
		scale_fact = 1.0 / ( 4 * np.pi * np.pi )
		S11 = ( ru * ru + zu * zu)
		S21 = ( ru * rv + zu * zv)
		S12 = ( S12 * ( 1.0 + lu ) - S11 * lv )
		S22 = ( ( rv * rv + zv * zv + r * r ) * ( 1.0 + lu ) - S21 * lv )
		S11 = np.trapz(S11 / g, x=zeta, axis=2)
		S12 = np.trapz(S12 / g, x=zeta, axis=2)
		S21 = np.trapz(S21 / g, x=zeta, axis=2)
		S22 = np.trapz(S22 / g, x=zeta, axis=2)
		S11 = np.trapz(S11, x=theta, axis=1)*scale_fact
		S12 = np.trapz(S12, x=theta, axis=1)*scale_fact
		S21 = np.trapz(S21, x=theta, axis=1)*scale_fact
		S22 = np.trapz(S22, x=theta, axis=1)*scale_fact
		return S11,S12,S21,S22

	def getSpline(self,*args,**kwargs):
		"""Returns a profile in the AUX_S/F form
		This routine returns the pressure, current or rotational
		transform profile in the form AUX form used by the VMEC input
		spline routines.

		Parameters
		----------
		ns : int (optional)
			Number of gridpoints in toroidal flux to use (default ns<=100)
		s : list (optional)
			List of points in normalized toroidal flux to use
		profile : string (optional)
			Which field to use (presf, jcurv, iotaf; default: iotaf)
		lprint : boolean (optional)
			Print the arrays to the screen (default: False)
		Returns
		----------
		aux_s : list
			Array of normalized toroidal flux knots
		aux_f : list
			Array of pressure values [Pa]
		"""
		import numpy as np
		ns_out  = kwargs.get('ns',min(self.ns,100))
		aux_s  = kwargs.get('s',None)
		f_type  = kwargs.get('profile','iotaf')
		lprint  = kwargs.get('lprint',False)
		if type(aux_s) is type(None):
			aux_s = np.linspace(0,1,ns_out)
		if len(aux_s) > 100:
			print(' ERROR: number of toroidal points must be <=100')
			return [],[]
		s = np.linspace(0,1,self.ns)
		f = getattr(self,f_type).flatten()
		aux_f = np.interp(aux_s,s,f)
		if lprint:
			print(rf'  AUX_S = {aux_s}')
			print(rf'  AUX_F = {aux_f}')
		return aux_s,aux_f

	def getCurrentPoloidal(self):
		"""Returns the poloidal total current
		This routine returns the total poloidal current as used by the
		BNORM code.

		Returns
		----------
		curpol : float
			Total poloidal current B_v*2*pi/nfp (m=0,n=0)
		"""
		import numpy as np
		curpol = 1.0
		for mn in range(self.mnmax_nyq):
			if (self.xm_nyq[mn]==0 and self.xn_nyq[mn]==0):
				curpol = 2.0*self.bsubvmnc[self.ns-1,mn]*np.pi/self.nfp 
		return curpol

	def getBcyl(self,R,phi,Z):
		"""Wrapper to the GetBcyl_WOUT function

		This routine wrappers the GetBcyl_WOUT function found in
		vmec_utils.  It takes R, phi, and Z as inputs and returns
		the Br, Bphi, Bz, s, and u values at that point.A status flag
		is also returned (info) which indicates
			 0: successfully find s,u point
			-1: did not converge
			-3: sflux > 1, probably

		Parameters
		----------
		R : real
			Cylindical R coordinate [m].
		phi : real
			Cylindical phi coordinate [rad].
		Z : real
			Cylindical Z coordinate [m].
		Returns
		-------
		br : real
			Magnetic field in cylindrical R direction [T].
		bphi : real
			Magnetic field in cylindrical phi direction [T].
		bz : real
			Magnetic field in cylindrical Z direction [T].
		s : real
			Normalized toroidal flux coordinate [arb].
		u : real
			Poloidal angle coordinate (VMEC angle) [rad].
		info: int
			Status of inverse lookup.
		"""
		return self.libStell.vmec_getBcyl_wout(R,phi,Z)

	def get_flxcoord(self,s,u,v):
		"""Wrapper to the get_flxcoord function

		This routine wrappers the get_flxcoord function found in
		vmec_utils.  It takes s, u, and v as inputs and returns
		the R, phi, Z, dRds, dZds, dRdu, and dZdu values at that
		point.

		Parameters
		----------
		s : real
			Normalized toroidal flux (VMEC).
		u : real
			Poloidal angle [rad] (VMEC)
		v : real
			Toroidal angle [rad] [VMEC]

		Returns
		-------
		R : real
			Cylindical R coordinate [m].
		phi : real
			Cylindical phi coordinate [rad].
		Z : real
			Cylindical Z coordinate [m].
		dRds : real
			Derivative of R coordiante with respect to s (dR/ds)
		dZds : real
			Derivative of Z coordiante with respect to s (dZ/ds)
		dRdu : real
			Derivative of R coordiante with respect to u (dR/du)
		dZdu : real
			Derivative of Z coordiante with respect to u (dZ/du)
		"""
		return self.libStell.vmec_get_flxcoord(s,u,v)

	def extrapSurface(self,surf=None,dist=0.1):
		"""Returns an extrapolated surface.
		This routine extrapolates a surface a given distance using the
		VMEC scaling for modes. The surface to extrapolate (default ns)
		and distance to extrapolate (default 0.1 [m]) are optional
		inputs.

		Parameters
		----------
		surf : int (optional)
			Surface to extrapolate (default: ns)
		dist : float (optional)
			Distance to extrapolate [m] (default 0.1)

		Returns
		----------
		rmnc : ndarray
			R cosine harmonics of extrapolated surface
		zmns : ndarray
			Z sine harmonics of extrapolated surface
		zmns : ndarray
			R sine harmonics of extrapolated surface
		rmnc : ndarray
			Z cosine harmonics of extrapolated surface
		"""
		import numpy as np
		# Handle the surface to use
		if surf:
			k=surf-1
		else:
			k=self.ns-1
		rho = np.sqrt(float(k)/(self.ns-1))
		# Extract surface and remove axis component
		rmnc = np.zeros((self.mnmax))
		zmns = np.zeros((self.mnmax))
		rmns = np.zeros((self.mnmax))
		zmnc = np.zeros((self.mnmax))
		for mn in range(self.mnmax):
			if self.xm[mn]==0:
				rmnc[mn] = self.rmnc[k,mn] - self.rmnc[0,mn]
				zmns[mn] = self.zmns[k,mn] - self.zmns[0,mn]
			else:
				rmnc[mn] = self.rmnc[k,mn]
				zmns[mn] = self.zmns[k,mn]
		if self.iasym == 1:
			zmnc = self.zmnc[k,:]
			rmns = self.rmns[k,:]
			for mn in range(self.mnmax):
				if self.xm[mn]==0:
					zmnc[mn] = zmnc[mn] - self.zmnc[0,mn]
					rmns[mn] = rmns[mn] - self.rmns[0,mn]
		# Scale by a factor
		scale = (self.aminor+dist)/self.aminor
		scale = scale*scale
		scalemn = np.squeeze(np.where(self.xm%2==1, rho*scale, scale))
		rmnc = rmnc * scalemn
		zmns = zmns * scalemn
		if self.iasym == 1:
			zmnc = zmnc * scalemn
			rmns = rmns * scalemn
		# Add axis back in
		for mn in range(self.mnmax):
			if self.xm[mn]==0:
				rmnc[mn] = rmnc[mn] + self.rmnc[0,mn]
				zmns[mn] = zmns[mn] + self.zmns[0,mn]
		if self.iasym == 1:
			for mn in range(self.mnmax):
				if self.xm[mn]==0:
					zmnc[mn] = zmnc[mn] + self.zmnc[0,mn]
					rmns[mn] = rmns[mn] + self.rmns[0,mn]
		return rmnc, zmns, rmns, zmnc

	def offsetCurve(self, R, Z, distance=0.1):
		"""Returns points an offset distance from a curve
		The routine takes a set of points in R and Z and calculates a
		set of offset points.

		Parameters
		----------
		R : ndarray
			Array of points in R [m]
		Z : ndarray
			Array of points in Z [m]
		dist : float (optional)
			Distance to extrapolate [m] (default 0.1)

		Returns
		----------
		R : ndarray
			Offset points in R [m]
		Z : ndarray
			Offset points in Z [m]
		"""
		import numpy as np

		# Calculate the tangents
		dx = np.diff(R, append=R[0])
		dy = np.diff(Z, append=Z[0])

		# Calculate the normals
		lengths = np.sqrt(dx**2 + dy**2)
		normals_x = -dy / lengths
		normals_y = dx / lengths

		# Offset the points
		offset_x = R + distance * normals_x
		offset_y = Z + distance * normals_y
		return offset_x, offset_y

	def fitSurface(self,surf=None,dist=0.1):
		"""Returns Fourier Harmonics of a fitted surface
		This routine calculates a fit of Fourier coefficients to
		a supplied surface or a surface a uniform distance from the
		VMEC equilibrium boundary. If a surface is provided it must be
		provided as a set of ordered points x,y,z [3,nu,nv].

		Parameters
		----------
		surf : ndarray [3,nu,nv] (optional)
			Surface to fit (default: use dist)
		dist : float (optional)
			Distance to extrapolate [m] (default 0.1)

		Returns
		----------
		rmnc : ndarray
			R cosine harmonics of extrapolated surface
		zmns : ndarray
			Z sine harmonics of extrapolated surface
		zmns : ndarray
			R sine harmonics of extrapolated surface
		rmnc : ndarray
			Z cosine harmonics of extrapolated surface
		"""
		import numpy as np
		from scipy.optimize import minimize
		ns1 = self.ns-1
		if surf:
			theta  = surf[0]
			phi    = surf[1]
			boundr = surf[3]
			boundz = surf[4]
			nu     = len(theta)
			nv     = len(phi)
		else:
			nu = 64
			nv = 16
			theta = np.ndarray((nu,1))
			phi   = np.ndarray((nv,1))
			for j in range(nu): theta[j]=2.0*np.pi*j/float(nu)
			if self.iasym == 1:
				for j in range(nv):  phi[j] = 2.0*np.pi*j/float(nv-1)
			else:
				for j in range(nv):  phi[j] =    np.pi*j/float(nv-1)
			phi = phi / self.nfp
			r = self.cfunct(theta,phi,self.rmnc,self.xm,self.xn)
			z = self.sfunct(theta,phi,self.zmns,self.xm,self.xn)
			if self.iasym==1:
				r = r + self.sfunct(theta,phi,self.rmns,self.xm,self.xn)
				z = z + self.cfunct(theta,phi,self.zmnc,self.xm,self.xn)
			boundr = np.squeeze(r[ns1,:,:])
			boundz = np.squeeze(z[ns1,:,:])
			for v in range(nv):
				x = boundr[:,v]
				y = boundz[:,v]
				x,y = self.offsetCurve(x,y,dist)
				boundr[:,v] = x
				boundz[:,v] = y
		opts = (nu,nv,self.xm,self.xn,theta,phi,boundr,boundz)
		self.fitFactor = (1.0+np.squeeze(self.xm))
		r_temp = np.squeeze(self.rmnc[ns1,:])*self.fitFactor
		z_temp = np.delete(self.zmns[ns1,:]*self.fitFactor,self.mn00)
		x0 = np.concatenate((r_temp,z_temp))
		if self.iasym == 1:
			r1_temp = np.delete(self.rmns[ns1,:]*self.fitFactor,self.mn00)
			z1_temp = np.squeeze(self.zmnc[ns1,:],self.mn00)*self.fitFactor
			x0 = np.concatenate((r_temp,z_temp,r1_temp,z1_temp))
		self.Nfeval = 0
		self.Rfit = boundr
		self.Zfit = boundz
		res = minimize(self.fit_func, x0, method='CG',
               args=opts, options={'disp': True},
               callback=self.callbackF,tol=1.0E-1)
		xf = res.x
		i1 = 0; i2 = i1+self.mnmax
		rmnc = np.broadcast_to(xf[i1:i2],(1,self.mnmax))/self.fitFactor
		i1 = i2; i2 = i1+self.mnmax-1
		xtemp = xf[i1:i2]
		zmns = np.broadcast_to(np.insert(xtemp,self.mn00,0.0),(1,self.mnmax))/self.fitFactor
		if self.iasym==1:
			i1 = i2; i2 = i1+self.mnmax-1
			xtemp = xf[i1:i2]
			rmns = np.broadcast_to(np.insert(xtemp,self.mn00,0.0),(1,self.mnmax))/self.fitFactor
			i1 = i2; i2 = i1+self.mnmax
			zmnc = np.broadcast_to(xf[i1:i2],(1,self.mnmax))/self.fitFactor
		else:
			rmns = np.zeros((self.mnmax))
			zmnc = np.zeros((self.mnmax))
		return np.squeeze(rmnc),np.squeeze(zmns),np.squeeze(rmns),np.squeeze(zmnc)

	def fit_func(self,x,*args):
		import numpy as np
		nu = args[0]
		nv = args[1]
		xm = args[2]
		xn = args[3]
		theta = args[4]
		phi = args[5]
		boundr = args[6]
		boundz = args[7]
		mnmax = xm.shape[0]
		i1 = 0; i2 = i1+mnmax
		rmnc = np.broadcast_to(x[i1:i2],(1,mnmax))/self.fitFactor
		i1 = i2; i2 = i1+mnmax-1
		xtemp = x[i1:i2]
		zmns = np.broadcast_to(np.insert(xtemp,self.mn00,0.0),(1,mnmax))/self.fitFactor
		r = self.cfunct(theta,phi,rmnc,xm,xn)
		z = self.sfunct(theta,phi,zmns,xm,xn)
		if x.shape[0] > 2*mnmax:
			i1 = i2; i2 = i1+mnmax-1
			xtemp = x[i1:i2]
			rmns = np.broadcast_to(np.insert(xtemp,self.mn00,0.0),(1,mnmax))/self.fitFactor
			i1 = i2; i2 = i1+mnmax
			zmnc = np.broadcast_to(x[i1:i2],(1,mnmax))/self.fitFactor
			r = r + self.sfunct(theta,phi,rmns,xm,xn)
			z = z + self.cfunct(theta,phi,zmnc,xm,xn)
		# Calc distance
		d = 0
		for u in range(nu):
			for v in range(nv):
				dr = boundr[u,v] - np.squeeze(r[0,:,v])
				dz = boundz[u,v] - np.squeeze(z[0,:,v])
				dl2 = dr*dr + dz*dz
				d = d + min(dl2)
		return d
		
	def callbackF(self,intermediate_result):
		print(f'ITER: {self.Nfeval} -- dval: {intermediate_result.fun}')
		self.Nfeval += 1

	def callbackF_plot(self,Xi):
		import numpy as np
		import matplotlib.pyplot as pyplot
		mnmax = self.mnmax
		nu = 32
		nv = 3
		theta = np.ndarray((nu,1))
		phi   = np.ndarray((nv,1))
		for j in range(nu): theta[j]=2.0*np.pi*j/float(nu-1)
		for j in range(nv):  phi[j] =    np.pi*j/float(nv-1)
		i1 = 0; i2 = i1+mnmax
		rmnc = np.broadcast_to(Xi[i1:i2],(1,mnmax))/self.fitFactor
		i1 = i2; i2 = i1+mnmax
		xtemp = Xi[i1:i2]
		zmns = np.broadcast_to(np.insert(xtemp,self.mn00,0.0),(1,mnmax))/self.fitFactor
		r = self.cfunct(theta,phi,rmnc,self.xm,self.xn/self.nfp)
		z = self.sfunct(theta,phi,zmns,self.xm,self.xn/self.nfp)
		px = 1/pyplot.rcParams['figure.dpi']
		fig,ax = pyplot.subplots(1,1,figsize=(1024*px,768*px))
		ax.plot(self.Rfit[:,0],self.Zfit[:,0],'ko')
		ax.plot(self.Rfit[:,-1],self.Zfit[:,-1],'ko')
		ax.plot(r[0,:,0],z[0,:,0],'r')
		ax.plot(r[0,:,1],z[0,:,1],'g')
		ax.plot(r[0,:,2],z[0,:,2],'b')
		pyplot.show(block=False)
		pyplot.pause(0.1)
		print('{0:4d}   {1: 3.6f}   {2: 3.6f}   {3: 3.6f}'.format(self.Nfeval, Xi[0], Xi[1], Xi[2]))
		self.Nfeval += 1
		pyplot.close(fig)

# VMEC INDATA Class
class VMEC_INDATA():
	"""Class for working with VMEC equilibria

	"""
	def __init__(self, parent=None):
		self.nfp = None
		self.libStell = LIBSTELL()

	def read_indata(self,filename):
		"""Reads INDATA namelist from a file

		This routine wrappers the vmec_input module reading routine.
		Parameters
		----------
		filename : string
			Input file name with INDATA namelist
		"""
		indata_dict = self.libStell.read_indata(filename)
		for key in indata_dict:
			setattr(self, key, indata_dict[key])
		# generate helpers
		#print(self.rbc.shape)

	def write_indata(self,filename):
		"""Writes INDATA namelist to a file

		This routine wrappers the vmec_input module writing routine.
		Parameters
		----------
		filename : string
			Input file name to write INDATA namelist to
		"""
		#print(self.__dict__)
		out_dict = vars(self)
		#del out_dict['libStell']
		#print(d(self))
		self.libStell.write_indata(filename,out_dict)

	def pmass(self,x):
		"""Wrapper to the PMASS function

		This routine wrappers the PMASS function which
		returns the mass(s) pressure function.

		Parameters
		----------
		s : real
			Value of normalized toroidal flux
		Returns
		-------
		val : real
			Value of mass(s)
		"""
		return self.libStell.pmass(x)

	def piota(self,x):
		"""Wrapper to the PIOTA function

		This routine wrappers the PIOTA function which
		returns the iota(s) function.

		Parameters
		----------
		s : real
			Value of normalized toroidal flux
		Returns
		-------
		val : real
			Value of iota
		"""
		return self.libStell.piota(x)

	def pcurr(self,x):
		"""Wrapper to the PCURR function

		This routine wrappers the PCURR function which
		returns the current profile.

		Parameters
		----------
		s : real
			Value of normalized toroidal flux
		Returns
		-------
		val : real
			Value of current profile
		"""
		return self.libStell.pcurr(x)



# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)