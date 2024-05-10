#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for reading and handling VMEC
equilibrium data.
"""

# Libraries
import sys
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



# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)