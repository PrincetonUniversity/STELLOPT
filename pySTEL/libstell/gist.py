##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for reading and handling 
GIST data.
"""

# Libraries
from libstell.libstell import LIBSTELL

# Constants

# GIST Class
class GIST():
	"""Class for working with GIST data

	"""
	def __init__(self):
		self.g11 = None

	def read_gist(self,filename):
		"""Reads a GIST geometry file

		This routine reads the GIST geometry files.

		Parameters
		----------
		file : str
			Path to GIST file.
		"""
		import numpy as np
		f = open(filename,'r')
		lines = f.readlines()
		f.close()
		i = 0
		temp_dict={}
		while (i < len(lines)):
			if '!s0, alpha0' in lines[i]:
				[txt1,txt2,txteq,s0_txt,alpha0_txt] = lines[i].split()
				self.s0 = float(s0_txt)
				self.alpha0 = float(alpha0_txt)
			elif '!major, minor radius[m]=' in lines[i]:
				[txt1,txt2,txt3,Rmajor_txt,Aminor_txt] = lines[i].split()
				self.Rmajor = float(Rmajor_txt)
				self.Aminor = float(Aminor_txt)
			elif 'my_dpdx =' in lines[i]:
				[txt1,dpdx_txt] = lines[i].split('=')
				self.dpdx = float(dpdx_txt)
			elif 'q0 =' in lines[i]:
				[txt1,q0_txt] = lines[i].split('=')
				self.q0 = float(q0_txt)
			elif 'shat =' in lines[i]:
				[txt1,shat_txt] = lines[i].split('=')
				self.shat = float(shat_txt)
			elif 'gridpoints =' in lines[i]:
				[txt1,gridpoints_txt] = lines[i].split('=')
				self.gridpoints = int(gridpoints_txt)
			elif 'n_pol =' in lines[i]:
				[txt1,n_pol_txt] = lines[i].split('=')
				self.n_pol = int(n_pol_txt)
			elif '/' in lines[i]:
				i0 = i + 1
				self.g11     = np.zeros(self.gridpoints)
				self.g12     = np.zeros(self.gridpoints)
				self.g22     = np.zeros(self.gridpoints)
				self.Bhat    = np.zeros(self.gridpoints)
				self.abs_jac = np.zeros(self.gridpoints)
				self.L2      = np.zeros(self.gridpoints)
				self.L1      = np.zeros(self.gridpoints)
				self.dBdt    = np.zeros(self.gridpoints)
				for j in range(self.gridpoints):
					k = i0 + j
					txt        = lines[k].split()
					self.g11[j]     = float(txt[0])
					self.g12[j]     = float(txt[1])
					self.g22[j]     = float(txt[2])
					self.Bhat[j]   = float(txt[3])
					self.abs_jac[j] = float(txt[4])
					self.L2[j]      = float(txt[5])
					self.L1[j]      = float(txt[6])
					self.dBdt[j]    = float(txt[7])
				self.kp1 = self.L2 - self.dpdx/2.0/self.Bhat
			i = i + 1

	def calcGist(self,alpha,maxpnt,vmec_data,s):
		"""Compute the GIST quantities from inputs

		This routine computes the various GIST quantities based on a 
		set of user inputs. It mimics the behavior of the STELLOPT
		stellopt_txport code.

		"""
		import numpy as np
		self.s0 = s
		self.alpha0 = alpha
		self.Rmajor = vmec_data.Rmajor
		self.Aminor = vmec_data.Aminor
		Ba = ABS(vmec_data.phi[-1]/(np.pi*vmec_data.Aminor*vmec_data.Aminor))
		iota = vmec_data.getiota(s)
		iotap = vmec_data.getiotaprime(s)
		pressp = vmec_data.getpressureprime(s)
		mu0 = 4E-7*np.pi
		q = 1.0/iota
		qprime = -iota*q*q
		shat   = 2*s/q*qprime
		dpdx   = -4.0*np.sqrt(s)/Ba**2 * pressp*mu0
		theta = np.linspace(-np.pi,np.pi,maxpnt)
		self.dpdx = dpdx
		self.q0   = q
		self.shat = shat
		self.gridpoints = maxpnt
		self.n_pol = 1
		self.g11     = np.zeros(self.gridpoints)
		self.g12     = np.zeros(self.gridpoints)
		self.g22     = np.zeros(self.gridpoints)
		self.Bhat    = np.zeros(self.gridpoints)
		self.abs_jac = np.zeros(self.gridpoints)
		self.L2      = np.zeros(self.gridpoints)
		self.L1      = np.zeros(self.gridpoints)
		self.dBdt    = np.zeros(self.gridpoints)
		for u in range(maxpnt):
			zeta = alpha + q*(theta[u]-0.0)
			thetastar = vmec_data.getTheta(s,theta[u],zeta)
			R,phi,Z,dRds,dZds,dRdu,dZdu,dRdv,dZdv = vmec_data.get_flxcoord(s,thetastar,zeta)
			# Calc covariant vectors
			esubs = [dRds,0.0,dZds]
			esubu = [dRdu,0.0,dZdu]
			esubv = [dRdv*vmec_data.nfp,R,dZdv*vmec_data.nfp]
			# Calc Jacobian
			sqrtg = R*(dRdu*dZds-dRds*dZdu)
			# Calc contravaiant vectors
			es = np.array([esubu[1]*esubv[2]-esubu[2]*esubv[1],
						   esubu[2]*esubv[0]-esubu[0]*esubv[2],
			    		   esubu[0]*esubv[1]-esubu[1]*esubv[0]])/sqrtg
			eu = np.array([esubv[1]*esubs[2]-esubv[2]*esubs[1],
						   esubv[2]*esubs[0]-esubv[0]*esubs[2],
						   esubv[0]*esubs[1]-esubv[1]*esubs[0]])/sqrtg
			ev = np.array([esubs[1]*esubu[2]-esubs[2]*esubu[1],
						   esubs[2]*esubu[0]-esubs[0]*esubu[2],
						   esubs[0]*esubu[1]-esubs[1]*esubu[0]])/sqrtg
			# Calc Grad(B)
			b = vmec_data.cfunct([thetastar],[zeta],vmec_data.bmnc,vmec_data.xm_nyq,vmec_data.xn_nyq)
			bumns = -vmec_data.bmnc*np.tile(vmec_data.xm_nyq,(1,vmec_data.ns)).T
			bvmns =  vmec_data.bmnc*np.tile(vmec_data.xn_nyq,(1,vmec_data.ns)).T
			m = np.tile(vmec_data.xm_nyq,[1,vmec_data.ns])
			x = np.linspace(0,1,vmec_data.ns)
			f = np.diff(b,prepend=0)*(self.ns-1)
			modb = np.interp(s,x,b)
			bs = np.interp(s,x,f)
			bu = vmec_data.sfunct([thetastar],[zeta],bumns,vmec_data.xm_nyq,vmec_data.xn_nyq)
			bv = vmec_data.sfunct([thetastar],[zeta],bvmns,vmec_data.xm_nyq,vmec_data.xn_nyq)
			gradb = bs*es + bu*eu + bv*ev
			# Adjust eu to include lambda factor
			lam = vmec_data.cfunct([thetastar],[zeta],vmec_data.bmnc,vmec_data.xm_nyq,vmec_data.xn_nyq)
			lumns = -vmec_data.lmnc*np.tile(vmec_data.xm,(1,vmec_data.ns)).T
			lvmns =  vmec_data.lmnc*np.tile(vmec_data.xn,(1,vmec_data.ns)).T
			x = np.linspace(0,1,vmec_data.ns)
			f = np.diff(lam,prepend=0)*(self.ns-1)
			ls = np.interp(s,x,f)
			lu = vmec_data.sfunct([thetastar],[zeta],lumns,vmec_data.xm,vmec_data.xn)
			lv = vmec_data.sfunct([thetastar],[zeta],lvmns,vmec_data.xm,vmec_data.xn)
			eu = eu + ls*es + lu*eu + lv*ev
			# Calc metric elments
			gradA = thetastar*qprime*es + q*eu - ev
			jac1 = 1.0/np.array([es[1]*gradA[2]-es[2]*gradA[1],
						   es[2]*gradA[0]-es[0]*gradA[2],
						   es[0]*gradA[1]-es[1]*gradA[0]])
			gss = np.sum(es*es)
			gsa = np.sum(es*gradA)
			gst = np.sum(es*eu)
			gaa = np.sum(gradA*gradA)
			gat = np.sum(gradA*eu)
			alpha = q*thetaStar - v
			# Now output the values
			self.Bhat[u] = modb/Ba
			self.g11[u]  = gss*self.Aminor*self.Aminor*0.25/s
			self.g12[u]  = gsa*self.Aminor*self.Aminor*iota*0.5
			self.g22[u] = (self.Bhat[u]*self.Bhat[u] + self.g12[u]*self.g12[u])/self.g11[u]
			self.abs_jac[u] = ABS(jac1*2*q/self.Aminor**3)
			# Reuse some variables (order matters)
			ea = np.array([eu[1]*es[2]-eu[2]*es[1],
						   eu[2]*es[0]-eu[0]*es[2],
						   eu[0]*es[1]-eu[1]*es[0]])*jac1
			et = np.array([es[1]*gradA[2]-es[2]*gradA[1],
						   es[2]*gradA[0]-es[0]*gradA[2],
						   es[0]*gradA[1]-es[1]*gradA[0]])*jac1
			es = np.array([gradA[1]*eu[2]-gradA[2]*eu[1],
						   gradA[2]*eu[0]-gradA[0]*eu[2],
						   gradA[0]*eu[1]-gradA[1]*eu[0]])*jac1
			gradB = gradB/Ba
			dBds  = np.sum(gradB,es)
			dBda  = np.sum(gradB,ea)
			self.dBdt[u] = np.sum(gradB,et)
			c     = iota*iota*self.Aminor**4
			self.L1[u] = q/np.sqrt(s)*(dBda + c*(gss*gat-gsa*gst)*self.dBdt[u]/(4*self.Bhat[u]**2))
			self.L2[u] = 2*np.sqrt(s)*(dBds + c*(gaa*gst-gsa*gat)*self.dBdt[u]/(4*self.Bhat[u]**2))
		self.kp1 = self.L2 - self.dpdx/2.0/self.Bhat

	def calcProxG11(self):
		"""Computes the G11 proxy

		This routine computes the G11 proxy

		Returns
		-------
		g11 : float
			The G11 proxy
		"""
		return self.g11

	def calcProx1(self):
		"""Computes the Prox1 proxy

		This routine computes the Prox1 proxy

		Returns
		-------
		Prox1 : float
			The Prox1 proxy
		"""
		import numpy as np
		from scipy.signal import savgol_filter
		dk     = 0.9596
		rhsrkp = 0.002
		tau_s  = 1.1267
		rkp_p  = 3.0
		rkp_cr = 0.05335
		tau    = 1.0
		vrbl     = self.g12/self.g11
		sloc     = np.zeros(self.gridpoints)
		sloc[1:] = (vrbl[1:] - vrbl[0:-1])*(self.gridpoints-1)/(2.0*np.pi)
		sloc[0] = sloc[1]
		slocav = savgol_filter(sloc, 8, 3)
		rkp1av = savgol_filter(self.kp1, 8, 3)
		qqfac = (1.0 + 1.0 / (rhsrkp * (1.0 + (tau_s*slocav)**2)))
		dkp = np.zeros(self.gridpoints)
		dkp[:] = rkp_p-rkp_cr
		dkpfac = np.zeros(self.gridpoints)
		dkpfac[dkp > 0.0] = dkp[dkp>0.0]
		vkp1fac = np.zeros(self.gridpoints)
		vkp1fac[rkp1av < 0.0] = -rkp1av[rkp1av<0.0]
		return (dk/rkp_p)*np.sqrt(tau*vkp1fac*dkpfac)*qqfac

	def calcProx1b(self):
		"""Computes the Prox1b proxy

		This routine computes the Prox1b proxy

		Returns
		-------
		Prox1b : float
			The Prox1b proxy
		"""
		import numpy as np
		from scipy.signal import savgol_filter
		dk     = 0.9596
		rhsrkp = 0.002
		tau_s  = 1.1267
		rkp_p  = 3.0
		rkp_cr = 0.05335
		tau    = 1.0
		rkp1av = savgol_filter(self.kp1, 8, 3)
		qqfac = (1.0 + 1.0 / (rhsrkp * (1.0)))
		dkp = np.zeros(self.gridpoints)
		dkp[:] = rkp_p-rkp_cr
		dkpfac = np.zeros(self.gridpoints)
		dkpfac[dkp > 0.0] = dkp[dkp>0.0]
		vkp1fac = np.zeros(self.gridpoints)
		vkp1fac[rkp1av < 0.0] = -rkp1av[rkp1av<0.0]
		return (dk/rkp_p)*np.sqrt(tau*vkp1fac*dkpfac)*qqfac

	def calcProx1d(self):
		"""Computes the Prox1d proxy

		This routine computes the Prox1d proxy

		Returns
		-------
		Prox1d : float
			The Prox1d proxy
		"""
		import numpy as np
		from scipy.signal import savgol_filter
		dk     = 0.9596
		rhsrkp = 0.002
		tau_s  = 1.1267
		rkp_p  = 3.0
		rkp_cr = 0.05335
		tau    = 1.0
		rkp1av = savgol_filter(self.kp1, 8, 3)
		qqfac = self.g11
		dkp = np.zeros(self.gridpoints)
		dkp[:] = rkp_p-rkp_cr
		dkpfac = np.zeros(self.gridpoints)
		dkpfac[dkp > 0.0] = dkp[dkp>0.0]
		vkp1fac = np.zeros(self.gridpoints)
		vkp1fac[rkp1av < 0.0] = -rkp1av[rkp1av<0.0]
		return (dk/rkp_p)*np.sqrt(tau*vkp1fac*dkpfac)*qqfac

	def calcProxL2(self):
		"""Computes the ProxL2 proxy

		This routine computes the ProxL2 proxy

		Returns
		-------
		ProxL2 : float
			The ProxL2 proxy
		"""
		vqqprox = self.L2
		vqqprox[vqqprox>0] = 0.01*vqqprox[vqqprox>0]
		return vqqprox

	def calcProxTEMOverlap(self):
		"""Computes the TEMOverlap proxy

		This routine computes the TEMOverlap proxy

		Returns
		-------
		Proxy : float
			The TEMOverlap proxy
		"""
		bmax = max(self.Bhat)
		bmin = min(self.Bhat)
		b0   = 0.5*(bmin+bmax)
		bamp = 0.5*(bmax-bmin)
		vqqprox = self.kp1*(self.Bhat-b0)/bamp
		return vqqprox

	def calcProxTEMBounce(self):
		"""Computes the TEMBounce proxy

		This routine computes the TEMBounce proxy

		Returns
		-------
		Proxy : float
			The TEMBounce proxy
		"""
		import numpy as np
		from scipy.interpolate import CubicSpline
		from scipy.integrate import quad
		bmax = max(self.Bhat)
		bmin = min(self.Bhat)
		x = np.linspace(0.0,1.0,self.gridpoints)
		self.Bhat_spl = CubicSpline(x,self.Bhat)
		self.L2_spl = CubicSpline(x,self.L2)
		A = 1.0 / bmax
		B = 1.0 / bmin
		f = lambda x: self.TEMfunc_Bounce(x)
		y,yerr = quad(f,A,B,epsabs=1.0E-9, epsrel=1.0E-3)
		return -y

	def TEMfunc_Bounce(self,lam_TEM):
		"""TEMBounce proxy function (may not need)

		This routine computes the TEMBounce proxy function

		Returns
		-------
		val : float
			The TEMBounce proxy function
		"""
		import numpy as np
		from scipy.integrate import quad
		f = lambda x,a: self.TEMsubfunc_Bounce(x,a)
		A = 0.0
		B = np.pi*2
		y,yerr = quad(f,A,B,args=(lam_TEM,),epsabs=1.0E-9, epsrel=1.0E-3)
		return y

	def TEMsubfunc_Bounce(self,x,lam_TEM):
		"""TEMBounce proxy subfunction

		This routine computes the TEMBounce proxy subfunction

		Returns
		-------
		val : float
			The TEMBounce proxy subfunction
		"""
		import numpy as np
		val = 0.0
		B0 = self.Bhat_spl(x)
		if (0.9999/lam_TEM - B0) >= 0:
			L20 = self.L2_spl(x)
			val = L20 * (1.0 - 0.5 * lam_TEM*B0)/np.sqrt(1.0 - lam_TEM*B0)
		return val

	def calcProxTEMTau(self):
		"""Computes the TEMTau proxy

		This routine computes the TEMTau proxy

		Returns
		-------
		Proxy : float
			The TEMTau proxy
		"""
		import numpy as np
		from scipy.interpolate import CubicSpline
		from scipy.integrate import quad
		bmax = max(self.Bhat)
		bmin = min(self.Bhat)
		x = np.linspace(0.0,1.0,self.gridpoints)
		self.Bhat_spl = CubicSpline(x,self.Bhat)
		self.L2_spl = CubicSpline(x,self.L2)
		A = 1.0 / bmax
		B = 1.0 / bmin
		f = lambda x: self.TEMfunc_Tau(x)
		y,yerr = quad(f,A,B,epsabs=1.0E-9, epsrel=1.0E-3)
		return -y

	def TEMsubfunc_Tau(self,x,lam_TEM):
		"""TEMTau proxy subfunction

		This routine computes the TEMTau proxy subfunction

		Returns
		-------
		val : float
			The TEMTau proxy subfunction
		"""
		import numpy as np
		val = 0.0
		B0 = self.Bhat_spl(x)
		if (0.9999/lam_TEM - B0) >= 0:
			val = 1.0/np.sqrt(1.0 - lam_TEM*B0)
		return val

	def TEMfunc_Tau(self,lam_TEM):
		"""TEMTau proxy function

		This routine computes the TEMTau proxy function

		Returns
		-------
		val : float
			The TEMTau proxy function
		"""
		import numpy as np
		from scipy.integrate import quad
		f = lambda x,a: self.TEMsubfunc_Tau(x,a)
		A = 0.0
		B = np.pi*2
		val_w, val_werr = quad(f,A,B,args=(lam_TEM,))
		if val_w == 0.0:
			return 0.0
		A = 0.0
		B = np.pi*2
		val_t,val_terr = quad(f,A,B,args=(lam_TEM,))
		return val_w/val_t

# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)

