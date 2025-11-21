##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for reading and handling BOOZER
transform data.
"""

# Libraries
from libstell.libstell import LIBSTELL, FourierRep

# Constants

# VMEC Class
class BOOZER(FourierRep):
	"""Class for working with VMEC equilibria

	"""
	def __init__(self):
		super().__init__()
		self.nfp = None
		self.libStell = LIBSTELL()

	def read_boozer(self,filename):
		"""Reads a BOOZER boozmn file

		This routine reads and initilizes the BOOZER class
		with variable information from a BOOZER boozmn file.

		Parameters
		----------
		file : str
			Path to wout file.
		"""
		import numpy as np
		import copy
		boozmn_dict = copy.deepcopy(self.libStell.read_boozer(filename))
		for key in boozmn_dict:
			setattr(self, key, boozmn_dict[key])
		self.mboz_b = int(max(np.squeeze(self.ixm_b)))
		nmax = int(max(np.squeeze(self.ixn_b))/self.nfp_b)

	def plotBmnSpectrum(self,sval,ax=None):
		"""Plots the boozer spectrum for a surface

		This routine plots the boozer spectrum for a given
		surface.

		Parameters
		----------
		sval : int
			Surface to plot
		ax : axes (optional)
			Matplotlib axes object to plot to.
		"""
		import numpy as np
		import matplotlib.pyplot as pyplot
		lplotnow = False
		if not ax:
			ax = pyplot.axes()
			lplotnow = True
		# Array extents
		mmax = int(max(np.squeeze(self.ixm_b)))
		nmax = int(max(np.squeeze(self.ixn_b))/self.nfp_b)
		# Sort BMN into array
		bmn = np.zeros((mmax+1,2*nmax+1))
		for mn in range(self.mnboz_b):
			m = self.ixm_b[mn,0]
			n = int(self.ixn_b[mn,0]/self.nfp_b) + nmax
			bmn[m,n] = self.bmnc_b[sval,mn]
		#Plot
		x = np.linspace(0,mmax,mmax+1)
		y = np.linspace(-nmax,nmax,2*nmax+1)
		hmesh=ax.pcolormesh(x,y,np.log10(np.abs(bmn.T)),cmap='jet',shading='gouraud')
		ax.set_xlabel('Poloidal Modes (m)')
		ax.set_ylabel('Toroidal Modes (n)')
		ax.set_title(rf'BOOZER $s$ = {sval/self.ns_b:5.2f}')
		pyplot.colorbar(hmesh,label='$log_{10}$[T]',ax=ax)
		if lplotnow: pyplot.show()

	def plotBsurf(self,sval,ax=None,cmap='jet'):
		"""Plots the boozer |B| on a surface

		This routine plots the boozer |B| for a given
		surface.

		Parameters
		----------
		sval : int
			Surface to plot
		ax : axes (optional)
			Matplotlib axes object to plot to.
		"""
		import numpy as np
		import matplotlib.pyplot as pyplot
		lplotnow = False
		if not ax:
			ax = pyplot.axes()
			lplotnow = True
		theta = np.linspace([0],[np.pi*2.0],360)
		zeta  = np.linspace([0],[np.pi*2.0],256)
		b = self.cfunct(theta,zeta,self.bmnc_b,self.ixm_b,self.ixn_b/self.nfp_b)
		hmesh=ax.pcolormesh(np.squeeze(zeta),np.squeeze(theta),np.squeeze(b[sval,:,:]),cmap=cmap,shading='gouraud')
		ax.plot(zeta,zeta*self.iota_b[sval],'w')
		ax.set_xlabel(r'Toroidal Angle ($\phi$) [rad]')
		ax.set_ylabel(r'Poloidal Angle ($\theta$) [rad]')
		ax.set_title(rf'BOOZER $s$ = {sval/self.ns_b:5.2f}')
		pyplot.colorbar(hmesh,label='[T]',ax=ax)
		if lplotnow: pyplot.show()

	def plot_fieldline(self,nlines=4,*args,**kwargs):
		"""Plots the boozer fieldline in 3D

		This routine plots the boozer fieldline in 3D.

		Parameters
		----------
		sval : int
			Surface to plot
		nlines : int
			Number of lines to plot on surface (default=4)
		plot3D : plot3D object (optional)
			Plotting object to render to.
		"""
		import numpy as np
		import vtk
		from libstell.plot3D import PLOT3D 
		plot3D  = kwargs.get('plot3D',None)
		nphi = kwargs.get('nphi',360)
		ntheta = kwargs.get('ntheta',4)
		phimin = kwargs.get('phimin',0.0)
		phimax = kwargs.get('phimax',2*np.pi)
		sdex   = kwargs.get('sdex',[self.ns_b-1])
		lrender = False
		if not plot3D:
			plt = PLOT3D()
			lrender = True
		theta0 = np.linspace([0],[np.pi*2.0],ntheta+1)
		theta0 = theta0[0:-1]
		phi   = np.linspace([phimin],[phimax],nphi)
		points_array = np.zeros((nphi,3))
		scalar       = np.zeros((nphi,1))
		xout         = np.zeros((self.ns_b,ntheta,nphi))
		yout         = np.zeros((self.ns_b,ntheta,nphi))
		zout         = np.zeros((self.ns_b,ntheta,nphi))
		bout         = np.zeros((self.ns_b,ntheta,nphi))
		ph           = np.zeros((ntheta,1))
		for i in sdex:
			for k,pht in enumerate(phi):
				ph[:,0] = pht
				th = theta0 + self.iota_b[i,0]*pht
				r = self.cfunct(th,ph,self.rmnc_b,self.ixm_b,self.ixn_b)
				z = self.sfunct(th,ph,self.zmns_b,self.ixm_b,self.ixn_b)
				b = self.cfunct(th,ph,self.bmnc_b,self.ixm_b,self.ixn_b)
				p = self.sfunct(th,ph,self.pmns_b,self.ixm_b,self.ixn_b)
				phi_cart = p+ph
				xout[i,:,k] = r[i,:,0]*np.cos(phi_cart[i,:,0])
				yout[i,:,k] = r[i,:,0]*np.sin(phi_cart[i,:,0])
				zout[i,:,k] = z[i,:,0]
				bout[i,:,k] = b[i,:,0]
		for i in sdex:
			for j in range(ntheta):
				points_array[:,0] = xout[i,j,:]
				points_array[:,1] = yout[i,j,:]
				points_array[:,2] = zout[i,j,:]
				scalar            = bout[i,j,:]
				points = vtk.vtkPoints()
				for point in points_array:
					points.InsertNextPoint(point)
				scalar = plt.valuesToScalar(scalar)
				plt.add3Dline(points,linewidth=2,scalars=scalar)
		# In case it isn't set by user.
		plt.setBGcolor()
		if lrender: plt.render()

	def calcQuasiError(self,m,n):
		"""Calculates the quasi-symmetry error for each surface

		This routine computes the quasi-symmetry error for each
		surface in the datastructure for which the boozer
		transformation has been performed. For QAS error (m=1,n=0),
		for QPS error (m=0,n=1), for helical symmetry error both
		m and n should be finite.

		Parameters
		----------
		m : int
			Poloidal spectrum symmetry
		n : int
			Toroidal spectrum symmetry

		Returns
		-------
		error : ndarray
			Quasi-symmetry error
		"""
		from copy import deepcopy
		import numpy as np
		error = np.zeros((self.ns_b))
		if m == 0:
			maskdex = self.ixn_b == 0
		elif n==0:
			maskdex = self.ixm_b == 0
		else:
			maskdex = (self.ixm_b == m) & (self.ixn_b == n)
		maskdex = np.squeeze(maskdex)
		bmnc = deepcopy(self.bmnc_b)
		b00  = self.bmnc_b[:,0]
		bmnc[:,maskdex] = 0.0
		for i in range(self.ns_b):
			if b00[i] == 0: continue
			error[i] = np.sqrt(np.sum(bmnc[i,:]*bmnc[i,:]))/b00[i]
		return error

	def calcB10B11(self):
		"""Calculates the B10/B11 for each surface

		This routine computes B(m,n) B10/B11 ratio which is
		a proxy for the bootstrap current in the design of W7-X.

		Returns
		-------
		B10/B11 : float
			Ratio of B10/B11
		"""
		import numpy as np
		mask01 = (self.ixm_b == 1) & (self.ixn_b == 0)
		mask11 = (self.ixm_b == 1) & (self.ixn_b == self.nfp_b)
		mask01 = np.squeeze(mask01)
		mask11 = np.squeeze(mask11)
		b01    = self.bmnc_b[:,mask01]
		b11    = self.bmnc_b[:,mask11]
		b11    = np.where(b11==0.0,1.0,b11)
		return np.abs(b01/b11)

	def calcQuasiIsodynamic(self,k,nalpha=65,ntheta0=5,nlambda=4,lplot=False):
		"""Calculates the quasi-isodynamic error

		This routine computes B(m,n) B10/B11 ratio which is
		a proxy for the bootstrap current in the design of W7-X.

		Parameters
		----------
		k : int
			Radial grid index (python indexing)
		nalpha : int
			Gridpoints along magnetic field (default: 65)
		ntheta0 : int
			Number of fieldlines considered (default: 5)
		nlambda : int
			Number of well depths considered (default: 4)
		lplot : boolean
			Make diagnostic plots

		Returns
		-------
		QIerror : float
			Quasi-isodynamic error
		"""
		import numpy as np
		if lplot:
			import matplotlib.pyplot as pyplot
			xplt = np.zeros((nalpha))
			px = 1/pyplot.rcParams['figure.dpi']
			fig,ax = pyplot.subplots(ntheta0,2,figsize=(1024*px,768*px))
		iota0 = self.iota_b[k]
		bvco  = self.bvco_b[k]
		# Helpers
		modb = np.zeros((nalpha))
		modbs = np.zeros((nalpha))
		dl = np.zeros((nalpha))
		integral_A = np.zeros((nalpha))
		integral_B = np.zeros((nalpha))
		J_C = np.zeros((nlambda,ntheta0))
		J_I = np.zeros((nlambda,ntheta0))
		# Loop over thetas
		for i in range(ntheta0):
			deltaphi = np.pi*2.0/iota0
			modbs[:] = 0.0
			modb[:] = 0.0
			theta0  = np.pi*2.0*i/ntheta0
			for l in range(nalpha):
				phi = deltaphi*float(l)/float(nalpha-1)
				theta = theta0 + iota0*phi
				if lplot: xplt[l] = phi
				for mn in range(self.mnboz_b):
					modb[l] = modb[l]+self.bmnc_b[k,mn]*np.cos(self.ixm_b[mn]*theta+self.ixn_b[mn]*phi/self.nfp_b)
			if lplot: 
				ax[i,0].plot(xplt,modb,'k')
				ax[i,0].set_ylabel('|B| [T]')
			# Find min/max values of |B|
			lmin = np.argmin(modb)
			lmax = np.argmax(modb)
			# Recalculate lenght of field line
			phimin = deltaphi*float(lmin)/float(nalpha-1)
			phimax = deltaphi*float(lmax)/float(nalpha-1)
			deltaphi = np.abs(phimax-phimin)
			# Now compute |B| centered around lmin
			modb[:] = 0.0
			for l in range(nalpha):
				phi = phimin - deltaphi + 2*deltaphi*float(l)/float(nalpha-1)
				theta = theta0 + iota0*phi
				if lplot: xplt[l] = phi
				for mn in range(self.mnboz_b):
					modb[l] = modb[l]+self.bmnc_b[k,mn]*np.cos(self.ixm_b[mn]*theta+self.ixn_b[mn]*phi/self.nfp_b)
			if lplot: ax[i,0].plot(xplt,modb,'r')
			# Find the Bmin and half point
			lh = int(np.round(nalpha*0.5))
			lmin = np.argmin(modb)
			lmin = min(max(lmin,1),nalpha-2)
			# Now shift the lmin to the middle
			modb = np.roll(modb,-(lmin-lh))
			# Now recalc lh as lmin
			lmin = np.argmin(modb)
			lmin = min(max(lmin,1),nalpha-2)
			# Squash the array
			modbs[:] = modb[:]
			for l in range(lmin,0,-1):
				if (modbs[l]<modbs[l+1]): modbs[l] = modbs[l+1]
			for l in range(lmin,nalpha,1):
				if (modbs[l]<modbs[l-1]): modbs[l] = modbs[l-1]
			# Stretch the array
			Bmin = np.min(modb)
			Bmax = np.max(modb)
			for l in range(0,lmin+1):
				modbs[l] = Bmin + (Bmax-Bmin)*(modbs[l]-Bmin)/(modbs[0]-Bmin)
			for l in range(lmin,nalpha):
				modbs[l] = Bmin + (Bmax-Bmin)*(modbs[l]-Bmin)/(modbs[nalpha-1]-Bmin)
			if lplot: 
				ax[i,1].plot(xplt,modbs,'k')
				ax[i,1].set_ylabel('|B| [T]')
			# Compute dl
			dl[:] = modb[:]*deltaphi/((nalpha-1)*bvco)
			# Now we evaluate integrals for differnet vlaues of lambda
			for m in range(nlambda):
				Bmir = Bmin + 0.9 * (Bmax-Bmin)*float(m+1)/float(nlambda)
				lam  = 1.0/Bmir
				i1   = np.count_nonzero(modbs[0:lmin+1]>Bmir)
				i2   = np.count_nonzero(modbs[lmin:]<Bmir) + lmin - 1
				i1   = max(i1,0)
				i2   = min(i2,nalpha-1)
				if lplot:
					y1 = np.max(modbs[[i1,i2]])
					ax[i,1].plot(xplt[[i1,i2]],[y1,y1],'r')
				integral_A[:] = 1.0 - lam * modb[:]
				signiA = np.sign(integral_A)
				integral_A[:] = signiA * np.sqrt(np.abs(integral_A[:]))*dl[:]
				#integral_A[:] = np.sqrt(np.abs(1.0 - lam *  modb[:]))*dl[:]
				integral_B[:] = np.sqrt(np.abs(1.0 - lam * modbs[:]))*dl[:]
				J_I[m,i] = np.sum(integral_A[i1:i2])
				J_C[m,i] = np.sum(integral_B[i1:i2])
		if lplot: 
			ax[0,0].set_title('Original')
			ax[0,1].set_title('Approximate')
			ax[-1,0].set_xlabel(r'$\phi_{Boozer}$ [rad]')
			ax[-1,1].set_xlabel(r'$\phi_{Boozer}$ [rad]')
			pyplot.show()
		# Construct Error
		ftemp = 0.0
		for m in range(nlambda):
			for i1 in range(ntheta0):
				for i2 in range(ntheta0):
					ftemp = ftemp + J_I[m,i1] - J_C[m,i2]
		norm = np.sum(J_C+J_I)/float(nlambda*ntheta0)
		return ftemp/norm










# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)









