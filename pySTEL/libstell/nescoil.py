#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for reading and handling NESCOIL
surface potential data.
"""

# Libraries
from libstell.libstell import LIBSTELL, FourierRep

# Constants

# NESCOIL Class
class NESCOIL(FourierRep):
	"""Class for working with NESCOIL output data.

	"""
	def __init__(self):
		super().__init__()
		self.libStell = LIBSTELL()

	def read_nescout(self,filename):
		"""Reads a NESCOIL output file.

		This routine reads and initilizes the NESCOIL
		class with information from a nescout file.

		Parameters
		----------
		file : str
			Path to nescout file.
		"""
		nescout_dict = self.libStell.read_nescout(filename)
		for key in nescout_dict:
			setattr(self, key, nescout_dict[key])

	def bfield(self,x,y,z):
		"""Evaluates the magnetic field from a NESCOIL surface current

		This routine evaluates the magnetic field at a point in space
		from a NESCOIL surface current.

		Parameters
		----------
		x : float
			X point to evaluate [m]
		y : float
			Y point to evaluate [m]
		z : float
			Z point to evaluate [m]
		Returns
		----------
		bx : float
			X-component of magnetic field [T]
		by : float
			Y-component of magnetic field [T]
		bZ : float
			Z-component of magnetic field [T]
		"""
		return self.libStell.nescout_bfield(x,y,z,istat=-327)

	def generatePotential(self,theta,zeta):
		"""Computes the potential on a grid

		This routine computes the potential on a grid

		Parameters
		----------
		theta : ndarray
			Poloidal Angle [rad]
		zeta : ndarray
			Field Period Angle [rad]
		Returns
		-------
		pot : ndarray
			NESCOIL potential []
		"""
		pot = self.sfunct(theta,zeta,self.potmns_surface.T,self.xm_pot,self.xn_pot)
		return pot

	def generateTotalPotential(self,theta,zeta):
		"""Computes the potential on a grid

		This routine computes the potential on a grid

		Parameters
		----------
		theta : ndarray
			Poloidal Angle [rad]
		zeta : ndarray
			Field Period Angle [rad]
		Returns
		-------
		pot : ndarray
			NESCOIL potential []
		"""
		import numpy as np
		pot = self.generatePotential(theta,zeta)
		nu = len(theta)
		nv = len(zeta)
		for j in range(self.nu): pot[0,j,:] = pot[0,j,:] - self.cut*0.5*theta[j]/np.pi
		for j in range(self.nv): pot[0,:,j] = pot[0,:,j] - self.cup*0.5*zeta[j]/(np.pi*self.np)
		return pot


	def plotpotential(self,ax=None):
		"""Plots the NESCOIL Potential

		This routine plots the NESCOIL code surface potential

		Parameters
		----------
		ax : axes (optional)
			Matplotlib axes object to plot to.
		"""
		import numpy as np
		import matplotlib.pyplot as pyplot
		lplotnow = False
		if not ax:
			ax = pyplot.axes()
			lplotnow = True
		theta = np.ndarray((self.nu,1))
		zeta  = np.ndarray((self.nv,1))
		for j in range(self.nu): theta[j]=2.0*np.pi*j/float(self.nu-1)
		for j in range(self.nv):  zeta[j]=    np.pi*j/float(self.nv-1)
		pot = self.generatePotential(theta,zeta)
		#pot = self.sfunct(theta,zeta,self.potmns_surface.T,self.xm_pot,self.xn_pot)
		hmesh=ax.pcolormesh(np.squeeze(zeta),np.squeeze(theta),np.squeeze(pot[0,:,:]),cmap='jet',shading='gouraud')
		ax.set_xlabel('Toroidal angle [rad]')
		ax.set_ylabel('Poloidal angle [rad]')
		ax.set_title(r'NESCOIL $\Phi$ Potential')
		pyplot.colorbar(hmesh,label='$Pot$ [arb]',ax=ax)

	def plottotalpotential(self,ax=None):
		"""Plots the NESCOIL Total Potential

		This routine plots the NESCOIL code surface potential

		Parameters
		----------
		ax : axes (optional)
			Matplotlib axes object to plot to.
		"""
		import numpy as np
		import matplotlib.pyplot as pyplot
		lplotnow = False
		if not ax:
			ax = pyplot.axes()
			lplotnow = True
		theta = np.ndarray((self.nu,1))
		zeta  = np.ndarray((self.nv,1))
		for j in range(self.nu): theta[j]=2.0*np.pi*j/float(self.nu-1)
		for j in range(self.nv):  zeta[j]=    np.pi*j/float(self.nv-1)
		pot = self.generateTotalPotential(theta,zeta)
		#pot = self.sfunct(theta,zeta,self.potmns_surface.T,self.xm_pot,self.xn_pot)
		#for j in range(self.nu): pot[0,j,:] = pot[0,j,:] - self.cut*1.0*j/float(self.nv-1)
		#for j in range(self.nv): pot[0,:,j] = pot[0,:,j] - self.cup*0.5*j/float(self.nv-1)
		hmesh=ax.pcolormesh(np.squeeze(zeta),np.squeeze(theta),np.squeeze(pot[0,:,:]),cmap='jet',shading='gouraud')
		#ax.semilogy(abscissa, data, **kwargs)
		ax.set_xlabel('Toroidal angle [rad]')
		ax.set_ylabel('Poloidal angle [rad]')
		ax.set_title(r'NESCOIL Total $\Phi$ Potential')
		pyplot.colorbar(hmesh,label='$Pot$ [arb]',ax=ax)
		if lplotnow: pyplot.show()

	def plotsurfaces(self,ax=None):
		"""Plots the NESCOIL Surfaces

		This routine plots the NESCOIL current potential surface
		and the plasma surface over a half field period.

		Parameters
		----------
		ax : axes (optional)
			Matplotlib axes object to plot to.
		"""
		import numpy as np
		import matplotlib.pyplot as pyplot
		lplotnow = False
		if not ax:
			#fig = pyplot.figure()
			#ax = fig.add_subplot(111,projection='3d')
			lplotnow = True
		theta = np.ndarray((self.nu,1))
		zeta  = np.ndarray((self.nv,1))
		for j in range(self.nu): theta[j]=2.0*np.pi*j/float(self.nu-1)
		for j in range(self.nv):  zeta[j]=np.pi*j/float(self.nv-1)
		rp = self.cfunct(theta,zeta,self.rmnc_plasma.T,self.xm_plasma,self.xn_plasma)
		zp = self.sfunct(theta,zeta,self.zmns_plasma.T,self.xm_plasma,self.xn_plasma)
		rc = self.cfunct(theta,zeta,self.rmnc_surface.T,self.xm_surface,self.xn_surface)
		zc = self.sfunct(theta,zeta,self.zmns_surface.T,self.xm_surface,self.xn_surface)
		ax = self.isotoro(rp,zp,zeta,0,color='red')
		hc = self.isotoro(rc,zc,zeta,0,axes=ax,color='green')
		#ax.set_aspect('equal', 'box')
		#ax.semilogy(abscissa, data, **kwargs)
		if lplotnow: pyplot.show()

	def cutcoils(self,ncoils_per_halfperiod):
		"""Cut coils from the NESCOIL potential

		This routine cuts coils from the NESCOIL potential.
		It allows the user to specify the number of coils per half 
		period.

		Parameters
		----------
		ncoils_per_halfperiod : integer
			Number of coils per half period (suggest 5)
		"""
		import numpy as np
		from libstell.coils import COILSET
		# First calculate the potential on a grid



# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)