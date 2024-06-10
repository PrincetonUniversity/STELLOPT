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
		return self.libStell.nescout_bfield(x,y,z)

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
		for j in range(self.nv):  zeta[j]=2.0*np.pi*j/float(self.nv-1)
		print(theta.shape)
		pot = self.cfunct(theta,zeta,self.potmnc_surface,self.xm_surface,self.xn_surface)
		#hmesh=ax.pcolormesh(theta,zeta,np.squeeze(pot[1,:,:]),cmap='jet',shading='gouraud')
		hmesh=ax.pcolormesh(np.squeeze(pot[1,:,:]),cmap='jet',shading='gouraud')
		#ax.semilogy(abscissa, data, **kwargs)
		ax.set_xlabel('Toroidal angle [rad]')
		ax.set_ylabel('Poloidal angle [rad]')
		ax.set_title('NESCOIL Potential')
		pyplot.colorbar(hmesh,label='$Pot$ [arb]',ax=ax)
		if lplotnow: pyplot.show()



# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)