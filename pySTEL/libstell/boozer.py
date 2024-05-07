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
		boozmn_dict = self.libStell.read_boozer(filename)
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

	def plotBsurf(self,sval,ax=None):
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
		theta = np.ndarray((360,1))
		zeta  = np.ndarray((256,1))
		for j in range(360): theta[j]=2.0*np.pi*j/359.0
		for j in range(256):  zeta[j]=2.0*np.pi*j/256.0
		b = self.cfunct(theta,zeta,self.bmnc_b,self.ixm_b,self.ixn_b/self.nfp_b)
		hmesh=ax.pcolormesh(np.squeeze(zeta),np.squeeze(theta),np.squeeze(b[sval,:,:]),cmap='jet',shading='gouraud')
		ax.plot(zeta,zeta*self.iota_b[sval],'w')
		ax.set_xlabel(r'Toroidal Angle ($\phi$) [rad]')
		ax.set_ylabel(r'Poloidal Angle ($\theta$) [rad]')
		ax.set_title(rf'BOOZER $s$ = {sval/self.ns_b:5.2f}')
		pyplot.colorbar(hmesh,label='[T]',ax=ax)
		if lplotnow: pyplot.show()



# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)









