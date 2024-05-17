##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for reading and handling BNORM
data.
"""

# Libraries
from libstell.libstell import LIBSTELL, FourierRep

# Constants

# VMEC Class
class BNORM(FourierRep):
	"""Class for working with VMEC equilibria

	"""
	def __init__(self):
		super().__init__()
		self.libStell = LIBSTELL()

	def read_bnorm(self,filename):
		"""Reads a BNORM file

		This routine reads and initilizes the BNORM class
		with variable information from a BNORM file.

		Parameters
		----------
		file : str
			Path to wout file.
		"""
		import numpy as np
		f = open(filename,'r')
		lines = f.readlines()
		f.close()
		self.mnmax = len(lines)
		self.xm = np.zeros((self.mnmax))
		self.xn = np.zeros((self.mnmax))
		self.bnmnc = np.zeros((1,self.mnmax))
		self.bnmns = np.zeros((1,self.mnmax))
		mn = 0
		for line in lines:
			(txt1,txt2,txt3) = line.split()
			self.xm[mn] = int(txt1)
			self.xn[mn] = int(txt2)
			self.bnmnc[0,mn] = float(txt3)
			mn = mn + 1

	def plotBnmnSpectrum(self,ax=None):
		"""Plots the Bnormal spectrum for a surface

		This routine plots the bnormal spectrum for a given
		surface.

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
		# Array extents
		mmax = int(max(np.squeeze(self.xm)))
		nmax = int(max(np.squeeze(self.xn)))
		# Sort BMN into array
		bmn = np.zeros((mmax+1,2*nmax+1))
		for mn in range(self.mnmax):
			m = self.xm[mn]
			n = int(self.xn[mn]) + nmax
			bmn[m,n] = self.bmnc_b[1,mn]
		#Plot
		x = np.linspace(0,mmax,mmax+1)
		y = np.linspace(-nmax,nmax,2*nmax+1)
		hmesh=ax.pcolormesh(x,y,np.log10(np.abs(bmn.T)),cmap='jet',shading='gouraud')
		ax.set_xlabel('Poloidal Modes (m)')
		ax.set_ylabel('Toroidal Modes (n)')
		ax.set_title(rf'BNORM Normal Field')
		pyplot.colorbar(hmesh,label='$log_{10}$[arb]',ax=ax)
		if lplotnow: pyplot.show()

	def plotBsurf(self,ax=None):
		"""Plots the Bnormal on a surface

		This routine plots the bnormal.

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
		theta = np.ndarray((360,1))
		zeta  = np.ndarray((256,1))
		for j in range(360): theta[j]=2.0*np.pi*j/359.0
		for j in range(256):  zeta[j]=2.0*np.pi*j/256.0
		b = self.cfunct(theta,zeta,self.bnmnc,self.xm,self.xn)
		hmesh=ax.pcolormesh(np.squeeze(zeta),np.squeeze(theta),np.squeeze(b[1,:,:]),cmap='jet',shading='gouraud')
		ax.set_xlabel(r'Toroidal Angle ($\phi$) [rad]')
		ax.set_ylabel(r'Poloidal Angle ($\theta$) [rad]')
		ax.set_title(rf'BNORM')
		pyplot.colorbar(hmesh,label=r'$B_{normal}$ [arb]',ax=ax)
		if lplotnow: pyplot.show()



# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)









