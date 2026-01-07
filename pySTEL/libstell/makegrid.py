##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for reading and handling MAKEGRID
data.
"""

# Libraries
from libstell.libstell import LIBSTELL

# Constants

# MGRID Class
class MAKEGRID():
	"""Class for working with MGRID vacuum field data.

	"""
	def __init__(self):
		super().__init__()
		self.libStell = LIBSTELL()

	def read_mgrid(self,filename,extcur,nv,nfp):
		"""Reads a MAKEGRID file

		This routine reads and initilizes the MAKEGRID class
		with variable information from a MAKEGRID file.

		Parameters
		----------
		file : str
			Path to makegrid file.
		extcur : extcur
			List of external currents ([A] or scale factor)
		nv : int
			Number of toroidal planes
		nfp : int
			Field periodicty
		"""
		import copy
		import numpy as np
		mgrid_dict = self.libStell.read_mgrid(filename,extcur,nv,nfp)
		for key in mgrid_dict:
			setattr(self, key, mgrid_dict[key])
		self.brvac = np.transpose(self.brvac,(2,1,0))
		self.bpvac = np.transpose(self.bpvac,(2,1,0))
		self.bzvac = np.transpose(self.bzvac,(2,1,0))

	def plot_bfield(self,nv=1,ax=None):
		"""Plot a toroidal cut of the magnetic field

		This routine plots a toroidal cut of the magnetic field.

		Parameters
		----------
		nv : int
			Toroidal index of plot.
		ax : axes (optional)
			Matplotlib axes object to plot to.
		"""
		import numpy as np
		import matplotlib.pyplot as pyplot
		lplotnow = False
		if not ax:
			ax = pyplot.axes()
			lplotnow = True
		#Plot
		x = np.linspace(self.rminb,self.rmaxb,self.nr0b)
		y = np.linspace(self.zminb,self.zmaxb,self.nz0b)
		b = np.sqrt(self.brvac**2+self.bzvac**2+self.bpvac**2)
		hmesh=ax.pcolormesh(x,y,np.squeeze(b[:,:,nv]).T,cmap='jet',shading='gouraud')
		ax.set_xlabel('R [m]')
		ax.set_ylabel('Z [m]')
		ax.set_title(rf'MAKEGRID Magnetic Field')
		pyplot.colorbar(hmesh,label='|B|',ax=ax)
		if lplotnow: pyplot.show()

# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)
