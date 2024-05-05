##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for reading and handling BOOZER
transform data.
"""

# Libraries
#from libstell import libstell
import libstell

# Constants

# VMEC Class
class BOOZER(libstell.FourierRep):
	"""Class for working with VMEC equilibria

	"""
	def __init__(self):
		super().__init__()
		self.nfp = None
		self.libStell = libstell.LIBSTELL()

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
		fig.colorbar(hmesh,label='$log_{10}$[T]')
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
		fig.colorbar(hmesh,label='[T]')
		if lplotnow: pyplot.show()



# Main routine
if __name__=="__main__":
	import sys
	from argparse import ArgumentParser
	import numpy as np
	import matplotlib.pyplot as pyplot
	parser = ArgumentParser(description= 
		'''Provides class for accessing boozer data also servers as a
		   simple tool for assessing boozer boozmn files.''')
	parser.add_argument("-b", "--boozer", dest="booz_ext",
		help="BOOZER boozmn file extension", default = None)
	parser.add_argument("-p", "--plot", dest="lplot", action='store_true',
		help="Plot the boozer file.", default = False)
	parser.add_argument("-s", "--surf", dest="sval",
		help="Equilibrium surface to plot.", default = None)
	boozmn_data = BOOZER()
	args = parser.parse_args()
	if args.booz_ext:
		#try:
		boozmn_data.read_boozer(args.booz_ext)
		#except:
		#	print(f'Could not file boozmn file: {args.booz_ext}')
		#	sys.exit(-1)
		if args.lplot:
			if not args.sval:
				sval = int(0.25*boozmn_data.ns_b)
			else:
				sval = int(args.sval)
			px = 1/pyplot.rcParams['figure.dpi']
			fig,(ax1,ax2) = pyplot.subplots(1,2,figsize=(1024*px,768*px))
			pyplot.subplots_adjust(hspace=0.1,wspace=0.3)
			boozmn_data.plotBmnSpectrum(sval,ax=ax1)
			boozmn_data.plotBsurf(sval,ax=ax2)
			pyplot.show()








