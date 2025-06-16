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
			self.bnmns[0,mn] = float(txt3)
			mn = mn + 1

	def read_bnorm_real(self,filename):
		"""Reads a BNORM_REAL file

		This routine reads and initilizes the BNORM class
		with variable information from a BNORM_REAL file.

		Parameters
		----------
		file : str
			Path to wout file.
		"""
		import numpy as np
		f = open(filename,'r')
		lines = f.readlines()
		f.close()
		self.nuv=int(lines[0])
		u             = np.zeros(self.nuv,dtype=np.int64)
		v             = np.zeros(self.nuv,dtype=np.int64)
		theta         = np.zeros(self.nuv)
		zeta          = np.zeros(self.nuv)
		phi           = np.zeros(self.nuv)
		rreal         = np.zeros(self.nuv)
		zreal         = np.zeros(self.nuv)
		Nx            = np.zeros(self.nuv)
		Ny            = np.zeros(self.nuv)
		Nz            = np.zeros(self.nuv)
		bnreal        = np.zeros(self.nuv)
		bcreal        = np.zeros(self.nuv)
		bnormal_total = np.zeros(self.nuv)
		for j in range(self.nuv):
			txt                   = lines[j+1].split()
			u[j]             = int(txt[1])
			v[j]             = int(txt[2])
			theta[j]         = float(txt[3])
			zeta[j]          = float(txt[4])
			phi[j]           = float(txt[5])
			rreal[j]         = float(txt[6])
			zreal[j]         = float(txt[7])
			Nx[j]            = float(txt[8])
			Ny[j]            = float(txt[9])
			Nz[j]            = float(txt[10])
			bnreal[j]        = float(txt[11])
			bcreal[j]        = float(txt[12])
			bnormal_total[j] = float(txt[13])
		nu = max(u)
		nv = max(v)
		self.u             = u.reshape(nu,nv)
		self.v             = v.reshape(nu,nv)
		self.theta          = theta.reshape(nu,nv)
		self.zeta          = zeta.reshape(nu,nv)
		self.phi           = phi.reshape(nu,nv)
		self.rreal         = rreal.reshape(nu,nv)
		self.zreal         = zreal.reshape(nu,nv)
		self.Nx            = Nx.reshape(nu,nv)
		self.Ny            = Ny.reshape(nu,nv)
		self.Nz            = Nz.reshape(nu,nv)
		self.bnreal        = bnreal.reshape(nu,nv)
		self.bcreal        = bcreal.reshape(nu,nv)
		self.bnormal_total = bnormal_total.reshape(nu,nv)

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
			bmn[m,n] = self.bnmns[1,mn]
		#Plot
		x = np.linspace(0,mmax,mmax+1)
		y = np.linspace(-nmax,nmax,2*nmax+1)
		hmesh=ax.pcolormesh(x,y,np.log10(np.abs(bmn.T)),cmap='jet',shading='gouraud')
		ax.set_xlabel('Poloidal Modes (m)')
		ax.set_ylabel('Toroidal Modes (n)')
		ax.set_title(rf'BNORM Normal Field')
		pyplot.colorbar(hmesh,label='$log_{10}$[arb]',ax=ax)
		if lplotnow: pyplot.show()

	def plot_bnorm_real_total(self,ax=None):
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
		x = self.theta[:,0]
		y = self.phi[0,:]
		hmesh=ax.pcolormesh(y,x,self.bnormal_total.T)
		ax.set_xlabel('Toroidal Angle (phi) [rad]')
		ax.set_ylabel('Poloidal Angle (phi) [rad]')
		ax.set_title(rf'Total B-Normal Field')
		#pyplot.colorbar(hmesh,label='$log_{10}$[arb]',ax=ax)
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
		b = self.sfunct(theta,zeta,self.bnmns,self.xm,self.xn)
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









