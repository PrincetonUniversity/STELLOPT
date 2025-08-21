##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for reading and handling COBRAVMEC
data.
"""

# Libraries
from libstell.libstell import LIBSTELL, FourierRep

# Constants

# COBRA Class
class COBRA(FourierRep):
	"""Class for working with COBRA equilibria

	"""
	def __init__(self):
		super().__init__()
		self.libStell = LIBSTELL()

	def read_cobra(self,filename):
		"""Reads a COBRA file

		This routine reads and initilizes the COBRA class
		with variable information from a COBRA_GRATE file.

		Parameters
		----------
		file : str
			Path to wout file.
		"""
		import numpy as np
		f = open(filename,'r')
		lines = f.readlines()
		f.close()
		th = []
		ze = []
		s  = []
		k  = []
		grate = []
		i = 0
		while (i < len(lines)):
			(ze_txt,th_txt,nsurf_txt) = lines[i].split()
			th_flt = float(th_txt)
			ze_flt = float(ze_txt)
			nsurf_int = int(nsurf_txt)
			for j in range(nsurf_int):
				i = i + 1
				(k_txt,s_txt,grate_txt) = lines[i].split()
				th.append(th_flt)
				ze.append(ze_flt)
				k.append(int(k_txt))
				s.append(float(s_txt))
				grate.append(float(grate_txt))
			i = i + 1
		theta = np.unique(th)
		zeta  = np.unique(ze)
		ks    = np.unique(k)
		st    = np.unique(s)
		nth   = len(theta)
		nze   = len(zeta)
		ns    = len(ks)
		self.theta = theta
		self.zeta = zeta
		self.s = st
		self.k = ks
		self.grate = np.reshape(grate,(nth,nze,ns))


	def plot_grate(self,ntheta=None,nzeta=None,ax=None):
		""" Plot the ballooning growth rate

		This routine plots the ballooning growth rate as a function of the radial grid.


		Parameters
		----------
		ntheta : list (int)
			Poloidal indices to plot (default: all)
		nzeta : list (int)
			Poloidal indices to plot (default: all)
		ax : axes (optional)
			Matplotlib axes object to plot to.
		"""
		import numpy as np
		import matplotlib.pyplot as pyplot
		lplotnow = False
		# Handles axes
		if not ax:
			ax = pyplot.axes()
			lplotnow = True
		# Handle ntheta
		if type(ntheta) == type(None):
			th_vec = np.linspace(0,len(self.theta)-1,len(self.theta),dtype=int)
		else:
			th_vmec = np.flatten([ntheta])
		# Handle ntheta
		if type(nzeta) == type(None):
			ze_vec = np.linspace(0,len(self.zeta)-1,len(self.zeta),dtype=int)
		else:
			ze_vmec = np.flatten([ntheta])
		# Plot
		for i in th_vec:
			for j in ze_vec:
				ax.plot(self.s,self.grate[i,j,:],label=rf'$\theta$: {self.theta[i]:5.2f}, $\zeta$: {self.zeta[j]:5.2f},')
		ax.set_xlabel('Norm. Toridal Flux (s)')
		ax.set_ylabel(r'Growth Rate $(-\gamma^2)$')
		ax.set_title(rf'COBRA Ballooning Stability')
		ax.legend()
		if lplotnow: pyplot.show()

# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)