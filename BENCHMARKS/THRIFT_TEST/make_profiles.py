#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a class to help setup a THRIFT
profiles file. The output file itself has arrays
raxis_prof nrho size (r/a)
taxis_prof ntime size (time)
te_prof ntime,nrho size (Te)
ne_prof ntime,nrho size (ne)
Z nion size (Z)
M nion size (mass)
ti_prof ntime,nrho,nion size (Ti)
ni_prof ntime,nrho,nion size (ni)
"""
# Libraries
import sys, os
#sys.path.insert(0, '../../pySTEL/')
import numpy as np                    #For Arrays
#from math import pi
#from libstell.beams3d import read_beams3d



class THRIFT_PROFILES():
	"""Class which handles performing calorimetry of the NBI system

	"""
	def __init__(self):
		self.amu = 1.66053906660E-27
		self.rho = None
		self.time = None
		self.te = None
		self.ne = None
		self.Z  = None
		self.M  = None
		self.ni = None
		self.ti = None

	def set_rho(self,rho):
		"""Set the rho profile

		Sets the rho profile array.  Here rho is the normalized
		radial coordinates <r/a> where the normalized toroidal
		flux is s = Phi/Phi_edge = rho^2

		Paramters
		---------
		rho : real
			rho <r/a>
		"""
		self.rho = rho

	def set_time(self,time):
		"""Set the time profile

		Sets the time profile array

		Paramters
		---------
		time : real
			time data [s]
		"""
		self.time = time

	def set_te(self,te):
		"""Set the Te profile

		Sets the Te profile first dimenion must be the same size as
		the time array, while the second must be the same size as
		the rho array.

		Paramters
		---------
		te : real
			te data [ev]
		"""
		self.te = te

	def set_ne(self,ne):
		"""Set the ne profile

		Sets the ne profile first dimenion must be the same size as
		the time array, while the second must be the same size as
		the rho array.

		Paramters
		---------
		ne : real
			ne data [m^-3]
		"""
		self.ne = ne

	def set_Z(self,Z):
		"""Set the ion Z array

		Sets ion Z array which is the Z for each ion species.

		Paramters
		---------
		Z : real
			Z (electron charge)
		"""
		self.Z = Z

	def set_mass(self,M):
		"""Set the ion mass array

		Sets ion mass array which is the mass for each ion species.

		Paramters
		---------
		M : real
			M [amu]
		"""
		self.M = M*self.amu

	def set_ti(self,ti):
		"""Set the Ti profile

		Sets the Ti profile first dimenion must be the same size as
		the Z array and M array, second must be the same size as the
		time array and third must be the same size as the rho array.

		Paramters
		---------
		ti : real
			ti data [eV]
		"""
		self.ti = ti

	def set_ni(self,ni):
		"""Set the ni profile

		Sets the ni profile first dimenion must be the same size as
		the Z array and M array, second must be the same size as the
		time array and third must be the same size as the rho array.

		Paramters
		---------
		ni : real
			ni data [m^-3]
		"""
		self.ni = ni

	def write_profile(self,name='profiles_example.h5'):
		"""Writes the profiles HDF5 file

		Writes the profiles HDF5 file.

		Paramters
		---------
		name : string
			File name
		"""
		import h5py
		nrho = len(self.rho)
		nt   = len(self.time)
		nZ   = len(self.Z) 
		nM   = len(self.M)
		nni  = self.ni.shape[2]
		nti  = self.ti.shape[2]
		if not(nZ == nM) or not(nZ == nni) or not(nZ == nti):
			print(' Check ion dimensionality')
			return 
		hf = h5py.File(name, 'w')
		hf.create_dataset('nrho', data=nrho)
		hf.create_dataset('nt', data=nt)
		hf.create_dataset('nion', data=nZ)
		hf.create_dataset('raxis_prof', data=self.rho)
		hf.create_dataset('taxis_prof', data=self.time)
		hf.create_dataset('Z_prof', data=self.Z)
		hf.create_dataset('mass_prof', data=self.M)
		hf.create_dataset('ne_prof', data=self.ne)
		hf.create_dataset('te_prof', data=self.te)
		hf.create_dataset('ni_prof', data=self.ni)
		hf.create_dataset('ti_prof', data=self.ti)
		hf.close()

	def generate_profiles_example(self):
		"""Example of how to generate profiles
		"""
		f_ne = lambda x,nemax,nemin: nemax - (nemax-nemin)*x**3
		f_te = lambda x,temax,temin: temax - (temax-temin)*x
		f_ti = lambda x,timax,timin: timax - (timax-timin)*x
		f_time = lambda x,tau: 1-np.exp(-x/tau)
		mass_array = np.array([1, 4, 12])
		Z_array    = np.array([1, 2,  6])
		frac_array = np.array([0, 0.01, 0.01])
		time = np.linspace(0,10.0,1000)
		rho  = np.linspace(0,1.0,64)
		ne   = np.zeros((1000,64))
		te   = np.zeros((1000,64))
		ni   = np.zeros((1000,64,3))
		ti   = np.zeros((1000,64,3))
		for i,t in enumerate(time):
			for j,r in enumerate(rho):
				ne[i,j] = f_time(t,2)*f_ne(r,5E19,1E19)
				te[i,j] = f_time(t,2)*f_te(r,3000,100)
				ti[i,j,0] = f_time(t,2)*f_ti(r,1500,100)
				ni[i,j,0] = ne[i,j]
		if ni.shape[0] > 1:
			for k,Z in enumerate(Z_array):
				if k == 0: continue
				for i,t in enumerate(time):
					for j,r in enumerate(rho):
						ti[i,j,k] = ti[i,j,0]
						ni[i,j,k] = ne[i,j]*frac_array[k]
						ni[i,j,k] = ni[i,j,0] - ni[i,j,k]*Z
		return rho, time, mass_array,Z_array, ne, te, ni, ti


if __name__ == "__main__":
	from argparse import ArgumentParser
	profs = THRIFT_PROFILES()
	[rho, time, mass_array,Z_array, ne, te, ni, ti] = profs.generate_profiles_example()
	profs.set_rho(rho)
	profs.set_time(time)
	profs.set_Z(Z_array)
	profs.set_mass(mass_array)
	profs.set_ne(ne)
	profs.set_te(te)
	profs.set_ti(ti)
	profs.set_ni(ni)
	profs.write_profile('profiles_test.h5')
	sys.exit(0)