#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for reading and handling 
BEAMS3D Energetic particle data.
"""

# Libraries

# Constants

# FIELDLINES Class
class BEAMS3D():
	"""Class for working with BEAMS3D data

	"""
	def __init__(self):
		test = 0

	def read_beams3d(self,file):
		"""Reads a BEAMS3D HDF5 file

		This routine reads and initilizes the BEAMS3D
		class with variable information from an HDF5 file.

		Parameters
		----------
		file : str
		Path to HDF5 file.
		"""
		import h5py
		import numpy as np
		beams_data = {}
		#read file
		with h5py.File(file,'r') as f:
			# Logicals
			for temp in ['lascot', 'lbeam', 'lbeam_simple', \
						 'lcoil', 'lcollision', 'ldepo', 'lflux', \
						 'lmgrid', 'lpies', 'lspec', 'lvac', \
						 'lvessel', 'lvmec', 'leqdsk', 'lhint', \
						 'lfusion', 'lboxsim', 'lhitonly']:
				if temp in f:
					setattr(self, temp, np.int64(f[temp][0])==1)
			# Integers
			for temp in ['nbeams', 'nface', 'nparticles', 'nphi', \
						 'npoinc', 'nz', 'nsteps', 'nvertex', 'nz', \
						 'ns_prof1', 'ns_prof2', 'ns_prof3', \
						 'ns_prof4', 'ns_prof5']:
				if temp in f:
					setattr(self, temp, np.int64(f[temp][0]))
			# Floats
			for temp in ['VERSION','plasma_mass','plasma_Zmean', \
						 'B_kick_min', 'B_kick_max', 'freq_kick', \
						 'E_kick', 'partvmax']:
				if temp in f:
					setattr(self, temp, float(f[temp][0]))
			# Arrays
			for temp in ['raxis', 'phiaxis', 'zaxis', 'B_R', \
						 'B_PHI', 'B_Z', 'S_ARR', 'U_ARR', 'POT_ARR', \
						 'GFactor', 'TE', 'NE', 'NI', 'TI', \
						 'ZEFF_ARR', 'wall_vertex', 'wall_faces', \
						 'mass', 'charge', 'Weight', 'Beam', 'Zatom', \
						 'end_state', 'Energy', 'wall_strikes', \
						 'wall_load', 'wall_shine', 'beam_density', \
						 't_end', 'R_lines', 'Z_lines', 'PHI_lines', \
						 'vll_lines', 'neut_lines', 'moment_lines', \
						 'S_lines', 'U_lines', 'B_lines', \
						 'dist_rhoaxis', 'dist_uaxis', 'dist_paxis', \
						 'dist_Vaxis', 'dist_Waxis', 'ndot_prof', \
						 'epower_prof', 'ipower_prof', 'j_prof', \
						 'dense_prof', 'dist_prof', 'Shinethrough', \
						 'Shineport', 'NEUTRON_RATE', 'E_NEUTRONS']:
				if temp in f:
					setattr(self, temp, np.array(f[temp][:]))
		# Make derived arrays
		self.X_lines = self.R_lines*np.cos(self.PHI_lines)
		self.Y_lines = self.R_lines*np.sin(self.PHI_lines)
		self.MODB    = np.sqrt(self.B_R**2 + self.B_PHI**2 + self.B_Z**2)
		return

	def calcBaxis(self,kphi=1):
		"""Calculates the magnetic field on axis

		This routine calcualtes the magnetic field on axis
		given a toroidal cut in the background grid. Default is
		to use the first toroidald gridpoint (phi=0).

		Parameters
		----------
		kphi : int (optional)
			Toroidal plane to evaluate
		
		Returns
		-------
		Baxis : float
			Magnetic field on axis.
		"""
		import numpy as np
		S2D = np.squeeze(self.S_ARR[:,kphi,:])
		B2D = np.squeeze(self.MODB[:,kphi,:])
		dex = np.argwhere(S2D == np.min(S2D))
		return B2D[dex[0],dex[1]]

	def calcAminor(self):
		"""Calculates the minor radius

		This routine calculates the minor radius using the
		normalized toroidal flux label (S).
		
		Returns
		-------
		Aminor : float
			Equilibrium minor radius
		"""
		import numpy as np
		from contourpy import contour_generator, LineType
		Aminor = 0
		for k in range(self.nphi): 
			S2D = np.squeeze(self.S_ARR[:,k,:])
			S2D = np.where(S2D>1.2,1.2,S2D)
			dex = np.argwhere(S2D == np.min(S2D))
			r0  = self.raxis[dex[0]]
			z0  = self.zaxis[dex[1]]
			cont_gen = contour_generator(x=self.raxis,y=self.zaxis,z=S2D, line_type=LineType.Separate)
			lines = cont_gen.lines(1.0)
			Aminor = Aminor + a
		return Aminor/self.nphi

	def calcVolume(self):
		"""Calculates the differential volume

		This routine calcualtes the differential volume and the volume
		as a function of normalized toroidal flux (S). 
		
		Returns
		-------
		S : ndarray
			Normalized toroidal flux array
		V : ndarray
			Volume [m^3]
		dVds : ndarray
			Differential volume dV/ds [m^3]
		"""
		import numpy as np
		dr = self.raxis[1] - self.raxis[0]
		dz = self.zaxis[1] - self.zaxis[0]
		dp = self.phiaxis[1] - self.phiaxis[0]
		nfp = np.round(2*np.pi/self.phiaxis[-1])

		area = dr*dz
		vol  = self.raxis*dphi*area
		# TBD
		return -1,-1,-1



