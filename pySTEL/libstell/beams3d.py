#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for reading and handling 
BEAMS3D Energetic particle data.
"""

# Libraries
from libstell.libstell import LIBSTELL

# Constants

# BEAMS3D Class
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
						 'npoinc', 'nr', 'nsteps', 'nvertex', 'nz', \
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
			# Arrays (1D)
			for temp in ['raxis', 'phiaxis', 'zaxis', \
						 'mass', 'charge', 'Weight', 'Beam', 'Zatom', \
						 'end_state', 'Energy', 'wall_strikes', \
						 'GFactor', 't_end',  \
						 'wall_strikes', \
						 'dist_rhoaxis', 'dist_uaxis', 'dist_paxis', \
						 'dist_Vaxis', 'dist_Waxis', 'Shinethrough', \
						 'Shineport', 'Energy']:
				if temp in f:
					setattr(self, temp, np.array(f[temp][:]))
			# Arrays (2D)
			for temp in ['wall_vertex', 'wall_faces', \
						 'wall_load', 'wall_shine', 'beam_density', \
						 'R_lines', 'Z_lines', 'PHI_lines', \
						 'vll_lines', 'neut_lines', 'moment_lines', \
						 'S_lines', 'U_lines', 'B_lines', \
						 'ndot_prof', \
						 'epower_prof', 'ipower_prof', 'j_prof', \
						 'dense_prof', 'E_NEUTRONS']:
				if temp in f:
					array = np.transpose(f[temp][:],(1,0))
					setattr(self, temp, np.array(array))
			# Arrays (3D)
			for temp in ['B_R', \
						 'B_PHI', 'B_Z', 'S_ARR', 'U_ARR', 'POT_ARR', \
						 'TE', 'NE', 'TI', \
						 'ZEFF_ARR']:
				if temp in f:
					array = np.transpose(f[temp][:],(2,1,0))
					setattr(self, temp, np.array(array))
			# Arrays (4D)
			for temp in ['NEUTRON_RATE','beam_density', 'NI']:
				if temp in f:
					array = np.transpose(f[temp][:],(3,2,1,0))
					setattr(self, temp, np.array(array))
			# Arrays (6D)
			for temp in ['NEUTRON_RATE','beam_density']:
				if temp in f:
					array = np.transpose(f[temp][:],(5,4,3,2,1,0))
					setattr(self, temp, np.array(array))
		# Make derived arrays
		self.X_lines = self.R_lines*np.cos(self.PHI_lines)
		self.Y_lines = self.R_lines*np.sin(self.PHI_lines)
		self.MODB    = np.sqrt(self.B_R**2 + self.B_PHI**2 + self.B_Z**2)
		if hasattr(self,'wall_faces'): self.wall_faces = self.wall_faces - 1
		return

	def calcVperp(self):
		"""Calculates the perpendicular marker velocity

		This routine calcualtes the perpendicular marker velocity
		from the magnetic moment.
		
		Returns
		-------
		Vperp : float
			Perpendicular velocity [m/s]
		"""
		mass2D = np.broadcast_to(self.mass,(self.nsteps+1,self.nparticles))
		vperp  = np.sqrt(2.0*self.moment_lines*self.B_lines/mass2D)
		vperp  = np.where(self.B_lines < 0,0,vperp)
		return vperp


	def calcBaxis(self,kphi=0):
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
		dex = np.argwhere(S2D == np.min(S2D)).flatten()
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
			dex = np.argwhere(S2D == np.min(S2D)).flatten()
			r0  = self.raxis[dex[0]]
			z0  = self.zaxis[dex[1]]
			cont_gen = contour_generator(x=self.raxis,y=self.zaxis,z=S2D.T, line_type=LineType.Separate)
			lines = cont_gen.lines(1.0)
			lines = lines[0]
			x   = lines[:,0] - r0
			y   = lines[:,1] - z0
			Aminor = Aminor + np.mean(np.sqrt(x*x+y*y))
		return Aminor/self.nphi

	def calcRmajor(self):
		"""Calculates the major radius

		This routine calculates the major radius using the
		normalized toroidal flux label (S).
		
		Returns
		-------
		Rmajor : float
			Equilibrium major radius
		"""
		import numpy as np
		from contourpy import contour_generator, LineType
		Rmajor = 0
		for k in range(self.nphi): 
			S2D = np.squeeze(self.S_ARR[:,k,:])
			S2D = np.where(S2D>1.2,1.2,S2D)
			dex = np.argwhere(S2D == np.min(S2D)).flatten()
			r0  = self.raxis[dex[0]]
			z0  = self.zaxis[dex[1]]
			cont_gen = contour_generator(x=self.raxis,y=self.zaxis,z=S2D.T, line_type=LineType.Separate)
			lines = cont_gen.lines(1.0)
			lines = lines[0]
			Rmajor = Rmajor + np.mean(lines[:,0])
		return Rmajor/self.nphi

	def calcVolume(self,ns=None):
		"""Calculates the differential volume

		This routine calcualtes the differential volume and the volume
		as a function of normalized toroidal flux (S). 
		
		Parameters
		----------
		ns : int (optional)
			Number of radial gridpoints (default ns_prof1)

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
		from scipy.signal import savgol_filter

		if not ns:
			ns = self.ns_prof1
		dr = self.raxis[1] - self.raxis[0]
		dz = self.zaxis[1] - self.zaxis[0]
		dp = self.phiaxis[1] - self.phiaxis[0]
		nfp = np.round(2*np.pi/self.phiaxis[-1])

		area = dr*dz
		vol  = self.raxis*dp*area
		vol2d = np.broadcast_to(vol,(self.nz,self.nr)).T
		dV = np.zeros((self.nr,self.nphi,self.nz))
		for j in range(self.nphi-1):
			grid = np.squeeze(self.S_ARR[:,j,:])
			dV[:,j,:] = np.where(grid<=1,vol2d,0.0)

		edges = np.linspace(0.0,1.0,ns+1)
		ds = edges[1]-edges[0]
		s  = (edges[1:]+edges[0:-1])*0.5

		plasma_dvolds = np.zeros((ns))
		for i in range(ns):
			dtemp = np.where( self.S_ARR > edges[i], \
							 dV,0.0)
			dtemp = np.where( self.S_ARR <= edges[i+1], \
							 dtemp,0.0)
			plasma_dvolds[i] = np.sum(dtemp)

		plasma_dvolds = savgol_filter(plasma_dvolds,5,2)*nfp/ds
		plasma_vol    = np.cumsum(plasma_dvolds)*ds

		z = np.polyfit(s,plasma_vol,3)
		p = np.poly1d(z)
		pp = np.polyder(p)
		plasma_dvolds = pp(s)
		plasma_vol    = np.cumsum(plasma_dvolds)*ds
		
		return s, plasma_vol, plasma_dvolds

	def calcProfiles(self,ns=None):
		"""Calculates the kinetic plasma profiles used

		This routine claculates the kinetic profiles used in the
		BEAMS3D simulation based on the background grid quantities.
		
		Parameters
		----------
		ns : int (optional)
			Number of radial gridpoints (default ns_prof1)

		Returns
		-------
		S : ndarray
			Normalized toroidal flux array
		ne : ndarray
			Electron Density [m^-3]
		ni : ndarray
			Ion Density [m^-3]
		te : ndarray
			Electron Temperature [eV]
		ti : ndarray
			Ion Temperature [eV]
		zeff : ndarray
			Effective ion charge
		"""
		import numpy as np
		from scipy.interpolate import PchipInterpolator

		# Radial grid
		if not ns:
			ns = self.ns_prof1
		s  = np.linspace(0.0,1.0,ns)

		# Extract data
		[C, IA] = np.unique(self.S_ARR,return_index=True)
		NE_FLAT = self.NE.flatten()
		TE_FLAT = self.TE.flatten()
		TI_FLAT = self.TI.flatten()
		ZEFF_FLAT = self.ZEFF_ARR.flatten()
		ne_temp = NE_FLAT[IA]
		te_temp = TE_FLAT[IA]
		ti_temp = TI_FLAT[IA]
		zeff_temp = ZEFF_FLAT[IA]

		# Extract Ion density
		nni = np.size(self.NI,0)
		ni_temp = np.zeros((nni,len(IA)))
		for i in range(nni):
			NI_FLAT = self.NI[i,:,:,:].flatten()
			ni_temp[i,:] = NI_FLAT[IA]

		# Make Mirror
		C = np.concatenate((-C[::-1], C))
		ne_temp = np.concatenate((ne_temp[::-1], ne_temp))
		te_temp = np.concatenate((te_temp[::-1], te_temp))
		ti_temp = np.concatenate((ti_temp[::-1], ti_temp))
		zeff_temp = np.concatenate((zeff_temp[::-1], zeff_temp))
		ni_temp = np.concatenate((ni_temp[:,::-1],ni_temp),axis=1)


		# Fit
		[C, IA] = np.unique(C,return_index=True)
		ne_spl = PchipInterpolator(C,ne_temp[IA])
		te_spl = PchipInterpolator(C,te_temp[IA])
		ti_spl = PchipInterpolator(C,ti_temp[IA])
		zeff_spl = PchipInterpolator(C,zeff_temp[IA])
		ne = ne_spl(s)
		te = te_spl(s)
		ti = ti_spl(s)
		zeff = zeff_spl(s)
		ni = np.zeros((nni,ns))
		for i in range(nni):
			ni_spl = PchipInterpolator(C,ni_temp[i,IA])
			ni[i,:] = ni_spl(s)

		return s,ne,ni,te,ti,zeff

	def calcEr(self,ns=None):
		"""Calculates the radial electric field

		This routine calcualtes the radial electric field and returns
		the radial grid (s), electrostatic scalar potential (V), and 
		the radial derivative of the electrostatic scalar potential.
		Note that Er = dV/ds * ds/dr
					 = dV/ds * ds/drho * drho/dr
					 = dV/ds * 2.0 * rho * (1.0/Aminor)
			s = rho*rho -> ds/drho = 2 * rho
			rho = r/Aminor -> drho/dr = 1.0 / Aminor 
		
		Parameters
		----------
		ns : int (optional)
			Number of radial gridpoints (default ns_prof1)

		Returns
		-------
		S : ndarray
			Normalized toroidal flux array
		V : ndarray
			Electrostatic scalar potential [V]
		dVds : ndarray
			Radial derivative of the ES potential (dV/ds) [V]
		"""
		import numpy as np
		from scipy.interpolate import PchipInterpolator

		# Radial grid
		if not ns:
			ns = self.ns_prof1
		s  = np.linspace(0.0,1.0,ns)

		# Extract data
		[C, IA] = np.unique(self.S_ARR,return_index=True)
		POT_FLAT = self.POT_ARR.flatten()
		pot_temp = POT_FLAT[IA]

		# Make Mirror
		C = np.concatenate((-C[::-1], C))
		pot_temp = np.concatenate((pot_temp[::-1], pot_temp))

		# Fit
		[C, IA] = np.unique(C,return_index=True)
		p = PchipInterpolator(C,pot_temp[IA])
		V = p(s)
		x = np.concatenate((-s[:0:-1],s))
		f = np.concatenate((V[:0:-1],V))
		dVds = np.gradient(f,x)
		return s, V, dVds

	def calcMagaxis(self):
		"""Calculates the magetic axis location

		This routine calcualtes the magnetic axis location in each
		toroidal cut of the background grid. The magnetic axis is
		returned in an R array and Z array for each toridal cut.
		
		Returns
		-------
		R : ndarray
			Magnetic axis radial position (R) [m]
		Z : ndarray
			Magnetic axis vertical position (R) [m]
		"""
		import numpy as np
		raxis = []
		zaxis = []
		for k in range(self.nphi-1):
			S2D = np.squeeze(self.S_ARR[:,k,:])
			S2D = np.where(S2D>1.2,1.2,S2D)
			dex = np.argwhere(S2D == np.min(S2D)).flatten()
			raxis.append(self.raxis[dex[0]])
			zaxis.append(self.zaxis[dex[1]])
		raxis.append(raxis[0])
		zaxis.append(zaxis[0])
		return np.array(raxis),np.array(zaxis)


	def calcDepo(self,ns=None):
		"""Calculates the deposition profile

		This routine calcualtes the radial birth profile in
		units of [part/m^3/s] 
		
		Parameters
		----------
		ns : int (optional)
			Number of radial gridpoints (default ns_prof1)

		Returns
		-------
		rho : ndarray
			Normalized minor radius array
		birth : ndarray
			Birth Rate [part/m^3/s]
		"""
		import numpy as np
		from scipy.interpolate import PchipInterpolator

		# Setup rho on centered grid
		if not ns:
			ns = self.ns_prof1
		edges = np.linspace(0.0,1.0,ns+1)
		rho  = (edges[1:]+edges[0:-1])/2.0

		# Calc volume elements
		[s,_,dVds]=self.calcVolume(ns=ns)
		rho_s = np.sqrt(s)
		p = PchipInterpolator(rho_s,2.0*rho_s*dVds)
		dVdrho = p(rho)

		# Select start index
		dex_start = 0
		if self.lbeam:
			dex_start = 1

		# Calc births
		births    = np.zeros((self.nbeams,ns))
		rho_lines = np.sqrt(self.S_lines)
		for b in range(self.nbeams):
			#dexb = np.nonzero(self.Beam == (b+1))
			dexb = self.Beam == (b+1)
			rho_temp = rho_lines[dex_start,dexb]
			for i in range(ns):
				w_temp = self.Weight[dexb]
				w_temp = np.where(rho_temp >= edges[i],w_temp,0.0)
				w_temp = np.where(rho_temp < edges[i+1],w_temp,0.0)
				births[b,i] = np.sum(w_temp)
			births[b,:] = ns * births[b,:] / dVdrho

		return rho,births

	def calcLoss(self):
		"""Calculates the losses as a function of time

		This routine calcualtes the cumulative losses as a function of
		time for each beam line.  Losses are of particles not
		markers.

		Returns
		-------
		time : ndarray
			Time array [s]
		loss : ndarray
			Losses [particles]
		"""
		import numpy as np

		# Setup time arrays
		time_ns  = np.linspace(0,0.999E-6,1000)
		time_mus = np.linspace(1E-6,0.999E-3,1000)
		time_ms  = np.linspace(1E-3,1,1000)
		time     = np.append(time_ns,[time_mus,time_ms])

		# Calc losses
		loss = np.zeros((self.nbeams,len(time)))
		for b in range(self.nbeams):
			ldex1 = self.end_state==2
			ldex2 = self.Beam == (b+1)
			lostdex = ldex1 & ldex2
			t = self.t_end[lostdex]
			w = self.Weight[lostdex]
			[counts,_] = np.histogram(t,bins=time,weights=w)
			nlost  = np.cumsum(counts)
			loss[b,1:] = nlost

		return time,loss

	def calcJrad(self,ns=None):
		"""Calculates radial fast ion current

		This routine calculates the radial fast ion current from the
		initial birth profile and the steady-state fast ion density.

		
		Parameters
		----------
		ns : int (optional)
			Number of radial gridpoints (default ns_prof1)

		Returns
		-------
		rho : ndarray
			Normalized minor radius array
		jrad : ndarray
			Radial current density [A/m^2]
		"""
		import numpy as np
		# Setup rho on centered grid
		if not ns:
			ns = self.ns_prof1
		edges = np.linspace(0.0,1.0,ns+1)
		# Determine subset of particles for initial distribution
		tdex = 1
		if self.lbeam: tdex = 3
		# Get plasma volume / area
		[s,V,dVds] = self.calcVolume(ns=ns)
		# Geta Aminor
		Aminor = self.calcAminor()
		area = np.broadcast_to(2*V/(Aminor*np.sqrt(s)),(self.nbeams,ns))
		# Calculate the Thermalized and Lost particles
		Itherm = np.zeros(self.nbeams,ns)
		Ilost  = np.zeros(self.nbeams,ns)
		for k in range(self.nbeams):
			thermdex = np.argwhere((self.end_state==1 and self.Beam==k+1))
			lostdex  = np.argwhere((self.end_state==2 and self.Beam==k+1))
			for i in thermdex:
				slines = self.S_lines[:,i]
				s1     = slines[tdex]
				temp   = np.argwhere(self.R_lines[:,i]>0)
				j1     = temp[-1]
				s2     = slines[j1]
				j1     = max(sum(edges>s1),1)
				j2     = min(sum(edges>s2),self.ns_prof1)
				Itherm[k,j1:j2] = Itherm[k,j1:j2] + (self.Weight[i]*beam_data.charge[i]*sign(j2-j1))
			for i in lostdex:
				slines = self.S_lines[:,i]
				s1     = slines[tdex]
				j1     = max(sum(edges>s1),1) 
				Ilost[k,j1:] = Ilost[k,j1:] + (self.Weight[i]*beam_data.charge[i])
		# Calc losses
		return (Itherm+Ilost)/area

	def plotorbit(self,markers=None,plot3D=None):
		"""Plots traces of the orbits in 3D

		This routine plots traces of the particle orbits in 3D.

		Parameters
		----------
		markers : list (optional)
			List of marker indices to plot (default: all)
		plot3D : plot3D object (optional)
			Plotting object to render to.
		"""
		import numpy as np
		import vtk
		from libstell.plot3D import PLOT3D
		# Handle optionals
		if plot3D: 
			lplotnow=False
			plt = plot3D
		else:
			lplotnow = True
			plt = PLOT3D()
		# Handle markers
		if type(markers) == type(None):
			markers_in = np.linspace(0,self.nparticles-1,dtype=int)
		else:
			markers_in = markers
		# Plot markers
		for i in markers_in:
			j = np.argwhere(np.squeeze(self.R_lines[i,:])>0)
			k = j[-1][0]
			points_array = np.zeros((k,3))
			points_array[:,0] = self.X_lines[i,0:k]
			points_array[:,1] = self.Y_lines[i,0:k]
			points_array[:,2] = self.Z_lines[i,0:k]
			# Convert numpy array to VTK points
			points = vtk.vtkPoints()
			for point in points_array:
				points.InsertNextPoint(point)
			plt.add3Dline(points,linewidth=2)
		# In case it isn't set by user.
		plt.setBGcolor()
		# Render if requested
		if lplotnow: plt.render()

	def plot_heatflux(self,factor=1.0,colormap='hot',beams=None,load_type='heatflux',plot3D=None):
		"""Plots the BEAMS3D wall heat flux

		This makes a 3D plot of the first wall heat flux. The user
		can specify the type of quantity plotted using load type
			'heatflux'	: First wall heat flux
			'shine'		: Shinethrough flux
			'strikes' 	: Wall strikes

		Parameters
		----------
		factor : float (optional)
			Scaleing factor to apply (default: 1.0)
		colormap : str (optional)
			Color map to plot points (default: hot)
		beams : list (optional)
			List of beams to consider in plot (default: all)
		load_type : str (optional)
			Type of plot (default: heatflux)
		plot3D : plot3D object (optional)
			Plotting object to render to.
		"""
		import numpy as np
		from libstell.plot3D import PLOT3D
		# Handle optionals
		if plot3D: 
			lplotnow=False
			plt = plot3D
		else:
			lplotnow = True
			plt = PLOT3D()
		if type(beams) is type(None):
			beams_use = list(range(self.nbeams))
		else:
			beams_use = [x - 1 for x in beams] 
		# Which quantitity to plot
		if load_type == 'heatflux':
			val = np.sum(self.wall_load[:,beams_use],axis=1)
		elif load_type == 'shine':
			val = np.sum(self.wall_shine[:,beams_use],axis=1)
		elif load_type == 'strikes':
			val = np.sum(self.wall_strikes[:,beams_use],axis=1)
		else:
			print(f'ERROR: plot_heatflux load_type must be heatflux, shine, or strikes. load_type={load_type} ')
			return
		# Make points
		points,triangles = plt.facemeshTo3Dmesh(self.wall_vertex.T,self.wall_faces.T)
		scalar = plt.valuesToScalar(val*factor)
		# Add to Render
		plt.add3Dmesh(points,triangles,FaceScalars=scalar,opacity=1.0,color=colormap)
		# In case it isn't set by user.
		plt.setBGcolor()
		# Render if requested
		if lplotnow: plt.render()

# BEASM3D Input Class
class BEAMS3D_INPUT():
	"""Class for working with BEAMS3D INPUT data

	"""
	def __init__(self, parent=None):
		self.libStell = LIBSTELL()

	def read_input(self,filename):
		"""Reads BEAMS3D_INPUT namelist from a file

		This routine wrappers the beams3d_input_mod module reading routine.
		Parameters
		----------
		filename : string
			Input file name with BEASM3D_INPUT namelist
		"""
		indata_dict = self.libStell.read_beams3d_input(filename)
		for key in indata_dict:
			setattr(self, key, indata_dict[key])

	def write_input(self,filename):
		"""Writes BEASM3D_INPUT namelist to a file

		This routine wrappers the beams3d_input_mod module writing routine.
		Parameters
		----------
		filename : string
			Input file name to write BEASM3D_INPUT namelist to
		"""
		out_dict = vars(self)
		self.libStell.write_beams3d_input(filename,out_dict)

# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)



