#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for reading and handling 
FIELDLINES fieldline data data.
"""

# Libraries
from libstell.libstell import LIBSTELL

# Constants

# FIELDLINES Class
class FIELDLINES():
	"""Class for working with FIELDLINES data

	"""
	def __init__(self):
		test = 0

	def read_fieldlines(self,filename):
		"""Reads a FIELDLINES HDF5 file

		This routine reads and initilizes the FIELDLINES
		class with variable information from an HDF5 file.

		Parameters
		----------
		file : str
			Path to HDF5 file.
		"""
		import h5py
		import numpy as np
		fieldline_data = {}
		#read file
		with h5py.File(filename,'r') as f:
			# Logicals
			for temp in ['ladvanced', 'laxis_i', 'lcoil', 'lmgrid', \
				'lmu', 'lpies', 'lreverse', 'lspec', 'lvac', 'lvessel',\
				'lvmec']:
				if temp in f:
					setattr(self, temp, np.int64(f[temp][0])==1)
			# Integers
			for temp in ['nlines', 'nphi', 'npoinc', 'nr', 'nsteps', 'nz']:
				if temp in f:
					setattr(self, temp, np.int64(f[temp][0]))
            # Floats
			for temp in ['VERSION','iota0']:
				if temp in f:
					setattr(self, temp, float(f[temp][0]))
			# Arrays (1D)
			for temp in ['raxis', 'phiaxis', 'zaxis', \
						 'wall_strikes','L_lines']:
				if temp in f:
					setattr(self, temp, np.array(f[temp][:]))
			# Arrays (2D)
			for temp in ['wall_vertex', 'wall_faces', \
						 'R_lines', 'Z_lines', 'PHI_lines', \
						 'B_lines', 'Rhc_lines', 'Zhc_lines']:
				if temp in f:
					array = np.transpose(f[temp][:],(1,0))
					setattr(self, temp, np.array(array))
			# Arrays (3D)
			for temp in ['B_R', 'B_PHI', 'B_Z', \
						 'A_R', 'A_PHI', 'A_Z']:
				if temp in f:
					array = np.transpose(f[temp][:],(2,1,0))
					setattr(self, temp, np.array(array))
		# Make derived arrays
		self.nfp     = round(2*np.pi/self.phiaxis[-1])
		self.X_lines = self.R_lines*np.cos(self.PHI_lines)
		self.Y_lines = self.R_lines*np.sin(self.PHI_lines)
		# Fix B_R and B_Z
		for i in range(self.nr):
			self.B_R[i,:,:] = self.B_R[i,:,:]*self.B_PHI[i,:,:]/self.raxis[i]
			self.B_Z[i,:,:] = self.B_Z[i,:,:]*self.B_PHI[i,:,:]/self.raxis[i]
		if hasattr(self,'wall_faces'): self.wall_faces = self.wall_faces - 1

	def calc_reff(self):
		"""Calculates the effective radius

		Using the first traced fieldline the routine calculates
		the effective minor radius [m].

		Returns
		----------
		reff : ndarray
			Effective minor radius [m].
		"""
		import numpy as np
		x = self.R_lines
		y = self.Z_lines
		for i in range(self.nsteps):
			x[:,i] = x[:,i] - x[0,i]
			y[:,i] = y[:,i] - y[0,i]
		theta = np.arctan2(y,x)
		theta  = np.arctan2(y,x)
		dtheta = np.diff(theta,axis=1)
		dtheta = np.where(dtheta<-np.pi,dtheta+2*np.pi,dtheta)
		dtheta = np.where(dtheta>np.pi,dtheta-2*np.pi,dtheta)
		theta  = np.cumsum(dtheta,axis=1)
		reff   = np.mean(np.sqrt(x*x+y*y),axis=1)
		return reff

	def calc_iota(self):
		"""Calculates the rotational transform

		Using the first traced fieldline the routine calculates
		the effective minor radius [m], rotational transform,
		and error in rotational transform.

		Returns
		----------
		reff : ndarray
			Effective minor radius [m]
		iota : ndarray
			Rotational transform
		iota_err : ndarray
			Error in rotational transform
		"""
		import numpy as np
		x = self.R_lines
		y = self.Z_lines
		for i in range(self.nsteps):
			x[:,i] = x[:,i] - x[0,i]
			y[:,i] = y[:,i] - y[0,i]
		theta  = np.arctan2(y,x)
		dtheta = np.diff(theta,axis=1)
		dtheta = np.where(dtheta<-np.pi,dtheta+2*np.pi,dtheta)
		dtheta = np.where(dtheta>np.pi,dtheta-2*np.pi,dtheta)
		theta  = np.cumsum(dtheta,axis=1)
		reff       = np.mean(np.sqrt(x*x+y*y),axis=1)
		iota       = np.zeros((self.nlines))
		iota_err   = np.zeros((self.nlines))
		for i in range(self.nlines):
			p, residuals, rank, singular_values, rcond = np.polyfit(self.PHI_lines[i,0:self.nsteps-1],theta[i,:],1,full=True)
			iota[i] = p[0]
			iota_err[i] = np.sqrt(residuals)
		iota[0] = 2.0 * iota[1] - iota[2]
		return reff, iota, iota_err

	def plot_poincare(self,phi,nskip=1,ax=None,color_data=None):
		"""Creates a basic Poincare plot

		Plots a Poincare plot given a toroidal angle in radians (phi).
		A number of fieldlines to skip can also be provided (nskip).
		The user may also provide an axes (ax) to plot to.

		Returns
		----------
		phi : float
			Toroidal angle to plot. [radians]
		nskip : int (optional)
			Number of fieldlines to skip.
		ax : axes (optional)
			Matplotlib axes object to plot to.
		"""
		import numpy as np
		import matplotlib.pyplot as pyplot
		lplotnow = False
		lcdata = False
		if not ax:
			ax = pyplot.axes()
			lplotnow = True
		if color_data is not None:
			lcdata = True
		k = int(self.npoinc*phi/self.phiaxis[-1])
		rmin = np.amin(self.raxis)
		rmax = np.amax(self.raxis)
		x = self.R_lines[0:self.nlines:nskip,k:self.nsteps-1:self.npoinc]
		y = self.Z_lines[0:self.nlines:nskip,k:self.nsteps-1:self.npoinc]
		if lcdata:
			c = color_data[0:self.nlines:nskip,k:self.nsteps-1:self.npoinc]
			ax.scatter(x,y,s=0.1,c=c,marker='.')
		else:
			ax.plot(x,y,'.k',markersize=0.1)
		ax.set_xlabel('R [m]')
		ax.set_ylabel('Z [m]')
		ax.set_title(rf'FIELDLINES $\phi$ = {np.rad2deg(phi):3.1f}')
		ax.set_aspect('equal')
		ax.set_xlim(rmin,rmax)
		if lplotnow: pyplot.show()

	def plot_cloud(self,k,pointsize=0.01,color='red',plot3D=None):
		"""Plots the FIELDILNES Poin in 3D

		This routine plots the FIELDLINES poincare points in 3D

		Parameters
		----------
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
		vertices = []
		scalar   = []
		for i in range(self.nsteps):
			if (self.R_lines[k,i] > 0):
				vertices.append([self.X_lines[k,i],self.Y_lines[k,i],self.Z_lines[k,i]])
				scalar.append(self.B_lines[k,i])
		vertices = np.array(vertices)
		scalar = plt.valuesToScalar(np.array(scalar))
		points = plt.vertexToPoints(vertices)
		# Add to Render
		if float(scalar.GetValueRange()[1]) > 0:
			plt.add3Dpoints(points,scalars=scalar,pointsize=pointsize)
			# Colorbar
			plt.colorbar()
		else:
			plt.add3Dpoints(points,pointsize=pointsize)
		# In case it isn't set by user.
		plt.setBGcolor()
		# Render if requested
		if lplotnow: plt.render()

	def plot_poincare3D(self,k=0,pointsize=0.01,color='red',plot3D=None):
		"""Plots the FIELDILNES Poincare cross section in 3D

		This routine makes a 3D Poincare plot.

		Parameters
		----------
		k : int (optional)
			Poincare cross section to plot (default: 0)
		pointsize : float (optional)
			Size of points (default: 0.01)
		color : str (optional)
			Color to plot points (default: red)
		plot3D : plot3D object (optional)
			Plotting object to render to.
		"""
		import numpy as np
		#import vtk
		from libstell.plot3D import PLOT3D
		# Handle optionals
		if plot3D: 
			lplotnow=False
			plt = plot3D
		else:
			lplotnow = True
			plt = PLOT3D()
		vertices = []
		scalar   = []
		P = np.mod(self.PHI_lines,self.phiaxis[-1])
		X = self.R_lines * np.cos(P)
		Y = self.R_lines * np.sin(P)
		for i in range(k,self.nsteps,self.npoinc):
			for j in range(self.nlines):
				if (self.R_lines[j,i] > 0):
					vertices.append([X[j,i],Y[j,i],self.Z_lines[j,i]])
		vertices = np.array(vertices)
		scalar = plt.valuesToScalar(np.array(scalar))
		points = plt.vertexToPoints(vertices)
		plt.add3Dpoints(points,pointsize=pointsize,color=color)
		# In case it isn't set by user.
		plt.setBGcolor()
		# Render if requested
		if lplotnow: plt.render()

	def plot_heatflux(self,factor=1.0,colormap='hot',plot3D=None):
		"""Plots the BEAMS3D wall heat flux

		This makes a 3D plot of the first wall heat flux.

		Parameters
		----------
		factor : float (optional)
			Scaleing factor to apply (default: 1.0)
		colormap : str (optional)
			Color map to plot points (default: hot)
		plot3D : plot3D object (optional)
			Plotting object to render to.
		"""
		from libstell.plot3D import PLOT3D
		# Handle optionals
		if plot3D: 
			lplotnow=False
			plt = plot3D
		else:
			lplotnow = True
			plt = PLOT3D()
		# Make points

		points,triangles = plt.facemeshTo3Dmesh(self.wall_vertex,self.wall_faces)
		scalar = plt.valuesToScalar(self.wall_strikes*factor)
		# Add to Render
		plt.add3Dmesh(points,triangles,FaceScalars=scalar,opacity=1.0,color=colormap)
		# In case it isn't set by user.
		plt.setBGcolor()
		# Render if requested
		if lplotnow: plt.render()

# FIELDLINES Input Class
class FIELDLINES_INPUT():
	"""Class for working with FIELDLINES INPUT data

	"""
	def __init__(self, parent=None):
		self.libStell = LIBSTELL()

	def read_input(self,filename):
		"""Reads FIELDLINES_INPUT namelist from a file

		This routine wrappers the fieldlines_input_mod module reading routine.
		Parameters
		----------
		filename : string
			Input file name with FIELDLINES_INPUT namelist
		"""
		indata_dict = self.libStell.read_fieldlines_input(filename)
		for key in indata_dict:
			setattr(self, key, indata_dict[key])

	def write_input(self,filename):
		"""Writes FIELDLINES_INPUT namelist to a file

		This routine wrappers the fieldlines_input_mod module writing routine.
		Parameters
		----------
		filename : string
			Input file name to write FIELDLINES_INPUT namelist to
		"""
		out_dict = vars(self)
		self.libStell.write_fieldlines_input(filename,out_dict)

# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)
