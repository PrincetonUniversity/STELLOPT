##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for reading and handling 
FIELDLINES fieldline data data.
"""

# Libraries
#from libstell import libstell

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
			# Arrays
			for temp in ['phiaxis', 'raxis', 'zaxis', 'B_lines', \
				'PHI_lines', 'R_lines', 'Z_lines', 'B_PHI', 'B_R', \
				'B_Z','wall_vertex', 'wall_faces', 'wall_strikes', \
				'A_R', 'A_PHI', 'A_Z', 'L_lines', 'Rhc_lines',  \
				'Zhc_lines']:
				if temp in f:
					setattr(self, temp, np.array(f[temp][:]))
		# Make derived arrays
		self.X_lines = self.R_lines*np.cos(self.PHI_lines)
		self.Y_lines = self.R_lines*np.sin(self.PHI_lines)
		# Fix B_R and B_Z
		for i in range(self.nr):
			self.B_R[i,:,:] = self.B_R[i,:,:]*self.B_PHI[i,:,:]/self.raxis[i]
			self.B_Z[i,:,:] = self.B_Z[i,:,:]*self.B_PHI[i,:,:]/self.raxis[i]

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
			x[i,:] = x[i,:] - x[i,0]
			y[i,:] = y[i,:] - y[i,0]
		theta = np.arctan2(y,x)
		theta  = np.arctan2(y,x)
		dtheta = np.diff(theta,axis=0)
		dtheta = np.where(dtheta<-np.pi,dtheta+2*np.pi,dtheta)
		dtheta = np.where(dtheta>np.pi,dtheta-2*np.pi,dtheta)
		theta  = np.cumsum(dtheta,axis=0)
		reff   = np.mean(np.sqrt(x*x+y*y),axis=0)
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
			x[i,:] = x[i,:] - x[i,0]
			y[i,:] = y[i,:] - y[i,0]
		theta = np.arctan2(y,x)
		theta  = np.arctan2(y,x)
		dtheta = np.diff(theta,axis=0)
		dtheta = np.where(dtheta<-np.pi,dtheta+2*np.pi,dtheta)
		dtheta = np.where(dtheta>np.pi,dtheta-2*np.pi,dtheta)
		theta  = np.cumsum(dtheta,axis=0)
		reff       = np.mean(np.sqrt(x*x+y*y),axis=0)
		iota       = np.zeros((self.nlines))
		iota_err   = np.zeros((self.nlines))
		for i in range(self.nlines):
			p, residuals, rank, singular_values, rcond = np.polyfit(self.PHI_lines[0:nstep-1,i],theta[:,i],1,full=True)
			iota[i] = p[0]
			iota_err[i] = np.sqrt(residuals)
		iota[0] = 2.0 * iota[1] - iota[2]
		return reff, iota, iota_err

	def plot_poincare(self,phi,nskip=1,ax=None):
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
		if not ax:
			ax = pyplot.axes()
			lplotnow = True
		k = int(self.npoinc*phi/self.phiaxis[-1])
		rmin = np.amin(self.raxis)
		rmax = np.amax(self.raxis)
		x = self.R_lines[k:self.nsteps-1:self.npoinc,0:self.nlines:nskip]
		y = self.Z_lines[k:self.nsteps-1:self.npoinc,0:self.nlines:nskip]
		ax.plot(x,y,'.k',markersize=0.1)
		ax.set_xlabel('R [m]')
		ax.set_ylabel('Z [m]')
		ax.set_title(rf'FIELDLINES $\phi$ = {np.rad2deg(phi):3.1f}')
		ax.set_aspect('equal')
		ax.set_xlim(rmin,rmax)
		if lplotnow: pyplot.show()

# Main routine
if __name__=="__main__":
	import sys
	from argparse import ArgumentParser
	import matplotlib.pyplot as pyplot
	import vmec
	import numpy as np
	parser = ArgumentParser(description= 
		'''Provides class for accessing fieldlines data also servers as a
		   simple tool for assessing fieldlines output.''')
	parser.add_argument("-f", "--fieldlines", dest="fieldlines_ext",
		help="FIELDLINES file extension", default = None)
	parser.add_argument("-p", "--plot", dest="lplot", action='store_true',
		help="Plot the fieldlines file.", default = False)
	parser.add_argument("-v", "--vmec", dest="vmec_ext", 
		help="Add VMEC equilbrium to plot", default = None)
	args = parser.parse_args()
	field_data = FIELDLINES()
	if args.fieldlines_ext:
		field_data.read_fieldlines('fieldlines_'+args.fieldlines_ext+'.h5')
		if args.lplot:
			px = 1/pyplot.rcParams['figure.dpi']
			fig,(ax1,ax2,ax3) = pyplot.subplots(1,3,sharey=True,figsize=(1024*px,512*px))
			pyplot.subplots_adjust(hspace=0.1,wspace=0.15)
			k = 0
			field_data.plot_poincare(k,6,ax=ax1)
			k = field_data.phiaxis[int(field_data.nphi/4)]
			field_data.plot_poincare(k,6,ax=ax2)
			k = field_data.phiaxis[int(field_data.nphi/2)]
			field_data.plot_poincare(k,6,ax=ax3)
			if args.vmec_ext:
				vmec_wout = vmec.VMEC()
				vmec_wout.read_wout(args.vmec_ext)
				theta = np.ndarray((360,1))
				zeta  = np.ndarray((3,1))
				for j in range(360): theta[j]=2.0*np.pi*j/359.0
				for j in range(3):   zeta[j]=     np.pi*j/2.0
				r = vmec_wout.cfunct(theta,zeta,vmec_wout.rmnc,vmec_wout.xm,vmec_wout.xn/vmec_wout.nfp)
				z = vmec_wout.sfunct(theta,zeta,vmec_wout.zmns,vmec_wout.xm,vmec_wout.xn/vmec_wout.nfp)
				j = vmec_wout.ns-1
				ax1.plot(r[j,:,0],z[j,:,0],'r')
				ax2.plot(r[j,:,1],z[j,:,1],'r')
				ax3.plot(r[j,:,2],z[j,:,2],'r')
			pyplot.show()
	sys.exit(0)
