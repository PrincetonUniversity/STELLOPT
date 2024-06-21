##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for reading and handling FOCUS
data.
"""

# Libraries
#from libstell.libstell import LIBSTELL, FourierRep

# Constants

# VMEC Class
class FOCUS():
	"""Class for working with VMEC equilibria

	"""
	def __init__(self):
		return

	def read_focus_HDF5(self,filename):
		"""Reads a FOCUS HDF5 file

		This routine reads and initilizes the FIELDLINES
		class with variable information from an HDF5 file.

		Parameters
		----------
		file : str
			Path to HDF5 file.
		"""
		import h5py
		import numpy as np
		focus_data = {}
		if '.h5' not in filename:
			filename = filename + '.h5'
		if 'focus_' not in filename:
			filename = 'focus_' + filename
		with h5py.File(filename,'r') as f:
			# Integers
			for temp in ['IsQuiet','IsSymmetric','case_surface',\
				'knotsurf','Nteta','Nzeta','case_init','case_coils',\
				'Ncoils','IsVaryCurrent','IsVaryGeometry','NFcoil',\
				'Nseg','case_optimize','IsNormalize','ISNormWeight',\
				'case_bnormal','case_length','case_curv','curv_alpha',\
				'DF_maxiter','CG_maxiter','HN_maxiter','TN_maxiter',\
				'TN_reorder','case_postproc','save_freq','save_coils',\
				'save_harmonics','save_filaments','update_plasma',\
				'pp_ns','pp_maxiter','Nfp','iout','ibnorm','mbnorm',\
				'ibharm','mbharm','itflux','mtflux','ittlen','mttlen',\
				'icssep','mcssep','icurv','mcurv']:
				if temp in f:
					setattr(self, temp, np.int64(f[temp][0]))
			# Reals
			for temp in ['knotsurf','init_current','init_radius',\
				'exit_xtol','weight_bnorm','weight_bharm','weight_tflux',\
				'target_tflux','weight_ttlen','target_length','k0',\
				'weight_specw','weight_cssep','weight_gnorm','weight_inorm',\
				'weight_mnorm','weight_curv','DF_tausta','DF_tauend','DF_xtol',\
				'CG_xtol','CG_wolfe_c1','CG_wolfe_c2','HN_xtol','HN_factor',\
				'TN_xtol','TN_cr','pp_phi','pp_raxis','pp_zaxis','pp_rmax',\
				'pp_zmax','pp_xtol','surf_vol','Inorm','Gnorm','Mnorm',\
				'overlap','time_initialize','time_optimize','time_postproc']:
				if temp in f:
					setattr(self, temp, float(f[temp][0]))
			# Arrays
			for temp in ['xsurf','ysurf','zsurf','nx','ny','nz','nn',\
				'plas_Bn','Bn','Bx','By','Bz','evolution','coilspace',\
				'deriv','Bmnin','Bmnim','initial_Bmnc','initial_Bmns',\
				'target_Bmnc','target_Bmns','Bmnc','Bmns','coil_importance',\
				'LM_fvec','LM_fjac','ppr','ppz','iota','XYZB','booz_mnc',\
				'booz_mns','bmim','bmin']:
				if temp in f:
					setattr(self, temp, np.array(f[temp][:]))

	def write_focus_plasma(self,nfp,xm,xn,rmnc,zmns,rmns=None,zmnc=None,xm_b=None,\
		xn_b=None,bmnc=None,bmns=None,filename='plasma.boundary'):
		"""Writes a focus boundary file

		This routine writes the FOCUS plasma boundary file.

		Parameters
		----------
		nfp : int
			Field periodicity
		xm : ndarray
			Poloidal mode array
		xn : ndarray
			Toroidal mode array
		rmnc : ndarray
			R cosine boundary harmonics
		zmns : ndarray
			Z sine boundary harmonics
		rmns : ndarray (optional)
			R sine boundary harmonics
		zmnc : ndarray (optional)
			Z cosine boundary harmonics
		xm_b : ndarray (optional)
			Poloidal mode array (B-normal)
		xn_b : ndarray (optional)
			Toroidal mode array (B-normal)
		bmnc : ndarray (optional)
			B-normal cosine boundary harmonics
		bmns : ndarray (optional)
			B-normal sine boundary harmonics
		filename : string (optional)
			Boundary file name (default: plasma.boundary)
		"""
		mnmax = len(xm)
		mnmax_b = 1
		if (xm_b and xn_b and (bmnc)):
			mnmax_b = len(xm_b)
		else:
			xm_b = [0]
			xn_b = [0]
			bmnc = [0]
			bmns = [0]
		if not (rmns and zmnc):
			rmns = np.zeros((mnmax))
			zmnc = np.zeros((mnmax))
		f=open(filename,'w')
		f.write(f'#Nfou Nfp NBnf\n')
		f.write(f'{int(mnmax)} {int(nfp)} {int(mnmax_b)}\n')
		f.write(f'#plasma boundary\n')
		f.write(f'# n m Rbc Rbs Zbc Zbs\n')
		for mn in range(mnmax):
			f.write(f'{int(xn[mn])} {int(xm[mn])} {rmnc[mn]:10.9f} {rmns[mn]:10.9f} {zmnc[mn]:10.9f} {zmns[mn]:10.9f}\n')
		f.write(f'#Bn harmonics\n')
		f.write(f'# n m bnc bns\n')
		for mn in range(mnmax_b):
			f.write(f'{int(xm_b[mn])} {int(xn_b[mn])} {bmnc[mn]:10.9f} {bmns[mn]:10.9f}\n')
		f.close()

	def plotConvergence(self,ax=None):
		"""Plots the FOCUS Convergence

		This routine plots the FOCUS code convergence for a given
		run.

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
		ax.semilogy(self.evolution[0],self.evolution[1],label=r'$\chi^2$')
		ax.semilogy(self.evolution[0],self.evolution[2],label=r"$|d \chi^2 / d {\bf X}|$")
		ax.semilogy(self.evolution[0],self.evolution[3],label=r"$f_{B_n}$")
		ax.semilogy(self.evolution[0],self.evolution[4],label=r"$f_{B_{mn}}$")
		ax.semilogy(self.evolution[0],self.evolution[5],label=r"$f_{\Psi}$")
		ax.semilogy(self.evolution[0],self.evolution[6],label=r"$f_L$")
		ax.semilogy(self.evolution[0],self.evolution[7],label=r"$f_{COIL-SURF}$")
		ax.semilogy(self.evolution[0],self.evolution[8],label=r"$f_{curv}$")
		ax.set_xlabel('Time [s]')
		ax.set_ylabel('Cost Function')
		ax.set_title('FOCUS Convergence')
		ax.legend()
		if lplotnow: pyplot.show()

	def plotBNormal(self,ax=None):
		"""Plots the FOCUS B-normal

		This routine plots the FOCUS code B-normal

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
		if self.IsSymmetric == 0:
			x = np.linspace(0,2*np.pi,int(self.Nzeta))
		elif self.IsSymmetric == 1:
			x = np.linspace(0,2*np.pi/self.Nfp,int(self.Nzeta))
		elif self.IsSymmetric == 2:
			x = np.linspace(0,np.pi/self.Nfp,int(self.Nzeta))
		else:
			return
		y = np.linspace(0,2.*np.pi,int(self.Nteta))
		hmesh=ax.pcolormesh(x,y,self.Bn.T,cmap='jet',shading='gouraud')
		#ax.semilogy(abscissa, data, **kwargs)
		ax.set_xlabel('Toroidal angle [rad]')
		ax.set_ylabel('Poloidal angle [rad]')
		ax.set_title('FOCUS B-Normal')
		pyplot.colorbar(hmesh,label='$B_n$ [T]',ax=ax)
		if lplotnow: pyplot.show()

	def vtkLUTHelper(self,ctable='jet'):
		"""Helper for VTK Lookup tables based on matplotlib

		The routine helps the user create a VTK lookup table based
		on the Matplotlib color tables.

		Parameters
		----------
		ctable : string
			Matplotlib colortable name
		Returns
		----------
		lut : VTKLookupTable
			VTK style lookup table class
		"""
		import vtk
		from matplotlib import cm
		lut = vtk.vtkLookupTable()
		lut.SetNumberOfTableValues(256)
		lut.Build()
		# Use matplotlib to generate the jet colormap
		jet = cm.get_cmap(ctable, 256)
		for i in range(256):
			rgba = jet(i / 255.0)
			lut.SetTableValue(i, rgba[0], rgba[1], rgba[2], rgba[3])
		return lut

	def plotBN3D(self,renderer=None,render_window=None):
		"""Plots the FOCUS B-normal in 3D

		This routine plots the FOCUS code B-normal in 3D

		Parameters
		----------
		renderer : vtkRenderer (optional)
			Renderer for plotting with VTK
		render_window : vtkRnderWindow (optional)
			Render window for plotting with VTK
		"""
		import numpy as np
		import matplotlib.pyplot as pyplot
		import vtk
		# Handle optionals
		lplotnow = True
		if renderer or render_window: lplotnow=False
		if not renderer: renderer = vtk.vtkRenderer()
		if not render_window: 
			render_window = vtk.vtkRenderWindow()
			render_window.AddRenderer(renderer)
			render_window_interactor = vtk.vtkRenderWindowInteractor()
			render_window_interactor.SetRenderWindow(render_window)
			render_window.SetSize(1024, 768)
		# Define VTK elements
		points = vtk.vtkPoints()
		triangles = vtk.vtkCellArray()
		polydata = vtk.vtkPolyData()
		actor = vtk.vtkActor()
		lut = self.vtkLUTHelper('jet')
		# Generate Mesh (v,u) indexing not (u,v)
		nu = self.Nteta # One short
		nv = self.Nzeta
		for v in range(nv):
			for u in range(nu):
				points.InsertNextPoint([self.xsurf[v,u],self.ysurf[v,u],self.zsurf[v,u]])
				if (v == nv -1): continue
				# Now do faces
				i1 = u + v * nu
				i2 = i1 + nu
				i3 = i2 + 1
				i4 = i1 + 1
				if u == nu-1:
					i3 = i1 + 1
					i4 = i1 - nu + 1
				triangle = vtk.vtkTriangle()
				triangle.GetPointIds().SetId(0, i1)
				triangle.GetPointIds().SetId(1, i2)
				triangle.GetPointIds().SetId(2, i4)
				triangles.InsertNextCell(triangle)
				triangle = vtk.vtkTriangle()
				triangle.GetPointIds().SetId(0, i2)
				triangle.GetPointIds().SetId(1, i3)
				triangle.GetPointIds().SetId(2, i4)
				triangles.InsertNextCell(triangle)
		polydata = vtk.vtkPolyData()
		polydata.SetPoints(points)
		polydata.SetPolys(triangles)
		vals = self.Bn.flatten()
		scalars = vtk.vtkFloatArray()
		scalars.SetNumberOfComponents(1)
		for value in vals:
			scalars.InsertNextValue(value)
		polydata.GetPointData().SetScalars(scalars)
		# Create a scalar bar (color bar) actor
		scalar_bar = vtk.vtkScalarBarActor()
		scalar_bar.SetLookupTable(lut)
		scalar_bar.SetTitle("")
		scalar_bar.SetNumberOfLabels(5)
		# Create a mapper and set the scalar range to the scalar values range
		mapper = vtk.vtkPolyDataMapper()
		mapper.SetInputData(polydata)
		mapper.SetLookupTable(lut)
		mapper.SetScalarRange(scalars.GetRange())
		actor.SetMapper(mapper)
		# Add actor to the scene
		renderer.AddActor(actor)
		# Render and interact
		renderer.AddActor2D(scalar_bar)
		renderer.SetBackground(0.1, 0.2, 0.3)
		if lplotnow:
			render_window.Render()
			render_window_interactor.Start()

	def plotPoincare(self,ax=None):
		"""Plots the FOCUS Poincare Plot

		This routine plots the FOCUS Poincare plots

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
		ax.plot(self.ppr,self.ppz,'.',markersize=1)
		rsurf = np.sqrt(self.xsurf**2+self.ysurf**2)
		rmin = min(rsurf[0,:])
		rmax = max(rsurf[0,:])
		zmin = min(self.zsurf[0,:])
		zmax = max(self.zsurf[0,:])
		ax.plot(rsurf[0,:],self.zsurf[0,:],'k')
		ax.set_xlabel('R [m]')
		ax.set_ylabel('Z [m]')
		ax.set_aspect('equal')
		ax.set_xlim(rmin*0.9,rmax*1.1)
		ax.set_ylim(zmin*1.1,zmax*1.1)
		if lplotnow: pyplot.show()

	def plotIota(self,ax=None):
		"""Plots the FOCUS Iota

		This routine plots the FOCUS rotational transform

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
		ax.plot(self.ppr[0]-self.ppr[0,1],self.iota)
		ax.set_xlabel('r [m]')
		ax.set_ylabel(r'$\iota$')
		ax.set_ylabel('Rotational Transform')
		if lplotnow: pyplot.show()

# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)

