#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for reading and handling NESCOIL
surface potential data.
"""

# Libraries
from libstell.libstell import LIBSTELL, FourierRep

# Constants

# NESCOIL Class
class NESCOIL(FourierRep):
	"""Class for working with NESCOIL output data.

	"""
	def __init__(self):
		super().__init__()
		self.libStell = LIBSTELL()

	def read_nescin(self,filename):
		"""Reads a NESCOIL input file.

		This routine reads and initilizes the NESCOIL
		class with information from a nescin file.

		Parameters
		----------
		file : str
			Path to nescin file.
		"""
		nescin_dict = self.libStell.read_nescoil_input(filename)
		for key in nescin_dict:
			setattr(self, key, nescin_dict[key])

	def write_nescin(self,filename):
		"""Writes a NESCOIL input file.

		This routine wrappers the write_nescoil_input module writing routine.
		Parameters
		----------
		filename : string
			Path to nescin file.
		"""
		out_dict = vars(self)
		self.libStell.write_nescoil_input(filename,out_dict)

	def read_nescout(self,filename):
		"""Reads a NESCOIL output file.

		This routine reads and initilizes the NESCOIL
		class with information from a nescout file.

		Parameters
		----------
		file : str
			Path to nescout file.
		"""
		nescout_dict = self.libStell.read_nescout(filename)
		for key in nescout_dict:
			setattr(self, key, nescout_dict[key])

	def bfield_init(self,nu=128,nv=128):
		"""Initilaizes the magnetic field from a NESCOIL surface current

		This routine initializes the magnetic field 
		from a NESCOIL surface current.

		Parameters
		----------
		nu : int (optional)
			Number of poloidal gridpoints (default: 128)
		nv : int (optional)
			Number of toroidal gridpoints (default: 128)
		"""
		return self.libStell.nescout_bfield_init(nu,nv)

	def bfield(self,x,y,z,istat=0):
		"""Evaluates the magnetic field from a NESCOIL surface current

		This routine evaluates the magnetic field at a point in space
		from a NESCOIL surface current. The status flag can be used to
		control the integration method.  If istat<0 then a discrete
		integral will be used, otherwise an adaptive integration is
		used.

		Parameters
		----------
		x : float
			X point to evaluate [m]
		y : float
			Y point to evaluate [m]
		z : float
			Z point to evaluate [m]
		istat : int (optional)
			Status flag
		Returns
		----------
		bx : float
			X-component of magnetic field [T]
		by : float
			Y-component of magnetic field [T]
		bZ : float
			Z-component of magnetic field [T]
		"""
		return self.libStell.nescout_bfield(x,y,z,istat=istat)

	def generatePotential(self,theta,zeta):
		"""Computes the potential on a grid

		This routine computes the potential on a grid

		Parameters
		----------
		theta : ndarray
			Poloidal Angle [rad]
		zeta : ndarray
			Field Period Angle [rad]
		Returns
		-------
		pot : ndarray
			NESCOIL potential []
		"""
		pot = self.sfunct(theta,zeta,self.potmns_surface.T,self.xm_pot,self.xn_pot)
		return pot

	def generateTotalPotential(self,theta,zeta):
		"""Computes the potential on a grid

		This routine computes the potential on a grid

		Parameters
		----------
		theta : ndarray
			Poloidal Angle [rad]
		zeta : ndarray
			Field Period Angle [rad]
		Returns
		-------
		pot : ndarray
			NESCOIL potential []
		"""
		import numpy as np
		pot = self.generatePotential(theta,zeta)
		nu = len(theta)
		nv = len(zeta)
		for j in range(nu): pot[0,j,:] = pot[0,j,:] - self.cut*0.5*theta[j]/np.pi
		for j in range(nv): 
			v = zeta[j]/(np.pi*2)
			pot[0,:,j] = pot[0,:,j] - self.cup*v
		return pot


	def plotpotential(self,ax=None):
		"""Plots the NESCOIL Potential

		This routine plots the NESCOIL code surface potential

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
		theta = np.ndarray((self.nu,1))
		zeta  = np.ndarray((self.nv,1))
		for j in range(self.nu): theta[j]=2.0*np.pi*j/float(self.nu-1)
		for j in range(self.nv):  zeta[j]=    np.pi*j/float(self.nv-1)
		pot = self.generatePotential(theta,zeta)
		#pot = self.sfunct(theta,zeta,self.potmns_surface.T,self.xm_pot,self.xn_pot)
		hmesh=ax.pcolormesh(np.squeeze(zeta),np.squeeze(theta),np.squeeze(pot[0,:,:]),cmap='jet',shading='gouraud')
		ax.set_xlabel('Toroidal angle [rad]')
		ax.set_ylabel('Poloidal angle [rad]')
		ax.set_title(r'NESCOIL $\Phi$ Potential')
		pyplot.colorbar(hmesh,label='$Pot$ [arb]',ax=ax)

	def plottotalpotential(self,ax=None):
		"""Plots the NESCOIL Total Potential

		This routine plots the NESCOIL code surface potential

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
		theta = np.ndarray((self.nu,1))
		zeta  = np.ndarray((self.nv,1))
		for j in range(self.nu): theta[j]=2.0*np.pi*j/float(self.nu-1)
		for j in range(self.nv):  zeta[j]=    np.pi*j/float(self.nv-1)
		pot = self.generateTotalPotential(theta,zeta)
		#pot = self.sfunct(theta,zeta,self.potmns_surface.T,self.xm_pot,self.xn_pot)
		#for j in range(self.nu): pot[0,j,:] = pot[0,j,:] - self.cut*1.0*j/float(self.nv-1)
		#for j in range(self.nv): pot[0,:,j] = pot[0,:,j] - self.cup*0.5*j/float(self.nv-1)
		hmesh=ax.pcolormesh(np.squeeze(zeta),np.squeeze(theta),np.squeeze(pot[0,:,:]),cmap='jet',shading='gouraud')
		#ax.semilogy(abscissa, data, **kwargs)
		ax.set_xlabel('Toroidal angle [rad]')
		ax.set_ylabel('Poloidal angle [rad]')
		ax.set_title(r'NESCOIL Total $\Phi$ Potential')
		pyplot.colorbar(hmesh,label='$Pot$ [arb]',ax=ax)
		if lplotnow: pyplot.show()

	def plotsurfaces(self,plot3D=None):
		"""Plots the NESCOIL Surfaces

		This routine plots the NESCOIL current potential surface
		and the plasma surface over a half field period.

		Parameters
		----------
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
		# Generate VTK objects
		theta = np.ndarray((self.nu,1))
		zeta  = np.ndarray((self.nv,1))
		for j in range(self.nu): theta[j]=2.0*np.pi*j/float(self.nu)
		for j in range(self.nv):  zeta[j]=np.pi*j/float(self.nv-1)
		rp = self.cfunct(theta,zeta,self.rmnc_plasma.T,self.xm_plasma,self.xn_plasma)
		zp = self.sfunct(theta,zeta,self.zmns_plasma.T,self.xm_plasma,self.xn_plasma)
		rc = self.cfunct(theta,zeta,self.rmnc_surface.T,self.xm_surface,self.xn_surface)
		zc = self.sfunct(theta,zeta,self.zmns_surface.T,self.xm_surface,self.xn_surface)
		self.isotoro(rp,zp,zeta/self.np,-1,plot3D=plt,lclosev=False,color='red')
		self.isotoro(rc,zc,zeta/self.np,-1,plot3D=plt,lclosev=False,color='green')
		# Render if requested
		if lplotnow: plt.render()

	def plotsurfaces_old(self,renderer=None,render_window=None):
		"""Plots the NESCOIL Surfaces

		This routine plots the NESCOIL current potential surface
		and the plasma surface over a half field period.

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
		theta = np.ndarray((self.nu,1))
		zeta  = np.ndarray((self.nv,1))
		for j in range(self.nu): theta[j]=2.0*np.pi*j/float(self.nu-1)
		for j in range(self.nv):  zeta[j]=np.pi*j/float(self.nv-1)
		rp = self.cfunct(theta,zeta,self.rmnc_plasma.T,self.xm_plasma,self.xn_plasma)
		zp = self.sfunct(theta,zeta,self.zmns_plasma.T,self.xm_plasma,self.xn_plasma)
		rc = self.cfunct(theta,zeta,self.rmnc_surface.T,self.xm_surface,self.xn_surface)
		zc = self.sfunct(theta,zeta,self.zmns_surface.T,self.xm_surface,self.xn_surface)
		self.isotoro(rp,zp,zeta/self.np,-1,renderer=renderer,render_window=render_window)
		self.isotoro(rc,zc,zeta/self.np,-1,renderer=renderer,render_window=render_window,color='green')
		if lplotnow:
			render_window.Render()
			render_window_interactor.Start()

	def cutcoils(self,ncoils_per_halfperiod,lplot=False):
		"""Cut coils from the NESCOIL potential

		This routine cuts coils from the NESCOIL potential.
		It allows the user to specify the number of coils per half 
		period.

		Parameters
		----------
		ncoils_per_halfperiod : integer
			Number of coils per half period (suggest 5)
		lplot : boolean (optional)
			Plot the potential and potential lines. (default: False)
		"""
		import numpy as np
		from libstell.coils import COILSET, COILGROUP, COIL
		from contourpy import contour_generator, LineType
		import matplotlib.pyplot as pyplot
		# Generate coilset
		coils = COILSET()
		coils.nfp = self.np
		coils.ngroups = ncoils_per_halfperiod
		coils.xmin = 1E9; coils.xmax=-1E9
		coils.ymin = 1E9; coils.ymax=-1E9
		coils.zmin = 1E9; coils.zmax=-1E9
		# First generate potential to determine contours
		theta = np.reshape( np.linspace(0,2*np.pi,self.nu),(self.nu,1))
		zeta  = np.reshape( np.linspace(0,np.pi,self.nv),(self.nv,1))
		pot = self.generateTotalPotential(theta,zeta)
		cont_vals = np.zeros((ncoils_per_halfperiod))
		for k in range(ncoils_per_halfperiod):
			u = 0
			v = round((k+0.5)*self.nv/(ncoils_per_halfperiod))
			cont_vals[k] = pot[0,u,v]
		# Now calculate a larger potential map so coils can span periods
		theta = np.reshape( np.linspace(0,2*np.pi,self.nu),(self.nu,1))
		zeta_min = (-1.0/ncoils_per_halfperiod)*np.pi
		zeta_max = (1.0+1.0/ncoils_per_halfperiod)*np.pi
		zeta  = np.reshape( np.linspace(zeta_min,zeta_max,self.nv),(self.nv,1))
		pot = self.generateTotalPotential(theta,zeta)
		# Now generate contours
		potmin = np.min(pot)
		potmax = np.max(pot)
		delta  = 2*(potmax-potmin)/float(2*ncoils_per_halfperiod+1.5)
		cont_gen = contour_generator(x=np.squeeze(zeta),y=np.squeeze(theta),z=np.squeeze(pot), line_type=LineType.Separate)
		# Make plot if requested
		if lplot:
			px = 1/pyplot.rcParams['figure.dpi']
			fig=pyplot.figure(figsize=(1024*px,768*px))
			ax=fig.add_subplot(111)
			hmesh=ax.contourf(np.squeeze(zeta),np.squeeze(theta),np.squeeze(pot),np.sort(cont_vals),extend='both',cmap='Greens')
			ax.contour(np.squeeze(zeta),np.squeeze(theta),np.squeeze(pot),np.sort(cont_vals),colors='black')
			ax.set_xlabel('Toroidal angle [rad]')
			ax.set_ylabel('Poloidal angle [rad]')
			ax.set_title(r'NESCOIL Coil Cutting')
			pyplot.colorbar(hmesh,label=r'Potential $\Phi$ [arb]',ax=ax)
			pyplot.show()
		# Now loop over contours
		for k in range(ncoils_per_halfperiod):
			level = cont_gen.lines(cont_vals[k])
			th = np.array([]); ze = np.array([])
			for temp in level:
				th = np.append(th,temp[:,1])
				ze = np.append(ze,temp[:,0])
			# Fourier transform the coil
			npts = len(th)
			r = np.zeros((npts)); z = np.zeros((npts))
			for mn in range(self.mnmax_surface):
				mtheta = th*self.xm_surface[mn]
				nzeta  = ze*self.xn_surface[mn]
				r  = r + np.cos(mtheta+nzeta)*self.rmnc_surface[mn]
				z  = z + np.sin(mtheta+nzeta)*self.zmns_surface[mn]
			# Convert to XYZ and make current/group
			ph = ze/float(self.np)
			x = r * np.cos(ph)
			y = r * np.sin(ph)
			c = np.ones((npts))*self.curpol/(np.pi*4E-7*ncoils_per_halfperiod*2)
			g = np.ones((npts))*(k+1)
			c[-1] = 0.0
			# Create stellarator symmetric coil
			ph = (2.0*np.pi - ze)/self.np
			xo = np.append(x,r[::-1]*np.cos(ph[::-1]))
			yo = np.append(y,r[::-1]*np.sin(ph[::-1]))
			zo = np.append(z,-z[::-1])
			co = np.append(c,c)
			go = np.append(g,g)
			x  = xo; y = yo; z = zo; c = co; g =go
			# Now make all field periods
			for mn in range(1,self.np):
				cop = np.cos(mn*self.alp)
				sip = np.sin(mn*self.alp)
				x = np.append(x,xo*cop - yo*sip)
				y = np.append(y,xo*sip + yo*cop)
				z = np.append(z,zo)
				c = np.append(c,co)
				g = np.append(g,go)
			coils.xmin = np.minimum(coils.xmin,np.min(x))
			coils.ymin = np.minimum(coils.ymin,np.min(y))
			coils.zmin = np.minimum(coils.zmin,np.min(z))
			coils.xmax = np.maximum(coils.xmax,np.max(x))
			coils.ymax = np.maximum(coils.ymax,np.max(y))
			coils.zmax = np.maximum(coils.zmax,np.max(z))
			# Now create group
			coil_name=f'MOD{k+1}'
			coils.groups.extend([COILGROUP(x,y,z,c,coil_name)])
		# Return a coil object
		return coils

# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)