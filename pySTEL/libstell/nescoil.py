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
		import copy
		nescin_dict = copy.deepcopy(self.libStell.read_nescoil_input(filename))
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
		import copy
		nescout_dict = copy.deepcopy(self.libStell.read_nescout(filename))
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

		This routine computes the normalised potential on a grid

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

		This routine computes the normalised potential on a grid

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
		# note that while techincally this should be pot - u and pot - v,
		# NESCOIL says dpot/dv = dphi/dv + v in surfcur_diag....so we use that
		for j in range(nu): pot[0,j,:] = pot[0,j,:] + self.cut*0.5*theta[j]/np.pi
		for j in range(nv): pot[0,:,j] = pot[0,:,j] + self.cup*0.5*zeta[j]/np.pi
		return pot

	def plotpotential(self,ax=None,cmap='jet'):
		"""Plots the NESCOIL Potential

		This routine plots the NESCOIL code surface potential

		Parameters
		----------
		ax : axes (optional)
			Matplotlib axes object to plot to.
		cmap : string (optional)
			Matplotlib colormap.

		Returns
		-------
		quadmesh : matplotlib.collections.Quadmesh
			Quadmesh as produced by pcolormesh
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
		quadmesh=ax.pcolormesh(np.squeeze(zeta),np.squeeze(theta),np.squeeze(pot[0,:,:]),cmap=cmap,shading='gouraud')
		ax.set_xlabel('Toroidal angle [rad]')
		ax.set_ylabel('Poloidal angle [rad]')
		ax.set_title(r'NESCOIL $\Phi$ Potential')
		pyplot.colorbar(quadmesh,label='$Pot$ [arb]',ax=ax)
		if lplotnow: pyplot.show()
		return quadmesh

	def plottotalpotential_old(self,ax=None,cmap='jet'):
		"""Plots the NESCOIL Total Potential (old)

		This routine plots the NESCOIL code surface potential

		Parameters
		----------
		ax : axes (optional)
			Matplotlib axes object to plot to.

		Returns
		-------
		quadmesh : matplotlib.collections.Quadmesh
			Quadmesh as produced by pcolormesh
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
		quadmesh=ax.pcolormesh(np.squeeze(zeta),np.squeeze(theta),np.squeeze(pot[0,:,:]),cmap=cmap,shading='gouraud')
		ax.set_xlabel('Toroidal angle [rad]')
		ax.set_ylabel('Poloidal angle [rad]')
		ax.set_title(r'NESCOIL Total $\Phi$ Potential')
		pyplot.colorbar(quadmesh,label='$Pot$ [arb]',ax=ax)
		if lplotnow: pyplot.show()
		return quadmesh

	def plottotalpotential(self,ax=None,nlevels=5,cmap='jet'):
		"""Plots the NESCOIL Total Potential

		This routine plots the NESCOIL code surface potential

		Parameters
		----------
		ax : axes (optional)
			Matplotlib axes object to plot to.
		nlevels : int (optional)
			Number of contour levels (default: 5)
		cmap : string (optional)
			Colormap (default: jet)

		Returns
		-------
		quadmesh : matplotlib.collections.Quadmesh
			Quadmesh as produced by pcolormesh
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
		cont_vals = np.zeros((nlevels))
		for k in range(nlevels):
			u = round(0.0*self.nu)
			v = round((k+0.5)*self.nv/(nlevels))
			cont_vals[k] = pot[0,u,v]
		hmesh=ax.contourf(np.squeeze(zeta),np.squeeze(theta),np.squeeze(pot),np.sort(cont_vals),extend='both',cmap='Greens')
		ax.contour(np.squeeze(zeta),np.squeeze(theta),np.squeeze(pot),np.sort(cont_vals),colors='black')
		#quadmesh=ax.pcolormesh(np.squeeze(zeta),np.squeeze(theta),np.squeeze(pot[0,:,:]),cmap=cmap,shading='gouraud')
		ax.set_xlabel('Toroidal angle [rad]')
		ax.set_ylabel('Poloidal angle [rad]')
		ax.set_title(r'NESCOIL Total $\Phi$ Potential')
		pyplot.colorbar(hmesh,label='$Pot$ [arb]',ax=ax)
		if lplotnow: pyplot.show()
		return hmesh

	def computesurfaces(self):
		"""Mesh the NESCOIL Surfaces in real space over half field period."""
		import numpy as np
		self.theta = np.ndarray((self.nu,1))
		self.zeta  = np.ndarray((self.nv,1))
		for j in range(self.nu): self.theta[j]=2.0*np.pi*j/float(self.nu-1)
		for j in range(self.nv): self.zeta[j]=np.pi*j/float(self.nv-1)   ## this is the toroidal angle \varphi/nfp
		self.rp = self.cfunct(self.theta,self.zeta,self.rmnc_plasma.T,self.xm_plasma,self.xn_plasma)
		self.zp = self.sfunct(self.theta,self.zeta,self.zmns_plasma.T,self.xm_plasma,self.xn_plasma)
		self.rc = self.cfunct(self.theta,self.zeta,self.rmnc_surface.T,self.xm_surface,self.xn_surface)
		self.zc = self.sfunct(self.theta,self.zeta,self.zmns_surface.T,self.xm_surface,self.xn_surface)
		return self

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
		self.computesurfaces()
		self.isotoro(self.rp,self.zp,self.zeta/self.np,-1,plot3D=plt,lclosev=False,color='red')
		self.isotoro(self.rc,self.zc,self.zeta/self.np,-1,plot3D=plt,lclosev=False,color='green')
		# Render if requested
		if lplotnow: plt.render()

	def cutcoils_old(self,ncoils_per_halfperiod,npts=128,lplot=False):
		"""Cut coils from the NESCOIL potential

		This routine cuts coils from the NESCOIL potential.
		It allows the user to specify the number of coils per half 
		period.

		Parameters
		----------
		ncoils_per_halfperiod : integer
			Number of coils per half period (suggest 5)
		npts : int
			Number of points in coil (default: 128)
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
		# Compute total current
		Ipol = self.curpol*self.np/(4.0E-7*np.pi)
		# First generate potential to determine contours
		theta = np.reshape( np.linspace(0,2*np.pi,self.nu),(self.nu,1))
		zeta  = np.reshape( np.linspace(0,np.pi,self.nv),(self.nv,1))
		pot = self.generateTotalPotential(theta,zeta)
		cont_vals = np.zeros((ncoils_per_halfperiod))
		for k in range(ncoils_per_halfperiod):
			u = round(0.0*self.nu)
			v = round((k+0.5)*self.nv/(ncoils_per_halfperiod))
			cont_vals[k] = pot[0,u,v]
		# Now calculate a larger potential map so coils can span periods
		theta = np.reshape( np.linspace(0,2*np.pi,self.nu),(self.nu,1))
		zeta_min = (-2.0/ncoils_per_halfperiod)*np.pi
		zeta_max = (1.0+2.0/ncoils_per_halfperiod)*np.pi
		zeta  = np.reshape( np.linspace(zeta_min,zeta_max,self.nv),(self.nv,1))
		pot = self.generateTotalPotential(theta,zeta)
		# Now generate contours
		cont_gen = contour_generator(x=np.squeeze(zeta),y=np.squeeze(theta),z=np.squeeze(pot), line_type=LineType.Separate, chunk_size=0)
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
			# One contour per level
			for temp in level:
				th_t = temp[:,1]
				ze_t = temp[:,0]
				if th_t[0] == th_t[-1]:
					if th_t[0] > 0:
						th_t = th_t - np.pi*2.0
				th = np.append(th_t[0:-1],th)
				ze = np.append(ze_t[0:-1],ze)
			#print('========')
			#print(temp)
			#print(th)
			#print(ze)
			# Wrap the coil so that poitive current is positive field (counterclockwise from top)
			if (th[16]-th[0] > 0):
				th = th[::-1]
				ze = ze[::-1]
				print(rf'Flipping coil {k}')
			# Now we need to interpolate the coil onto the interval [0,2*pi] in theta.
			l_in   = np.linspace(0.0,1.0,len(th))
			l_out  = np.linspace(0.0,1.0,npts)
			th_out = np.interp(l_out,l_in,th)
			ph_out = np.interp(l_out,l_in,ze)/self.np
			# Fourier transform the coil
			r = np.zeros((npts)); z = np.zeros((npts))
			for mn in range(self.mnmax_surface):
				mtheta = th_out*self.xm_surface[mn]
				nzeta  = ph_out*self.xn_surface[mn]*self.np
				r  = r + np.cos(mtheta+nzeta)*self.rmnc_surface[mn]
				z  = z + np.sin(mtheta+nzeta)*self.zmns_surface[mn]
			# Convert to XYZ and make current/group
			x = r * np.cos(ph_out)
			y = r * np.sin(ph_out)
			c = np.ones((npts))*Ipol/(self.np*ncoils_per_halfperiod*2)
			g = np.ones((npts))*(k+1)
			c[-1] = 0.0
			# Create stellarator symmetric coil
			#phn = (2.0*np.pi/self.np - ph_out)
			phn = -ph_out
			xo = np.append(x,r[::-1]*np.cos(phn[::-1]))
			yo = np.append(y,r[::-1]*np.sin(phn[::-1]))
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


	def cutcoils(self,ncoils_per_halfperiod,npts=128,lplot=False):
		"""Cut coils from the NESCOIL potential

		This routine cuts coils from the NESCOIL potential.
		It allows the user to specify the number of coils per half 
		period.

		Parameters
		----------
		ncoils_per_halfperiod : integer
			Number of coils per half period (suggest 5)
		npts : int
			Number of points in coil (default: 128)
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
		# Compute total current
		Ipol = self.curpol*self.np/(4.0E-7*np.pi)
		# Calculate the potential map over full field period
		theta = np.linspace([0],[2.0*np.pi],self.nu)
		zeta = np.linspace([-np.pi],[np.pi],self.nv*2)
		pot = np.squeeze(self.generateTotalPotential(theta,zeta))
		# Recompute theta and zeta to match format
		theta = np.squeeze(theta)
		zeta = np.squeeze(zeta)
		# Now generate contours
		# Now loop over contours
		for k in range(ncoils_per_halfperiod):
			u = 0.0
			v = np.pi*(k+0.5)/ncoils_per_halfperiod
			#print(k,u,v)
			th,ze = self.trace_isocontour(np.squeeze(theta),np.squeeze(zeta),np.squeeze(pot), u, v, num_points=npts, period_x=True, period_y=True)
			#print(th)
			#print(ze)
			# Wrap the coil so that poitive current is positive field (counterclockwise from top)
			if (th[16]-th[0] > 0):
				th = th[::-1]
				ze = ze[::-1]
				print(rf'Flipping coil {k}')
			# Now we need to interpolate the coil onto the interval [0,2*pi] in theta.
			l_in   = np.linspace(0.0,1.0,len(th))
			l_out  = np.linspace(0.0,1.0,npts)
			th_out = np.interp(l_out,l_in,th)
			ph_out = np.interp(l_out,l_in,ze)/self.np
			# Fourier transform the coil
			r = np.zeros((npts)); z = np.zeros((npts))
			for mn in range(self.mnmax_surface):
				mtheta = th_out*self.xm_surface[mn]
				nzeta  = ph_out*self.xn_surface[mn]*self.np
				r  = r + np.cos(mtheta+nzeta)*self.rmnc_surface[mn]
				z  = z + np.sin(mtheta+nzeta)*self.zmns_surface[mn]
			# Convert to XYZ and make current/group
			x = r * np.cos(ph_out)
			y = r * np.sin(ph_out)
			c = np.ones((npts))*Ipol/(self.np*ncoils_per_halfperiod*2)
			g = np.ones((npts))*(k+1)
			c[-1] = 0.0
			# Create stellarator symmetric coil
			#phn = (2.0*np.pi/self.np - ph_out)
			phn = -ph_out
			xo = np.append(x,r[::-1]*np.cos(phn[::-1]))
			yo = np.append(y,r[::-1]*np.sin(phn[::-1]))
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

                # Make plot if requested
		if lplot:
			px = 1/pyplot.rcParams['figure.dpi']
			fig=pyplot.figure(figsize=(1024*px,768*px))
			ax=fig.add_subplot(111)
			hmesh=ax.contourf(np.squeeze(zeta),np.squeeze(theta),np.squeeze(pot),levels=2*ncoils_per_halfperiod+1,extend='both',cmap='Greens')
			ax.contour(np.squeeze(zeta),np.squeeze(theta),np.squeeze(pot),levels=2*ncoils_per_halfperiod+1,colors='black')
			ax.set_xlabel('Toroidal angle [rad]')
			ax.set_ylabel('Poloidal angle [rad]')
			ax.set_title(r'NESCOIL Coil Cutting')
			pyplot.colorbar(hmesh,label=r'Potential $\Phi$ [arb]',ax=ax)
			pyplot.show()
                        
		return coils

	def cutcoils_helical(self,nhelical_coils,lplot=False):
		"""Cut coils from the NESCOIL potential

		This routine cuts coils from the NESCOIL potential.
		It allows the user to specify the number of coils per half 
		period.

		Parameters
		----------
		nhelical_coils : integer
			Number of helical coils (suggest 2)
		lplot : boolean (optional)
			Plot the potential and potential lines. (default: False)
		"""
		import numpy as np
		from libstell.coils import COILSET, COILGROUP, COIL
		from contourpy import contour_generator, LineType
		import matplotlib.pyplot as pyplot
		# Generate coilset
		coils = COILSET()
		print('!!!!!!!!!!!!!!!!!!!!!!')
		print('!!  NOT IMPLEMENTED !!')
		print('!!!!!!!!!!!!!!!!!!!!!!')
		return coils


	def trace_isocontour(self, x, y, f, x0, y0, num_points=64, period_x=False, period_y=False):
		"""
		Traces an isocontour line from a given starting point on a 2D grid.
		
		Parameters:
			x (1D array): Grid coordinates along the first axis, shape (M,)
			y (1D array): Grid coordinates along the second axis, shape (N,)
			f (2D array): Evaluated functional values, shape (M, N)
			x0, y0 (float): Starting coordinate for the trace
			num_points (int): Exact number of points to return along the trajectory
			period_x (bool or float): Periodicity in x. If True, inferred from grid.
			period_y (bool or float): Periodicity in y. If True, inferred from grid.
			
		Returns:
			resampled_x (1D array): X coordinates of the trace, length `num_points`
			resampled_y (1D array): Y coordinates of the trace, length `num_points`
		"""
		import numpy as np
		from scipy.interpolate import RegularGridInterpolator
		from scipy.integrate import solve_ivp

		# 1. Parse Periodicity, Offsets, and Grid Spacings
		dx_grid = np.abs(x[1] - x[0])
		dy_grid = np.abs(y[1] - y[0])
		min_spacing = min(dx_grid, dy_grid)
		
		x_min, y_min = x[0], y[0]

		# PERIOD DEFINITION: 
		# Change to (x[-1] - x[0] + dx_grid) ONLY if your grid stops short of the repeating boundary.
		# If x[0]=-1 and x[-1]=1 are the exact same physical point, leave it as (x[-1] - x[0]).
		px = (x[-1] - x[0]) if period_x is True else (period_x if period_x else None)
		py = (y[-1] - y[0]) if period_y is True else (period_y if period_y else None)
		
		# Domain-shifted wrap function
		def wrap(val, p, offset):
			return offset + ((val - offset) % p) if p else val

		# 2. Setup Grid Interpolator
		interp = RegularGridInterpolator((x, y), f, method='linear', bounds_error=False, fill_value=None)
		
		# 3. Numerical Gradient Function (with proper domain shifts)
		def get_grad(pt):
			cx, cy = pt
			eps = 1e-5
			cx_p, cx_m = wrap(cx + eps, px, x_min), wrap(cx - eps, px, x_min)
			cy_p, cy_m = wrap(cy + eps, py, y_min), wrap(cy - eps, py, y_min)
			
			df_dx = (interp((cx_p, wrap(cy, py, y_min))) - interp((cx_m, wrap(cy, py, y_min)))) / (2 * eps)
			df_dy = (interp((wrap(cx, px, x_min), cy_p)) - interp((wrap(cx, px, x_min), cy_m))) / (2 * eps)
			return np.array([np.asarray(df_dx).item(), np.asarray(df_dy).item()])

		# 4. Define the ODE System
		def ode_func(s, pt):
			grad = get_grad(pt)
			norm = np.linalg.norm(grad)
			if norm < 1e-9:
				return np.array([0.0, 0.0])
			return np.array([-grad[1] / norm, grad[0] / norm])

		# 5. Scale-Independent Event Detection
		tol = 0.5 * min_spacing  
		lockout_distance = 4.0 * tol  

		def close_to_start(s, pt):
			if s < lockout_distance: 
				return 1.0  
			cx, cy = pt
			dx = cx - x0
			if px: dx = (dx + px/2) % px - px/2
			dy = cy - y0
			if py: dy = (dy + py/2) % py - py/2
			return np.sqrt(dx**2 + dy**2) - tol
		
		def out_of_bounds(s, pt):
			cx, cy = pt
			if not period_x and (cx < x[0] or cx > x[-1]): return -1.0
			if not period_y and (cy < y[0] or cy > y[-1]): return -1.0
			return 1.0

		close_to_start.terminal = True
		out_of_bounds.terminal = True

		max_len = 1.0 * ((x[-1] - x[0]) + (y[-1] - y[0]))
		max_step_val = 0.2 * tol

		# 6. Execute Integration
		sol_f = solve_ivp(
			ode_func, t_span=(0, max_len), y0=[x0, y0],
			events=[close_to_start, out_of_bounds], 
			rtol=1e-5, atol=1e-5, max_step=max_step_val
		)
		
		closed_loop = len(sol_f.t_events[0]) > 0
		
		if closed_loop:
			path_x, path_y = sol_f.y[0], sol_f.y[1]
		else:
			sol_b = solve_ivp(
				ode_func, t_span=(0, -max_len), y0=[x0, y0],
				events=[out_of_bounds], 
				rtol=1e-5, atol=1e-5, max_step=max_step_val
			)
			bx, by = sol_b.y[0][::-1], sol_b.y[1][::-1]
			path_x = np.concatenate((bx[:-1], sol_f.y[0]))
			path_y = np.concatenate((by[:-1], sol_f.y[1]))

		# 7. Uniform Resampling & SEAM-FREE Final Wrapping
		dx_pts, dy_pts = np.diff(path_x), np.diff(path_y)
		step_lens = np.sqrt(dx_pts**2 + dy_pts**2)
		arc_lengths = np.concatenate(([0], np.cumsum(step_lens)))
		
		if arc_lengths[-1] == 0:
			return np.full(num_points, x0), np.full(num_points, y0)
			
		target_arcs = np.linspace(0, arc_lengths[-1], num_points)
		resampled_x = np.interp(target_arcs, arc_lengths, path_x)
		resampled_y = np.interp(target_arcs, arc_lengths, path_y)
		
		# FIX: Wrap relative to the start point (x0, y0) instead of the grid minimum.
		# This keeps the curve visually continuous and eliminates boundary grazing spikes.
		if px: resampled_x = x0 + ((resampled_x - x0 + px/2) % px - px/2)
		if py: resampled_y = y0 + ((resampled_y - y0 + py/2) % py - py/2)
			
		return resampled_x, resampled_y

# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)
