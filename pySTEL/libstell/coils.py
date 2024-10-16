##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for working with coils files
"""

# Libraries

# Constants

# VMEC Class
class COILSET():
	"""Class for working with coils files

	"""
	def __init__(self):
		from collections import deque
		super().__init__()
		self.nfp = None
		self.ngroups = None
		self.groups = []
		self.xmin=None; self.xmax=None;
		self.ymin=None; self.ymax=None;
		self.zmin=None; self.zmax=None;
		self.color_cycle = deque(['g', 'b', 'c', 'm', 'y', 'k'])

	def read_coils_file(self,filename):
		"""Directly reads a coils file

		This routine reads a coils file into the class.

		Parameters
		----------
		filename : str
			Path to coils file.
		"""
		import numpy as np
		f = open(filename,'r')
		lines = f.readlines()
		f.close()
		if  'periods' in lines[0]:
			self.nfp = int(lines[0][8:])
		else:
			print("Bad Synatx line 1 in coils file")
		if 'begin filament' not in lines[1]:
			print("Bad Synatx line 2 in coils file")
		if 'mirror' not in lines[2]:
			print("Bad Synatx line 3 in coils file")
		npts = int(len(lines))-3
		if 'end' in lines[-1]: npts = npts - 1
		coords = np.zeros((3,npts))
		group = np.zeros((npts))
		current = np.zeros((npts))
		coilnames = []
		last_dex = 0
		self.ngroups = 0
		for i in range(0,npts):
			line = lines[i+3].split()
			if len(line) < 4: continue
			#print(line)
			coords[0,i] = float(line[0])
			coords[1,i] = float(line[1])
			coords[2,i] = float(line[2])
			current[i]  = float(line[3])
			if len(line) == 6:
				if int(line[4]) not in group:
					coilnames.extend([line[5]])
				self.ngroups = max([self.ngroups,int(line[4])])
				group[last_dex:i+1] = int(line[4])
				last_dex=i+1
		# Set extents
		self.xmin = min(coords[0]); self.xmax = max(coords[0])
		self.ymin = min(coords[1]); self.ymax = max(coords[1])
		self.zmin = min(coords[2]); self.zmax = max(coords[2])
		# Create the coil object
		for i in range(self.ngroups):
			x = coords[0,group==(i+1)]
			y = coords[1,group==(i+1)]
			z = coords[2,group==(i+1)]
			c = current[group==(i+1)]
			self.groups.extend([COILGROUP(x,y,z,c,coilnames[i])])

	def rescalecoils(self,npts_new):
		"""Changes coil resolution

		This routine changes the number of points defining the
		coils by splineing to a new resolution.

		Parameters
		----------
		npts_new : integer
			Number of new points to spline to.
		"""
		import numpy as np
		from scipy.interpolate import CubicSpline
		npts_new_array = np.full(self.ngroups,npts_new,dtype=int)
		for i in range(self.ngroups):
			s_new = np.linspace(0.0,1.0,npts_new_array[i])
			for j in range(self.groups[i].ncoils):
				x = self.groups[i].coils[j].x
				y = self.groups[i].coils[j].y
				z = self.groups[i].coils[j].z
				s = np.linspace(0.0,1.0,self.groups[i].coils[j].npts)
				cx = CubicSpline(s,x,bc_type='periodic')
				cy = CubicSpline(s,y,bc_type='periodic')
				cz = CubicSpline(s,z,bc_type='periodic')
				self.groups[i].coils[j].x = cx(s_new)
				self.groups[i].coils[j].y = cy(s_new)
				self.groups[i].coils[j].z = cz(s_new)
				self.groups[i].coils[j].npts = npts_new_array[i]

	def plotcoils(self,plot3D=None):
		"""Plots a coilset in 3D using VTK

		This routine plots coils in 3D using VTK

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
		# Setup color array
		color_txt=['red','green','blue','yellow','magenta','cyan','aqua']
		# Plot coils
		for i in range(self.ngroups):
			for j in range(self.groups[i].ncoils):
				points_array = np.zeros((self.groups[i].coils[j].npts,3))
				points_array[:,0] =self.groups[i].coils[j].x
				points_array[:,1] =self.groups[i].coils[j].y
				points_array[:,2] =self.groups[i].coils[j].z
				# Convert numpy array to VTK points
				points = vtk.vtkPoints()
				for point in points_array:
					points.InsertNextPoint(point)
				# Add to render
				plt.add3Dline(points,color=color_txt[i % len(color_txt)],linewidth=5)
		# In case it isn't set by user.
		plt.setBGcolor()
		# Render if requested
		if lplotnow: plt.render()

	def plotcoilsHalfFP(self,plot3D=None):
		"""Plots a half field period of a coilset in 3D using VTK

		This routine plots a half field period of a coilset in 3D using VTK

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
		# Setup color array
		color_txt=['red','green','blue','yellow','magenta','cyan','aqua']
		# Plot coils
		for i in range(self.ngroups):
			j=0
			points_array = np.zeros((self.groups[i].coils[j].npts,3))
			points_array[:,0] =self.groups[i].coils[j].x
			points_array[:,1] =self.groups[i].coils[j].y
			points_array[:,2] =self.groups[i].coils[j].z
			# Convert numpy array to VTK points
			points = vtk.vtkPoints()
			for point in points_array:
				points.InsertNextPoint(point)
			# Add to render
			plt.add3Dline(points,color=color_txt[i % len(color_txt)],linewidth=5)
		# In case it isn't set by user.
		plt.setBGcolor()
		# Render if requested
		if lplotnow: plt.render()

	def plotcoilplasmaDist(self,plot3D=None):
		"""Plots coil with coil-plasma distance

		This routine plots a half field period of a coilset in 3D using
		VTK. The coil is color coded by the coil-plasma distance.

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
		# Get the min and max values
		cmin = 1E20; cmax=-1E20;
		for i in range(self.ngroups):
			j = 0
			cmin = min(cmin,min(self.groups[i].coils[j].dist_surf)) 
			cmax = max(cmin,max(self.groups[i].coils[j].dist_surf)) 
		# Plot coils
		for i in range(self.ngroups):
			j=0
			points_array = np.zeros((self.groups[i].coils[j].npts,3))
			points_array[:,0] =self.groups[i].coils[j].x
			points_array[:,1] =self.groups[i].coils[j].y
			points_array[:,2] =self.groups[i].coils[j].z
			# Handle Color
			vals = np.array(self.groups[i].coils[j].dist_surf, dtype=float)
			scalar = plt.valuesToScalar(vals)
			# Convert numpy array to VTK points
			points = vtk.vtkPoints()
			for point in points_array:
				points.InsertNextPoint(point)
			# Add to render
			plt.add3Dline(points,scalars=scalar,linewidth=5)
		# Set color limits
		plt.setClim(cmin,cmax)
		# In case it isn't set by user.
		plt.setBGcolor()
		# Colorbar
		plt.colorbar(title='Distance (coil-plasma) [m]')
		# Render if requested
		if lplotnow: plt.render()

	def plotcoilcoilDist(self,plot3D=None):
		"""Plots coil with coil-coil distance

		This routine plots a half field period of a coilset in 3D using
		VTK. The coil is color coded by the minimum coil-coil distance.

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
		# Get the min and max values
		cmin = 1E20; cmax=-1E20;
		for i in range(self.ngroups):
			j = 0
			cmin = min(cmin,min(self.groups[i].coils[j].dist_coil)) 
			cmax = max(cmin,max(self.groups[i].coils[j].dist_coil)) 
		# Plot coils
		for i in range(self.ngroups):
			j=0
			points_array = np.zeros((self.groups[i].coils[j].npts,3))
			points_array[:,0] =self.groups[i].coils[j].x
			points_array[:,1] =self.groups[i].coils[j].y
			points_array[:,2] =self.groups[i].coils[j].z
			# Handle Color
			vals = np.array(self.groups[i].coils[j].dist_coil, dtype=float)
			scalar = plt.valuesToScalar(vals)
			# Convert numpy array to VTK points
			points = vtk.vtkPoints()
			for point in points_array:
				points.InsertNextPoint(point)
			# Add to render
			plt.add3Dline(points,scalars=scalar,linewidth=5)
		# Set color limits
		plt.setClim(cmin,cmax)
		# In case it isn't set by user.
		plt.setBGcolor()
		# Colorbar
		plt.colorbar(title='Distance (coil-coil) [m]')
		# Render if requested
		if lplotnow: plt.render()

	def plotcoilsRZ(self,*args,**kwargs):
		"""Plots each coil in the RZ plot

		This routine plots coils in 2D R,Z projection.
		"""
		import numpy as np
		import matplotlib.pyplot as pyplot
		c_temp = self.color_cycle[0]
		for i in range(self.ngroups):
			fig=kwargs.pop('fig',pyplot.figure())
			ax=kwargs.pop('axes',fig.add_subplot(111))
			rmax = -1E7; rmin = 1E7
			zmax = -1E7; zmin = 1E7
			for j in range(self.groups[i].ncoils):
				if j == 0:
					r = np.sqrt(self.groups[i].coils[j].x**2 + self.groups[i].coils[j].y**2)
					phi = np.arctan2(self.groups[i].coils[j].y,self.groups[i].coils[j].x)
					ax.plot(r,self.groups[i].coils[j].z, c=c_temp)
					rmin = min(rmin,min(r))
					rmax = max(rmin,max(r))
					zmin = min(zmin,min(self.groups[i].coils[j].z))
					zmax = max(zmin,max(self.groups[i].coils[j].z))
			ax.set_xlim(rmin*0.95,rmax*1.05); ax.set_xlabel('R [m]')
			ax.set_ylim(zmin*1.05,zmax*1.05); ax.set_ylabel('Z [m]')
			ax.set_title(f"Coil - {self.groups[i].name}")
			pyplot.show()
			self.color_cycle.rotate(1)
			c_temp = self.color_cycle[0]

	def write_coils_file(self,filename):
		"""Writes a coils file

		This routine writes a coils file into a file.

		Parameters
		----------
		filename : str
			Path to coils file.
		"""
		import numpy as np
		f = open(filename,'w')
		f.write(f"periods {self.nfp}\n")
		f.write(f"begin filament\n")
		f.write(f"mirror NIL\n")
		for i in range(self.ngroups):
			for j in range(self.groups[i].ncoils):
				current = np.ones((self.groups[i].coils[j].npts))*self.groups[i].current
				current[-1] = 0
				npts_write = self.groups[i].coils[j].npts
				for k in range(npts_write-1):
					f.write(f"{self.groups[i].coils[j].x[k]:.10E} {self.groups[i].coils[j].y[k]:.10E} {self.groups[i].coils[j].z[k]:.10E} {current[k]:.10E}\n")
				k = self.groups[i].coils[j].npts-1
				f.write(f"{self.groups[i].coils[j].x[k]:.10E} {self.groups[i].coils[j].y[k]:.10E} {self.groups[i].coils[j].z[k]:.10E} {current[k]:.10E} {i+1} {self.groups[i].name}\n")
		f.close()

	def coilbiot(self,x,y,z,extcur=None):
		"""Calculates field at point in space

		This routine calculates the magnetic field at a point in space
		given the point and external current array.

		Parameters
		----------
		x : real
			Cartesian x value [m].
		y : real
			Cartesian y value [m].
		z : real
			Cartesian z value [m].
		extcur : list
			Array of currents in coil groups [A]
		Returns
		----------
		bx : real
			Magnetic field in cartesian x direction [T]
		by : real
			Magnetic field in cartesian y direction [T]
		bz : real
			Magnetic field in cartesian z direction [T]
		"""
		import numpy as np
		bx = 0; by = 0; bz = 0
		for i in range(self.ngroups):
			for j in range(self.groups[i].ncoils):
				if extcur:
					bxt,byt,bzt = self.groups[i].coils[j].bfield(x,y,z,extcur[i])
				else:
					bxt,byt,bzt = self.groups[i].coils[j].bfield(x,y,z,self.groups[i].current)
				bx = bx + bxt
				by = by + byt
				bz = bz + bzt
		return bx,by,bz

	def coilvecpot(self,x,y,z,extcur=None):
		"""Calculates vector potential at point in space

		This routine calculates the vector potential at a point in space
		given the point and external current array.

		Parameters
		----------
		x : real
			Cartesian x value [m].
		y : real
			Cartesian y value [m].
		z : real
			Cartesian z value [m].
		extcur : list
			Array of currents in coil groups [A]
		Returns
		----------
		ax : real
			Vector potential in cartesian x direction []
		ay : real
			Vector potential in cartesian y direction []
		az : real
			Vector potential in cartesian z direction []
		"""
		import numpy as np
		ax = 0; ay = 0; az = 0
		for i in range(self.ngroups):
			for j in range(self.groups[i].ncoils):
				if extcur:
					axt,ayt,azt = self.groups[i].coils[j].vecpot(x,y,z,extcur[i])
				else:
					axt,ayt,azt = self.groups[i].coils[j].vecpot(x,y,z,self.groups[i].current)
				ax = ax + axt
				ay = ay + ayt
				az = az + azt
		return ax,ay,az

	def coiloffset(self,dist=0.0):
		"""Calculates offset from coils

		This routine calculates an offset position based on the
		geometric mean of the coil. It returns an odered list of points.

		Parameters
		----------
		dist : real
			Offset distance in [m] (- is toward geometric mean)
		Returns
		----------
		vertex : ndarray [3,ncoils,npts]
			Order set of points for each coil.
		"""
		import numpy as np
		from scipy import interpolate
		ntheta = 64
		l_new = np.linspace(0,1,ntheta)
		ncoils_total = 0
		for i in range(self.ngroups):
			for j in range(self.groups[i].ncoils):
				ncoils_total = ncoils_total + 1
		k = 0
		vertex = np.zeros((3,ncoils_total,ntheta-1))
		P_order = np.zeros(ncoils_total)
		for i in range(self.ngroups):
			for j in range(self.groups[i].ncoils):
				[x0,y0,z0] = self.groups[i].coils[j].geomCenter()
				R0 = np.sqrt(x0*x0+y0*y0)
				P0 = np.arctan2(y0,x0)
				P_order[k] = P0
				l = np.linspace(0,1,self.groups[i].coils[j].npts)
				x = np.interp(l_new,l,self.groups[i].coils[j].x,period=1)
				y = np.interp(l_new,l,self.groups[i].coils[j].y,period=1)
				z = np.interp(l_new,l,self.groups[i].coils[j].z,period=1)
				R = np.sqrt(x*x+y*y)
				p = np.arctan2(y,x)
				dr = R-R0
				dz = z-z0
				d  = np.sqrt(dr*dr+dz*dz)
				r2 = R + dr*dist/d
				z2 = z + dz*dist/d
				x2 = r2 * np.cos(p)
				y2 = r2 * np.sin(p)
				vertex[0][k][:] = x2[0:-1]
				vertex[1][k][:] = y2[0:-1]
				vertex[2][k][:] = z2[0:-1]
				k = k + 1
		# Reorder in phi
		Pdex = P_order.argsort()
		vertex = vertex[:,Pdex,:]
		P_order = P_order[Pdex]
		# Now Smooth
		#nphi = 180
		#phi_new = np.linspace(-np.pi,np.pi,nphi)
		#vertex2 = np.zeros((3,nphi,ntheta))
		#for i in range(ntheta):
		#	x_spl = interpolate.splrep(P_order, np.squeeze(vertex[0,:,i]),per=True)
		#	y_spl = interpolate.splrep(P_order, np.squeeze(vertex[1,:,i]),per=True)
		#	z_spl = interpolate.splrep(P_order, np.squeeze(vertex[2,:,i]),per=True)
		#	vertex2[0,:,i] = interpolate.splev(phi_new,x_spl)
		#	vertex2[1,:,i] = interpolate.splev(phi_new,y_spl)
		#	vertex2[2,:,i] = interpolate.splev(phi_new,z_spl)
		return vertex

	def coilSurfDist(self,xs,ys,zs):
		"""Calculates coil-surface distance

		This routine calculates the distance between a coil and a 
		surface defined by points in cartesian coordiantes (x,y,z).
		Values are stored in the coil atribute dist_coil.

		Parameters
		----------
		xs : ndarray
			X points defining surface [m]
		ys : ndarray
			Y points defining surface [m]
		zs : ndarray
			Z points defining surface [m]
		"""
		for i in range(self.ngroups):
			for j in range(self.groups[i].ncoils):
				self.groups[i].coils[j].surfDist(xs,ys,zs)

	def coilCoilDist(self):
		"""Calculates coil-coil distance

		This routine calculates the distance between a coil and
		another coil. Values are stored in the coil atribute coil_coil.
		"""
		for i in range(self.ngroups):
			for j in range(self.groups[i].ncoils):
				for k in range(self.ngroups):
					for l in range(self.groups[i].ncoils):
						if (i == k) and (j == l): continue
						xs = self.groups[k].coils[l].x
						ys = self.groups[k].coils[l].y
						zs = self.groups[k].coils[l].z
						self.groups[i].coils[j].coilDist(xs,ys,zs)

	def singleToMultiFilament(self,width=0.2,height=0.2, nwidth=4, nheight=4, lskip=None):
		"""Generates a multi-filament coil from a single filament coil set

		This routine utilizes the finiteBuildCoil routine to generate a
		multi-filament model for a coil. It returns a coil object. The
		coil width, height, and number of filaments can be specified in
		each dimmension.

		Parameters
		----------
		width : float
			Finite build coil width [m]
		height : float
			Finite build coil height [m]
		nwidth : int
			Number of coils in the width direction.
		nheight : int
			Number of coils in height direction.
		lskip : list
			List of coil groups to skip

		Returns
		----------
		mulit_coil : Object
			CoilGroup Object of the multi-filament coil
		"""
		import numpy as np
		# Check lskip
		if type(lskip) is type(None):
			lskip_array = [False] * self.ngroups
		else:
			lskip_array = lskip
		# Loop over coils
		vertices = []
		faces = []
		l = int(0)
		mcoil = COILSET()
		mcoil.nfp = self.nfp
		mcoil.xmin = self.xmin; mcoil.xmax = self.xmax
		mcoil.ymin = self.ymin; mcoil.ymax = self.ymax
		mcoil.zmin = self.zmin; mcoil.zmax = self.zmax
		mcoil.ngroups = self.ngroups
		for i in range(self.ngroups):
			if lskip_array[i]:
				mcoil.groups.extend(self.groups[i])
			else:
				ncoil_max = self.groups[i].ncoils
				xc = []; yc = []; zc = []; cc = []
				current = self.groups[i].current
				coilname = self.groups[i].name
				for j in range(ncoil_max):
					# Generate bounding boxes
					xx,yy,zz = self.groups[i].coils[j].finiteBuildCoil(width=width,height=height)
					# Get vectors
					Nwx = np.squeeze(xx[1,:] - xx[0,:])
					Nwy = np.squeeze(yy[1,:] - yy[0,:])
					Nwz = np.squeeze(zz[1,:] - zz[0,:])
					Nhx = np.squeeze(xx[2,:] - xx[0,:])
					Nhy = np.squeeze(yy[2,:] - yy[0,:])
					Nhz = np.squeeze(zz[2,:] - zz[0,:])
					x0  = np.squeeze(xx[0,:]); y0 = np.squeeze(yy[0,:]); z0 = np.squeeze(zz[0,:])
					c0  = np.ones((self.groups[i].coils[j].npts)) * current
					c0[-1] = 0.0
					for l in range(nwidth):
						for m in range(nheight):
							xx = x0 + l*Nwx/(nwidth-1) + m*Nhx/(nheight-1)
							yy = y0 + l*Nwy/(nwidth-1) + m*Nhy/(nheight-1)
							zz = z0 + l*Nwz/(nwidth-1) + m*Nhz/(nheight-1)
							xc.extend(xx); yc.extend(yy); zc.extend(zz); cc.extend(c0)
				mcoil.groups.extend([COILGROUP(np.array(xc),np.array(yc),np.array(zc),np.array(cc),coilname)])
		return mcoil

	def fitSurface(self):
		"""Returns best fit surface to coil set

		This routine returns the fourier harmonics which are a best fit
		to a point cloud generated from the coil set.

		Returns
		----------
		xm : ndarray
			Poloidal mode array
		xn : ndarray
			Toroidal mode array
		rmnc : ndarray
			R cos(m*th+n*phi) harmonics [m]
		zmns : ndarray
			Z sin(m*th+n*phi) harmonics [m]
		"""
		from libstell.plot3D import PLOT3D
		from libstell.libstell import FourierRep
		import numpy as np
		from scipy.optimize import minimize
		from pathlib import Path
		# Setup initial condition
		self.Nfeval = 0
		mmax = 1; nmax = 3
		nu = max(mmax*4,32); nv = max(nmax*2,32)
		# First extract points
		xt = []
		yt = []
		zt = []
		for i in range(self.ngroups):
			for j in range(self.groups[i].ncoils):
				xc = self.groups[i].coils[j].x
				yc = self.groups[i].coils[j].y
				zc = self.groups[i].coils[j].z
				xt.extend(xc)
				yt.extend(yc)
				zt.extend(zc)
		xt = np.array(xt)
		yt = np.array(yt)
		zt = np.array(zt)
		# Only one field period (and stellsym)
		pt = np.arctan2(yt,xt)
		xt = xt[pt>=0.0]
		yt = yt[pt>=0.0]
		zt = zt[pt>=0.0]
		pt = np.arctan2(yt,xt)
		xt = xt[pt<=np.pi/self.nfp]
		yt = yt[pt<=np.pi/self.nfp]
		zt = zt[pt<=np.pi/self.nfp]
		pt = np.arctan2(yt,xt)
		rt = np.sqrt(xt*xt+yt*yt)
		theta = np.ndarray((nu,1))
		phi   = np.ndarray((nv,1))
		theta = np.linspace([0],[np.pi*2],nu+1)
		phi   = np.linspace([0],[np.pi/self.nfp],nv)
		theta = theta[0:-1]
		xm = []; xn = []
		m=0; mn10 = 0; mn20 = 0; mn21 = 0; mn21m=0; mn01 = 0;
		for n in range(0,nmax+1):
			xm.append(m)
			xn.append(n)
		l = nmax
		for m in range(1,mmax+1):
			for n in range(-nmax,nmax+1):
				l = l + 1
				xm.append(m)
				xn.append(n)
				if (m==1 and n==0):
					mn10 = l
				if (m==0 and n==1):
					mn01 = l
				if (m==2 and n==0):
					mn20 = l
				if (m==2 and n==1):
					mn21 = l
				if (m==2 and n==-1):
					mn21m = l
		mnmax = len(xm)
		xm = np.array(xm,ndmin=2).T
		xn = np.array(xn,ndmin=2).T*self.nfp
		mn00 = 0
		#fit_factor = 1.0*(np.squeeze(xm)+abs(np.squeeze(xn)))
		#fit_factor = 1.0*(np.squeeze(xm)+1.0)
		fit_factor[fit_factor<1.0] = 1.0
		# Set options
		opts = (nu,nv,xm,xn,theta,phi,rt,zt,fit_factor,mn00)
		# Set initial condition
		rmnc_temp = np.zeros((1,mnmax))
		zmns_temp = np.zeros((1,mnmax))
		# Estimate from coils
		#rmnc_temp[0,:] = 1.0E-3
		#zmns_temp[0,:] = -1.0E-3
		rmnc_temp[0,mn00] = np.mean(rt)
		rmnc_temp[0,mn10] = (max(rt)-min(rt))*0.2
		#rmnc_temp[0,mn20] = rmnc_temp[0,mn10]*0.3
		#rmnc_temp[0,mn21] = rmnc_temp[0,mn20]*0.33
		#rmnc_temp[0,mn21m] = rmnc_temp[0,mn20]*0.33
		zmns_temp[0,mn10] = -max(zt)*1.0
		# Check for previous run
		my_file = Path("./coil_surf_harmonics.csv")
		if my_file.is_file():
			f = open("./coil_surf_harmonics.csv",'r')
			lines = f.readlines()
			f.close()
			for line in lines:
				if '#' in line: continue
				n_txt,m_txt,rbc_txt,rbs_txt,zbc_txt,zbs_txt = line.split(',')
				n = int(n_txt); m = int(m_txt)
				l = 0
				for l in range(mnmax):
					if (xn[l,0] == n) and (xm[l,0] == m):
							rmnc_temp[0,l] = float(rbc_txt)
							zmns_temp[0,l] = float(zbs_txt)
			for l in range(mnmax):
				print(f'{xn[l,0]} {xm[l,0]} {rmnc_temp[0,l]} {zmns_temp[0,l]}')
			if False:
				rmnc_temp = rmnc_temp*1.05
				zmns_temp = zmns_temp*0.95
		# Apply factor and remove Z00
		rmnc_temp = rmnc_temp*fit_factor
		zmns_temp = zmns_temp*fit_factor
		#rmnc_temp = np.sign(rmnc_temp[0,:])*np.power(abs(rmnc_temp[0,:]),fit_factor)
		#zmns_temp = np.sign(zmns_temp[0,:])*np.power(abs(zmns_temp[0,:]),fit_factor)
		zmns_temp = np.delete(zmns_temp,mn00)
		x0 = np.concatenate((np.squeeze(rmnc_temp),np.squeeze(zmns_temp)))
		#print(x0)
		# Call fit function
		if True:
			res = minimize(self.surf_fit_func, x0, method='CG',
					args=opts, options={'disp': True,'maxiter':100},
					callback=self.surf_fit_callbackF,tol=1.0E-1)
			xf = res.x
		else:
			xf = x0
		i1 = 0; i2 = i1+mnmax
		rmnc = np.broadcast_to(xf[i1:i2],(1,mnmax))
		rmnc = rmnc/fit_factor
		#rmnc = np.sign(rmnc)*np.power(abs(rmnc),1.0/fit_factor)
		i1 = i2; i2 = i1+mnmax-1
		xtemp = xf[i1:i2]
		#zmns = np.broadcast_to(np.insert(xtemp,mn00,0.0),(1,mnmax))**(1.0/fit_factor)
		zmns = np.broadcast_to(np.insert(xtemp,mn00,0.0),(1,mnmax))
		zmns = zmns/fit_factor
		#zmns = np.sign(zmns)*np.power(abs(zmns),1.0/fit_factor)
		#
		plt3d = PLOT3D()
		print('  Fitting Surface')
		theta = np.linspace([0],[np.pi*2],64)
		phi   = np.linspace([0],[np.pi*2],64)/self.nfp
		FR = FourierRep()
		r = FR.cfunct(theta,phi,rmnc,xm,xn)
		z = FR.sfunct(theta,phi,zmns,xm,xn)
		FR.isotoro(r,z,phi,0,plot3D=plt3d,lclosev=False)
		vertices = np.array([xt,yt,zt]).T
		points = plt3d.vertexToPoints(vertices)
		plt3d.add3Dpoints(points,pointsize=0.1)
		plt3d.render()
		return xm,xn,rmnc,zmns

	def surf_fit_func(self,x,*args):
		import numpy as np
		from libstell.libstell import FourierRep
		FR = FourierRep()
		nu = args[0]
		nv = args[1]
		xm = args[2]
		xn = args[3]
		theta = args[4]
		phi = args[5]
		boundr = args[6]
		boundz = args[7]
		fit_factor = args[8]
		mn00 = args[9]
		mnmax = xm.shape[0]
		i1 = 0; i2 = i1+mnmax
		#rmnc = np.broadcast_to(x[i1:i2],(1,mnmax))**(1.0/fit_factor)
		rmnc = np.broadcast_to(x[i1:i2],(1,mnmax))
		rmnc = rmnc/fit_factor
		#rmnc = np.sign(rmnc)*np.power(abs(rmnc),1.0/fit_factor)
		i1 = i2; i2 = i1+mnmax-1
		xtemp = x[i1:i2]
		#zmns = np.broadcast_to(np.insert(xtemp,mn00,0.0),(1,mnmax))**(1.0/fit_factor)
		zmns = np.broadcast_to(np.insert(xtemp,mn00,0.0),(1,mnmax))
		zmns = zmns/fit_factor
		#zmns = np.sign(zmns)*np.power(abs(zmns),1.0/fit_factor)
		r = FR.cfunct(theta,phi,rmnc,xm,xn)
		z = FR.sfunct(theta,phi,zmns,xm,xn)
		d = 0
		for u in range(nu):
			for v in range(nv):
				dr = boundr - np.squeeze(r[0,u,v])
				dz = boundz - np.squeeze(z[0,u,v])
				dl2 = dr*dr + dz*dz
				#dl2[dl2<0.001] = 100.0
				d = d + min(dl2)
		#print(d)
		return d
		
	def surf_fit_callbackF(self,intermediate_result):
		print(f'ITER: {self.Nfeval} -- dval: {intermediate_result.fun}')
		self.Nfeval += 1

	def blenderCoil(self,width=0.2,height=0.2,lfield_period=False):
		"""Generates the lists Blender needs to render a coilset

		This routine generates the verticies and faces lists which
		Blender needs to render a coil.

		Parameters
		----------
		width : float
			Finite build coil width [m]
		height : float
			Finite build coil height [m]
		lfield_period : boolean
			Return coilset over one field period (default: False)

		Returns
		----------
		vertices : list
			List of tuples defining verticies
		faces: list
			List of tubles defining faces
		"""
		import numpy as np
		# Generate volumetric coil rendering
		vertices = []
		faces = []
		l = int(0)
		for i in range(self.ngroups):
			ncoil_max = self.groups[i].ncoils
			if (lfield_period): ncoil_max = 1
			for j in range(ncoil_max):
				xx,yy,zz = self.groups[i].coils[j].finiteBuildCoil(width=width,height=height)
				for k in range(xx.shape[1]-1):
					vertices.append((xx[0,k],yy[0,k],zz[0,k]))
					vertices.append((xx[1,k],yy[1,k],zz[1,k]))
					vertices.append((xx[2,k],yy[2,k],zz[2,k]))
					vertices.append((xx[3,k],yy[3,k],zz[3,k]))
					faces.append((l,l+1,l+5))
					faces.append((l,l+5,l+4))
					l=l+1
					faces.append((l,l+1,l+5))
					faces.append((l,l+5,l+4))
					l=l+1
					faces.append((l,l+1,l+5))
					faces.append((l,l+5,l+4))
					l=l+1
					faces.append((l,l-3,l+1))
					faces.append((l,l+1,l+4))
					l=l+1
				vertices.append((xx[0, 0],yy[0, 0],zz[0, 0]))
				vertices.append((xx[1, 0],yy[1, 0],zz[1, 0]))
				vertices.append((xx[2, 0],yy[2, 0],zz[2, 0]))
				vertices.append((xx[3, 0],yy[3, 0],zz[3, 0]))
				l = l + 4
		return vertices,faces

class COILGROUP():
	"""Class which defines a coil group

	"""
	def __init__(self, x, y, z, current,name):
		self.name = name
		self.current = current[0]
		self.coils = []
		idex = [index for index, value in enumerate(current) if value == 0]
		self.ncoils = len(idex)
		i = 0
		for j in idex:
			# Reverse order the coil if current sign changes
			if current[j-1] != current[0]:
				x2 = x[i:j+1]
				y2 = y[i:j+1]
				z2 = z[i:j+1]
				self.coils.extend([COIL(x2[::-1],y2[::-1],z2[::-1])])
			else:
				self.coils.extend([COIL(x[i:j+1],y[i:j+1],z[i:j+1])])
			i = j+1

class COIL():
	"""Class which defines a coil

	"""
	def __init__(self,x,y,z):
		self.npts = len(x)
		self.x = x
		self.y = y
		self.z = z
		self.dx = x[1:]-x[0:-1]
		self.dy = y[1:]-y[0:-1]
		self.dz = z[1:]-z[0:-1]
		self.vx = self.y[0:-1]*self.dz - self.z[0:-1]*self.dy
		self.vy = self.z[0:-1]*self.dx - self.x[0:-1]*self.dz
		self.vz = self.x[0:-1]*self.dy - self.y[0:-1]*self.dx
		self.xt = None
		self.yt = None
		self.zt = None
		self.dist_surf = None
		self.dist_coil = None

	def vecpot(self,x,y,z,current):
		"""Calculates Vector potential

		This routine calculates the vector potential for a given coil
		a position in space and a current in said coil.

		Parameters
		----------
		x : real
			Cartesian x value [m].
		y : real
			Cartesian y value [m].
		z : real
			Cartesian z value [m].
		current : real
			Current in coil [A]
		Returns
		----------
		ax : real
			Vector potential in cartesian x direction [A/m]
		ay : real
			Vector potential in cartesian y direction [A/m]
		az : real
			Vector potential in cartesian z direction [A/m]
		"""
		import numpy as np
		x1 = x - self.x
		y1 = y - self.y
		z1 = z - self.z
		rw = np.sqrt(x1*x1+y1*y1+z1*z1)
		fa = ( rw[1:] + rw[0:-1] ) / \
			 ( rw[1:] * rw[0:-1]   * \
			 ( rw[1:] * rw[0:-1] + x1[1:] * x1[0:-1] + \
			   y1[1:] * y1[0:-1] + z1[1:] * z1[0:-1] ) )
		ax = sum( fa * self.dx ) * current
		ay = sum( fa * self.dy ) * current
		az = sum( fa * self.dz ) * current
		return ax, ay, az

	def bfield(self,x,y,z,current):
		"""Calculates magnetic field

		This routine calculates the magnetic field for a given coil
		a position in space and a current in said coil.

		Parameters
		----------
		x : real
			Cartesian x value [m].
		y : real
			Cartesian y value [m].
		z : real
			Cartesian z value [m].
		current : real
			Current in coil [A]
		Returns
		----------
		bx : real
			Magnetic field in cartesian x direction [T]
		by : real
			Magnetic field in cartesian y direction [T]
		bz : real
			Magnetic field in cartesian z direction [T]
		"""
		import numpy as np
		fac = 1.0E-7
		x1 = x - self.x
		y1 = y - self.y
		z1 = z - self.z
		rw = np.sqrt(x1*x1+y1*y1+z1*z1)
		fa = ( rw[1:] + rw[0:-1] ) / \
			 ( rw[1:] * rw[0:-1]   * \
			 ( rw[1:] * rw[0:-1] + x1[1:] * x1[0:-1] + \
			   y1[1:] * y1[0:-1] + z1[1:] * z1[0:-1] ) )
		ax = sum(fa*self.dx)*current
		ay = sum(fa*self.dy)*current
		az = sum(fa*self.dz)*current
		bx = sum( fa * self.vx * current ) - y * az + z * ay
		by = sum( fa * self.vy * current ) - z * ax + x * az
		bz = sum( fa * self.vz * current ) - x * ay + y * ax
		return fac*bx, fac*by, fac*bz

	def geomCenter(self):
		"""Calculates geometric center of the coil

		This routine calculates the geometric center of the coil.

		Parameters
		----------
		Returns
		----------
		x : real
			Center in x coordinate [m]
		y : real
			Center in y coordinate [m]
		z : real
			Center in z coordinate [m]
		"""
		import numpy as np
		return np.mean(self.x),np.mean(self.y),np.mean(self.z)

	def surfDist(self,xs,ys,zs):
		"""Calculates coil-surface distance

		This routine calculates the distance between a coil and a 
		surface defined by points in cartesian coordiantes (x,y,z).

		Parameters
		----------
		xs : ndarray
			X points defining surface [m]
		ys : ndarray
			Y points defining surface [m]
		zs : ndarray
			Z points defining surface [m]
		"""
		import numpy as np
		nsurf = len(xs)
		xc = np.broadcast_to(self.x,(nsurf,self.npts)).T
		yc = np.broadcast_to(self.y,(nsurf,self.npts)).T
		zc = np.broadcast_to(self.z,(nsurf,self.npts)).T
		x  = np.broadcast_to(xs,(self.npts,nsurf))
		y  = np.broadcast_to(ys,(self.npts,nsurf))
		z  = np.broadcast_to(zs,(self.npts,nsurf))
		dx = xc - x
		dy = yc - y
		dz = zc - z
		dl2 = dx*dx + dy * dy + dz * dz
		self.dist_surf = np.sqrt(np.min(dl2,axis=1))
		return

	def coilDist(self,xs,ys,zs):
		"""Calculates coil-coil distance

		This routine calculates the distance between a coil and
		another coil defined by points in cartesian coordiantes 
		(x,y,z).

		Parameters
		----------
		xs : ndarray
			X points defining the other coil [m]
		ys : ndarray
			Y points defining the other coil [m]
		zs : ndarray
			Z points defining the other coil [m]
		"""
		import numpy as np
		nsurf = len(xs)
		xc = np.broadcast_to(self.x,(nsurf,self.npts)).T
		yc = np.broadcast_to(self.y,(nsurf,self.npts)).T
		zc = np.broadcast_to(self.z,(nsurf,self.npts)).T
		x  = np.broadcast_to(xs,(self.npts,nsurf))
		y  = np.broadcast_to(ys,(self.npts,nsurf))
		z  = np.broadcast_to(zs,(self.npts,nsurf))
		dx = xc - x
		dy = yc - y
		dz = zc - z
		dl2 = dx*dx + dy * dy + dz * dz
		if type(self.dist_coil) is type(None):
			self.dist_coil = np.sqrt(np.min(dl2,axis=1))
		else:
			dist_coil = np.sqrt(np.min(dl2,axis=1))
			self.dist_coil = np.minimum(self.dist_coil,dist_coil)
		return

	def spline_tangent(self, order=3, der=1):
		"""Calculate the tangent of coil using spline interpolation

		This routine calculates the tangent of a coil using spline
		interpolation. Order and derivative level can be set by user.

		Parameters
		----------
		order : int (optional)
			Order of spline (default: 3)
		der : int
			Derivative order (default:1)
		"""
		import numpy as np
		from scipy import interpolate
		t = np.linspace(0, 1.0, len(self.x), endpoint=True)
		fx = interpolate.splrep(t, self.x, s=0, k=order, per=1)
		fy = interpolate.splrep(t, self.y, s=0, k=order, per=1)
		fz = interpolate.splrep(t, self.z, s=0, k=order, per=1)
		self.xt = interpolate.splev(t, fx, der=1)
		self.yt = interpolate.splev(t, fy, der=1)
		self.zt = interpolate.splev(t, fz, der=1)
		if der == 2:
			self.xa = interpolate.splev(t, fx, der=2)
			self.ya = interpolate.splev(t, fy, der=2)
			self.za = interpolate.splev(t, fz, der=2)
		return

	def finiteBuildCoil(self, width=0.1, height=0.1, frame="centroid", **kwargs):
		"""Expand single coil filament to a finite-build coil.

		This routine expands a coil filament of a finite-build coil.
		The coil width, height and build frame can be set by the user.

		Parameters
		----------
		width : float (optional)
			Toroidal width of the coil [m] (default: 0.1)
		height : float (optional)
			Radial height of the coil [m] (default: 0.1)
		frame : string (optional)
			Finite build frame "centroid", "frenet", "parallel" (default:centroid)

		Returns
		----------
		x : ndarry
			X-coordiante for plotting as a mesh [m]
		y : ndarry
			Y-coordiante for plotting as a mesh [m]
		z : ndarry
			Z-coordiante for plotting as a mesh [m]
		"""
		import numpy as np
		n = self.npts
		# calculate the tangent
		if self.xt is None:
			self.spline_tangent(order=1)
		xt = self.xt
		yt = self.yt
		zt = self.zt
		tt = np.sqrt(xt * xt + yt * yt + zt * zt)
		xt = xt / tt
		yt = yt / tt
		zt = zt / tt

		# use surface normal if needed
		if frame == "centroid":
			# use the geometry center is a good idea
			[center_x,center_y,center_z]=self.geomCenter()
			xn = self.x - center_x
			yn = self.y - center_y
			zn = self.z - center_z
			nt = xn * xt + yn * yt + zn * zt
			xn = xn - nt * xt
			yn = yn - nt * yt
			zn = zn - nt * zt
		if frame == "centroid_cyl":
			# use the geometry center but define in cylindrical coords
			[center_x,center_y,center_z]=self.geomCenter()
			phi = np.arctan2(self.y,self.x)
			r  = np.sqrt(self.x**2+self.y**2)
			center_r = np.sqrt(center_x**2+center_y**2)
			rn = r - center_r
			xn = rn*np.cos(phi)
			yn = rn*np.sin(phi)
			zn = self.z - center_z
			nt = xn * xt + yn * yt + zn * zt
			xn = xn - nt * xt
			yn = yn - nt * yt
			zn = zn - nt * zt
		elif frame == "frenet":
			self.spline_tangent(der=2)
			xn = self.xa
			yn = self.ya
			zn = self.za
		elif frame == "parallel":
			# parallel transport frame
			# Hanson & Ma, Parallel Transp ort Approach to Curve Framing, 1995
			def rotate(x, ang):
				c = np.cos(ang)
				s = np.sin(ang)
				return [
					[
						c + x[0] ** 2 * (1 - c),
						x[0] * x[1] * (1 - c) - s * x[2],
						x[2] * x[0] * (1 - c) + s * x[1],
					],
					[
						x[0] * x[1] * (1 - c) + s * x[2],
						c + x[1] ** 2 * (1 - c),
						x[2] * x[1] * (1 - c) - s * x[0],
					],
					[
						x[0] * x[2] * (1 - c) - s * x[1],
						x[1] * x[2] * (1 - c) + s * x[0],
						c + x[2] ** 2 * (1 - c),
					],
				]

			T = np.transpose([self.xt, self.yt, self.zt])
			T = T / np.linalg.norm(T, axis=1)[:, np.newaxis]
			B = np.cross(T[:-1], T[1:], axis=1)
			B = B / np.linalg.norm(B, axis=1)[:, np.newaxis]
			theta = np.arccos(np.sum(T[:-1] * T[1:], axis=1))
			V = np.zeros_like(T)
			kwargs.setdefault("vx", self.x[0] - np.average(self.x[0:-1]))
			kwargs.setdefault("vy", self.y[0] - np.average(self.y[0:-1]))
			vx = kwargs["vx"]
			vy = kwargs["vy"]
			vz = -(vx * T[0, 0] + vy * T[0, 1]) / T[0, 2]
			vv = np.linalg.norm([vx, vy, vz])
			V[0, :] = [vx / vv, vy / vv, vz / vv]
			print(np.dot(V[0, :], T[0, :]))
			for i in range(len(theta)):
				V[i + 1, :] = rotate(B[i, :], theta[i]) @ V[i, :]
			xn = V[:, 0]
			yn = V[:, 1]
			zn = V[:, 2]
		else:
			assert True, "not finished"

		nn = np.sqrt(xn * xn + yn * yn + zn * zn)
		xn = xn / nn
		yn = yn / nn
		zn = zn / nn
		# calculate the bi-normal
		xb = yt * zn - yn * zt
		yb = zt * xn - zn * xt
		zb = xt * yn - xn * yt
		bb = np.sqrt(xb * xb + yb * yb + zb * zb)
		xb = xb / bb
		yb = yb / bb
		zb = zb / bb
		# get the boundary lines
		z1 = self.z - width / 2 * zb + height / 2 * zn
		x1 = self.x - width / 2 * xb + height / 2 * xn
		x2 = self.x + width / 2 * xb + height / 2 * xn
		y2 = self.y + width / 2 * yb + height / 2 * yn
		z2 = self.z + width / 2 * zb + height / 2 * zn
		x3 = self.x + width / 2 * xb - height / 2 * xn
		y3 = self.y + width / 2 * yb - height / 2 * yn
		z3 = self.z + width / 2 * zb - height / 2 * zn
		x4 = self.x - width / 2 * xb - height / 2 * xn
		y4 = self.y - width / 2 * yb - height / 2 * yn
		z4 = self.z - width / 2 * zb - height / 2 * zn
		y1 = self.y - width / 2 * yb + height / 2 * yn
		# assemble
		xx = np.array([x1, x2, x3, x4, x1])
		yy = np.array([y1, y2, y3, y4, y1])
		zz = np.array([z1, z2, z3, z4, z1])
		return xx, yy, zz

	def curveBuildCoil(self,nsegments=6):
		"""Builds a coil fitting curved segments.

		This routine builds a coil via a set of segmented curves.
		To do this the code calculates the coil curvature and
		divides the existing coil into N curves where the endpoints
		of each curve coincide with a maximum in curvature. Each
		sector is fit to the original coils segments. The new coil
		x,y,z is returned.

		Parameters
		----------
		nsegments : int (optional)
			Number of coil segments(default: 6)

		Returns
		----------
		x : ndarry
			X-coordiantes [m]
		y : ndarry
			Y-coordiantes [m]
		z : ndarry
			Z-coordiantes [m]
		"""
		# First calculate the coil curvature
		self.spline_tangent(der=2)
		xn = self.xa
		yn = self.ya
		zn = self.za
		kappa = np.sqrt(xn**2+yn**2+zn**2)
		# Now determine points of maximum curvature
		# Loop over segments
		for k in range(nsegments):
			# Define starting/ending points
			# Define points to fit
			# Estimate initial guess
			# Fit curve
			# Add to total curve
			print('test')


if __name__=="__main__":
	import sys
	sys.exit(0)

