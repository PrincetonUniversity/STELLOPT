##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for working with coils files
"""

# Libraries
import libstell

# Constants

# VMEC Class
class COILSET(libstell.LIBSTELL):
	"""Class for working with coils files

	"""
	def __init__(self):
		from collections import deque
		super().__init__()
		self.libStell = libstell.LIBSTELL()
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

	def plotcoils(self,*args,**kwargs):
		"""Plots a coilset in 3D

		This routine plots coils in 3D
		"""
		import numpy as np
		import matplotlib.pyplot as pyplot
		import mpl_toolkits.mplot3d as mplot3d
		fig=kwargs.pop('fig',pyplot.figure())
		ax=kwargs.pop('axes',fig.add_subplot(111,projection='3d'))
		c_temp = self.color_cycle[0]
		for i in range(self.ngroups):
			for j in range(self.groups[i].ncoils):
				if j == 0:
					ax.plot(self.groups[i].coils[j].x, \
						   self.groups[i].coils[j].y, \
						   self.groups[i].coils[j].z, c=c_temp, label=self.groups[i].name)
				else:
					ax.plot(self.groups[i].coils[j].x, \
						   self.groups[i].coils[j].y, \
						   self.groups[i].coils[j].z, c=c_temp)
			self.color_cycle.rotate(1)
			c_temp = self.color_cycle[0]
		ax.set_xlim(self.xmin*1.05,self.xmax*1.05); ax.set_xlabel('X [m]')
		ax.set_ylim(self.ymin*1.05,self.ymax*1.05); ax.set_ylabel('Y [m]')
		ax.set_zlim(self.zmin*1.05,self.zmax*1.05); ax.set_zlabel('Z [m]')
		ax.set_title('COILS')
		ax.set_aspect('equal', adjustable='box')
		pyplot.legend(loc="upper left")
		pyplot.show()

	def plotcoilsRZ(self,*args,**kwargs):
		"""Plots each coil in the RZ plot

		This routine plots coils in 3D
		"""
		import numpy as np
		import matplotlib.pyplot as pyplot
		import mpl_toolkits.mplot3d as mplot3d
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
				offset = 0
				if j == self.groups[i].ncoils-1: offset = 1
				current = np.ones((self.groups[i].coils[j].npts))*self.groups[i].current
				current[-1] = 0
				for k in range(self.groups[i].coils[j].npts-offset):
					f.write(f"{self.groups[i].coils[j].x[k]:.10E} {self.groups[i].coils[j].y[k]:.10E} {self.groups[i].coils[j].z[k]:.10E} {current[k]:.10E}\n")
			k = self.groups[i].coils[j].npts-offset-1
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

if __name__=="__main__":
	import sys
	from argparse import ArgumentParser
	parser = ArgumentParser(description= 
		'''Provides class for accessing coils files also servers as a
		   simple tool for assessing coils or coils files.''')
	parser.add_argument("-c", "--coil", dest="coils_file",
		help="Coils file for input", default = None)
	parser.add_argument("-p", "--plot", dest="lplot", action='store_true',
		help="Plot the coils file.", default = False)
	parser.add_argument("-prz", "--plotRZ", dest="lplotRZ", action='store_true',
		help="Plot each coil group in RZ.", default = False)
	parser.add_argument("-b", "--bfield", dest="bxyz",
		help="Output B field at x,y,z", default = None)
	parser.add_argument("-a", "--afield", dest="axyz",
		help="Output A field at x,y,z", default = None)	
	parser.add_argument("-o", "--output", dest="loutput", action='store_true',
		help="Output the coil", default = False)
	args = parser.parse_args()
	coils = COILSET()
	if args.coils_file: 
		coils.read_coils_file(args.coils_file)
		if args.lplot: coils.plotcoils()
		if args.lplotRZ: coils.plotcoilsRZ()
		if args.loutput: coils.write_coils_file(args.coils_file+'_new')
		if args.axyz:
			x,y,z = args.axyz.split(',')
			ax,ay,az = coils.coilvecpot(float(x),float(y),float(z))
			print(f"Vector Potential ({x},{y},{z}) : {ax}, {ay}, {az} ")
		if args.bxyz:
			x,y,z = args.bxyz.split(',')
			bx,by,bz = coils.coilbiot(float(x),float(y),float(z))
			print(f"B-Field ({x},{y},{z}) : {bx}, {by}, {bz} [T]")
	sys.exit(0)

