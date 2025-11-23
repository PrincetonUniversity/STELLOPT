#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for reading and handling 
DIAGNO synthetic magnetic diagnostic data.
"""

# Libraries
from libstell.libstell import LIBSTELL

# Constants

# DIAGNO Class
class DIAGNO():
	"""Class for working with DIAGNO data

	"""
	def __init__(self):
		test = 0

	def read_diagno_flux(self,filename):
		"""Reads a DIAGNO fluxloop data file.

		This routine reads the diagno fluxloop response data.

		Parameters
		----------
		filename : str
			Path to diagno_flux file.

		Returns
		-------
		names : list
			List of flux loop names
		flux : ndarray
			Flux [Wb]
		"""
		import numpy as np
		f = open(filename,'r')
		line = f.readline()
		nels = int(line)
		self.flux = np.zeros((nels))
		for i in range(nels):
			line = f.readline()
			self.flux[i] = float(line)
		#self.names = [f.readline()]
		#print(self.names)
		self.names = []
		for i in range(nels):
			self.names.append(f.readline().strip())
		f.close()

	def read_diagno_segrog(self,filename):
		"""Reads a DIAGNO segmented Rogowski data file.

		This routine reads the DIAGNO segmented Rogowski data. Note
		that the computed signal is in T-m^3. This is because DIAGNO
		multiplies the integral B.dl by the area of the Rogowski loops.
		This area should be A*N/L so that the returned value is in V.

		Parameters
		----------
		filename : str
			Path to diagno_segrog file.

		Returns
		-------
		names : list
			List of flux loop names
		"""
		import numpy as np
		f = open(filename,'r')
		line = f.readline()
		nels = int(line)
		self.signal = np.zeros((nels))
		for i in range(nels):
			line = f.readline()
			self.signal[i] = float(line)
		self.names = []
		for i in range(nels):
			self.names.append(f.readline().strip())
		f.close()
		return

	def read_diagno_bprobes(self,filename):
		"""Reads a DIAGNO bprobes data file.

		This routine reads the DIAGNO probe data file. Note
		that the computed signal is in T-m^3. This is because DIAGNO
		multiplies the integral B.dl by the area of the Rogowski loops.
		This area should be A*N/L so that the returned value is in V.

		Parameters
		----------
		filename : str
			Path to diagno_segrog file.

		Returns
		-------
		names : list
			List of flux loop names
		"""
		import numpy as np
		f = open(filename,'r')
		line = f.readline()
		xp=[]; yp=[]; zp=[]; modb=[]; bflux=[]
		while '#' not in line:
			(t1,t2,t3,t4,t5,t6)=line.strip()
			xp.append(float(t2))
			yp.append(float(t3))
			zp.append(float(t4))
			modb.append(float(t5))
			bflux.append(float(t6))
			line = f.readline()
		f.close()
		self.x_probe = np.array(xp)
		self.y_probe = np.array(yp)
		self.z_probe = np.array(zp)
		self.modb_probe = np.array(modb)
		self.bflux_probe = np.array(t6)
		return

# DIAGNO Class
class DIAGNO_DIAG():
	"""Class for working with DIAGNO DIAGNOSTIC files

	"""
	def __init__(self):
		test = 0

	def read_fluxloops(self,filename,scale_factor=1.0):
		"""Reads a DIAGNO fluxloop definition file.

		This routine reads the diagno fluxloop deffinition file.

		Parameters
		----------
		filename : str
			Path to diagno_fluxloops file.
		scale_factor : float
			Scale factor to apply to data (default = 1.0)
		"""
		import numpy as np
		f = open(filename,'r')
		line = f.readline()
		nloops = int(line)
		self.floops = {}
		for i in range(nloops):
			line = f.readline()
			[txt1,txt2,txt3,txt4] = line.split()
			nels = int(txt1)
			name = txt4
			xyz = np.ndarray((3,nels))
			for j in range(nels):
				line = f.readline()
				xyz[:,j]=[float(x) for x in line.split()]
			self.floops[name] = xyz
		f.close()

	def read_segrog(self,filename,scale_factor=1.0):
		"""Reads a DIAGNO segrog definition file.

		This routine reads the diagno segmented rogowski deffinition
		file.

		Parameters
		----------
		filename : str
			Path to diagno_segrog file.
		scale_factor : float
			Scale factor to apply to data (default = 1.0)
		"""
		import numpy as np
		f = open(filename,'r')
		line = f.readline()
		nloops = int(line)
		self.rog = {}
		for i in range(nloops):
			line = f.readline()
			[txt1,txt2,txt3,txt4] = line.split()
			nels = int(txt1)
			name = txt4
			xyz = np.ndarray((4,nels))
			for j in range(nels):
				line = f.readline()
				xyz[:,j]=[float(x) for x in line.split()]
			self.rog[name] = xyz[0:3,:]
		f.close()

	def plot_diagno_fluxloops(self,plot3D=None,color='purple'):
		"""Plots a diagno fluxloop in 3D using VTK

		This routine plots diagno fluxloops in 3D using VTK

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
		# Plot coils
		for key, value in self.floops.items():
			xyz0 = value[:,0]
			xyz = np.insert(value,value.shape[1],xyz0,axis=1) # add first to last
			points_array = np.swapaxes(xyz,0,1)
			# Convert numpy array to VTK points
			points = vtk.vtkPoints()
			for point in points_array:
				points.InsertNextPoint(point)
			# Add to render
			plt.add3Dline(points,color=color,linewidth=5)
		# In case it isn't set by user.
		plt.setBGcolor()
		# Render if requested
		if lplotnow: plt.render()

	def plot_diagno_segrog(self,plot3D=None,color='green'):
		"""Plots a diagno segrog in 3D using VTK

		This routine plots diagno segmented Rogowski in 3D using VTK

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
		# Plot coils
		for key, value in self.rog.items():
			xyz0 = value[:,0]
			xyz = np.insert(value,value.shape[1],xyz0,axis=1) # add first to last
			points_array = np.swapaxes(xyz,0,1)
			# Convert numpy array to VTK points
			points = vtk.vtkPoints()
			for point in points_array:
				points.InsertNextPoint(point)
			# Add to render
			plt.add3Dline(points,color=color,linewidth=5)
		# In case it isn't set by user.
		plt.setBGcolor()
		# Render if requested
		if lplotnow: plt.render()



# DIAGNO Input Class
class DIAGNO_IN():
	"""Class for working with DIAGNO_IN data

	"""
	def __init__(self, parent=None):
		self.libStell = LIBSTELL()
		self.libStell.read_diagno_in('CALLED_FROM_PYTHON')

	def read_input(self,filename):
		"""Reads DIAGNO_IN namelist from a file

		This routine wrappers the diagno_input_mod module reading routine.
		Parameters
		----------
		filename : string
			Input file name with DIAGNO_IN namelist
		"""
		diagno_in_dict = self.libStell.read_diagno_in(filename)
		for key in diagno_in_dict:
			setattr(self, key, diagno_in_dict[key])

	def write_input(self,filename):
		"""Writes DIAGNO_IN namelist to a file

		This routine wrappers the diagno_input_mod module writing routine.
		Parameters
		----------
		filename : string
			Input file name to write DIAGNO_IN namelist to
		"""
		out_dict = vars(self)
		self.libStell.write_diagno_in(filename,out_dict)

# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)