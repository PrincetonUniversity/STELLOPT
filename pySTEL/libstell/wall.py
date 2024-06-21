##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for working with wall_data
"""

# Libraries
from libstell.libstell import LIBSTELL

# Constants

# VMEC Class
class WALL(LIBSTELL):
	"""Class for working with wall files

	"""
	def __init__(self):
		self.name = None
		self.date = None
		self.nfaces = None
		self.nvertex = None
		self.vertex = None
		self.faces = None
		self.laccel = False

	def read_wall(self,filename):
		"""Directly reads a wall file

		This routine reads a wall file into the class.

		Parameters
		----------
		filename : str
			Path to wall file.
		"""
		import numpy as np
		f = open(filename,'r')
		lines = f.readlines()
		f.close()
		if  'MACHINE:' in lines[0]:
			self.name = lines[0][8:]
		else:
			print("Bad Synatx line 1 in wall file")
		if  'DATE:' in lines[1]:
			self.date = lines[1][6:]
		else:
			print("Bad Synatx line 2 in wall file")
		n1, n2 = lines[2].split()
		i1 = 3
		if (n1 == 0 and n2 == 0):
			self.laccel = True
			n1, n2 = lines[3].split()
			i1 = 4
		self.nvertex = int(n1)
		self.nfaces  = int(n2)
		self.vertex  = np.zeros((self.nvertex,3), dtype=float)
		self.faces   = np.zeros((self.nfaces,3), dtype=int)
		for i in range(self.nvertex):
			line = lines[i+i1].split()
			#print(line)
			self.vertex[i,0] = float(line[0])
			self.vertex[i,1] = float(line[1])
			self.vertex[i,2] = float(line[2])
		i1 = i1 + self.nvertex
		for i in range(self.nfaces):
			line = lines[i+i1].split()
			# note we convert to python indexing
			self.faces[i,0] = int(line[0])-1
			self.faces[i,1] = int(line[1])-1
			self.faces[i,2] = int(line[2])-1

	def write_wall(self,filename):
		"""Directly writes a wall file

		This routine writes a wall file from the class.

		Parameters
		----------
		filename : str
			Path to wall file.
		"""
		f = open(filename,'w')
		f.write(f"MACHINE: {self.name}\n")
		f.write(f"DATE: {self.date}\n")
		f.write(f"{self.nvertex} {self.nfaces}\n")
		for i in range(self.nvertex):
			f.write(f"{self.vertex[i,0]:20.10E} {self.vertex[i,1]:20.10E} {self.vertex[i,2]:20.10E}\n")
		# Note we convert back to matlab indexing
		for i in range(self.nfaces):
			f.write(f"{int(self.faces[i,0])+1} {int(self.faces[i,1])+1} {int(self.faces[i,2])+1}\n")
		f.close()

	def plot_wall_cloud(self,ax=None):
		"""Plots the vertices of the wall

		This routine plots the vertices of a wall at a point cloud.
		It takes an axis as an optional argument.

		Parameters
		----------
		ax : axis object (optional)
			Axis onto which to plot
		"""
		import matplotlib.pyplot as pyplot
		lplotnow = False
		if not ax:
			ax = pyplot.axes(projection='3d')
			lplotnow = True
		ax.scatter(self.vertex[:,0],self.vertex[:,1],self.vertex[:,2],marker='.')
		if lplotnow: pyplot.show()

	def plot_wall_3D(self,wallcolor=None,renderer=None,render_window=None):
		"""Plots a wall in 3D using VTK

		This routine plots walls in 3D using VTK

		Parameters
		----------
		wallcolor : ndarray (optional)
			Array of values to color code wall.
		renderer : vtkRenderer (optional)
			Renderer for plotting with VTK
		render_window : vtkRnderWindow (optional)
			Render window for plotting with VTK
		"""
		import numpy as np
		import matplotlib.pyplot as pyplot
		import vtk
		from vtkmodules.vtkCommonColor import vtkNamedColors
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
		# Convert numpy arrays to VTK arrays
		points = vtk.vtkPoints()
		for vertex in self.vertex:
			points.InsertNextPoint(vertex.tolist())
		triangles = vtk.vtkCellArray()
		for index in self.faces:
			triangle = vtk.vtkTriangle()
			triangle.GetPointIds().SetId(0, index[0])
			triangle.GetPointIds().SetId(1, index[1])
			triangle.GetPointIds().SetId(2, index[2])
			triangles.InsertNextCell(triangle)
		# Create a polydata object
		polydata = vtk.vtkPolyData()
		polydata.SetPoints(points)
		polydata.SetPolys(triangles)
		# Create an actor
		actor = vtk.vtkActor()
		# Add scalar values to the polydata (or make red)
		scalars=None
		if (wallcolor): 
			vals = args[0][s[k],:,:].T.flatten()
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
		else:
			colors = vtkNamedColors()
			actor.GetProperty().SetColor(0.5,0.5,0.5)
		# Create a mapper and set the scalar range to the scalar values range
		mapper = vtk.vtkPolyDataMapper()
		mapper.SetInputData(polydata)
		if scalars:
			mapper.SetLookupTable(lut)
			mapper.SetScalarRange(scalars.GetRange())
		actor.SetMapper(mapper)
		# Add actor to the scene
		renderer.AddActor(actor)
		# Render and interact
		if (scalars): renderer.AddActor2D(scalar_bar)
		renderer.SetBackground(0.1, 0.2, 0.3)
		if lplotnow:
			render_window.Render()
			render_window_interactor.Start()

	def blenderWall(self):
		"""Generates the lists Blender needs to render a wall

		This routine generates the verticies and faces lists which
		Blender needs to render a wall.

		Returns
		----------
		vertices : list
			List of tuples defining verticies
		faces: list
			List of tubles defining faces
		"""
		vertices = []
		faces = []
		for i in range(self.nvertex):
			vertices.append((self.vertex[i,0],self.vertex[i,1],self.vertex[i,2]))
		for i in range(self.nfaces):
			faces.append((int(self.faces[i,0]),int(self.faces[i,1]),int(self.faces[i,2])))
		return vertices,faces

	def genWallfromOrdered(self,vertex):
		"""Generates a wall from an orderer array

		This routine generates a wall from an ordered set of points
		where the array size is (3,nphi,ntheta). It assumes
		that datapoints are not repeated in the theta  or
		phi direction, but those coordinates are periodic.

		Parameters
		----------
		vertex : ndarray
			Ordered set of verticies
		"""
		import numpy as np
		from datetime import datetime
		temp,nphi,ntheta = vertex.shape
		nvertex = nphi*ntheta
		nfaces  = 2*nphi*ntheta
		x = np.zeros((nvertex))
		y = np.zeros((nvertex))
		z = np.zeros((nvertex))
		faces  = np.zeros((nfaces,3))
		# Do Vertices
		k=0
		for v in range(nphi):
			for u in range(ntheta):
				x[k] = vertex[0,v,u]
				y[k] = vertex[1,v,u]
				z[k] = vertex[2,v,u]
				k = k + 1
		# Do faces
		n = 0
		for v in range(nphi):
			for u in range(ntheta):
				i1 = u + v * ntheta
				# Catch special case #1
				if u == ntheta-1:
					i2 = i1 + ntheta
					i3 = i1 + 1
					i4 = i1 - ntheta + 1
					if v == nphi - 1:
						i2 = u
						i3 = 0
				elif u < ntheta-1:
					i2 = i1 + ntheta
					i3 = i1 + ntheta +1
					i4 = i1 + 1
					if v == nphi -1:
						i2 = u
						i3 = u + 1
				faces[n] = [i1, i2, i4]
				n = n + 1
				faces[n] = [i2, i3, i4]
				n = n + 1
		self.name = f"Generated using genWallfromOrdered in Python."
		self.date = datetime.today().strftime('%Y-%m-%d')
		self.nvertex = nvertex
		self.nfaces  = nfaces
		self.vertex = np.column_stack((x,y,z))
		self.faces = faces

if __name__=="__main__":
	import sys
	wall_data = WALL()
	wall_data.read_wall('QSM.dat')
	wall_data.plot_wall_3D()
	sys.exit(0)
