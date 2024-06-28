##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class simplifiying the VTK interface
"""

# Libraries
import vtk
from vtkmodules.vtkCommonColor import vtkNamedColors

# Constants
_WPIX_ = 1024
_HPIX_ = 768
_R_RGB_ = 0.1
_G_RGB_ = 0.2
_B_RGB_ = 0.3
_CMAP_  = 'jet'

# VMEC Class
class PLOT3D():
	"""Class for working with wall files

	"""
	def __init__(self,lwindow=True):
		self.cmap = _CMAP_
		self.renderer = vtk.vtkRenderer()
		self.render_window = vtk.vtkRenderWindow()
		self.scalar_bar = vtk.vtkScalarBarActor()
		self.lookupTable = None
		if lwindow:
			self.render_window_interactor = vtk.vtkRenderWindowInteractor()
			self.setRendererWindow(self.render_window)
			self.render_window.SetSize(_WPIX_, _HPIX_)
		else:
			self.render_window_interactor = None


	def setRenderer(self,renderer):
		"""Set the renderer

		The routine allows one to set the renderer object for VTK.

		Parameters
		----------
		r : VTK Renderer
			Renderer Object from VTK
		"""
		self.renderer = renderer

	def setRendererWindow(self,render_window):
		"""Set the renderer window

		The routine allows one to set the render window object for VTK.

		Parameters
		----------
		r : VTK Render Window
			Renderer Window Object from VTK
		"""
		self.render_window = render_window
		self.render_window.AddRenderer(self.renderer)
		self.render_window_interactor.SetRenderWindow(self.render_window)

	def setLookupTable(self,lut):
		"""Set the color lookupTable

		The routine allows one to set the renderer lookupTable

		Parameters
		----------
		lut : VTK LookupTable
			Lookup Table object from VTK
		"""
		self.lookupTable = lut
		actors = self.renderer.GetActors()
		for actor in actors:
			mapper = actor.GetMapper()
			if hasattr(mapper,'LookupTable'):
				mapper.SetLookupTable(lut)


	def setBGcolor(self,r=_R_RGB_,g=_G_RGB_,b=_B_RGB_):
		"""Set the background color for the render

		The routine takes the R, G, and B color channels as
		values from 0 to 1 and sets the background color
		of the render.  Calling with no parameters sets
		the value to a default.

		Parameters
		----------
		r : float
			Red color channel [0.0,1.0]
		g : float
			Green color channel [0.0,1.0]
		b : float
			Blue color channel [0.0,1.0]
		"""
		self.renderer.SetBackground(r,g,b)

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

	def vtkColor(self,j):
		"""Helper to vary colors

		This routine returns a VTK color object based on a number.
		This is done to mimic the behavior of changing colors when
		plotting lines in matplotlib.

		Parameters
		----------
		j : int
			Index into predefined colors
		"""
		from vtkmodules.vtkCommonColor import vtkNamedColors
		colors = vtkNamedColors()
		color_txt=['red','green','blue','yellow','magenta','cyan','aqua']
		j = j % len(color_txt)
		return colors.GetColor3d(color_txt[j])

	def setClim(self,cmin,cmax,llast=False):
		"""Set the color limits for all actors

		Sets the color limits for all VTK actors using common limits

		Parameters
		----------
		cmin : float
			Min value for color limits
		cmax : float
			Max value for color limits
		llast : boolean
			Only for last actor (default: False)
		"""
		if llast:
			actor = actors.GetLastActor()
			mapper = actor.GetMapper()
			mapper.SetScalarRange(cmin,cmax)
		else:
			actors = self.renderer.GetActors()
			for actor in actors:
				mapper = actor.GetMapper()
				mapper.SetScalarRange(cmin,cmax)

	def torusvertexTo3Dmesh(self,x,y,z,lcloseu=True,lclosev=True):
		"""Generate points and triangle objects from x,y,z data for a torus

		This routine returns a vtkPoints object and vtkCellArray
		given a vertex array and indicies array. There is an assumption that
		f[0,:] /= f[nu,:] and f[:,0] /= f[:,nv] for x, y, and z. The lcloseu
		and lclosev arguments can be used to tell the routine that faces
		closing the torus in the u and v direction should not be generated.

		Parameters
		----------
		x : ndarray
			X torus points to plot
		y : ndarray
			Y torus points to plot
		z : ndarray
			Z torus points to plot
		lcloseu : boolean (optional)
			Close surface in poloidal direction (default: True)
		lclosev : boolean (optional)
			Close surface in toroidal direction (default: True)
		Returns
		-------
		points : VTK Points object
			Points object for 3D meshes.
		triangle : VTK Cell Array object
			Triangle object for 3D meshes.
		"""
		import numpy as np
		# Create objects
		points = vtk.vtkPoints()
		triangles = vtk.vtkCellArray()
		# Get sizes
		nu = np.size(x,0)
		nv = np.size(x,1)
		# Loop
		for v in range(nv):
			for u in range(nu):
				points.InsertNextPoint([x[u,v],y[u,v],z[u,v]])
				if (u == nu - 1) and (not lcloseu): continue
				if (v == nv - 1) and (not lclosev): continue
				i1 = u + v * nu
				i2 = i1 + nu
				i3 = i1 + nu + 1
				i4 = i1 + 1
				if (v == (nv - 1)):
					i2 = u
					i3 = u + 1
				if (u == (nu - 1)):
					i3 = i1 + 1
					i4 = i1 - nu + 1
					if (v == nv - 1): i3 = 0 
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
		return points, triangles

	def valuesToScalar(self,vals):
		"""Generate scalar object from values

		This routine returns a scalar VTK object from a set
		of values for use in coloring objects.

		Parameters
		----------
		vals : ndarray
			Values to use for coloring.
		Returns
		-------
		scalars : VTK FloatArray
			Scalar values for coloring objects.
		"""
		scalars = vtk.vtkFloatArray()
		scalars.SetNumberOfComponents(1)
		for value in vals.T.flatten():
			scalars.InsertNextValue(value)
		return scalars

	def facemeshTo3Dmesh(self,vertices,indices):
		"""Generate points and triangle objects from a facemesh

		This routine returns a vtkPoints object and vtkCellArray
		given a vertex array and indicies array.

		Parameters
		----------
		vertex : ndarray
			Points to plot (npts,3)
		indices : ndarray
			Face list (nfaces,3)
		Returns
		-------
		points : VTK Points object
			Points object for 3D meshes.
		triangle : VTK Cell Array object
			Triangel object for 3D meshes.
		"""
		# Create objects
		points = vtk.vtkPoints()
		triangles = vtk.vtkCellArray()
		# Convert numpy arrays to VTK arrays
		for vertex in vertices:
			points.InsertNextPoint(vertex.tolist())
		for index in indices:
			triangle = vtk.vtkTriangle()
			triangle.GetPointIds().SetId(0, index[0])
			triangle.GetPointIds().SetId(1, index[1])
			triangle.GetPointIds().SetId(2, index[2])
			triangles.InsertNextCell(triangle)
		return points,triangles

	def add3Dline(self,points,scalars=None,linewidth=None,color='red'):
		"""Add a 3D line to a render

		This routine adds a line using VTK where points is an object
		as returned by vtk.vtkPoints().

		Parameters
		----------
		points : VTK Points object
			Points to plot
		scalars : VTK FloatArray object (optional)
			Color values to plot
		linewidth : float (optional)
			Width of line to plot.
		color : string (optional)
			Line color name, see VTK (scalars overrides)
		"""
		# Create actor/mapper
		actor = vtk.vtkActor()
		mapper = vtk.vtkPolyDataMapper()
		# Create a polyline to connect the points
		polyline = vtk.vtkPolyLine()
		npts=points.GetNumberOfPoints()
		polyline.GetPointIds().SetNumberOfIds(npts)
		for k in range(npts):
			polyline.GetPointIds().SetId(k, k)
		# Create a cell array to store the polyline
		cells = vtk.vtkCellArray()
		cells.InsertNextCell(polyline)
		# Create a polydata object and add the points and the polyline to it
		polydata = vtk.vtkPolyData()
		polydata.SetPoints(points)
		polydata.SetLines(cells)
		# Setup the mapper
		mapper.SetInputData(polydata)
		# Handle scalars or use line color
		if type(scalars) != type(None): 
			polydata.GetPointData().SetScalars(scalars)
			self.lookupTable = self.vtkLUTHelper(self.cmap)
			mapper.SetLookupTable(self.lookupTable)
			mapper.SetScalarRange(scalars.GetRange())
		else:
			colors = vtkNamedColors()
			actor.GetProperty().SetColor(colors.GetColor3d(color))
		# Set line thickness
		if linewidth: actor.GetProperty().SetLineWidth(linewidth)
		# Set Mapper
		actor.SetMapper(mapper)
		# Add actor
		self.renderer.AddActor(actor)

	def add3Dmesh(self,points,triangles,scalars=None,opacity=1.0,color='red'):
		"""Add a 3D mesh to a render

		This routine adds a mesh using VTK where points is an object
		as returned by vtk.vtkPoints() and triangles is an object as
		returned by vtk.CellArray()

		Parameters
		----------
		points : VTK Points object
			Points to plot
		triangles : VTK CellArray object
			Face list
		scalars : VTK FloatArray object (optional)
			Color values to plot
		opacity : float (optional)
			Surface opacity [0.0,1.0] (default:1.0)
		color : string (optional)
			Surface color name, see VTK (scalars overrides)
		"""
		from vtkmodules.vtkCommonColor import vtkNamedColors
		# Create actor/mapper
		actor = vtk.vtkActor()
		mapper = vtk.vtkPolyDataMapper()
		# Create polydata
		polydata = vtk.vtkPolyData()
		polydata.SetPoints(points)
		polydata.SetPolys(triangles)
		# Link Mapper to polydata
		mapper.SetInputData(polydata)
		# Handle scalars or make red
		if type(scalars) != type(None): 
			polydata.GetPointData().SetScalars(scalars)
			if not self.lookupTable: self.lookupTable = self.vtkLUTHelper(self.cmap)
			mapper.SetLookupTable(self.lookupTable)
			mapper.SetScalarRange(scalars.GetRange())
		else:
			colors = vtkNamedColors()
			actor.GetProperty().SetColor(colors.GetColor3d(color))
		# Set Opacity
		actor.GetProperty().SetOpacity(opacity)
		# Set Mapper
		actor.SetMapper(mapper)
		# Add actor to the scene
		self.renderer.AddActor(actor)

	def colorbar(self,show=True,title=""):
		"""Add a colorbar to a render

		This routine adds a colorbar to the render.

		Parameters
		----------
		show : logical (optional)
			Show the colorbar (default: True)
		title : string
			Colorbar title
		"""
		if show:
			self.scalar_bar.SetLookupTable(self.lookupTable)
			self.scalar_bar.SetTitle(title)
			self.scalar_bar.SetNumberOfLabels(5)
			self.renderer.AddActor2D(self.scalar_bar)

	def render(self):
		"""Render the window

		This routine renders the window.
		"""
		self.render_window.Render()
		self.render_window_interactor.Start()

if __name__=="__main__":
	import sys
	sys.exit(0)