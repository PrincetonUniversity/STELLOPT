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
_R_RGB_ = 1.0
_G_RGB_ = 1.0
_B_RGB_ = 1.0
_CMAP_  = 'jet'
_FONTSIZE_ = 24

# VMEC Class
class PLOT3D():
	"""Class for working with wall files

	"""
	def __init__(self,lwindow=True):
		import yaml
		import os
		from pathlib import Path
		self.cmap = None
		self.renderer = vtk.vtkRenderer()
		self.render_window = vtk.vtkRenderWindow()
		self.scalar_bar = vtk.vtkScalarBarActor()
		self.camera = vtk.vtkCamera()
		self.lookupTable = None
		self.colordex = 0
		self.fontsize = _FONTSIZE_
		if lwindow:
			self.render_window_interactor = vtk.vtkRenderWindowInteractor()
			self.setRendererWindow(self.render_window)
			self.render_window.SetSize(_WPIX_, _HPIX_)
		else:
			self.render_window_interactor = None
		cfg_file=os.path.join(Path.home(),'.pystelrc')
		if os.path.isfile(cfg_file):
			with open(cfg_file, "r") as ymlfile:
				cfg = yaml.safe_load(ymlfile)
			if 'mycolors' in cfg.keys(): self.colororder = cfg['mycolors']
			if 'cmaps' in cfg.keys(): self.cmap = cfg['cmaps']
			if 'fontsize' in cfg.keys(): self.fontsize = cfg['fontsize']
		else:
			self.colororder = ['red','green','blue','yellow','magenta','cyan','aqua']

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

	def vtkLUTHelper(self,ctable=None):
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
		#from matplotlib import cm
		import matplotlib
		lut = vtk.vtkLookupTable()
		# Check for custom colormaps	
		if ctable in list(matplotlib.colormaps.keys()):
			lut.SetNumberOfTableValues(256)
			lut.Build()
			# Use matplotlib to generate the jet colormap
			mlibmap = matplotlib.colormaps[ctable]
			for i in range(256):
				rgba = mlibmap(i / 255.0)
				lut.SetTableValue(i, rgba[0], rgba[1], rgba[2], rgba[3])
		elif type(self.cmap) is dict:
			if ctable in list(self.cmap.keys()):
				cmap = self.cmap[ctable]
			else:
				# Use first colormap
				res = list(self.cmap.keys())[0]
				cmap = self.cmap[res]
			# Now construct
			ncolors = len(cmap)
			lut.SetNumberOfTableValues((ncolors-1)*256)
			lut.Build()
			k = 0
			for i in range(ncolors-1):
				c1 = cmap[i]
				c2 = cmap[i+1]
				for j in range(256):
					r = (c2[0]-c1[0])*float(j)/255.0 + c1[0]
					g = (c2[1]-c1[1])*float(j)/255.0 + c1[1]
					b = (c2[2]-c1[2])*float(j)/255.0 + c1[2]
					lut.SetTableValue(k,r,g,b,1.0)
					k = k + 1
		else:
			ctable = _CMAP_ # default from above
			lut.SetNumberOfTableValues(256)
			lut.Build()
			# Use matplotlib to generate the jet colormap
			mlibmap = matplotlib.colormaps[ctable]
			for i in range(256):
				rgba = mlibmap(i / 255.0)
				lut.SetTableValue(i, rgba[0], rgba[1], rgba[2], rgba[3])
		return lut

	def setActorColor(self,actor,color=None):
		"""Set the color property of the actor

		This routine helps set the color of an actor. It does so by
		cycling through the colororder property of this class. If the
		user passes a specific color then the named color is used.

		Parameters
		----------
		actor : VTK Actor Object
			Actor to set color property of
		color : str (optional)
			Use the named color instead
		"""
		cNamed = vtkNamedColors()
		if not color:
			ctemp = self.colororder[self.colordex]
			self.colordex = self.colordex + 1
			if self.colordex >= len(self.colororder):
				self.colordex = 0
			if type(ctemp) == type('str'):
				vtkcolor = cNamed.GetColor3d(ctemp)
			else:
				vtkcolor = vtk.vtkColor3d(ctemp[0],ctemp[1],ctemp[2])
		elif type(color) == type('str'):
			vtkcolor = cNamed.GetColor3d(color)
		else:
			vtkcolor = cNamed.GetColor3d('red')
		actor.GetProperty().SetColor(vtkcolor)

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

	def setCamera(self,pos=None,focus=None,camup=None,angle=None):
		"""Set the camera properties

		The routine allows the user to set the camera properties of a
		VTK render.

		Parameters
		----------
		pos : list (optional)
			Set position of camera x,y,z
		focus : list (optional)
			Set focus of camera nx,ny,nz
		camup : list (optional)
			Set up-vector of camera ux,uy,uz
		angle : float (optional)
			Set the camera viewing angle. [deg]
		"""
		camera = self.renderer.GetActiveCamera()
		if type(pos) is not type(None): camera.SetPosition(pos[0], pos[1], pos[2])  # Set Camera position
		if type(focus) is not type(None): camera.SetFocalPoint(focus[0], focus[1], focus[2])  # Set Camera focal point
		if type(camup) is not type(None): camera.SetViewUp(camup[0], camup[1], camup[2])  # Set Camera up vector
		if type(angle) is not type(None): camera.SetViewAngle(angle)  # Set Camera viewing angle [deg]

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

	def vertexToPoints(self,vertices):
		"""Generate points objects from an array

		This routine returns a vtkPoints object given a vertex array.

		Parameters
		----------
		vertex : ndarray
			Points to plot (npts,3)
		Returns
		-------
		points : VTK Points object
			Points object for 3D meshes.
		"""
		# Create objects
		points = vtk.vtkPoints()
		# Convert numpy arrays to VTK arrays
		for vertex in vertices:
			points.InsertNextPoint(vertex.tolist())
		return points


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
		# Convert numpy arrays to VTK arrays
		points = self.vertexToPoints(vertices)
		# Create objects
		triangles = vtk.vtkCellArray()
		# Calculate indices
		for index in indices:
			triangle = vtk.vtkTriangle()
			triangle.GetPointIds().SetId(0, index[0])
			triangle.GetPointIds().SetId(1, index[1])
			triangle.GetPointIds().SetId(2, index[2])
			triangles.InsertNextCell(triangle)
		return points,triangles

	def tetrameshTo3DTetra(self,vertices,indices):
		"""Generate points and tetrahedron objects from a tetramesh

		This routine returns a vtkPoints object and vtkCellArray
		given a vertex array and indicies array.

		Parameters
		----------
		vertex : ndarray
			Points to plot (npts,3)
		indices : ndarray
			Face list (nfaces,4)
		Returns
		-------
		points : VTK Points object
			Points object for 3D meshes.
		tetra : VTK Cell Array object
			Tetrahedron object for 3D meshes.
		"""
		# Convert numpy arrays to VTK arrays
		points = self.vertexToPoints(vertices)
		# Create objects
		tetra = vtk.vtkCellArray()
		# Convert numpy arrays to VTK arrays
		for index in indices:
			tet = vtk.vtkTetra()
			tet.GetPointIds().SetId(0, index[0])
			tet.GetPointIds().SetId(1, index[1])
			tet.GetPointIds().SetId(2, index[2])
			tet.GetPointIds().SetId(3, index[3])
			tetra.InsertNextCell(tet)
		return points,tetra

	def add3Dpoints(self,points,pointsize=0.01,scalars=None,color=None):
		"""Adds 3D point cloud to a render

		This routine adds a point cloud using VTK where points is an object
		as returned by vtk.vtkPoints().

		Parameters
		----------
		points : VTK Points object
			Points to plot
		scalars : VTK FloatArray object (optional)
			Color values to plot
		pointsize : float (optional)
			Point size (default=1)
		color : string (optional)
			Line color name, see VTK (scalars overrides)
		"""
		# Create actor,mapper,sphere,glyph
		actor = vtk.vtkActor()
		mapper = vtk.vtkPolyDataMapper()
		sphere = vtk.vtkSphereSource()
		glyph = vtk.vtkGlyph3D()
		# Create a polydata object and add the points and the polyline to it
		polydata = vtk.vtkPolyData()
		polydata.SetPoints(points)
		# Setup the mapper
		mapper.SetInputData(polydata)
		# Create the glyph
		sphere.SetRadius(pointsize)
		glyph.SetSourceConnection(sphere.GetOutputPort())
		glyph.SetInputData(polydata)
		glyph.SetColorModeToColorByScalar()
		glyph.SetScaleModeToDataScalingOff()  # Disable scaling by scalar values
		mapper.SetInputConnection(glyph.GetOutputPort())
		# Handle scalars or use line color
		if type(scalars) != type(None): 
			polydata.GetPointData().SetScalars(scalars)
			self.lookupTable = self.vtkLUTHelper(self.cmap)
			mapper.SetLookupTable(self.lookupTable)
			mapper.SetScalarRange(scalars.GetRange())
		else:
			self.setActorColor(actor,color)
		# Set Mapper
		actor.SetMapper(mapper)
		# Add actor
		self.renderer.AddActor(actor)

	def add3Dline(self,points,scalars=None,linewidth=None,color=None):
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
			self.setActorColor(actor,color)
		# Set line thickness
		if linewidth: actor.GetProperty().SetLineWidth(linewidth)
		# Set Mapper
		actor.SetMapper(mapper)
		# Add actor
		self.renderer.AddActor(actor)

	def add3Dmesh(self,points,triangles,scalars=None,opacity=1.0,color=None,FaceScalars=None):
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
			if not self.lookupTable: self.lookupTable = self.vtkLUTHelper(color)
			mapper.SetLookupTable(self.lookupTable)
			mapper.SetScalarRange(scalars.GetRange())
		elif type(FaceScalars) != type(None):
			FaceScalars.SetName("FaceValues")
			polydata.GetCellData().SetScalars(FaceScalars)
			if not self.lookupTable: self.lookupTable = self.vtkLUTHelper(color)
			mapper.SetLookupTable(self.lookupTable)
			mapper.SetScalarRange(FaceScalars.GetRange())
			mapper.SetScalarModeToUseCellData()
		else:
			self.setActorColor(actor,color)
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
			# Set Label Text
			self.scalar_bar.GetLabelTextProperty().SetColor(0,0,0)
			self.scalar_bar.GetLabelTextProperty().ShadowOff()
			self.scalar_bar.GetLabelTextProperty().SetFontSize(self.fontsize)
			# Set Titel Text
			self.scalar_bar.GetTitleTextProperty().SetColor(0,0,0)
			self.scalar_bar.GetTitleTextProperty().ShadowOff()
			self.scalar_bar.GetTitleTextProperty().SetFontSize(self.fontsize)
			self.scalar_bar.UnconstrainedFontSizeOn()
			self.scalar_bar.SetLookupTable(self.lookupTable)
			self.scalar_bar.SetTitle(title)
			self.scalar_bar.SetNumberOfLabels(5)
			self.scalar_bar.SetBarRatio(0.25)
			self.renderer.AddActor2D(self.scalar_bar)

	def setView(self,az=None,el=None):
		"""Set the view based on azimuth and elevation

		This routine sets the camera position using position and focal
		point.

		Parameters
		----------
		az : float (optional)
			Azimuth about focal point (degrees)
		el : float (optional)
			Elevation about focal point (degrees)
		"""
		if az: self.camera.Azimuth(az)
		if el: self.camera.Elevation(el)
		self.camera.OrthogonalizeViewUp()
		self.renderer.SetActiveCamera(self.camera)

	def setZoom(self,zoom=None):
		"""Set the view zoom

		This routine sets the camera zoom using a factor.

		Parameters
		----------
		zoom : float (optional)
			Zoom factor (degrees)
		"""
		if zoom: self.camera.Zoom(zoom)
		self.renderer.SetActiveCamera(self.camera)

	def render(self):
		"""Render the window

		This routine renders the window.
		"""
		self.render_window.Render()
		self.render_window_interactor.Start()

	def saveImage(self,filename='vtkImage.png'):
		"""Save VTK Render as image

		This routine saves the current VTK render as an image based
		on the filename provided.

		Parameters
		----------
		filename : string (optional)
			Name of file (bmp,jpg,pnm,ps,tiff,png are supported)
		"""
		from vtkmodules.vtkIOImage import vtkBMPWriter,vtkJPEGWriter,vtkPNGWriter,vtkPNMWriter,vtkPostScriptWriter,vtkTIFFWriter
		from vtkmodules.vtkRenderingCore import vtkWindowToImageFilter
		rgba = True
		if '.bmp' in filename:
			writer = vtkBMPWriter()
		elif '.jpg' in filename:
			writer = vtkJPEGWriter()
		elif '.pnm' in filename:
			writer = vtkPNMWriter()
		elif '.ps' in filename:
			if rgba:
				rgba = False
			writer = vtkPostScriptWriter()
		elif '.tiff' in filename:
			writer = vtkTIFFWriter()
		elif '.png' in filename:
			writer = vtkPNGWriter()
		else:
			print(f'Unsupported file extension: {filename}')
			return
		windowto_image_filter = vtkWindowToImageFilter()
		windowto_image_filter.SetInput(self.render_window)
		windowto_image_filter.SetScale(1)  # image quality
		windowto_image_filter.Update()
		if rgba:
			windowto_image_filter.SetInputBufferTypeToRGBA()
		else:
			windowto_image_filter.SetInputBufferTypeToRGB()
			# Read from the front buffer.
			windowto_image_filter.ReadFrontBufferOff()
			windowto_image_filter.Update()
		writer.SetFileName(filename)
		writer.SetInputConnection(windowto_image_filter.GetOutputPort())
		writer.Write()

if __name__=="__main__":
	import sys
	sys.exit(0)