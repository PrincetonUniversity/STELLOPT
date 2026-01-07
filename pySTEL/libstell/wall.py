##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for working with wall_data
"""

# Libraries

# Constants

# WALL Class
class WALL():
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

		This routine reads a wall file into the class. Note if the
		file ends in STL then the numpy-stl routine is used to
		read the file.

		Parameters
		----------
		filename : str
			Path to wall file.
		"""
		import numpy as np
		from stl import mesh
		import re
		from datetime import datetime
		if '.stl' in filename:
			mesh_data = mesh.Mesh.from_file(filename)
			self.vertex = mesh_data.vectors.reshape((-1, 3))
			self.faces = np.arange(len(self.vertex)).reshape((-1, 3))
			byte_string = mesh_data.name
			string = byte_string.decode('utf-8')
			match = re.search(r'\d{4}-\d{2}-\d{2}', string)
			self.name  = string
			if match: 
				self.date  = match.group()
			else:
				self.date = datetime.today().strftime('%Y-%m-%d')
			self.nvertex = self.vertex.shape[0]
			self.nfaces = self.faces.shape[0]
			return
		f = open(filename,'r')
		lines = f.readlines()
		f.close()
		if  'MACHINE:' in lines[0]:
			self.name = lines[0][8:].strip()
		else:
			print("Bad Synatx line 1 in wall file")
		if  'DATE:' in lines[1]:
			self.date = lines[1][6:].strip()
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

	def write_wall_stl(self,filename):
		"""Directly writes a wall STL file.

		This routine writes a wall STL file from the class.

		Parameters
		----------
		filename : str
			Path to wall file.
		"""
		import numpy as np
		from stl import mesh
		wall_mesh = mesh.Mesh(np.zeros(self.nfaces, dtype=mesh.Mesh.dtype))
		for i, f in enumerate(self.faces):
			for j in range(3):
				wall_mesh.vectors[i][j] = self.vertex[f[j],:]
		wall_mesh.save(filename)

	def write_wall_kisslinger(self,filename,nphi,nfp=1):
		"""Creates a wall file in Kisslinger format

		This routine makes use of the meshcut library to produce
		a Kisslinger format wall from the existing wall. Such a format
		is also used by TRAVIS.

		Parameters
		----------
		filename : str
			Name of file to output
		nphi : int
			Number of toroidal cuts per field period
		nfp : int
			Number of field period (default: 1)
		"""
		import meshcut
		import numpy as np
		npts = 128
		plane_orig = (0.0,0.0,0.0)
		phiarr = np.linspace(0,np.pi*2/nfp,int(nphi))
		sout = np.linspace(0.0,1.0,npts)
		f = open(filename,'w')
		f.write(self.name+" <This_is_the_vessel_file>\n")
		rshift = 0.0
		zshift = 0.0
		dshift = 0.0
		icol   = 1
		scalf  = 0.01
		isyt   = 1
		isyp   = 1
		# Note not sure what last two values in kisslinger format are
		f.write(f"{int(nphi)} {int(npts)} {int(nfp)} {rshift} {zshift} {dshift} {icol} {scalf} {isyt} {isyp}\n")
		for phi in phiarr:
			f.write(f"{180.0*phi/np.pi}\n")
			nx = -np.sin(phi)
			ny = np.cos(phi)
			plane_normal = (nx,ny,0.0)
			mesh = meshcut.cross_section(self.vertex,self.faces, \
				plane_orig=plane_orig,plane_normal=plane_normal)
			submesh = mesh[0]
			submesh = np.concatenate((submesh,[submesh[0,:]]))
			s = np.linspace(0.0,1.0,submesh.shape[0])
			x = np.interp(sout,s,submesh[:,0])
			y = np.interp(sout,s,submesh[:,1])
			z = np.interp(sout,s,submesh[:,2])
			x[-1] = x[0]
			y[-1] = y[0]
			z[-1] = z[0]
			r = np.sqrt(x*x+y*y)
			for i in range(npts):
				#f.write(f"{x[i]*100.0} {y[i]*100.0} {z[i]*100.0}\n")
				f.write(f"{r[i]/scalf} {z[i]/scalf}\n")

	def wallAdd(self,wall_in):
		"""Add a wall to this wall

		This routine adds a wall structure to this wall structure.

		Parameters
		----------
		wall_in : WALL obj
			Wall to add to this wall.
		"""
		import numpy as np
		verts1 = self.vertex.tolist()
		verts2 = wall_in.vertex.tolist()
		verts1.extend(verts2)
		faces1 = self.faces.tolist()
		faces2 = (wall_in.faces+self.nvertex).tolist()
		faces1.extend(faces2)
		self.faces = np.array(faces1, dtype=int)
		self.vertex = np.array(verts1)
		self.nfaces = self.faces.shape[0]
		self.nvertex = self.vertex.shape[0]

	def wallClean(self):
		"""Clean wall of zero area elements

		This routine removes any zero area elements in the wall mesh.
		"""
		import numpy as np
		i0 = self.faces[:,0]
		i1 = self.faces[:,1]
		i2 = self.faces[:,2]
		V1 = self.vertex[i1,:]-self.vertex[i0,:]
		V2 = self.vertex[i2,:]-self.vertex[i0,:]
		N  = np.zeros((self.nfaces,3), dtype=float)
		N[:,0]  = V1[:,1]*V2[:,2] - V1[:,2]*V2[:,1]
		N[:,1]  = V1[:,2]*V2[:,0] - V1[:,0]*V2[:,2]
		N[:,2]  = V1[:,0]*V2[:,1] - V1[:,1]*V2[:,0]
		A = np.sum(N*N,axis=1)
		mask = A<=0
		print(f'Removing {np.sum(mask)} bad triangles.')
		self.faces = np.delete(self.faces,mask,axis=0)
		self.nfaces = self.faces.shape[0]

	def cut_wall_RZ(self,phi=0):
		"""Creates an R/Z cut of the wall

		This routine plots a cut of the wall at constant phi angle.

		Parameters
		----------
		phi : float
			Toroidal angle [rad] (default = 0)

		Returns
		----------
		R : float
			Major radius points along cut [m]
		Z : float
			Vertical points along cut [m]
		"""
		import meshcut
		import numpy as np
		# Make cut
		nx = -np.sin(phi)
		ny = np.cos(phi)
		plane_orig = (0.0,0.0,0.0)
		plane_normal = (nx,ny,0.0)
		mesh = meshcut.cross_section(self.vertex,self.faces, \
				plane_orig=plane_orig,plane_normal=plane_normal)
		R_out = []
		Z_out = []
		for submesh in mesh:
			R = np.sqrt(submesh[:,0]**2+submesh[:,1]**2)
			Z = submesh[:,2]
			if all(R > 0.0):
				R_out.append(R)
				Z_out.append(Z)
		return R_out,Z_out

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

	def plot_wall_2D(self,phi=0,ax=None):
		"""Plots an RZ cut of a wall

		This routine plots a cut of the wall at constant phi angle.

		Parameters
		----------
		phi : float
			Toroidal angle [rad] (default = 0)
		ax : axis object (optional)
			Axis onto which to plot
		"""
		import matplotlib.pyplot as pyplot
		import meshcut
		import numpy as np
		lplotnow = False
		if not ax:
			ax = pyplot.axes()
			lplotnow = True
		# Make cut
		nx = -np.sin(phi)
		ny = np.cos(phi)
		plane_orig = (0.0,0.0,0.0)
		plane_normal = (nx,ny,0.0)
		mesh = meshcut.cross_section(self.vertex,self.faces, \
				plane_orig=plane_orig,plane_normal=plane_normal)
		for submesh in mesh:
			R = np.sqrt(submesh[:,0]**2+submesh[:,1]**2)
			Z = submesh[:,2]
			if all(R > 0.0):
				ax.plot(R,Z,'k')
		ax.set_aspect('equal')
		if lplotnow: pyplot.show()

	def plot_wall_3D(self,wallcolor=None,plot3D=None):
		"""Plots a wall in 3D using VTK

		This routine plots walls in 3D using VTK

		Parameters
		----------
		wallcolor : ndarray (optional)
			Array of values to color code wall.
		plot3D : plot3D object (optional)
			Plotting object to render to.
		"""
		from libstell.plot3D import PLOT3D
		# Handle optionals
		if plot3D: 
			lplotnow=False
			plt = plot3D
		else:
			lplotnow = True
			plt = PLOT3D()
		# Generate VTK objects
		[points, triangles]=plt.facemeshTo3Dmesh(self.vertex,self.faces)
		# Generate Wall colors
		scalar = None
		if type(wallcolor) != type(None): 
			if isinstance(wallcolor, str):
				plt.add3Dmesh(points,triangles,color=wallcolor)
			else:
				scalar = plt.valuesToScalar(wallcolor)
				plt.add3Dmesh(points,triangles,FaceScalars=wallcolor)
		else:
			plt.add3Dmesh(points,triangles,color='gray')
		# Render if requested
		if lplotnow: plt.render()

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
		faces  = np.zeros((nfaces,3),dtype=int)
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
		self.laccel = False

	def refineWall(self,dlmin=0.001,dlmax=0.10,info=1):
		"""Remeshes a wall using GMSH

		This routine remeshes a wall using GMSH based on a minimum
		and maximum grid size. Screen output can be controled via the
		infor parameter:
			0: No screen output
			1: Errors only (default)
			2: Errors and warnings
			3: Info, Errors and warnings

		Parameters
		----------
		dlmin : float (optional)
			Minimum mesh size [m] (default=0.001)
		dlmax : float (optional)
			Maximum mesh size [m] (default=0.100)
		info : int (optional)
			GMSH screen output (default=1)
		"""
		import numpy as np
		import gmsh
		gmsh.initialize()
		gmsh.option.setNumber("General.Terminal", 1)
		gmsh.model.add("custom_mesh")
		# Step 1: Add vertices (need to recode here)
		node_tags = []
		for i in range(self.nvertex):
			tag = gmsh.model.geo.addPoint(self.vertex[i,0],self.vertex[i,1],self.vertex[i,2])
			node_tags.append(tag)
		# Step 2: Add faces as triangles
		for i in range(self.nfaces):
			v1, v2, v3 = self.faces[i,:] + 1 # GMSH uses 1 based indexing
			l1 = gmsh.model.geo.addLine(node_tags[v1 - 1], node_tags[v2 - 1])
			l2 = gmsh.model.geo.addLine(node_tags[v2 - 1], node_tags[v3 - 1])
			l3 = gmsh.model.geo.addLine(node_tags[v3 - 1], node_tags[v1 - 1])
			loop = gmsh.model.geo.addCurveLoop([l1, l2, l3])
			gmsh.model.geo.addPlaneSurface([loop])
		# Step 3: Set Mesh min/max lengths
		gmsh.option.setNumber("Mesh.CharacteristicLengthMin", dlmin)
		gmsh.option.setNumber("Mesh.CharacteristicLengthMax", dlmax)
		# Step 4: Synchronize and mesh
		gmsh.model.geo.synchronize()
		gmsh.model.mesh.generate(2)
		# Step 5: get mesh
		node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
		vertices = [(node_coords[i], node_coords[i+1], node_coords[i+2]) for i in range(0, len(node_coords), 3)]
		# Get all faces (triangles) in the mesh
		faces = []
		for dim, tag in gmsh.model.getEntities(2):  # Dimension 2 entities are surfaces
			element_types, element_tags, node_tags = gmsh.model.mesh.getElements(dim, tag)
			for elem_type, elems, nodes in zip(element_types, element_tags, node_tags):
				if elem_type == 2:  # '2' corresponds to triangular elements
					# Gmsh returns a flat list of node indices for triangles
					for i in range(0, len(nodes), 3):
						faces.append((nodes[i], nodes[i+1], nodes[i+2]))
		# Check
		print(rf"  Vertices remeshed from {self.nvertex} to {np.array(vertices).shape[0]}")
		print(rf"  Faces remeshed from {self.nfaces} to {np.array(faces).shape[0]}")
		# Store new mesh
		self.vertex = np.array(vertices)
		self.faces  = np.array(faces)-1
		self.nvertex = self.vertex.shape[0]
		self.nfaces = self.faces.shape[0]



# LINESEG Class
class LINESEG():
	"""Class for linesegment

	"""
	def __init__(self,R=None,phi=None,Z=None,RHat=None,ZHat=None,L=None):
		from numpy import sqrt
		self.Phi = phi
		self.R   = R
		self.Z   = Z
		self.RHat = RHat
		self.ZHat = ZHat
		self.L    = L
		if (RHat and ZHat):
			norm = sqrt(RHat**2+ZHat**2)
			self.RHat = RHat/norm
			self.ZHat = ZHat/norm

	def getEndpoints(self):
		"""Return cartesian coordinates of the LINESEG

		This routine returns the cartesian coordinates of the LINESEG
		end points.

		Returns
		----------
		points : list
			The [2,x,y,z] coordiantes of the three vertices. [m]
		"""
		from numpy import sin,cos
		t1 = [self.R-0.5*self.L*self.ZHat,self.Z+0.5*self.L*self.RHat]
		t2 = [self.R+0.5*self.L*self.ZHat,self.Z-0.5*self.L*self.RHat]
		p1 = [t1[0]*cos(self.Phi),t1[0]*sin(self.Phi),t1[1]]
		p2 = [t2[0]*cos(self.Phi),t2[0]*sin(self.Phi),t2[1]]
		return [p1,p2]

# WEDGE Class
class WEDGE(LINESEG):
	"""Class for wedge

	"""
	def __init__(self,R=None,phi=None,Z=None,RHat=None,ZHat=None,L=None,alpha=None):
		from numpy import sin,cos
		super().__init__(R,phi,Z,RHat,ZHat,L)
		self.alpha = alpha
		if alpha:
			self.sinah = sin(alpha*0.5)
			self.cosah = cos(alpha*0.5)

	def getEndpoints(self):
		"""Return cartesian coordinates of the WEDGE

		This routine returns the cartesian coordinates of the WEDGE
		end points.

		Returns
		----------
		points : list
			The [3,x,y,z] coordiantes of the three vertices. [m]
		"""
		from numpy import sin,cos
		t1 = [self.R-self.L*self.sinah*self.ZHat,self.Z+self.L*self.sinah*self.RHat]
		t2 = [self.R+self.L*self.cosah*self.RHat,self.Z+self.L*self.cosah*self.ZHat]
		t3 = [self.R+self.L*self.sinah*self.ZHat,self.Z-self.L*self.sinah*self.RHat]
		p1 = [t1[0]*cos(self.Phi),t1[0]*sin(self.Phi),t1[1]]
		p2 = [t2[0]*cos(self.Phi),t2[0]*sin(self.Phi),t2[1]]
		p3 = [t3[0]*cos(self.Phi),t3[0]*sin(self.Phi),t3[1]]
		return [p1,p2,p3]

# CIRCLE Class
class CIRCLE(LINESEG):
	"""Class for wedge

	"""
	def __init__(self,R=None,phi=None,Z=None,RHat=None,ZHat=None,L=None,N=None):
		from numpy import sin,cos
		super().__init__(R,phi,Z,RHat,ZHat,L)
		self.N = N

	def getEndpoints(self):
		"""Return cartesian coordinates of the CIRCLE

		This routine returns the cartesian coordinates of the CIRCLE
		end points. Note the circle does not close.

		Returns
		----------
		outarr : list
			The [n,x,y,z] coordiantes of the vertices. [m]
		"""
		import numpy as np
		outarr = []
		for k in range(self.N):
			t1 = [self.R+self.L*cos(2*np.pi*k/self.N)*self.RHat,self.Z+self.L*sin(2*np.pi*k/self.N)*self.ZHat]
			p1 = [t1[0]*cos(self.Phi),t1[0]*sin(self.Phi),t1[1]]
			temp.append(t1)
		return outarr

# Parameterized wall model
class PARAM_WALL():
	"""Class for defining parameterized walls

	"""
	def __init__(self):
		self.elements = None
		pass

	def addElement(self,subset):
		"""Add a subset to the elements

		This routine appends a list of primative shapes to the
		elements list

		Parameters
		----------
		subset : list
			A list of primatives
		"""
		# Sort by phi
		subset.sort(key=lambda x: x.Phi, reverse=False)
		# Add to elements
		if type(self.elements) is type(None):
			self.elements=[subset]
		else:
			self.elements.extend([subset])

	def getWall(self):
		"""Return a wall using the elements

		This routine computes the vertices and faces and returns a wall
		object based on the elements. Each item of elements is a list.
		This subset list is composed of a group of similar shapes which
		are toroidally linked. The notion being that these elements
		define a toroidal shape. Thus each item of the subset must be
		the same shape.

		Returns
		----------
		wall : WALL class
			Returns a wall object.
		"""
		import numpy as np
		from datetime import datetime
		out_wall = WALL()
		out_wall.faces=[]
		out_wall.vertex=[]
		nvertex = 0
		nfaces  = 0
		k       = 0
		for subset in self.elements:
			# Append the vertex information
			for item in subset:
				p=item.getEndpoints()
				for ps in p:
					out_wall.vertex.append(ps)
			# Determine how many points are in the shape
			npoints = len(subset[0].getEndpoints())
			# Determine number of toroidal elements
			nelements = len(subset)
			# Append the face information
			for i in range(nelements-1):
				for j in range(npoints-1):
					faces = [k, k+1, k+npoints]
					out_wall.faces.append(faces)
					faces = [k+1, k+npoints+1, k+npoints]
					out_wall.faces.append(faces)
					k = k + 1
				# Close the circle
				if npoints > 3:
					faces = [k-1, k-npoints, k+npoints-1]
					out_wall.faces.append(faces)
					faces = [k-npoints, k, k+npoints-1]
				k = k + 1 # skip endpoint
			k = k + npoints # skip to next subset
		# Setup wall object
		out_wall.nvertex = len(out_wall.vertex)
		out_wall.nfaces  = len(out_wall.faces)
		out_wall.vertex  = np.array(out_wall.vertex)
		out_wall.faces   = np.array(out_wall.faces)
		out_wall.name = f"Generated using simpilfied wall elements in Python."
		out_wall.date = datetime.today().strftime('%Y-%m-%d')
		return out_wall

# SOLID WALL MODEL
class SOLIDWALL():
	"""Class for defining parameterized solid walls like in parastell

	"""
	def __init__(self,radial_build_dict,poloidal_angles,toroidal_angles):
		self.radial_build_dict = radial_build_dict
		self.poloidal_angles = poloidal_angles
		self.toroidal_angles = toroidal_angles
		self.nfp = int(360.0/max(toroidal_angles))
		self.spline_order = 3

	def setInitialSurface(self,r,z,nr,nz):
		"""Define the inital R,Z surface

		This routine is used to define the initial surface in terms of
		R, Z, n_R, and n_Z where each quantity has the dimension of the
		poloidal and toroidal angles.

		Parameters
		----------
		r : array
			Cylindrical R values of surface (npoloidal,ntoroidal)
		z : array
			Cylindrical Z values of surface (npoloidal,ntoroidal)
		nr : array
			Cylindrical R normal values of surface (npoloidal,ntoroidal)
		nz : array
			Cylindrical Z normal values of surface (npoloidal,ntoroidal)
		"""
		self.r = r
		self.z = z
		self.nr = nr
		self.nz = nz

	def plotWalls2D(self,phi=0.0,ax=None):
		"""Plot 2D cuts of the wall

		This routine is used to define the initial surface in terms of
		R, Z, n_R, and n_Z where each quantity has the dimension of the
		poloidal and toroidal angles.

		Parameters
		----------
		phi : float
			Toroidal angle to plot
		"""
		import numpy as np
		import matplotlib.pyplot as plt
		from scipy.interpolate import make_interp_spline

		lplotnow = False
		if not ax:
			ax = plt.axes()
			lplotnow = True
		# Helpers
		nradial = len(self.radial_build_dict)+1
		ntheta  = len(self.poloidal_angles)
		nphi    = len(self.toroidal_angles)
		ntheta_out   = 90
		theta_plot   = np.linspace(0.0,360.0,ntheta_out)
		# Generate the shells
		total_thickness = np.zeros((ntheta,nphi))
		R_shells        = np.zeros((nradial,ntheta,nphi))
		Z_shells        = np.zeros((nradial,ntheta,nphi))
		R_shells[0,:,:] = self.r
		Z_shells[0,:,:] = self.z
		i = 1
		for name,properties in self.radial_build_dict.items():
			print(f'Working on: {name}')
			thick = properties['thickness_matrix']
			total_thickness = total_thickness + thick
			R_shells[i,:,:] = self.r + total_thickness*self.nr
			Z_shells[i,:,:] = self.z + total_thickness*self.nz
			i = i + 1
		# Force periodic
		R_shells[:,-1,:] = R_shells[:,0,:]
		R_shells[:,:,-1] = R_shells[:,:,0]
		Z_shells[:,-1,:] = Z_shells[:,0,:]
		Z_shells[:,:,-1] = Z_shells[:,:,0]
		Rspl=make_interp_spline(np.squeeze(self.toroidal_angles), R_shells, k=self.spline_order, bc_type='periodic',axis=2)
		Zspl=make_interp_spline(np.squeeze(self.toroidal_angles), Z_shells, k=self.spline_order, bc_type='periodic',axis=2)
		Rtmp = Rspl(phi)
		Ztmp = Zspl(phi)
		Rspl=make_interp_spline(np.squeeze(self.poloidal_angles), Rtmp, k=self.spline_order, bc_type='periodic',axis=1)
		Zspl=make_interp_spline(np.squeeze(self.poloidal_angles), Ztmp, k=self.spline_order, bc_type='periodic',axis=1)
		Rplt = Rspl(theta_plot)
		Zplt = Zspl(theta_plot)
		names = list(self.radial_build_dict.keys())
		for i in range(nradial-1,0,-1):
			ax.fill(Rplt[i,:],Zplt[i,:],label=names[i-1])
		ax.fill(Rplt[0,:],Zplt[0,:],'r',label='plasma')
		ax.legend()
		ax.axis('equal')
		ax.set_xlabel('R [m]')
		ax.set_ylabel('Z [m]')
		ax.set_title(rf'Wall Layers ($\phi={np.rad2deg(phi)}$)')
		if lplotnow: plt.show()

	def generateWalls(self,npol=360,ntor=90):
		"""Creates the solid walls

		The routine creates the solid walls with an output grid
		resolution equal to that of the npol and ntor values.
		This routine produces both wall.dat files and wall.stl
		files for each layer.

		Parameters
		----------
		npol : int
			Number of poloidal points to use for wall generation
		ntor : int
			Number of toroidla points to use for wall generation
		"""
		import numpy as np
		from scipy.interpolate import make_interp_spline
		from copy import deepcopy
		from datetime import datetime
		# Define output arrays
		poloidal_angles_out = np.linspace(0.0,360.0,npol)
		toroidal_angles_out = np.linspace(0.0,max(self.toroidal_angles),ntor)
		# Surface helpers
		nradial = len(self.radial_build_dict)+1
		ntheta  = len(self.poloidal_angles)
		nphi    = len(self.toroidal_angles)
		# Construct layers
		R_shells        = np.zeros((nradial,ntheta,nphi))
		Z_shells        = np.zeros((nradial,ntheta,nphi))
		total_thickness = np.zeros((ntheta,nphi))
		R_shells[0,:,:] = self.r
		Z_shells[0,:,:] = self.z
		i = 1
		for name,properties in self.radial_build_dict.items():
			print(f'Working on: {name}')
			thick = properties['thickness_matrix']
			total_thickness = total_thickness + thick
			R_shells[i,:,:] = R_shells[0,:,:] + total_thickness*self.nr
			Z_shells[i,:,:] = Z_shells[0,:,:] + total_thickness*self.nz
			i = i + 1
		# Force periodic
		R_shells[:,-1,:] = R_shells[:,0,:]
		R_shells[:,:,-1] = R_shells[:,:,0]
		Z_shells[:,-1,:] = Z_shells[:,0,:]
		Z_shells[:,:,-1] = Z_shells[:,:,0]
		# Spline over the poloidal direction
		splr = make_interp_spline(np.squeeze(self.poloidal_angles), R_shells, k=self.spline_order, bc_type='periodic',axis=1)
		splz = make_interp_spline(np.squeeze(self.poloidal_angles), Z_shells, k=self.spline_order, bc_type='periodic',axis=1)
		Rp_shells = splr(poloidal_angles_out)
		Zp_shells = splz(poloidal_angles_out)
		# Spline over the toroidal direction
		splr = make_interp_spline(np.squeeze(self.toroidal_angles), Rp_shells, k=self.spline_order, bc_type='periodic',axis=2)
		splz = make_interp_spline(np.squeeze(self.toroidal_angles), Zp_shells, k=self.spline_order, bc_type='periodic',axis=2)
		Rpt_shells = splr(toroidal_angles_out)
		Zpt_shells = splz(toroidal_angles_out)
		# Create a geometry helper for completing the torus
		phirad = np.tile(np.deg2rad(toroidal_angles_out),(nradial,npol,1))
		xtarr  = Rpt_shells*np.cos(phirad)
		ytarr  = Rpt_shells*np.sin(phirad)
		ztarr  = Zpt_shells
		# Extend to full torus
		xarr = xtarr; yarr = ytarr; zarr=ztarr
		for v in range(1,self.nfp):
			cop = np.cos(v*2.0*np.pi/self.nfp)
			sip = np.sin(v*2.0*np.pi/self.nfp)
			xarr = np.append(xarr,xtarr*cop-ytarr*sip,2)
			yarr = np.append(yarr,ytarr*cop+xtarr*sip,2)
			zarr = np.append(zarr,ztarr,2)
		# Adjust ordering (s,th,ph) -> (s,ph,th)
		xarr = np.swapaxes(xarr,1,2)
		yarr = np.swapaxes(yarr,1,2)
		zarr = np.swapaxes(zarr,1,2)
		# Now generate the shells as walls
		shells=[]
		for k in range(nradial):
			shells.append(WALL())
			shells[k].genWallfromOrdered(np.array([xarr[k,:,:],yarr[k,:,:],zarr[k,:,:]]))
			shells[k].wallClean()
		# Now generate the soild wall from the shell models
		#names = list(self.radial_build_dict.keys())
		#for k in range(1,nradial):
		#	wall_out = WALL()
		#	wall_out = deepcopy(shells[k-1])
		#	wall_out.faces = shells[k-1].faces[:,[2,1,0]] # flips normal direction for inner layer
		#	wall_outer = deepcopy(shells[k])
		#	wall_out.wallAdd(wall_outer)
		#	wall_out.name = names[k-1] + "Generated by pySTEL SOLIDWALL"
		#	wall_out.date = datetime.today().strftime('%Y-%m-%d')
		#	wall_out.write_wall(names[k-1]+'.dat')
		#	wall_out.writeSTL(names[k-1]+'.stl')
		# For now just output the inner surface
		names = list(self.radial_build_dict.keys())
		for k in range(1,nradial):
			wall_out = WALL()
			wall_out = deepcopy(shells[k-1])
			wall_out.faces = shells[k-1].faces[:,[2,1,0]] # flips normal direction for inner layer
			#wall_outer = deepcopy(shells[k])
			#wall_out.wallAdd(wall_outer)
			wall_out.name = names[k-1] + "Generated by pySTEL SOLIDWALL"
			wall_out.date = datetime.today().strftime('%Y-%m-%d')
			wall_out.write_wall(names[k-1]+'.dat')
			wall_out.writeSTL(names[k-1]+'.stl')
		# Do final outer surface
		wall_out = WALL()
		wall_out = deepcopy(shells[nradial-1])
		wall_out.faces = shells[nradial-1].faces[:,[2,1,0]] # flips normal direction for inner layer
		wall_out.name = "Outermost surface Generated by pySTEL SOLIDWALL"
		wall_out.date = datetime.today().strftime('%Y-%m-%d')
		wall_out.write_wall('outermost_surface.dat')
		wall_out.writeSTL('outermost_surface.stl')

if __name__=="__main__":
	import sys
	sys.exit(0)
