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

	def writeSTL(self,filename='wall.stl'):
		"""Outputs an STL object from wall object

		This routine generates an stereolithography (STL) file from
		the wall object.

		Parameters
		----------
		filename : str (optional)
			Filename for output file
		"""
		import numpy as np
		from stl import mesh
		stlobj = mesh.Mesh(np.zeros(self.faces.shape[0],dtype=mesh.Mesh.dtype))
		for i, f in enumerate(self.faces):
			for j in range(3):
				stlobj.vectors[i][j] = self.vertex[f[j],:]
		stlobj.save(filename)

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













if __name__=="__main__":
	import sys
	sys.exit(0)
