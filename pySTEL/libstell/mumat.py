#!/usr/bin/env python3
"""
This library provides a python class for reading and handling 
MUMATERIAL magnetic material data.
"""

# Libraries
from libstell.libstell import LIBSTELL
import numpy as np

# Constants

# MUMAT Class
class MUMAT():
	"""Class for working with MUMAT data

	"""
	def __init__(self):
		self.machine_string=''
		self.date=''
		self.ntet=-1
		self.nvertex=-1
		self.nstate=0
		self.vertex=None
		self.tet=None
		self.state_dex = np.array([],dtype=int)
		self.state_type = np.array([],dtype=int)
		self.constant_mu = np.array([])
		self.constant_mu_o = np.array([])
		self.Mrem = np.array([])
		self.state_function_H = []
		self.state_function_M = []

	def read_mumat(self,filename):
		"""Reads a MUMAT data file.

		This routine reads the mumat data file.

		Parameters
		----------
		filename : str
			Path to mumat file.

		"""
		import numpy as np
		f = open(filename,'r')
		self.machine_string = f.readline()
		self.date = f.readline()
		(self.nvertex, self.ntet, self.nstate) = [int(x) for x in f.readline().split()]
		self.vertex = np.zeros((3,self.nvertex))
		self.tet = np.zeros((4,self.ntet),dtype=int)
		self.state_dex = np.zeros((self.ntet),dtype=int)
		self.state_type = np.zeros((self.nstate),dtype=int)
		self.constant_mu = np.zeros((self.nstate))
		self.constant_mu_o = np.zeros((self.nstate))
		self.Mrem = np.zeros((3,self.nstate))
		self.state_function_H = [None] * self.nstate
		self.state_function_M = [None] * self.nstate
		for i in range(self.nvertex):
			self.vertex[:,i] = [float(x) for x in f.readline().split()]
		for i in range(self.ntet):
			(n0,n1,n2,n3,n4) = [int(x) for x in f.readline().split()]
			self.tet[:,i] = [n0,n1,n2,n3]
			self.state_dex[i] = n4
		for i in range(self.nstate):
			self.state_type[i] = int(f.readline())
			if self.state_type[i] == 1:
				(self.constant_mu[i],self.constant_mu_o[i]) = [float(x) for x in f.readline().split()]
				self.Mrem[:,i] = [float(x) for x in f.readline().split()]
			elif self.state_type[i] == 2:
				nMH = int(f.readline())
				self.state_function_H[i] = [float(x) for x in f.readline().split()]
				self.state_function_M[i] = [float(x) for x in f.readline().split()]
			elif self.state_type[i] ==3:
				self.constant_mu[i] = float(f.readline())
			else:
				print(f'!!! Unknown state_type == {self.state_type[i]}')
		f.close()

	def write_mumat(self,filename):
		"""Writes a MUMAT data file.

		This routine writes the mumat data file.

		Parameters
		----------
		filename : str
			Path to mumat file.

		"""
		f = open(filename,'w')
		f.write(f'{self.machine_string}\n')
		f.write(f'{self.date}\n')
		f.write(f'{self.nvertex:d} {self.ntet:d} {self.nstate:d}\n')
		for i in range(self.nvertex):
			f.write(f'{self.vertex[0,i]:.10E} {self.vertex[1,i]:.10E} {self.vertex[2,i]:.10E}\n')
		for i in range(self.ntet):
			f.write(f'{self.tet[0,i]:d} {self.tet[1,i]:d} {self.tet[2,i]:d} {self.tet[3,i]:d} {self.state_dex[i]:d}\n')
		for i in range(self.nstate):
			f.write(f'{self.state_type[i]:d}\n')
			if self.state_type[i] == 1:
				f.write(f'{self.constant_mu[i]:.10E} {self.constant_mu_o[i]:.10E}\n')
				f.write(f'{self.Mrem[0,i]:.10E} {self.Mrem[1,i]:.10E} {self.Mrem[2,i]:.10E}\n')
			elif self.state_type[i] == 2:
				nMH = len(self.state_function_H[i])
				f.write(f'{nMH:d}\n')
				for j in range(nMH):
					f.write(f'{self.state_function_H[i][j]:.10E} ')
				f.write('\n')
				for j in range(nMH):
					f.write(f'{self.state_function_M[i][j]:.10E} ')
				f.write('\n')
			elif self.state_type[i] ==3:
				f.write(f'{self.constant_mu[i]}\n')
		f.close()

	def load_gmsh(self,gmsh):
		"""Loads a mumat data structure from gmsh object

		This routine loads the mumat geometry data from a gmsh
		object. The vertices are loaded from gmsh.model.mesh.getNodes()
		and the 4-node tetrahedrons are loaded from
		gmsh.model.mesh.getElements(dim=3). Other element types are
		ignored, but an error message is printed. The state function
		is set to 0 for all tets.  This needs to be set seperately by
		invoking add_state and set_state.  The example usage of this
		routine is as follows:

			import gmsh
			gmsh.initialize()
			gmsh.open(file)
			gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 5.0)
			gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 50.000)
			gmsh.model.geo.synchronize()
			gmsh.model.mesh.generate(3)
			mumat_data.load_gmsh(gmsh)

		Parameters
		----------
		gmsh : GMSH object
			The main GMSH object
		"""
		from datetime import datetime
		if not gmsh.isInitialized():
			print('GMSH object not initialized!')
			return
		(a,b,c) = gmsh.model.mesh.getNodes()
		self.machine_string = f'GMSH Mesh of {gmsh.model.get_file_name()}'
		self.date = datetime.today().strftime('%Y-%m-%d')
		self.ntet = -1
		self.vertex= b.reshape((-1,3)).T
		self.nvertex = self.vertex.shape[1]
		elemTypes, elemTags, nodeTags = gmsh.model.mesh.getElements(dim=3)
		tets = None
		for i in range(len(elemTypes)):
			if elemTypes[i] == 4: # Standard 4-node tetrahedron
				if type(tets) == type(None):
					tets = nodeTags[i]
				else:
					tets = np.append(tets,nodeTags[i])
			else:
				print(f'Found non-tetrahedron(4-node) elemType == {lemTypes[i]:d} skipping.')
		self.tet = tets.reshape((-1,4)).T
		self.ntet = self.tet.shape[1]
		self.state_dex = np.zeros((self.ntet),dtype=int)

	def add_state(self,state_type,mu=None,mu_o=None,H=None,M=None,Mrem=None):
		"""Adds a state function to the mumat object

		This routine adds a state function to the mumat object.
		Currently three state_types are supported:
			state_type = 1
				mu:   Constant permeability
				mu_o: Constant permeability 0
				Mrem: Remnant magnetization (vector)
			state_type = 2
				H:    H function knots [A/m]
				M:    M function values [T]
			state_type = 3
				mu:   Constant permeability

		Parameters
		----------
		state_type : int
			State function type (1,2,3)
		mu : float
			Magnetic permeability
		mu_o : float
			Magnetic permeability
		Mrem : numpy array
			Magnetization vector
		H : numpy array
			Magnetic field strength [A/m]
		M : numpy array
			Magnetic Flux array [T]
		"""
		if state_type <= 0 or state_type > 3:
			print(f'!!! Unknown state_type == {self.state_type[i]}')
			return
		self.nstate = self.nstate + 1
		self.state_type = np.append(self.state_type,state_type)
		if state_type == 1:
			self.constant_mu = np.append(self.constant_mu,mu)
			self.constant_mu_o = np.append(self.constant_mu_o,mu_o)
			self.Mrem = np.append(self.Mrem,Mrem)
		elif state_type == 2:
			self.state_function_H.append(H)
			self.state_function_M.append(M)
		elif state_type == 3:
			self.constant_mu = np.append(self.constant_mu,mu)

	def set_state(self,state_dex,tet_dex):
		"""Set the state function for a given set of tetrahedrons

		This subroutine sets the state function for a given set of
		tetrahedrons. The state index and an array of tetrahedron
		indices are provided by the user. These may both be
		a single value. Tet_dex may be a list of python indices
		and state_dex a single value. Tet_dex and state_dex may
		both be lists of the same size.

		Parameters
		----------
		state_dex : int or int-list
			Index of state function array to use
		tet_dex : int or int-list
			Index of tetrahedron to set
		"""
		self.state_dex[tet_dex] = state_dex

	def plot_mesh(self,plot3D=None):
		"""Plots a mumat mesh in 3D

		This routine plots the mumaterial mesh in 3D

		Parameters
		----------
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
		[points, tetra]=plt.tetrameshTo3DTetra(self.vertex.T,self.tet.T-1)
		# Generate Wall colors
		plt.add3Dwireframe(points,tetra,color='red')
		plt.setBGcolor(1.0,1.0,1.0)
		# Render if requested
		if lplotnow: plt.render()

	def plot_state(self,state_dex=None,ax=None):
		"""Plots a mumat state function

		This routine plots the mumaterial state functions. If given a
		state_dex it plots that state function. Otherwise it will plot
		all state functions.

		Parameters
		----------
		state_dex : int
			State to plot (default: plots all)
		ax : axes (optional)
			Matplotlib axes object to plot to.
		"""
		import matplotlib.pyplot as pyplot
		lplotnow = False
		# Handles axes
		if not ax: lplotnow = True
		if type(state_dex) == type(None):
			states = list(range(self.nstate))
		else:
			states = state_dex
		for i in states:
			if self.state_type==2:
				if lplotnow: ax = pyplot.axes()
				ax.plot(self.state_function_H[i],self.state_function_M[i],color='black',linewidth=2.0)
				ax.set_xlabel('H [A/m]')
				ax.set_ylabel('M [T]')
				ax.set_title(f'Mumaterial State Function ({i:02d})')
				if lplotnow: pyplot.show()

# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)
