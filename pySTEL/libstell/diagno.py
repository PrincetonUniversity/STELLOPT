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

	def read_fluxloops(self,filename,scale_factor=1.0):
		"""Reads a DIAGNO fluxloop definition file.

		This routine reads the diagno fluxloop deffinition file.

		Parameters
		----------
		filename : str
			Path to diagno_fluxloops file.
		scale_factor : float
			Scale factor to apply to data (default = 1.0)

		Returns
		-------
		floops : dict
			Dictionary of fluxloop names and array of x,y,z data.
		"""
		import numpy as np
		f = open(filename,'r')
		line = f.readline()
		nloops = int(line)
		floops = {}
		for i in range(nloops):
			line = f.readline()
			[txt1,txt2,txt3,txt4] = line.split()
			nels = int(txt1)
			name = txt4
			xyz = np.ndarray((3,nels))
			for j in range(loops):
				xyz[:,j]=[float(x) for x in line.split()]
			floops[name] = xyz
		f.close()
		return floops

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
		line = f.readline()
		flux = [float(x) for x in line.split()]
		for line in f:
			names = names.append(line)
		f.close()
		return names,np.array(flux)



# BEASM3D Input Class
class DIAGNO_IN():
	"""Class for working with DIAGNO_IN data

	"""
	def __init__(self, parent=None):
		self.libStell = LIBSTELL()

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