#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for reading and handling MUMATERIALS
data.
"""

# Libraries
from libstell.libstell import LIBSTELL

# Constants

# MUMAT Class
class MUMAT():
	"""Class for working with MUMATERIAL data

	"""
	def __init__(self, parent=None):
		super().__init__()
		self.libStell = LIBSTELL()

	def read_mumat_file(self,filename):
		"""Reads a MUMATERIAL tetrahedron file

		This routine reads the MUMATERIAL tetrahedron file.

		Parameters
		----------
		file : str
			Path to tetrahedron file.
		"""
		import numpy as np
		import copy
		mumat_dict = copy.deepcopy(self.libStell.read_mumat_file(filename))
		for key in mumat_dict:
			setattr(self, key, mumat_dict[key])


# MUMAT Input Class
class MUMAT_INPUT():
	"""Class for working with MUMAT INPUT data

	"""
	def __init__(self, parent=None):
		self.libStell = LIBSTELL()

	def read_input(self,filename):
		"""Reads MUMAT_INPUT namelist from a file

		This routine wrappers the mumaterial_mod module reading routine.
		Parameters
		----------
		filename : string
			Input file name with MUMAT_INPUT namelist
		"""
		mumat_dict = self.libStell.read_mumat_input(filename)
		for key in mumat_dict:
			setattr(self, key, mumat_dict[key])

	def write_input(self,filename):
		"""Writes MUMAT_INPUT namelist to a file

		This routine wrappers the mumaterial_mod module writing routine.
		Parameters
		----------
		filename : string
			Input file name to write MUMAT_INPUT namelist to
		"""
		out_dict = vars(self)
		self.libStell.write_mumat_input(filename,out_dict)

# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)









