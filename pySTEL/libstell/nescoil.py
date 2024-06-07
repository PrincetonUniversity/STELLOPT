#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for reading and handling NESCOIL
surface potential data.
"""

# Libraries
from libstell.libstell import LIBSTELL, FourierRep

# Constants

# NESCOIL Class
class NESCOIL(FourierRep):
	"""Class for working with NESCOIL output data.

	"""
	def __init__(self):
		super().__init__()
		self.libStell = LIBSTELL()

	def read_nescout(self,filename):
		"""Reads a NESCOIL output file.

		This routine reads and initilizes the NESCOIL
		class with information from a nescout file.

		Parameters
		----------
		file : str
			Path to nescout file.
		"""
		nescout_dict = self.libStell.read_nescout(filename)
		for key in nescout_dict:
			setattr(self, key, nescout_dict[key])



# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)