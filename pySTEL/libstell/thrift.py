#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for reading and handling 
BEAMS3D Energetic particle data.
"""

# Libraries
from libstell.libstell import LIBSTELL

# Constants

# BEAMS3D Class
class THRIFT():
	"""Class for working with THRIFT data

	"""
	def __init__(self):
		test = 0

	def read_thrift(self,file):
		"""Reads a THRIFT HDF5 file

		This routine reads and initilizes the THRIFT
		class with variable information from an HDF5 file.

		Parameters
		----------
		file : str
		Path to HDF5 file.
		"""
		import h5py
		import numpy as np
		beams_data = {}
		#read file
		with h5py.File(file,'r') as f:
			# Logicals
			for temp in ['leccd','lohmic','lvmec','lnbcd']:
				if temp in f:
					setattr(self, temp, np.int64(f[temp][0])==1)
			# Integers
			for temp in ['npicard','nrho','nssize','ntimesteps']:
				if temp in f:
					setattr(self, temp, np.int64(f[temp][0]))
			# Floats
			for temp in ['VERSION','jtol','picard_factor']:
				if temp in f:
					setattr(self, temp, float(f[temp][0]))
			# Arrays
			for temp in ['THRIFT_ALPHA1','THRIFT_ALPHA2','THRIFT_ALPHA3','THRIFT_ALPHA4','THRIFT_AMINOR',\
				'THRIFT_BAV','THRIFT_BSQAV','THRIFT_BVAV','THRIFT_COEFF_A','THRIFT_COEFF_B','THRIFT_COEFF_BP',\
				'THRIFT_COEFF_C','THRIFT_COEFF_CP','THRIFT_COEFF_D','THRIFT_COEFF_DP','THRIFT_ETAPARA','THRIFT_I',\
				'THRIFT_IBOOT','THRIFT_IECCD','THRIFT_INBCD','THRIFT_IOHMIC','THRIFT_IOTA','THRIFT_IPLASMA','THRIFT_ISOURCE',\
				'THRIFT_J','THRIFT_JBOOT','THRIFT_JECCD','THRIFT_JNBCD','THRIFT_JOHMIC','THRIFT_JPLASMA','THRIFT_JSOURCE',\
				'THRIFT_MATLD','THRIFT_MATMD','THRIFT_MATRHS','THRIFT_MATUD','THRIFT_P','THRIFT_PHIEDGE','THRIFT_PPRIME',\
				'THRIFT_RHO','THRIFT_RHOFULL','THRIFT_RMAJOR','THRIFT_S','THRIFT_S11','THRIFT_S12','THRIFT_SNOB','THRIFT_T','\
				THRIFT_UGRID','THRIFT_VP']:
				if temp in f:
					setattr(self, temp, np.array(f[temp][:]))
		# Make derived arrays
		return

# THRIFT Input Class
class THRIFT_INPUT():
	"""Class for working with THRIFT INPUT data

	"""
	def __init__(self, parent=None):
		self.libStell = LIBSTELL()

	def read_input(self,filename):
		"""Reads THRIFT_INPUT namelist from a file

		This routine wrappers the thrift_input_mod module reading routine.
		Parameters
		----------
		filename : string
			Input file name with THRIFT_INPUT namelist
		"""
		# Not yet implemented
		#indata_dict = self.libStell.read_thrift_input(filename)
		#for key in indata_dict:
		#	setattr(self, key, indata_dict[key])

	def write_input(self,filename):
		"""Writes THRIFT_INPUT namelist to a file

		This routine wrappers the thrift_input_mod module writing routine.
		Parameters
		----------
		filename : string
			Input file name to write THRIFT_INPUT namelist to
		"""
		# Not yet implemented
		#out_dict = vars(self)
		#self.libStell.write_thrift_input(filename,out_dict)

# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)
