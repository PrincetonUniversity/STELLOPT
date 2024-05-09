##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for calculating
quantities related to nuclear fusion
"""

# Libraries
#from libstell import libstell

# Constants
EC = 1.602176634E-19 # Electron charge

# VMEC Class
class FUSION():
	"""Class for working with FUSION data

	"""
	def __init__(self):
		self.C_DICT = { "DT"    : [ 1.17302E-09,  1.51361E-02,  7.51886E-02, \
						 	        4.60643E-03,  1.35000E-02, -1.06750E-04, \
							        1.36600E-05] , \
						"DDT"   : [ 5.65718E-12,  3.41267E-03,  1.99167E-03, \
								    0.00000E+00,  1.05060E-05,  0.00000E+00, \
								    0.00000E+00] , \
						"DDHe3" : [ 5.43360E-12,  5.85778E-03,  7.68222E-03, \
								    0.00000E+00, -2.96400E-06,  0.00000E+00, \
								    0.00000E+00] , \
						"DHe3"  : [ 5.51036E-10,  6.41918E-03, -2.02896E-03, \
								   -1.91080E-05,  1.35776E-04,  0.00000E+00, \
								    0.00000E+00] }
		self.MRC2_DICT = { "DT"    : 1124656, \
						   "DDT"   :  937814 ,\
						   "DDHe3" :  937814 ,\
						   "DHe3"  : 1124572 }
		self.BG_DICT = { "DT"    : 34.3827, \
						 "DDT"   : 31.3970 ,\
						 "DDHe3" : 31.3970 ,\
						 "DHe3"  : 68.7508 }
		self.E_DT_He    = 3.52E6*EC
		self.E_DD_T     = 1.01E6*EC
		self.E_DD_p     = 3.02E6*EC
		self.E_DD_He3   = 0.82E6*EC
		self.E_DHe3_He4 = 3.60E6*EC
		self.E_DHe3_p   = 1.47E7*EC

	def sigmaBH(self, ti, reaction='DT'):
		"""Computs the Bosch Hale Fusion cross sections

		This routine calculates the fusion cross sections given
		density and temperature using the Bosch Hale model.
		H.-S. Bosch and G. M. Hale 1992 Nucl. Fusion 32 611
		https://doi.org/10.1088/0029-5515/32/4/I07

		Parameters
		----------
		ti : real
			Ion Temperature in [eV]
		reaction : string
			Reaction type 'DT','DDT','DDHe3','DHe3'
		Returns
		----------
		sigma : real
			Reaction rate [m^3/s]
		"""
		import math
		C = self.C_DICT[reaction]
		BG = self.BG_DICT[reaction]
		MRC2 = self.MRC2_DICT[reaction]
		ti_kev = ti * 1E-3
		zeta   = ( ( ( C[5] * ti_kev ) + C[3] ) * ti_kev + C[1] ) * ti_kev
		zeta   = zeta / ( ( ( ( C[6] * ti_kev ) + C[4] ) * ti_kev + C[2] ) * ti_kev + 1.0 )
		zeta   = 1.0 - zeta
		theta  = ti_kev / zeta
		eta    = ( 0.25 * BG * BG / theta ) ** (1.0/3.0)
		return 1.0E-6 * C[0] * theta * math.sqrt( eta / ( MRC2 * ti_kev * ti_kev * ti_kev ) ) * math.exp( -3 * eta )


if __name__=="__main__":
	import sys
	temp = FUSION()
	print(temp.sigmaBH(10000,'DDT') * 1E20 * 1E20)
	#temp = FIELDLINES()
	#val=temp.read_indata('input.ORBITS')
	#temp.read_wout('fieldlines_test.nc')
	sys.exit(0)
