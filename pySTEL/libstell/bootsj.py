##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for reading and handling 
BOOTSJ data.
"""

# Libraries
#from libstell import libstell

# Constants

# BOOTSJ Class
class BOOTSJ():
	"""Class for working with BOOTSJ data

	"""
	def __init__(self):
		test = 0

	def read_jBbs(self,filename):
		"""Reads a BOOTSJ jBbs file

		This routine reads a BOOTSJ jBbs file.

		Parameters
		----------
		file : str
			Path to jBbs file.
		"""
		import numpy as np
		f = open(filename,'r')
		lines = f.readlines()
		f.close()
		#check header
		if  'surface' not in lines[0]:
			print("Bad Synatx line 1 in jBbs file")
		npts = int(len(lines))
		self.jBbs_s = np.zeros((npts))
		self.jBbs_JbsdotB = np.zeros((npts))
		self.jBbs_tokfrac = np.zeros((npts))
		for i in range(npts):
			line = lines[i+1].split()
			self.jBbs_s[i]       = float(line[0])
			self.jBbs_JbsdotB[i] = float(line[1])
			self.jBbs_tokfrac[i] = float(line[2])

	def read_answers_plot(self,filename):
		"""Reads a BOOTSJ answers_plot file

		This routine reads a BOOTSJ answers_plot file.

		Parameters
		----------
		file : str
			Path to answers_plot file.
		"""
		import numpy as np
		f = open(filename,'r')
		lines = f.readlines()
		f.close()
		# No header
		npts = int(len(lines))
		self.rhoar = np.zeros((npts))
		self.gbsnorm = np.zeros((npts))
		self.amain = np.zeros((npts))
		self.aiterm1 = np.zeros((npts))
		self.other1 = np.zeros((npts))
		self.dibs = np.zeros((npts))
		self.bsdense = np.zeros((npts))
		self.bsdensi = np.zeros((npts))
		self.bstempe = np.zeros((npts))
		self.bstempi = np.zeros((npts))
		self.qsafety = np.zeros((npts))
		self.ftrapped = np.zeros((npts))
		self.bsnorm = np.zeros((npts))
		self.tempe1 = np.zeros((npts))
		self.tempi1 = np.zeros((npts))
		self.dense = np.zeros((npts))
		self.densi = np.zeros((npts))
		self.betar = np.zeros((npts))
		self.ajBbs = np.zeros((npts))
		for i in range(npts):
			line = lines[i].split()
			self.rhoar[i]    = float(line[0])
			self.gbsnorm[i]  = float(line[1])
			self.amain[i]    = float(line[2])
			self.aiterm1[i]  = float(line[3])
			self.other1[i]   = float(line[4])
			self.dibs[i]     = float(line[5])
			self.bsdense[i]  = float(line[6])
			self.bsdensi[i]  = float(line[7])
			self.bstempe[i]  = float(line[8])
			self.bstempi[i]  = float(line[9])
			self.qsafety[i]  = float(line[10])
			self.ftrapped[i] = float(line[11])
			self.bsnorm[i]   = float(line[12])
			self.tempe1[i]   = float(line[13])
			self.tempi1[i]   = float(line[14])
			self.dense[i]    = float(line[15])
			self.densi[i]    = float(line[16])
			self.betar[i]    = float(line[17])
			self.ajBbs[i]    = float(line[18])
		self.Itotal = np.trapz(self.dibs,self.rhoar)*1.0E6

	def calc_Itotal(self,order=6):
		"""Computes the total current

		This routine computes the total current density from the dibs
		and rhoar arrays.  Note that rhoar is actually s.
		"""
		import numpy as np
		from scipy import integrate
		self.Ibs = integrate.cumulative_trapezoid( \
			self.dibs,self.rhoar, initial=0)
		# Note we fit in rho no s (rhoar is s)
		rho = np.sqrt(self.rhoar)
		x = np.append(-rho[::-1],rho)
		y = np.append(self.Ibs[::-1],self.Ibs)
		z = np.polyfit(x,y,order)
		# Set I(0) = 0.0
		z[-1] = 0.0
		# Set I(1.0) = Itotal
		z = self.Itotal*z/np.sum(z)
		# Save the fit
		p = np.poly1d(z)
		self.Ibs_fit = p(np.sqrt(self.rhoar))
		# Compute the dI/ds = dI/drho/(2*rho) from the fit.
		dz = z[0:-1]
		dp = np.poly1d(dz)
		self.dibs_fit = 0.5*dp(np.sqrt(self.rhoar))/rho
		

if __name__=="__main__":
	import sys
	sys.exit(0)




