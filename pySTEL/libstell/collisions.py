##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for calculating
quantities related to collisions
"""

# Libraries
#from libstell import libstell

# Constants
EC = 1.602176634E-19 # Electron charge [C]
HBAR = 1.054571817E-34 # Planks reduced constant [J.s]
DA = 1.66053906660E-27 # Dalton
ME = 9.1093837E-31 # Electron mass [kg]
C  = 299792458 # Speed of light [m/s]

# VMEC Class
class COLLISIONS():
	"""Class for working with FUSION data

	"""
	def __init__(self):
		# Electron Mass
		self.ME = ME
		# Proton Mass
		self.MP = 1.007276466621*DA 
		# Deuturium Mass
		self.MD = 2.01410177811*DA
		# Tritium Mass
		self.MT = 3.01604928*DA
		# Helium Mass
		self.MHe3 = 3.0160293*DA
		self.MHe4 = 4.002603254*DA

	def coullog_ee(self,ne,te):
		"""Computes the Coulomb e-e logarithm

		This routine calculates the Coulomb Logarithm acording to the
		NRL plasma formulary general deffinition.

		Parameters
		----------
		ne : real
			Electron density [m^-3]
		te : real
			Electron Temperature [eV]
		Returns
		----------
		clog : real
			Coulomb logarithm
		"""
		import numpy as np
		ne_cm = ne*1E-6
		a1 = -np.log(np.sqrt(ne_cm)*te**(-1.25))
		a2 = -np.sqrt(1E-5 + ((np.log(te)-2)**2)/16.0)
		return 23.5 + a1 + a2

	def coullog_ei(self,ne,te,mi,Z,ni,ti):
		"""Computes the Coulomb e-i logarithm

		This routine calculates the Coulomb Logarithm acording to the
		NRL plasma formulary general deffinition.

		Parameters
		----------
		ne : real
			Electron density [m^-3]
		te : real
			Electron Temperature [eV]
		mi : real
			Ion mass [kg]
		Z : real
			Ion Charge number
		ni : real
			Ion density [m^-3]
		ti : real
			Ion Temperature [eV]
		Returns
		----------
		clog : real
			Coulomb logarithm
		"""
		import numpy as np
		ne_cm = ne*1E-6
		mu = mi/self.MP
		# No need to convert masses from kg to g
		if ti/mi < (te/self.ME):
			if (te < (10 * Z * Z)):
				clog = 23 - np.log(np.sqrt(ne_cm)*Z*te**(-1.5))
			else:
				clog = 24 - np.log(np.sqrt(ne_cm)/te)
		else:
			mu = (self.ME*mi)/(self.ME+mi)
			clog  = 16 - np.log(mu*Z*Z*np.sqrt(ni)*ti**-1.5)
		return clog

	def coullog_ii(self,mi1,Z1,ni1,ti1,mi2,Z2,ni2,ti2):
		"""Computes the Coulomb i-i logarithm

		This routine calculates the Coulomb Logarithm acording to the
		NRL plasma formulary general deffinition.

		Parameters
		----------
		mi1 : real
			Ion #1 mass [kg]
		Z1  : real
			Ion #1 Charge number
		ni1 : real
			Ion #1 density [m^-3]
		ti1 : real
			Ion #1 Temperature [eV]
		mi2 : real
			Ion #2 mass [kg]
		Z2  : real
			Ion #2 Charge number
		ni2 : real
			Ion #2 density [m^-3]
		ti2 : real
			Ion #2 Temperature [eV]
		Returns
		----------
		clog : real
			Coulomb logarithm
		"""
		import numpy as np
		ni1_cm = ni1*1E-6
		ni2_cm = ni2*1E-6
		mu1 = mi1/self.MP
		mu2 = mi2/self.MP
		# No need to convert masses from kg to g
		clog = (ni1_cm*Zi1*Zi1/ti1) + (ni2_cm*Zi2*Zi2/ti2)
		clog = Zi1*Zi2*(mu1+mu2)*np.sqrt(clog)/(mu1*ti2+mu2*ti1)
		clog = 23 - np.log(clog)
		return clog

	def coullog_iifast(self,ne,te,mi1,Z1,mi2,Z2,vd):
		"""Computes the Coulomb counterstreaming i-i logarithm

		This routine calculates the Coulomb Logarithm acording to the
		NRL plasma formulary general deffinition.

		Parameters
		----------
		ne : real
			Electron density [m^-3]
		te : real
			Electron Temperature [eV]
		mi1 : real
			Ion #1 mass [kg]
		Z1  : real
			Ion #1 Charge number
		mi2 : real
			Ion #2 mass [kg]
		Z2  : real
			Ion #2 Charge number
		vd  : real
			Ion relative velocity [m/s]
		Returns
		----------
		clog : real
			Coulomb logarithm
		"""
		import numpy as np
		ne_cm = ne*1E-6
		mu1 = mi1/self.MP
		mu2 = mi2/self.MP
		# No need to convert masses from kg to g
		beta = vd/C
		clog = ne_cm/te
		clog = Z1*Z2*(mu1+mu2) * np.sqrt(clog) / (mu1*mu2*beta*beta)
		clog = 43-np.log(clog)
		return clog

	def collisionfreq(self,m1,Z1,T1,m2,Z2,n2,T2,clog):
		"""Computes the collision frequency for two species

		This routine calculates the collision frequency acording to the
		NRL plasma formulary general deffinition.

		Parameters
		----------
		m1 : real
			Particle #1 mass [kg]
		Z1 : real
			Particle #1 Charge number
		T1 : real
			Particle #1 Temperature [eV]
		m2 : real
			Particle #2 mass [kg]
		Z2 : real
			Particle #2 Charge number
		n2 : real
			Particle #2 Density [m^-3]
		T2 : real
			Particle #2 Temperature [eV]
		clog : real
			Coulomb logarithm
		Returns
		----------
		freq : real
			Collision frequency
		"""
		import numpy as np
		n2_cm = n2*1E-6
		# factor of 1000 comes from sqrt(mass in kg to g)
		freq = 1.8E-19*np.sqrt(m1*m2)*Z1*Z1*Z2*Z2*n2_cm*clog*1000
		freq = freq / (m1*T2+m2*T1)**1.5
		return freq

	def tauspitzer(self,mi,Zi,ne,Te,clog):
		"""Computes the Spitzer ion-electron exchange time

		This routine calculates the Spitzer ion-electron momentum
		exchange time.

		Parameters
		----------
		mi : real
			Ion Mass [kg]
		Zi : real
			Ion Charge number
		ne : real
			Electron Density [m^-3]
		Te : real
			Electron Temperature [eV]
		clog : real
			Coulomb logarithm
		Returns
		----------
		tau : real
			Spitzer ion-electron momentum exchange time
		"""
		import numpy as np
		return 3.777183E41*mi*np.sqrt(Te*Te*Te)/(ne*Zi*Zi*clog)

	def criticalvelocity(self,mp,Te):
		"""Computes the critical velocity

		This routine calculates the critical velocity for an energetic
		particle.

		Parameters
		----------
		mi : real
			Plasma Ion Mass [kg]
		Te : real
			Electron Temperature [eV]
		Returns
		----------
		vcrit : real
			Critical velocity [m/s]
		"""
		import numpy as np
		return (0.75*np.sqrt(np.pi*mp/self.ME))**(1.0/3.0) * np.sqrt(2*EC*Te/mp)

if __name__=="__main__":
	import sys
	sys.exit(0)
