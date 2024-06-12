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

	def bfield(self,x,y,z):
		"""Evaluates the magnetic field from a NESCOIL surface current

		This routine evaluates the magnetic field at a point in space
		from a NESCOIL surface current.

		Parameters
		----------
		x : float
			X point to evaluate [m]
		y : float
			Y point to evaluate [m]
		z : float
			Z point to evaluate [m]
		Returns
		----------
		bx : float
			X-component of magnetic field [T]
		by : float
			Y-component of magnetic field [T]
		bZ : float
			Z-component of magnetic field [T]
		"""
		return self.libStell.nescout_bfield(x,y,z,istat=-327)

	def generatePotential(self,theta,zeta):
		"""Computes the potential on a grid

		This routine computes the potential on a grid

		Parameters
		----------
		theta : ndarray
			Poloidal Angle [rad]
		zeta : ndarray
			Field Period Angle [rad]
		Returns
		-------
		pot : ndarray
			NESCOIL potential []
		"""
		pot = self.sfunct(theta,zeta,self.potmns_surface.T,self.xm_pot,self.xn_pot)
		return pot

	def generateTotalPotential(self,theta,zeta):
		"""Computes the potential on a grid

		This routine computes the potential on a grid

		Parameters
		----------
		theta : ndarray
			Poloidal Angle [rad]
		zeta : ndarray
			Field Period Angle [rad]
		Returns
		-------
		pot : ndarray
			NESCOIL potential []
		"""
		import numpy as np
		pot = self.generatePotential(theta,zeta)
		nu = len(theta)
		nv = len(zeta)
		for j in range(self.nu): pot[0,j,:] = pot[0,j,:] - self.cut*0.5*theta[j]/np.pi
		for j in range(self.nv): pot[0,:,j] = pot[0,:,j] - self.cup*0.5*zeta[j]/(np.pi*self.np)
		return pot


	def plotpotential(self,ax=None):
		"""Plots the NESCOIL Potential

		This routine plots the NESCOIL code surface potential

		Parameters
		----------
		ax : axes (optional)
			Matplotlib axes object to plot to.
		"""
		import numpy as np
		import matplotlib.pyplot as pyplot
		lplotnow = False
		if not ax:
			ax = pyplot.axes()
			lplotnow = True
		theta = np.ndarray((self.nu,1))
		zeta  = np.ndarray((self.nv,1))
		for j in range(self.nu): theta[j]=2.0*np.pi*j/float(self.nu-1)
		for j in range(self.nv):  zeta[j]=    np.pi*j/float(self.nv-1)
		pot = self.generatePotential(theta,zeta)
		#pot = self.sfunct(theta,zeta,self.potmns_surface.T,self.xm_pot,self.xn_pot)
		hmesh=ax.pcolormesh(np.squeeze(zeta),np.squeeze(theta),np.squeeze(pot[0,:,:]),cmap='jet',shading='gouraud')
		ax.set_xlabel('Toroidal angle [rad]')
		ax.set_ylabel('Poloidal angle [rad]')
		ax.set_title(r'NESCOIL $\Phi$ Potential')
		pyplot.colorbar(hmesh,label='$Pot$ [arb]',ax=ax)

	def plottotalpotential(self,ax=None):
		"""Plots the NESCOIL Total Potential

		This routine plots the NESCOIL code surface potential

		Parameters
		----------
		ax : axes (optional)
			Matplotlib axes object to plot to.
		"""
		import numpy as np
		import matplotlib.pyplot as pyplot
		lplotnow = False
		if not ax:
			ax = pyplot.axes()
			lplotnow = True
		theta = np.ndarray((self.nu,1))
		zeta  = np.ndarray((self.nv,1))
		for j in range(self.nu): theta[j]=2.0*np.pi*j/float(self.nu-1)
		for j in range(self.nv):  zeta[j]=    np.pi*j/float(self.nv-1)
		pot = self.generateTotalPotential(theta,zeta)
		#pot = self.sfunct(theta,zeta,self.potmns_surface.T,self.xm_pot,self.xn_pot)
		#for j in range(self.nu): pot[0,j,:] = pot[0,j,:] - self.cut*1.0*j/float(self.nv-1)
		#for j in range(self.nv): pot[0,:,j] = pot[0,:,j] - self.cup*0.5*j/float(self.nv-1)
		hmesh=ax.pcolormesh(np.squeeze(zeta),np.squeeze(theta),np.squeeze(pot[0,:,:]),cmap='jet',shading='gouraud')
		#ax.semilogy(abscissa, data, **kwargs)
		ax.set_xlabel('Toroidal angle [rad]')
		ax.set_ylabel('Poloidal angle [rad]')
		ax.set_title(r'NESCOIL Total $\Phi$ Potential')
		pyplot.colorbar(hmesh,label='$Pot$ [arb]',ax=ax)
		if lplotnow: pyplot.show()

	def plotsurfaces(self,ax=None):
		"""Plots the NESCOIL Surfaces

		This routine plots the NESCOIL current potential surface
		and the plasma surface over a half field period.

		Parameters
		----------
		ax : axes (optional)
			Matplotlib axes object to plot to.
		"""
		import numpy as np
		import matplotlib.pyplot as pyplot
		lplotnow = False
		if not ax:
			#fig = pyplot.figure()
			#ax = fig.add_subplot(111,projection='3d')
			lplotnow = True
		theta = np.ndarray((self.nu,1))
		zeta  = np.ndarray((self.nv,1))
		for j in range(self.nu): theta[j]=2.0*np.pi*j/float(self.nu-1)
		for j in range(self.nv):  zeta[j]=np.pi*j/float(self.nv-1)
		rp = self.cfunct(theta,zeta,self.rmnc_plasma.T,self.xm_plasma,self.xn_plasma)
		zp = self.sfunct(theta,zeta,self.zmns_plasma.T,self.xm_plasma,self.xn_plasma)
		rc = self.cfunct(theta,zeta,self.rmnc_surface.T,self.xm_surface,self.xn_surface)
		zc = self.sfunct(theta,zeta,self.zmns_surface.T,self.xm_surface,self.xn_surface)
		ax = self.isotoro(rp,zp,zeta,0,color='red')
		hc = self.isotoro(rc,zc,zeta,0,axes=ax,color='green')
		#ax.set_aspect('equal', 'box')
		#ax.semilogy(abscissa, data, **kwargs)
		if lplotnow: pyplot.show()

	def cutcoils(self,ncoils_per_halfperiod):
		"""Cut coils from the NESCOIL potential

		This routine cuts coils from the NESCOIL potential.
		It allows the user to specify the number of coils per half 
		period.

		Parameters
		----------
		ncoils_per_halfperiod : integer
			Number of coils per half period (suggest 5)
		"""
		import numpy as np
		from libstell.coils import COILSET, COILGROUP, COIL
		from contourpy import contour_generator, LineType
		# Generate coilset
		coils = COILSET()
		coils.nfp = self.np
		coils.ngroups = ncoils_per_halfperiod
		coils.xmin = 1E9; coils.xmax=-1E9
		coils.ymin = 1E9; coils.ymax=-1E9
		coils.zmin = 1E9; coils.zmax=-1E9
		#coils.groups = []
		# First calculate the potential on a grid
		theta = np.ndarray((self.nu,1))
		zeta  = np.ndarray((self.nv,1))
		for j in range(self.nu): theta[j]=2.0*np.pi*j/float(self.nu-1)
		for j in range(self.nv):  zeta[j]=np.pi*j/float(self.nv-1)
		pot = self.generateTotalPotential(theta,zeta)
		# Now generate contours
		potmin = np.min(pot)
		potmax = np.max(pot)
		delta  = 2*(potmax-potmin)/float(2*ncoils_per_halfperiod+1.5)
		cont_gen = contour_generator(x=np.squeeze(zeta),y=np.squeeze(theta),z=np.squeeze(pot), line_type=LineType.Separate)
		# Now loop over contours
		for k in range(ncoils_per_halfperiod):
			u = round(self.nu/2)
			v = round((k+0.5)*self.nv/(ncoils_per_halfperiod))
			level = cont_gen.lines(pot[0,u,v])
			level = level[0]
			th = level[:,1]
			ze = level[:,0]
			# Fourier transform the coil
			r = th*0.0 ;z = th*0.0
			for mn in range(self.mnmax_surface):
				mtheta = th*self.xm_surface[mn]
				nzeta  = ze*self.xn_surface[mn]
				r  = r + np.cos(mtheta+nzeta)*self.rmnc_surface[mn]
				z  = z + np.sin(mtheta+nzeta)*self.zmns_surface[mn]
			# Convert to XYZ and make current/group
			ph = ze/self.np
			x = r * np.cos(ph)
			y = r * np.sin(ph)
			c = x * 0.0 + self.curpol/(np.pi*4E-7)
			g = x * 0.0 + k+1
			c[-1] = 0.0
			# Create stellarator symmetric coil
			ph = (2.0*np.pi - ze)/self.np
			x2 = r * np.cos(ph)
			y2 = r * np.sin(ph)
			z2 =-z
			x2=x2[::-1]; y2=y2[::-1]; z2 = z2[::-1]
			g2 = g
			c2 = c
			xo = np.append(x,x2)
			yo = np.append(y,y2)
			zo = np.append(z,z2)
			co = np.append(c,c)
			go = np.append(g,g)
			co[-1] = 0.0
			print(zo)
			x  = xo; y = yo; z = zo; c = co; g =go
			# Now make all field periods
			for mn in range(1,self.np):
				cop = np.cos(mn*self.alp)
				sip = np.sin(mn*self.alp)
				x = np.append(x,xo*sip)
				y = np.append(y,yo*sip)
				z = np.append(y,zo)
				c = np.append(c,co)
				g = np.append(g,go)
			#print(c)
			coils.xmin = np.minimum(coils.xmin,np.min(x))
			coils.ymin = np.minimum(coils.ymin,np.min(y))
			coils.zmin = np.minimum(coils.zmin,np.min(z))
			coils.xmax = np.maximum(coils.xmax,np.max(x))
			coils.ymax = np.maximum(coils.ymax,np.max(y))
			coils.zmax = np.maximum(coils.zmax,np.max(z))
			# Now create group
			coil_name=f'MOD{k+1}'
			coils.groups.extend([COILGROUP(x,y,z,c,coil_name)])
		# Return a coil object
		return coils










# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)