##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for reading and handling VMEC
equilibrium data.
"""

# Libraries
#from libstell import libstell
import libstell

# Constants

# VMEC Class
class VMEC(libstell.FourierRep):
	"""Class for working with VMEC equilibria

	"""
	def __init__(self):
		super().__init__()
		self.nfp = None
		self.libStell = libstell.LIBSTELL()

	def read_wout(self,filename):
		"""Reads a VMEC wout_file

		This routine reads and initilizes the VMEC class
		with variable information from a VMEC wout file.

		Parameters
		----------
		file : str
			Path to wout file.
		"""
		import numpy as np
		wout_dict = self.libStell.read_wout(filename)
		for key in wout_dict:
			setattr(self, key, wout_dict[key])
		# (mu-nv) -> (mu+nv)
		self.xn = -self.xn
		# Half to full grid
		self.buco = self.h2f(self.buco)
		self.bvco = self.h2f(self.bvco)
		self.vp = self.h2f(self.vp)
		self.overr = self.h2f(self.overr)
		self.specw = self.h2f(self.specw)
		for mn in range(self.mnmax):
			self.lmns[:,mn] = self.h2f(self.lmns[:,mn])
		for mn in range(self.mnmax_nyq):
			self.bmnc[:,mn] = self.h2f(self.bmnc[:,mn])
			self.gmnc[:,mn] = self.h2f(self.gmnc[:,mn])
			self.bsupumnc[:,mn] = self.h2f(self.bsupumnc[:,mn])
			self.bsupvmnc[:,mn] = self.h2f(self.bsupvmnc[:,mn])
			self.bsubsmns[:,mn] = self.h2f(self.bsubsmns[:,mn])
			self.bsubumnc[:,mn] = self.h2f(self.bsubumnc[:,mn])
			self.bsubvmnc[:,mn] = self.h2f(self.bsubvmnc[:,mn])
		if self.iasym==1:
			for mn in range(self.mnmax):
				self.lmnc[:,mn] = self.h2f(self.lmnc[:,mn])
			for mn in range(self.mnmax_nyq):
				self.bmns[:,mn] = self.h2f(self.bmns[:,mn])
				self.gmns[:,mn] = self.h2f(self.gmns[:,mn])
				self.bsupumns[:,mn] = self.h2f(self.bsupumns[:,mn])
				self.bsupvmns[:,mn] = self.h2f(self.bsupvmns[:,mn])
				self.bsubsmnc[:,mn] = self.h2f(self.bsubsmnc[:,mn])
				self.bsubumns[:,mn] = self.h2f(self.bsubumns[:,mn])
				self.bsubvmns[:,mn] = self.h2f(self.bsubvmns[:,mn])
		# Calc Eplasma
		self.eplasma = 1.5*4*np.pi*np.pi*sum( self.vp * self.presf ) / self.ns

	def h2f(self,var_half):
		"""Half to full grid

		This routine takes a 1D field and interpolates it from the half
		to the full grid. For an ns sized array we assumes that the
		first index [0]=0 and is just a placeholder.

		Parameters
		----------
		var_half : list
			Variable on half grid
		Returns
		----------
		var_full : list
			Variable on full grid
		"""
		temp = var_half.copy()
		temp[0] = 1.5 * temp[1] - 0.5 * temp[2]
		temp[1:-1] = 0.5 * (temp[1:-1] + temp[2:])
		temp[-1] = 2.0 * temp[-1] - 1.0 * temp[-2]
		return temp

	def calc_jll(self, theta, phi ):
		"""Compute parallel current density

		This routine computes the paralle current density 

		Parameters
		----------
		theta : ndarray
			Polidal angle grid [rad]
		phi : ndarray
			Toroidal angle grid [rad]
		Returns
		----------
		jll : ndarray
			Parallel current density [A/m^2]
		"""
		b  = self.cfunct(theta,phi,self.bmnc,self.xm_nyq,self.xn_nyq)
		g  = self.cfunct(theta,phi,self.gmnc,self.xm_nyq,self.xn_nyq)
		bu = self.cfunct(theta,phi,self.bsubumnc,self.xm_nyq,self.xn_nyq)
		bv = self.cfunct(theta,phi,self.bsubvmnc,self.xm_nyq,self.xn_nyq)
		ju = self.cfunct(theta,phi,self.currumnc,self.xm_nyq,self.xn_nyq)
		jv = self.cfunct(theta,phi,self.currvmnc,self.xm_nyq,self.xn_nyq)

		if (self.iasym==1):
			b  = b + self.sfunct(theta,phi,self.bmns,self.xm_nyq,self.xn_nyq)
			g  = g + self.sfunct(theta,phi,self.gmns,self.xm_nyq,self.xn_nyq)
			bu = bu + self.sfunct(theta,phi,self.bsubumns,self.xm_nyq,self.xn_nyq)
			bv = bv + self.sfunct(theta,phi,self.bsubvmns,self.xm_nyq,self.xn_nyq)
			ju = ju + self.sfunct(theta,phi,self.currumns,self.xm_nyq,self.xn_nyq)
			jv = jv + self.sfunct(theta,phi,self.currvmns,self.xm_nyq,self.xn_nyq)
		jll = (bu*ju+bv*jv)/(g*b)
		return jll

	def calc_grad_rhosq(self):
		"""Compute <|grad(rho)|^2> 
		This routine flux surface average of |grad(rho)|^2 

		Parameters
		----------
		theta : ndarray
			Polidal angle grid [rad]
		phi : ndarray
			Toroidal angle grid [rad]
		Returns
		----------
		avgrho2 : ndarray
			Flux surface average of |grad(rho)|^2
		"""
		import numpy as np
		nu = 64
		nv = 128
		theta = np.linspace(0,2*np.pi,nu)
		phi   = np.linspace(0,2*np.pi,nv)
		# Create derivatives
		xm2d  = np.broadcast_to(self.xm,(self.ns,self.mnmax))
		xn2d  = np.broadcast_to(self.xn,(self.ns,self.mnmax))
		rumns = - xm2d * self.rmnc
		rvmns = - xn2d * self.rmnc
		zumnc =   xm2d * self.zmns
		zvmnc =   xn2d * self.zmns
		if self.iasym==1:
			rumnc =   xm2d * self.rmns
			rvmnc =   xn2d * self.rmns
			zumns = - xm2d * self.zmnc
			zvmns = - xn2d * self.zmnc
		r = self.cfunct(theta,zeta,self.rmnc,self.xm,self.xn)
		g = self.cfunct(theta,zeta,self.gmnc,self.xm_nyq,self.xn_nyq)
		ru = self.sfunct(theta,zeta,rumns,self.xm,self.xn)
		rv = self.sfunct(theta,zeta,rvmns,self.xm,self.xn)
		zu = self.cfunct(theta,zeta,zumnc,self.xm,self.xn)
		zv = self.cfunct(theta,zeta,zvmnc,self.xm,self.xn)
		if self.iasym==1:
			r  = r  + self.sfunct(theta,zeta,self.rmns,self.xm,self.xn)
			g  = g  + self.sfunct(theta,zeta,self.gmns,self.xm_nyq,self.xn_nyq)
			ru = ru + self.cfunct(theta,zeta,rumnc,self.xm,self.xn)
			rv = rv + self.cfunct(theta,zeta,rvmnc,self.xm,self.xn)
			zu = zu + self.sfunct(theta,zeta,zumns,self.xm,self.xn)
			zv = zv + self.sfunct(theta,zeta,zvmns,self.xm,self.xn)
		# Calc metrics
		gsr = - zu * r
		gsp = zu * rv - ru * zv
		gsz = ru * r
		gs  = ( gsr * gsr + gsp * gsp + gsz * gsz ) / ( g * g )
		for i in range(self.ns):
			gs[i,:,:] = 0.25 * gs[i,:,:] / ( rho[i] * rho[i] )
		vp = np.sum(g,axis=(1,2))
		avgrho2 = np.sum(gs * g,axis=(1,2)) / vp
		avgrho2[1] = 2.0 * avgrho2[1] - avgrho2[2]
		return avgrho2

	def calc_susceptance(self):
		"""Compute susceptance matrix elements 
		This routine calculates the susceptance matrix elements
		S11, S12, S21, S22.

		Returns
		----------
		S11 : ndarray
			Susceptance matrix elements S11
		S12 : ndarray
			Susceptance matrix elements S11
		S21 : ndarray
			Susceptance matrix elements S11
		S22 : ndarray
			Susceptance matrix elements S11
		"""
		import numpy as np
		nu = 64
		nv = 128
		theta = np.linspace(0,2*np.pi,nu)
		phi   = np.linspace(0,2*np.pi,nv)
		# Create derivatives
		xm2d  = np.broadcast_to(self.xm,(self.ns,self.mnmax))
		xn2d  = np.broadcast_to(self.xn,(self.ns,self.mnmax))
		rumns = - xm2d * self.rmnc
		rvmns = - xn2d * self.rmnc
		zumnc =   xm2d * self.zmns
		zvmnc =   xn2d * self.zmns
		lumnc =   xm2d * self.lmns
		lvmnc =   xn2d * self.lmns
		if self.iasym==1:
			rumnc =   xm2d * self.rmns
			rvmnc =   xn2d * self.rmns
			zumns = - xm2d * self.zmnc
			zvmns = - xn2d * self.zmnc
			lumns = - xm2d * self.lmnc
			lvmns = - xn2d * self.lmnc
		r = self.cfunct(theta,zeta,self.rmnc,self.xm,self.xn)
		g = self.cfunct(theta,zeta,self.gmnc,self.xm_nyq,self.xn_nyq)
		ru = self.sfunct(theta,zeta,rumns,self.xm,self.xn)
		rv = self.sfunct(theta,zeta,rvmns,self.xm,self.xn)
		zu = self.cfunct(theta,zeta,zumnc,self.xm,self.xn)
		zv = self.cfunct(theta,zeta,zvmnc,self.xm,self.xn)
		lu = self.cfunct(theta,zeta,lumnc,self.xm,self.xn)
		lv = self.cfunct(theta,zeta,lvmnc,self.xm,self.xn)
		if self.iasym==1:
			r  = r  + self.sfunct(theta,zeta,self.rmns,self.xm,self.xn)
			g  = g  + self.sfunct(theta,zeta,self.gmns,self.xm_nyq,self.xn_nyq)
			ru = ru + self.cfunct(theta,zeta,rumnc,self.xm,self.xn)
			rv = rv + self.cfunct(theta,zeta,rvmnc,self.xm,self.xn)
			zu = zu + self.sfunct(theta,zeta,zumns,self.xm,self.xn)
			zv = zv + self.sfunct(theta,zeta,zvmns,self.xm,self.xn)
			lu = lu + self.sfunct(theta,zeta,lumns,self.xm,self.xn)
			lv = lv + self.sfunct(theta,zeta,lvmns,self.xm,self.xn)
		# Calc suscpetance matrices
		scale_fact = 1.0 / ( 4 * np.pi * np.pi )
		S11 = ( ru * ru + zu * zu)
		S21 = ( ru * rv + zu * zv)
		S12 = ( S12 * ( 1.0 + lu ) - S11 * lv )
		S22 = ( ( rv * rv + zv * zv + r * r ) * ( 1.0 + lu ) - S21 * lv )
		S11 = np.trapz(S11 / g, x=zeta, axis=2)
		S12 = np.trapz(S12 / g, x=zeta, axis=2)
		S21 = np.trapz(S21 / g, x=zeta, axis=2)
		S22 = np.trapz(S22 / g, x=zeta, axis=2)
		S11 = np.trapz(S11, x=theta, axis=1)*scale_fact
		S12 = np.trapz(S12, x=theta, axis=1)*scale_fact
		S21 = np.trapz(S21, x=theta, axis=1)*scale_fact
		S22 = np.trapz(S22, x=theta, axis=1)*scale_fact
		return S11,S12,S21,S22


# VMEC INDATA Class
class VMEC_INDATA():
	"""Class for working with VMEC equilibria

	"""
	def __init__(self, parent=None):
		self.nfp = None
		self.libStell = libstell.LIBSTELL()

	def read_indata(self,filename):
		"""Reads INDATA namelist from a file

		This routine wrappers the vmec_input module reading routine.
		Parameters
		----------
		filename : string
			Input file name with INDATA namelist
		"""
		indata_dict = self.libStell.read_indata(filename)
		for key in indata_dict:
			setattr(self, key, indata_dict[key])
		# generate helpers
		#print(self.rbc.shape)



# Main routine
if __name__=="__main__":
	import sys
	from argparse import ArgumentParser
	import numpy as np
	import matplotlib.pyplot as pyplot
	parser = ArgumentParser(description= 
		'''Provides class for accessing vmec data also servers as a
		   simple tool for assessing vmec wout or input files.''')
	parser.add_argument("-v", "--vmec", dest="vmec_ext",
		help="VMEC file extension", default = None)
	parser.add_argument("-p", "--plot", dest="lplot", action='store_true',
		help="Plot the coils file.", default = False)
	args = parser.parse_args()
	vmec_wout = VMEC()
	vmec_input = VMEC_INDATA()
	libStell = libstell.LIBSTELL()
	if args.vmec_ext:
		linput = False
		loutput = False
		try:
			vmec_input.read_indata('input.'+args.vmec_ext)
			linput = True
		except:
			print(f'Could not file input file: input.{args.vmec_ext}')
		try:
			vmec_wout.read_wout(args.vmec_ext)
			loutput = True
		except:
			print(f'Could not file input file: wout_{args.vmec_ext}.nc or wout.{args.vmec_ext}')
		if not (linput or loutput): sys.exit(-1)
		# Do Input file plot
		if (args.lplot and linput):
			px = 1/pyplot.rcParams['figure.dpi']
			fig=pyplot.figure(figsize=(1024*px,768*px))
			ax=fig.add_subplot(221)
			pyplot.subplots_adjust(hspace=0.4,wspace=0.3)
			s = np.linspace(0.0,1.0,128)
			f = np.zeros(128)
			for i,x in enumerate(s):
				f[i] = libStell.pmass(x)
			ax.plot(s,f/1E3,'k')
			ax.text(0.02,0.19,rf'PHIEDGE={vmec_input.phiedge:4.3f} [Wb]', horizontalalignment='left',\
				verticalalignment='center', transform=ax.transAxes)
			ax.text(0.02,0.12,rf'PRESS_SCALE={vmec_input.pres_scale}', horizontalalignment='left',\
				verticalalignment='center', transform=ax.transAxes)
			ax.text(0.02,0.05,rf'MGRID_FILE: {vmec_input.mgrid_file}', horizontalalignment='left',\
				verticalalignment='center', transform=ax.transAxes)
			ax.text(0.98,0.95,f'PMASS_TYPE: {vmec_input.pmass_type}', horizontalalignment='right',\
				verticalalignment='center', transform=ax.transAxes)
			ax.set_xlim(0,1)
			ax.set_xlabel('Norm. Tor. Flux (s)')
			ax.set_ylabel("Pressure [kPa]")
			ax.set_title(f'VMEC Input: {args.vmec_ext}')
			ax=fig.add_subplot(222)
			if vmec_input.ncurr==1:
				for i,x in enumerate(s):
					f[i] = libStell.pcurr(x)
					ax.set_ylabel("Current dI/ds")
					temp_str='PCURR'
					temp_type = vmec_input.pcurr_type
			else:
				for i,x in enumerate(s):
					f[i] = libStell.piota(x)
					ax.set_ylabel(r"$\iota$")
					temp_str='PIOTA'
					temp_type = vmec_input.piota_type
			ax.plot(s,f,'k')
			ax.set_xlabel('Norm. Tor. Flux (s)')
			ax.text(0.98,0.95,f'{temp_str}_TYPE: {temp_type}', horizontalalignment='right',\
				verticalalignment='center', transform=ax.transAxes)
			ax.text(0.02,0.12,rf'CURTOR: {vmec_input.curtor:6.1f}', horizontalalignment='left',\
				verticalalignment='center', transform=ax.transAxes)
			ax.text(0.02,0.05,rf'NCURR: {vmec_input.ncurr}', horizontalalignment='left',\
				verticalalignment='center', transform=ax.transAxes)
			pyplot.show()
			
		# Do wout file plot
		if (args.lplot and loutput):
			px = 1/pyplot.rcParams['figure.dpi']
			fig=pyplot.figure(figsize=(1024*px,768*px))
			ax=fig.add_subplot(221)
			pyplot.subplots_adjust(hspace=0.4,wspace=0.3)
			ax.plot(np.linspace(0.0,1.0,vmec_wout.ns),vmec_wout.presf/1E3,'k')
			ax.text(0.02,0.47,rf'$B_0$={vmec_wout.b0:4.3f} [T]', horizontalalignment='left',\
				verticalalignment='center', transform=ax.transAxes)
			ax.text(0.02,0.40,rf'$R/a$={vmec_wout.aspect:4.3f}', horizontalalignment='left',\
				verticalalignment='center', transform=ax.transAxes)
			ax.text(0.02,0.33,rf'$R$={vmec_wout.rmajor:4.3f}', horizontalalignment='left',\
				verticalalignment='center', transform=ax.transAxes)
			ax.text(0.02,0.26,rf'$a$={vmec_wout.aminor:4.3f}', horizontalalignment='left',\
				verticalalignment='center', transform=ax.transAxes)
			ax.text(0.02,0.19,rf'$W_p$={vmec_wout.eplasma[0]/1000:4.3f} [kJ]', horizontalalignment='left',\
				verticalalignment='center', transform=ax.transAxes)
			ax.text(0.02,0.12,rf'$<\beta>$={vmec_wout.betatot*100:3.2f} %', horizontalalignment='left',\
				verticalalignment='center', transform=ax.transAxes)
			ax.text(0.02,0.05,rf'MGRID_FILE: {vmec_wout.mgrid_file}', horizontalalignment='left',\
				verticalalignment='center', transform=ax.transAxes)
			ax.text(0.98,0.95,f'PMASS_TYPE: {vmec_wout.pmass_type}', horizontalalignment='right',\
				verticalalignment='center', transform=ax.transAxes)
			ax.set_xlim(0,1)
			ax.set_xlabel('Norm. Tor. Flux (s)')
			ax.set_ylabel("Pressure [kPa]")
			ax.set_title(f'VMEC Ouput: {args.vmec_ext}')
			ax=fig.add_subplot(222)
			ax2 = ax.twinx()
			ax.plot(np.linspace(0.0,1.0,vmec_wout.ns),vmec_wout.iotaf,'k')
			ax2.plot(np.linspace(0.0,1.0,vmec_wout.ns),vmec_wout.jcurv/1E3,'r')
			ax.text(0.05,0.05,rf'$I$={vmec_wout.itor/1000:4.2f} [kA]', horizontalalignment='left',\
				verticalalignment='center', transform=ax.transAxes)
			ax.text(0.98,0.95,f'PCURR_TYPE: {vmec_wout.pcurr_type}', horizontalalignment='right',\
				verticalalignment='center', transform=ax.transAxes)
			ax.text(0.98,0.88,f'PIOTA_TYPE: {vmec_wout.piota_type}', horizontalalignment='right',\
				verticalalignment='center', transform=ax.transAxes)
			ax.set_xlim(0,1)
			ax.set_xlabel('Norm. Tor. Flux (s)')
			ax.set_ylabel(r"$\iota$",color='k')
			ax2.set_ylabel(r"$<j_v>$ [$kA/m^2$]",color='r')
			ax=fig.add_subplot(223)
			theta = np.ndarray((360,1))
			zeta  = np.ndarray((3,1))
			for j in range(360): theta[j]=2.0*np.pi*j/359.0
			for j in range(3):   zeta[j]=     np.pi*j/2.0
			r = vmec_wout.cfunct(theta,zeta,vmec_wout.rmnc,vmec_wout.xm,vmec_wout.xn/vmec_wout.nfp)
			z = vmec_wout.sfunct(theta,zeta,vmec_wout.zmns,vmec_wout.xm,vmec_wout.xn/vmec_wout.nfp)
			ax.plot(r[1,1,0],z[1,1,0],'+r')
			ax.plot(r[1,1,1],z[1,1,1],'+g')
			ax.plot(r[1,1,2],z[1,1,2],'+b')
			j = vmec_wout.ns-1
			ax.plot(r[j,:,0],z[j,:,0],'r')
			ax.plot(r[j,:,1],z[j,:,1],'g')
			ax.plot(r[j,:,2],z[j,:,2],'b')
			j = int(vmec_wout.ns/4)
			ax.plot(r[j,:,0],z[j,:,0],'--r')
			ax.plot(r[j,:,1],z[j,:,1],'--g')
			ax.plot(r[j,:,2],z[j,:,2],'--b')
			ax.set_aspect('equal', adjustable='box')
			ax.text(0.02,0.05,rf'NFP: {vmec_wout.nfp}', horizontalalignment='left',\
				verticalalignment='center', transform=ax.transAxes)
			ax=fig.add_subplot(224)
			theta = np.ndarray((256,1))
			zeta  = np.ndarray((256,1))
			for j in range(256): theta[j]=2.0*np.pi*j/255.0
			for j in range(256):  zeta[j]=2.0*np.pi*j/255.0
			b = vmec_wout.cfunct(theta,zeta,vmec_wout.bmnc,vmec_wout.xm_nyq,vmec_wout.xn_nyq/vmec_wout.nfp)
			j = int(vmec_wout.ns/4)
			h=ax.pcolormesh(np.squeeze(b[j,:,:]),cmap='jet',shading='gouraud')
			ax.set_xlabel(r"$\zeta [rad]$")
			ax.set_ylabel(r"$\theta_{VMEC}$ [rad]")
			ax.set_title("|B| at mid radius")
			fig.colorbar(h,label='[T]')
			pyplot.show()
	sys.exit(0)