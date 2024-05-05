##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for interfacing to libstell
"""

# Libraries

# Constants

# VMEC Class
class LIBSTELL():
	"""Class for working with VMEC equilibria

	"""
	def __init__(self, parent=None):
		import os,sys
		import ctypes as ct
		from subprocess import Popen, PIPE
		self.STELLOPT_PATH = os.environ["STELLOPT_PATH"]
		self.PATH_TO_LIBSTELL = os.path.join(self.STELLOPT_PATH,'LIBSTELL','Release','libstell.so')
		try:
			self.libstell = ct.cdll.LoadLibrary(self.PATH_TO_LIBSTELL)
		except:
			print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
			print("!!  Could not load shared libraray libstell.so    !!")
			print(f"!!  PATH: {self.PATH_TO_LIBSTELL}    !!")
			print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
		# Figure out uderscoring
		out = Popen(args="nm "+self.PATH_TO_LIBSTELL, 
			shell=True, 
			stdout=PIPE).communicate()[0].decode("utf-8")
		attrs = [ i.split(" ")[-1].replace("\r", "") \
			for i in out.split("\n") if " T " in i]
		func = '_write_wout_file'
		module = 'read_wout_mod_' 
		names = [ s for s in attrs if module in s and func in s]
		name = names[0].replace(module, ',')
		name = name.replace(func, ',')
		self.s1, self.s2, self.s3 = name.split(',')
		# Weird OSX behavior
		if self.s1=='___':
			self.s1='__'

	def safe_open(self,iunit,istat,filename,filestat,fileform,record_in=None,access_in='SEQUENTIAL',delim_in='APOSTROPHE'):
		"""Wrapper to safe_open in safe_open_mod

		This routine wrappers safe_open from safe_open_mod

		Parameters
		----------
		iunit : int
			Fortran file unit number
		istat : int
			Status value
		filename : str
			Filename string
		filestat : str
			Fortran open file status ('OLD' 'NEW' 'UNKNOWN' 'SCRATCH')
		fileform : str
			Fortran open file form ('FORMATTED' 'UNFORMATTED' 'PRINT')
		record_in : int
			Fortran open record size
		access_in : str
			Fortran open access ('APPEND' 'DIRECT' 'SEQUENTIAL')
		delim_in : str
			Fortran open delimter ()
		"""
		import ctypes as ct
		module_name = self.s1+'safe_open_mod_'+self.s2
		safe_open_h = getattr(self.libstell,module_name+'_safe_open'+self.s3)
		# SUBROUTINE safe_open(int iunit, int istat, char filename, char filestat, char fileform, int record_in, char access_in, char delim_in)
		safe_open_h.argtypes= [ ct.POINTER(ct.c_int), ct.POINTER(ct.c_int), ct.c_char_p, ct.c_char_p, ct.c_char_p, \
			ct.POINTER(ct.c_int), ct.c_char_p, ct.c_char_p, \
			ct.c_long, ct.c_long, ct.c_long, ct.c_long, ct.c_long]
		safe_open_h.restype=None
		iunit_temp = ct.c_int(iunit)
		istat_temp = ct.c_int(istat)
		opt1 = ct.c_bool(True)
		opt2 = ct.c_bool(True)
		opt3 = ct.c_bool(True)
		#if record_in:
		if record_in:
			record_in_temp = ct.c_int(record_in)
		else:
			record_in_temp = ct.c_int(0)
		safe_open_h(ct.byref(iunit_temp),ct.byref(istat_temp), \
			filename.encode('UTF-8'),filestat.encode('UTF-8'),fileform.encode('UTF-8'),\
			ct.byref(record_in_temp),access_in.encode('UTF-8'),delim_in.encode('UTF-8'), \
			len(filename),len(filestat),len(fileform),len(access_in),len(delim_in))
		istat = istat_temp
		iunit = iunit_temp
		return

	def safe_close(self,iunit):
		"""Wrapper to safe_close in safe_open_mod

		This routine wrappers safe_close from safe_open_mod

		Parameters
		----------
		iunit : int
			Fortran file unit number
		"""
		import ctypes as ct
		module_name = self.s1+'safe_open_mod_'+self.s2
		safe_close_h = getattr(self.libstell,module_name+'_safe_close'+self.s3)
		safe_close_h.restype=None
		iunit_temp = ct.c_int(iunit)
		safe_close_h(ct.byref(iunit_temp))
		return

	def read_indata(self,filename):
		"""Reads a VMEC INDATA namelist

		This routine wrappers read_indata_namelist function in
		vmec_input. To do this we wrapper safe_open_mod

		Parameters
		----------
		file : str
			Path to wout file.
		"""
		import ctypes as ct
		# These are defined in vparams.f but for some reason they're no in to .so
		ntord  = 101
		mpol1d = 101
		ndatafmax = 101
		# We use an added routine as a helper
		module_name = self.s1+'vmec_input_'+self.s2
		read_indata_namelist = getattr(self.libstell,module_name+'_read_indata_namelist_byfile'+self.s3)
		read_indata_namelist.argtypes = [ct.c_char_p,ct.c_long]
		read_indata_namelist.restype=None
		read_indata_namelist(filename.encode('UTF-8'),len(filename))
		# Get vars
		booList=['lpofr','lmac','lfreeb','lrecon','loldout','ledge_dump','lasym','lforbal','lrfp',\
			'lmovie','lmove_axis','lwouttxt','ldiagno','lmoreiter','lfull3d1out','l_v3fit',\
			'lspectrum_dump','loptim','lgiveup','lbsubs','lgiveup']
		booLen=[1]*len(booList)
		intList=['nfp','ncurr','nsin','niter','nstep','nvacskip','mpol','ntor','ntheta','nzeta', \
			'mfilter_fbdy','nfilter_fbdy','max_main_iterations','omp_num_threads',\
			'imse','isnodes','itse','ipnodes','iopt_raxis','imatch_phiedge','nflxs']
		intLen=[1]*len(intList)
		intList.extend(['ns_array','niter_array'])
		intLen.extend([(100,1)]*2)
		realList=['time_slice','curtor','delt','ftol','tcon0','gamma','phiedge','phidiam',\
			'sigma_current','sigma_delphid','tensi','tensp','tensi2','fpolyi','presfac',\
			'mseangle_offset','pres_offset','mseangle_offsetm','spres_ped','bloat',\
			'pres_scale','prec2d_threshold','bcrit','fgiveup']
		realLen=[1]*len(realList)
		realList.extend(['rbc','zbs','rbs','zbc'])
		realLen.extend([(2*ntord+1,mpol1d+1)]*4)
		realList.extend(['am_aux_s','am_aux_f','ai_aux_s','ai_aux_f','ac_aux_s','ac_aux_f'])
		realLen.extend([(ndatafmax,1)]*6)
		realList.extend(['am','ai','ac','aphi','ah','at'])
		realLen.extend([(21,1)]*6)
		realList.extend(['ah_aux_s','ah_aux_f','at_aux_s','at_aux_f'])
		realLen.extend([(ndatafmax,1)]*4)
		realList.extend(['raxis','zaxis','raxis_cc','raxis_cs','zaxis_cc','zaxis_cs'])
		realLen.extend([(ntord+1,1)]*6)
		realList.extend(['ftol_array','extcur'])
		realLen.extend([(100,1),(300,1)])
		charList=['pcurr_type','piota_type','pmass_type','pt_type','ph_type']
		charLen=[(20,1)]*6
		charList.extend(['mgrid_file','input_extension'])
		charLen.extend([(200,1),(200,1)])
		out_data = self.get_module_vars(module_name,booList,booLen,intList,intLen,realList,realLen,charList,charLen,ldefined_size_arrays=True)
		return out_data

	def write_indata(self,filename):
		"""Wrappers writing of the VMEC INDATA namelist

		This routine wrappers write_indata_namelist in LIBSTELL

		Parameters
		----------
		file : str
			Path to input file.
		"""
		import ctypes as ct
		module_name = self.s1+'vmec_input_'+self.s2
		write_indata_namelist = getattr(self.libstell,module_name+'_write_indata_namelist_byfile'+self.s3)
		write_indata_namelist.argtypes = [ct.c_char_p,ct.c_long]
		write_indata_namelist.restype=None
		write_indata_namelist(filename.encode('UTF-8'),len(filename))

	def read_wout(self,file):
		"""Reads a wout file and returns a dictionary

		This routine wrappers read_wout in LIBSTELL and returns
		a dictionary of values

		Parameters
		----------
		file : str
			Path to wout file.
		Returns
		----------
		vars : dict
			Dictionary of module variables
		"""
		import ctypes as ct
		module_name = self.s1+'read_wout_mod_'+self.s2
		read_wout = getattr(self.libstell,module_name+'_readw_and_open'+self.s3)
		read_wout.argtypes=[ct.c_char_p, ct.POINTER(ct.c_int), ct.POINTER(ct.c_int), ct.c_long]
		read_wout.restype=None
		ierr = ct.c_int(0)
		iopen = ct.c_int(0)
		read_wout(file.encode('UTF-8'), ct.byref(ierr), ct.byref(iopen), len(file))
		if not (ierr.value == 0):
			return None
		# Setup Arrays
		out_data={}
		# Get Scalars
		booList  = ['lasym','lthreed','lwout_opened']
		booLen   = [1]*len(booList)
		intList  = ['ns','nfp','mpol','ntor','mnmax','mnmax_nyq','iasym','ierr_vmec']
		intLen   = [1]*len(intList)
		realList = ['wb','wp','gamma','pfac','rmax_surf','rmin_surf','zmax_surf',\
			'aspect','betatot','betapol','betator','betaxis','b0','version_',\
			'ionlarmor','volavgb','fsql','fsqr','fsqz','ftolv','aminor','rmajor',\
			'volume','rbtor','rbtor0','itor','machsq']
		realLen = [1]*len(realList)
		scalar_data = self.get_module_vars(module_name,booList,booLen,intList,intLen,realList,realLen)
		ns = scalar_data['ns']
		# Get 1D Real Arrays
		realList = ['iotas','iotaf','presf','phipf','chipf','chi','phi','mass',\
			'pres','beta_vol','phip','buco','bvco','vp','overr','jcuru',\
			'jcurv','specw','jdotb','dmerc','dwell','dcurr','dgeod','equif']
		realLen = [(scalar_data['ns'],1)]*len(realList)
		realList.extend(['xm','xn','xm_nyq','xn_nyq'])
		realLen.extend([(scalar_data['mnmax'],1)]*2)
		realLen.extend([(scalar_data['mnmax_nyq'],1)]*2)
		# Add 2D Arrays
		realList.extend(['rmnc','zmns','lmns'])
		realLen.extend([(scalar_data['ns'],scalar_data['mnmax'])]*3)
		realList.extend([ 'bmnc','gmnc','bsupumnc','bsupvmnc',\
			'bsubsmns','bsubumnc','bsubvmnc',\
			'currumnc','currvmnc'])
		realLen.extend([(scalar_data['ns'],scalar_data['mnmax_nyq'])]*9)
		if scalar_data['iasym']:
			realList.extend(['rmns','zmnc','lmnc'])
			realLen.extend([(scalar_data['ns'],scalar_data['mnmax'])]*3)
			realList.extend([ 'bmns','gmns','bsupumns','bsupvmns',\
				'bsubsmnc','bsubumns','bsubvmns',\
				'currumns','currvmns'])
			realLen.extend([(scalar_data['ns'],scalar_data['mnmax_nyq'])]*9)
		array_data = self.get_module_vars(module_name,realVar=realList,realLen=realLen)
		# Try reading strings
		charVar=['mgrid_file','input_extension','pmass_type','pcurr_type','piota_type']
		charLen=[(200,1),(100,1),(20,1),(20,1),(20,1)]
		string_data = self.get_module_vars(module_name,charVar=charVar,charLen=charLen,ldefined_size_arrays=True)
		# Return
		return scalar_data | array_data | string_data

	def get_module_vars(self,modName,booVar=None,booLen=None,\
		intVar=None,intLen=None,realVar=None,realLen=None,\
		charVar=None,charLen=None,\
		ldefined_size_arrays=False):
		"""Reads module variables into a python dictionary

		This rountine helps to streamline reading of module variables
		into a python dictionary. The variable lists (booVar, intVar,
		and realVar) are lists of strings referencing the variable names.
		The length lists (booLen,intLen,realLen) are a list of tupules
		defining the size of the variable.

		Parameters
		----------
		modName : str
			Name of the module eg 'read_wout_mod'
		booVar : list (optional)
			List of strings referencing boolean variables to pull
		booLen : list (optional)
			List of tuples defining array size
		intVar : list (optional)
			List of strings referencing integer variables to pull
		intLen : list (optional)
			List of tuples defining array size
		realVar : list (optional)
			List of strings referencing integer variables to pull
		realLen : list (optional)
			List of tuples defining array size
		ldefined_size_arrays : logical (optional)
			Set to true when reading a predefined size arrays temp(0:20)
		Returns
		-------
		out_data : dict
			Dictionary of module variables
		"""
		import ctypes as ct
		import numpy.ctypeslib as npct
		from math import prod
		out_data={}
		# Booleans
		if booVar:
			ftemp = ct.POINTER(ct.c_bool)
			for i,temp in enumerate(booVar):
				if booLen[i]==1:
					out_data[temp]=ct.c_bool.in_dll(self.libstell,modName+'_'+temp+self.s3).value
				else:
					if ldefined_size_arrays : ftemp=ct.c_bool*prod(booLen[i])
					out_data[temp]=npct.as_array(ftemp.in_dll(self.libstell,modName+'_'+temp+self.s3),booLen[i])
		# Integers
		if intVar:
			ftemp = ct.POINTER(ct.c_int)
			for i,temp in enumerate(intVar):
				if intLen[i]==1:
					out_data[temp]=ct.c_int.in_dll(self.libstell,modName+'_'+temp+self.s3).value
				else:
					if ldefined_size_arrays : ftemp=ct.c_int*prod(intLen[i])
					out_data[temp]=npct.as_array(ftemp.in_dll(self.libstell,modName+'_'+temp+self.s3),intLen[i])
		# Reals
		if realVar:
			ftemp = ct.POINTER(ct.c_double)
			for i,temp in enumerate(realVar):
				if realLen[i]==1:
					out_data[temp]=ct.c_double.in_dll(self.libstell,modName+'_'+temp+self.s3).value
				else:
					if ldefined_size_arrays : ftemp=ct.c_double*prod(realLen[i])
					out_data[temp]=npct.as_array(ftemp.in_dll(self.libstell,modName+'_'+temp+self.s3),realLen[i])
		# Reals
		if charVar:
			ftemp = ct.POINTER(ct.c_char)
			for i,temp in enumerate(charVar):
				if charLen[i]==1:
					out_data[temp]=ct.c_char.in_dll(self.libstell,modName+'_'+temp+self.s3).value.decode('UTF-8')
				else:
					if ldefined_size_arrays : ftemp=ct.c_char*prod(charLen[i])
					out_data[temp]=ftemp.in_dll(self.libstell,modName+'_'+temp+self.s3).value.decode('UTF-8')
		return out_data

	def set_module_var(self,modName,var,val):
		import ctypes as ct
		import numpy as np
		import numpy.ctypeslib as npct
		if type(val) == bool:
		    f = ct.c_bool
		elif type(val) == int:
		    f = ct.c_int
		elif type(val) == float:
		    f = ct.c_double
		elif type(val) == str:
		    n = len(val)
		    f = ct.c_char*n
		elif type(val) in (np.ndarray,list):
		    if type(val[0]) == bool:
		        tt = ct.c_bool
		    elif type(val[0]) == np.int32:
		        tt = ct.c_int
		    elif type(val[0]) == np.float64:
		        tt = ct.c_double
		    else:
		        print('   Unrecognized type:',type(val[0]))
		        return
		    n = val.ndim
		    f = tt*val.size
		else:
		    print('   Unrecognized type:',type(val))
		    return
		temp=f.in_dll(self.libstell,self.s1+''+modName+'_'+self.s2+'_'+var+self.s3)
		if type(val) == np.ndarray:
		    if n==1:
		        for i,col in enumerate(val):
		            temp[i].value = val[i]
		elif type(val) == str:
		    temp.value = val.encode('UTF-8')
		else:
		    temp.value = val
		return

	def pcurr(self,s):
		"""Wrapper to the PCURR function

		This routine wrappers the PCURR function which
		returns either the I'(s) or I(s) depending on the
		choice of pcurr_type.

		Parameters
		----------
		s : real
			Value of normalized toroidal flux.
		Returns
		-------
		val : real
			Value of I'(s) or I(s)
		"""
		import ctypes as ct
		# Load Libraries
		pcurr_func = getattr(self.libstell,'pcurr_')
		pcurr_func.argtypes = [ct.POINTER(ct.c_double)]
		pcurr_func.restype=ct.c_double
		s_temp = ct.c_double(s)
		val = pcurr_func(ct.byref(s_temp))
		return val

	def piota(self,s):
		"""Wrapper to the PIOTA function

		This routine wrappers the PIOTA function which
		returns the iota(s) rotational transform.

		Parameters
		----------
		s : real
			Value of normalized toroidal flux
		Returns
		-------
		val : real
			Value of iota(s)
		"""
		import ctypes as ct
		# Load Libraries
		pcurr_func = getattr(self.libstell,'piota_')
		pcurr_func.argtypes = [ct.POINTER(ct.c_double)]
		pcurr_func.restype=ct.c_double
		s_temp = ct.c_double(s)
		val = pcurr_func(ct.byref(s_temp))
		return val

	def pmass(self,s):
		"""Wrapper to the PMASS function

		This routine wrappers the PMASS function which
		returns the mass(s) pressure function.

		Parameters
		----------
		s : real
			Value of normalized toroidal flux
		Returns
		-------
		val : real
			Value of mass(s)
		"""
		import ctypes as ct
		# Load Libraries
		pcurr_func = getattr(self.libstell,'pmass_')
		pcurr_func.argtypes = [ct.POINTER(ct.c_double)]
		pcurr_func.restype=ct.c_double
		s_temp = ct.c_double(s)
		val = pcurr_func(ct.byref(s_temp))
		return val

	def parse_coils_file(self,filename):
		"""Parses a coils file

		This routine wrappers parsecoils in LIBSTELL biotsavart and returns
		a dictionary of values

		Parameters
		----------
		file : str
			Path to coils file.
		Returns
		----------
		vars : dict
			Dictionary of module variables
		"""
		import ctypes as ct
		module_name = self.s1+'biotsavart_'+self.s2
		parse_coils_file = getattr(self.libstell,module_name+'_parse_coils_file'+self.s3)
		parse_coils_file.argtypes=[ct.c_char_p, ct.c_bool, ct.c_long]
		parse_coils_file.restype=None
		lgrps = ct.c_bool(False)
		parse_coils_file(file.encode('UTF-8'), ct.byref(lgrps), len(file))
		if not (ierr.value == 0):
			return None

class FourierRep():
	def __init__(self, parent=None):
		self.test = None

	def cfunct(self,theta,phi,fmnc,xm,xn):
		"""Cos transformation

		This routine performs the cosine transformation of the
		kernel (m*theta+n*phi).

		Parameters
		----------
		theta : ndarray
			Poloidal grid in radians.
		phi : ndarray
			Toroidal grid in radians.
		fmnc : ndarray
			Array to transform (radial,fourier)
		xm : ndarray
			Poloidal harmonic array.
		xn : ndarray
			Toroidal harmonic array.
		Returns
		----------
		f : ndarray
			Fourier transformed quantity (radial,poloidal,toroidal)
		"""
		import numpy as np
		f=0
		(ns,mn)=fmnc.shape
		lt = len(theta)
		lz = len(phi)
		mt=np.matmul(xm,theta.T)
		nz=np.matmul(xn,phi.T)
		cosmt=np.cos(mt)
		sinmt=np.sin(mt)
		cosnz=np.cos(nz)
		sinnz=np.sin(nz)
		f = np.zeros((ns,lt,lz))
		fmn = np.ndarray((mn,lt))
		for k in range(ns):
			fmn = np.broadcast_to(fmnc[k,:],(lt,mn)).T
			f[k,:,:]=np.matmul((fmn*cosmt).T, cosnz)-np.matmul((fmn*sinmt).T, sinnz)
		return f

	def sfunct(self,theta,phi,fmnc,xm,xn):
		"""Sine transformation

		This routine performs the sine transformation of the
		kernel (m*theta+n*phi).

		Parameters
		----------
		theta : ndarray
			Poloidal grid in radians.
		phi : ndarray
			Toroidal grid in radians.
		fmnc : ndarray
			Array to transform (radial,fourier)
		xm : ndarray
			Poloidal harmonic array.
		xn : ndarray
			Toroidal harmonic array.
		Returns
		----------
		f : ndarray
			Fourier transformed quantity (radial,poloidal,toroidal)
		"""
		import numpy as np
		f=0
		(ns,mn)=fmnc.shape
		lt = len(theta)
		lz = len(phi)
		mt=np.matmul(xm,theta.T)
		nz=np.matmul(xn,phi.T)
		cosmt=np.cos(mt)
		sinmt=np.sin(mt)
		cosnz=np.cos(nz)
		sinnz=np.sin(nz)
		f = np.zeros((ns,lt,lz))
		fmn = np.ndarray((mn,lt))
		for k in range(ns):
			fmn = np.broadcast_to(fmnc[k,:],(lt,mn)).T
			f[k,:,:]=np.matmul((fmn*sinmt).T, cosnz)-np.matmul((fmn*cosmt).T, sinnz)
		return f
		
	def isotoro(self,r,z,zeta,svals,*args,**kwargs):
		import numpy as np
		import matplotlib.pyplot as pyplot
		import mpl_toolkits.mplot3d as mplot3d
		import math as math
		import matplotlib.tri as mtri
		#from mayavi import mlab
		nr = np.size(svals)
		if (nr == 1):
			s= [svals]
			nr = 1
		else:
			s=svals
		nt = np.size(r,1)
		nz = np.size(r,2)
		vertex = np.zeros((nt*nz,3,nr))
		for k in range(0,nr):
			ivertex = 0
			ifaces = 0
			for j in range(0,nz):
				for i in range(0,nt):
					vertex[ivertex,0,k]=r[s[k],i,j]*math.cos(zeta[j])
					vertex[ivertex,1,k]=r[s[k],i,j]*math.sin(zeta[j])
					vertex[ivertex,2,k]=z[s[k],i,j]
					ivertex = ivertex + 1
		u = np.linspace(0, 1, endpoint=True, num=nt)
		v = np.linspace(0, 1, endpoint=True, num=nz)
		u, v = np.meshgrid(u, v)
		u, v = u.flatten(), v.flatten()
		tri = mtri.Triangulation(u, v)
		test=len(kwargs)
		fig=kwargs.pop('fig',pyplot.figure())
		h=kwargs.pop('axes',fig.add_subplot(111,projection='3d'))
		for k in range(0,nr):
			if (len(args)==0):
				tsurf=h.plot_trisurf(vertex[:,0,k],vertex[:,1,k],vertex[:,2,k], triangles=tri.triangles,color='red',shade='yes',linewidth=0.0,alpha=1)
				#tsurf=mlab.triangular_mesh(vertex[:,0,k],vertex[:,1,k],vertex[:,2,k], tri.triangless)
			else:
				# Matplotlib way (SLOW)
				vals = args[0][s[k],:,:].T.flatten()
				colors = np.mean(vals[tri.triangles], axis=1)
				tsurf=h.plot_trisurf(vertex[:,0,k],vertex[:,1,k],vertex[:,2,k], triangles=tri.triangles,cmap='jet',shade=False,linewidth=0.0,alpha=1)
				tsurf.set_array(colors)
				tsurf.autoscale()
				#MAYAVI Way (need to figure out how to embed)
				#h    = mlab.figure()
				#vals = args[0][s[k],:,:].T.flatten()
				#tsurf=mlab.triangular_mesh(vertex[:,0,k],vertex[:,1,k],vertex[:,2,k], tri.triangles, scalars=vals, colormap='jet',figure=h)
				#print(type(tsurf))
		if (test==0):
			pyplot.show()
		return h



# Main routine
if __name__=="__main__":
	import sys
	temp = LIBSTELL()
	val=temp.read_indata('input.ORBITS')
	temp.set_module_var('vmec_input','mpol',12)
	temp.set_module_var('vmec_input','niter_array',[1,2,3,4,5])
	temp.write_indata('input.test')
	#wout=temp.read_wout('wout_W7X_AIM_n04_e30_i15_8SH2_slow.nc')
	sys.exit(0)