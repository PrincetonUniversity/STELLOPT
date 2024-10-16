##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for interfacing to libstell
"""

# Libraries

# Constants

# LIBSTELL Class
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

	def vtkColor(self,j):
		"""Helper to vary colors

		This routine returns a VTK color object based on a number.
		This is done to mimic the behavior of changing colors when
		plotting lines in matplotlib.

		Parameters
		----------
		j : int
			Index into predefined colors
		"""
		from vtkmodules.vtkCommonColor import vtkNamedColors
		colors = vtkNamedColors()
		color_txt=['red','green','blue','yellow','magenta','cyan','aqua']
		j = j % len(color_txt)
		return colors.GetColor3d(color_txt[j])


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
		import numpy as np
		# Get constants
		module_name = self.s1+'vsvd0_'+self.s2
		get_constant = getattr(self.libstell,module_name+'_getnigroup'+self.s3)
		get_constant.argtypes = None
		get_constant.restype=ct.c_int
		nigroup = get_constant()
		module_name = self.s1+'vparams_'+self.s2
		get_constant = getattr(self.libstell,module_name+'_getndatafmax'+self.s3)
		get_constant.argtypes = None
		get_constant.restype=ct.c_int
		ndatafmax = get_constant()
		get_constant = getattr(self.libstell,module_name+'_getmpol1d'+self.s3)
		get_constant.argtypes = None
		get_constant.restype=ct.c_int
		mpol1d = get_constant()
		get_constant = getattr(self.libstell,module_name+'_getntord'+self.s3)
		get_constant.argtypes = None
		get_constant.restype=ct.c_int
		ntord = get_constant()
		get_constant = getattr(self.libstell,module_name+'_getnsd'+self.s3)
		get_constant.argtypes = None
		get_constant.restype=ct.c_int
		nsd = get_constant()
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
			'imse','isnodes','itse','ipnodes','iopt_raxis','imatch_phiedge','nflxs','pre_niter']
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
		realLen.extend([(100,1),(nigroup,1)])
		charList=['pcurr_type','piota_type','pmass_type','pt_type','ph_type']
		charLen=[(20,1)]*6
		charList.extend(['mgrid_file','input_extension'])
		charLen.extend([(200,1),(200,1)])
		out_data = self.get_module_vars(module_name,booList,booLen,intList,intLen,realList,realLen,charList,charLen,ldefined_size_arrays=True)
		# Note we need to reshape the rbc/s zbs/c arrays.
		out_data['rbc'] = np.reshape(out_data['rbc'],(mpol1d+1,2*ntord+1))
		out_data['zbc'] = np.reshape(out_data['zbc'],(mpol1d+1,2*ntord+1))
		out_data['rbs'] = np.reshape(out_data['rbs'],(mpol1d+1,2*ntord+1))
		out_data['zbs'] = np.reshape(out_data['zbs'],(mpol1d+1,2*ntord+1))
		return out_data

	def write_indata(self,filename,out_dict=None):
		"""Wrappers writing of the VMEC INDATA namelist

		This routine wrappers write_indata_namelist in LIBSTELL

		Parameters
		----------
		file : str
			Path to input file.
		out_dict : dict (optional)
			Dictionary of items to change.
		"""
		import ctypes as ct
		module_name = self.s1+'vmec_input_'+self.s2
		# Check if we want to update values
		if out_dict:
			for key in out_dict:
				self.set_module_var(module_name,key,out_dict[key])
		write_indata_namelist = getattr(self.libstell,module_name+'_write_indata_namelist_byfile'+self.s3)
		write_indata_namelist.argtypes = [ct.c_char_p,ct.c_long]
		write_indata_namelist.restype=None
		write_indata_namelist(filename.encode('UTF-8'),len(filename))

	def read_bootin(self,filename):
		"""Reads a BOOTSJ BOOTIN namelist

		This routine wrappers read_boot_namelist function in
		bootsj_input. To do this we wrapper safe_open_mod

		Parameters
		----------
		file : str
			Path to wout file.
		"""
		import ctypes as ct
		# We use an added routine as a helper
		module_name = self.s1+'bootsj_input_'+self.s2
		read_bootin_namelist = getattr(self.libstell,module_name+'_read_boot_namelist_byfile'+self.s3)
		read_bootin_namelist.argtypes = [ct.c_char_p,ct.c_long]
		read_bootin_namelist.restype=None
		read_bootin_namelist(filename.encode('UTF-8'),len(filename))
		# Get vars
		intList=['nrho','mbuse','nbuse','isymm0']
		intLen=[1]*len(intList)
		realList=['tempres', 'zeff1', 'dens0', 'teti', 'damp', 'damp_bs']
		realLen=[1]*len(realList)
		realList.extend(['ate','ati'])
		realLen.extend([(12,1)]*2)
		out_data = self.get_module_vars(module_name,intVar=intList,intLen=intLen,realVar=realList,realLen=realLen,ldefined_size_arrays=True)
		return out_data

	def set_bootin(self,in_dict):
		"""Set the bootin namelist variables via a passed dict

		This routine passes variable information back to the module
		for writing using a dict
		'rho':		number of rho values to use (depricated)
		'mbuse':	number of poloidal modes in B-field to use
		'nbuse':	number of toroidal modes in B-field to use
		'zeff1':	Effective ion charge
		'dens0':	Central electron density in 1E20 m^-3 units
		'teti':		Electron to ion temeprature ratio
		'tempres':	Te(s) = P(s)**tempres or if -1 Te(s) = sqrt(p)
		'damp':     Superceeded by damp_bs
		'damp_bs':  Resonance damping factor
		'isymm0':   if !=0 force symmetry
		'ate':      Polynomial electron temperature
		'ati':      Polynomial ion temperature

		Parameters
		----------
		in_dict : dict
			Input dictionary
		"""
		import ctypes as ct
		module_name = self.s1+'bootsj_input_'+self.s2
		for key in in_dict.keys():
			self.set_module_var(module_name,key,in_dict[key])


	def write_bootin(self,filename):
		"""Wrappers writing of the BOOTSJ BOOTIN namelist

		This routine wrappers write_bootsj_input in LIBSTELL

		Parameters
		----------
		file : str
			Path to input file.
		"""
		import ctypes as ct
		module_name = self.s1+'bootsj_input_'+self.s2
		write_bootin_namelist = getattr(self.libstell,module_name+'_write_boot_namelist_byfile'+self.s3)
		write_bootin_namelist.argtypes = [ct.c_char_p,ct.c_long]
		write_bootin_namelist.restype=None
		write_bootin_namelist(filename.encode('UTF-8'),len(filename))

	def read_beams3d_input(self,filename):
		"""Reads a BEAMS3D BEAMS3D_INPUT namelist

		This routine wrappers read_beams3d_input function in
		beams3d_input_mod.

		Parameters
		----------
		file : str
			Path to wout file.
		"""
		import ctypes as ct
		# A few constants defined in globals
		module_name = self.s1+'beams3d_globals_'+self.s2
		get_constant = getattr(self.libstell,module_name+'_getmaxparticles'+self.s3)
		get_constant.argtypes = None
		get_constant.restype=ct.c_int
		maxparticles = get_constant()
		get_constant = getattr(self.libstell,module_name+'_getmaxbeams'+self.s3)
		get_constant.argtypes = None
		get_constant.restype=ct.c_int
		maxbeams = get_constant()
		get_constant = getattr(self.libstell,module_name+'_getmaxproflen'+self.s3)
		get_constant.argtypes = None
		get_constant.restype=ct.c_int
		maxproflen = get_constant()
		get_constant = getattr(self.libstell,module_name+'_getnion'+self.s3)
		get_constant.argtypes = None
		get_constant.restype=ct.c_int
		nion = get_constant()
		# Call the initialization routine
		module_name = self.s1+'beams3d_input_mod_'+self.s2
		init_beams3d_input = getattr(self.libstell,module_name+'_init_beams3d_input'+self.s3)
		init_beams3d_input.argtypes = None
		init_beams3d_input.restype = None
		init_beams3d_input()
		# Only read a file if we didn't pass an empty string.
		if filename != '':
			# We use an added routine as a helper
			module_name = self.s1+'beams3d_input_mod_'+self.s2
			read_beams3d_input = getattr(self.libstell,module_name+'_read_beams3d_input'+self.s3)
			read_beams3d_input.argtypes = [ct.c_char_p,ct.POINTER(ct.c_int),ct.c_long]
			read_beams3d_input.restype = None
			istat = ct.c_int(0)
			read_beams3d_input(filename.encode('UTF-8'),ct.byref(istat),len(filename))
			if not (istat.value == 0):
				return None
		# Get vars
		intList=['nr','nphi','nz','nparticles_start','npoinc', \
				'duplicate_factor', 'ns_prof1','ns_prof2', \
				'ns_prof3','ns_prof4','ns_prof5', 'nr_fida', \
				'nphi_fida', 'nz_fida', 'nenergy_fida', 'npitch_fida']
		intLen=[1]*len(intList)
		intList.extend(['dex_beams','ni_aux_z'])
		intLen.extend([(maxbeams,1),(nion,1)])
		realList=['rmin', 'rmax', 'zmin', 'zmax', 'phimin', 'phimax', \
				  'follow_tol', 'vc_adapt_tol', \
				  'ne_scale', 'te_scale', 'ti_scale', 'zeff_scale', \
				  'fusion_scale', 'plasma_zmean', 'plasma_mass', \
				  'therm_factor', 'lendt_m', 'te_col_min', \
				  'b_kick_min', 'b_kick_max', 'freq_kick', 'e_kick', \
				  'rho_fullorbit', 'partvmax', 'rmin_fida', \
				  'rmax_fida', 'zmin_fida', 'zmax_fida', 'phimin_fida', \
				  'phimax_fida', 't_fida', 'rho_max_dist']
		realLen=[1]*len(realList)
		realList.extend(['r_start_in', 'phi_start_in', 'z_start_in',\
						 'vll_start_in', 'mu_start_in', 't_end_in', \
						 'charge_in', 'mass_in', 'zatom_in', \
						 'weight_in', 'vr_start_in', 'vphi_start_in', \
						 'vz_start_in'])
		realLen.extend([(maxparticles,1)]*13)
		realList.extend(['mass_beams', 'charge_beams', 'e_beams',\
						 'p_beams', 'div_beams', 'asize_beams', \
						 'adist_beams'])
		realLen.extend([(maxbeams,1)]*7)
		realList.extend(['r_beams', 'z_beams', 'phi_beams'])
		realLen.extend([(maxbeams,2)]*3)
		realList.extend(['ne_aux_s', 'te_aux_s', 'ti_aux_s', \
						 'ne_aux_f', 'te_aux_f', 'ti_aux_f', \
						 'pot_aux_s', 'pot_aux_f', 'zeff_aux_s', \
						 'zeff_aux_f'])
		realLen.extend([(maxproflen,1)]*10)
		realList.extend(['ni_aux_s', 'ni_aux_f', 'ni_aux_m'])
		realLen.extend([(maxproflen,1),(nion,maxproflen),(nion,1)])
		module_name = self.s1+'beams3d_globals_'+self.s2
		out_data = self.get_module_vars(module_name,intVar=intList,intLen=intLen,realVar=realList,realLen=realLen,ldefined_size_arrays=True)
		return out_data

	def write_beams3d_input(self,filename,out_dict=None):
		"""Wrappers writing of the BEAMS3D_INPUT namelist

		This routine wrappers write_beams3d_input in LIBSTELL

		Parameters
		----------
		file : str
			Path to input file.
		out_dict : dict (optional)
			Dictionary of items to change.
		"""
		import ctypes as ct
		module_name = self.s1+'beams3d_globals_'+self.s2
		# Check if we want to update values
		if out_dict:
			for key in out_dict:
				self.set_module_var(module_name,key,out_dict[key])
		module_name = self.s1+'beams3d_input_mod_'+self.s2
		write_beams3d_namelist = getattr(self.libstell,module_name+'_write_beams3d_namelist_byfile'+self.s3)
		write_beams3d_namelist.argtypes = [ct.c_char_p,ct.c_long]
		write_beams3d_namelist.restype=None
		write_beams3d_namelist(filename.encode('UTF-8'),len(filename))

	def read_diagno_in(self,filename):
		"""Reads a DIAGNO_IN namelist

		This routine wrappers read_diagno_input function in
		diagno_input_mod.

		Parameters
		----------
		filename : str
			Path to input file.
		Returns
		-------
		out_data : dict
			Dictionary of items.
		"""
		import ctypes as ct
		# We use an added routine as a helper
		module_name = self.s1+'diagno_input_mod_'+self.s2
		read_diagno_input = getattr(self.libstell,module_name+'_read_diagno_input'+self.s3)
		read_diagno_input.argtypes = [ct.c_char_p, ct.POINTER(ct.c_int), ct.c_long]
		read_diagno_input.restype=None
		ierr = ct.c_int(0)
		read_diagno_input(filename.encode('UTF-8'),ct.byref(ierr),len(filename))
		# Get vars
		booList=['lrphiz','luse_mut','lvc_field','lapoints_accurate_output','lbpoints_accurate_output']
		booLen=[1]*len(booList)
		booList.extend(['luse_extcur'])
		booLen.extend([(512,1)])
		intList=['nu','nv']
		intLen=[1]*len(intList)
		realList=['units','vc_adapt_tol','vc_adapt_rel','int_step']
		realLen=[1]*len(realList)
		realList.extend(['bprobe_turns','flux_turns','segrog_turns'])
		realLen.extend([(2048,1),(512,1),(256,1)])
		charList=['int_type','afield_points_file','bfield_points_file',\
				'bprobes_file','mirnov_file', 'seg_rog_file','flux_diag_file',\
				'bprobes_mut_file', 'mir_mut_file', 'rog_mut_file', 'flux_mut_file']
		charLen=[(256,1)]*11
		module_name = self.s1+'diagno_runtime_'+self.s2
		out_data = self.get_module_vars(module_name,booList,booLen,intList,intLen,realList,realLen,charList,charLen,ldefined_size_arrays=True)
		return out_data

	def write_diagno_in(self,filename,out_dict=None):
		"""Wrappers writing of the DIAGNO_IN namelist

		This routine wrappers write_diagno_input in LIBSTELL

		Parameters
		----------
		filename : str
			Path to input file.
		out_dict : dict (optional)
			Dictionary of items to change.
		"""
		import ctypes as ct
		module_name = self.s1+'diagno_runtime_'+self.s2
		# Check if we want to update values
		if out_dict:
			for key in out_dict:
				self.set_module_var(module_name,key,out_dict[key])
		module_name = self.s1+'diagno_input_mod_'+self.s2
		write_diagno_input = getattr(self.libstell,module_name+'_write_diagno_input_byfile'+self.s3)
		write_diagno_input.argtypes = [ct.c_char_p, ct.c_long]
		write_diagno_input.restype=None
		write_diagno_input(filename.encode('UTF-8'),len(filename))

	def read_fieldlines_input(self,filename):
		"""Reads a FIELDLINES_INPUT namelist

		This routine wrappers read_fieldlines_input function in
		fieldlines_input_mod.

		Parameters
		----------
		filename : str
			Path to input file.
		Returns
		-------
		out_data : dict
			Dictionary of items.
		"""
		import ctypes as ct
		# A few constants defined in globals
		module_name = self.s1+'fieldlines_globals_'+self.s2
		get_constant = getattr(self.libstell,module_name+'_getmaxlines'+self.s3)
		get_constant.argtypes = None
		get_constant.restype=ct.c_int
		maxlines = get_constant()
		# We use an added routine as a helper
		module_name = self.s1+'fieldlines_input_mod_'+self.s2
		read_fieldlines_input = getattr(self.libstell,module_name+'_read_fieldlines_input'+self.s3)
		read_fieldlines_input.argtypes = [ct.c_char_p, ct.POINTER(ct.c_int), ct.c_long]
		read_fieldlines_input.restype=None
		ierr = ct.c_int(0)
		read_fieldlines_input(filename.encode('UTF-8'),ct.byref(ierr),len(filename))
		# Get vars
		intList=['nr','nphi','nz','npoinc','num_hcp']
		intLen=[1]*len(intList)
		realList=['rmin','rmax','zmin','zmax','phimin','phimax','vc_adapt_tol','follow_tol','mu','delta_hc']
		realLen=[1]*len(realList)
		realList.extend(['r_start','z_start','phi_start','phi_end','r_hc','z_hc','phi_hc'])
		realLen.extend([(maxlines,1)]*7)
		realList.extend(['errorfield_amp','errorfield_phase'])
		realLen.extend([(20,1)]*2)
		charList=['int_type']
		charLen=[(256,1)]
		module_name = self.s1+'fieldlines_globals_'+self.s2
		out_data = self.get_module_vars(module_name,None,None,intList,intLen,realList,realLen,charList,charLen,ldefined_size_arrays=True)
		return out_data

	def write_fieldlines_input(self,filename,out_dict=None):
		"""Wrappers writing of the FIELDLINES_INPUT namelist

		This routine wrappers write_fieldlines_input in LIBSTELL

		Parameters
		----------
		filename : str
			Path to input file.
		out_dict : dict (optional)
			Dictionary of items to change.
		"""
		import ctypes as ct
		module_name = self.s1+'fieldlines_globals_'+self.s2
		# Check if we want to update values
		if out_dict:
			for key in out_dict:
				self.set_module_var(module_name,key,out_dict[key])
		module_name = self.s1+'fieldlines_input_mod_'+self.s2
		write_fieldlines_input = getattr(self.libstell,module_name+'_write_fieldlines_namelist_byfile'+self.s3)
		write_fieldlines_input.argtypes = [ct.c_char_p, ct.c_long]
		write_fieldlines_input.restype=None
		write_fieldlines_input(filename.encode('UTF-8'),len(filename))

	def read_nescoil_input(self,filename):
		"""Reads a NESCOIL input file

		This routine wrappers read_nescin function in
		read_nescoil_mod.

		Parameters
		----------
		filename : str
			Path to input file.
		Returns
		-------
		out_data : dict
			Dictionary of items.
		"""
		import ctypes as ct
		# We use an added routine as a helper
		module_name = self.s1+'read_nescoil_mod_'+self.s2
		read_nescin = getattr(self.libstell,module_name+'_read_nescin'+self.s3)
		read_nescin.argtypes=[ct.c_char_p, ct.POINTER(ct.c_int), ct.c_long]
		read_nescin.restype=None
		ierr = ct.c_int(0)
		read_nescin(filename.encode('UTF-8'), ct.byref(ierr), len(filename))
		# Setup Arrays
		out_data={}
		# Get Scalars
		booList  = ['lasym']
		booLen   = [1]*len(booList)
		intList  = ['nu', 'nv', 'nu1', 'nv1', 'mpol', 'ntor', 'mf', 'nf', 'md', 'nd', 'np', \
					'ibex', 'mstrt', 'mstep', 'mkeep', 'mdspw', 'w_psurf', 'w_csurf', \
					'w_bnuv', 'w_jsurf', 'w_xerr', 'w_svd', 'mnmax_plasma', \
					'mnmax_surface', 'nmax', 'mnd', 'nuv', 'nuv1', 'nuvh', 'nuvh1']
		intLen   = [1]*len(intList)
		realList = ['iota_edge', 'phip_edge', 'curpol', 'cut', 'cup', 'curwt', 'trgwt','alp']
		realLen = [1]*len(realList)
		scalar_data = self.get_module_vars(module_name,booList,booLen,intList,intLen,realList,realLen)
		# Get 1D Int Arrays
		intList= ['xm_plasma', 'xn_plasma','xm_surface', 'xn_surface']
		intLen = [(scalar_data['mnmax_plasma'],1),(scalar_data['mnmax_plasma'],1),\
			(scalar_data['mnmax_surface'],1),(scalar_data['mnmax_surface'],1)]
		# Get 1D Real Arrays
		realList = ['rmnc_plasma', 'zmns_plasma', 'rmns_plasma', \
			'zmnc_plasma', 'lmnc_plasma', 'lmns_plasma' ]
		realLen = [(scalar_data['mnmax_plasma'],1)]*len(realList)
		realList.extend(['rmnc_surface', 'zmns_surface', 'rmns_surface', \
			'zmnc_surface'])
		realLen.extend([(scalar_data['mnmax_surface'],1)]*4)
		array_data = self.get_module_vars(module_name,intVar=intList,intLen=intLen,realVar=realList,realLen=realLen)
		# Return
		return scalar_data | array_data

	def write_nescoil_input(self,filename,out_dict=None):
		"""Wrappers writing of the NESCOIL nescin file

		This routine wrappers write_nescin in LIBSTELL

		Parameters
		----------
		filename : str
			Path to input file.
		out_dict : dict (optional)
			Dictionary of items to change.
		"""
		import ctypes as ct
		module_name = self.s1+'read_nescoil_mod_'+self.s2
		# Check if we want to update values
		if out_dict:
			for key in out_dict:
				self.set_module_var(module_name,key,out_dict[key])
		write_nescin = getattr(self.libstell,module_name+'_write_nescin'+self.s3)
		write_nescin.argtypes = [ct.c_char_p, ct.POINTER(ct.c_int), ct.c_long]
		write_nescin.restype=None
		ierr = ct.c_int(0)
		write_nescin(filename.encode('UTF-8'),ct.byref(ierr),len(filename))

	def read_stellopt_input(self,filename):
		"""Reads a STELLOPT OPTIMUM namelist

		This routine wrappers read_stellopt_input function in
		stellopt_input_mod.

		Parameters
		----------
		filename : str
			Path to input file.
		Returns
		-------
		out_data : dict
			Dictionary of items.
		"""
		import ctypes as ct
		import numpy as np
		# A few constants defined in globals
		module_name = self.s1+'vsvd0_'+self.s2
		get_constant = getattr(self.libstell,module_name+'_getnigroup'+self.s3)
		get_constant.argtypes = None
		get_constant.restype=ct.c_int
		nigroup = get_constant()
		module_name = self.s1+'vparams_'+self.s2
		get_constant = getattr(self.libstell,module_name+'_getndatafmax'+self.s3)
		get_constant.argtypes = None
		get_constant.restype=ct.c_int
		ndatafmax = get_constant()
		get_constant = getattr(self.libstell,module_name+'_getmpol1d'+self.s3)
		get_constant.argtypes = None
		get_constant.restype=ct.c_int
		mpol1d = get_constant()
		get_constant = getattr(self.libstell,module_name+'_getntord'+self.s3)
		get_constant.argtypes = None
		get_constant.restype=ct.c_int
		ntord = get_constant()
		get_constant = getattr(self.libstell,module_name+'_getnsd'+self.s3)
		get_constant.argtypes = None
		get_constant.restype=ct.c_int
		nsd = get_constant()
		module_name = self.s1+'stellopt_globals_'+self.s2
		get_constant = getattr(self.libstell,module_name+'_getmaxwindsurf'+self.s3)
		get_constant.argtypes = None
		get_constant.restype=ct.c_int
		maxwindsurf = get_constant()
		get_constant = getattr(self.libstell,module_name+'_getbigno'+self.s3)
		get_constant.argtypes = None
		get_constant.restype=ct.c_int
		bigno = get_constant()
		# Call the initialization routine
		module_name = self.s1+'stellopt_input_mod_'+self.s2
		init_stellopt_input = getattr(self.libstell,module_name+'_init_stellopt_input'+self.s3)
		init_stellopt_input.argtypes = None
		init_stellopt_input.restype = None
		init_stellopt_input() 
		# We use an added routine as a helper
		module_name = self.s1+'stellopt_input_mod_'+self.s2
		read_stellopt_input = getattr(self.libstell,module_name+'_read_stellopt_input'+self.s3)
		read_stellopt_input.argtypes = [ct.c_char_p,ct.POINTER(ct.c_int),ct.c_long]
		read_stellopt_input.restype = None
		istat = ct.c_int(0)
		read_stellopt_input(filename.encode('UTF-8'),ct.byref(istat),len(filename))
		if not (istat.value == 0):
			return None
		# Get vars Globals
		module_name = self.s1+'stellopt_globals_'+self.s2
		booList=['lcentered_differences', 'lkeep_mins', 'lrefit', 'lcoil_geom', 'lno_restart', 'ltriangulate']
		booLen=[1]*len(booList)
		intList=['nfunc_max','cr_strategy', 'npopulation', 'noptimizers', 'mode', 'rho_exp']
		intLen=[1]*len(intList)
		realList=['ftol', 'xtol', 'gtol', 'epsfcn', 'factor', 'refit_param']
		realLen=[1]*len(realList)
		charList=['opt_type', 'axis_init_option']
		charLen=[(256,1),(256,1)]
		global_data = self.get_module_vars(module_name,booList,booLen,intList,intLen,realList,realLen,charList,charLen,ldefined_size_arrays=True)
		# Get VARS
		module_name = self.s1+'stellopt_vars_'+self.s2
		booList=['lphiedge_opt', 'lcurtor_opt', 'lpscale_opt', \
			'lbcrit_opt', 'lmix_ece_opt', 'lregcoil_winding_surface_separation_opt',\
			 'lregcoil_current_density_opt', 'lxval_opt', 'lyval_opt', \
			 'lxics_v0_opt','mango_bound_constraints']
		booLen=[1]*len(booList)
		booList.extend(['lextcur_opt','laphi_opt', 'lam_opt', \
					'lac_opt', 'lai_opt','lah_opt', 'lat_opt','lne_opt', \
					'lte_opt', 'lti_opt','lth_opt', 'lzeff_opt'])
		booLen.extend([(nigroup,1),(20,1)])
		booLen.extend([(21,1)]*10)
		booList.extend(['lam_s_opt', 'lam_f_opt', 'lac_s_opt', 'lac_f_opt', 'lai_s_opt', 'lai_f_opt',\
					'lne_f_opt', 'lte_f_opt', 'lti_f_opt', 'lth_f_opt', 'lphi_s_opt', 'lphi_f_opt', \
					'lzeff_f_opt', 'lemis_xics_f_opt', 'lbootj_f_opt', 'lbeamj_f_opt', 'lah_f_opt', 'lat_f_opt'])
		booLen.extend([(ndatafmax,1)]*18)
		booList.extend(['laxis_opt','lbound_opt','lrho_opt','lmode_opt','ldeltamn_opt'])
		booLen.extend([(ntord+1,1),(2*ntord+1,mpol1d+1),(2*ntord+1,mpol1d+1),(2*ntord+1,mpol1d+1),(2*ntord+1,2*mpol1d+1)])
		booList.extend(['lcoil_spline','lwindsurf'])
		booLen.extend([(nigroup,40),(maxwindsurf,1)])
		booList.extend(['lregcoil_rcws_rbound_c_opt','lregcoil_rcws_rbound_s_opt',\
			'lregcoil_rcws_zbound_c_opt','lregcoil_rcws_zbound_s_opt'])
		booLen.extend([(65,65)]*4)
		intList=['regcoil_nlambda', 'regcoil_num_field_periods', \
			'sfincs_min_procs', 'vboot_max_iterations']
		booList.extend(['lrosenbrock_x_opt'])
		booLen.extend([(20,1)])
		intLen=[1]*len(intList)
		intList.extend(['coil_nctrl'])
		intLen.extend([(nigroup,1)])
		realList=['dphiedge_opt', 'dcurtor_opt', 'dbcrit_opt', \
			'dpscale_opt', 'dmix_ece_opt', 'dxval_opt', 'dyval_opt', \
			'dregcoil_winding_surface_separation_opt', \
			'dregcoil_current_density_opt', 'dxics_v0_opt', \
			'phiedge_min', 'curtor_min', 'bcrit_min', \
			'pscale_min', 'mix_ece_min', 'xval_min', 'yval_min', 
			'regcoil_winding_surface_separation_min', \
			'regcoil_current_density_min', 'xics_v0_min', \
			'phiedge_max', 'curtor_max', 'bcrit_max', \
			'pscale_max', 'mix_ece_max', 'xval_max', 'yval_max', \
			'regcoil_winding_surface_separation_max', \
			'regcoil_current_density_max', 'xics_v0_max', \
			'mix_ece', 'xval', 'yval', 'xics_v0', \
			'regcoil_winding_surface_separation', \
			'regcoil_current_density','vboot_tolerance']
		realLen=[1]*len(realList)
		realList.extend(['dextcur_opt','extcur_min','extcur_max'])
		realLen.extend([(nigroup,1)]*3)
		realList.extend(['daphi_opt', 'aphi_min', 'aphi_max'])
		realLen.extend([(20,1)]*3)
		realList.extend(['dam_opt', 'dac_opt', 'dai_opt', 'dah_opt', 'dat_opt','dte_opt', 'dne_opt', 'dti_opt', 'dth_opt','dzeff_opt', \
			'am_min', 'ac_min', 'ai_min', 'ah_min', 'at_min', 'am_max', 'ac_max', 'ai_max', 'ah_max', 'at_max',\
			'te_min', 'ne_min', 'ti_min', 'th_min', 'te_max', 'ne_max', 'ti_max', 'th_max', 'zeff_max', 'zeff_min'])
		realLen.extend([(21,1)]*30)
		realList.extend(['bnfou'])
		realLen.extend([(25,41)])
		realList.extend(['dregcoil_rcws_rbound_c_opt','dregcoil_rcws_rbound_s_opt','dregcoil_rcws_zbound_c_opt','dregcoil_rcws_zbound_s_opt'])
		realLen.extend([(65,65)]*4)
		realList.extend(['te_opt','ti_opt','ne_opt','th_opt','zeff_opt'])
		realLen.extend([(21,1)]*5)
		realList.extend(['ne_aux_s', 'te_aux_s', 'ti_aux_s', 'th_aux_s', 'zeff_aux_s', \
			'phi_aux_s', 'beamj_aux_s', 'bootj_aux_s', 'sfincs_s', 'emis_xics_s',\
			'ne_aux_f', 'te_aux_f', 'ti_aux_f', 'th_aux_f', 'zeff_aux_f', \
			'phi_aux_f', 'beamj_aux_f', 'bootj_aux_f', 'emis_xics_f'])
		realLen.extend([(ndatafmax,1)]*19)
		realList.extend(['dam_s_opt', 'dam_f_opt', 'dac_s_opt', 'dac_f_opt', 'dai_s_opt', 'dai_f_opt', 'dphi_s_opt', 'dphi_f_opt', \
			'dne_f_opt', 'dte_f_opt', 'dti_f_opt', 'dth_f_opt', 'dzeff_f_opt', 'dbeamj_f_opt', 'dbootj_f_opt','dat_f_opt', \
			'dah_f_opt', 'demis_xics_f_opt'])
		realLen.extend([(ndatafmax,1)]*18)
		realList.extend(['am_f_min', 'ac_f_min', 'ai_f_min', 'phi_f_min', 'ne_f_min', \
			'te_f_min', 'ti_f_min', 'th_f_min', 'zeff_f_min', 'emis_xics_f_min', \
			'beamj_f_min', 'bootj_f_min', 'ah_f_min', 'at_f_min'])
		realLen.extend([(ndatafmax,1)]*14)
		realList.extend(['am_f_max', 'ac_f_max', 'ai_f_max', 'phi_f_max', 'ne_f_max', \
			'te_f_max', 'ti_f_max', 'th_f_max', 'zeff_f_max',  'emis_xics_f_max', \
			'beamj_f_max', 'bootj_f_max', 'ah_f_max', 'at_f_max'])
		realLen.extend([(ndatafmax,1)]*14)
		realList.extend(['daxis_opt','raxis_min','raxis_max','zaxis_min','zaxis_max'])
		realLen.extend([(ntord+1,1)]*5)
		realList.extend(['rhobc', 'modemn', 'dbound_opt', 'drho_opt', 'rbc_min', \
			'rbc_max', 'zbs_min', 'zbs_max', 'bound_min', 'bound_max', 'zbc_min',\
			'zbc_max', 'rbs_min', 'rbs_max'])
		realLen.extend([(2*ntord+1,1+mpol1d)]*14)
		realList.extend(['deltamn','ddeltamn_opt','delta_min','delta_max'])
		realLen.extend([(2*ntord+1,2*mpol1d+1)]*4)
		realList.extend(['coil_splinesx','coil_splinesy','coil_splinesz'])
		realLen.extend([(nigroup,44)]*3)
		realList.extend(['coil_splinefx','coil_splinefy','coil_splinefz','dcoil_spline', \
			'coil_splinefx_min','coil_splinefy_min','coil_splinefz_min','coil_splinefx_max','coil_splinefy_max','coil_splinefz_max'])
		realLen.extend([(nigroup,40)]*10)
		realList.extend(['regcoil_rcws_rbound_c', 'regcoil_rcws_rbound_s','regcoil_rcws_rbound_c_min',\
			'regcoil_rcws_rbound_s_min','regcoil_rcws_rbound_c_max', 'regcoil_rcws_rbound_s_max',\
			'regcoil_rcws_zbound_c', 'regcoil_rcws_zbound_s','regcoil_rcws_zbound_c_min', \
			'regcoil_rcws_zbound_s_min','regcoil_rcws_zbound_c_max', 'regcoil_rcws_zbound_s_max'])
		realLen.extend([(65,65)]*12)
		realList.extend(['drosenbrock_x_opt','rosenbrock_x','rosenbrock_x_min','rosenbrock_x_max'])
		realLen.extend([(20,1)]*4)
		charList=['sfincs_er_option', 'equil_type', 'te_type', 'ne_type', \
			'ti_type', 'th_type', 'beamj_type','bootj_type','zeff_type','emis_xics_type',\
			'fixedcoilname','regcoil_nescin_filename','bootcalc_type','phi_type','coil_type']
		charLen=[(256,1),(256,1),(256,1),(256,1),(256,1),(256,1),(256,1),\
			(256,1),(256,1),(256,1),(256,1),(256,1),(256,1),(256,1),(nigroup,1)]
		var_data = self.get_module_vars(module_name,booList,booLen,intList,intLen,realList,realLen,charList,charLen,ldefined_size_arrays=True)
		# Get target
		module_name = self.s1+'stellopt_targets_'+self.s2
		booList=['lneed_magdiag', 'lglobal_txport']
		booLen=[1]*len(booList)
		booList.extend(['lbooz','lmse_extcur'])
		booLen.extend([(nsd,1),(512,1)])
		intList=['mboz', 'nboz', 'numjstar', 'nz_txport', 'nalpha_txport', 'nruns_dkes',\
			'nu_orbit', 'nv_orbit', 'np_orbit', 'mlmnb_kink', 'ivac_kink', 'mmaxdf_kink', \
			'nmaxdf_kink', 'nra_ece', 'nphi_ece', 'numws', 'nu_bnorm', 'nv_bnorm', \
			'npts_biot', 'npts_clen', 'npts_torx', 'npts_curv', 'npts_csep', 'npts_cself', \
			'npts_crect', 'npts_cpoly']
		intLen=[1]*len(intList)
		intList.extend(['mlmns_kink', 'lssl_kink', 'lssd_kink','nj_kink','nk_kink'])
		intLen.extend([(16,1)]*5)
		realList=['target_phiedge', 'sigma_phiedge', 'target_curtor', 'sigma_curtor', 'target_curtor_max', 'sigma_curtor_max', \
			'target_volume', 'sigma_volume', 'target_beta', 'sigma_beta', 'target_betator', 'sigma_betator', \
			'target_betapol', 'sigma_betapol', 'target_wp', 'sigma_wp', 'target_aspect', 'sigma_aspect', \
			'target_rbtor', 'sigma_rbtor', 'target_r0', 'sigma_r0', 'target_z0', 'sigma_z0', \
			'target_b0', 'sigma_b0', 'target_aspect_max', 'sigma_aspect_max', 'width_aspect_max', 'target_gradp_max', \
			'sigma_gradp_max', 'width_gradp_max', 'target_pmin', 'sigma_pmin', 'width_pmin', 'target_curvature', \
			'sigma_curvature', 'target_kappa', 'sigma_kappa', 'phi_kappa', 'target_kappa_box', 'sigma_kappa_box', \
			'phi_kappa_box', 'target_kappa_avg', 'sigma_kappa_avg', 'target_x', 'sigma_x' ,'target_y', \
			'sigma_y', 'qm_ratio', 'cutoff_te_line', 'target_vessel', 'sigma_vessel', 'alpha_start_txport', \
			'alpha_end_txport', 'nu_dkes_erdiff', 'ep_dkes_erdiff', 'em_dkes_erdiff', 'mass_orbit', 'z_orbit', \
			'target_coil_bnorm', 'sigma_coil_bnorm', 'target_regcoil_winding_surface_separation', 'sigma_regcoil_winding_surface_separation',\
			'target_regcoil_current_density', 'sigma_regcoil_current_density', 'target_curvature_p2', 'sigma_curvature_p2', \
			'target_coilsep',  'sigma_coilsep', 'coilrectpfw']
		realLen=[1]*len(realList)
		realList.extend(['target_rosenbrock_f','sigma_rosenbrock_f'])
		realLen.extend([(20,1),(20,1)])
		realList.extend(['target_press', 'sigma_press', 'r_press', 'z_press', 'phi_press', \
			's_press', 'target_pressprime', 'sigma_pressprime', 'r_pressprime', 'z_pressprime', \
			'phi_pressprime', 's_pressprime', 'target_te', 'sigma_te', 'r_te', \
			'z_te', 'phi_te', 's_te', 'target_ne', 'sigma_ne', \
			'r_ne', 'z_ne', 'phi_ne', 's_ne', 'target_ne_line', \
			'sigma_ne_line', 'r0_ne_line', 'phi0_ne_line', 'z0_ne_line', 'r1_ne_line', \
			'phi1_ne_line', 'z1_ne_line', 'target_te_line', 'sigma_te_line', 'r0_te_line', \
			'phi0_te_line', 'z0_te_line', 'r1_te_line', 'phi1_te_line', 'z1_te_line', \
			'target_ti_line','sigma_ti_line', 'r0_ti_line', 'phi0_ti_line', 'z0_ti_line', \
			'r1_ti_line', 'phi1_ti_line', 'z1_ti_line', 'target_zeff_line', 'sigma_zeff_line', \
			'r0_zeff_line', 'phi0_zeff_line', 'z0_zeff_line', 'r1_zeff_line', 'phi1_zeff_line', \
			'z1_zeff_line', 'target_visbrem_line', 'sigma_visbrem_line', 'lambda_visbrem_line', 'r0_visbrem_line', \
			'phi0_visbrem_line', 'z0_visbrem_line', 'r1_visbrem_line', 'phi1_visbrem_line', 'z1_visbrem_line', \
			'calib_visbrem_line', 'target_xics', 'sigma_xics', 'target_xics_bright', 'sigma_xics_bright', \
			'target_xics_w3', 'sigma_xics_w3', 'target_xics_v', 'sigma_xics_v', 'r0_xics', \
			'phi0_xics', 'z0_xics', 'r1_xics', 'phi1_xics', 'z1_xics', \
			'target_faraday', 'sigma_faraday', 'r0_faraday', 'phi0_faraday', 'z0_faraday', \
			'r1_faraday', 'phi1_faraday', 'z1_faraday', 'target_sxr','sigma_sxr', \
			'r0_sxr', 'phi0_sxr', 'z0_sxr', 'r1_sxr', 'phi1_sxr', \
			'z1_sxr', 'target_ti', 'sigma_ti', 'r_ti', 'z_ti', \
			'phi_ti', 's_ti', 'target_vphi', 'sigma_vphi', 'r_vphi', \
			'z_vphi', 'phi_vphi', 's_vphi', 'target_iota', 'sigma_iota', \
			'r_iota', 'z_iota', 'phi_iota', 's_iota', 'target_vaciota', \
			'sigma_vaciota', 'r_vaciota', 'z_vaciota', 'phi_vaciota', 's_vaciota', \
			'target_mse', 'sigma_mse', 'r_mse', 'z_mse', 'phi_mse', \
			's_mse', 'a1_mse', 'a2_mse', 'a3_mse', 'a4_mse', \
			'a5_mse', 'a6_mse', 'a7_mse', 'vac_mse', 'target_segrog', \
			'sigma_segrog', 'target_fluxloop', 'sigma_fluxloop', 'target_extcur', 'sigma_extcur', \
			'balloon_theta', 'balloon_zeta', 'e_dkes', 'nu_dkes', 'nup_dkes_alpha', \
			'num_dkes_alpha', 'ep_dkes_alpha', 'em_dkes_alpha'])
		realLen.extend([(512,1)]*148)
		realList.extend(['target_bmin', 'sigma_bmin', 'target_bmax', 'sigma_bmax', 'target_jcurv', \
			'sigma_jcurv', 'target_jdotb', 'sigma_jdotb', 'target_balloon', 'sigma_balloon', \
			'target_bootstrap', 'sigma_bootstrap', 'target_neo', 'sigma_neo', 'target_jstar', \
			'sigma_jstar', 'target_magwell', 'sigma_magwell', 'target_helicity', 'sigma_helicity', \
			'target_helicity_old', 'sigma_helicity_old', 'target_resjac', 'sigma_resjac', 'xm_resjac', \
			'xn_resjac', 'target_txport', 'sigma_txport', 's_txport', 'target_dkes', \
			'sigma_dkes', 'target_dkes_erdiff', 'sigma_dkes_erdiff', 'target_dkes_alpha', 'sigma_dkes_alpha', \
			'target_gamma_c', 'sigma_gamma_c', 'target_orbit', 'sigma_orbit'])
		realLen.extend([(nsd,1)]*39)
		realList.extend(['target_bprobe', 'sigma_bprobe'])
		realLen.extend([(2048,1)]*2)
		realList.extend(['target_separatrix', 'sigma_separatrix', 'r_separatrix', 'z_separatrix', 'phi_separatrix', \
			'target_limiter', 'sigma_limiter', 'r_limiter', 'z_limiter', 'phi_limiter'])
		realLen.extend([(361,361)]*10)
		realList.extend(['vll_orbit', 'mu_orbit', 'vperp_orbit'])
		realLen.extend([(16384,1)]*3)
		realList.extend(['target_kink', 'sigma_kink'])
		realLen.extend([(16,1)]*2)
		realList.extend(['target_ece', 'sigma_ece', 'freq_ece'])
		realLen.extend([(16,512)]*3)
		realList.extend(['antennaposition_ece', 'targetposition_ece', 'rbeam_ece', 'rfocus_ece'])
		realLen.extend([(16,3)]*4)
		realList.extend(['target_regcoil_chi2_b', 'sigma_regcoil_chi2_b'])
		realLen.extend([(16900,1)]*4)
		realList.extend(['target_coillen', 'sigma_coillen', 'target_coilsegvar', 'sigma_coilsegvar', 'target_coilcrv',  \
			'sigma_coilcrv', 'target_coilself', 'sigma_coilself', 'target_coiltorvar', 'sigma_coiltorvar', \
			'thwt_coiltorvar', 'coilrectvmin', 'coilrectvmax', 'coilrectduu', 'coilrectdul', \
			'target_coilrect', 'sigma_coilrect', 'target_coilpoly', 'sigma_coilpoly'])
		realLen.extend([(nigroup,1)]*19)
		realList.extend(['kopolyu', 'kopolyv'])
		realLen.extend([(128,16)]*2)
		charList=['magdiag_coil', 'vessel_string', 'txport_proxy', 'vessel_ece', 'mirror_ece', 'targettype_ece', 'antennatype_ece']
		charLen=[(256,1)]*7
		target_data = self.get_module_vars(module_name,booList,booLen,intList,intLen,realList,realLen,charList,charLen,ldefined_size_arrays=True)
		var_data['lbound_opt'] = np.reshape(var_data['lbound_opt'],(mpol1d+1,2*ntord+1))
		var_data['dbound_opt'] = np.reshape(var_data['dbound_opt'],(mpol1d+1,2*ntord+1))
		var_data['bound_min'] = np.reshape(var_data['bound_min'],(mpol1d+1,2*ntord+1))
		var_data['bound_max'] = np.reshape(var_data['bound_max'],(mpol1d+1,2*ntord+1))
		var_data['lrho_opt'] = np.reshape(var_data['lrho_opt'],(mpol1d+1,2*ntord+1))
		var_data['drho_opt'] = np.reshape(var_data['drho_opt'],(mpol1d+1,2*ntord+1))
		var_data['ldeltamn_opt'] = np.reshape(var_data['ldeltamn_opt'],(2*mpol1d+1,2*ntord+1))
		var_data['ddeltamn_opt'] = np.reshape(var_data['ddeltamn_opt'],(2*mpol1d+1,2*ntord+1))
		var_data['delta_min'] = np.reshape(var_data['delta_min'],(2*mpol1d+1,2*ntord+1))
		var_data['delta_max'] = np.reshape(var_data['delta_max'],(2*mpol1d+1,2*ntord+1))
		return global_data, var_data, target_data

	def write_stellopt_input(self,filename,global_dict=None,var_dict=None,target_dict=None):
		"""Wrappers writing of the STELLOPT OPTIMUM namelist

		This routine wrappers write_stellopt_input in LIBSTELL

		Parameters
		----------
		filename : str
			Path to input file.
		out_dict : dict (optional)
			Dictionary of items to change.
		"""
		import ctypes as ct
		# Global module
		if global_dict:
			module_name = self.s1+'stellopt_globals_'+self.s2
			for key in global_dict:
				self.set_module_var(module_name,key,global_dict[key])
		# var module
		if var_dict:
			module_name = self.s1+'stellopt_vars_'+self.s2
			for key in var_dict:
				self.set_module_var(module_name,key,var_dict[key])
		# target module
		if target_dict:
			module_name = self.s1+'stellopt_targets_'+self.s2
			for key in target_dict:
				self.set_module_var(module_name,key,target_dict[key])
		# Now write namelist
		module_name = self.s1+'stellopt_input_mod_'+self.s2
		write_stellopt_input = getattr(self.libstell,module_name+'_write_optimum_namelist_byfile'+self.s3)
		write_stellopt_input.argtypes = [ct.c_char_p, ct.c_long]
		write_stellopt_input.restype=None
		write_stellopt_input(filename.encode('UTF-8'),len(filename))

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
			'jcurv','specw','jdotb','bdotb','dmerc','dshear','dwell','dcurr','dgeod','equif',\
			'bdotgradv']
		realLen = [(scalar_data['ns'],1)]*len(realList)
		realList.extend(['xm','xn','xm_nyq','xn_nyq'])
		realLen.extend([(scalar_data['mnmax'],1)]*2)
		realLen.extend([(scalar_data['mnmax_nyq'],1)]*2)
		realList.extend(['am','ac','ai'])
		realLen.extend([(21,1)]*3)
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

	def read_boozer(self,file):
		"""Reads a boozmn file and returns a dictionary

		This routine wrappers read_boozer in LIBSTELL and returns
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
		module_name = self.s1+'read_boozer_mod_'+self.s2
		read_boozer = getattr(self.libstell,module_name+'_read_boozer_file'+self.s3)
		read_boozer.argtypes=[ct.c_char_p, ct.POINTER(ct.c_int), ct.POINTER(ct.c_int), ct.c_long]
		read_boozer.restype=None
		ierr = ct.c_int(0)
		iopen = ct.c_int(0)
		read_boozer(file.encode('UTF-8'), ct.byref(ierr), ct.byref(iopen), len(file))
		if not (ierr.value == 0):
			return None
		# Setup Arrays
		out_data={}
		# Get Scalars
		booList  = ['lasym_b']
		booLen   = [1]*len(booList)
		intList  = ['mnboz_b', 'mboz_b', 'nboz_b', 'nfp_b', 'ns_b']
		intLen   = [1]*len(intList)
		realList  = ['aspect_b', 'rmax_b', 'rmin_b', 'betaxis_b']
		realLen   = [1]*len(realList)
		scalar_data = self.get_module_vars(module_name,booList,booLen,intList,intLen,realList,realLen)
		ns = scalar_data['ns_b']
		mnmax = scalar_data['mnboz_b']
		# Get 1D Real Arrays
		intList  = ['idx_b', 'ixm_b', 'ixn_b']
		intLen   = [(ns,1),(mnmax,1),(mnmax,1)]
		realList = ['iota_b','pres_b','phip_b','phi_b','beta_b','buco_b','bvco_b']
		realLen = [(ns,1)]*len(realList)
		# Add 2D Arrays
		realList.extend(['bmnc_b','rmnc_b','zmns_b','pmns_b', 'gmnc_b'])
		realLen.extend([(ns,mnmax)]*5)
		if scalar_data['lasym_b']:
			realList.extend(['bmns_b', 'rmns_b', 'zmnc_b', 'pmnc_b', 'gmns_b'])
			realLen.extend([(ns,mnmax)]*5)
		array_data = self.get_module_vars(module_name,intVar=intList,intLen=intLen,realVar=realList,realLen=realLen)
		# Return
		return scalar_data | array_data

	def read_nescout(self,file):
		"""Reads a nescout file and returns a dictionary

		This routine wrappers read_nescout in LIBSTELL and returns
		a dictionary of values

		Parameters
		----------
		file : str
			Path to nescout file.
		Returns
		----------
		vars : dict
			Dictionary of module variables
		"""
		import ctypes as ct
		module_name = self.s1+'read_nescoil_mod_'+self.s2
		read_wout = getattr(self.libstell,module_name+'_read_nescout'+self.s3)
		read_wout.argtypes=[ct.c_char_p, ct.POINTER(ct.c_int), ct.c_long]
		read_wout.restype=None
		ierr = ct.c_int(0)
		read_wout(file.encode('UTF-8'), ct.byref(ierr), len(file))
		if not (ierr.value == 0):
			return None
		# Setup Arrays
		out_data={}
		# Get Scalars
		booList  = ['lasym']
		booLen   = [1]*len(booList)
		intList  = ['nu', 'nv', 'nu1', 'nv1', 'mpol', 'ntor', 'mf', 'nf', 'md', 'nd', 'np', \
					'ibex', 'mstrt', 'mstep', 'mkeep', 'mdspw', 'w_psurf', 'w_csurf', \
					'w_bnuv', 'w_jsurf', 'w_xerr', 'w_svd', 'mnmax_plasma', \
					'mnmax_surface', 'nmax', 'mnd', 'nuv', 'nuv1', 'nuvh', 'nuvh1',\
					'mnmax_pot']
		intLen   = [1]*len(intList)
		realList = ['iota_edge', 'phip_edge', 'curpol', 'cut', 'cup', 'curwt', 'trgwt','alp']
		realLen = [1]*len(realList)
		scalar_data = self.get_module_vars(module_name,booList,booLen,intList,intLen,realList,realLen)
		# Get 1D Int Arrays
		intList= ['xm_plasma', 'xn_plasma','xm_surface', 'xn_surface', 'xm_pot', 'xn_pot']
		intLen = [(scalar_data['mnmax_plasma'],1),(scalar_data['mnmax_plasma'],1),\
			(scalar_data['mnmax_surface'],1),(scalar_data['mnmax_surface'],1),\
			(scalar_data['mnmax_pot'],1),(scalar_data['mnmax_pot'],1)]
		# Get 1D Real Arrays
		realList = ['rmnc_plasma', 'zmns_plasma', 'rmns_plasma', \
			'zmnc_plasma', 'lmnc_plasma', 'lmns_plasma' ]
		realLen = [(scalar_data['mnmax_plasma'],1)]*len(realList)
		realList.extend(['rmnc_surface', 'zmns_surface', 'rmns_surface', \
			'zmnc_surface'])
		realLen.extend([(scalar_data['mnmax_surface'],1)]*4)
		realList.extend(['potmns_surface'])
		realLen.extend([(scalar_data['mnmax_pot'],1)]*1)
		if scalar_data['w_psurf'] > 1:
			realList.extend(['x_plasma', 'y_plasma', 'z_plasma', 'r_plasma'])
			realLen.extend([(scalar_data['nuvh1'],1)]*4)
		if scalar_data['w_psurf'] > 2:
			realList.extend(['dsur_plasma', 'nx_plasma', 'ny_plasma', 'nz_plasma'])
			realLen.extend([(scalar_data['nuvh1'],1)]*4)
		if scalar_data['w_psurf'] > 3:
			realList.extend(['dxdu_plasma', 'dydu_plasma', 'dxdv_plasma', 'dydv_plasma'])
			realLen.extend([(scalar_data['nuvh1'],1)]*4)
		if scalar_data['w_bnuv'] > 0:
			realList.extend(['db_normal', 'babs'])
			realLen.extend([(scalar_data['nuvh1'],1)]*2)
		if scalar_data['w_bnuv'] > 1:
			realList.extend(['bn_plasma'])
			realLen.extend([(scalar_data['nuvh1'],1)])
		if scalar_data['w_csurf'] > 1:
			realList.extend(['x_surface', 'y_surface', 'z_surface', 'r_surface'])
			realLen.extend([(scalar_data['nuvh'],1)]*4)
		if scalar_data['w_csurf'] > 2:
			realList.extend(['dsur_surface','nx_surface', 'ny_surface', 'nz_surface'])
			realLen.extend([(scalar_data['nuvh'],1)]*4)
		if scalar_data['w_csurf'] > 3:
			realList.extend(['xcur_surface', 'ycur_surface', 'zcur_surface'])
			realLen.extend([(scalar_data['nuvh'],1)]*3)
		array_data = self.get_module_vars(module_name,intVar=intList,intLen=intLen,realVar=realList,realLen=realLen)
		# Return
		return scalar_data | array_data

	def nescout_bfield_init(self,nu=128,nv=128):
		"""Initializes the magnetic field from a NESCOIL surface current

		This routine evaluates the magnetic field at a point in space
		from a NESCOIL surface current.

		Parameters
		----------
		nu : int (optional)
			Number of poloidal gridpoints (default: 128)
		nv : int (optional)
			Number of toroidal gridpoints (default: 128)
		"""
		import ctypes as ct
		module_name = self.s1+'read_nescoil_mod_'+self.s2
		# initialize values
		bfield_init = getattr(self.libstell,module_name+'_nescoil_bfield_init_ctypes'+self.s3)
		bfield_init.argtypes=[ct.POINTER(ct.c_int),ct.POINTER(ct.c_int)]
		bfield_init.restype=None
		nu_ctype = ct.c_int(nu)
		nv_ctype = ct.c_int(nv)
		bfield_init(nu_ctype,nv_ctype)

	def nescout_bfield(self,x,y,z,istat=0):
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
		istat : integer
			Status optional (set to -327 to use explicit integration)
		Returns
		----------
		bx : float
			X-component of magnetic field [T]
		by : float
			Y-component of magnetic field [T]
		bZ : float
			Z-component of magnetic field [T]
		"""
		import ctypes as ct
		module_name = self.s1+'read_nescoil_mod_'+self.s2
		bfield = getattr(self.libstell,module_name+'_nescoil_bfield_adapt_dbl'+self.s3)
		bfield.argtypes=[ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), \
						 ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), \
						 ct.POINTER(ct.c_int)]
		bfield.restype=None
		x_ctype = ct.c_double(x)
		y_ctype = ct.c_double(y)
		z_ctype = ct.c_double(z)
		bx_ctype = ct.c_double(0.0)
		by_ctype = ct.c_double(0.0)
		bz_ctype = ct.c_double(0.0)
		istat_ctype = ct.c_int(istat)
		bfield(ct.byref(x_ctype),ct.byref(y_ctype),ct.byref(z_ctype),\
				 ct.byref(bx_ctype),ct.byref(by_ctype),ct.byref(bz_ctype),ct.byref(istat_ctype))
		return bx_ctype.value,by_ctype.value,bz_ctype.value

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
				#print(temp,booLen[i])
				if booLen[i]==1:
					out_data[temp]=ct.c_bool.in_dll(self.libstell,modName+'_'+temp+self.s3).value
				else:
					# This works because fortran has 4 byte sized booleans
					if ldefined_size_arrays : ftemp=ct.c_int*prod(booLen[i])
					out_data[temp]=npct.as_array(ftemp.in_dll(self.libstell,modName+'_'+temp+self.s3),booLen[i])>0
		# Integers
		if intVar:
			ftemp = ct.POINTER(ct.c_int)
			for i,temp in enumerate(intVar):
				#print(temp,intLen[i])
				if intLen[i]==1:
					out_data[temp]=ct.c_int.in_dll(self.libstell,modName+'_'+temp+self.s3).value
				else:
					if ldefined_size_arrays : ftemp=ct.c_int*prod(intLen[i])
					out_data[temp]=npct.as_array(ftemp.in_dll(self.libstell,modName+'_'+temp+self.s3),intLen[i])
		# Reals
		if realVar:
			ftemp = ct.POINTER(ct.c_double)
			for i,temp in enumerate(realVar):
				#print(temp,realLen[i])
				if realLen[i]==1:
					out_data[temp]=ct.c_double.in_dll(self.libstell,modName+'_'+temp+self.s3).value
				else:
					if ldefined_size_arrays : ftemp=ct.c_double*prod(realLen[i])
					out_data[temp]=npct.as_array(ftemp.in_dll(self.libstell,modName+'_'+temp+self.s3),realLen[i])
		# Characters
		if charVar:
			ftemp = ct.POINTER(ct.c_char)
			for i,temp in enumerate(charVar):
				#print(temp,charLen[i])
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
		elif type(val) == list:
			if type(val[0]) == bool:
				tt = ct.c_bool
			else:
				print(f'   Unrecognized list type ({var}):',type(val[0]))
				return
			n = val.ndim
			f = tt*val.size
		elif type(val) == np.ndarray:
			if val.dtype.type == np.bool_:
				tt = ct.c_bool
			elif val.dtype.type == np.int32:
				tt = ct.c_int
			elif val.dtype.type == np.int64:
				tt = ct.c_int
			elif val.dtype.type == np.float64:
				tt = ct.c_double
			else:
				print(f'   Unrecognized ndarray type ({var}):',val.dtype.type)
				return
			n = val.ndim
			f = tt*val.size
		elif type(val) == type(self):
			return # Do nothing
		else:
			print(f'   Unrecognized type ({var}):',type(val))
			return
		temp=f.in_dll(self.libstell,modName+'_'+var+self.s3)
		if type(val) == np.ndarray:
			if n==1:
				for i,col in enumerate(val):
					temp[i] = val[i]
		elif type(val) == str:
			temp.value = val.encode('UTF-8')
		else:
			temp.value = val
		return

	def vmec_get_flxcoord(self,s,u,v):
		"""Wrapper to the get_flxcoord function

		This routine wrappers the get_flxcoord function found in
		vmec_utils.  It takes s, u, and v as inputs and returns
		the R, phi, Z, dRds, dZds, dRdu, and dZdu values at that
		point.

		Parameters
		----------
		s : real
			Normalized toroidal flux (VMEC).
		u : real
			Poloidal angle [rad] (VMEC)
		v : real
			Toroidal angle [rad] [VMEC]

		Returns
		-------
		R : real
			Cylindical R coordinate [m].
		phi : real
			Cylindical phi coordinate [rad].
		Z : real
			Cylindical Z coordinate [m].
		dRds : real
			Derivative of R coordiante with respect to s (dR/ds)
		dZds : real
			Derivative of Z coordiante with respect to s (dZ/ds)
		dRdu : real
			Derivative of R coordiante with respect to u (dR/du)
		dZdu : real
			Derivative of Z coordiante with respect to u (dZ/du)
		"""
		import ctypes as ct
		module_name = self.s1+'vmec_utils_'+self.s2
		get_flxcoord = getattr(self.libstell,module_name+'_get_flxcoord_python'+self.s3)
		get_flxcoord.argtypes = [ct.POINTER(ct.c_double),ct.POINTER(ct.c_double),
			ct.POINTER(ct.c_double),ct.POINTER(ct.c_double),ct.POINTER(ct.c_double),ct.POINTER(ct.c_double), \
			ct.c_long,ct.c_long]
		get_flxcoord.restype=None
		x1 = (ct.c_double*3)(0,0,0)
		c_flx = (ct.c_double*3)(s,u,v)
		rs = ct.c_double(0)
		zs = ct.c_double(0)
		ru = ct.c_double(0)
		zu = ct.c_double(0)
		get_flxcoord(x1,c_flx,\
			ct.byref(rs),ct.byref(zs),ct.byref(ru),ct.byref(zu),len(x1),len(c_flx))
		R = x1[0]
		v = x1[1]
		Z = x1[2]
		rs1 = rs.value
		zs1 = zs.value
		ru1 = ru.value
		zu1 = zu.value
		return R,v,Z,rs1,zs1,ru1,zu1

	def vmec_getBcyl_wout(self,R,phi,Z):
		"""Wrapper to the GetBcyl_WOUT function

		This routine wrappers the GetBcyl_WOUT function found in
		vmec_utils.  It takes R, phi, and Z as inputs and returns
		the Br, Bphi, Bz, s, and u values at that point. A status flag
		is also returned (info) which indicates
			 0: successfully find s,u point
			-1: did not converge
			-3: sflux > 1, probably

		Parameters
		----------
		R : real
			Cylindical R coordinate [m].
		phi : real
			Cylindical phi coordinate [rad].
		Z : real
			Cylindical Z coordinate [m].
		Returns
		-------
		br : real
			Magnetic field in cylindrical R direction [T].
		bphi : real
			Magnetic field in cylindrical phi direction [T].
		bz : real
			Magnetic field in cylindrical Z direction [T].
		s : real
			Normalized toroidal flux coordinate [arb].
		u : real
			Poloidal angle coordinate (VMEC angle) [rad].
		info: int
			Status of inverse lookup.
		"""
		import ctypes as ct
		module_name = self.s1+'vmec_utils_'+self.s2
		getBcyl = getattr(self.libstell,module_name+'_getbcyl_wout'+self.s3)
		getBcyl.argtypes = [ct.POINTER(ct.c_double),ct.POINTER(ct.c_double),ct.POINTER(ct.c_double), \
			ct.POINTER(ct.c_double),ct.POINTER(ct.c_double),ct.POINTER(ct.c_double), \
			ct.POINTER(ct.c_double),ct.POINTER(ct.c_double),ct.POINTER(ct.c_int)]
		getBcyl.restype=None
		r_temp = ct.c_double(R)
		phi_temp = ct.c_double(phi)
		z_temp = ct.c_double(Z)
		br_temp = ct.c_double(0)
		bphi_temp = ct.c_double(0)
		bz_temp = ct.c_double(0)
		s_temp = ct.c_double(0)
		u_temp = ct.c_double(0)
		info_temp = ct.c_int(0)
		getBcyl(ct.byref(r_temp),ct.byref(phi_temp),ct.byref(z_temp), \
			ct.byref(br_temp),ct.byref(bphi_temp),ct.byref(bz_temp), \
			ct.byref(s_temp),ct.byref(u_temp),ct.byref(info_temp))
		Br = br_temp.value
		Bphi = bphi_temp.value
		Bz = bz_temp.value
		s = s_temp.value
		u = u_temp.value
		info = info_temp.value
		return Br,Bphi,Bz,s,u,info

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
		test = None

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
			f[k,:,:]=np.matmul((fmn*sinmt).T, cosnz)+np.matmul((fmn*cosmt).T, sinnz)
		return f

	def isotoro(self,r,z,phi,svals,*args,**kwargs):
		"""Plot a surface in 3D using VTK

		This routine plots a 3D surface using the VTK library when
		passed a r [m], z [m], and phi [rad] arrays as produced by
		the sfunct and cfunct functions. The user may supply a surface
		to plot as an index or array of indices. Pass svals=-1 if the
		arrays only have one radial gridpoint to plot.

		Parameters
		----------
		r : ndarray
			Ordered list of R verticies [m] (ns,nu,nv)
		z : ndarray
			Ordered list of Z verticies [m] (ns,nu,nv)
		phi : ndarray
			Ordered list of phi coordiantes [rad] (nv)
		surface : int
			Surface to generate in ns
		vals : ndarray (optional)
			Ordered list of vertex values for coloring (ns,nu,nv)
		plot3D : plot3D object (optional)
			Plotting object to render to.
		lcloseu : boolean (optional)
			Close the grid in the poloidal direction. (default: True)
		lclosev : boolean (optional)
			Close the grid in the toroidal direction. (default: True)
		lcolorbar : boolean (optional)
			Turn the colorbar on (default: False)
		color : string (optional)
			Surface color name, overriden by vals (default: 'red')
		"""
		import numpy as np
		from libstell.plot3D import PLOT3D 
		# Handle input arguments
		plt  = kwargs.get('plot3D',None)
		lcu  = kwargs.get('lcloseu',True)
		lcv  = kwargs.get('lclosev',True)
		vals = kwargs.get('vals',None)
		lbar = kwargs.get('lcolorbar',False)
		color = kwargs.get('color','red')
		lrender = False
		if not plt:
			plt = PLOT3D()
			lrender = True
		# Figure out number of surfaces to plot
		if type(svals) is list:
			s = svals
		else:
			# Aviod plotting axis
			if svals == 0: svals = 1
			# Flag for plotting single surface array
			if r.shape[0] == 1: svals = 0
			s= [svals]
		nr = np.size(s)
		# Setup x,y,z helpers
		nu = np.size(r,1)
		nv = np.size(r,2)
		x_s = np.zeros((nu,nv))
		y_s = np.zeros((nu,nv))
		z_s = np.zeros((nu,nv))
		# Loop over radial values
		for k in range(nr):
			# Generate Coords
			x_s = r[s[k],:,:]*np.cos(np.broadcast_to(phi,(nv,nu))).T
			y_s = r[s[k],:,:]*np.sin(np.broadcast_to(phi,(nv,nu))).T
			z_s = z[s[k],:,:]
			[points,triangles] = plt.torusvertexTo3Dmesh(x_s,y_s,z_s,lcloseu=lcu,lclosev=lcv)
			# Handle Color
			if type(vals) != type(None): 
				scalar = plt.valuesToScalar(vals[s[k],:,:])
				# Add to Render
				plt.add3Dmesh(points,triangles,scalars=scalar,opacity=1.0/nr)
			else:
				plt.add3Dmesh(points,triangles,color=color,opacity=1.0/nr)
		# In case it isn't set by user.
		plt.setBGcolor()
		# Render if requested
		if lrender: plt.render()

	def generateSurface(self,r,z,phi,surface=None):
		"""Generates the vertex and indices arrays to render a surface

		This routine generates the verticies and faces lists which
		are used by many routines to render a surface in 3D.

		Parameters
		----------
		r : ndarray
			Ordered list of R verticies [m] (ns,nu,nv)
		z : ndarray
			Ordered list of Z verticies [m] (ns,nu,nv)
		phi : ndarray
			Ordered list of phi coordiantes [rad] (nv)
		surface : int (optional)
			Surface to generate in ns (default: outermost)

		Returns
		----------
		vertices : ndarray
			Vertex values [m]
		faces: ndarray
			Indices defining triangles
		"""
		import numpy as np
		[vertices_list,faces_list] = self.blenderSurface(r,z,phi,surface=surface)
		ndex = len(vertices_list)
		vertex = np.zeros((ndex,3))
		for i in range(ndex):
			vertex[i,0] = vertices_list[i][0]
			vertex[i,1] = vertices_list[i][1]
			vertex[i,2] = vertices_list[i][2]
		# For indices we need to get it out of tupple format
		ndex = len(faces_list)
		faces = np.zeros((ndex,3),dtype=int)
		for i in range(ndex):
			faces[i,0] = faces_list[i][0]
			faces[i,1] = faces_list[i][1]
			faces[i,2] = faces_list[i][2]
		return vertex,faces

	def surfaceSTL(self,r,z,phi,surface=None,filename='surface.stl'):
		"""Outputs an STL file from a surface

		This routine outputs an STL file for a given surface.

		Parameters
		----------
		r : ndarray
			Ordered list of R verticies [m] (ns,nu,nv)
		z : ndarray
			Ordered list of Z verticies [m] (ns,nu,nv)
		phi : ndarray
			Ordered list of phi coordiantes [rad] (nv)
		surface : int (optional)
			Surface to generate in ns (default: outermost)
		filename : str (optional)
			Filename for output file
		"""
		import numpy as np
		from stl import mesh
		vertex,faces = self.generateSurface(r,z,phi,surface)
		nfaces = faces.shape[0]
		wall_mesh = mesh.Mesh(np.zeros(nfaces, dtype=mesh.Mesh.dtype))
		for i, f in enumerate(faces):
			for j in range(3):
				wall_mesh.vectors[i][j] = vertex[f[j],:]
		wall_mesh.save(filename)

	def blenderSurface(self,r,z,phi,surface=None):
		"""Generates the lists Blender needs to render a flux surface

		This routine generates the verticies and faces lists which
		Blender needs to render a surface. It assumes
		that datapoints are not repeated in the theta  or
		phi direction, but those coordinates are periodic.

		Parameters
		----------
		r : ndarray
			Ordered list of R verticies [m] (ns,nu,nv)
		z : ndarray
			Ordered list of Z verticies [m] (ns,nu,nv)
		phi : ndarray
			Ordered list of phi coordiantes [rad] (nv)
		surface : int (optional)
			Surface to generate in ns (default: outermost)

		Returns
		----------
		vertices : list
			List of tuples defining verticies
		faces: list
			List of tubles defining faces
		"""
		import numpy as np
		# Generate volumetric coil rendering
		vertices = []
		faces = []
		if surface:
			k = surface
		else:
			k = r.shape[0]-1
		nu = np.size(r,1)
		nv = np.size(r,2)
		# Handle only a partial period
		if phi[nv-1]%(2.0*np.pi) == phi[0]:
			nv = nv - 1
			lseem = True
		else:
			lseem = False
		# Loop and construct
		for v in range(nv):
			for u in range(nu):
				x = r[k,u,v] * np.cos(phi[v])
				y = r[k,u,v] * np.sin(phi[v])
				vertices.append((x,y,z[k,u,v]))
				if (v == nv -1) and (not lseem): continue
				# Now do faces
				i1 = u + v * nu
				# Catch special case #1
				if u == nu-1:
					i2 = i1 + nu
					i3 = i1 + 1
					i4 = i1 - nu + 1
					if (v == nv - 1):
						i2 = u
						i3 = 0
				elif u < nu-1:
					i2 = i1 + nu
					i3 = i1 + nu +1
					i4 = i1 + 1
					if (v == nv -1):
						i2 = u
						i3 = u + 1
				faces.append((i1,i2,i4))
				faces.append((i2,i3,i4))
		return vertices,faces


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