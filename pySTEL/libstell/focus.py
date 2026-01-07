##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for reading and handling FOCUS
data.
"""

# Libraries
#from libstell.libstell import LIBSTELL, FourierRep

# Constants

# VMEC Class
class FOCUS():
	"""Class for working with VMEC equilibria

	"""
	def __init__(self):
		return

	def read_focus_HDF5(self,filename):
		"""Reads a FOCUS HDF5 file

		This routine reads and initilizes the FIELDLINES
		class with variable information from an HDF5 file.

		Parameters
		----------
		file : str
			Path to HDF5 file.
		"""
		import h5py
		import numpy as np
		focus_data = {}
		if '.h5' not in filename:
			filename = filename + '.h5'
		if 'focus_' not in filename:
			filename = 'focus_' + filename
		with h5py.File(filename,'r') as f:
			# Integers
			for temp in ['IsQuiet','IsSymmetric','case_surface',\
				'knotsurf','Nteta','Nzeta','case_init','case_coils',\
				'Ncoils','IsVaryCurrent','IsVaryGeometry','NFcoil',\
				'Nseg','case_optimize','IsNormalize','ISNormWeight',\
				'case_bnormal','case_length','case_curv','case_straight','curv_alpha',\
				'DF_maxiter','CG_maxiter','HN_maxiter','TN_maxiter',\
				'TN_reorder','case_postproc','save_freq','save_coils',\
				'save_harmonics','save_filaments','update_plasma',\
				'pp_ns','pp_maxiter','Nfp','iout','ibnorm','mbnorm',\
				'ibharm','mbharm','itflux','mtflux','ittlen','mttlen',\
				'icssep','mcssep','icurv','mcurv','axis_npoints',\
				'istr','mstr','iccsep','mccsep','itors','mtors','inissin','mnissin']:
				if temp in f:
					setattr(self, temp, np.int64(f[temp][0]))
			# Reals
			for temp in ['knotsurf','init_current','init_radius',\
				'exit_xtol','weight_bnorm','weight_sbnorm','weight_bharm','weight_tflux',\
				'target_tflux','weight_isum','target_isum','weight_ttlen','target_length','curv_k0',\
				'weight_specw','weight_cssep','weight_gnorm','weight_inorm',\
				'weight_mnorm','weight_curv','weight_straight','weight_ccsep','weight_tors',\
				'ccsep_alpha','ccsep_beta','ccsep_skip','tors_alpha','case_tors','tors0',\
				'nissin_alpha','nissin_beta','penfun_nissin','nissin0','nissin_sigma','nissin_gamm',\
				'DF_tausta','DF_tauend','DF_xtol','DF_maxiter',\
				'CG_maxiter','CG_xtol','CG_wolfe_c1','CG_wolfe_c2',\
				'HN_maxiter','HN_xtol','HN_factor',\
				'TN_maxiter','TN_xtol','TN_cr',\
				'pp_phi','pp_raxis','pp_zaxis','pp_rmax',\
				'pp_zmax','pp_xtol','surf_vol','surf_vol_lim','Inorm','Gnorm','Mnorm',\
				'overlap','time_initialize','time_optimize','time_postproc']:
				if temp in f:
					setattr(self, temp, float(f[temp][0]))
			# Arrays
			for temp in ['xsurf','ysurf','zsurf','nx','ny','nz','nn',\
				'xsurf_lim','ysurf_lim','zsurf_lim','nx_lim','ny_lim','nz_lim','nn_lim'\
				'plas_Bn','Bn','Bx','By','Bz',\
				'limiter_Bn','Bn_limiter','Bx_limiter','By_limiter','Bz_limiter',\
				'evolution','coilspace',\
				'deriv','Bmnin','Bmnim','initial_Bmnc','initial_Bmns',\
				'target_Bmnc','target_Bmns','Bmnc','Bmns','coil_importance',\
				'LM_fvec','LM_fjac','ppr','ppz','iota','axis_phi','axis_r','axis_z',\
				'XYZB','booz_mnc','booz_mns','bmim','bmin']:
				if temp in f:
					setattr(self, temp, np.array(f[temp][:]))

	def read_focus_focus(self,filename):
		"""Reads a .focus file

		This routine reads the .focus files containg the various representations
		of coils and magnetic field sources.

		Parameters
		----------
		file : str
			Path to .focus file.
		"""
		import numpy as np
		f = open(filename,'r')
		lines = f.readlines()
		f.close()
		cline = 1
		ncoils = int(lines[cline])
		self.ncoils = ncoils
		self.coil_type = np.zeros((ncoils),dtype=int)
		self.symmetry_type = np.zeros((ncoils),dtype=int)
		self.name = []
		self.nseg = np.zeros((ncoils),dtype=int)
		self.current = np.zeros((ncoils))
		self.ifree = np.zeros((ncoils),dtype=int)
		self.length = np.zeros((ncoils))
		self.lfree = np.zeros((ncoils),dtype=int)
		self.lz = np.zeros((ncoils),dtype=int)
		self.target_length = np.zeros((ncoils))
		self.ncoef = np.zeros((ncoils),dtype=int)
		self.bz = np.zeros((ncoils))
		self.coef1x = []
		self.coef2x = []
		self.coef1y = []
		self.coef2y = []
		self.coef1z = []
		self.coef2z = []
		self.knots = []
		for i in range(ncoils):
			cline = cline + 3
			(ctype_txt,symm_txt,name_txt) = lines[cline].split()
			self.coil_type[i] = int(ctype_txt)
			self.symmetry_type[i] = int(symm_txt)
			self.name.append(name_txt)
			if self.coil_type[i] == 1: # Fourier
				cline = cline + 2
				(nseg_txt,current_txt,ifree_txt,length_txt,lfree_txt,target_length_txt) = lines[cline].split()
				self.nseg[i] = int(nseg_txt)
				self.ifree[i] = int(ifree_txt)
				self.lfree[i] = int(lfree_txt)
				self.current[i] = float(current_txt)
				self.length[i] = float(length_txt)
				self.target_length[i] = float(target_length_txt)
				cline = cline + 2
				self.ncoef[i] = int(lines[cline])
				cline = cline + 2
				self.coef1x.append([float(s) for s in lines[cline].split()])
				cline = cline + 1
				self.coef2x.append([float(s) for s in lines[cline].split()])
				cline = cline + 1
				self.coef1y.append([float(s) for s in lines[cline].split()])
				cline = cline + 1
				self.coef2y.append([float(s) for s in lines[cline].split()])
				cline = cline + 1
				self.coef1z.append([float(s) for s in lines[cline].split()])
				cline = cline + 1
				self.coef2z.append([float(s) for s in lines[cline].split()])
			elif self.coil_type[i] == 2: # Magnet
				cline = cline + 2
				(lc_txt,ox_txt,oy_txt,oz_txt,ic_txt,I_txt,mt_txt,mp_txt) = lines[cline].split()
				print('WARNING: Permanent Magnets not implemented')
			elif self.coil_type[i] == 3: # Background Btor, Bz
				cline = cline + 2
				(ifree_txt,current_txt,lz_txt,bz_txt) = lines[cline].split()
				self.ifree[i] = int(ifree_txt)
				self.current[i] = float(current_txt)
				self.lz[i] = int(lz_txt)
				self.bz[i] = float(bz_txt)
			elif self.coil_type[i] == 5: # Spline
				cline = cline + 2
				(nseg_txt,current_txt,ifree_txt,length_txt,lfree_txt,target_length_txt) = lines[cline].split()
				self.nseg[i] = int(nseg_txt)
				self.ifree[i] = int(ifree_txt)
				self.lfree[i] = int(lfree_txt)
				self.current[i] = float(current_txt)
				self.length[i] = float(length_txt)
				self.target_length[i] = float(target_length_txt)
				cline = cline + 2
				self.ncoef[i] = int(lines[cline])
				cline = cline + 2
				self.knots.append([float(s) for s in lines[cline].split()])
				cline = cline + 2
				self.coef1x.append([float(s) for s in lines[cline].split()])
				cline = cline + 1
				self.coef1y.append([float(s) for s in lines[cline].split()])
				cline = cline + 1
				self.coef1z.append([float(s) for s in lines[cline].split()])
			else:
				print(rf'Unkown coil type {self.coil_type[i]}, {self.name[i]}')
		return

	def write_focus_focus(self,filename='new.focus'):
		"""Writes a .focus source deffinition file

		This routine writes the .focus magnetic field source deffinition
		file.

		Parameters
		----------
		filename : string (optional)
			Source deffinition file name (default: new.focus)
		"""
		f=open(filename,'w')
		f.write(f" # Total number of coils\n")
		f.write(f"        {self.ncoils:6d}\n")
		for i in range(self.ncoils):
			f.write(f" #----------------- {i+1:d} ---------------------------\n")
			f.write(f" #coil_type   coil_symm  coil_name\n")
			f.write(f"   {self.coil_type[i]:3d}    {self.symmetry_type[i]:3d}    {self.name[i]}\n")
			if self.coil_type[i] == 1:
				f.write(f" #Nseg        current         Ifree         Length         Lfree  target_length\n")
				f.write(f"  {self.nseg[i]:4d}{self.current[i]:23.15E}   {self.ifree[i]:3d}{self.length[i]:23.15E}   {self.lfree[i]:3d}{self.target_length[i]:23.15E}\n")
				f.write(f' #NFcoil\n')
				f.write(f'  {self.ncoef[i]:3d}\n')
				f.write(f' #Fourier harmonics for coils ( xc; xs; yc; ys; zc; zs)\n')
				for j in range(self.ncoef[i]): f.write(f"{self.coef1x[i][j]:23.15E}")
				for j in range(self.ncoef[i]): f.write(f"{self.coef2x[i][j]:23.15E}")
				f.write('\n')
				for j in range(self.ncoef[i]): f.write(f"{self.coef1y[i][j]:23.15E}")
				for j in range(self.ncoef[i]): f.write(f"{self.coef2y[i][j]:23.15E}")
				f.write('\n')
				for j in range(self.ncoef[i]): f.write(f"{self.coef1z[i][j]:23.15E}")
				for j in range(self.ncoef[i]): f.write(f"{self.coef2z[i][j]:23.15E}")
				f.write('\n')
			elif self.coil_type[i] ==2:
				f.write(f' #  Lc  ox   oy   oz  Ic  I  mt  mp (note yet implemented in python interface)\n')
				f.write(f'   1   0.0  0.0  0.0  1 1.0E6  0.0  0.0\n')
			elif self.coil_type[i] == 3:
				f.write(f' # Ic     I    Lc  Bz  (Ic control I; Lc control Bz)\n')
				f.write(f'  {self.ifree[i]}  {self.current[i]:21.15E}  {self.lz[i]}  {self.bz[i]:21.15E}\n')
			elif self.coil_type[i] == 5:
				f.write(f' #NS          current         Ifree         Length         Lfree  target_length\n')
				f.write(f"  {self.nseg[i]:4d}{self.current[i]:23.15E}   {self.ifree[i]:3d}{self.length[i]:23.15E}   {self.lfree[i]:3d}{self.target_length[i]:23.15E}\n")
				f.write(f' #NCP\n')
				f.write(f'  {self.ncoef[i]:3d}\n')
				f.write(f' #Vector of knots\n')
				for j in range(self.ncoef[i]+4): f.write(f"{self.knots[i][j]:23.15E}")
				f.write(f'\n #Control Points Coordinates for coils ( x; y; z)  \n')
				for j in range(self.ncoef[i]): f.write(f"{self.coef1x[i][j]:23.15E}")
				f.write('\n')
				for j in range(self.ncoef[i]): f.write(f"{self.coef1y[i][j]:23.15E}")
				f.write('\n')
				for j in range(self.ncoef[i]): f.write(f"{self.coef1z[i][j]:23.15E}")
				f.write('\n')
			else:
				print(rf'Unkown coil type {self.ctype[i]}, {self.name[i]}')
		f.close()
		return

	def write_focus_plasma(self,nfp,xm,xn,rmnc,zmns,rmns=None,zmnc=None,xm_b=None,\
		xn_b=None,bmnc=None,bmns=None,filename='plasma.boundary'):
		"""Writes a focus boundary file

		This routine writes the FOCUS plasma boundary file.

		Parameters
		----------
		nfp : int
			Field periodicity
		xm : ndarray
			Poloidal mode array
		xn : ndarray
			Toroidal mode array
		rmnc : ndarray
			R cosine boundary harmonics
		zmns : ndarray
			Z sine boundary harmonics
		rmns : ndarray (optional)
			R sine boundary harmonics
		zmnc : ndarray (optional)
			Z cosine boundary harmonics
		xm_b : ndarray (optional)
			Poloidal mode array (B-normal)
		xn_b : ndarray (optional)
			Toroidal mode array (B-normal)
		bmnc : ndarray (optional)
			B-normal cosine boundary harmonics
		bmns : ndarray (optional)
			B-normal sine boundary harmonics
		filename : string (optional)
			Boundary file name (default: plasma.boundary)
		"""
		import numpy as np
		mnmax = len(xm)
		mnmax_b = 1
		if (type(xm_b) is not type(None)) and \
		   (type(xn_b) is not type(None)) and \
		   ((type(bmns) is not type(None)) or \
		   	(type(bmnc) is not type(None))):
			mnmax_b = len(xm_b)
		else:
			xm_b = [0]
			xn_b = [0]
			bmnc = [0]
			bmns = [0]
		if not (rmns and zmnc):
			rmns = np.zeros((mnmax))
			zmnc = np.zeros((mnmax))
		f=open(filename,'w')
		f.write(f'#Nfou Nfp NBnf\n')
		f.write(f'{int(mnmax)} {int(nfp)} {int(mnmax_b)}\n')
		f.write(f'#plasma boundary\n')
		f.write(f'# n m Rbc Rbs Zbc Zbs\n')
		for mn in range(mnmax):
			f.write(f'{int(xn[mn])} {int(xm[mn])} {rmnc[mn]:10.9e} {rmns[mn]:10.9e} {zmnc[mn]:10.9e} {zmns[mn]:10.9e}\n')
		f.write(f'#Bn harmonics\n')
		f.write(f'# n m bnc bns\n')
		for mn in range(mnmax_b):
			f.write(f'{int(xn_b[mn]):d} {int(xm_b[mn]):d} {bmnc[mn]:10.9e} {bmns[mn]:10.9e}\n')
		f.close()

	def write_focus_plasma_booz(self,nfp,xm,xn,rmnc,zmns,pmns,rmns=None,zmnc=None,pmnc=None,xm_b=None,\
		xn_b=None,bmnc=None,bmns=None,filename='plasma.boundary'):
		"""Writes a focus boundary file

		This routine writes the FOCUS plasma boundary file.

		Parameters
		----------
		nfp : int
			Field periodicity
		xm : ndarray
			Poloidal mode array
		xn : ndarray
			Toroidal mode array
		rmnc : ndarray
			R cosine boundary harmonics
		zmns : ndarray
			Z sine boundary harmonics
		pmns : ndarray
			P sine toridal angle harmonics
		rmns : ndarray (optional)
			R sine boundary harmonics
		zmnc : ndarray (optional)
			Z cosine boundary harmonics
		pmnc : ndarray
			P cosine toridal angle harmonics
		xm_b : ndarray (optional)
			Poloidal mode array (B-normal)
		xn_b : ndarray (optional)
			Toroidal mode array (B-normal)
		bmnc : ndarray (optional)
			B-normal cosine boundary harmonics
		bmns : ndarray (optional)
			B-normal sine boundary harmonics
		filename : string (optional)
			Boundary file name (default: plasma.boundary)
		"""
		import numpy as np
		mnmax = len(xm)
		mnmax_b = 1
		if (type(xm_b) is not type(None)) and \
		   (type(xn_b) is not type(None)) and \
		   ((type(bmns) is not type(None)) or \
		   	(type(bmnc) is not type(None))):
			mnmax_b = len(xm_b)
		else:
			xm_b = [0]
			xn_b = [0]
			bmnc = [0]
			bmns = [0]
		if not (rmns and zmnc):
			rmns = np.zeros((mnmax))
			zmnc = np.zeros((mnmax))
			pmnc = np.zeros((mnmax))
		f=open(filename,'w')
		f.write(f'#Nfou Nfp NBnf\n')
		f.write(f'{int(mnmax)} {int(nfp)} {int(mnmax_b)}\n')
		f.write(f'#plasma boundary\n')
		f.write(f'# n m Rbc Rbs Zbc Zbs\n')
		for mn in range(mnmax):
			f.write(f'{int(xn[mn]):5d} {int(xm[mn]):5d} {rmnc[mn]:10.9e} {rmns[mn]:10.9e} {zmnc[mn]:10.9e} {zmns[mn]:10.9e} {pmnc[mn]:10.9e} {pmns[mn]:10.9e}\n')
		f.write(f'#Bn harmonics\n')
		f.write(f'# n m bnc bns\n')
		for mn in range(mnmax_b):
			f.write(f'{int(xn_b[mn]):5d} {int(xm_b[mn]):5d} {bmnc[mn]:10.9e} {bmns[mn]:10.9e}\n')
		f.close()

	def write_focus_harmonics(self,tol,nfp,xm,xn,rmnc,zmns,pmns=None,rmns=None,zmnc=None,pmnc=None,xm_b=None,\
		xn_b=None,bmnc=None,bmns=None,filename='target.harmonics'):
		"""Writes a focus boundary file

		This routine writes the FOCUS plasma boundary file.

		Parameters
		----------
		tol : float
			Tolerance for boundary Harmonics targeting
		nfp : int
			Field periodicity
		xm : ndarray
			Poloidal mode array
		xn : ndarray
			Toroidal mode array
		rmnc : ndarray
			R cosine boundary harmonics
		zmns : ndarray
			Z sine boundary harmonics
		pmns : ndarray
			P sine toridal angle harmonics
		rmns : ndarray (optional)
			R sine boundary harmonics
		zmnc : ndarray (optional)
			Z cosine boundary harmonics
		pmnc : ndarray
			P cosine toridal angle harmonics
		xm_b : ndarray (optional)
			Poloidal mode array (B-normal)
		xn_b : ndarray (optional)
			Toroidal mode array (B-normal)
		bmnc : ndarray (optional)
			B-normal cosine boundary harmonics
		bmns : ndarray (optional)
			B-normal sine boundary harmonics
		filename : string (optional)
			Boundary file name (default: target.harmonics)
		"""
		import numpy as np
		mnmax = len(xm)
		mnmax_b = 1
		if (type(xm_b) is not type(None)) and \
		   (type(xn_b) is not type(None)) and \
		   ((type(bmns) is not type(None)) or \
		   	(type(bmnc) is not type(None))):
			mnmax_b = len(xm_b)
		else:
			xm_b = [0]
			xn_b = [0]
			bmnc = [0]
			bmns = [0]
		if type(pmns) is not type(None):
			pmns = np.zeros((mnmax))
		if not (rmns and zmnc and pmnc):
			rmns = np.zeros((mnmax))
			zmnc = np.zeros((mnmax))
			pmnc = np.zeros((mnmax))
		xm_out = []
		xn_out = []
		for mn in range(mnmax):
			if (abs(rmnc[mn])>=tol) or ((abs(zmns[mn])>=tol)) or (abs(rmns[mn])>=tol) or (abs(zmnc[mn])>=tol) or ((abs(pmnc[mn])>=tol)) or ((abs(pmns[mn])>=tol)):
				xm_out.append(xm[mn])
				xn_out.append(xn[mn])
		nharm = len(xm_out)
		weight = np.zeros((nharm))
		for i in range(nharm):
			if xm_out[i] > 0 and xn_out[i]!=0:
				weight[i] = float(abs(xn_out[i])+abs(xm_out[i]))/float(abs(xm_out[i]*xn_out[i]))
		#weight = float(abs(xn_out)+abs(xm_out))/float(abs(xm_out*xn_out))
		weight = np.where(weight==0,1.0,weight)
		weight = np.where(weight>1.0,1.0,weight)
		f=open(filename,'w')
		f.write(f'# Number of Harmonics\n')
		f.write(f'{int(nharm)}\n')
		f.write(f'# N,M,BMNC,BMNS,Weight\n')
		for mn in range(nharm):
			f.write(f'{int(xn_out[mn]):5d} {int(xm_out[mn]):5d} {0.0:10.2f} {0.0:10.2f} {weight[mn]:10.3f}\n')
		f.close()

	def plotConvergence(self,ax=None):
		"""Plots the FOCUS Convergence

		This routine plots the FOCUS code convergence for a given
		run.

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
		plotlabels = [r'Time [s]',r'$\chi^2$',r"$|d \chi^2 / d {\bf X}|$",r"$f_{B_n}$",
						r"$f_{B_{mn}}$",r"$f_{\Psi}$",r"f_{I~total}",r"$f_L$",
						r"$f_{COIL-SURF}$",r"$f_{curv}$",r"$f_{COIL-COIL}$",
						r"$f_{Torsion}$",r"$f_{nissin}$",r"$f_{B_{navg}}$",
						r"$f_{Straight}$",r"$f_{B_n}~Limiter$",r"$f_{B_{mn}~Limiter}$"]
		nele = len(self.evolution)
		print(nele)
		for i in range(nele):
			if i == 0: continue
			if any(self.evolution[i] > 0): ax.semilogy(self.evolution[0],self.evolution[i],label=plotlabels[i])
		ax.set_xlabel('Time [s]')
		ax.set_ylabel('Cost Function')
		ax.set_title('FOCUS Convergence')
		ax.legend()
		if lplotnow: pyplot.show()

	def plotBNormal(self,ax=None):
		"""Plots the FOCUS B-normal

		This routine plots the FOCUS code B-normal

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
		if self.IsSymmetric == 0:
			x = np.linspace(0,2*np.pi,int(self.Nzeta))
		elif self.IsSymmetric == 1:
			x = np.linspace(0,2*np.pi/self.Nfp,int(self.Nzeta))
		elif self.IsSymmetric == 2:
			x = np.linspace(0,np.pi/self.Nfp,int(self.Nzeta))
		else:
			return
		y = np.linspace(0,2.*np.pi,int(self.Nteta))
		hmesh=ax.pcolormesh(x,y,self.Bn.T,cmap='jet',shading='gouraud')
		#ax.semilogy(abscissa, data, **kwargs)
		ax.set_xlabel('Toroidal angle [rad]')
		ax.set_ylabel('Poloidal angle [rad]')
		ax.set_title('FOCUS B-Normal')
		pyplot.colorbar(hmesh,label='$B_n$ [T]',ax=ax)
		if lplotnow: pyplot.show()

	def plotBN3D(self,plot3D=None):
		"""Plots the FOCUS B-normal in 3D

		This routine plots the FOCUS code B-normal in 3D

		Parameters
		----------
		plot3D : plot3D object (optional)
			Plotting object to render to.
		"""
		import numpy as np
		import vtk
		from libstell.plot3D import PLOT3D
		# Handle optionals
		if plot3D: 
			lplotnow=False
			plt = plot3D
		else:
			lplotnow = True
			plt = PLOT3D()
		[points,triangles] = plt.torusvertexTo3Dmesh(self.xsurf.T,self.ysurf.T,self.zsurf.T,lcloseu=True,lclosev=False)
		# Handle Bn
		if hasattr(self,'Bn'):
			scalar = plt.valuesToScalar(self.Bn.flatten())
			# Add to Render
			plt.add3Dmesh(points,triangles,scalars=scalar)
		else:
			plt.add3Dmesh(points,triangles)
		# In case it isn't set by user.
		plt.setBGcolor()
		# Colorbar
		plt.colorbar()
		# Render if requested
		if lplotnow: plt.render()

	def plotLimiter3D(self,plot3D=None):
		"""Plots the FOCUS Limiter surface in 3D

		This routine plots the FOCUS code limiter surface in 3D

		Parameters
		----------
		plot3D : plot3D object (optional)
			Plotting object to render to.
		"""
		import numpy as np
		import vtk
		from libstell.plot3D import PLOT3D
		# Check to see if a limiter exists
		if not hasattr(self,'xsurf_lim'):
			return
		# Handle optionals
		if plot3D: 
			lplotnow=False
			plt = plot3D
		else:
			lplotnow = True
			plt = PLOT3D()
		[points,triangles] = plt.torusvertexTo3Dmesh(self.xsurf_lim.T,self.ysurf_lim.T,self.zsurf_lim.T,lcloseu=True,lclosev=False)
		# Handle Bn
		#scalar = plt.valuesToScalar(self.Bn.flatten())
		# Add to Render
		plt.add3Dmesh(points,triangles)
		# In case it isn't set by user.
		plt.setBGcolor()
		# Colorbar
		#plt.colorbar()
		# Render if requested
		if lplotnow: plt.render()

	def plotPoincare(self,ax=None):
		"""Plots the FOCUS Poincare Plot

		This routine plots the FOCUS Poincare plots

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
		ax.plot(self.ppr,self.ppz,'.',markersize=1)
		rsurf = np.sqrt(self.xsurf**2+self.ysurf**2)
		rmin = min(rsurf[0,:])
		rmax = max(rsurf[0,:])
		zmin = min(self.zsurf[0,:])
		zmax = max(self.zsurf[0,:])
		ax.plot(rsurf[0,:],self.zsurf[0,:],'r')
		if hasattr(self,'xsurf_lim'):			
			rsurf_lim = np.sqrt(self.xsurf_lim**2+self.ysurf_lim**2)
			rmin = min(rsurf_lim[0,:])
			rmax = max(rsurf_lim[0,:])
			zmin = min(self.zsurf_lim[0,:])
			zmax = max(self.zsurf_lim[0,:])
			ax.plot(rsurf_lim[0,:],self.zsurf_lim[0,:],'k')
		ax.set_xlabel('R [m]')
		ax.set_ylabel('Z [m]')
		ax.set_aspect('equal')
		ax.set_xlim(rmin*0.9,rmax*1.1)
		ax.set_ylim(zmin*1.1,zmax*1.1)
		if lplotnow: pyplot.show()

	def plotIota(self,ax=None):
		"""Plots the FOCUS Iota

		This routine plots the FOCUS rotational transform

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
		ax.plot(self.ppr[0]-self.ppr[0,1],self.iota)
		ax.set_xlabel('r [m]')
		ax.set_ylabel(r'$\iota$')
		ax.set_ylabel('Rotational Transform')
		if lplotnow: pyplot.show()

# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)

