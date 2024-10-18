##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for reading and handling 
STELLOPT data.
"""

# Libraries
from libstell.libstell import LIBSTELL

# Constants

# STELLOPT Class
class STELLOPT():
	"""Class for working with STELLOPT data

	"""
	def __init__(self):
		test = 0
		self.target_names=['TEST_X',      \
			'TEST_Y', 'ROSENBROCK_X', 'PHIEDGE', 'CURTOR',             \
			'CURTOR_MAX', 'RBTOR', 'R0', 'Z0', 'B0', 'VOLUME', 'BETA', \
			'BETAPOL', 'BETATOR', 'WP', 'ASPECT', 'CURVATURE', 'KAPPA',\
			'KAPPA_BOX', 'KAPPA_AVG', 'ASPECT_MAX', 'PMIN', 'EXTCUR',  \
			'LINE_NE', 'LINE_TE', 'LINE_TI', 'LINE_ZEFF',              \
			'XICS_BRIGHT', 'XICS', 'XICS_W3', 'XICS_V', 'SXR',         \
			'FARADAY', 'LINE_VISBREM', 'ECE', 'PRESS', 'PRESSPRIME',   \
			'NE', 'TE', 'TI', 'VPHI', 'IOTA', 'VACIOTA', 'JDOTB',      \
			'MAGWELL', 'BMIN', 'BMAX', 'MSE', 'JCURV', 'COILLEN',      \
			'COILSEP', 'COILCRV', 'COILSELF', 'COILTORVAR',            \
			'COILRECT', 'COILPOLY', 'BPROBES', 'FLUXLOOPS', 'SEGROG',  \
			'VESSEL', 'SPEARATRIX', 'LIMITER', 'BALLOON', 'BOOTSTRAP', \
			'NEO', 'DKES', 'DKES_ERDIFF', 'DKES_ALPHA', 'TXPORT',      \
			'ORBIT', 'HELICITY', 'HELICITY_FULL', 'JSTAR', 'RESJAC',   \
			'COIL_BNORM', 'REGCOIL_CHI2_B', 'CURVATURE_P2', 'GAMMA_C', \
			'KINK']

	def read_stellopt_map(self,filename='map.dat'):
		"""Reads a STELLOPT MAP output file

		This routine reads the STELLOPT map.dat  output file.

		Parameters
		----------
		filename : str
			Path to STELLOPT file.
		"""
		import numpy as np
		import re
		f = open(filename,'r')
		content = f.read()
		f.close()
		numbers = re.findall(r'-?\d+\.?\d*(?:[eE][+-]?\d+)?', content)
		numbers = [float(num) for num in numbers]
		mtargets  = int(numbers[0])
		nvars     = int(numbers[1])
		ndiv      = int(numbers[2])
		numsearch = int(numbers[3])
		nnext     = nvars * numsearch + 4
		x         = numbers[4:nnext]
		fval      = numbers[nnext:]
		x2d       = np.reshape(x,(numsearch,nvars)).T
		f2d       = np.reshape(fval,(numsearch,mtargets)).T
		self.x_map = x2d
		self.f_map = f2d

	def read_stellopt_profile(self,filename):
		"""Reads a STELLOPT tprof output file

		This routine reads the STELLOPT tprof output file.

		Parameters
		----------
		file : str
			Path to tprof file.
		Returns
		-------
		s : ndarray
			Array of normalized toroidal flux (s).
		ne : ndarray
			Array of electron density [m^-3]
		te : ndarray
			Array of electron temperature [eV]
		ti : ndarray
			Array of ion temperatures [eV]
		zeff : ndarray
			Effective ion charge.
		p : ndarray
			Pressure [Pa]
		"""
		import numpy as np
		f = open(filename,'r')
		line = f.readline() # header
		s    = []
		ne   = []
		te   = []
		ti   = []
		zeff = []
		p    = []
		for line in f: # read rest of lines
			[txt1,txt2,txt3,txt4,txt5,txt6] = line.split()
			s.append(float(txt1))
			ne.append(float(txt2))
			te.append(float(txt3))
			ti.append(float(txt4))
			zeff.append(float(txt5))
			p.append(float(txt6))
		f.close()
		return np.array(s),np.array(ne),np.array(te),np.array(ti),np.array(zeff),np.array(p)

	def read_stellopt_varlabels(self,filename='var_labels'):
		"""Reads a STELLOPT var_labels output file

		This routine reads the STELLOPT var_labels output file.

		Parameters
		----------
		file : str
			Path to var_labels file. (default 'var_labels')
		"""
		f = open(filename,'r')
		line = f.readline()
		nvars = int(line)
		var      = []
		varnames = []
		for i in range(nvars):
			line = f.readline()
			line.replace('\n','')
			[var1,var2] = line.split(':')
			var.append(var1.strip())
			varnames.append(var2.strip())
		line = f.readline()
		mtargets = int(line)
		targetnames = []
		for i in range(mtargets):
			line = f.readline()
			line.replace('\n','')
			targetnames.append(line.strip())
		f.close()
		self.varnames = varnames
		self.var = var
		self.targetnames = targetnames

	def read_stellopt_jacobian(self,filename):
		"""Reads a STELLOPT jacobian output file

		This routine reads the STELLOPT jacobian output file.

		Parameters
		----------
		file : str
			Path to jacobian file.
		"""
		import numpy as np
		import re
		f = open(filename,'r')
		content = f.read()
		f.close()
		numbers = re.findall(r'-?\d+\.?\d*(?:[eE][+-]?\d+)?', content)
		numbers = [float(num) for num in numbers]
		mtargets  = int(numbers[0])
		nvars     = int(numbers[1])
		jac       = numbers[2:]
		jac2d       = np.reshape(jac,(nvars,mtargets)).T
		self.jac2d  = jac2d

	def read_stellopt_output(self,filename):
		"""Reads a STELLOPT output file

		This routine reads the STELLOPT output file.

		Parameters
		----------
		file : str
			Path to STELLOPT file.
		"""
		import numpy as np
		f = open(filename,'r')
		lines = f.readlines()
		f.close()
		i = 0
		temp_dict={}
		self.niter = len([line for line in lines if 'ITER' in line])
		iter_val = -1
		while (i < len(lines)):
			if 'VERSION' in lines[i]:
				[temp,version_txt] = lines[i].split()
				self.stellopt_version = float(version_txt)
				i = i +1
				continue
			elif 'ITER' in lines[i]:
				# Note this fixes issue where ITER is not necessarily
				# sequential
				[temp,iter_txt] = lines[i].split()
				if 'ITER' not in temp_dict.keys():
					temp_dict['ITER'] = np.zeros((self.niter,1))
				iter_val = iter_val + 1
				temp_dict['ITER'][iter_val] = int(iter_txt)
				i = i +1
				continue
			elif 'TARGETS' in lines[i]:
				[targ_name,nrow_txt,ncol_txt] = lines[i].split()
				nrow = int(nrow_txt)
				ncol = int(ncol_txt)
				i1 = i + 2 # skip header
				if 'TARGETS' not in temp_dict.keys():
					temp_dict['TARGETS'] = np.zeros((self.niter,nrow))
				for j in range(nrow):
					temp_txt = lines[i1+j]
					temp_dict['TARGETS'][iter_val,j] = float(temp_txt)
				i = i + 2 + nrow - 1
			elif 'SIGMAS' in lines[i]:
				[targ_name,nrow_txt,ncol_txt] = lines[i].split()
				nrow = int(nrow_txt)
				ncol = int(ncol_txt)
				i1 = i + 2 # skip header
				if 'SIGMAS' not in temp_dict.keys():
					temp_dict['SIGMAS'] = np.zeros((self.niter,nrow))
				for j in range(nrow):
					temp_txt = lines[i1+j]
					temp_dict['SIGMAS'][iter_val,j] = float(temp_txt)
				i = i + 2 + nrow - 1
			elif 'VALS' in lines[i]:
				[targ_name,nrow_txt,ncol_txt] = lines[i].split()
				nrow = int(nrow_txt)
				ncol = int(ncol_txt)
				i1 = i + 2 # skip header
				if 'VALS' not in temp_dict.keys():
					temp_dict['VALS'] = np.zeros((self.niter,nrow))
				for j in range(nrow):
					temp_txt = lines[i1+j]
					temp_dict['VALS'][iter_val,j] = float(temp_txt)
				i = i + 2 + nrow - 1
			if any(x in lines[i] for x in self.target_names):
				[targ_name,nrow_txt,ncol_txt] = lines[i].split()
				nrow = int(nrow_txt)
				ncol = int(ncol_txt)
				# Header line i+1
				header = lines[i+1].split()
				# Get data
				i1 = i+2
				#print(f'---- TARGET: {targ_name} {nrow} {ncol}')
				for j in range(nrow):
					temp_txt = lines[i1+j].split()
					for k in range(ncol):
						# Fix the header name if bad.
						header_fix = header[k].replace('#','K')
						header_fix = header_fix.replace('|B|','MODB')
						header_fix = header_fix.replace('EPS_EFF^(3/2)','EPS_EFF32')
						header_fix = header_fix.replace('<B**2>','BSQAVG')
						header_fix = header_fix.replace('+','p')
						header_fix = header_fix.replace('-','m')
						header_fix = header_fix.replace('/','')
						#print(f'------ HEADER: {header_fix}')
						targ_name_full = targ_name+'_'+header_fix
						if targ_name_full not in temp_dict.keys():
							temp_dict[targ_name_full] = np.zeros((self.niter,nrow))
						temp_dict[targ_name_full][iter_val,j] = float(temp_txt[k])
				# Add VAL if not explicitly named
				if 'VAL' not in header:
					if targ_name+'_VAL' not in temp_dict.keys():
						temp_dict[targ_name+'_VAL'] = np.zeros((self.niter,nrow))
					k = header.index('SIGMA')+1
					for j in range(nrow):
						temp_txt = lines[i1+j].split()
						temp_dict[targ_name+'_VAL'][iter_val,j] = float(temp_txt[k])
				i = i + 2 + nrow - 1
			i = i + 1
		#print(temp_dict.keys())
		# Convert to attributes
		for key in temp_dict:
			setattr(self, key, temp_dict[key])
		# Calculate Chisq for each value
		for targ_name in self.target_names:
			if hasattr(self,targ_name+'_TARGET'):
				targ  = getattr(self,targ_name+'_TARGET')
				sigma = getattr(self,targ_name+'_SIGMA')
				val   = getattr(self,targ_name+'_VAL')
				chisq = (targ-val)/sigma
				setattr(self,targ_name+'_CHISQ',chisq*chisq)

# STELLOPT Input Class
class STELLOPT_INPUT():
	"""Class for working with STELLOPT INPUT data

	"""
	def __init__(self, parent=None):
		self.libStell = LIBSTELL()
		self.global_data = STELLOPT_INPUT_GLOBAL()
		self.var_data = STELLOPT_INPUT_VAR()
		self.target_data = STELLOPT_INPUT_TARGET()

	def read_input(self,filename):
		"""Reads STELLOPT_INPUT namelist from a file

		This routine wrappers the stellopt_input_mod module reading routine.
		Parameters
		----------
		filename : string
			Input file name with OPTIMUM namelist
		"""
		# there are three separate module we need to deal with
		global_dict, var_dict, target_dict = self.libStell.read_stellopt_input(filename)
		for key in global_dict:
			setattr(self.global_data, key, global_dict[key])
		for key in var_dict:
			setattr(self.var_data, key, var_dict[key])
		for key in target_dict:
			setattr(self.target_data, key, target_dict[key])

	def write_input(self,filename):
		"""Writes STELLOPT_INPUT namelist to a file

		This routine wrappers the stellopt_input_mod module writing routine.
		Parameters
		----------
		filename : string
			Input file name to write OPTIMUM namelist to
		"""
		# there are three separate module we need to deal with
		global_dict = vars(self.global_data)
		var_dict = vars(self.var_data)
		target_dict = vars(self.target_data)
		self.libStell.write_stellopt_input(filename,global_dict,var_dict,target_dict)

class STELLOPT_INPUT_GLOBAL():
	def __init__(self):
		pass

class STELLOPT_INPUT_VAR():
	def __init__(self):
		pass

class STELLOPT_INPUT_TARGET():
	def __init__(self):
		pass




# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)



