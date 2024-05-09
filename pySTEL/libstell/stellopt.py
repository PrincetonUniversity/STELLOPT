##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for reading and handling 
STELLOPT data.
"""

# Libraries
#from libstell import libstell

# Constants

# STELLOPT Class
class STELLOPT():
	"""Class for working with BOOTSJ data

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
		while (i < len(lines)):
			if 'VERSION' in lines[i]:
				[temp,version_txt] = lines[i].split()
				self.stellopt_version = float(version_txt)
				i = i +1
				continue
			elif 'ITER' in lines[i]:
				[temp,iter_txt] = lines[i].split()
				iter_val = int(iter_txt)
				#print(f'-- ITERATION: {lines[i]}')
				i = i +1
				continue
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


