#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for reading and handling 
THRIFT data.
"""

# Libraries
import numpy as np
import matplotlib.pyplot as plt

plt.rc('font', size=18)

# Constants

# THRIFT Class
class THRIFT():
    """" Class for working with THRIFT data
    
    """
    def __init__(self):
        test = 0
        
    def read_thrift(self,file):
        """Reads a THRIFT HDF5 file

		This routine reads and initilizes the THRIFT
		class with variable information from an HDF5 file.

		Parameters
		----------
		file : str
		Path to HDF5 file.
		"""

        import h5py
        import numpy as np

		#read file
        with h5py.File(file,'r') as f:
            # Logicals
            for temp in ['leccd','lohmic','lvmec','lnbcd']:
                if temp in f:
                    setattr(self, temp, np.int64(f[temp][0])==1)
            # Integers
            for temp in ['npicard','nrho','nssize','ntimesteps']:
                if temp in f:
                    setattr(self, temp, np.int64(f[temp][0]))
            # Floats
            for temp in ['VERSION','jtol','picard_factor']:
                if temp in f:
                    a = 1
                    setattr(self,temp, float(f[temp][0]))
            # Arrays                
            for temp in ['THRIFT_ALPHA1','THRIFT_ALPHA2','THRIFT_ALPHA3','THRIFT_ALPHA4','THRIFT_AMINOR',\
       			'THRIFT_BAV','THRIFT_BSQAV','THRIFT_BVAV','THRIFT_COEFF_A','THRIFT_COEFF_B','THRIFT_COEFF_BP',\
				'THRIFT_COEFF_C','THRIFT_COEFF_CP','THRIFT_COEFF_D','THRIFT_COEFF_DP','THRIFT_EPARB','THRIFT_ETAPARA','THRIFT_I',\
				'THRIFT_IBOOT','THRIFT_IECCD','THRIFT_INBCD','THRIFT_IOHMIC','THRIFT_IOTA','THRIFT_IPLASMA','THRIFT_ISOURCE',\
				'THRIFT_J','THRIFT_JBOOT','THRIFT_JECCD','THRIFT_JNBCD','THRIFT_JOHMIC','THRIFT_JPLASMA','THRIFT_JSOURCE',\
				'THRIFT_MATLD','THRIFT_MATMD','THRIFT_MATRHS','THRIFT_MATUD','THRIFT_P','THRIFT_PHIEDGE','THRIFT_PPRIME',\
				'THRIFT_RHO','THRIFT_RHOFULL','THRIFT_RMAJOR','THRIFT_S','THRIFT_S11','THRIFT_S12','THRIFT_SNOB','THRIFT_T',\
				'THRIFT_UGRID','THRIFT_VP']:
                if temp in f:
                    setattr(self, temp, np.array(f[temp][:]))
                    
        self.units_dictionary = {}
        for current in ['THRIFT_I','THRIFT_IBOOT','THRIFT_IECCD','THRIFT_INBCD','THRIFT_IOHMIC','THRIFT_IPLASMA','THRIFT_ISOURCE']:
            self.units_dictionary[current] = r'[A]'
        for current_density in ['THRIFT_J','THRIFT_JBOOT','THRIFT_JECCD','THRIFT_JNBCD','THRIFT_JOHMIC','THRIFT_JPLASMA','THRIFT_JSOURCE']:
            self.units_dictionary[current_density] = r'[A/m$^2]$'
        self.units_dictionary['THRIFT_ETAPARA'] = r'$[\Omega\,$m]'
                    
    def plot_vars_in_time(self,*vars,time_array=[0,1/4,1/2,3/4,1]):
        # plots var as a funciton of roa at different times
        # time_array is in fractions of t_end
        # vars is any variable of the type THRIFT_## with dimension (ntimesteps,nssize)
        
        for var in vars:
            plot_var = getattr(self,var)
            # check dimension of var is (ntimesteps,nssize)
            self.check_var_shape(plot_var,self.ntimesteps,self.nssize)
            t_end = self.THRIFT_T[-1]
            times = np.array(time_array)*t_end
            
            idx = [np.argmin(np.abs(self.THRIFT_T-t)) for t in times]
            times = self.THRIFT_T[idx]
            plot_var = plot_var[idx,:]
            
            _, ax = plt.subplots(figsize=(11,8))
            for it,time in enumerate(times):
                ax.plot(np.sqrt(self.THRIFT_S),plot_var[it,:],label=f't={time}s')
                ax.set_xlabel('r/a') 
                ax.set_title(var)   
            ax.grid()    
            ax.legend()
            try:
                ax.set_ylabel(self.units_dictionary[var])
            except:
                ax.set_ylabel('')
            # ax.set_yscale('log')
            plt.show()
            
    def check_var_shape(self,var, nt, nrho):
        
        if isinstance(var, np.ndarray):
            if var.ndim == 2:
                if var.shape != (nt, nrho):
                    print(f"{var} is a 2D array but has the wrong shape: {var.shape}")
                    exit(0)
            else:
                print(f"{var} is not a 2D array, it has {var.ndim} dimensions.")
                exit(0)
        else:
            print(f"{var} is not a numpy array.")
            exit(0)
            
    def get_vars(self,var,time=0.0):
        # return an array with var(roa) at t=time
        # var is any variable of the type THRIFT_## with dimension (ntimesteps,nssize)
        
        plot_var = getattr(self,var)
        
        idx = np.argmin(np.abs(self.THRIFT_T-time))
        
        real_time = self.THRIFT_T[idx]
        
        print(f'Returning variable {var} at t={real_time}s')
        
        return plot_var[idx,:]
        
        
        
            
    def plot_vars_vs_iota(self,*vars,time_array=[0,1/4,1/2,3/4,1]):
        # plots var as a funciton of iota at different times
        # time_array is in fractions of t_end
        # vars is any variable of the type THRIFT_## with dimension (ntimesteps,nssize)
        
        for var in vars:
            plot_var = getattr(self,var)
            # check dimension of var is (ntimesteps,nssize)
            self.check_var_shape(plot_var,self.ntimesteps,self.nssize)
            t_end = self.THRIFT_T[-1]
            times = np.array(time_array)*t_end
            
            idx = [np.argmin(np.abs(self.THRIFT_T-t)) for t in times]
            times = self.THRIFT_T[idx]
            plot_var = plot_var[idx,:]
            
            iota = self.THRIFT_IOTA[idx,:]
            
            _, ax = plt.subplots(figsize=(11,8))
            for it,time in enumerate(times): 
                ax.plot(iota[it,:],plot_var[it,:],label=f't={time}s')
                ax.set_xlabel('iota') 
                ax.set_title(var)   
            ax.grid()    
            ax.legend()
            try:
                ax.set_ylabel(self.units_dictionary[var])
            except:
                ax.set_ylabel('')
            plt.show()

            
    # def compare_J(self,time_array=[0,1/4,1/2,3/4,1]):
        
    #     J1 = getattr(self,'THRIFT_J')
    #     J2B = getattr(self,'THRIFT_JDOTB')
    #     B = getattr(self,'THRIFT_BAV')
    #     J2 = J2B/B
        
    #     t_end = self.THRIFT_T[-1]
    #     times = np.array(time_array)*t_end
        
    #     idx = [np.argmin(np.abs(self.THRIFT_T-t)) for t in times]
    #     times = self.THRIFT_T[idx]
        
    #     J1 = J1[idx,:]
    #     J2 = J2[idx,:]
        
    #     print(np.max(J1-J2))
        
    #     _, ax = plt.subplots(figsize=(11,8))
        
    #     for it,time in enumerate(times):
    #         ax.plot(np.sqrt(self.THRIFT_S),J1[it,:],'-',label=f't={time}s')
    #         ax.plot(np.sqrt(self.THRIFT_S),J2[it,:],'.',label=f't={time}s')
            
    #     ax.set_xlabel('r/a')
    #     ax.set_ylabel('')
    #     #ax.set_title(var)   
    #     ax.grid()    
    #     ax.legend()
    #     plt.show()
        
        
        
        
# # THRIFT Input Class
# class THRIFT_INPUT():

# 	def __init__(self, parent=None):
# 		self.libStell = LIBSTELL()

# 	def read_input(self,filename):
		
# 		# Not yet implemented
# 		#indata_dict = self.libStell.read_thrift_input(filename)
# 		#for key in indata_dict:
# 		#	setattr(self, key, indata_dict[key])

# 	def write_input(self,filename):
		
# 		# Not yet implemented
# 		#out_dict = vars(self)
# 		#self.libStell.write_thrift_input(filename,out_dict)

# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)