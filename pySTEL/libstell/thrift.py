#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for reading and handling 
THRIFT data.
"""

# Libraries
import numpy as np
import matplotlib.pyplot as plt
import h5py

plt.rc('font', size=18)
default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
custom_colors = ['#5faf30', '#1D2258', '#004817', '#a1cdc8']
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=custom_colors+default_colors)
plt.rcParams['lines.linewidth'] = 2.5
# plt.rcParams['axes.prop_cycle'] = plt.cycler(color=['#5faf30','#1D2258','#004817','#a1cdc8'])

# Constants

# THRIFT Class
class THRIFT():
    """" Class for working with THRIFT data
    
    """
    def __init__(self):
        test = 0
        
    def read_thrift(self,*files):
        """Reads THRIFT HDF5 files

		This routine reads and initilizes the THRIFT
		class with variable information from one or more HDF5 files.

		Parameters
		----------
		files : str
		Path to HDF5 files.
		"""
        
        # checks order of time is correct
        self.check_time_order(*files)
        
        # set integers attributes
        self.set_integers_attribute(*files)
        
        # set arrays attributes
        self.set_arrays_attribute(*files) 
        
        # set units of quantities
        self.set_units()       
    
    def check_time_order(self,*files):

        time = []
        
        for file in files:
            with h5py.File(file,'r') as f:
                time.append( f['THRIFT_T'][:] )
        time = np.concatenate(time)
        
        #check ordering
        if(not np.all(np.diff(time) > 0) ):
            print('ERROR: the given THRIFT output files are not in the correct order...')
            exit(0)
            
    def set_integers_attribute(self,*files):
        # checks that nssize didn't change among different ouput files
        # computes the total number of ntimesteps
        
        npicard = []
        nssize = []
        ntimesteps = []
        
        for file in files: 
            with h5py.File(file,'r') as f:
                npicard.append(np.int64(f['npicard'][0]))
                nssize.append(np.int64(f['nssize'][0]))
                ntimesteps.append(np.int64(f['ntimesteps'][0]))

        if(not np.all(nssize==nssize[0])):
            print('ERROR: nssize changes bewteen output files!')
            exit(0)
            
        setattr(self,'ntimesteps',np.sum(ntimesteps))
        setattr(self,'nssize',nssize[0])
        setattr(self,'npicard',npicard)
            
    def set_arrays_attribute(self,*files):
        
        for file in files: 
            with h5py.File(file,'r') as f:
                # Arrays
                for temp in ['THRIFT_ALPHA1','THRIFT_ALPHA2','THRIFT_ALPHA3','THRIFT_ALPHA4','THRIFT_AMINOR',\
       			    'THRIFT_BAV','THRIFT_BETATOT','THRIFT_BSQAV','THRIFT_BVAV','THRIFT_COEFF_A','THRIFT_COEFF_B','THRIFT_COEFF_BP',\
				    'THRIFT_COEFF_C','THRIFT_COEFF_CP','THRIFT_COEFF_D','THRIFT_COEFF_DP','THRIFT_EPARB','THRIFT_ER','THRIFT_ETAPARA','THRIFT_GNEO',\
                    'THRIFT_I','THRIFT_IBOOT','THRIFT_IECCD','THRIFT_INBCD','THRIFT_IOHMIC','THRIFT_IOTA','THRIFT_IPLASMA','THRIFT_ISOURCE',\
				    'THRIFT_J','THRIFT_JBOOT','THRIFT_JECCD','THRIFT_JNBCD','THRIFT_JOHMIC','THRIFT_JPLASMA','THRIFT_JSOURCE',\
				    'THRIFT_MATLD','THRIFT_MATMD','THRIFT_MATRHS','THRIFT_MATUD','THRIFT_P','THRIFT_PHIEDGE','THRIFT_PPRIME',\
				    'THRIFT_QNEO','THRIFT_RMAJOR','THRIFT_S11','THRIFT_S12','THRIFT_T', 'THRIFT_UGRID','THRIFT_VP',\
                    'THRIFT_DENS', 'THRIFT_TEMP', 'THRIFT_PRESS']:
                    if temp in f:
                        # Get the data from the file
                        data = np.array(f[temp][:])

                        # Check if the attribute exists; if not, initialize it
                        if not hasattr(self, temp):
                            setattr(self, temp, data)
                        else:
                            # Concatenate the new data to the existing attribute
                            existing_data = getattr(self, temp)
                            setattr(self, temp, np.concatenate((existing_data, data)))
                        
        # set THRIFT_S array (no concatenation needed); in set_integers_attribute already checked nssize is the same for ALL files
        with h5py.File(file,'r') as f:
            setattr(self, 'THRIFT_S', np.array(f['THRIFT_S'][:]))
            setattr(self, 'THRIFT_SNOB', np.array(f['THRIFT_SNOB'][:]))
                        
    def set_units(self):
        # sets units of different variables (useful when plotting)
        
        self.units_dictionary = {}
        for current in ['THRIFT_I','THRIFT_IBOOT','THRIFT_IECCD','THRIFT_INBCD','THRIFT_IOHMIC','THRIFT_IPLASMA','THRIFT_ISOURCE']:
            self.units_dictionary[current] = r'[A]'
        for current_density in ['THRIFT_J','THRIFT_JBOOT','THRIFT_JECCD','THRIFT_JNBCD','THRIFT_JOHMIC','THRIFT_JPLASMA','THRIFT_JSOURCE']:
            self.units_dictionary[current_density] = r'[A/m$^2]$'
        self.units_dictionary['THRIFT_ETAPARA'] = r'$[\Omega\,$m]'
        self.units_dictionary['THRIFT_ER'] = r'$[V/$m]'
        self.units_dictionary['THRIFT_GNEO'] = r'[m$^{-2}s$^{-1}$]'
        self.units_dictionary['THRIFT_QNEO'] = r'[$\text{eV}~\text{m}^{-2}s$^{-1}$]'
             
    def plot_vars_in_time(self,*vars,time_slice=None,time_array=None):
        # plots var as a funciton of roa at different times
        # the times can be given as time_slices (fractions of t_end)
        # or as time_array
        # if both given, time_slice prevails
        # vars is any variable of the type THRIFT_## with dimension (ntimesteps,nssize)
        
        if( (time_array is None) or (time_slice is not None and time_array is not None)):
            t_end = self.THRIFT_T[-1]
            time_slice = [0,1/4,1/2,3/4,1]
            times = np.array(time_slice)*t_end
        else:
            times = np.atleast_1d(time_array)
        
        for var in vars:
            plot_var = getattr(self,var)
            # check dimension of var is (ntimesteps,nssize)
            self.check_var_shape(plot_var,self.ntimesteps,self.nssize)

            idx = [np.argmin(np.abs(self.THRIFT_T-t)) for t in times]
            times = self.THRIFT_T[idx]
            plot_var = plot_var[idx,:]
            
            _, ax = plt.subplots(figsize=(11,8))
            for it,time in enumerate(times):
                try:
                    # ax.plot(np.sqrt(self.THRIFT_S),plot_var[it,:],label=f't={time}s')
                    ax.plot(np.sqrt(self.THRIFT_S),plot_var[it,:],label=f't={time:.1f}s'+r', $\beta=$'+f'{self.THRIFT_BETATOT[idx[it]]*100:.2f}%')   
                except:
                    ax.plot(np.sqrt(self.THRIFT_SNOB),plot_var[it,:],label=f't={time}s')
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
        # return an array with var(roa) at t=time; time can be an array
        # var is any variable of the type THRIFT_## with dimension (ntimesteps,nssize)
        
        plot_var = getattr(self,var)
        
        time = np.atleast_1d(time)
        
        idx = [np.argmin(np.abs(self.THRIFT_T-t)) for t in time]
        #idx = np.argmin(np.abs(self.THRIFT_T-time))
        
        real_time = self.THRIFT_T[idx]
        
        print(f'Returning variable {var} at t={real_time}s')
        
        if(plot_var.ndim == 1):
            return plot_var[idx]
        elif(plot_var.ndim == 2):
            return plot_var[idx,:]
        elif( plot_var.ndim ==3):
            return plot_var[idx,:,:]
        else:
            print('ERROR: What are you trying to get?!')
            exit(0)

    def plot_plasma_profile(self,plasma_file):
        
        import h5py
        
        hf = h5py.File(plasma_file, 'r')

        #read taxis and raxis
        nt = np.array( hf['nt'] )
        raxis = np.array( hf['raxis_prof'][:] )
        taxis = np.array( hf['taxis_prof'][:] )
        
        # select some t's for plotting
        num_values = 1000
        step = max(1, len(taxis) // num_values)
        selected_indices = np.arange(0, len(taxis), step)[:num_values]
        t_selected = taxis[selected_indices]
        
        ne = np.array( hf['ne_prof'][:] )
        Te = np.array( hf['te_prof'][:] )
        
        ni = np.array( hf['ni_prof'][:] )
        Ti = np.array( hf['ti_prof'][:] )
        
        nion = np.int64( hf['nion'] )
        
        hf.close()
        
        _, ax_n = plt.subplots(figsize=(11,8))
        _, ax_T = plt.subplots(figsize=(11,8))    
        
        ax_n.plot(raxis,ne[:,selected_indices]/1e20)
        ax_n.set_xlabel('r/a') 
        ax_n.set_title('ne [1E20 m^-3]')   
        ax_n.grid()   
        
        ax_T.plot(raxis,Te[:,selected_indices]/1e3)
        ax_T.set_xlabel('r/a') 
        ax_T.set_title('Te [keV]')   
        ax_T.grid()  
        
        plt.show()

        for i in range(nion):
            
            _, ax_n = plt.subplots(figsize=(11,8))
            _, ax_T = plt.subplots(figsize=(11,8))
            
            ax_n.plot(raxis,ni[:,selected_indices,i]/1e20)
            ax_n.set_xlabel('r/a') 
            ax_n.set_title(f'ni [1E20 m^-3], ion={i+1}')   
            ax_n.grid()
            
            ax_T.plot(raxis,Ti[:,selected_indices,i]/1e3)
            ax_T.set_xlabel('r/a') 
            ax_T.set_title(f'Ti [keV], ion={i+1}')   
            ax_T.grid()
            
            plt.show()
            
    def plot_plasma_current_decay(self,tstart=20.0):
        # plots total plasma current as a function of time
        # estimates decay time with LR circuit eqvalent time-scale
        
        Iplasma = self.THRIFT_IPLASMA[:,-1]
        t = self.THRIFT_T
        
        _, ax = plt.subplots(figsize=(11,8))
        ax.plot(t,Iplasma)
        ax.set_xlabel('t [s]') 
        ax.set_title('Total Plasma Current')   
        ax.grid()    
        # ax.legend()
        ax.set_ylabel('[A]')
        # plt.show()
        
        _, ax = plt.subplots(figsize=(13,8))
        ax.plot(t,np.abs(Iplasma))
        ax.set_xlabel('t [s]')    
        ax.grid()    
        ax.set_ylabel('[A]')
        ax.set_yscale('log')
        # plt.show()
        
        t_fit = t[t>tstart]
        I_fit = np.abs(Iplasma[t>tstart])
        
        p1,p0 = np.polyfit(t_fit,np.log(I_fit),1)
        
        tau = -1/p1
        
        ax.plot(t_fit,np.exp(p0+p1*t_fit),'--',label=r'$\tau_{\text{fit}}=$'+f'{tau:.1f}s')
        ax.legend()
        # plt.show()
        
        mu0 = 1.256e-6
        R0 = np.mean(self.THRIFT_RMAJOR)
        a  = np.mean(self.THRIFT_AMINOR[:,-1])

        Lext = mu0*R0*(np.log(8*R0/a)-2.0)
        
        A = np.pi*a*a
        L = 2*np.pi*R0
        
        integrated_cond = np.trapz(1/self.THRIFT_ETAPARA[-1,:],self.THRIFT_S)
        Rohm = (L/A) * 1/integrated_cond
        
        tau_LR = Lext / Rohm
        
        # ax.text(60,4e3,r'$\tau_{L/R}=$'+f'{tau_LR:.1f}s')
        ax.set_title('|Total Plasma Current|, '+r'$\tau_{L/R}=$'+f'{tau_LR:.1f}s')
        plt.show()
        
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
            
    def plot_all_ambipolar_roots(self,*ambipolar_files,make_plot=True):
        # plots all roots of the electric field 
        # each file correspond to a different instant of time
        
        from itertools import groupby
        from collections import defaultdict

        for file in ambipolar_files:
            
            try:
                time = np.loadtxt(file,skiprows=2,max_rows=1)
            except:
                time = -1
                
            penta = np.loadtxt(file,skiprows=4)
        
            roa = penta[:,0]
            Er = penta[:,1]
            JBS = penta[:,2]
            
            Gamma_e = penta[:,3]
            
            # self.roa = defaultdict(list)
            # self.Er = defaultdict(list)
            
            num_roots = []
            roa_all = []
            Er_all = []
            JBS_all = []
            Gamma_e_all = []
            
            paired  = zip(roa,Er,JBS,Gamma_e)
            
            # Group by the first element (roa)
            for _, group in groupby(paired, key=lambda x: x[0]):
                
                group_list = list(group)  # Convert the group to a list
                
                num_roots.append(len(group_list))
                
                roa_all.append( [x[0] for x in group_list])  # Extract the roa part of the group
                Er_all.append(  [x[1] for x in group_list])  # Extract the other_array part of the group
                JBS_all.append(  [x[2] for x in group_list])  # Extract the other_array part of the group
                Gamma_e_all.append(  [x[3] for x in group_list])  # Extract the other_array part of the group
            
            # plots    
            _, ax = plt.subplots(figsize=(11,8))

            # Plot each root
            for j in range(max(num_roots)):
                x = []  # roa values
                y = []  # Er values
                for k, (r_vals, er_vals) in enumerate(zip(roa_all, Er_all)):
                    if j < len(er_vals):  # Only include if the j-th value exists in Er[k]
                        x.append(r_vals[0])  
                        y.append(er_vals[j])
                plt.plot(x, y, marker='o')  # Plot the j-th curve

            
            ax.set_xlabel('r/a')
            ax.set_ylabel(r'$E_r$ [V/cm]')
            ax.set_title(f't={time}s')
            
            #JBS plot
            _, ax = plt.subplots(figsize=(11,8))

            # Plot each root
            for j in range(max(num_roots)):
                x = []  # roa values
                y = []  # JBS values
                for k, (r_vals, jbs_vals) in enumerate(zip(roa_all, JBS_all)):
                    if j < len(jbs_vals):  # Only include if the j-th value exists in Er[k]
                        x.append(r_vals[0])  
                        y.append(jbs_vals[j])
                plt.plot(x, y, marker='o')  # Plot the j-th curve

            
            ax.set_xlabel('r/a')
            ax.set_ylabel(r'$J_{BS}~[A/m^2]$')
            ax.set_title(f't={time}s')
            
            #Fluxes plot
            _, ax = plt.subplots(figsize=(11,8))

            # Plot each root
            for j in range(max(num_roots)):
                x = []  # roa values
                y = []  # gamma_e values
                for k, (r_vals, ge_vals) in enumerate(zip(roa_all, Gamma_e_all)):
                    if j < len(ge_vals):  # Only include if the j-th value exists in Er[k]
                        x.append(r_vals[0])  
                        y.append(ge_vals[j])
                plt.plot(x, y, marker='o')  # Plot the j-th curve
  
            ax.set_xlabel('r/a')
            ax.set_ylabel(r'$\Gamma_e~[m^{-2}~s^{-1}]$')
            ax.set_title(f't={time}s')
            
        if(make_plot): plt.show()
        
        return roa_all, Er_all, JBS_all, Gamma_e_all
        
    def Maxwell_construction(self,fluxes_vs_Er_file,Zions):
        # plots fluxes*Z as function of Er
        
        import matplotlib.pyplot as plt
        from collections import defaultdict
        
        time = np.loadtxt(fluxes_vs_Er_file,skiprows=2,max_rows=1)
            
        penta = np.loadtxt(fluxes_vs_Er_file,skiprows=5)
        
        roa = penta[:,0]
        Er = penta[:,1]
        gamma_e = penta[:,2]      
        
        #check dimension of Zion equals that of file 
        if (len(penta[0,3:]) != len(Zions) ):
            print('ERROR: length of Zions does not match file size')
            exit(0)
            
        gamma_i_tot = np.sum(penta[:,3:]*Zions,axis=1)
        
        gamma_i = penta[:,3:]
        
        Jr = gamma_i_tot - gamma_e
        
        Er_dict = defaultdict(list)
        Jr_dict = defaultdict(list)
        Ge_dict = defaultdict(list)
        Gi_dict = defaultdict(lambda: defaultdict(list))
        
        for r,er,jr,ge in zip(roa,Er,Jr,gamma_e):
            Er_dict[r].append(er)
            Jr_dict[r].append(jr)
            Ge_dict[r].append(ge)
        
        for k,_ in enumerate(Zions):
            for r, gi in zip(roa,gamma_i[:,k]):
                Gi_dict[k][r].append(gi)
        
        roa_unique = np.unique(roa)
        # check size of dictionaries correspond to size of roa_unique
        if( (len(roa_unique) != len(Er_dict)) or (len(roa_unique) != len(Jr_dict)) ):
            print('ERROR" Sizes of roa_unique and of dictionaries should match...')
            exit(0)
              
        Er_ambipolar = []
            
        # set the ambipolar root for each r/a using Maxwell construction criterium
        for roa in roa_unique:
            sign_changes = np.diff(np.sign(Jr_dict[roa])) != 0
            change_indices = np.where(sign_changes)[0]
            num_sign_changes = len(change_indices)
            
            if(num_sign_changes==1):
                Er_val = Er_dict[roa][change_indices[0]]
            elif( num_sign_changes>1):
                # integrate between first and last root
                first_change = change_indices[0]
                last_change = change_indices[-1]
                integral = np.trapz(Jr_dict[roa][first_change:last_change+1],Er_dict[roa][first_change:last_change+1])
                if(integral>0):
                    Er_val = Er_dict[roa][change_indices[0]]
                else:
                    Er_val = Er_dict[roa][change_indices[-1]]  
            else:
                print('ERROR: THis else should not be possible....')
                exit(0)
            
            Er_ambipolar.append(Er_val)
        
        Er_ambipolar = np.array(Er_ambipolar)
            
        # plot
        # _, ax = plt.subplots(figsize=(11,8))
        plt.plot(roa_unique,Er_ambipolar,label='Maxwell')
        # plt.set_xlabel('r/a')
        # plt.set_ylabel(r'$E_r$ [V/cm]')
        # plt.set_title(f't={time}s')
        plt.legend()
        plt.show()
        
        _, ax = plt.subplots(figsize=(11,8))
        ax.plot(Er_dict[roa_unique[4]],Jr_dict[roa_unique[4]],'.-')
        ax.set_xlabel('Er [V/cm]')
        ax.set_ylabel(r'$\Sigma Z_i\Gamma_i-\Gamma_e$')
        ax.set_title(f'r/a={roa_unique[4]}')
        ax.grid()
        
        _, ax = plt.subplots(figsize=(11,8))
        ax.plot(Er_dict[roa_unique[4]],Ge_dict[roa_unique[4]],'.-',label='Gamma_e')
        for k,_ in enumerate(Zions):
            ax.plot(Er_dict[roa_unique[4]],Gi_dict[k][roa_unique[4]],'.-',label=f'Gamma_i_{k}')
        ax.set_xlabel('Er [V/cm]')
        ax.set_ylabel(r'$\Gamma$')
        ax.set_yscale('symlog',linthresh=0.1)
        ax.set_title(f'r/a={roa_unique[4]}')
        ax.grid()
        plt.legend()
        
        plt.show()
        
# THRIFT Class
class THRIFT_plasma_solver():
    """" Class for working with plasma solver implemented in THRIFT
    
    """
    
    def __init__(self, plasma=None, list_of_species=None):
        
        import sys 
        sys.path.insert(1,'/home/antonio/STELLOPT/pySTEL/libstell')
        from plasma import PLASMA
        
        # should give plasma OR list_of_species. If both, plasma prevails
        if plasma is None and list_of_species is None:
            print('Could not create class: Need to provide a plasma class or a list with name of species')
            exit(1)
        elif plasma is None and list_of_species is not None:
            self.plasma_class = PLASMA(list_of_species=list_of_species)
        else:
            self.plasma_class = plasma
        
    def read_thrift_plasma_solver(self,file):
        """Reads plasma_solver THRIFT HDF5 file

		Parameters
		----------
		file : str
		Path to HDF5 files.
		"""
        import h5py

        with h5py.File(file,'r') as f:
                self.Nt = f['Nt_plasma_grid']
                self.Nr = f['Nr_plasma_grid']
                #
                self.time_grid = f['time_plasma_grid'][:]
                self.rho_grid  = f['rho_plasma_grid'][:]
                #
                self.plasma_N = f['plasma_N'][:,:,:]
                self.plasma_T = f['plasma_T'][:,:,:]
                self.plasma_P = f['plasma_P'][:,:,:]
                #
                Zions = np.array( f['Zions'][:], dtype=int)
                     
                # check if Zions coincides with that in plasma class
                Zcharge = np.array( [self.plasma_class.Zcharge[ion] for ion in self.plasma_class.ion_species], dtype=int )
                if not np.array_equal(Zions, Zcharge):
                    raise ValueError("Zcharge from file does not match that of plasma class!")
                
                # transpose plasma_N and plasma_T
                self.plasma_N = np.transpose(self.plasma_N, axes=[2,1,0])
                self.plasma_T = np.transpose(self.plasma_T, axes=[2,1,0])
                self.plasma_P = np.transpose(self.plasma_P, axes=[2,1,0])
                
    def create_input_sources_file(self,filename,nt,nrho,tfin):
        # creates 
        
        Zcharge_ions = np.array( [self.plasma_class.Zcharge[ion] for ion in self.plasma_class.ion_species], dtype=int )
        mass_ions    = [self.plasma_class.mass[ion] for ion in self.plasma_class.ion_species]
        nZ   = self.plasma_class.num_ion_species
        
        self.rho_grid_source = np.linspace(0,1,nrho)
        self.t_grid_source = np.linspace(0,tfin,nt)
        
        SE_out = np.zeros((nrho,nt,nZ+1))
        Sn_out = np.zeros((nrho,nt,nZ+1))

        hf = h5py.File(filename, 'w')
                    
        hf.create_dataset('nrho', data=nrho)
        hf.create_dataset('nt', data=nt)
        #
        hf.create_dataset('nion', data=nZ)
        #
        hf.create_dataset('raxis_source', data=self.rho_grid_source)
        hf.create_dataset('taxis_source', data=self.t_grid_source)
        #
        hf.create_dataset('Z_prof', data=Zcharge_ions)
        hf.create_dataset('mass_prof', data=mass_ions)
        #
        hf.create_dataset('S_energy', data=SE_out)
        hf.create_dataset('S_particle', data=Sn_out)
        #
        hf.close()
        
    def add_energy_source_to_input_file(self,filename,which_species,source_type,dVdrho=None,total_power=None,sigma_rho=None,rho_0=None,time_dependent_factor=None,cte_source=None):
        
        from scipy.integrate import quad
        
        try:
            species_id = self.plasma_class.list_of_species.index(which_species) 
        except:
            raise ValueError(f'{which_species} not possible!')
        
        match source_type:
            case 'external_gaussian':
                if((total_power is None) or (sigma_rho is None) or (rho_0 is None) or (dVdrho is None)):
                    print('ERROR: Need to provide total_power [W], dVdrho, sigma_rho and rho_0 for gaussian external source')
                    exit(1) 
                else:
                    # integrand = lambda rho: np.exp(-(rho-rho_0)**2/sigma_rho**2) * dVdrho(rho)
                    # #
                    # cte = total_power / quad(integrand,0,1)[0]
                    # #
                    # source = lambda t,rho: cte * np.exp(-(rho-rho_0)**2/sigma_rho**2)
                    integrand = np.exp(-(self.rho_grid_source - rho_0)**2/sigma_rho**2) * dVdrho(self.rho_grid_source)
                    integrand = integrand.flatten()
                    #
                    cte = total_power / np.trapz(integrand,self.rho_grid_source)
                    #
                    source = lambda t: cte * np.exp(-(self.rho_grid_source - rho_0)**2/sigma_rho**2)
            
            case 'time_dependent_gaussian':
                if((total_power is None) or (sigma_rho is None) or (rho_0 is None) or (dVdrho is None) or (time_dependent_factor is None)):
                    print('ERROR: Need to provide total_power [W], time_Dependent_factor, dVdrho, sigma_rho and rho_0 for gaussian external source')
                    exit(1) 
                else:
                    integrand = np.exp(-(self.rho_grid_source - rho_0)**2/sigma_rho**2) * dVdrho(self.rho_grid_source)
                    integrand = integrand.flatten()
                    #
                    cte = total_power / np.trapz(integrand,self.rho_grid_source)
                    #
                    source = lambda t: time_dependent_factor(t) * cte * np.exp(-(self.rho_grid_source - rho_0)**2/sigma_rho**2)
                    
            case 'constant':
                if(cte_source is None):
                    print('ERROR: cte_source is needed in order to generate a constant source.')
                    exit(0)
                else:
                    source = lambda t: cte_source
                    
            case _:
                print(f'ERROR: Source type {source_type} is NOT possible')
                exit(0)
            
        # add source to which_species without changing the others           
        with h5py.File(filename, 'r+') as f:
            dset = f['S_energy']
            for it,t in enumerate(self.t_grid_source): 
                dset[:,it,species_id] = source(t)#,self.rho_grid_source)
                
    def add_particle_source_to_input_file(self,filename,which_species,source_type,dVdrho=None,injected_particles_per_sec=None,sigma_rho=None,rho_0=None,time_dependent_factor=None,cte_source=None):
        
        from scipy.integrate import quad
        
        try:
            species_id = self.plasma_class.list_of_species.index(which_species) 
        except:
            raise ValueError(f'{which_species} not possible!')
            
        match source_type:
            case 'external_gaussian':
                if((injected_particles_per_sec is None) or (sigma_rho is None) or (rho_0 is None) or (dVdrho is None)):
                    print('ERROR: Need to provide injected_particles_per_sec, dVdrho, sigma_rho and rho_0 for gaussian external source')
                    exit(1) 
                else:
                    # integrand = lambda rho: np.exp(-(rho-rho_0)**2/sigma_rho**2) * dVdrho(rho)
                    # #
                    # cte = injected_particles_per_sec / quad(integrand,0,1)[0]
                    # #
                    # source = lambda t,rho: cte * np.exp(-(rho-rho_0)**2/sigma_rho**2)
                    
                    integrand = np.exp(-(self.rho_grid_source - rho_0)**2/sigma_rho**2) * dVdrho(self.rho_grid_source)
                    integrand = integrand.flatten()
                    #
                    cte = injected_particles_per_sec / np.trapz(integrand,self.rho_grid_source)
                    #
                    source = lambda t: cte * np.exp(-(self.rho_grid_source - rho_0)**2/sigma_rho**2)
            
            case 'time_dependent_gaussian':
                if((injected_particles_per_sec is None) or (sigma_rho is None) or (rho_0 is None) or (dVdrho is None) or (time_dependent_factor is None)):
                    print('ERROR: Need to provide injected_particles_per_sec, time_Dependent_factor, dVdrho, sigma_rho and rho_0 for gaussian external source')
                    exit(1) 
                else:
                    integrand = np.exp(-(self.rho_grid_source - rho_0)**2/sigma_rho**2) * dVdrho(self.rho_grid_source)
                    integrand = integrand.flatten()
                    #
                    cte = injected_particles_per_sec / np.trapz(integrand,self.rho_grid_source)
                    #
                    source = lambda t: time_dependent_factor(t) * cte * np.exp(-(self.rho_grid_source - rho_0)**2/sigma_rho**2)
                    
            case 'constant':
                if(cte_source is None):
                    print('ERROR: cte_source is needed in order to generate a constant source.')
                    exit(0)
                else:
                    source = lambda t: cte_source
                    
            case _:
                print(f'ERROR: Source type {source_type} is NOT possible')
                exit(0)
                
        # add source to which_species without changing the others           
        with h5py.File(filename, 'r+') as f:
            dset = f['S_particle']
            for it,t in enumerate(self.t_grid_source): 
                dset[:,it,species_id] = source(t) #,self.rho_grid_source)
            
                   
# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)