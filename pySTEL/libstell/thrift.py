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

plt.rc('font', size=20)
plt.rcParams['lines.linewidth'] = 4
plt.rcParams['lines.markersize'] = 16
custom_colors = color=['#5FAF30','#1D2258','#A1CDC8','#8b3843','#014817','#cdcd15']
default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=custom_colors+default_colors)

EC = 1.602176634E-19 # Electron charge [C]

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
        
    def read_thrift_folder(self, folder_path):
        """Reads THRIFT HDF5 files inside folder_path

        Files should be named: output_<n>.h5, where <n> is a positive integer.
        The function finds the minimum n and checks that files are consecutively
        numbered with no gaps.

        Parameters
        ----------
        folder_path : str
            Path to folder containing output_#.h5 files.
        """
        import os
        import re

        # Match files like output_123.h5
        pattern = re.compile(r"output_(\d+)\.h5$")
        numbered_files = []

        for file_name in os.listdir(folder_path):
            match = pattern.match(file_name)
            if match:
                number = int(match.group(1))
                full_path = os.path.join(folder_path, file_name)
                if os.path.isfile(full_path):
                    numbered_files.append((number, full_path))

        if not numbered_files:
            raise FileNotFoundError(f"No output_*.h5 files found in {folder_path}")

        # Sort by output number
        numbered_files.sort()
        numbers, files = zip(*numbered_files)

        # Check for sequential numbering
        expected_numbers = list(range(min(numbers), max(numbers) + 1))
        if list(numbers) != expected_numbers:
            missing = sorted(set(expected_numbers) - set(numbers))
            raise FileNotFoundError(f"Missing expected files: {', '.join(f'output_{n}.h5' for n in missing)}")

        # read files
        self.read_thrift(*files)
    
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
				    'THRIFT_COEFF_C','THRIFT_COEFF_CP','THRIFT_COEFF_D','THRIFT_COEFF_DP','THRIFT_DPECRHDV','THRIFT_EPARB','THRIFT_ER','THRIFT_ETAPARA','THRIFT_GNEO',\
                    'THRIFT_I','THRIFT_IBOOT','THRIFT_IECCD','THRIFT_INBCD','THRIFT_IOHMIC','THRIFT_IOTA','THRIFT_IPLASMA','THRIFT_ISOURCE',\
				    'THRIFT_J','THRIFT_JBOOT','THRIFT_JECCD','THRIFT_JNBCD','THRIFT_JOHMIC','THRIFT_JPLASMA','THRIFT_JSOURCE',\
				    'THRIFT_MATLD','THRIFT_MATMD','THRIFT_MATRHS','THRIFT_MATUD','THRIFT_P','THRIFT_PECRH','THRIFT_PHIEDGE','THRIFT_PPRIME',\
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
        self.units_dictionary['THRIFT_QNEO'] = r'[$\text{W}~\text{m}^{-2}$]'
             
    def plot_vars_in_time(self,*vars,time_slice=None,time_array=None,make_plot=True):
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
                    ax.plot(np.sqrt(self.THRIFT_S),plot_var[it,:],label=f't={time:.2f}s'+r', $\beta=$'+f'{self.THRIFT_BETATOT[idx[it]]*100:.2f}%')   
                except:
                    ax.plot(np.sqrt(self.THRIFT_SNOB),plot_var[it,:],label=f't={time:.2f}s')
                ax.set_xlabel('r/a') 
                ax.set_title(var)   
            ax.grid()    
            ax.legend()
            try:
                ax.set_ylabel(self.units_dictionary[var])
            except:
                ax.set_ylabel('')
            # ax.set_yscale('log')
            if(make_plot): plt.show()
            
    def check_var_shape(self,var, nt, nrho):
        
        if isinstance(var, np.ndarray):
            if var.ndim == 2:
                if var.shape != (nt, nrho) and var.shape != (nt, nrho-2):
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
        
        ni = np.atleast_3d( np.array( hf['ni_prof'][:] ) )
        Ti = np.atleast_3d( np.array( hf['ti_prof'][:] ) ) 
        
        nion = np.int64( hf['nion'] )[0]
        
        hf.close()
        
        _, ax_n = plt.subplots(figsize=(11,8))
        _, ax_T = plt.subplots(figsize=(11,8))    
        
        ax_n.plot(raxis,ne[:,selected_indices]/1e20,'.-')
        ax_n.set_xlabel('r/a') 
        ax_n.set_title('ne [1E20 m^-3]')   
        ax_n.grid()   
        
        ax_T.plot(raxis,Te[:,selected_indices]/1e3,'.-')
        ax_T.set_xlabel('r/a') 
        ax_T.set_title('Te [keV]')   
        ax_T.grid()  
        
        plt.show()

        for i in range(nion):
            
            _, ax_n = plt.subplots(figsize=(11,8))
            _, ax_T = plt.subplots(figsize=(11,8))
            
            ax_n.plot(raxis,ni[:,selected_indices,i]/1e20,'.-')
            ax_n.set_xlabel('r/a') 
            ax_n.set_title(f'ni [1E20 m^-3], ion={i+1}')   
            ax_n.grid()
            
            ax_T.plot(raxis,Ti[:,selected_indices,i]/1e3,'.-')
            ax_T.set_xlabel('r/a') 
            ax_T.set_title(f'Ti [keV], ion={i+1}')   
            ax_T.grid()
            ax_T.legend()
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
            if(make_plot):
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
                
                plt.show()
        
        return roa_all, Er_all, JBS_all, Gamma_e_all, num_roots
        
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
        
        for ROA in roa_unique:
        
            _, ax = plt.subplots(figsize=(11,8))
            ax.plot(Er_dict[ROA],Jr_dict[ROA],'.-')
            ax.set_xlabel('Er [V/cm]')
            ax.set_ylabel(r'$\Sigma Z_i\Gamma_i-\Gamma_e$')
            ax.set_title(f'r/a={ROA}')
            ax.grid()
            
            # _, ax = plt.subplots(figsize=(11,8))
            # ax.plot(Er_dict[ROA],Ge_dict[ROA],'.-',label='Gamma_e')
            # for k,_ in enumerate(Zions):
            #     ax.plot(Er_dict[ROA],Gi_dict[k][ROA],'.-',label=f'Gamma_i_{k}')
            # ax.set_xlabel('Er [V/cm]')
            # ax.set_ylabel(r'$\Gamma$')
            # ax.set_yscale('symlog',linthresh=0.1)
            # ax.set_title(f'r/a={ROA}')
            # ax.grid()
            # plt.legend()
            
            plt.show()
    
    def process_wout_folders(self,file_pattern: str, *folders):
        """
        Processes folders containing the wout files outputed by THRIFT, named like:
            file_pattern.N1_N2.nc
        or already processed:
            file_pattern.N1.nc

        For each N1, keeps only the file with the largest N2, renames it to file_pattern.N1.nc,
        and removes the rest. Fully processed folders are handled gracefully.

        Returns a list with all final file Paths across all folders 
        (in the order as given by *folders)
        """
        import re
        from pathlib import Path
        final_files = []

        # Example: file_pattern.001_004.nc
        pattern_full = re.compile(rf"^{re.escape(file_pattern)}\.(\d{{3}})_(\d{{3}})\.(.+)$")
        # Example: file_pattern.001.nc
        pattern_final = re.compile(rf"^{re.escape(file_pattern)}\.(\d{{3}})\.(.+)$")

        for folder in map(Path, folders):
            if not folder.is_dir():
                raise ValueError(f"Folder does not exist: {folder}")

            grouped = {}       # N1 → list of (N2, Path)
            preprocessed = {}  # N1 → Path

            for f in folder.iterdir():
                if not f.is_file():
                    continue

                m_full = pattern_full.match(f.name)
                m_final = pattern_final.match(f.name)

                if m_full:
                    n1, n2, ext = m_full.groups()
                    grouped.setdefault(n1, []).append((int(n2), f))

                elif m_final:
                    n1, ext = m_final.groups()
                    preprocessed[n1] = f

            # If folder already contains only final files
            if grouped == {}:
                final_files.extend(sorted(preprocessed.values()))
                continue

            folder_results = []

            # Process N1 groups
            for n1, files in grouped.items():

                # Already final → keep final, delete full versions
                if n1 in preprocessed:
                    folder_results.append(preprocessed[n1])
                    for _, f in files:
                        f.unlink()
                    continue

                # Select largest N2
                files.sort(key=lambda t: t[0])
                _, best_file = files[-1]

                ext = best_file.suffix  # includes dot

                # New name: file_pattern.N1.ext
                new_name = f"{file_pattern}.{n1}{ext}"
                new_path = best_file.with_name(new_name)

                best_file.rename(new_path)

                # Remove others
                for _, old_file in files[:-1]:
                    old_file.unlink()

                folder_results.append(new_path)

            # Add final files with N1 not appearing in grouped
            for n1, f in preprocessed.items():
                if n1 not in grouped:
                    folder_results.append(f)

            # Sort within folder for deterministic order
            folder_results = sorted(folder_results, key=lambda p: p.name)
            final_files.extend(folder_results)

        return final_files
    
    def get_I_total(self,time=None):
        """ Returns the total current

        This subroutine returns the total current in A.
        The user may provide a timeslice.

        Parameters
        ----------
        time : float (optional)
            Time at which to evaluate profile. (default: last timestamp)
        Returns
        ----------
        I : float
            Total current [A]
        """
        import numpy as np
        from scipy.interpolate import RegularGridInterpolator
        if type(time) == type(None):
            t0 = self.THRIFT_T[-1]
        else:
            t0 = max(time,self.THRIFT_T[-1])
        ftemp = RegularGridInterpolator((self.THRIFT_T,self.THRIFT_S),self.THRIFT_I)
        return float(ftemp([t0,1.0])[0])

    def get_temperature_prof(self,species=0,time=None,ns=64):
        """ Returns a species temperature array

        This subroutine returns the species tempterature.

        Parameters
        ----------
        species : int (optional)
            Species to return. (default: 0 - electrons)
        time : float (optional)
            Time at which to evaluate profile. (default: last timestamp)
        ns : int (optional)
            Number of points to use in evaluation (default: 64)
        Returns
        ----------
        sflx : ndarray
            Array of knots in normalized toroidal flux (s)
        T : ndarray
            Array of values of temperature [eV]
        """
        import numpy as np
        from scipy.interpolate import RegularGridInterpolator
        if type(time) == type(None):
            t0 = self.THRIFT_T[-1]
        else:
            t0 = max(time,self.THRIFT_T[-1])
        sflx = np.linspace(0,1.0,ns)
        tval = np.ones_like(sflx)*t0
        x    = np.vstack((tval,sflx))
        ftemp = RegularGridInterpolator((self.THRIFT_T,self.THRIFT_S),np.squeeze(self.THRIFT_TEMP[:,:,species]))
        return sflx,ftemp(x.T)

    def get_density_prof(self,species=0,time=None,ns=64):
        """ Returns a species density array

        This subroutine returns the species density.

        Parameters
        ----------
        species : int (optional)
            Species to return. (default: 0 - electrons)
        time : float (optional)
            Time at which to evaluate profile. (default: last timestamp)
        ns : int (optional)
            Number of points to use in evaluation (default: 64)
        Returns
        ----------
        sflx : ndarray
            Array of knots in normalized toroidal flux (s)
        N : ndarray
            Array of values of density [m^-3]
        """
        import numpy as np
        from scipy.interpolate import RegularGridInterpolator
        if type(time) == type(None):
            t0 = self.THRIFT_T[-1]
        else:
            t0 = max(time,self.THRIFT_T[-1])
        sflx = np.linspace(0,1.0,ns)
        tval = np.ones_like(sflx)*t0
        x    = np.vstack((tval,sflx))
        ftemp = RegularGridInterpolator((self.THRIFT_T,self.THRIFT_S),np.squeeze(self.THRIFT_DENS[:,:,species]))
        return sflx,ftemp(x.T)

    def get_j_prof(self,time=None,ns=64):
        """ Returns a current profile array

        This subroutine returns the total current density
        array in a 2D array where the first dimension are the points
        in s and the second dimension is the current density in A/m^2.
        The user may provide a timeslice or number of points in an 
        array.

        Parameters
        ----------
        time : float (optional)
            Time at which to evaluate profile. (default: last timestamp)
        ns : int (optional)
            Number of points to use in evaluation (default: 64)
        Returns
        ----------
        sflx : ndarray
            Array of knots in normalized toroidal flux (s)
        jtotal : ndarray
            Array of values of total current [A/m^2]
        """
        import numpy as np
        from scipy.interpolate import RegularGridInterpolator
        if type(time) == type(None):
            t0 = self.THRIFT_T[-1]
        else:
            t0 = max(time,self.THRIFT_T[-1])
        sflx = np.linspace(0,1.0,ns)
        tval = np.ones_like(sflx)*t0
        x    = np.vstack((tval,sflx))
        ftemp = RegularGridInterpolator((self.THRIFT_T,self.THRIFT_S),self.THRIFT_J)
        return sflx,ftemp(x.T)

    def get_Iboot_total(self,time=None):
        """ Returns the total boostrap current

        This subroutine returns the total bootstrap current in kA.
        The user may provide a timeslice.

        Parameters
        ----------
        time : float (optional)
            Time at which to evaluate profile. (default: last timestamp)
        Returns
        ----------
        I : float
            Total bootstrap current [A]
        """
        import numpy as np
        from scipy.interpolate import RegularGridInterpolator
        if type(time) == type(None):
            t0 = self.THRIFT_T[-1]
        else:
            t0 = max(time,self.THRIFT_T[-1])
        ftemp = RegularGridInterpolator((self.THRIFT_T,self.THRIFT_S),self.THRIFT_IBOOT)
        return float(ftemp([t0,1.0])[0])

    def get_jboot_prof(self,time=None,ns=64):
        """ Returns the bootstrap current density

        This subroutine returns the bootstrap current density [A/m^2].

        Parameters
        ----------
        time : float (optional)
            Time at which to evaluate profile. (default: last timestamp)
        ns : int (optional)
            Number of points to use in evaluation (default: 64)
        Returns
        ----------
        sflx : ndarray
            Array of knots in normalized toroidal flux (s)
        jboot : ndarray
            Array of values of bootstrap current [A/m^2]
        """
        import numpy as np
        from scipy.interpolate import RegularGridInterpolator
        if type(time) == type(None):
            t0 = self.THRIFT_T[-1]
        else:
            t0 = max(time,self.THRIFT_T[-1])
        sflx = np.linspace(0,1.0,ns)
        tval = np.ones_like(sflx)*t0
        x    = np.vstack((tval,sflx))
        ftemp = RegularGridInterpolator((self.THRIFT_T,self.THRIFT_S),self.THRIFT_JBOOT)
        return sflx,ftemp(x.T)

    def get_jsource_prof(self,time=None,ns=64):
        """ Returns the source current density

        This subroutine returns the source current density [A/m^2].

        Parameters
        ----------
        time : float (optional)
            Time at which to evaluate profile. (default: last timestamp)
        ns : int (optional)
            Number of points to use in evaluation (default: 64)
        Returns
        ----------
        sflx : ndarray
            Array of knots in normalized toroidal flux (s)
        jboot : ndarray
            Array of values of source current [A/m^2]
        """
        import numpy as np
        from scipy.interpolate import RegularGridInterpolator
        if type(time) == type(None):
            t0 = self.THRIFT_T[-1]
        else:
            t0 = max(time,self.THRIFT_T[-1])
        sflx = np.linspace(0,1.0,ns)
        tval = np.ones_like(sflx)*t0
        x    = np.vstack((tval,sflx))
        ftemp = RegularGridInterpolator((self.THRIFT_T,self.THRIFT_S),self.THRIFT_JSOURCE)
        return sflx,ftemp(x.T)

    def get_iota_prof(self,time=None,ns=64):
        """ Returns the rotational transform array

        This subroutine returns the rotational transform (iota)
        array in a 2D array where the first dimension are the points
        in s and the second dimension is the rotational transform.
        The user may provide a timeslice or number of points in an 
        array.

        Parameters
        ----------
        time : float (optional)
            Time at which to evaluate profile. (default: last timestamp)
        ns : int (optional)
            Number of points to use in evaluation (default: 64)
        Returns
        ----------
        sflx : ndarray
            Array of knots in normalized toroidal flux (s)
        iota : ndarray
            Array of values of rotational transform
        """
        import numpy as np
        from scipy.interpolate import RegularGridInterpolator
        if type(time) == type(None):
            t0 = self.THRIFT_T[-1]
        else:
            t0 = max(time,self.THRIFT_T[-1])
        sflx = np.linspace(0,1.0,ns)
        tval = np.ones_like(sflx)*t0
        x    = np.vstack((tval,sflx))
        ftemp = RegularGridInterpolator((self.THRIFT_T,self.THRIFT_S),self.THRIFT_IOTA)
        return sflx,ftemp(x.T)

    def get_pot_prof(self,time=None,ns=64):
        """ Returns the electrostatic potential array

        This subroutine returns the electrostatic potential
        array in a 2D array where the first dimension are the points
        in s and the second dimension is the electrostatic potential 
        in V. The user may provide a timeslice or number of points 
        in an array.

        Parameters
        ----------
        time : float (optional)
            Time at which to evaluate profile. (default: last timestamp)
        ns : int (optional)
            Number of points to use in evaluation (default: 64)
        Returns
        ----------
        sflx : ndarray
            Array of knots in normalized toroidal flux (s)
        phi : ndarray
            Array of values of electrostatic potential [V]
        """
        import numpy as np
        from scipy.interpolate import RegularGridInterpolator, interpn
        from scipy.integrate import cumulative_trapezoid
        if type(time) == type(None):
            t0 = self.THRIFT_T[-1]
        else:
            t0 = max(time,self.THRIFT_T[-1])
        s,A = self.get_Aminor(time=time,ns=ns)
        s,Er = self.get_er_prof(time=time,ns=ns)
        # Compute the estatic potential
        phi = -cumulative_trapezoid(Er,A,initial=0)
        return s,phi

    def get_er_prof(self,time=None,ns=64):
        """ Returns the radial electric field array

        This subroutine returns the radial electric field
        array in a 2D array where the first dimension are the points
        in s and the second dimension is Er in V/m.
        The user may provide a timeslice or number of points in an 
        array.

        Parameters
        ----------
        time : float (optional)
            Time at which to evaluate profile. (default: last timestamp)
        ns : int (optional)
            Number of points to use in evaluation (default: 64)
        Returns
        ----------
        sflx : ndarray
            Array of knots in normalized toroidal flux (s)
        phi : ndarray
            Array of values of radial electric field [V/m]
        """
        import numpy as np
        from scipy.interpolate import RegularGridInterpolator
        if type(time) == type(None):
            t0 = self.THRIFT_T[-1]
        else:
            t0 = max(time,self.THRIFT_T[-1])
        sflx = np.linspace(0,1.0,ns)
        tval = np.ones_like(sflx)*t0
        x    = np.vstack((tval,sflx))
        ftemp = RegularGridInterpolator((self.THRIFT_T,self.THRIFT_S),self.THRIFT_ER)
        return sflx,ftemp(x.T)

    def get_Aminor(self,time=None,ns=64):
        """ Returns the minor radius array

        This subroutine returns the minor radius
        array in a 2D array where the first dimension are the points
        in s and the second dimension is the minor radius in m.
        The user may provide a timeslice or number of points in an 
        array.

        Parameters
        ----------
        time : float (optional)
            Time at which to evaluate profile. (default: last timestamp)
        ns : int (optional)
            Number of points to use in evaluation (default: 64)
        Returns
        ----------
        sflx : ndarray
            Array of knots in normalized toroidal flux (s)
        Aminor : ndarray
            Array of values of minor radius [m]
        """
        import numpy as np
        from scipy.interpolate import RegularGridInterpolator
        if type(time) == type(None):
            t0 = self.THRIFT_T[-1]
        else:
            t0 = max(time,self.THRIFT_T[-1])
        sflx = np.linspace(0,1.0,ns)
        tval = np.ones_like(sflx)*t0
        x    = np.vstack((tval,sflx))
        ftemp = RegularGridInterpolator((self.THRIFT_T,self.THRIFT_S),self.THRIFT_AMINOR)
        return sflx,ftemp(x.T)
    
    def create_dkes_results_file_from_DKES_coeffs(self,DKES_coeffs_file,k):
        """This function creates a results.surface_k file similar to the results file outputed by DKES
        The DKES_coeffs_file is a file outputed by THRIFT with the DKES coefficients at all surfaces
        """
        # Containers for rows with the requested k
        selected_rows = []

        with open(DKES_coeffs_file, "r") as f:
            for line in f:
                line = line.strip()

                # Skip empty lines or header
                if not line or line.lower().startswith("dkes_k"):
                    continue

                parts = line.split()
                if len(parts) != 6:
                    continue  # defensive: ignore malformed lines

                k_val = int(parts[0])

                if k_val == k:
                    # Parse values
                    _, Er_v, nu_v, D11, D31, D33 = parts
                    selected_rows.append(
                        (
                            float(nu_v),   # cmul
                            float(Er_v),   # efield
                            float(D11),
                            float(D31),
                            float(D33),
                        )
                    )

        # Check that k was found
        if not selected_rows:
            raise ValueError("The given k is not in the input file!")

        # Output file
        outname = f"results.surface_{k}"

        with open(outname, "w") as fout:
            # Header
            fout.write("*\n")
            fout.write(
                "cmul\tefield\tweov\twtov\t"
                "L11m\tL11p\tL31m\tL31p\tL33m\tL33p\n"
            )

            # Write data preserving original order
            for nu_v, Er_v, D11, D31, D33 in selected_rows:
                weov = 0.0
                wtov = 0.0

                fout.write(
                    f"{nu_v:.8e}\t{Er_v:.8e}\t{weov:.1f}\t{wtov:.1f}\t"
                    f"{D11:.8e}\t{D11:.8e}\t"
                    f"{D31:.8e}\t{D31:.8e}\t"
                    f"{D33:.8e}\t{D33:.8e}\n"
                )

        return outname
    
    def solve_time_dependent_current_equation(self,L_inductance,edge_factor=None,add_neglected_terms=False):
        """ Solves the time dependent current diffusion equation. Profiles are time-dependent
            
            Equation is written in the conservative form:
            dI/dt = (S11/phia^2) * d/ds[ D(s,t)*dI/ds + P(s,t)*I ] + F(s)
            where:
            D(s) = etapar*V' * <B2>/mu0
            P(s) = etapar*V' * p'
            F(s) = -S11/phia^2 * d/ds[etapar*V' * phia/mu0 * <JNI*B>]  +  phia/mu0 * [dS12/dt + iota*dS11/dt]
            
            We solve this time-dependent PDE using an implicit backward Euler time scheme:
            (I_new - I_old)/dt = L*I_new + F
        """
        from scipy.sparse import diags, identity
        from scipy.sparse.linalg import spsolve
        mu0 = 4*np.pi*1e-7
        
        # Define time array
        t_solver = self.THRIFT_T
        Nt = len(t_solver)
        #
        Ns = len(self.THRIFT_S)
        ds = self.THRIFT_S[1] - self.THRIFT_S[0]
        
        J_SOURCE = self.THRIFT_JSOURCE.copy() # copy to avoid modifying original array
        if(edge_factor is not None):
            J_SOURCE[:,-1] *= edge_factor 
        
        phia = self.THRIFT_PHIEDGE
        #
        aux1 = self.THRIFT_S11/phia[:,np.newaxis]**2
        aux2 = self.THRIFT_ETAPARA*self.THRIFT_VP
        aux3 = aux2 * (phia[:,np.newaxis]/mu0) * J_SOURCE *self.THRIFT_BAV
        #
        D = aux2 * self.THRIFT_BSQAV/mu0
        P = aux2 * self.THRIFT_PPRIME
        F = -aux1 * np.gradient(aux3,self.THRIFT_S,axis=1)
        #
        if(add_neglected_terms):
            iota = self.THRIFT_IOTA
            S11 = self.THRIFT_S11
            S12 = self.THRIFT_S12
            dS11dt = np.gradient(S11,t_solver,axis=0)
            dS12dt = np.gradient(S12,t_solver,axis=0)
            #
            F = F + (phia[:,np.newaxis]/mu0) * (dS12dt + iota*dS11dt)
        
        I_solution = np.zeros((Nt,Ns)) # initial condition is I(t=0,s) = 0.0
        #
        for it in range(1,Nt):
            
            # Assemble L operator
            Dminus = 0.5 * (D[it,0:-2] + D[it,1:-1])
            Dplus  = 0.5 * (D[it,1:-1] + D[it,2:])
            #
            Pminus = 0.5 * (P[it,0:-2] + P[it,1:-1])
            Pplus  = 0.5 * (P[it,1:-1] + P[it,2:])
            #
            main_diag = np.zeros(Ns)
            upper_diag = np.zeros(Ns-1)
            lower_diag = np.zeros(Ns-1)
            #
            main_diag[1:-1] = -Dplus/ds - Pplus/2 - Dminus/ds - Pminus/2
            upper_diag[1:] = Dplus/ds + Pplus/2 
            lower_diag[0:-1] = Dminus/ds + Pminus/2
            #
            L = diags([lower_diag/ds, main_diag/ds, upper_diag/ds], offsets=[-1, 0, 1], format="csr")  
            
            # Multiply L operator by S11/phia^2 term
            G = diags(aux1[it,:], 0, format="csr")
            L = G @ L
            
            # Create LHS matrix
            dt = t_solver[it] - t_solver[it-1]
            Id = identity(Ns, format="csr")
            LHS = Id - dt*L
            
            # Create RHS matrix
            RHS_last = -phia[it]*L_inductance / (self.THRIFT_VP[it,-1]*self.THRIFT_ETAPARA[it,-1])
            # 
            diag = np.ones(Ns)
            diag[-1] = RHS_last
            #
            RHS = diags(diag, offsets=0, format="csr")
            
            # Source Vector with boundary conditions
            SOURCE = dt*F[it,:]
            SOURCE[0] = 0.0
            SOURCE[-1] = -J_SOURCE[it,-1]*self.THRIFT_BAV[it,-1]*dt
            
            # Apply edge BC on LHS
            LHS = LHS.tolil()  # Convert to LIL for easy row modification
            VpEtapar = self.THRIFT_VP[it,-1]*self.THRIFT_ETAPARA[it,-1]
            B2 = self.THRIFT_BSQAV[it,-1]
            pp = self.THRIFT_PPRIME[it,-1]
            LHS[-1,:] = 0.0 
            LHS[-1,-1] = -phia[it]*L_inductance/(VpEtapar) - 1.5*B2*dt/(phia[it]*ds) - mu0*pp*dt/phia[it]
            LHS[-1,-2] = 2.0*B2*dt / (phia[it]*ds)
            LHS[-1,-3] = -0.5*B2*dt / (phia[it]*ds)
            LHS = LHS.tocsr()
            
            # Solve    
            I_solution[it,:] = spsolve(LHS, RHS.dot(I_solution[it-1,:]) + SOURCE)
        
        return t_solver, I_solution
    
    def create_plasma_profiles_file_from_transport_solver_jolib(self,joblib_file,output_filename='plasma_profiles.h5'):
        """ This function creates a plasma profiles files, which is read by THRIFT, 
        from a joblib output file of transport solver solver
        """
        import joblib
        import h5py
        from libstell.plasma import PLASMA
        
        solver = joblib.load(joblib_file)
        
        time_array = solver.time
        nt = len(time_array)
        raxis_prof = solver.rho_grid
        nrho = len(raxis_prof)
        
        ne = solver.N['electrons'].T
        Te = solver.T['electrons'].T
        
        ions = [s for s in solver.list_of_species if s != 'electrons']
        num_ions = len(ions)
        
        plasma = PLASMA(list_of_species=solver.list_of_species)
        Zcharge_ions = np.array( [plasma.Zcharge[ion] for ion in plasma.ion_species], dtype=float )
        mass_ions    = [plasma.mass[ion] for ion in plasma.ion_species]
        
        ni = np.stack([solver.N[ion].T for ion in ions], axis=-1)
        Ti = np.stack([solver.T[ion].T for ion in ions], axis=-1)  

        hf = h5py.File(output_filename, 'w')
        #
        hf.create_dataset('nrho', data=nrho)
        hf.create_dataset('nt', data=nt)
        hf.create_dataset('nion', data=num_ions)
        hf.create_dataset('raxis_prof', data=raxis_prof)
        hf.create_dataset('taxis_prof', data=time_array)
        hf.create_dataset('Z_prof', data=Zcharge_ions)
        hf.create_dataset('mass_prof', data=mass_ions)
        hf.create_dataset('ne_prof', data=ne)
        hf.create_dataset('te_prof', data=Te)
        hf.create_dataset('ni_prof', data=ni)
        hf.create_dataset('ti_prof', data=Ti)
        #
        hf.close()
    
    def plot_thrift_vars_subiterations(self,folder, plot_var):
        """
        Plots arrays from files named:

            thrift_vars*.xxx_yyy

        where:
            xxx = iteration number
            yyy = subiteration number

        Parameters
        ----------
        folder : str
            Path to folder containing the files.

        plot_var : str
            One of:
                'THRIFT_J'
                'THRIFT_JBOOT'
                'THRIFT_JECCD'
        """
        import re
        import os
        import glob
        allowed_vars = ['THRIFT_J', 'THRIFT_JBOOT', 'THRIFT_JECCD']

        if plot_var not in allowed_vars:
            raise ValueError(
                f"plot_var must be one of {allowed_vars}"
            )

        # Column mapping
        col_map = {
            'THRIFT_J': 1,
            'THRIFT_JBOOT': 2,
            'THRIFT_JECCD': 3
        }

        # Regex for extracting iteration/subiteration
        pattern = re.compile(r'.*\.(\d+)_(\d+)$')

        # Find files
        files = glob.glob(os.path.join(folder, 'thrift_vars*.*_*'))

        if len(files) == 0:
            raise ValueError(f'No matching files found in {folder}')

        # Organize files by iteration
        iterations = {}

        for f in files:

            basename = os.path.basename(f)

            match = pattern.match(basename)

            if match is None:
                continue

            iteration = int(match.group(1))
            subiteration = int(match.group(2))

            if iteration not in iterations:
                iterations[iteration] = []

            iterations[iteration].append((subiteration, f))

        # Sort iterations
        sorted_iterations = sorted(iterations.keys())

        # Create one figure per iteration
        for iteration in sorted_iterations:

            _, ax = plt.subplots(figsize=(11,8))
            _, ax2 = plt.subplots(figsize=(11,8))

            # Sort subiterations
            subiter_files = sorted(iterations[iteration],
                                key=lambda x: x[0])

            data_all = []
            for subiteration, filepath in subiter_files:

                # Load data
                data = np.loadtxt(filepath, skiprows=1)

                roa = data[:,0]
                y = data[:, col_map[plot_var]]

                ax.plot(roa,y,label=f'{subiteration:03d}')
                data_all.append(y)
                
            ax.set_xlabel('r/a')
            ax.set_ylabel(plot_var)
            ax.set_title(f'{plot_var}, it={iteration}')
            ax.legend()
            ax.grid(True)
            #
            data_all = np.array(data_all)
            ax2.plot(data_all,'.-',markersize=10,linewidth=2.5,label=roa)
            # ax2.plot(data_all[:,0:-1:20],'.-',markersize=10,linewidth=2.5,label=roa[0:-1:20])
            ax2.set_xlabel('nsubiter')
            ax2.set_ylabel(plot_var)
            ax2.set_title(f'{plot_var}, it={iteration}')
            ax2.legend()
            ax2.grid(True)
            
        plt.show()
    
# THRIFT Class
class THRIFT_plasma_solver():
    """" Class for working with plasma solver implemented in THRIFT
    
    """
    
    def __init__(self, plasma=None, list_of_species=None):
        
        from libstell.plasma import PLASMA
        
        # should give plasma OR list_of_species. If both, plasma prevails
        if plasma is None and list_of_species is None:
            print('Could not create class: Need to provide a plasma class or a list with name of species')
            exit(1)
        elif plasma is None and list_of_species is not None:
            self.plasma_class = PLASMA(list_of_species=list_of_species)
        else:
            self.plasma_class = plasma
            
        self.list_of_species = list_of_species
        
    def read_thrift_plasma_solver_folder(self, folder_path):
        """Reads THRIFT PLASMA_SOLVER HDF5 files inside folder_path

        Files should be named: plasma_solver_<n>.h5, where <n> is a positive integer.
        The function finds the minimum n and checks that files are consecutively
        numbered with no gaps.

        Parameters
        ----------
        folder_path : str
            Path to folder containing plasma_solver_#.h5 files.
        """
        import os
        import re

        # Match files like plasma_solver_123.h5
        pattern = re.compile(r"plasma_solver_(\d+)\.h5$")
        numbered_files = []

        for file_name in os.listdir(folder_path):
            match = pattern.match(file_name)
            if match:
                number = int(match.group(1))
                full_path = os.path.join(folder_path, file_name)
                if os.path.isfile(full_path):
                    numbered_files.append((number, full_path))

        if not numbered_files:
            raise FileNotFoundError(f"No plasma_solver_*.h5 files found in {folder_path}")

        # Sort by output number
        numbered_files.sort()
        numbers, files = zip(*numbered_files)

        # Check for sequential numbering
        expected_numbers = list(range(min(numbers), max(numbers) + 1))
        if list(numbers) != expected_numbers:
            missing = sorted(set(expected_numbers) - set(numbers))
            raise FileNotFoundError(f"Missing expected files: {', '.join(f'plasma_solver_{n}.h5' for n in missing)}")

        # read files
        self.read_thrift_plasma_solver(*files)  
                
    def read_thrift_plasma_solver(self,*files):
        """Reads plasma_solver THRIFT HDF5 files

		Parameters
		----------
		files : str
		Path to HDF5 files.
		"""
        import h5py
        
        ############### CHECK TIME ORDER #######################
        time = []
        for file in files:
            with h5py.File(file,'r') as f:
                time.append( f['time_plasma_grid'][:] )
        time = np.concatenate(time)
        #check ordering
        if(not np.all(np.diff(time) >= -1e-10) ):
            print('ERROR: plasma_solver files are not in the correct order...')
            for diff in np.diff(time): 
                if diff<0: print(diff)
            exit(0)
        else:
            self.time_grid = time
            self.Nt = len(time)
            
        ######## CHECK Zions IS THE SAME IN ALL FILES and that it coincides with that in plasma class  ##########
        Z_plasma_class = np.array( [self.plasma_class.Zcharge[ion] for ion in self.plasma_class.ion_species], dtype=int )
        for file in files:
            with h5py.File(file,'r') as f:
                Zions = np.array( f['Zions'][:], dtype=int)
                if not np.array_equal(Zions, Z_plasma_class):
                    raise ValueError("Zcharge from file does not match that of plasma class!")
                
        ############ CHECK RHO_GRID IS THE SAME IN ALL FILES ################
        with h5py.File(files[0],'r') as f:
                self.rho_grid  = f['rho_plasma_grid'][:]
                self.Nr = f['Nr_plasma_grid']
        for file in files:
            with h5py.File(file,'r') as f:
                rho_grid  = f['rho_plasma_grid'][:]
                if (not np.array_equal(rho_grid, self.rho_grid)):
                    raise ValueError("rho_grid NOT EQUAL among all the files...")
                    
        ################ CONCATENATE DATA ##################################
        for file in files:
            with h5py.File(file,'r') as f:
                for temp in ['r_plasma_grid','plasma_N','plasma_T','N_fast_alphas','Dn_NEO','cn_NEO','Dp_NEO',\
                             'cp_NEO','G_NEO_complet','Q_NEO_complet','Dp_total','cp_total','Dn_total','cn_total',\
                             'S_radiated_power','S_alpha_power','S_energy_ext','S_particle_ext','dVdr','plasma_Er',\
                             'iota2o3','aminor','Baxis']:
                    data = np.array(f[temp])
                    # Check if the attribute exists; if not, initialize it
                    if not hasattr(self, temp):
                        setattr(self, temp, data)
                    else:
                        # Concatenate the new data to the existing attribute along the last axis
                        existing_data = getattr(self, temp)
                        if data.ndim == 1:
                            setattr(self, temp, np.concatenate((existing_data, data), axis=0))       
                        else:
                            setattr(self, temp, np.concatenate((existing_data, data), axis=1))                    
                    
        ##################### TRANSPOSE DATA #################################
        #### 2D arrays should be [time,rho]
        #### 3D arrays should be [species,time,rho]
        for name, value in vars(self).items():
            if isinstance(value, np.ndarray):
                if value.ndim == 2:
                    # Transpose 2D array
                    setattr(self, name, value.T)
                elif value.ndim == 3:
                    # Transpose 3D array
                    setattr(self, name, np.transpose(value, axes=[2, 1, 0]))

        # Rename r-grid
        self.r_grid = self.r_plasma_grid
        delattr(self,'r_plasma_grid')
        
        ##################  DELETE REPEATED TIMES ##########################       
        # Repeated times occur when sequent plasma_solver files are read
        mask = np.concatenate(([True],np.abs(np.diff(time)) > 1E-6 ))
        self.time_grid = self.time_grid[mask]
        
        # Apply mask to all 2D and 3D arrays
        for name, value in vars(self).items():
            if isinstance(value, np.ndarray):
                # 2D array
                if value.ndim == 2:
                    setattr(self,name,value[mask,:])
                # 3D array
                elif value.ndim == 3:
                    setattr(self,name,value[:,mask,:])        
        
    def create_input_sources_file(self,filename,nrho,tgrid):
        # creates 
        
        Zcharge_ions = np.array( [self.plasma_class.Zcharge[ion] for ion in self.plasma_class.ion_species], dtype=int )
        mass_ions    = [self.plasma_class.mass[ion] for ion in self.plasma_class.ion_species]
        nZ   = self.plasma_class.num_ion_species
        
        self.rho_grid_source = np.linspace(0,1,nrho)
        self.t_grid_source = tgrid
        nt = len(tgrid)
        
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
                dset[:,it,species_id] = source(t)
                
    def add_particle_source_to_input_file(self,filename,which_species,source_type,dVdrho=None,injected_particles_per_sec=None,
                                          sigma_rho=None,rho_0=None,time_dependent_factor=None,cte_source=None,source_2D=None):
        
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
                    
            case 'source_2D':
                if(source_2D is None):
                    raise ValueError('Need to provide the 2D arrays source_2D')
                # Check if source_2D has the expected dimensions: nt x nrho
                nt = len(self.t_grid_source)
                nrho = len(self.rho_grid_source)
                if(source_2D.shape != (nt,nrho)):
                    raise ValueError(f'source_2D array does not have the expected ({nt},{nrho}) shape!')
                # Save directly in file and leave
                with h5py.File(filename, 'r+') as f:
                    dset = f['S_particle']
                    dset[:,:,species_id] = source_2D.T
                return
                            
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
                dset[:,it,species_id] = source(t)
                
    def convert_to_joblib(self,dt_save=0.1,filename='thrift_transport_simul'):
        """ This function creates a joblib file with the transport simulation data
        We can the use the same post-processing tools we use to analyse transport simulations
        performed by pySTEL class plasma_solver
        """
        from types import SimpleNamespace
        from pathlib import Path
        import joblib
        from collections import defaultdict
        
        # check if extension of filename is .joblib; if not, add
        output_filename = str(Path(filename).with_suffix(".joblib"))
        
        saved_class = SimpleNamespace()
        saved_class.rho_grid = self.rho_grid
        saved_class.list_of_species = self.list_of_species
        
        # simulation dt
        dt = np.diff(self.time_grid)[-1]
        
        # saving frequency 
        freq = max(1, round(dt_save / dt))
        sl = slice(0, None, freq)  # defines the slice once

        saved_class.time = self.time_grid[sl]
        saved_class.Nt = len(self.time_grid[sl])
        
        saved_class.r_grid = self.r_grid[sl,:]
        saved_class.dVdr = self.dVdr[sl,:]
        
        saved_class.Er = self.plasma_Er[sl,:]
        
        for attr1,attr2 in zip(\
            ('N','T','Dp','cp','Dn','cn','Dn_NEO','Dp_NEO','cp_NEO','cn_NEO','Q_NEO_complet'),\
            ('plasma_N','plasma_T','Dp_total','cp_total','Dn_total','cn_total','Dn_NEO','Dp_NEO','cp_NEO','cn_NEO','Q_NEO_complet')):
            setattr(saved_class, attr1, {})
            for ispecies,species in enumerate(self.list_of_species):
                getattr(saved_class, attr1)[species] = getattr(self, attr2)[ispecies,sl,:]
        
        saved_class.N['alphas_fast'] = self.N_fast_alphas[sl,:]
        saved_class.aminor = self.aminor[sl]
        saved_class.iota2o3 = self.iota2o3[sl]
        saved_class.Baxis = self.Baxis[sl]

        saved_class.explicit_energy_sources   = defaultdict(dict)
        saved_class.explicit_particle_sources = defaultdict(dict)
        saved_class.Q_total = defaultdict(dict)
        saved_class.Gamma_total = defaultdict(dict)
        saved_class.Q_NEO = defaultdict(dict)
        saved_class.Gamma_NEO = defaultdict(dict)
        saved_class.Q_turb = defaultdict(dict)
        saved_class.Gamma_turb = defaultdict(dict)

        saved_class.explicit_energy_sources['electrons']['Bremsstrahlung'] = -self.S_radiated_power[sl,:]
        
        for ispecies,species in enumerate(self.list_of_species):
            saved_class.explicit_energy_sources[species]['time_dependent_gaussian'] = self.S_energy_ext[ispecies,sl,:]
            saved_class.explicit_particle_sources[species]['time_dependent_gaussian'] = self.S_particle_ext[ispecies,sl,:]
            saved_class.explicit_energy_sources[species]['alpha_heating'] = self.S_alpha_power[ispecies,sl,:]
            
            # Reconstruct Fluxes and save
            r_grid   = self.r_grid[sl,:]
            #
            p_r = self.plasma_N[ispecies,sl,:]*self.plasma_T[ispecies,sl,:]*EC
            dpdr = akima_derivative(r_grid,p_r,axis=1)
            #
            n_r = self.plasma_N[ispecies,sl,:]
            dndr = akima_derivative(r_grid,n_r,axis=1)
            #
            saved_class.Q_total[species] = -self.Dp_total[ispecies,sl,:]*dpdr + self.cp_total[ispecies,sl,:]*p_r
            saved_class.Gamma_total[species] = -self.Dn_total[ispecies,sl,:]*dndr + self.cn_total[ispecies,sl,:]*n_r
            #
            saved_class.Q_NEO[species] = -self.Dp_NEO[ispecies,sl,:]*dpdr + self.cp_NEO[ispecies,sl,:]*p_r
            saved_class.Gamma_NEO[species] = -self.Dn_NEO[ispecies,sl,:]*dndr + self.cn_NEO[ispecies,sl,:]*n_r
            #
            saved_class.Q_turb[species] = saved_class.Q_total[species] - saved_class.Q_NEO[species]
            saved_class.Gamma_turb[species] = saved_class.Gamma_total[species] - saved_class.Gamma_NEO[species]

        joblib.dump(saved_class, output_filename)
        
def _akima_derivative_1d(x, y):
    """
    1D Akima derivative routine.
    """
    n = x.size
    dy = np.zeros_like(y)
    if n < 2:
        raise ValueError("Need at least 2 points")

    # First divided differences
    m = (y[1:] - y[:-1]) / (x[1:] - x[:-1])

    if n == 2:
        dy[0] = m[0]
        dy[1] = m[0]
        return dy

    # Boundary slopes
    cxp, cxpp = m[0], m[1]
    cxm, cxmm = m[-1], m[-2]

    dy[0] = 1.5 * cxp - 0.5 * cxpp
    dy[-1] = 1.5 * cxm - 0.5 * cxmm

    # Ghost slopes
    cxtrap0 = 2.0 * dy[0] - cxp
    cxtrap1 = 2.0 * dy[-1] - cxm

    # Interior points
    for i in range(1, n - 1):
        if i == 1:
            cxmm = cxtrap0
        else:
            cxmm = (y[i - 1] - y[i - 2]) / (x[i - 1] - x[i - 2])

        cxm = (y[i] - y[i - 1]) / (x[i] - x[i - 1])
        cxp = (y[i + 1] - y[i]) / (x[i + 1] - x[i])

        if i == n - 2:
            cxpp = cxtrap1
        else:
            cxpp = (y[i + 2] - y[i + 1]) / (x[i + 2] - x[i + 1])

        w1 = abs(cxp - cxpp)
        w2 = abs(cxm - cxmm)

        if (w1 + w2) == 0.0:
            dy[i] = 0.5 * (cxm + cxp)
        else:
            dy[i] = (w1 * cxm + w2 * cxp) / (w1 + w2)

    return dy


def akima_derivative(x, y, axis=None):
    """
    General multidimensional Akima derivative.

    Parameters
    ----------
    x : ndarray, 1D or same shape as y
        Grid values along the differentiation axis.
    y : ndarray
        Values to differentiate.
    axis : int
        Axis along which the derivative is taken.

    Returns
    -------
    dy : ndarray
        Derivative of y along the chosen axis.
    """
    x = np.asarray(x)
    y = np.asarray(y)

    if axis is None:
        axis = y.ndim - 1
    axis = np.core.numeric.normalize_axis_index(axis, y.ndim)

    # Move the target axis to the last dimension
    y_m = np.moveaxis(y, axis, -1)

    # Broadcast x to match y_m.shape
    if x.ndim == 1:
        # x is 1D: must match the size of the differentiation axis
        if x.size != y_m.shape[-1]:
            raise ValueError(
                f"1D x has length {x.size}, but y has size {y_m.shape[-1]} along axis {axis}."
            )
        # Expand x to the same shape as y_m
        # This mirrors numpy.trapz behavior
        shape = (1,) * (y_m.ndim - 1) + (x.size,)
        x_m = np.broadcast_to(x.reshape(shape), y_m.shape)
    else:
        # x must have the same full shape as y
        if x.shape != y.shape:
            raise ValueError("If x is not 1D, it must have the same shape as y.")
        x_m = np.moveaxis(x, axis, -1)

    # Flatten all dimensions except the last
    leading_shape = y_m.shape[:-1]
    N = y_m.shape[-1]

    y_flat = y_m.reshape(-1, N)
    x_flat = x_m.reshape(-1, N)

    dy_flat = np.empty_like(y_flat)

    # Apply the 1D Akima routine along the last axis for each slice
    for i in range(y_flat.shape[0]):
        dy_flat[i] = _akima_derivative_1d(x_flat[i], y_flat[i])

    # Restore multidimensional shape
    dy_m = dy_flat.reshape(y_m.shape)

    # Move axis back to original position
    dy = np.moveaxis(dy_m, -1, axis)

    return dy

def animate_time_series(
    time_axis,
    panels,
    legends=None,
    ylabels=None,
    xlabel=None,
    title=None,
    **FuncAnimation_kwargs):
    """
    Time-series animation with any number of panels.

    Parameters
    ----------
    time_axis : 1D array
        The common x-axis for all panels.
    panels : list
        List of panels. Each panel may be:
            * a single 1D array
            * a list/tuple of arrays (multiple signals on same panel)
    legends : list (same length as panels)
        Each element is:
            * None (auto labels)
            * a single string (for single-signal panel)
            * a list of strings (for multi-signal panel)
    ylabels : list (same length as panels)
        Y-axis labels for each panel (None allowed).
    xlabel : str
        Label for the bottom x-axis.

    Returns
    -------
    anim : matplotlib FuncAnimation
    """
    from matplotlib.animation import FuncAnimation
    
    time_axis = np.asarray(time_axis)
    if time_axis.ndim != 1:
        raise ValueError("time_axis must be a 1D array-like.")

    # ---- Normalize panel structure ----
    normalized_panels = []
    for p in panels:
        if p is None:
            normalized_panels.append([])
            continue
        if isinstance(p, (list, tuple)):
            arrs = [np.asarray(a) for a in p]
        else:
            arrs = [np.asarray(p)]
        for a in arrs:
            if a.shape != time_axis.shape:
                raise ValueError("All signals must match time_axis shape.")
        normalized_panels.append(arrs)

    # keep only non-empty panels
    non_empty_indices = [i for i, p in enumerate(normalized_panels) if len(p) > 0]
    panels_non_empty = [normalized_panels[i] for i in non_empty_indices]
    Npanels = len(panels_non_empty)
    if Npanels == 0:
        raise ValueError("No non-empty panels provided.")

    # ---- Normalize legends ----
    if legends is None:
        legends = [None] * len(normalized_panels)
    if len(legends) != len(normalized_panels):
        raise ValueError("legends must have same length as panels.")

    legends_non_empty = []
    for i in non_empty_indices:
        lg = legends[i]
        sigs = normalized_panels[i]
        n = len(sigs)

        if lg is None:
            legends_non_empty.append([f"Signal {j+1}" for j in range(n)])
        elif isinstance(lg, str):
            if n == 1:
                legends_non_empty.append([lg])
            else:
                raise ValueError("Multi-signal panel requires a list of labels.")
        else:
            lg_list = list(lg)
            if len(lg_list) != n:
                raise ValueError("Legend length mismatch.")
            legends_non_empty.append(lg_list)

    # ---- Normalize ylabels ----
    if ylabels is None:
        ylabels = [None] * len(normalized_panels)
    if len(ylabels) != len(normalized_panels):
        raise ValueError("ylabels must match panels length.")

    ylabels_non_empty = [ylabels[i] for i in non_empty_indices]

    # ---- Create figure ----
    fig, axes = plt.subplots(Npanels, 1, sharex=True, figsize=(10, 3*Npanels))
    if Npanels == 1:
        axes = [axes]

    all_line_objs = []
    vlines = []

    # ---- Prepare each panel ----
    for ax, sig_list, lg_list, ylabel in zip(axes, panels_non_empty, legends_non_empty, ylabels_non_empty):

        # lines
        lines = []
        for _ in sig_list:
            ln, = ax.plot([], [])
            lines.append(ln)
        all_line_objs.append(lines)

        # y-label
        if ylabel is not None:
            ax.set_ylabel(ylabel)

        # x-limits fixed
        ax.set_xlim(time_axis[0], time_axis[-1])

        # y-limits
        mins = [np.min(s) for s in sig_list]
        maxs = [np.max(s) for s in sig_list]
        ymin, ymax = min(mins), max(maxs)
        if np.isclose(ymin, ymax):
            span = abs(ymin) if ymin != 0 else 1.0
            ymin -= 0.1 * span
            ymax += 0.1 * span
        else:
            pad = 0.05*(ymax - ymin)
            ymin -= pad
            ymax += pad
        ax.set_ylim(ymin, ymax)

        # legend
        #ax.legend(lg_list, loc="upper right")  
        ax.legend(lg_list) ## this way legend location is updated automatically at each frame
        
        # vertical time marker
        vlines.append(ax.axvline(time_axis[0], ls="--", color="k"))

    # bottom xlabel
    if xlabel is not None:
        axes[-1].set_xlabel(xlabel)
        
    # title
    if title is not None:
        axes[0].set_title(title)

    # ---- Animation update ----
    def update(frame):
        xnow = time_axis[frame]

        for p_idx, lines in enumerate(all_line_objs):
            sigs = panels_non_empty[p_idx]
            for s_idx, ln in enumerate(lines):
                ln.set_data(time_axis[:frame+1], sigs[s_idx][:frame+1])
            vlines[p_idx].set_xdata([xnow])

        return [artist for sub in all_line_objs for artist in sub] + vlines

    anim = FuncAnimation(
        fig,
        update,
        frames=len(time_axis),
        blit=False,
        **FuncAnimation_kwargs
    )

    return anim
                   
# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)