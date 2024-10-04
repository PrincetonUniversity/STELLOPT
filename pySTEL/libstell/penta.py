"""
This library provides a python class for reading and handling PENTA3
results data.
"""

import numpy as np
import sys

# Constants
EC = 1.602176634E-19 # Electron charge [C]

# PENTA Class
class PENTA:
    
    def __init__(self, folder_path, plasma=None, Zions=None):
        #folder_path is a path to the folder containnig the following PENTA3 results files:
        # - flows_vs_Er
        # - flows_vs_roa
        # - fluxes_vs_Er
        # - fluxes_vs_roa
        # - Jprl_vs_roa
        # - contra_vs_roa
        # - sigmas_vs_roa
        
        print('\nPENTA class being created...')
        
        self.folder_path = folder_path
        
        #check if a plasma class is given. If not check if Zions is given
        if plasma is None and Zions is None:
            print('Could not create PENTA class: Need to provide a plasma class or an array with char number (Z) of the ions')
            exit(1)
        elif(plasma is not None and Zions is not None):
            print('Both plasma class and Zions provided. Considering plasma class and discarding Zions array')
            self.list_of_species = plasma.list_of_species
            print(f'List of species: {self.list_of_species}')
            self.Zions = [plasma.Zcharge[species] for species in self.list_of_species]
            #remove electrons
            self.Zions = self.Zions[1:]
            print(f'Zions={self.Zions}')
        elif(plasma is not None and Zions is None):
            self.list_of_species = plasma.list_of_species
            print(f'List of species: {self.list_of_species}')
            self.Zions = [plasma.Zcharge[species] for species in self.list_of_species]
            #remove electrons
            self.Zions = self.Zions[1:]
            print(f'Zions={self.Zions}')
        elif(plasma is None and Zions is not None):
            #check if size of Zions is compatible with number of ions in results files
            self.check_size_Zions(Zions)
            self.list_of_species = ['e'] + [f'i{k}' for k in range(1,len(Zions)+1)]
            print(f'List of species: {self.list_of_species}')
            self.Zions = Zions
            print(f'Zions={self.Zions}')
            
        #sets the arrays self.roa_unique and self.Er_search
        self.set_independent_variables()
        
        # sets the dictionaries self.##[root], self.##[root]], self.##[root]],
        # with ## being: roa, Er, Jprl_total, J_BS
        self.set_variables_by_root()
        
        # sets the dictionaries:
        # - self.uprl[species,root]
        # - self.Jprl[species,root]
        # - self.Gamma[species,root]
        # species is one of the species in self.list_of_species and root should be 
        # 'ion_root', 'electron_root' or 'unstable_root'
        self.set_fluxes_flows_by_root()
            
            
    def check_size_Zions(self,Zions):
        #checks if len(Zions) is the same as the number of ions in the file fluxes_vs_Er
        
        filename = self.folder_path + '/fluxes_vs_Er'
            
        penta = np.loadtxt(filename,skiprows=2)        
        num_ion_species = len(penta[0,:]) - 3
        
        if( len(Zions) is not num_ion_species):
            print(f'ERROR: There are {num_ion_species} species, but only provided Zcharge of {len(Zions)} !!')
            exit(0)
    
    def set_independent_variables(self):
        # sets the arrays self.roa_unique, self.Er_search
        # and sets the integer self.Smax
        # note that Er_search is in V/cm
        
        filename = self.folder_path + '/fluxes_vs_Er'
            
        penta = np.loadtxt(filename,skiprows=2) 
        
        self.roa_unique = np.unique(penta[:,0])
        self.Er_search = np.unique(penta[1,:])
        
        filename = self.folder_path + '/flows_vs_roa'
            
        penta = np.loadtxt(filename,skiprows=2)
        
        number_flows_per_species = len(penta[0,3:]) / len(self.list_of_species)
        
        self.Smax = int(number_flows_per_species - 1)
        
        print(f'Smax={self.Smax}')
        
    def set_variables_by_root(self):
        # sets the arrays self.##_ion_root, self.##_electron_root, self.##_unstable_root,
        # with ## being: roa, Er, Jprl_total, JBS
        
        from itertools import groupby
        from collections import defaultdict
        
        filename = self.folder_path + '/Jprl_vs_roa'
            
        penta = np.loadtxt(filename,skiprows=2)
        
        roa = penta[:,0]
        Er = penta[:,1]
        Jprl_total = penta[:,-2]
        JBS = penta[:,-1]
        
        self.roa = defaultdict(list)
        self.Er = defaultdict(list)
        self.Jprl_total = defaultdict(list)
        self.JBS = defaultdict(list)

        i=0
        for _, group in groupby(roa):
            num_roots = len( list(group) )
            
            if(num_roots == 1):
                self.roa['ion_root'].append(roa[i])
                self.Er['ion_root'].append(Er[i])
                self.Jprl_total['ion_root'].append(Jprl_total[i])
                self.JBS['ion_root'].append(JBS[i])
            elif(num_roots ==3):
                # ion root
                self.roa['ion_root'].append(roa[i])
                self.Er['ion_root'].append(Er[i])
                self.Jprl_total['ion_root'].append(Jprl_total[i])
                self.JBS['ion_root'].append(JBS[i])
                # unstble root
                self.roa['unstable_root'].append(roa[i+1])
                self.Er['unstable_root'].append(Er[i+1])
                self.Jprl_total['unstable_root'].append(Jprl_total[i+1])
                self.JBS['unstable_root'].append(JBS[i+1])
                # electron root
                self.roa['electron_root'].append(roa[i+2])
                self.Er['electron_root'].append(Er[i+2])
                self.Jprl_total['electron_root'].append(Jprl_total[i+2])
                self.JBS['electron_root'].append(JBS[i+2])
            else:
                raise ValueError(f"How come you have {num_roots} roots ??")   
            i += num_roots        
  
    def set_fluxes_flows_by_root(self):
        # sets the dictionaries self.uprl[species,root], self.Jprl[species,root] and self.Gamma[species,root]
        # they dictionaries return arrays
        # species is one of the species in self.list_of_species and root should be 'ion_root', 'electron_root' or 'unstable_root'
        
        from itertools import groupby
        from collections import defaultdict
        
        #get uprl0 flows for all species
        filename = self.folder_path + '/flows_vs_roa'   
        penta = np.loadtxt(filename,skiprows=2)
        roa = penta[:,0]
        uprl0_penta = penta[:,3:-1:(self.Smax+1)]
        
        #get Jprl flows for all species
        filename = self.folder_path + '/Jprl_vs_roa'   
        penta = np.loadtxt(filename,skiprows=2)
        Jprl_penta = penta[:,3:(3+len(self.list_of_species))]
        
        #get particle fluxes for all species
        filename = self.folder_path + '/fluxes_vs_roa'
        penta = np.loadtxt(filename,skiprows=2)
        Gamma_penta = penta[:,np.r_[3,5:(5+len(self.Zions))]]
        
        self.uprl = defaultdict(list)
        self.Jprl = defaultdict(list)
        self.Gamma = defaultdict(list)
        
        for k,species in enumerate(self.list_of_species):
        
            i=0
            for _, group in groupby(roa):
                num_roots = len( list(group) )
                
                if(num_roots == 1):
                    self.uprl[species,'ion_root'].append(uprl0_penta[i,k])
                    self.Jprl[species,'ion_root'].append(Jprl_penta[i,k])
                    self.Gamma[species,'ion_root'].append(Gamma_penta[i,k])
                elif(num_roots ==3):
                    # ion root
                    self.uprl[species,'ion_root'].append(uprl0_penta[i,k])
                    self.Jprl[species,'ion_root'].append(Jprl_penta[i,k])
                    self.Gamma[species,'ion_root'].append(Gamma_penta[i,k])
                    # unstable root
                    self.uprl[species,'unstable_root'].append(uprl0_penta[i+1,k])
                    self.Jprl[species,'unstable_root'].append(Jprl_penta[i+1,k])
                    self.Gamma[species,'unstable_root'].append(Gamma_penta[i+1,k])
                    # electron root
                    self.uprl[species,'electron_root'].append(uprl0_penta[i+2,k])
                    self.Jprl[species,'electron_root'].append(Jprl_penta[i+2,k])
                    self.Gamma[species,'electron_root'].append(Gamma_penta[i+2,k])
                else:
                    raise ValueError(f"How come you have {num_roots} roots ??")   
                i += num_roots
                
    def plot_Er_vs_roa(self,which_root='all',plot=True):
        # plots ambipolar Er vs roa
        # which_root is: 'ion_root', 'unstable_root', 'electron_root' or 'all'
        # if 'all', plots all roots in the same figure
        
        import matplotlib.pyplot as plt
        
        plt.rc('font', size=18)
        fig, ax = plt.subplots(figsize=(11,8))
        
        if which_root is ('ion_root' or 'electron_root' or 'unstable_root'):
            ax.plot(self.roa[which_root],self.Er[which_root],'.-',label=which_root)
        elif which_root is 'all':
            ax.plot(self.roa['ion_root'],self.Er['ion_root'],'.-',label='ion_root')
            ax.plot(self.roa['unstable_root'],self.Er['unstable_root'],'.-',label='unstable_root')
            ax.plot(self.roa['electron_root'],self.Er['electron_root'],'.-',label='electron_root')
        else:
            print('ERROR: which_root can only be ion_root, electron_root, unstable_root or all')
            exit(0)
                    
        ax.set_ylabel(r'Er [V/cm]')
        ax.set_xlabel(r'r/a')
        ax.legend(fontsize=12)
        ax.grid()
        plt.legend()
        if plot:
            plt.show()
            
    def plot_Gamma_vs_roa(self,which_root='all',which_species='all',plot=True):
        # plots ambipolar particle flux vs roa
        # which_root is: 'ion_root', 'unstable_root', 'electron_root' or 'all'
        # if 'all', plots all roots in the same figure
        # which_species is any species in self.list_of_species
        
        import matplotlib.pyplot as plt
        
        if which_species is 'all': 
            plotting_species = self.list_of_species
        elif which_species not in self.list_of_species: 
            print(f'Error: species {which_species} is not valid. Pick species from: {self.list_of_species}')
        else: 
            plotting_species = [which_species]
        
        plt.rc('font', size=18)
        fig, ax = plt.subplots(figsize=(11,8))
        for species in plotting_species:
        
            if which_root is ('ion_root' or 'electron_root' or 'unstable_root'):
                ax.plot(self.roa[which_root],self.Gamma[species,which_root],'.-',label=which_root+', '+species)
            elif which_root is 'all':
                ax.plot(self.roa['ion_root'],self.Gamma[species,'ion_root'],'.-',label='ion_root'+', '+species)
                ax.plot(self.roa['unstable_root'],self.Gamma[species,'unstable_root'],'.-',label='unstable_root'+', '+species)
                ax.plot(self.roa['electron_root'],self.Gamma[species,'electron_root'],'.-',label='electron_root'+', '+species)
            else:
                print('ERROR: which_root can only be ion_root, electron_root, unstable_root or all')
                exit(0)
                    
        ax.set_ylabel(r'$\Gamma~[m^{-2}s^{-1}]$')
        ax.set_xlabel(r'r/a')
        ax.set_title('ambipolar particle fluxes')
        ax.legend(fontsize=12)
        ax.grid()
        plt.legend()
        if plot:
            plt.show()
            
    def plot_flows_vs_roa(self,which_root='all',which_species='all',plot=True):
        # plots ambipolar flows vs roa
        # which_root is: 'ion_root', 'unstable_root', 'electron_root' or 'all'
        # if 'all', plots all roots in the same figure
        # which_species is any species in self.list_of_species
        
        import matplotlib.pyplot as plt
        
        if which_species is 'all': 
            plotting_species = self.list_of_species
        elif which_species not in self.list_of_species: 
            print(f'Error: species {which_species} is not valid. Pick species from: {self.list_of_species}')
        else: 
            plotting_species = [which_species]
        
        plt.rc('font', size=18)
        fig, ax = plt.subplots(figsize=(11,8))
        for species in plotting_species:
        
            if which_root is ('ion_root' or 'electron_root' or 'unstable_root'):
                ax.plot(self.roa[which_root],self.uprl[species,which_root],'.-',label=which_root+', '+species)
            elif which_root is 'all':
                ax.plot(self.roa['ion_root'],self.uprl[species,'ion_root'],'.-',label='ion_root'+', '+species)
                ax.plot(self.roa['unstable_root'],self.uprl[species,'unstable_root'],'.-',label='unstable_root'+', '+species)
                ax.plot(self.roa['electron_root'],self.uprl[species,'electron_root'],'.-',label='electron_root'+', '+species)
            else:
                print('ERROR: which_root can only be ion_root, electron_root, unstable_root or all')
                exit(0)
                    
        ax.set_ylabel(r'$\left<u_{\parallel}B\right>/\left<B^2\right>~[m~s^{-1}~T^{-1}]$')
        ax.set_xlabel(r'r/a')
        ax.set_title('ambipolar parallel flows')
        ax.legend(fontsize=12)
        ax.grid()
        plt.legend()
        if plot:
            plt.show()
            
    def plot_currents_vs_roa(self,which_root='all',plot=True):
        # plots ambipolar parallel currents vs roa
        # which_root is: 'ion_root', 'unstable_root', 'electron_root' or 'all'
        # if 'all', plots all roots in the same figure
        
        import matplotlib.pyplot as plt
        
        plt.rc('font', size=18)
        fig, ax = plt.subplots(figsize=(11,8))
        for species in self.list_of_species:
        
            if which_root is ('ion_root' or 'electron_root' or 'unstable_root'):
                ax.plot(self.roa[which_root],np.array(self.Jprl[species,which_root])/1000,'.-',label=which_root+', '+species)
            elif which_root is 'all':
                ax.plot(self.roa['ion_root'],np.array(self.Jprl[species,'ion_root'])/1000,'.-',label='ion_root'+', '+species)
                ax.plot(self.roa['unstable_root'],np.array(self.Jprl[species,'unstable_root'])/1000,'.-',label='unstable_root'+', '+species)
                ax.plot(self.roa['electron_root'],np.array(self.Jprl[species,'electron_root'])/1000,'.-',label='electron_root'+', '+species)
            else:
                print('ERROR: which_root can only be ion_root, electron_root, unstable_root or all')
                exit(0)
                
        # now plot Jprl_total and JBS (when 'all', it doesn't plot, otherwise the plot gets super messy)
        if which_root is ('ion_root' or 'electron_root' or 'unstable_root'):
                ax.plot(self.roa[which_root],np.array(self.Jprl_total[which_root])/1000,'.-',label='total, '+which_root)
                ax.plot(self.roa[which_root],np.array(self.JBS[which_root])/1000,'.-',label='BS, '+which_root)   
       
        ax.set_ylabel(r'$\left<J_{\parallel}\right>~[kA~m^{-2}]$')
        ax.set_xlabel(r'r/a')
        ax.set_title('ambipolar parallel currents')
        ax.legend(fontsize=12)
        ax.grid()
        plt.legend()
        if plot:
            plt.show()
            
    def plot_fluxes_vs_Er(self,roa_user,plot=True):
        # plots fluxes*Z as function of Er
        
        import matplotlib.pyplot as plt
        from collections import defaultdict
        
        filename = self.folder_path + '/fluxes_vs_Er'
            
        penta = np.loadtxt(filename,skiprows=2)
        
        roa = penta[:,0]
        Er = penta[:,1]
        gamma_e = penta[:,2]       
        
        gamma_i_tot = np.sum(penta[:,3:]*self.Zions,axis=1)
        
        Er_dict = defaultdict(list)
        Gamma_e_dict = defaultdict(list)
        Gamma_i_dict = defaultdict(list)
        
        for r,er,ge,gi in zip(roa,Er,gamma_e,gamma_i_tot):
            Er_dict[r].append(er)
            Gamma_e_dict[r].append(ge)
            Gamma_i_dict[r].append(gi)
        
        #conver roa_user to array in case it is not yet
        roa_user = np.array(roa_user)
        
        for r_user in roa_user:
            
            #check that roa_user is inside the limits
            #if yes, plot the closest to roa
            if(r_user <= np.min(roa) or r_user>= np.max(roa)):
                print(f'ERROR!! The provided roa={r_user} is outside the interval of PENTA data: [{np.min(roa)},{np.max(roa)}]')
                exit(0)        
            else:
                roa_closest = self.roa_unique[ np.argmin(np.abs(self.roa_unique-r_user)) ]
            
            plt.rc('font', size=16)
            fig=plt.figure(figsize=(8,6))
            ax = fig.add_subplot(111)
            ax.plot(Er_dict[roa_closest],Gamma_e_dict[roa_closest],label=r'$\Gamma_e$')
            ax.plot(Er_dict[roa_closest],Gamma_i_dict[roa_closest],label=r'$\Sigma~Z_i\Gamma_i$')
            ax.set_xlabel(r'Er [V/cm]')
            ax.set_ylabel(r'$\Gamma~~[\text{m}^{-2}\,\text{s}^{-1}]$')
            ax.set_title(f'particle fluxes, r/a={roa_closest}')
            ax.set_yscale('log')
            ax.legend(fontsize=12)
            ax.grid()
            plt.legend()
        if plot:
            plt.show()
            
    def plot_plasma_profiles_check(self,plot=True):
        # plots plasma profiles as in plasma_profiles_check file provided by PENTA3
        # this checks if splines of profiles and their derivatives were done correctly
        
        import matplotlib.pyplot as plt
        import numpy as np
        
        filename = self.folder_path+'/plasma_profiles_check'
            
        penta = np.loadtxt(filename,skiprows=2)
        
        roa = penta[:,0]
        
        num_ion_species = len(self.Zions)
        
        species_list = ['e'] + [f'i{j}' for j in range(1, num_ion_species + 1)]
        
        #densities: electrons, i1, i2, ...
        densities = penta[:,np.r_[2,(5+num_ion_species):(5+2*num_ion_species)]]
        
        #temperatures: electrons, i1, i2, ...
        temperatures = penta[:,np.r_[1,5:(5+num_ion_species)]]
        
        #grad_densitites: electrons, i1, i2, ...
        grad_densitites = penta[:,np.r_[3,(5+2*num_ion_species):(5+3*num_ion_species)]]
        
        #grad_tempratures: electrons, i1, i2, ...
        grad_temperatures = penta[:,np.r_[4,(5+3*num_ion_species):(5+4*num_ion_species)]]
        
        plt.rc('font', size=18)
        fig, ax = plt.subplots(2,2,figsize=(19,12))
        for k,s in enumerate(species_list): ax[0,0].plot(roa,densities[:,k],label=s)
        ax[0,0].set_ylabel(r'$n~~[m^{-3}]$')
        for k,s in enumerate(species_list): ax[0,1].plot(roa,temperatures[:,k]/1000,label=s)
        ax[0,1].set_ylabel(r'$T$  [keV]')
        for k,s in enumerate(species_list): ax[1,0].plot(roa,grad_densitites[:,k],label=s)
        ax[1,0].set_ylabel(r'$dn/dr~~[m^{-4}]$')
        for k,s in enumerate(species_list): ax[1,1].plot(roa,grad_temperatures[:,k]/1000,label=s)
        ax[1,1].set_ylabel(r'$dT/dr$  [keV/m]')
        
        ax[1,0].set_xlabel(r'r/a')
        ax[1,1].set_xlabel(r'r/a')
        
        for a in ax.flat:
            a.grid(True)
            a.legend()
            #a.legend(fontsize=12)
        if plot:
            plt.show()
            
    def plot_conductivity(self,aspect_ratio=None,plot=True):
        # plots parallel conducitivity as given by PENTA (when ran in 'SN' mode)
        # if aspect_ratio is given, spitzer-NEO is computed
        
        import matplotlib.pyplot as plt
        import numpy as np
        
        filename = self.folder_path+'/sigmas_vs_roa'
        
        try: 
            penta = np.loadtxt(filename,skiprows=2)
        except:
            print('Could not read file sigma_vs_roa. Probably PENTA was not run in "DKES" mode...?')
            exit(1)
        
        roa = penta[:,0]
        sigma_par = penta[:,2]
        sigma_par_Spitzer = penta[:,3]

        plt.rc('font', size=18)
        fig=plt.figure(figsize=(11,8))
        ax = fig.add_subplot(111)
        ax.plot(roa,sigma_par,label='PENTA')
        ax.plot(roa,sigma_par_Spitzer,label='Spitzer')
        if aspect_ratio is not None:
            sigma_par_Spitzer_NEO = sigma_par_Spitzer*(1-np.sqrt(roa/aspect_ratio))**2
            ax.plot(roa,sigma_par_Spitzer_NEO,'--',label='Spitzer-NEO')
        ax.set_ylabel(r'$\sigma_{\parallel}~~[\Omega^{-1}~m^{-1}]$')
        ax.set_xlabel(r'r/a')
        ax.grid()
        plt.legend()
        
        plt.rc('font', size=18)
        fig=plt.figure(figsize=(11,8))
        ax = fig.add_subplot(111)
        ax.plot(roa,1/sigma_par,label='PENTA')
        ax.plot(roa,1/sigma_par_Spitzer,label='Spitzer')
        if aspect_ratio is not None:
            ax.plot(roa,1/sigma_par_Spitzer_NEO,'--',label='Spitzer-NEO')
        ax.set_ylabel(r'$\eta_{\parallel}~~[\Omega~m]$')
        ax.set_xlabel(r'r/a')
        ax.set_yscale('log')
        ax.grid()
        plt.legend()
        
        if plot:
            plt.show()

# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)