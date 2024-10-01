"""
This library provides a python class for reading and handling DKES
results data. It also prepares the input PENTA coefficients
"""

import numpy as np
import sys

# Constants
EC = 1.602176634E-19 # Electron charge [C]

# DKES Class
class DKES:
    
    def __init__(self, surface, wout_file, eps_rel=0.03):
        self.eps_rel = eps_rel
        
        #check surface is an integer
        if not isinstance(surface,int):
            print('ERROR: surface must be an integer')
            exit(1)
        else:
            self.surface = surface
            
        #read wout file and get self.Bsq and self.roa
        self.get_VMEC_quantities(wout_file)
        
    def read_DKES_results(self,filename):
        
        dkes = np.loadtxt(filename,skiprows=2)
        
        self.cmul = dkes[:,0]
        self.efield = dkes[:,1]
        self.L11m = dkes[:,4]
        self.L11p = dkes[:,5]
        self.L31m = dkes[:,6]
        self.L31p = dkes[:,7]
        self.L33m = dkes[:,8]
        self.L33p = dkes[:,9]
        
        #Put Lijm and Lijp in dictionary
        self.Lm = {
            'L11' : self.L11m,
            'L31' : self.L31m,
            'L33' : self.L33m
        }
        
        self.Lp = {
            'L11' : self.L11p,
            'L31' : self.L31p,
            'L33' : self.L33p
        }
        
        #assuming that data is regular:
        # get ncmul, nefield
        self.ncmul = np.sum(self.efield == self.efield[0])
        self.nefield = int(np.round(len(self.cmul)/self.ncmul,decimals=0))
        
        #Check convergence of coefficients
        self.check_convergence()
        
        #Take average of coefficients
        L11  = 0.5*(self.L11m+self.L11p)
        L31  = 0.5*(self.L31m+self.L31p)
        L33  = 0.5*(self.L33m+self.L33p)
        
        #compute Dij_star correcting D13 and D33 with the B factors (see J.Lore documentation on PENTA)
        self.D11_star = L11
        self.D31_star = L31 * np.sqrt(self.Bsq)
        self.D33_star = L33 * self.Bsq
        
        #according to J. Lore documentation and also C. Beidler...
        # this is also what is done internally in PENTA
        self.D13_star = -self.D31_star
        
    def check_convergence(self):
    
        self.rel_error_L11 = abs(self.L11m-self.L11p) / abs(self.L11p)
        self.rel_error_L31 = abs(self.L31m-self.L31p) / abs(self.L31p)
        self.rel_error_L33 = abs(self.L33m-self.L33p) / abs(self.L33p)
        
        self.rel_error = {
            'L11' : self.rel_error_L11,
            'L31' : self.rel_error_L31,
            'L33' : self.rel_error_L33
        }
        
        for Lvar in self.rel_error:
            #give indexes where rel error larger than epsilon
            idx_eps = np.where(self.rel_error[Lvar] > self.eps_rel)[0]
            
            Lvarm = Lvar+'m'
            Lvarp = Lvar+'p'
            
            #in case there is any index, print in the command line where this happens
            if idx_eps.size > 0:
                print('\nRELATIVE ERROR OF '+Lvar+' >',self.eps_rel*100,'%')
                for idx in idx_eps:
                    print(f'{"cmul":<12} {"efield":<12} {Lvarm:<12} {Lvarp:<12} {"rel-error (%)":<15}')
                    print(f'{self.cmul[idx]:<12.3E} {self.efield[idx]:<12.3E} {self.Lm[Lvar][idx]:<12.4E} {self.Lp[Lvar][idx]:<12.4E} {self.rel_error[Lvar][idx]*100:<15.2f}')
                print(f'Making the average anyway: {Lvar} = 0.5*({Lvarp} + {Lvarm})')
            
    def plot_DKES_coeffs(self):
        # plots the species-independent L11, L13 and L33 
        
        import matplotlib.pyplot as pyplot
        
        var_names = {
                'D11_star': '$D_{11}^*~~[m^{-1}~T^{-2}]$',
                'D31_star': '$D_{31}^*$',
                'D33_star': '$D_{33}^*~~[m~T^2]$'
                }
    
        for plot_var in ['D11_star', 'D31_star', 'D33_star']:
            
            yplot = getattr(self,plot_var)
                
            px = 1/pyplot.rcParams['figure.dpi']
            pyplot.rc('font', size=20)
            #pyplot.rc('legend', fontsize=24)
            fig=pyplot.figure(figsize=(1024*px,768*px))
            ax = fig.add_subplot(111)
            for i in range(self.nefield):
                i1 = i*self.ncmul
                i2 = i1 + self.ncmul
                # plot without error bar
                #ax.plot(self.cmul[i1:i2],yplot[i1:i2],marker='+',label=rf'$E_s/v$={self.efield[i1]:3.1E}',linewidth=4,markersize=18)
                # plot with error bar
                [yerr_lower, yerr_upper] = self.compute_yerr(yplot[i1:i2],self.Lm[plot_var][i1:i2],self.Lp[plot_var][i1:i2])
                ax.errorbar(self.cmul[i1:i2],yplot[i1:i2],yerr=[yerr_lower,yerr_upper],fmt='-o',label=rf'$E_s/v$={self.efield[i1]:3.1E}',capsize=5, elinewidth=2, markeredgewidth=2)
                    
            ax.set_xlabel(r'$\nu/v\,\,[\text{m}^{-1}]$')
            ax.set_ylabel(f'{var_names[plot_var]}')
            ax.set_xscale('log')
            if(plot_var=='D11_star' or plot_var=='D33_star'):
                ax.set_yscale('log')
            if(plot_var=='D31_star'):
                ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
            ax.set_title(f'r/a={self.roa:.2f}')
            ax.legend(fontsize=12)
            ax.grid()
    
        pyplot.show()
        
    def compute_yerr(self,y,ym,yp):
    # Computes the distance between the central value 'y' and its upper and lower limits
    # A priori we don't know which one of these is the upper and lower since they change
    # depending on their sign
    # The output of this function is used for the argument of the maplotlib errorbar function

        # Check if arrays have different sizes
        if len(y) != len(yp) or len(y) != len(ym):
            print("Error: Arrays 'y', 'yp', and 'ym' must have the same size.")
            sys.exit(1)  # Exit the program with a status code of 1 (indicates error)
            
        yerr_lower = []
        yerr_upper = []
        
        for yi, ypi, ymi in zip(y, yp, ym):
            if(ymi>yi and ypi<yi):
                yerr_lower.append(yi-ypi)
                yerr_upper.append(ymi-yi)
                
            elif(ypi>yi and ymi<yi):
                yerr_lower.append(yi-ymi)
                yerr_upper.append(ypi-yi)
            else:
                #print('\nym and yp are equal!')
                yerr_lower.append(yi-ymi)
                yerr_upper.append(ypi-yi)
        
        return [yerr_lower, yerr_upper]
            
            
    # def compute_PENTA_coeffs(self,wout_filename):
    #     # Computes PENTA input coefficients lstar, mstar and nstar
    #     # size of lstar, mstar, nstar is n_cmul x n_efield
    #     # Check DKES/PENTA documentation to see the definition of lstar, mstar, nstar
    #     # In the documentation, D_ij^* corresponds to self.Lij, which are the species-independent DKES coefficientes
    #     print('\n#############################################################################')
    #     print('###################   Computing PENTA coefficients  #########################')
    #     print('#############################################################################')
        
    #     ######  WARNING: as of now this assumes an hydrogen plasma, qa=e_charge  #####
    #     print('\nWARNING: THIS ASSUMES A PLASMA WITH Z=1, qa=echarge')
        
    #     # Read Pfirsch-Schluter flow from external file
    #     # To do later...
    #     print('\nFailed to read Pfirsch-Schluter flow, <U^2>, from external file')
    #     print('Assuming <U^2>=0\n')
    #     self.Usq = 0.0
        
    #     #compute PENTA lstar
    #     aux = 1 - 1.5*self.cmul*self.L33/self.Bsq
    #     self.lstar = self.L11 - (2./3.)*self.cmul*self.Usq + (1.5*self.cmul*self.L31*self.L31/self.Bsq)/aux
    #     self.lstar = self.lstar / (EC*EC)
        
    #     #compute PENTA mstar
    #     self.mstar = self.cmul*self.cmul*self.L33 / aux
        
    #     #compute PENTA nstar
    #     self.nstar = self.cmul*self.L31 / aux
    #     self.nstar = self.nstar / EC
    
    def compute_PENTA1_coeffs(self):
        # Computes PENTA input coefficients lstar, mstar and nstar
        # size of lstar, mstar, nstar is n_cmul x n_efield
        # Check DKES/PENTA documentation to see the definition of lstar, mstar, nstar
        # In the documentation, D_ij^* corresponds to self.Lij, which are the species-independent DKES coefficientes
        print('\n#############################################################################')
        print('###################   Computing coefficients for PENTA1/v2.0 ##################')
        print('###############################################################################')
        
        ######  WARNING: as of now this assumes an hydrogen plasma, qa=e_charge  #####
        print('\nWARNING: THIS ASSUMES A PLASMA WITH Z=1, qa=echarge')
        
        # Read Pfirsch-Schluter flow from external file
        # To do later...
        print('\nFailed to read Pfirsch-Schluter flow, <U^2>, from external file')
        print('Assuming <U^2>=0\n')
        self.Usq = 0.0
        
        aux = 1 - 1.5*self.cmul*self.D33_star/self.Bsq
        
        #compute PENTA lstar
        self.lstar = self.D11_star - (2./3.)*self.cmul*self.Usq + (1.5*self.cmul*self.D13_star*self.D13_star/self.Bsq)/aux
        self.lstar = self.lstar / (EC*EC)
        
        #compute PENTA mstar
        self.mstar = self.cmul*self.cmul*self.D33_star / aux

        #compute PENTA nstar
        self.nstar = self.cmul*self.D13_star / aux
        self.nstar = self.nstar / EC
    
    def get_VMEC_quantities(self,wout_file):
        
        ##########################################################
        import netCDF4 as nc
        try:
            dataset = nc.Dataset(wout_file, 'r')
            Bsq = dataset.variables['bdotb'][:]
            phi = dataset.variables['phi'][:]
            bmnc = dataset.variables['bmnc'][:]
            dataset.close()
        except:
            print('\nERROR: Could not read wout file')
            sys.exit(0)
        ##########################################################
        
        self.Bsq = Bsq[self.surface-1]
        self.roa = np.sqrt(phi[self.surface-1]/phi[-1])
        # <|B|> = B00
        self.B0 = bmnc[self.surface-1,0]
        
        print(f'r/a = sqrt(PHI/PHIEDGE) = {self.roa}')
        print(f'Bsq = {self.Bsq}')
        print(f'B0 = {self.B0}')
        
     
    def plot_PENTA1_coeffs(self):
        # creates 3 graphs: lstar, mstar and nstar vs cmul (for each efield)
        
        import matplotlib.pyplot as pyplot
        
        for plot_var in ['lstar', 'mstar', 'nstar']:
            
            yplot = getattr(self,plot_var)
                
            px = 1/pyplot.rcParams['figure.dpi']
            pyplot.rc('font', size=24)
            pyplot.rc('legend', fontsize=24)
            fig=pyplot.figure(figsize=(1024*px,768*px))
            ax = fig.add_subplot(111)
            for i in range(self.nefield):
                i1 = i*self.ncmul
                i2 = i1 + self.ncmul
                ax.plot(self.cmul[i1:i2],yplot[i1:i2],marker='+',label=rf'$E_s/v$={self.efield[i1]:3.1E}',linewidth=4,markersize=18)
            ax.set_xlabel(r'$\nu/v [m^{-1}]$')
            ax.set_ylabel(f'PENTA {plot_var}')
            ax.set_xscale('log')
            if(plot_var=='lstar' or plot_var=='mstar'):
                ax.set_yscale('log')
            if(plot_var=='nstar'):
                ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
            ax.set_title(f'r/a={self.roa:.2f}')
            ax.grid()
            ax.legend(fontsize=12)
     
        pyplot.show()
        
    def h2f(self,var_half):
        """Half to full grid
        Copied from vmec.py

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
    
    def write_PENTA3_coeffs_to_files(self,where_to):
        # name of the files are 'D11_star_##' , 'D13_star_##', 'D33_star_##'
        
        #where_to save -- path should not have final /
        
        #checks that data is regular, i.e.: for each efield there are always the same cmul
        # if data is regular the following arrays are computed:
        # self.cmul_regular
        # self.efield_regular
        
        self.check_data_is_regular(self.cmul,self.efield)
        
        for var in ['D11_star', 'D31_star', 'D33_star']:
            
            filename = where_to + '/' + var + '_surface_' + f'{self.surface}'
            y = getattr(self,var)

            combined = np.concatenate((self.cmul_regular,self.efield_regular,y))
            
            #create file
            with open(filename, 'w') as file:
                file.write(f'{self.ncmul} {self.nefield}\n')
                for value in combined:
                    file.write(f'{value:.10e}\n')
    
    def write_PENTA1_coeffs_to_files(self,where_to):
        # name of the files are 'lstar_lijs_##' , 'mstar_lijs_##', 'nstar_lijs_##'
        
        #where_to save -- path should not have final /
        
        #checks that data is regular, i.e.: for each efield there are always the same cmul
        # if data is regular the following arrays are computed:
        # self.cmul_regular
        # self.efield_regular
        
        self.check_data_is_regular(self.cmul,self.efield)
        
        for var in ['lstar', 'mstar', 'nstar']:
            
            filename = where_to + '/' + var + '_lijs_' + 'surface_' + f'{self.surface}'
            y = getattr(self,var)

            combined = np.concatenate((self.cmul_regular,self.efield_regular,y))
            
            #create file
            with open(filename, 'w') as file:
                file.write(f'{self.ncmul} {self.nefield}\n')
                for value in combined:
                    file.write(f'{value:.10e}\n')            
                
    def check_data_is_regular(self,cmul,efield):
        # checks that data is regular, i.e.: for each efield there are always the same cmul
        # if data is regular the following arrays are computed:
        # self.cmul_regular
        # self.efield_regular
        
        # Get the unique values and counts of efield
        unique_efields, counts = np.unique(efield, return_counts=True)
        
        if not np.all(counts == counts[0]):
            print('ERROR: Each cmul does not have the same number of efields. Cannot proceed...')
            exit(1)
        else:
            N = counts[0]
            expected_cmul = cmul[:N]  # First block of cmul
              
        if(N != self.ncmul):
            print('Error: Disparity between counts and ncmul. Cannot proceed...')

        # Check each block of cmul corresponding to each unique efield
        for unique_value in unique_efields:
            cmul_block = cmul[efield == unique_value]  # Extract cmul values for current efield
            
            # Check if the cmul block matches the expected pattern
            if not np.array_equal(cmul_block, expected_cmul):
                print('ERROR: data is not regular. cmul and efield are not meshgrid-like. Cannot proceed')
                exit(1)
            else:
                self.efield_regular = unique_efields
                self.cmul_regular   = expected_cmul
                
    def set_cmul_species(self,K,make_plot=False):
        
        import matplotlib.pyplot as plt
        
        self.cmul_species = {}
        for species in self.plasma_class.list_of_species:
            cmul_temp = []
            for k in K:
                vth = self.plasma_class.get_thermal_speed(species,self.roa)
                vparticle = vth * np.sqrt(k)
                nu = self.plasma_class.get_collisionality(species,self.roa,vparticle)
                cmul_temp.append( nu / vparticle )
            
            self.cmul_species[species] = np.array(cmul_temp)
        
        if make_plot==True:
            # plot here cmul vs K for all species
            cmul_min = np.min(self.cmul)
            cmul_max = np.max(self.cmul)
            fig, ax = plt.subplots(figsize=(13,11))
            for species in self.plasma_class.list_of_species:
                plt.rc('font', size=18)
                plt.plot(K,self.cmul_species[f'{species}'],'o-',label=f'{species}')
                plt.plot(K,np.full_like(K,cmul_min),'-r',linewidth=2)
                plt.plot(K,np.full_like(K,cmul_max),'-r')
                ax.set_yscale('log')
                ax.set_ylabel(r'$\nu_D/v~~[m^{-1}]$')
                ax.set_xlabel(f'K')
                ax.set_title(f'r/a={self.roa:.2f}')
                ax.legend()      
                ax.grid()
            plt.show()
            
    def plot_Erv_K_species(self,Er_V_cm,plasma_class,K=None,Er_v_resonance=None):
        #makes plot of Er/v as a function of K for all species
        # Er_V_cm is in V/cm and it can wether be an array or a double
        
        import matplotlib.pyplot as plt
        
        #if K is not given, check if self.K exists
        if K is None:
            try:
                K = self.K
            except:
                print('ERROR!! K should be provided or defined before...')
                exit()
                
        # if Er_V_cm is a scalar, convert it to array
        if np.isscalar(Er_V_cm):
            Er_V_cm = np.array([Er_V_cm])
        else:
            Er_V_cm = np.array(Er_V_cm)              
        
        # convert Er to SI and take absolute value
        Er_V_m = np.abs(Er_V_cm * 100)
        
        for Er in Er_V_m:
            efield_min = np.min(self.efield)
            efield_max = np.max(self.efield)
            plt.rc('font', size=18)
            fig=plt.figure(figsize=(11,8))
            ax = fig.add_subplot(111)
            for species in plasma_class.list_of_species:
                
                #get thermal speed of species
                vth = plasma_class.get_thermal_speed(species,self.roa)
                Er_over_v = Er / (vth*np.sqrt(K))
                
                ax.plot(K,Er_over_v ,'o-',label=f'{species}',markersize=1) 
                ax.set_yscale('log')
                ax.set_ylabel(r'|Er/v|')
                ax.set_xlabel(f'K')
                ax.set_title(f'r/a={self.roa:.2f}, |Er|={Er/100} V/cm')      
                ax.grid()
            ax.plot(K,np.full_like(K,efield_min),'-r',linewidth=3,label='DKES lim')
            ax.plot(K,np.full_like(K,efield_max),'-r',linewidth=3)
            if Er_v_resonance is not None:
                ax.plot(K,np.full_like(K,Er_v_resonance),'--r',linewidth=3,label='Er/v res')
            plt.legend()
            plt.show()
    
    def set_PENTA1_integrands_energy_conv(self,intj,plasma_class,make_plots=True):
        # This function computes the integrand of the energy convolution as in PENTA for each efield
        # and plots it as function of K
        # Integrand = sqrt(K) * exp(-K) * (K-5/2)^{intj-1} * [lstar,m,star,nstar] * K^{3/2}
        # This function requires computing collisionality nu_D
        # We also spline interpolate lstar,mstar and star as function of cmul for each efield
        
        import matplotlib.pyplot as plt
        
        Kmin = 1e-4      #PENTA default value
        Kmax = 20        #PENTA default value
        numKsteps = 2000  #PENTA default value
        K = np.linspace(Kmin,Kmax,numKsteps)
        #tK = np.linspace(Kmin,0.6,80)
        #K = np.unique( np.concatenate( (K,tK) ) )
        #K = np.sort(K)
        self.K = K
        
        self.plasma_class = plasma_class
        
        # Computes dicitionary of arrays self.cmul_species. 
        # Contains cmul for each species for array K
        self.set_cmul_species(K,make_plot=make_plots)
        
        # Set integrands = l/m/n-star * fix func
        # this creates dictionary of arrays: self.lstar_integrand, self.nstar_integrand, self.mstar_integrand
        # for instance, self.lstar_integrand['electrons'][3] gives the arrays of integrand (as function of K) 
        # for l* for electrons for the 4th (3+1) electric field 
        self.set_integrands(intj,K,make_plot=make_plots)  
            
            
    def get_fix_func(self,K,intj,make_plot=False):
        
        import matplotlib.pyplot as plt
        
        fix_func = np.sqrt(K) * np.exp(-K) * (K-2.5)**(intj-1) * K**1.5
        
        if make_plot is True:
            fig, ax = plt.subplots(figsize=(8,6))
            ax.plot(K,fix_func,'o-')
            #ax.set_yscale('log')
            ax.set_ylabel(r'$\sqrt{K}e^{-K}\left(K-5/2\right)^{j-1}\,K^{3/2}$')
            ax.set_xlabel(f'K')   
            ax.set_title(f'j={intj}')
            ax.grid()    
            plt.show()
        
        return fix_func        
        
                
    def set_integrands(self,intj,K,make_plot=False):
        
        from scipy.interpolate import interp1d
        import matplotlib.pyplot as plt
        from collections import defaultdict
        
        print('\n ############################################################')
        print('############# COMPUTING INTEGRANDS AS IN PENTA #################')
        print('############################################################')
        
        # fix_func = f_j(K)*K^(3/2)
        fix_func = self.get_fix_func(K,intj,make_plot=False)
        
        # create dicionaries of lists
        self.lstar_integrand = defaultdict(list)
        self.mstar_integrand = defaultdict(list)
        self.nstar_integrand = defaultdict(list)
        
        for i in range(self.nefield):
            i1 = i*self.ncmul
            i2 = i1 + self.ncmul
            
            efield = self.efield[i1]
            
            x = self.cmul[i1:i2]
            yl = self.lstar[i1:i2]
            yn = self.nstar[i1:i2]
            ylogm = np.log( self.mstar[i1:i2] )
            
            # quadratic spline as in PENTA. Assuming log_interp = true
            xlog = np.log(x)
            lstar_interp = interp1d(xlog,yl,kind='quadratic',bounds_error=False,fill_value=0.0)
            nstar_interp = interp1d(xlog,yn,kind='quadratic',bounds_error=False,fill_value=0.0)
            logmstar_interp = interp1d(xlog,ylogm,kind='quadratic',bounds_error=False,fill_value=0.0)
            #this function is used to multiply exp(logmstar_interp1d), otherwise exp(0)=1 is taken outside the interpolating region
            filter_logmstar = interp1d(xlog,np.ones_like(xlog),bounds_error=False,fill_value=0.0)
            
            xspline = np.logspace(np.log10(x[0]),np.log10(x[-1]),100)
            
            fig, ax = plt.subplots(1,3,figsize=(17,6))
            ax[0].plot(x,yl,'ob')
            ax[0].plot(xspline, lstar_interp(np.log(xspline)),'red',label='spline')
            ax[0].set_yscale('log')
            ax[0].set_xscale('log')
            ax[0].set_ylabel(r'lstar')
            ax[0].set_xlabel(f'cmul')   
            ax[0].set_title(f'Er/v={efield}')
            ax[0].grid()
            ax[0].legend()
            
            ax[1].plot(x,yn,'ob')
            ax[1].plot(xspline, nstar_interp(np.log(xspline)),'red',label='spline')
            ax[1].set_xscale('log')
            ax[1].set_ylabel(r'nstar')
            ax[1].set_xlabel(f'cmul')   
            ax[1].set_title(f'Er/v={efield}')
            ax[1].grid()   
            #ax[1].legend()
            
            ax[2].plot(x,ylogm,'ob')
            ax[2].plot(xspline, logmstar_interp(np.log(xspline)),'red',label='spline')
            ax[2].set_xscale('log')
            ax[2].set_ylabel('ln(mstar)')
            ax[2].set_xlabel(f'cmul')   
            ax[2].set_title(f'Er/v={efield}')
            ax[2].grid()   
            #ax[1].legend()
            
            plt.tight_layout(pad=2)
            
            # full integrand of l*
            fig,ax = plt.subplots(figsize=(8,6))
            for species in self.plasma_class.list_of_species:
                integrand = lstar_interp(np.log(self.cmul_species[species]))*fix_func
                integral = self.get_integral(integrand,K)
                
                #save integrand
                self.lstar_integrand[species].append( integrand )
                
                #plot
                ax.plot(K,integrand,'o-',label=f'{species}, {integral:.3e}')       
                ax.set_xlabel('K')
                ax.set_ylabel(fr'$f_{intj}(K)~l^*(K)~K^{{3/2}}$')
                ax.set_title(f'Er/v={efield}')
                ax.grid()
            plt.legend()
                      
            # full integrand of n*
            fig,ax = plt.subplots(figsize=(10,6))
            for species in self.plasma_class.list_of_species:
                integrand = nstar_interp(np.log(self.cmul_species[species]))*fix_func
                integral = self.get_integral(integrand,K)
                
                #save integrand
                self.nstar_integrand[species].append( integrand )                
                
                #plot
                ax.plot(K,integrand,'o-',label=f'{species}, {integral:.3e}')       
                ax.set_xlabel('K')
                ax.set_ylabel(f'$f_{intj}(K)~n^*(K)~K^{{3/2}}$')
                ax.set_title(f'Er/v={efield}')
                ax.grid()
            plt.legend()
            
            # full integrand of m*
            fig,ax = plt.subplots(figsize=(10,6))
            for species in self.plasma_class.list_of_species:
                integrand = np.exp(logmstar_interp(np.log(self.cmul_species[species])))*filter_logmstar(np.log(self.cmul_species[species]))*fix_func
                integral = self.get_integral(integrand,K)
                
                #save integrand
                self.mstar_integrand[species].append( integrand )
                
                #plot
                ax.plot(K,integrand,'o-',label=f'{species}, {integral:.3e}')       
                ax.set_xlabel('K')
                ax.set_ylabel(f'$f_{intj}(K)~m^*(K)~K^{{3/2}}$')
                ax.set_title(f'Er/v={efield}')
                ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
                ax.grid()
            plt.legend()
            if make_plot:
                plt.show()
            else:
                plt.close('all')
            
    def get_integral(self,y,x,xmin=None,xmax=None,make_plot=False,plot_title=None):
        # this function computes the trapezoid integral of y=y(x)
        # if xmin or xmax are provided, the limits of the integral are changed
        # an error is raised if xmin or xmax fall outside the domain defined by x
        # if the array x does not contain xmin or xmax exactly, then the closest value is considered
        
        import matplotlib.pyplot as plt
        from scipy.integrate import trapezoid
        
        # Check if xmin and xmax are within the domain of x
        if xmin is not None and (xmin < x.min() or xmin > x.max()):
            raise ValueError(f"xmin ({xmin}) is outside the domain of x.")
        if xmax is not None and (xmax < x.min() or xmax > x.max()):
            raise ValueError(f"xmax ({xmax}) is outside the domain of x.")

        # If xmin or xmax are provided, adjust the limits
        if xmin is not None:
            xmin_index = np.argmin(np.abs(x - xmin))  # Find closest value to xmin in x
        else:
            xmin_index = 0  # If xmin is None, start from the beginning

        if xmax is not None:
            xmax_index = np.argmin(np.abs(x - xmax))  # Find closest value to xmax in x
        else:
            xmax_index = len(x) - 1  # If xmax is None, go to the end

        # Perform the trapezoidal integration
        x_selected = x[xmin_index:xmax_index+1]
        y_selected = y[xmin_index:xmax_index+1]
        integral = trapezoid(y_selected, x_selected)
        
        integral_exact = trapezoid(y,x)
        rel_error = np.abs(integral_exact-integral) / integral_exact

        # Optionally plot the function
        if make_plot:
            plt.rc('font', size=16)
            fig=plt.figure(figsize=(10,8))
            plt.plot(x, y, '.-')
            plt.fill_between(x_selected, y_selected, alpha=0.3, label=f'rel. error={rel_error*100:.1f}%')
            plt.xlabel('x')
            plt.legend()
            if plot_title is not None:
                plt.title(plot_title)
            plt.show()

        return integral
    
    def compute_energy_convolution(self,which_convol, cmin, cmax):
        #which_convol should be 'lstar', 'mstar or nstar'
        
        convol_type = ['lstar','mstar','nstar']
        
        if which_convol not in convol_type:
            print(f'ERROR: which_convol should take one of the following: {convol_type}')
            exit(1)
            
        #check if integrand exist
        if not hasattr(self,which_convol+'_integrand'):
            print('ERROR: integrand does not exist! Need to set it up first!!')
            exit(1)
        else:
            integrand = getattr(self,which_convol+'_integrand')
            
        #check if cmin and cmax are inside the cmul domain
        if cmin>cmax or cmin<np.min(self.cmul) or cmax>np.max(self.cmul):
            print('ERROR: limits of integral not correct. Cannot proceed')
            exit(1)
            
        #loop in field
        for i in range(self.nefield):
            i1 = i*self.ncmul
            #i2 = i1 + self.ncmul
            
            efield = self.efield[i1]
          
            for species in self.plasma_class.list_of_species:
                
                # compute Kmin and Kmax according to cmin and cmax
                # we take the values in self.cmul closest to cmin and cmax
                cmin_index = np.argmin(np.abs(self.cmul_species[species] - cmin))
                cmax_index = np.argmin(np.abs(self.cmul_species[species] - cmax))
            
                Kmin = self.K[np.min([cmin_index,cmax_index])]
                Kmax = self.K[np.max([cmin_index,cmax_index])]
                
                self.get_integral(integrand[species][i],self.K,xmin=Kmin,xmax=Kmax,make_plot=True,plot_title=which_convol+f', {species}, Er/v={efield}')
                       
    def plot_U2_estimate(self):
        
        import matplotlib.pyplot as plt
        
        #check what are the slopes of D11*K^3/2  for the smallest and largest electric field
        j = 0
        i1 = j*self.ncmul
        i2 = i1 + self.ncmul
        D11_Emin = self.D11_star[i1:i2]
        cmul_fit = self.cmul[i1:i2]
        
        j = self.nefield-1
        i1 = j*self.ncmul
        i2 = i1 + self.ncmul
        D11_Emax = self.D11_star[i1:i2]
        
        #take only the last 6 points of cmul for the fit
        cmul_fit = cmul_fit[-6:]
        D11_Emin = D11_Emin[-6:]
        D11_Emax = D11_Emax[-6:]
        
        p_Emin = np.polyfit(np.log10(cmul_fit),np.log10(D11_Emin),1)
        p_Emax = np.polyfit(np.log10(cmul_fit),np.log10(D11_Emax),1)
        
        x_fit = cmul_fit
        Emin_fit = 10**p_Emin[1] * x_fit**p_Emin[0]
        Emax_fit = 10**p_Emax[1] * x_fit**p_Emax[0]
        
        Emin_fit_plot = 10* 10**p_Emin[1] * x_fit**p_Emin[0] 
        Emax_fit_plot = 10* 10**p_Emax[1] * x_fit**p_Emax[0]
        
        print(f'slope Emin={p_Emin[0]}')
        print(f'slope Emax={p_Emax[0]}')
        
        fig=plt.figure(figsize=(11,8))
        ax = fig.add_subplot(111)
        for i in [0,self.nefield-1]:
            i1 = i*self.ncmul
            i2 = i1 + self.ncmul
            # from plot_DKES_coefficients
            [yerr_lower, yerr_upper] = self.compute_yerr(self.L11[i1:i2],self.Lm['L11'][i1:i2],self.Lp['L11'][i1:i2])
            ax.errorbar(self.cmul[i1:i2],self.L11[i1:i2],yerr=[yerr_lower,yerr_upper],fmt='-o',label=rf'$E_s/v$={self.efield[i1]:3.1E}',capsize=5, elinewidth=2, markeredgewidth=2)
            ax.plot(x_fit,Emin_fit_plot,'--r',linewidth=5)
            ax.plot(x_fit,Emax_fit_plot,'--r',linewidth=5)
                
        ax.set_xlabel(r'$\nu/v\,\,[\text{m}^{-1}]$')
        ax.set_ylabel(r'$D_{11}^*K^{3/2}~~[\text{m}^{-1}~\text{T}^{-2}]$')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_title(r'$\left<U^2\right>=$'+f'{1.5*10**p_Emin[1]:.2e}'+' T'+'$^{-2}$')
        ax.text(x_fit[-4],Emin_fit_plot[-2],f'slope={p_Emin[0]:.3f}')
        ax.text(x_fit[-4],Emax_fit_plot[-2],f'slope={p_Emax[0]:.3f}')
        ax.legend(fontsize=12)
        ax.grid()
        plt.show()   
            
    def write_U2_to_file(self,rho,filename=None):
        print(f'\n ########### ASSUMING U2 IS CONSTANT: <U^2>={self.Usq} ##############')
        
        if filename is None:
            filename = 'Utilde2_profile'
            
        with open(filename, 'w') as file:
            # Write the size of rho as the first line
            file.write(f'{rho.size}\n')
            
            # Write the data to the file
            for i in range(rho.size):
                row_data = [rho[i], self.Usq]

                # Write the row to the file, formatted as space-separated values
                file.write(" ".join(map(str, row_data)) + '\n')
                
    def read_PENTA_fluxes_vs_Er(self,plasma_class,folderpath=None):
        #reads the file fluxes_vs_Er creatd by PENTA code
        #folderpath should contain path to file without final '/'
        # if not given, it is assumed we are already inside the folder
        
        import matplotlib.pyplot as plt
        
        print('\n#############################################################')
        print('####################### PLOTTING FLUXES ######################')
        print('########## ASSUMES ALL SPECIES HAVE SAME CHARGE #############')
        print('#############################################################')
        print('#############################################################')
        
        if folderpath is None:
            filename = 'fluxes_vs_Er'
        else:
            filename = folderpath+'/fluxes_vs_Er'
            
        num_ions = plasma_class.num_ion_species
            
        penta = np.loadtxt(filename,skiprows=2)
        
        roa = penta[:,0]
        Er = penta[:,1]
        gamma_e = penta[:,2]
        
        gamma_i = {}
        j = 3
        for ion in plasma_class.ion_species:
            gamma_i[ion] = penta[:,j]
            j+=1
            
        gamma_i_tot = sum(gamma_i[ion] for ion in plasma_class.ion_species)
        
        plt.rc('font', size=18)
        fig=plt.figure(figsize=(11,8))
        ax = fig.add_subplot(111)
        ax.plot(Er,gamma_e,label='$\Gamma_e$')
        ax.plot(Er,gamma_i_tot,label='$\Gamma_{i,tot}$')
        for ion in plasma_class.ion_species:
            ax.plot(Er,gamma_i[ion],'--',label=f'$\Gamma({ion})$')
        ax.set_xlabel(r'Er [V/cm]')
        ax.set_ylabel(r'$\Gamma~~[\text{m}^{-2}\,\text{s}^{-1}]$')
        ax.set_title(f'r/a={roa[0]:.2f}')
        ax.set_yscale('log')
        ax.legend(fontsize=12)
        ax.grid()
        plt.legend()
        #plt.show()
        
        # plt.rc('font', size=18)
        # fig=plt.figure(figsize=(11,8))
        # ax = fig.add_subplot(111)
        # ax.plot(Er,gamma_e-gamma_i_tot,label='$\sum Z_j\Gamma_j$')
        # ax.set_xlabel(r'Er [V/cm]')
        # ax.set_ylabel(r'$\Gamma~~[\text{m}^{-2}\,\text{s}^{-1}]$')
        # ax.set_title(f'r/a={roa[0]:.2f}')
        # #ax.set_yscale('log')
        # ax.legend(fontsize=12)
        # ax.grid()
        # plt.legend()
        plt.show()

            
        
            
        
        
        
        
        
        
        
        
        
        
        
        
    
    # def perp_coll_freq(self,n,m,Z,T,v,species_num,electron_num):
    #     """Calc perpendicular collision frequency of species species_num

	# 	This routine takes a 1D arrays of n,m,Z,T,v which contains the
    #     density, mass, charge number, Temperature and velocity of all species
    #     in a plasma. Then it computes nu_D of species indicated by the 
    #     number species_num, which refers to the index in the arrays of 
    #     the species we are computing the collision frequency
        
    #     This is the collision frequency of the pitch-angle scatering
    #     operator and is defined in 
    #     S. P. Hirshman and D. J. Sigmar, Nucl. Fusion 21, 1079 (1981)

	# 	Parameters
	# 	----------
	# 	n : list
	# 		density in m^-3
    #     m : list
	# 		mass in kg
    #     Z : list
	# 		charge number; if electron, Z=-1
    #     T : list
	# 		temperature in eV
    #     v : list
	# 		velocity in m/s
    #     species num : integer
    #         identifies the species we want the collision freq of: nu_D_num_index
    #     electron_num : integer
    #         identifies the index of the electrons
		
    #     Returns
	# 	----------
	# 	nu_D : real
	# 		Perpendicular collision frequency [s^-1]
	# 	"""
  
    #     from scipy.special import erf
        
    #     EC = 1.602176634E-19 # Electron charge [C]
    #     EPS0 = 8.8541878188E-12 # Vacuum permittivity [F/m]
          
    #     Za = Z[species_num]
    #     ma = m[species_num]
    #     va = v[species_num]
        
    #     ne = n[electron_num]
    #     Te = T[electron_num]
        
    #     #compute loglambda as in PENTA
    #     if(Te>50):
    #         loglambda = 25.3 - 1.15*np.log10(ne/1e6) + 2.3*np.log10(Te)
    #     else:
    #         loglambda = 23.4 - 1.15*np.log10(ne/1e6) + 3.45*np.log10(Te)
        
    #     nu_D = 0.0
    #     for nb,mb,Zb,Tb,vb in zip(n,m,Z,T,v):
            
    #         vthb = np.sqrt(2*EC*Tb/mb)
    #         xb = vb/vthb
            
    #         Hb = (1-0.5/(xb*xb))*erf(xb) + np.exp(-xb*xb)/(xb*np.sqrt(np.pi))
            
    #         num = nb*(EC**4)*Za*Za*Zb*Zb*loglambda*Hb
    #         den = ma*ma*va*va*va*4*np.pi*EPS0*EPS0
            
    #         nu_D = nu_D + num/den
            
    #     return nu_D
            
            
        
        
        
        
        

      

     


  
# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)