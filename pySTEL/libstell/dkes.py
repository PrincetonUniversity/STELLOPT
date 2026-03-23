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
        
        #K array as in PENTA3
        Kmin = 1e-4      #PENTA default value
        Kmax = 10        #PENTA default value
        numKsteps = 1000  #PENTA default value
        #K_PENTA1 = np.linspace(Kmin,Kmax,numKsteps)
        # K in PENTA 3
        K = 10**np.linspace(np.log10(Kmin),np.log10(Kmax),numKsteps)
        self.K = K
        
    def read_DKES_results(self,filename):
        """
        Reads a DKES results file named 'filename'. Assumes that data is regular, i.e.
        cmul and efield form a rectangular meshgrid
        """
        
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
        self.L11  = 0.5*(self.L11m+self.L11p)
        self.L31  = 0.5*(self.L31m+self.L31p)
        self.L33  = 0.5*(self.L33m+self.L33p)
        
        #compute Dij_star, correcting D13 and D33 with the B factors (see J.Lore documentation on PENTA)
        self.D11_star = self.L11
        self.D31_star = self.L31 * np.sqrt(self.Bsq)
        self.D33_star = self.L33 * self.Bsq
        
        # We need to take D13* = +D31* in order to get the correct sign of JBS.B (this has been confirmed by comparing with experiments)
        # It has also been confirmed by running a tokamak case, where the bootstrap must always run in the same direction as the already
        # existing ohmic current (assuming that the density and temperature gradients are negative! See Freidberg)
        # Note that D13_star is used in this class ONLY to create PENTA Dstar input files
        self.D13_star = self.D31_star
        
        # Inside PENTA, D31_star is computed by taking the negative of D13_star, so we should make the same here otherwise when computing
        # bootstrap with self.get_BS_current() the sign is different from that given by PENTA
        self.D31_star = -self.D13_star
        
        # Compute <U2> as in PENTA: by looking at D11_star at maximum nu/v and minimum Er/v (ideally 0)
        self.U2 = 1.5 * self.D11_star[self.ncmul-1] / self.cmul[self.ncmul-1]
        
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
            
    def plot_DKES_coeffs(self,R0=None,iota=None,use_nustar=False):
        """
        Plots the species-independent DKES coefficients L11, L13 and L33 
        To use nu_star at the x-axis instead of cmul, R0 and iota need to be provided:
        nu_star = (nu/v) * (R0/iota) = cmul * (R0/iota)
        """
        # 
        import matplotlib.pyplot as pyplot
        
        var_names = {
                'L11': '$D_{11}^*~~[m^{-1}~T^{-2}]$',
                'L31': '$D_{31}^*$',
                'L33': '$D_{33}^*~~[m~T^2]$',
                } 
        
        for plot_var in ['L11', 'L31', 'L33']:
            
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
                if(use_nustar and (R0 is not None) and (iota is not None)):
                    nu_star = self.cmul[i1:i2] * (R0/iota)
                    ax.set_xlabel(r'$\nu^*$')
                    x_axis = nu_star
                else:
                    x_axis = self.cmul[i1:i2]
                    ax.set_xlabel(r'$\nu/v\,\,[\text{m}^{-1}]$')
                ax.errorbar(x_axis,yplot[i1:i2],yerr=[yerr_lower,yerr_upper],fmt='-o',label=rf'$E_s/v$={self.efield[i1]:3.1E}',capsize=5, elinewidth=2, markeredgewidth=2)
            ax.set_ylabel(f'{var_names[plot_var]}')
            ax.set_xscale('log')
            if(plot_var=='L11' or plot_var=='L33'):
                ax.set_yscale('log')
            if(plot_var=='L31'):
                ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
            ax.set_title(f'r/a={self.roa:.2f}')
            ax.legend(fontsize=12)
            ax.grid()
    
        pyplot.show()
        
    def plot_coeffs(self,which_coeff):
        """
        Plots coefficients a la carte. which_coeff can be:
        D11_star, D31_star, D31_over_D33_corrected or D33_star
        """
        import matplotlib.pyplot as pyplot
        
        if(which_coeff == 'D11_star'):
            yplot = self.D11_star
            var_name = r'$D_{11}^*~~[m^{-1}~T^{-2}]$'
            yscale_log = True
        elif(which_coeff == 'D31_star'):
            yplot = self.D31_star
            yscale_log = False
            var_name = r'$D_{31}^*$'
        elif(which_coeff == 'D31_over_D33_corrected'):
            D33_corrected = (2./3.)*self.Bsq/self.cmul - self.D33_star
            yplot = self.D31_star / D33_corrected
            var_name = r'$D_{31}^*\,\,/\,\,[(2/3)(\nu/v)^{-1}-D_{33}^*]$'
            yscale_log = False
        elif(which_coeff == 'D33_star'):
            yplot = self.D33_star
            var_name = r'$D_{33}^*$'
            yscale_log = False
        elif(which_coeff == 'D11_plus_D31sq_over_D33_corrected'):
            D33_corrected = (2./3.)*self.Bsq/self.cmul - self.D33_star
            yplot = self.D11_star + self.D31_star**2 / D33_corrected
            var_name = r'$D_{11}^*+[D_{31}^*]^2\,\,/\,\,[(2/3)(\nu/v)^{-1}-D_{33}^*]$'
            yscale_log = True
        else:
            print('Coeff not found...')
            exit(1)
                        
        px = 1/pyplot.rcParams['figure.dpi']
        pyplot.rc('font', size=20)
        #pyplot.rc('legend', fontsize=24)
        fig=pyplot.figure(figsize=(1024*px,768*px))
        ax = fig.add_subplot(111)
        for i in range(self.nefield):
            i1 = i*self.ncmul
            i2 = i1 + self.ncmul
            # plot without error bar
            ax.plot(self.cmul[i1:i2],yplot[i1:i2],marker='+',label=rf'$E_s/v$={self.efield[i1]:3.1E}',linewidth=4,markersize=18)
            # plot with error bar
            # [yerr_lower, yerr_upper] = self.compute_yerr(yplot[i1:i2],self.Lm[plot_var][i1:i2],self.Lp[plot_var][i1:i2])
            x_axis = self.cmul[i1:i2]
            ax.set_xlabel(r'$\nu/v\,\,[\text{m}^{-1}]$')
        ax.set_ylabel(var_name)
        ax.set_xscale('log')
        if(yscale_log): ax.set_yscale('log')
        ax.set_title(f'r/a={self.roa:.2f}')
        ax.legend(fontsize=12)
        ax.grid()
        pyplot.show()
    
    def plot_L33_test(self):
        
        import matplotlib.pyplot as pyplot
        
        L33 = 0.5*(self.L33m+self.L33p)
        
        L33_Spitzer = (2/3)/self.cmul - L33
        print(L33_Spitzer)
        
        #check what are the slopes of D11*K^3/2  for the smallest and largest electric field
        j = 0
        i1 = j*self.ncmul
        i2 = i1 + self.ncmul
        L33_fit = L33[i1:i2]
        cmul_fit = self.cmul[i1:i2]
        
        #fit to L33p
        #take only the last 6 points of cmul for the fit
        cmul_fit = cmul_fit[-6:]
        L33_fit = L33_fit[-6:]
        L33_Spitzer_fit = L33_Spitzer[-6:]
        
        p_L33 = np.polyfit(np.log10(cmul_fit),np.log10(L33_fit),1)
        p_L33_Spitzer = np.polyfit(np.log10(cmul_fit),np.log10(L33_Spitzer_fit),1)
        
        x_fit = cmul_fit
        #Emin_fit = 10**p_L33[1] * x_fit**p_L33[0]
        
        L33_fit_plot = 10* 10**p_L33[1] * x_fit**p_L33[0] 
        L33_Spitzer_fit_plot = 10* 10**p_L33_Spitzer[1] * x_fit**p_L33_Spitzer[0] 
        
        print(f'slope L33={p_L33[0]}')
        print(f'slope L33_Spitzer={p_L33_Spitzer[0]}')
                
        px = 1/pyplot.rcParams['figure.dpi']
        pyplot.rc('font', size=20)
        #pyplot.rc('legend', fontsize=24)
        fig=pyplot.figure(figsize=(1024*px,768*px))
        ax = fig.add_subplot(111)
        for i in range(self.nefield):
            i1 = i*self.ncmul
            i2 = i1 + self.ncmul
            ax.plot(self.cmul[i1:i2],self.L33m[i1:i2],marker='.',label=rf'L33m',linewidth=4,markersize=18)
            ax.plot(self.cmul[i1:i2],self.L33p[i1:i2],marker='.',label=rf'L33p',linewidth=4,markersize=18)
            ax.plot(self.cmul[i1:i2],L33_Spitzer[i1:i2],marker='.',label=rf'L33_with_Spitzer',linewidth=4,markersize=18)
            ax.plot(x_fit,L33_fit_plot,'--r',linewidth=5)
            ax.plot(x_fit,L33_Spitzer_fit_plot,'--r',linewidth=5)
                
        ax.set_xlabel(r'$\nu/v\,\,[\text{m}^{-1}]$')
        ax.set_ylabel(f'L33')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_title(f'r/a={self.roa:.2f}')
        ax.legend(fontsize=12)
        ax.grid()
        pyplot.show()
        
    def plot_Er_resonance(self,Er_v_resonance=None):
        """ Make plot of D11* as function of Er/v for each cmul"""
        import matplotlib.pyplot as pyplot
        
        px = 1/pyplot.rcParams['figure.dpi']
        pyplot.rc('font', size=20)
        fig=pyplot.figure(figsize=(1024*px,768*px))
        ax = fig.add_subplot(111)
        
        yplot = self.D11_star
        
        for i1 in range(self.ncmul):
            ax.plot(self.efield[i1::self.ncmul],yplot[i1::self.ncmul],marker='+',label=rf'$\nu/v$={self.cmul[i1]:3.1E}',linewidth=4,markersize=18)
            ax.set_xlabel(r'$E_r/v$')
            ax.set_ylabel(r'$D_{11}^*$')
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.legend(fontsize=12)
            ax.grid(True)
            ax.set_title(f'r/a={self.roa:.2f}')
        if Er_v_resonance is not None:
                ax.vlines(Er_v_resonance,np.min(yplot),np.max(yplot),'r','dashed',linewidth=3,label='Er/v res')
        pyplot.legend()
        pyplot.show()
        
    def compute_yerr(self,y,ym,yp):
        """ 
        Computes the distance between the central value 'y' and its upper and lower limits
        A priori we don't know which one of these is the upper and lower since they change
        depending on their sign
        The output of this function is used for the argument of the maplotlib errorbar function
        """

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
    
    def get_VMEC_quantities(self,wout_file):
        import netCDF4 as nc
        
        try:
            dataset = nc.Dataset(wout_file, 'r')
            Bsq = dataset.variables['bdotb'][:]
            phi = dataset.variables['phi'][:]
            bmnc = dataset.variables['bmnc'][:]
            self.aspect_ratio = dataset.variables['aspect'][:]
            self.aminor = dataset.variables['Aminor_p'][:]
            dataset.close()
        except:
            print('\nERROR: Could not read wout file')
            sys.exit(0)
        
        self.Bsq = Bsq[self.surface-1]
        self.roa = np.sqrt(phi[self.surface-1]/phi[-1])
        # <|B|> = B00
        self.B0 = bmnc[self.surface-1,0]
        
        print(f'r/a = sqrt(PHI/PHIEDGE) = {self.roa}')
        print(f'Bsq = {self.Bsq}')
        print(f'B0 = {self.B0}')
        print(f'R/a = {self.aspect_ratio}')
    
    def write_PENTA3_coeffs_to_files(self, where_to, cmul_list=None, efield_list=None):
        """ 
        Creates files to be read by PENTA3: 'D11_star_##' , 'D13_star_##', 'D33_star_##'
        where_to specifies ssaving directory (path should not have final /)
        
        If cmul_list and/or efield_list are None, then it uses ALL available values
        If not None, then looks for the closest available
        This is used to reduce the parameter space of DKES to make it faster (at the cost of loosing accuracy, of course)
        It's also used to make convergence tests in cmul/efield
        """
        
        #checks that data is regular, i.e.: for each efield there are always the same cmul
        # if data is regular the following arrays are computed:
        # self.cmul_regular
        # self.efield_regular
        self.check_data_is_regular(self.cmul,self.efield)
        
        # If no values are provided, use the full range of regular values
        if cmul_list is None:
            cmul_list = self.cmul_regular
        else:
            cmul_list = np.array(cmul_list)
            #check that cmul_list is inside the domain of self.cmul_regular
            if( np.min(cmul_list) < np.min(self.cmul_regular) or np.max(cmul_list) > np.max(self.cmul_regular)):
                print('ERROR!! the given cmul values are outisde the domain of cmul....')
            
            
        if efield_list is None:
            efield_list = self.efield_regular
        else:
            efield_list = np.array(efield_list)
            #check that efield_list is inside the domain of self.efield_regular
            if( np.min(efield_list) < np.min(self.efield_regular) or np.max(efield_list) > np.max(self.efield_regular)):
                print('ERROR!! the given efield values are outisde the domain of cmul....')
                
        
        for var in ['D11_star', 'D13_star', 'D33_star']:
            
            filename = where_to + '/' + var + '_surface_' + f'{self.surface}'
            y = getattr(self,var)
            
            # Get the closest indices for cmul and efield (if cmul_list and/or efield_list are not provided, this ends up with the same array)
            cmul_indices = self.get_closest_indices(cmul_list, self.cmul_regular)
            efield_indices = self.get_closest_indices(efield_list, self.efield_regular)
            
            # Get the new cmul_regular and efield_regular arrays using the closest indices
            new_cmul_regular = np.array([self.cmul_regular[i] for i in cmul_indices])
            new_efield_regular = np.array([self.efield_regular[i] for i in efield_indices])
            
            ncmul = len(new_cmul_regular)
            nefield = len(new_efield_regular)
            
            # Create an empty list to store the new y values
            new_y = []
            
            # Loop over the closest indices and extract the corresponding D11_star values
            for e_idx in efield_indices:
                for c_idx in cmul_indices:
                    # Compute the position in the y array
                    index = e_idx * len(self.cmul_regular) + c_idx
                    new_y.append(y[index])

            combined = np.concatenate((new_cmul_regular,new_efield_regular,new_y))
            
            #create file
            with open(filename, 'w') as file:
                file.write(f'{ncmul} {nefield}\n')
                for value in combined:
                    file.write(f'{value:.10e}\n')           
                
    def check_data_is_regular(self,cmul,efield):
        """ Checks that data is regular, i.e. for each efield there are always the same cmul.
        If data is regular the following arrays are computed:
        self.cmul_regular
        self.efield_regular
        """
        
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
                
    def get_closest_indices(self,values, regular_values):
        """For each value in 'values', find the index of the closest value in 'regular_values'."""
        indices = np.abs(np.array(regular_values)[:, np.newaxis] - values).argmin(axis=0)
        return indices
                
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
        """ 
        Makes a plot of Er/v as a function of K for all species
        Er_V_cm is in V/cm and it can be an array or a double
        """
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
                
                ax.plot(K,Er_over_v ,'.',label=f'{species}',markersize=2.5) 
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
    
    def get_PENTA3_energy_convolution(self,which_coeff,which_species,Er,plasma_class,K_exp=0,jval=0,log_interp_coeff=False,make_plot=True):
        """
        Plots the convolution integrand as in PENTA3:
        
        convolution = 2*na/sqrt(pi) * integral( dk * integrand )
        integrand = k**K_exp * sqrt(k) * exp(-k) * assoc_laguerre_{jval}^(3/2) * which_coeff
        
        which coeff can be:
        D11_star, D31_star, D33_star, D31_over_D33_corrected, D31sq_over_D33_corrected, LHS_SN_flow_eq
        
        Er must be provided in V/cm
        """
        import matplotlib.pyplot as plt
        from scipy.integrate import trapezoid
        from scipy.special import assoc_laguerre
        
        if(which_coeff == 'D11_star'):
            coeff = self.D11_star
        elif(which_coeff == 'D31_star'):
            coeff = self.D31_star
        elif(which_coeff == 'D33_star'):
            coeff = self.D33_star
        elif(which_coeff == 'D31_over_D33_corrected'):
            D33_corrected = (2./3.)*self.Bsq/self.cmul - self.D33_star
            coeff = self.D31_star / D33_corrected
        elif(which_coeff == 'D31sq_over_D33_corrected'):
            D33_corrected = (2./3.)*self.Bsq/self.cmul - self.D33_star
            coeff = self.D31_star**2 / D33_corrected
        elif(which_coeff == 'LHS_SN_flow_eq'):
            D33_corrected = (2./3.)*self.Bsq/self.cmul - self.D33_star
            coeff = (2/3)*self.Bsq/D33_corrected - self.cmul
        else:
            print('Coeff not found...')
            exit(1)
            
        self.check_data_is_regular(self.cmul,self.efield)
        
        cmul_log = np.log10( np.unique(self.cmul) )
        efield_log = np.log10( np.unique(self.efield) )
        
        # convert coeff to 2D array
        coeff_2d = coeff.reshape(len(efield_log), len(cmul_log)).transpose()
            
        #get thermal speed
        vth = plasma_class.get_thermal_speed(which_species,self.roa)

        #get cmul for self.K
        cmul_species = []
        Erv_species = []
        for k in self.K:
            vparticle = vth * np.sqrt(k)
            nu = plasma_class.get_collisionality(which_species,self.roa,vparticle)
            cmul_species.append( nu / vparticle )
            Erv_species.append(np.abs(Er*100/vparticle))
                         
        log_cmul_K = np.log10(cmul_species)
        log_efield_K = np.log10( np.abs(Er)*100/(vth*np.sqrt(self.K)) )
        
        # handle points that are out of range as in PENTA3
        log_cmul_K_non_clipped = log_cmul_K
        log_efield_K_non_clipped = log_efield_K
        log_cmul_K = np.clip(log_cmul_K, np.min(cmul_log), np.max(cmul_log))
        log_efield_K = np.clip(log_efield_K, np.min(efield_log), np.max(efield_log))    
        # get indexes where values were clipped 
        idx_clipped = (log_cmul_K_non_clipped!=log_cmul_K) + (log_efield_K_non_clipped!=log_efield_K)
        
        integrand = self.get_interpolated_coeff(which_coeff,cmul=cmul_species,efield=Erv_species,log_interp_coeff=log_interp_coeff)
            
        integrand = integrand * self.K**K_exp * np.sqrt(self.K) * np.exp(-self.K) * assoc_laguerre(self.K, jval, k=1.5)
        
        na = plasma_class.get_density(which_species,self.roa)
        norm =  2 * na / np.sqrt(np.pi)
        convolution = norm*trapezoid(integrand,self.K)
        
        print(f'Energy convolution = {convolution:.1e}')
        
        # Optionally plots the coefficient as function of cmul and the integrand on a 2nd axis
        # this allows to understand what are the most important points
        if make_plot:
            plt.rc('font', size=16)
            fig, ax1 = plt.subplots(figsize=(10,8))
            
            cmul_species = np.array(cmul_species)
            Erv_species  = np.array(Erv_species)
            integrand = np.array(integrand)
            
            # First plot the coefficient
            for ie in range(0, self.nefield):
                ax1.plot(np.unique(self.cmul), coeff_2d[:, ie], '.-', label=f'$E_r/v={np.unique(self.efield)[ie]:3.1E}$')

            # Create a second y-axis on the right
            ax2 = ax1.twinx()  # Create another axis that shares the same x-axis
            ax2.set_yscale('linear') 

            # Now plot the integrand
            ax2.plot(cmul_species[~idx_clipped], norm*integrand[~idx_clipped],'.-',color='black')
            ax2.plot(cmul_species[idx_clipped], norm*integrand[idx_clipped],'.-',color='red')
            ax2.fill_between(cmul_species, norm*integrand, alpha=0.3, label=f'|Er|={np.abs(Er):.2f} V/cm')  
            #
            #
            idx = [np.argmin(np.abs(cmul_species-cmul)) for cmul in np.unique(self.cmul)]
            ax2.plot(np.unique(self.cmul),norm*integrand[idx],'x',markersize=10,markeredgewidth=5)        

            if(log_interp_coeff): ax1.set_yscale('log')
            ax1.set_xscale('log')
            ax2.set_xscale('log')
            ax1.set_ylabel(which_coeff)
            ax2.set_ylabel(f'{which_species} ||{which_coeff} K^{K_exp} L_{jval}||')
            ax1.grid()
            ax1.set_xlabel(r'$\nu/v$')
            ax2.legend(loc='upper left')
            ax1.set_title(f'r/a={self.roa:.2f}')
            ax1.legend()
            #plt.show()
            
            _, ax4 = plt.figure(figsize=(10,8)), plt.gca()
            ax4.plot(Erv_species[~idx_clipped], norm*integrand[~idx_clipped],'.-',color='black')
            ax4.plot(Erv_species[idx_clipped], norm*integrand[idx_clipped],'.-',color='red')
            #
            idx = [np.argmin(np.abs(Erv_species-er)) for er in np.unique(self.efield)]
            ax4.plot(np.unique(self.efield),norm*integrand[idx],'x',markersize=10,markeredgewidth=5)
            ax4.set_xlabel(r'$Er/v$')
            ax4.grid()
            ax4.set_xscale('log')
            ax4.set_xlim(np.min(Erv_species),np.max(Erv_species))
            # plt.show()

            #now make a plot that shows how many integration points are out of range
            fig, ax3 = plt.figure(figsize=(10,8)), plt.gca()
            ax3.plot(self.cmul,self.efield,'x',label='DKES data')
            ax3.plot(10**log_cmul_K[idx_clipped],10**log_efield_K[idx_clipped],'.',color='r',label='integration points (out of range)')
            ax3.plot(10**log_cmul_K[~idx_clipped],10**log_efield_K[~idx_clipped],'.',color='k',label='integration points')
            ax3.set_yscale('log')
            ax3.set_xscale('log')
            ax3.set_xlabel(r'$\nu/v$')
            ax3.set_ylabel(r'$E_r/v$')
            ax3.set_title(f'{which_species},  |Er|={np.abs(Er)} V/cm')
            plt.legend()
            plt.show()
            
            # get the maxima...
            y = np.abs(norm*integrand)
            x1 = cmul_species
            x2 = Erv_species
            
            idx = np.argsort(y)[::-1]
            y_sorted = y[idx]
            x1_sorted = x1[idx]
            x2_sorted = x2[idx]
            
            def sci_round(val):
                signif=1
                if val == 0:
                    return 0.0
                exp = int(np.floor(np.log10(abs(val))))
                mant = round(val / 10**exp, signif)
                return mant * 10**exp   # “aEb” form
            
            x1_rounded = np.array([sci_round(v) for v in x1_sorted])
            x2_rounded = np.array([sci_round(v) for v in x2_sorted])

            for x1,x2 in zip(x1_rounded,x2_rounded):
                print(f'{x1:.1e}  ',f'  {x2:.1e}')
             
        return convolution
    
    def get_interpolated_coeff(self,which_coeff,cmul,efield,log_interp_coeff=False):
        """ Interpolates which_coeff at nu/v=cmul and Er/v=efield.
            cmul and efield can be floats or arrays with same size.
            
            which_coeff can be:
            D11_star, D31_star, D33_star, D31_over_D33_corrected, D31sq_over_D33_corrected, LHS_SN_flow_eq,
            capped_fluxes_coefficient
            
            If log_interp_coeff is True, then the log of the coefficient is used when producing the spline
            (this is particularly helpful for coeffs that vary by many orders of magnitude such as D11_star)
            The output is always of the quantity and not of its log
        """
        from scipy.interpolate import RectBivariateSpline
        
        if(which_coeff == 'D11_star'):
            coeff = self.D11_star
        elif(which_coeff == 'D31_star'):
            coeff = self.D31_star
        elif(which_coeff == 'D33_star'):
            coeff = self.D33_star
        elif(which_coeff == 'D31_over_D33_corrected'):
            D33_corrected = (2./3.)*self.Bsq/self.cmul - self.D33_star
            coeff = self.D31_star / D33_corrected
        elif(which_coeff == 'D31sq_over_D33_corrected'):
            D33_corrected = (2./3.)*self.Bsq/self.cmul - self.D33_star
            coeff = self.D31_star**2 / D33_corrected
        elif(which_coeff == 'LHS_SN_flow_eq'):
            D33_corrected = (2./3.)*self.Bsq/self.cmul - self.D33_star
            coeff = (2/3)*self.Bsq/D33_corrected - self.cmul
        elif(which_coeff == 'capped_fluxes_coefficient'):
            D33_corrected = (2./3.)*self.Bsq/self.cmul - self.D33_star
            coeff = self.D11_star - (2/3)*self.U2*self.cmul + self.D31_star**2 / D33_corrected
            coeff = np.clip(coeff,a_min=0.0,a_max=None)
        else:
            print('Coeff not found...')
            exit(1)
            
        # check that there are no Er/v=0
        if np.any(self.efield == 0):
            raise ValueError("Error: 'efield' array contains zero values, which is not allowed.")
        
        cmul_log = np.log10( np.unique(self.cmul) )
        efield_log = np.log10( np.unique(self.efield) )
        
        emin = np.min(efield_log)
        emax = np.max(efield_log)
        
        # Handle out of range as in PENTA
        cmul   = np.clip(cmul, np.min(np.unique(self.cmul)), np.max(np.unique(self.cmul)))
        efield = np.clip(efield, np.min(np.unique(self.efield)), np.max(np.unique(self.efield)))
        
        enrm_log_database = (efield_log-emin)/(emax-emin)
        enrm_log_in       = (np.log10(efield)-emin)/(emax-emin)
        
        # convert coeff to 2D array
        coeff_2d = coeff.reshape(len(efield_log), len(cmul_log)).transpose()

        # spline interpolate the log of the coeff
        # if the log is not taken, then for coefficients that span many orders of magnitude (as D11star)
        # it will give bad results...
        # k=1 corresponds to kord=2 in PENTA spline functions
        if(log_interp_coeff):
            interp_func = RectBivariateSpline(cmul_log, enrm_log_database, np.log10(coeff_2d), kx=1, ky=1)#kx=2, ky=2)
            out = 10**interp_func(np.log10(cmul),enrm_log_in,grid=False)
        else:
            interp_func = RectBivariateSpline(cmul_log, enrm_log_database, coeff_2d,  kx=1, ky=1)#kx=2, ky=2)
            out = interp_func(np.log10(cmul),enrm_log_in,grid=False)            
        return out
    
    def get_auxiliary_integrand(self,which_species,plasma_class,K_exp=0,jval=0,Er=None,make_plot=True):
        """
        Plots the "axuiliary" function of PENTA integrand:
        
        auxiliary_integrand = k**K_exp * sqrt(k) * exp(-k) * assoc_laguerre_{jval}^(3/2)
        
        which_species and a plasma_class must be provided in order to tranform K into nu/v. 
        If Er provided [in V/cm], then same plot is produced as function of Er/v
        
        Note that this function can be called without having read a dkes results file beforehand!
        """
        import matplotlib.pyplot as plt
        from scipy.special import assoc_laguerre
            
        #get thermal speed
        vth = plasma_class.get_thermal_speed(which_species,self.roa)

        #get cmul for self.K
        cmul_species = []
        Erv_species = []
        for k in self.K:
            vparticle = vth * np.sqrt(k)
            nu = plasma_class.get_collisionality(which_species,self.roa,vparticle)
            cmul_species.append( nu / vparticle )
            if(Er is not None):
                Erv_species.append(np.abs(Er*100/vparticle))
            
        auxiliary_integrand = self.K**K_exp * np.sqrt(self.K) * np.exp(-self.K) * assoc_laguerre(self.K, jval, k=1.5)
        
        plt.rc('font', size=16)
        _, ax = plt.subplots(figsize=(10,8))
        
        cmul_species = np.array(cmul_species)
        if(Er is not None):
            Erv_species  = np.array(Erv_species)

        if(make_plot):
            ax.plot(cmul_species, auxiliary_integrand,'.-',color='black')
            ax.set_title(f'r/a={self.roa:.2f}')
            ax.set_xscale('log')
            ax.set_ylabel(f'{which_species}: sqrt(K) e^(-K) K^{K_exp} L_{jval}')
            ax.grid()
            ax.set_xlabel(r'$\nu/v$')
            # ax.legend()
        
        if(Er is not None and make_plot):   
            _, ax = plt.subplots(figsize=(10,8))
            ax.plot(Erv_species, auxiliary_integrand,'.-',color='black',label=f'Er={Er:.1f} V/cm')
            ax.set_ylabel(f'{which_species}: sqrt(K) e^(-K) K^{K_exp} L_{jval}')
            ax.set_xlabel(r'$Er/v$')
            ax.set_title(f'r/a={self.roa:.2f}')
            ax.grid()
            ax.set_xscale('log')
            ax.legend()
        plt.show()
        
        if(Er is None):
            return cmul_species,auxiliary_integrand
        else:
            return cmul_species,Erv_species,auxiliary_integrand
        
    def get_BS_current(self,Er_Vcm,plasma_class,Smax,inspect=False,output=None):
        """ This function computes the SN Bootstrap Current <JBS.b> as PENTA3 does using the Sugama-Nishimura method. This amounts to 
        solve Eq. (3.2.1.1) in PENTA documentation by Jeremy Lore. The equation has one typo though: they are missing the elementary 
        charge multiplying Ta (since it was defined to be in eV)
        
        Note that the equation gives <u_parallel x B>/<B^2> for each species. The parallel BS current of each species is then:
        <JBS.b>_species = n[species] * q[species] * np.sqrt(self.Bsq) * <u_parallel x B>/<B^2>[species]
        
        And the total BS current is simply:
        <JBS.b> = sum_species ( <JBS.b>_species )
        
        !! Er is given in V/cm !!
        
        If inspect is True, then will plot inv(flow_mat)*A1 and inv(flow_mat)*A2 point-wise, to check which terms
        are contributing the most to the BS current
        
        By default, this function returns an array with <JBS.b>_species of each species
        BUT, if output=='flows', then this function returns <u_parallel x B>/<B^2> for each species and j=0,...,Smax
        """
        from scipy.special import assoc_laguerre
        
        Er = Er_Vcm
        # Check Er!=0
        if(np.abs(Er) < 1E-6):
            raise ValueError('Please Choose Er != 0, cause this brings problems when taking the log of Er')
        
        ############################################################################################################################
        ####################### Setup Convolutions in the SN Flow Equation (Eq. 3.2.1.1 in PENTA's notes)  #########################
        ############################################################################################################################
        vth = {}
        n = {}
        T = {}
        q = {}
        m = {}
        LHS_conv    = {species: np.zeros((Smax+1,Smax+1)) for species in plasma_class.list_of_species}
        RHS_A1_conv = {species: np.zeros((Smax+1)) for species in plasma_class.list_of_species}
        RHS_A2_conv = {species: np.zeros((Smax+1)) for species in plasma_class.list_of_species}
        #
        for species in plasma_class.list_of_species:
            vth[species] = plasma_class.get_thermal_speed(species,self.roa)
            n[species] = plasma_class.get_density(species,self.roa)
            T[species] = plasma_class.get_temperature(species,self.roa)
            q[species] = plasma_class.charge[species]
            m[species] = plasma_class.mass[species]
            
            #get cmul for self.K
            cmul_K = []
            for k in self.K:
                vparticle = vth[species] * np.sqrt(k)
                nu = plasma_class.get_collisionality(species,self.roa,vparticle)
                cmul_K.append( nu / vparticle )
                            
            cmul_K = np.array(cmul_K)
            efield_K = np.abs(Er)*100/(vth[species]*np.sqrt(self.K))
            
            D31_over_D33_corrected = self.get_interpolated_coeff('D31_over_D33_corrected',cmul_K,efield_K,log_interp_coeff=False)
            # compute: 2/3B^2/D33 - cmul
            fact_LHS_SN_flow = self.get_interpolated_coeff('LHS_SN_flow_eq',cmul=cmul_K,efield=efield_K,log_interp_coeff=False)
            
            #### RHS convolutions ####
            for jval in range(Smax+1):
                
                # A1
                integrand = D31_over_D33_corrected * self.K**1.5 * np.sqrt(self.K) * np.exp(-self.K) * assoc_laguerre(self.K, jval, k=1.5)
                integrand = integrand * n[species] * 2 / np.sqrt(np.pi)
                RHS_A1_conv[species][jval] =  np.trapz(integrand,x=self.K)
                
                # A2
                integrand = D31_over_D33_corrected * self.K**2.5 * np.sqrt(self.K) * np.exp(-self.K) * assoc_laguerre(self.K, jval, k=1.5)
                integrand = integrand * n[species] * 2 / np.sqrt(np.pi)
                RHS_A2_conv[species][jval] =  np.trapz(integrand,x=self.K)
                
            #### LHS convolution ####
            for jval in range(Smax+1):
                for kval in range(Smax+1):
                    integrand = fact_LHS_SN_flow * assoc_laguerre(self.K, kval, k=1.5)
                    integrand = integrand * self.K**1.5 * np.sqrt(self.K) * np.exp(-self.K) * assoc_laguerre(self.K, jval, k=1.5)
                    LHS_conv[species][jval,kval] = np.trapz(integrand,x=self.K) * n[species] * 2 / np.sqrt(np.pi)          
                               
        ############################################################################################################################
        ###############################################  Compute lmat ##############################################################
        ############################################################################################################################
        # Calculate Coulomb logarithm as in PENTA
        if ( T['electrons'] > 50 ):
            loglambda = 25.3 - 1.15*np.log10(n['electrons']/1.E6) + 2.3*np.log10(T['electrons'])
        else:
            loglambda = 23.4 - 1.15*np.log10(n['electrons']/1.E6) + 3.45*np.log10(T['electrons'])
        #
        lmat = define_friction_coeffs(masses=np.fromiter(m.values(), dtype=float), 
                                     charges=np.fromiter(q.values(), dtype=float), 
                                     v_ths = np.fromiter(vth.values(), dtype=float), 
                                     Temps = np.fromiter(T.values(), dtype=float), 
                                     dens  = np.fromiter(n.values(), dtype=float),
                                     loglambda = loglambda, 
                                     num_species=len(plasma_class.list_of_species), 
                                     Smax=Smax)
        
        ############################################################################################################################
        ###########################################  Assemble flow matrix ##########################################################
        ############################################################################################################################
        Nspecies = len(plasma_class.list_of_species)
        flow_mat = np.zeros(((Smax+1)*Nspecies,(Smax+1)*Nspecies))
        
        for ispec1,species1 in enumerate(plasma_class.list_of_species):
            fact1 = q[species1] / (EC*T[species1])
            fact2 = 1.5*q[species1] / (m[species1]*vth[species1]*EC*T[species1])
            
            for jval in range(Smax+1):
                
                for kval in range(Smax+1):
                    # ind1_LHS1 = (ispec1-1)*(Smax+1)
                    ind1_LHS1 = (ispec1)*(Smax+1)
                    flow_mat[ind1_LHS1+jval,ind1_LHS1+kval] = fact1*LHS_conv[species1][jval,kval]
                    
                    for ispec2,species2 in enumerate(plasma_class.list_of_species):
                        # ind1_LHS2 = ( ispec1 - 1 ) * ( Smax + 1 ) + jval +1
                        ind1_LHS2 = ( ispec1 ) * ( Smax + 1 ) + jval
                        # ind2_LHS2 = ( ispec2 - 1 ) * ( Smax + 1 ) + kval +1 
                        ind2_LHS2 = ( ispec2 ) * ( Smax + 1 ) + kval
                        flow_mat[ind1_LHS2,ind2_LHS2] += -fact2*lmat[ind1_LHS2,ind2_LHS2]
        
        ############################################################################################################################
        ##########################################  Compute A1 and A2 for each species #############################################
        ############################################################################################################################
        A1 = {}
        A2 = {}
        for species in plasma_class.list_of_species:
            na_prime_r = plasma_class.get_density_der(species,self.roa) / self.aminor
            Ta_prime_r = plasma_class.get_temperature_der(species,self.roa) / self.aminor
            #
            A1[species] = na_prime_r/n[species] - 1.5*Ta_prime_r/(T[species]) - q[species]*Er*100/(EC*T[species])
            A2[species] = Ta_prime_r/(T[species])
            
        ############################################################################################################################
        ###############################################  Assemble RHS vector #######################################################
        ############################################################################################################################
        RHS = []; RHS1 = []; RHS2 = []
        for species in plasma_class.list_of_species:
            for jval in range(Smax+1):
                rhs = -RHS_A1_conv[species][jval]*A1[species] - RHS_A2_conv[species][jval]*A2[species]
                RHS.append( rhs )
                if(inspect):
                    RHS1.append(-RHS_A1_conv[species][jval]*A1[species])
                    RHS2.append(-RHS_A2_conv[species][jval]*A2[species])
        RHS = np.array(RHS); RHS1 = np.array(RHS1); RHS2 = np.array(RHS2)
        
        ############################################################################################################################
        ###############################################  Solve Linear System #######################################################
        ############################################################################################################################
        
        uB_over_B2 = np.linalg.solve(flow_mat,RHS)
        if output=='flows':
            return uB_over_B2
        
        uB_over_B2_0 = uB_over_B2[0::Smax+1]
        
        JBS_species = np.zeros(Nspecies)
        for ispecies,species in enumerate(plasma_class.list_of_species):
            JBS_species[ispecies] = n[species]*q[species]*np.sqrt(self.Bsq)*uB_over_B2_0[ispecies]
            
        JBS_tot = np.sum(JBS_species)
        print(f'JBS = {JBS_tot/1E3:.1f} kA/m2')
        print('JBS_species = ',JBS_species)
        
        if(inspect):
            import matplotlib.pyplot as plt
            
            flow_mat_inv = np.linalg.inv(flow_mat)
            flow_mat_inv_point_A1 = flow_mat_inv * RHS1 #[:,None].T
            flow_mat_inv_point_A2 = flow_mat_inv * RHS2 #[:,None].T
            
            uB_over_B2_0_A1 = flow_mat_inv_point_A1[0::Smax+1,:]
            uB_over_B2_0_A2 = flow_mat_inv_point_A2[0::Smax+1,:]

            x_axis = []
            for species in plasma_class.list_of_species:
                for jval in range(Smax+1):
                    x = f'jval={jval:.0f},{species[0]}'
                    x_axis.append(x)
                    
            # nqB = []
            # for ispecies,species in enumerate(plasma_class.list_of_species):
            #     nqB.append( n[species]*q[species]*np.sqrt(self.Bsq) )

            _, ax = plt.subplots(figsize=(11,8))
            ax.plot(uB_over_B2_0_A1.T,'.-',label=plasma_class.list_of_species)
            ax.set_xticks(range(len(x_axis)))
            ax.grid()
            ax.set_xticklabels(x_axis)#, rotation=90, va='bottom')  # rotate bottom→up
            ax.set_title('inv(flow)_point_wise_RHS_A1')
            ax.legend()
            # plt.show()
            #
            _, ax = plt.subplots(figsize=(11,8))
            ax.plot(uB_over_B2_0_A2.T,'.-',label=plasma_class.list_of_species)
            ax.set_xticks(range(len(x_axis)))
            ax.grid()
            ax.set_xticklabels(x_axis)#, rotation=90, va='bottom')  # rotate bottom→up
            ax.set_title('inv(flow)_point_wise_RHS_A2')
            ax.legend()
            plt.show()
            
            # Sanity Check
            # summation = np.sum((uB_over_B2_0_A1+uB_over_B2_0_A2),axis=1)
            # for ispecies,species in enumerate(plasma_class.list_of_species):
            #     print('sanity =',n[species]*q[species]*np.sqrt(self.Bsq)*summation[ispecies] )
        
        return JBS_species
    
    def get_fluxes(self,Er_Vcm,plasma_class,Smax):
        """ This function computes the neoclassical particle and heat fluxes as PENTA3 does using the Sugama-Nishimura method. This amounts to 
        solve Eq. (3.2.3.1) in PENTA documentation by Jeremy Lore.
        
        !! Er is given in V/cm !!
        
        This function returns: Gamma, QoverT
        where Gammma and QoverT are arrays with the fluxes for each species
        """
        from scipy.special import assoc_laguerre
        
        Er = Er_Vcm
        # Check Er!=0
        if(np.abs(Er) < 1E-6):
            raise ValueError('Please Choose Er != 0, cause this brings problems when taking the log of Er')
        
        ############################################################################################################################
        ####################### Setup Convolutions in the SN Flow Equation (Eq. 3.2.3.1 in PENTA's notes)  #########################
        ############################################################################################################################
        vth = {}
        n = {}
        T = {}
        q = {}
        m = {}
        #
        A1_conv_Gamma = {}
        A1_conv_QoT = {}
        A2_conv_Gamma = {}
        A2_conv_QoT = {}
        flow_conv_Gamma = {species: np.zeros((Smax+1)) for species in plasma_class.list_of_species}
        flow_conv_QoT = {species: np.zeros((Smax+1)) for species in plasma_class.list_of_species}
        #
        for species in plasma_class.list_of_species:
            vth[species] = plasma_class.get_thermal_speed(species,self.roa)
            n[species] = plasma_class.get_density(species,self.roa)
            T[species] = plasma_class.get_temperature(species,self.roa)
            q[species] = plasma_class.charge[species]
            m[species] = plasma_class.mass[species]
            
            #get cmul for self.K
            cmul_K = []
            for k in self.K:
                vparticle = vth[species] * np.sqrt(k)
                nu = plasma_class.get_collisionality(species,self.roa,vparticle)
                cmul_K.append( nu / vparticle )
                            
            cmul_K = np.array(cmul_K)
            efield_K = np.abs(Er)*100/(vth[species]*np.sqrt(self.K))
            
            D31_over_D33_corrected = self.get_interpolated_coeff('D31_over_D33_corrected',cmul_K,efield_K,log_interp_coeff=False)
            capped_fluxes_coefficient = self.get_interpolated_coeff('capped_fluxes_coefficient',cmul_K,efield_K,log_interp_coeff=False)
            
            # Convolutions
            integrand = capped_fluxes_coefficient * self.K**1.5 * np.sqrt(self.K) * np.exp(-self.K) * assoc_laguerre(self.K, 0, k=1.5)
            A1_conv_Gamma[species] = np.trapezoid(integrand,x=self.K) * n[species] * 2 / np.sqrt(np.pi)
            
            integrand = capped_fluxes_coefficient * self.K**2.5 * np.sqrt(self.K) * np.exp(-self.K) * assoc_laguerre(self.K, 0, k=1.5)
            A2_conv_Gamma[species] = np.trapezoid(integrand,x=self.K) * n[species] * 2 / np.sqrt(np.pi)
            A1_conv_QoT[species]   = np.trapezoid(integrand,x=self.K) * n[species] * 2 / np.sqrt(np.pi)
            
            integrand = capped_fluxes_coefficient * self.K**3.5 * np.sqrt(self.K) * np.exp(-self.K) * assoc_laguerre(self.K, 0, k=1.5)
            A2_conv_QoT[species] = np.trapezoid(integrand,x=self.K) * n[species] * 2 / np.sqrt(np.pi)
            
            #### Convolutions that multiply flows ####
            for jval in range(Smax+1):

                integrand = D31_over_D33_corrected * self.K**1.5 * np.sqrt(self.K) * np.exp(-self.K) * assoc_laguerre(self.K, jval, k=1.5)
                flow_conv_Gamma[species][jval] =  np.trapz(integrand,x=self.K) * n[species] * 2 / np.sqrt(np.pi)
                
                integrand = D31_over_D33_corrected * self.K**2.5 * np.sqrt(self.K) * np.exp(-self.K) * assoc_laguerre(self.K, jval, k=1.5)
                flow_conv_QoT[species][jval] =  np.trapz(integrand,x=self.K) * n[species] * 2 / np.sqrt(np.pi)
                
        
        ############################################################################################################################
        ##########################################  Compute A1 and A2 for each species #############################################
        ############################################################################################################################
        A1 = {}
        A2 = {}
        for species in plasma_class.list_of_species:
            na_prime_r = plasma_class.get_density_der(species,self.roa) / self.aminor
            Ta_prime_r = plasma_class.get_temperature_der(species,self.roa) / self.aminor
            #
            A1[species] = na_prime_r/n[species] - 1.5*Ta_prime_r/(T[species]) - q[species]*Er*100/(EC*T[species])
            A2[species] = Ta_prime_r/(T[species])
            
        ############################################################################################################################
        ###############################################  Compute lmat ##############################################################
        ############################################################################################################################
        # Calculate Coulomb logarithm as in PENTA
        if ( T['electrons'] > 50 ):
            loglambda = 25.3 - 1.15*np.log10(n['electrons']/1.E6) + 2.3*np.log10(T['electrons'])
        else:
            loglambda = 23.4 - 1.15*np.log10(n['electrons']/1.E6) + 3.45*np.log10(T['electrons'])
        #
        lmat = define_friction_coeffs(masses=np.fromiter(m.values(), dtype=float), 
                                     charges=np.fromiter(q.values(), dtype=float), 
                                     v_ths = np.fromiter(vth.values(), dtype=float), 
                                     Temps = np.fromiter(T.values(), dtype=float), 
                                     dens  = np.fromiter(n.values(), dtype=float),
                                     loglambda = loglambda, 
                                     num_species=len(plasma_class.list_of_species), 
                                     Smax=Smax)
            
        ############################################################################################################################
        ######################################################  Compute PS term ####################################################
        ############################################################################################################################
        
        PS_term_Gamma = {}
        PS_term_QoT = {}
        
        for ispec1,speciesa in enumerate(plasma_class.list_of_species):
            
            PS_term_Gamma[speciesa] = 0.0
            PS_term_QoT[speciesa] = 0.0
            
            for ispec2,speciesb in enumerate(plasma_class.list_of_species):
                nb_prime_r = plasma_class.get_density_der(speciesb,self.roa) / self.aminor
                Tb_prime_r = plasma_class.get_temperature_der(speciesb,self.roa) / self.aminor
                
                # lmat_ind1 = ( ispec1 - 1 ) * ( Smax + 1 ) + 1
                lmat_ind1 = ( ispec1 ) * ( Smax + 1 )
                # lmat_ind2 = ( ispec2 - 1 ) * ( Smax + 1 ) + 1
                lmat_ind2 = ( ispec2 ) * ( Smax + 1 )
                
                l_ab_11 = lmat[lmat_ind1    , lmat_ind2]
                l_ab_21 = lmat[lmat_ind1+1  , lmat_ind2]
                l_ab_12 = lmat[lmat_ind1    , lmat_ind2+1]
                l_ab_22 = lmat[lmat_ind1+1  , lmat_ind2+1]
                
                PS_term_Gamma[speciesa] += (nb_prime_r*T[speciesb]*EC + n[speciesb]*Tb_prime_r*EC) / (q[speciesb]*n[speciesb]) * l_ab_11 \
                                            - Tb_prime_r*EC/q[speciesb] * l_ab_12
                
                PS_term_QoT[speciesa] += (nb_prime_r*T[speciesb]*EC + n[speciesb]*Tb_prime_r*EC) / (q[speciesb]*n[speciesb]) * \
                                            ( (5/2)*l_ab_11 - l_ab_21 ) \
                                        - Tb_prime_r*EC/q[speciesb] * ( (5/2)*l_ab_12 - l_ab_22 )
                
            PS_term_Gamma[speciesa] *= self.U2/q[speciesa]
            PS_term_QoT[speciesa] *= self.U2/q[speciesa]
            
        ############################################################################################################################
        ###################################################  Compute fluxes ########################################################
        ############################################################################################################################         
            
        uB_over_B2 = self.get_BS_current(Er,plasma_class,Smax,output='flows')
        
        num_species = len(plasma_class.list_of_species)
        
        Gamma = np.zeros(num_species)
        QoverT = np.zeros(num_species
                         )
        for ispecies,species in enumerate(plasma_class.list_of_species):
            
            mono_flux_1 = - (m[species]**2 * vth[species]**3) / (2*q[species]**2)*A1_conv_Gamma[species]*A1[species]
            mono_flux_2 = - (m[species]**2 * vth[species]**3) / (2*q[species]**2)*A2_conv_Gamma[species]*A2[species]
            flux_Ua = -(2/3)*self.Bsq*m[species]*vth[species]/q[species] * np.dot(flow_conv_Gamma[species][:],uB_over_B2[ispecies*(Smax+1):ispecies*(Smax+1)+(Smax+1)])
            
            Gamma[ispecies]  = mono_flux_1 + mono_flux_2 + flux_Ua + PS_term_Gamma[species]
            
            part1 = - (m[species]**2 * vth[species]**3) / (2*q[species]**2)*A1_conv_QoT[species]*A1[species]
            part2 = - (m[species]**2 * vth[species]**3) / (2*q[species]**2)*A2_conv_QoT[species]*A2[species]
            part3 = -(2/3)*self.Bsq*m[species]*vth[species]/q[species] * np.dot(flow_conv_QoT[species][:],uB_over_B2[ispecies*(Smax+1):ispecies*(Smax+1)+(Smax+1)])
            
            QoverT[ispecies] = part1 + part2 + part3 + PS_term_QoT[species]
            
        return Gamma,QoverT
                       
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
            
    # def write_U2_to_file(self,rho,filename=None):
    #     print(f'\n ########### ASSUMING U2 IS CONSTANT: <U^2>={self.Usq} ##############')
        
    #     if filename is None:
    #         filename = 'Utilde2_profile'
            
    #     with open(filename, 'w') as file:
    #         # Write the size of rho as the first line
    #         file.write(f'{rho.size}\n')
            
    #         # Write the data to the file
    #         for i in range(rho.size):
    #             row_data = [rho[i], self.Usq]

    #             # Write the row to the file, formatted as space-separated values
    #             file.write(" ".join(map(str, row_data)) + '\n')
    
    # def compute_PENTA1_coeffs(self):
    #     # Computes PENTA input coefficients lstar, mstar and nstar
    #     # size of lstar, mstar, nstar is n_cmul x n_efield
    #     # Check DKES/PENTA documentation to see the definition of lstar, mstar, nstar
    #     # In the documentation, D_ij^* corresponds to self.Lij, which are the species-independent DKES coefficientes
    #     print('\n#############################################################################')
    #     print('###################   Computing coefficients for PENTA1/v2.0 ##################')
    #     print('###############################################################################')
        
    #     ######  WARNING: as of now this assumes an hydrogen plasma, qa=e_charge  #####
    #     print('\nWARNING: THIS ASSUMES A PLASMA WITH Z=1, qa=echarge')
        
    #     # Read Pfirsch-Schluter flow from external file
    #     # To do later...
    #     print('\nFailed to read Pfirsch-Schluter flow, <U^2>, from external file')
    #     print('Assuming <U^2>=0\n')
    #     self.Usq = 0.0
        
    #     aux = 1 - 1.5*self.cmul*self.D33_star/self.Bsq
        
    #     #compute PENTA lstar
    #     self.lstar = self.D11_star - (2./3.)*self.cmul*self.Usq + (1.5*self.cmul*self.D13_star*self.D13_star/self.Bsq)/aux
    #     self.lstar = self.lstar / (EC*EC)
        
    #     #compute PENTA mstar
    #     self.mstar = self.cmul*self.cmul*self.D33_star / aux

    #     #compute PENTA nstar
    #     self.nstar = self.cmul*self.D13_star / aux
    #     self.nstar = self.nstar / EC
    
    # def plot_PENTA1_coeffs(self):
    #     # creates 3 graphs: lstar, mstar and nstar vs cmul (for each efield)
        
    #     import matplotlib.pyplot as pyplot
        
    #     for plot_var in ['lstar', 'mstar', 'nstar']:
            
    #         yplot = getattr(self,plot_var)
            
    #         if plot_var=='lstar':
    #             yplot1 = getattr(self,plot_var+'_1')
    #             yplot2 = getattr(self,plot_var+'_2')
    #             yplot3 = getattr(self,plot_var+'_3')
                
    #         px = 1/pyplot.rcParams['figure.dpi']
    #         pyplot.rc('font', size=24)
    #         pyplot.rc('legend', fontsize=24)
    #         fig=pyplot.figure(figsize=(1024*px,768*px))
    #         ax = fig.add_subplot(111)
    #         for i in range(self.nefield):
    #             i1 = i*self.ncmul
    #             i2 = i1 + self.ncmul
    #             ax.plot(self.cmul[i1:i2],yplot[i1:i2],marker='+',label=rf'$E_s/v$={self.efield[i1]:3.1E}',linewidth=4,markersize=18)
    #             if plot_var=='lstar' and i==1:
    #                 ax.plot(self.cmul[i1:i2],yplot1[i1:i2],'--',label='D11*',linewidth=4)
    #                 ax.plot(self.cmul[i1:i2],yplot2[i1:i2],'--',label='U2*',linewidth=4)
    #                 ax.plot(self.cmul[i1:i2],yplot3[i1:i2],'--',label='(D31*)^2/D33*',linewidth=4)
                    
    #         ax.set_xlabel(r'$\nu/v [m^{-1}]$')
    #         ax.set_ylabel(f'PENTA {plot_var}')
    #         ax.set_xscale('log')
    #         if(plot_var=='lstar' or plot_var=='mstar'):
    #             ax.set_yscale('log')
    #         if(plot_var=='nstar'):
    #             ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    #         ax.set_title(f'r/a={self.roa:.2f}')
    #         ax.grid()
    #         ax.legend(fontsize=12)
     
    #     pyplot.show()
    
    # def write_PENTA1_coeffs_to_files(self,where_to):
    #     # name of the files are 'lstar_lijs_##' , 'mstar_lijs_##', 'nstar_lijs_##'
        
    #     #where_to save -- path should not have final /
        
    #     #checks that data is regular, i.e.: for each efield there are always the same cmul
    #     # if data is regular the following arrays are computed:
    #     # self.cmul_regular
    #     # self.efield_regular
        
    #     self.check_data_is_regular(self.cmul,self.efield)
        
    #     for var in ['lstar', 'mstar', 'nstar']:
            
    #         filename = where_to + '/' + var + '_lijs_' + 'surface_' + f'{self.surface}'
    #         y = getattr(self,var)

    #         combined = np.concatenate((self.cmul_regular,self.efield_regular,y))
            
    #         #create file
    #         with open(filename, 'w') as file:
    #             file.write(f'{self.ncmul} {self.nefield}\n')
    #             for value in combined:
    #                 file.write(f'{value:.10e}\n') 
                
    # def plot_PENTA_Lcoeffs_vs_roa_ambi(self,folder_path=None,plot=True):
    #     #plots L11,L12 and L13 for each species if PENTA3 was run in DKES mode
        
    #     import matplotlib.pyplot as plt
    #     import numpy as np
    #     from itertools import groupby
        
    #     if folder_path is None:
    #         filename = 'thermalCoeffs_vs_roa'
    #     else:
    #         filename = folder_path+'/thermalCoeffs_vs_roa'
            
    #     penta = np.loadtxt(filename,skiprows=2)
        
    #     roa = penta[:,0]
        
    #     num_ion_species = int(len(penta[0,2:])/3) - 1
        
    #     species_list = ['e'] + [f'i{j}' for j in range(1, num_ion_species + 1)]        
        
    #     L11 = penta[:,np.r_[2:(3+num_ion_species)]] 
    #     L12 = penta[:,np.r_[(3+num_ion_species):(4+2*num_ion_species)]]
    #     L13 = penta[:,np.r_[(4+2*num_ion_species):(5+3*num_ion_species)]]

    #     plt.rc('font', size=18)
    #     fig11, ax11 = plt.subplots(1,1,figsize=(10,7))
    #     fig12, ax12 = plt.subplots(1,1,figsize=(10,7))
    #     fig13, ax13 = plt.subplots(1,1,figsize=(10,7))

    #     #loop in species
    #     for js,s in enumerate(species_list):
    #         L11_species = L11[:,js]
    #         L12_species = L12[:,js]
    #         L13_species = L13[:,js]
            
    #         L11_electron_root = []
    #         L11_ion_root = []
            
    #         L12_electron_root = []
    #         L12_ion_root = []
            
    #         L13_electron_root = []
    #         L13_ion_root = []
            
    #         roa_electron_root = []
    #         roa_ion_root = []
            
    #         #divide diferent roots (discard unstable one)
    #         i=0
    #         for _, group in groupby(roa):
    #             num_roots = len( list(group) )
                
    #             if(num_roots == 2):
    #                 raise ValueError(f"How come you have 2 roots??")
    #             elif(num_roots == 1):
    #                 L11_ion_root.append(L11_species[i])
    #                 L12_ion_root.append(L11_species[i])
    #                 L13_ion_root.append(L11_species[i])
    #                 roa_ion_root.append(roa[i])
    #             elif(num_roots ==3):
    #                 # ion root
    #                 L11_ion_root.append(L11_species[i])
    #                 L12_ion_root.append(L12_species[i])
    #                 L13_ion_root.append(L13_species[i])
    #                 roa_ion_root.append(roa[i])
                    
    #                 # electron root
    #                 L11_electron_root.append(L11_species[i+2])
    #                 L12_electron_root.append(L12_species[i+2])
    #                 L13_electron_root.append(L13_species[i+2])
    #                 roa_electron_root.append(roa[i+2])
    #             else:
    #                 raise ValueError(f"How come you have {num_roots} roots ??")   
    #             i += num_roots
            
    #         # plot L11
    #         ax11.plot(roa_ion_root,L11_ion_root,'.-',label=species_list[js]+', i-root')
    #         ax11.plot(roa_electron_root,L11_electron_root,'.-',label=species_list[js]+', e-root')
    #         ax11.set_title(r'$L_{11}$')
    #         # plot L12
    #         ax12.plot(roa_ion_root,L12_ion_root,'.-',label=species_list[js]+', i-root')
    #         ax12.plot(roa_electron_root,L12_electron_root,'.-',label=species_list[js]+', e-root')
    #         ax12.set_title(r'$L_{12}$')
    #         # plot L13
    #         ax13.plot(roa_ion_root,L13_ion_root,'.-',label=species_list[js]+', i-root')
    #         ax13.plot(roa_electron_root,L13_electron_root,'.-',label=species_list[js]+', e-root')
    #         ax13.set_title(r'$L_{13}$')
    #     ax11.set_xlabel(r'r/a')
    #     ax12.set_xlabel(r'r/a')
    #     ax13.set_xlabel(r'r/a')
    #     #ax.legend(fontsize=12)
    #     ax11.grid()
    #     ax12.grid()
    #     ax13.grid()
    #     plt.legend()
    #     if plot:
    #         plt.show()
            
    # def plot_PENTA_Lcoeffs_vs_Er(self,roa_user,folder_path=None,plot=True):
    #     #reads the file thermalCoeffs_vs_Er creatd by PENTA code 
    #     # when ran in DKES mode
    #     #folderpath should contain path to file without final '/'
    #     # if not given, it is assumed we are already inside the folder
        
    #     import matplotlib.pyplot as plt
    #     import numpy as np
    #     from collections import defaultdict
        
    #     if folder_path is None:
    #         filename = 'thermalCoeffs_vs_Er'
    #     else:
    #         filename = folder_path+'/thermalCoeffs_vs_Er'
            
    #     penta = np.loadtxt(filename,skiprows=2)
        
    #     roa = penta[:,0]
    #     Er = penta[:,1]
        
    #     num_ion_species = int(len(penta[0,2:])/3) - 1
        
    #     species_list = ['e'] + [f'i{j}' for j in range(1, num_ion_species + 1)]        
        
    #     L11_electrons = penta[:,2]
    #     L12_electrons = penta[:,3+num_ion_species]
    #     L13_electrons = penta[:,4+2*num_ion_species]
        
    #     Er_dict = defaultdict(list)
    #     L11_dict = defaultdict(list)
    #     L12_dict = defaultdict(list)
    #     L13_dict = defaultdict(list)
        
    #     for r,er,l11,l12,l13 in zip(roa,Er,L11_electrons,L12_electrons,L13_electrons):
    #         Er_dict[r].append(er)
    #         L11_dict[r].append(l11)
    #         L12_dict[r].append(l12)
    #         L13_dict[r].append(l13)

    #     roa_non_repeated = np.unique(roa)
        
    #     #convert roa_user to array in case it is not yet
    #     roa_user = np.array(roa_user)
        
    #     for r_user in roa_user:
            
    #         #check that roa_user is inside the limits
    #         #if yes, plot the closest to roa
    #         if(r_user <= np.min(roa) or r_user>= np.max(roa)):
    #             print(f'ERROR!! The provided roa={r_user} is outside the interval of PENTA data: [{np.min(roa)},{np.max(roa)}]')
    #             exit(0)        
    #         else:
    #             roa_closest = roa_non_repeated[ np.argmin(np.abs(roa_non_repeated-r_user)) ]
            
    #         plt.rc('font', size=16)
    #         fig=plt.figure(figsize=(8,6))
    #         ax = fig.add_subplot(111)
    #         ax.plot(Er_dict[roa_closest],L11_dict[roa_closest],label=r'L11, e')
    #         ax.plot(Er_dict[roa_closest],L12_dict[roa_closest],label=r'L12, e')
    #         ax.plot(Er_dict[roa_closest],L13_dict[roa_closest],label=r'L13, e')
    #         ax.set_xlabel(r'Er [V/cm]')
    #         ax.set_title(f'r/a={roa_closest}')
    #         #ax.set_yscale('log')
    #         ax.legend(fontsize=12)
    #         ax.grid()
    #         plt.legend()
    #     if plot:
    #         plt.show()
            
    # def plot_PENTA_Lcoeff_vs_Er(self,which_coeff,which_species: int,roa_user,folder_path=None,plot=True):
    #     #reads the file thermalCoeffs_vs_Er creatd by PENTA code 
    #     # when ran in DKES mode
    #     #folderpath should contain path to file without final '/'
    #     # if not given, it is assumed we are already inside the folder
    #     # Lcoeff should be 'L11','L12' or 'L13'
    #     # which_species should be an integer: 0 for electrons, 1 for ion1, 2 for ion2, ...
        
    #     import matplotlib.pyplot as plt
    #     import numpy as np
    #     from collections import defaultdict
        
    #     if folder_path is None:
    #         filename = 'thermalCoeffs_vs_Er'
    #     else:
    #         filename = folder_path+'/thermalCoeffs_vs_Er'
            
    #     penta = np.loadtxt(filename,skiprows=2)
        
    #     roa = penta[:,0]
    #     Er = penta[:,1]
        
    #     num_ion_species = int(len(penta[0,2:])/3) - 1
        
    #     species_list = ['e'] + [f'i{j}' for j in range(1, num_ion_species + 1)] 
        
    #     try:
    #         print(f'Plotting {which_coeff} coefficient for species {species_list[which_species]}')  
    #     except:
    #         print('ERROR: THE GIVEN SPECIES DOES NOT EXIST !!')
    #         exit(0)
        
    #     if(which_coeff=='L11'):
    #         L = penta[:,2+which_species]   
    #     elif(which_coeff=='L12'):
    #         L = penta[:,3+num_ion_species+which_species]
    #     elif(which_coeff=='L13'):
    #         L = penta[:,4+2*num_ion_species+which_species]
    #     else:
    #         print(f'ERROR!! The coefficient type must be L11, L12 or L13. Instead you provided {which_coeff}')
    #         exit(0)
                    
    #     Er_dict = defaultdict(list)
    #     L_dict = defaultdict(list)
        
    #     for r,er,l in zip(roa,Er,L):
    #         Er_dict[r].append(er)
    #         L_dict[r].append(l)

    #     roa_non_repeated = np.unique(roa)
        
    #     #convert roa_user to array in case it is not yet
    #     roa_user = np.array(roa_user)
        
    #     plt.rc('font', size=16)
    #     fig=plt.figure(figsize=(8,6))
    #     ax = fig.add_subplot(111)
    #     for r_user in roa_user:
            
    #         #check that roa_user is inside the limits
    #         #if yes, plot the closest to roa
    #         if(r_user <= np.min(roa) or r_user>= np.max(roa)):
    #             print(f'ERROR!! The provided roa={r_user} is outside the interval of PENTA data: [{np.min(roa)},{np.max(roa)}]')
    #             exit(0)        
    #         else:
    #             roa_closest = roa_non_repeated[ np.argmin(np.abs(roa_non_repeated-r_user)) ]
            
    #         ax.plot(Er_dict[roa_closest],L_dict[roa_closest],label=f'r/a={roa_closest}')
    #         ax.set_xlabel(r'Er [V/cm]')
    #         ax.set_title(f'{which_coeff} of {species_list[which_species]}')
    #         #ax.set_yscale('log')
    #         ax.legend(fontsize=12)
    #         ax.grid()
    #         plt.legend()
    #     if plot:
    #         plt.show()
            
    # def plot_PENTA_Lcoeff_vs_roa(self,which_coeff,which_species: int,Er_user,folder_path=None,plot=True):
    #     # plots PENTA Lcoeff as a afunction of r/a for a given Er=Er_user. Er_user can be an array
    #     # Er_user should come in V/cm
    #     # folder_path is the folder that contains tha thermlCoeffs_vs_Er file
    #     # note that this is not plotting the 'ambipolar coefficient', i.e., the coefficient at the Er
    #     # set b ambipolarity; here, instead, we choose at which Er we want to plot
        
    #     import matplotlib.pyplot as plt
    #     import numpy as np
    #     from collections import defaultdict
        
    #     if folder_path is None:
    #         filename = 'thermalCoeffs_vs_Er'
    #     else:
    #         filename = folder_path+'/thermalCoeffs_vs_Er'
            
    #     penta = np.loadtxt(filename,skiprows=2)
        
    #     roa = penta[:,0]
    #     Er = penta[:,1]
        
    #     num_ion_species = int(len(penta[0,2:])/3) - 1
        
    #     species_list = ['e'] + [f'i{j}' for j in range(1, num_ion_species + 1)] 
        
    #     try:
    #         print(f'Plotting {which_coeff} coefficient for species {species_list[which_species]}')  
    #     except:
    #         print('ERROR: THE GIVEN SPECIES DOES NOT EXIST !!')
    #         exit(0)
        
    #     if(which_coeff=='L11'):
    #         L = penta[:,2+which_species]   
    #     elif(which_coeff=='L12'):
    #         L = penta[:,3+num_ion_species+which_species]
    #     elif(which_coeff=='L13'):
    #         L = penta[:,4+2*num_ion_species+which_species]
    #     else:
    #         print(f'ERROR!! The coefficient type must be L11, L12 or L13. Instead you provided {which_coeff}')
    #         exit(0)
                    
    #     L_dict = defaultdict(list)
        
    #     for r,er,l in zip(roa,Er,L):
    #         L_dict[r,er].append(l)

    #     roa_unique = np.unique(roa)
    #     Er_unique = np.unique(Er)
        
    #     #convert roa_user to array in case it is not yet
    #     Er_user = np.array(Er_user)
        
    #     plt.rc('font', size=16)
    #     fig=plt.figure(figsize=(8,6))
    #     ax = fig.add_subplot(111)   
        
    #     for E_user in Er_user:
            
    #         #check that roa_user is inside the limits
    #         #if yes, plot the closest to roa
    #         if(E_user <= np.min(Er) or E_user>= np.max(Er)):
    #             print(f'ERROR!! The provided Er={Er_user} is outside the interval of PENTA data: [{np.min(Er)},{np.max(Er)}]')
    #             exit(0)        
    #         else:
    #             Er_closest = Er_unique[ np.argmin(np.abs(Er_unique-E_user)) ]
                
    #         L_plot = [l for r in roa_unique for l in L_dict[(r,Er_closest)]]
            
    #         ax.plot(roa_unique,L_plot,'.-',label=f'Er={Er_closest} V/cm')
    #         ax.set_xlabel(r'r/a')
    #         ax.set_title(f'{which_coeff} of {species_list[which_species]}')
    #         #ax.set_yscale('log')
    #         ax.legend(fontsize=12)
    #         ax.grid()
    #         plt.legend()
    #     if plot:
    #         plt.show()
                
    # def plot_PENTA1_Er_vs_roa(self,folder_path=None,plot=True):
        
    #     import matplotlib.pyplot as plt
    #     import numpy as np
    #     from collections import defaultdict
        
    #     if folder_path is None:
    #         filename = 'flows_vs_roa'
    #     else:
    #         filename = folder_path+'/flows_vs_roa'
            
    #     penta = np.loadtxt(filename,skiprows=2)
        
    #     roa = penta[:,0]
    #     Er = penta[:,1]
        
    #     #separate different roots
    #     # Initialize dictionaries to store the separated values
    #     roa_ion_root = []
    #     roa_unstable_root = []
    #     roa_electron_root = []
    #     Er_ion_root = []
    #     Er_unstable_root = []
    #     Er_electron_root = []

    #     # Create a dictionary to group all Er values by corresponding roa
    #     roa_to_Er = defaultdict(list)

    #     # Group Er values by roa
    #     for r, e in zip(roa, Er):
    #         roa_to_Er[r].append(e)

    #     # Iterate over the grouped roa values
    #     for r, er_values in roa_to_Er.items():
    #         # Sort the corresponding Er values for each roa
    #         sorted_er = sorted(er_values)
            
    #         if len(sorted_er) == 1:
    #             # If only one value of roa, assign it to ion root
    #             roa_ion_root.append(r)
    #             Er_ion_root.append(sorted_er[0])
            
    #         # elif len(sorted_er) == 2:
    #         #     print(f'ERROR: HOW COME YOU HAVE TWO ROOTS AT r/a={r} ???')
    #         #     exit()
            
    #         elif len(sorted_er) == 3:
    #             # If three Er values exist, assign smallest to ion root, middle to unstable root, and largest to electron root
    #             roa_ion_root.append(r)
    #             Er_ion_root.append(sorted_er[0])
    #             roa_unstable_root.append(r)
    #             Er_unstable_root.append(sorted_er[1])
    #             roa_electron_root.append(r)
    #             Er_electron_root.append(sorted_er[2])

    #     # Convert results to numpy arrays (optional, if you need them as arrays)
    #     roa_ion_root = np.array(roa_ion_root)
    #     Er_ion_root = np.array(Er_ion_root)
    #     roa_unstable_root = np.array(roa_unstable_root)
    #     Er_unstable_root = np.array(Er_unstable_root)
    #     roa_electron_root = np.array(roa_electron_root)
    #     Er_electron_root = np.array(Er_electron_root)
        
    #     plt.rc('font', size=18)
    #     fig=plt.figure(figsize=(11,8))
    #     ax = fig.add_subplot(111)
    #     ax.plot(roa_ion_root,Er_ion_root,'.-',label='ion root')
    #     ax.plot(roa_unstable_root,Er_unstable_root,'.-',label='unstable root')
    #     ax.plot(roa_electron_root,Er_electron_root,'.-',label='electron root')
    #     ax.set_ylabel(r'Er [V/cm]')
    #     ax.set_xlabel(r'r/a')
    #     ax.legend(fontsize=12)
    #     ax.grid()
    #     plt.legend()
    #     if plot:
    #         plt.show()
    
    # def set_PENTA1_integrands_energy_conv(self,intj,plasma_class,make_plots=True):
    #     # This function computes the integrand of the energy convolution as in PENTA for each efield
    #     # and plots it as function of K
    #     # Integrand = sqrt(K) * exp(-K) * (K-5/2)^{intj-1} * [lstar,m,star,nstar] * K^{3/2}
    #     # This function requires computing collisionality nu_D
    #     # We also spline interpolate lstar,mstar and star as function of cmul for each efield
        
    #     import matplotlib.pyplot as plt
        
    #     self.plasma_class = plasma_class
        
    #     # Computes dicitionary of arrays self.cmul_species. 
    #     # Contains cmul for each species for array K
    #     self.set_cmul_species(K,make_plot=make_plots)
        
    #     # Set integrands = l/m/n-star * fix func
    #     # this creates dictionary of arrays: self.lstar_integrand, self.nstar_integrand, self.mstar_integrand
    #     # for instance, self.lstar_integrand['electrons'][3] gives the arrays of integrand (as function of K) 
    #     # for l* for electrons for the 4th (3+1) electric field 
    #     self.set_integrands(intj,K,make_plot=make_plots)  
            
            
    # def get_fix_func(self,K,intj,make_plot=False):
        
    #     import matplotlib.pyplot as plt
        
    #     fix_func = np.sqrt(K) * np.exp(-K) * (K-2.5)**(intj-1) * K**1.5
        
    #     if make_plot is True:
    #         fig, ax = plt.subplots(figsize=(8,6))
    #         ax.plot(K,fix_func,'o-')
    #         #ax.set_yscale('log')
    #         ax.set_ylabel(r'$\sqrt{K}e^{-K}\left(K-5/2\right)^{j-1}\,K^{3/2}$')
    #         ax.set_xlabel(f'K')   
    #         ax.set_title(f'j={intj}')
    #         ax.grid()    
    #         plt.show()
        
    #     return fix_func
    
    # def set_integrands(self,intj,K,make_plot=False):
        
    #     from scipy.interpolate import interp1d
    #     import matplotlib.pyplot as plt
    #     from collections import defaultdict
        
    #     print('\n ############################################################')
    #     print('############# COMPUTING INTEGRANDS AS IN PENTA #################')
    #     print('############################################################')
        
    #     # fix_func = f_j(K)*K^(3/2)
    #     fix_func = self.get_fix_func(K,intj,make_plot=False)
        
    #     # create dicionaries of lists
    #     self.lstar_integrand = defaultdict(list)
    #     self.mstar_integrand = defaultdict(list)
    #     self.nstar_integrand = defaultdict(list)
        
    #     for i in range(self.nefield):
    #         i1 = i*self.ncmul
    #         i2 = i1 + self.ncmul
            
    #         efield = self.efield[i1]
            
    #         x = self.cmul[i1:i2]
    #         yl = self.lstar[i1:i2]
    #         yn = self.nstar[i1:i2]
    #         ylogm = np.log( self.mstar[i1:i2] )
            
    #         # quadratic spline as in PENTA. Assuming log_interp = true
    #         xlog = np.log(x)
    #         lstar_interp = interp1d(xlog,yl,kind='quadratic',bounds_error=False,fill_value=0.0)
    #         nstar_interp = interp1d(xlog,yn,kind='quadratic',bounds_error=False,fill_value=0.0)
    #         logmstar_interp = interp1d(xlog,ylogm,kind='quadratic',bounds_error=False,fill_value=0.0)
    #         #this function is used to multiply exp(logmstar_interp1d), otherwise exp(0)=1 is taken outside the interpolating region
    #         filter_logmstar = interp1d(xlog,np.ones_like(xlog),bounds_error=False,fill_value=0.0)
            
    #         xspline = np.logspace(np.log10(x[0]),np.log10(x[-1]),100)
            
    #         fig, ax = plt.subplots(1,3,figsize=(17,6))
    #         ax[0].plot(x,yl,'ob')
    #         ax[0].plot(xspline, lstar_interp(np.log(xspline)),'red',label='spline')
    #         ax[0].set_yscale('log')
    #         ax[0].set_xscale('log')
    #         ax[0].set_ylabel(r'lstar')
    #         ax[0].set_xlabel(f'cmul')   
    #         ax[0].set_title(f'Er/v={efield}')
    #         ax[0].grid()
    #         ax[0].legend()
            
    #         ax[1].plot(x,yn,'ob')
    #         ax[1].plot(xspline, nstar_interp(np.log(xspline)),'red',label='spline')
    #         ax[1].set_xscale('log')
    #         ax[1].set_ylabel(r'nstar')
    #         ax[1].set_xlabel(f'cmul')   
    #         ax[1].set_title(f'Er/v={efield}')
    #         ax[1].grid()   
    #         #ax[1].legend()
            
    #         ax[2].plot(x,ylogm,'ob')
    #         ax[2].plot(xspline, logmstar_interp(np.log(xspline)),'red',label='spline')
    #         ax[2].set_xscale('log')
    #         ax[2].set_ylabel('ln(mstar)')
    #         ax[2].set_xlabel(f'cmul')   
    #         ax[2].set_title(f'Er/v={efield}')
    #         ax[2].grid()   
    #         #ax[1].legend()
            
    #         plt.tight_layout(pad=2)
            
    #         # full integrand of l*
    #         fig,ax = plt.subplots(figsize=(8,6))
    #         for species in self.plasma_class.list_of_species:
    #             integrand = lstar_interp(np.log(self.cmul_species[species]))*fix_func
    #             integral = self.get_integral(integrand,K)
                
    #             #save integrand
    #             self.lstar_integrand[species].append( integrand )
                
    #             #plot
    #             ax.plot(K,integrand,'o-',label=f'{species}, {integral:.3e}')       
    #             ax.set_xlabel('K')
    #             ax.set_ylabel(fr'$f_{intj}(K)~l^*(K)~K^{{3/2}}$')
    #             ax.set_title(f'Er/v={efield}')
    #             ax.grid()
    #         plt.legend()
                      
    #         # full integrand of n*
    #         fig,ax = plt.subplots(figsize=(10,6))
    #         for species in self.plasma_class.list_of_species:
    #             integrand = nstar_interp(np.log(self.cmul_species[species]))*fix_func
    #             integral = self.get_integral(integrand,K)
                
    #             #save integrand
    #             self.nstar_integrand[species].append( integrand )                
                
    #             #plot
    #             ax.plot(K,integrand,'o-',label=f'{species}, {integral:.3e}')       
    #             ax.set_xlabel('K')
    #             ax.set_ylabel(f'$f_{intj}(K)~n^*(K)~K^{{3/2}}$')
    #             ax.set_title(f'Er/v={efield}')
    #             ax.grid()
    #         plt.legend()
            
    #         # full integrand of m*
    #         fig,ax = plt.subplots(figsize=(10,6))
    #         for species in self.plasma_class.list_of_species:
    #             integrand = np.exp(logmstar_interp(np.log(self.cmul_species[species])))*filter_logmstar(np.log(self.cmul_species[species]))*fix_func
    #             integral = self.get_integral(integrand,K)
                
    #             #save integrand
    #             self.mstar_integrand[species].append( integrand )
                
    #             #plot
    #             ax.plot(K,integrand,'o-',label=f'{species}, {integral:.3e}')       
    #             ax.set_xlabel('K')
    #             ax.set_ylabel(f'$f_{intj}(K)~m^*(K)~K^{{3/2}}$')
    #             ax.set_title(f'Er/v={efield}')
    #             ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    #             ax.grid()
    #         plt.legend()
    #         if make_plot:
    #             plt.show()
    #         else:
    #             plt.close('all')
    
    # def get_integral(self,y,x,xmin=None,xmax=None,make_plot=False,plot_title=None):
    #     # this function computes the trapezoid integral of y=y(x)
    #     # if xmin or xmax are provided, the limits of the integral are changed
    #     # an error is raised if xmin or xmax fall outside the domain defined by x
    #     # if the array x does not contain xmin or xmax exactly, then the closest value is considered
        
    #     import matplotlib.pyplot as plt
    #     from scipy.integrate import trapezoid
        
    #     # Check if xmin and xmax are within the domain of x
    #     if xmin is not None and (xmin < x.min() or xmin > x.max()):
    #         raise ValueError(f"xmin ({xmin}) is outside the domain of x.")
    #     if xmax is not None and (xmax < x.min() or xmax > x.max()):
    #         raise ValueError(f"xmax ({xmax}) is outside the domain of x.")

    #     # If xmin or xmax are provided, adjust the limits
    #     if xmin is not None:
    #         xmin_index = np.argmin(np.abs(x - xmin))  # Find closest value to xmin in x
    #     else:
    #         xmin_index = 0  # If xmin is None, start from the beginning

    #     if xmax is not None:
    #         xmax_index = np.argmin(np.abs(x - xmax))  # Find closest value to xmax in x
    #     else:
    #         xmax_index = len(x) - 1  # If xmax is None, go to the end

    #     # Perform the trapezoidal integration
    #     x_selected = x[xmin_index:xmax_index+1]
    #     y_selected = y[xmin_index:xmax_index+1]
    #     integral = trapezoid(y_selected, x_selected)
        
    #     integral_exact = trapezoid(y,x)
    #     rel_error = np.abs(integral_exact-integral) / integral_exact

    #     # Optionally plot the function
    #     if make_plot:
    #         plt.rc('font', size=16)
    #         fig=plt.figure(figsize=(10,8))
    #         plt.plot(x, y, '.-')
    #         plt.fill_between(x_selected, y_selected, alpha=0.3, label=f'rel. error={rel_error*100:.1f}%')
    #         plt.xlabel('x')
    #         plt.legend()
    #         if plot_title is not None:
    #             plt.title(plot_title)
    #         plt.show()

    #     return integral
    
    # def compute_energy_convolution(self,which_convol, cmin, cmax):
    #     #which_convol should be 'lstar', 'mstar or nstar'
        
    #     convol_type = ['lstar','mstar','nstar']
        
    #     if which_convol not in convol_type:
    #         print(f'ERROR: which_convol should take one of the following: {convol_type}')
    #         exit(1)
            
    #     #check if integrand exist
    #     if not hasattr(self,which_convol+'_integrand'):
    #         print('ERROR: integrand does not exist! Need to set it up first!!')
    #         exit(1)
    #     else:
    #         integrand = getattr(self,which_convol+'_integrand')
            
    #     #check if cmin and cmax are inside the cmul domain
    #     if cmin>cmax or cmin<np.min(self.cmul) or cmax>np.max(self.cmul):
    #         print('ERROR: limits of integral not correct. Cannot proceed')
    #         exit(1)
            
    #     #loop in field
    #     for i in range(self.nefield):
    #         i1 = i*self.ncmul
    #         #i2 = i1 + self.ncmul
            
    #         efield = self.efield[i1]
          
    #         for species in self.plasma_class.list_of_species:
                
    #             # compute Kmin and Kmax according to cmin and cmax
    #             # we take the values in self.cmul closest to cmin and cmax
    #             cmin_index = np.argmin(np.abs(self.cmul_species[species] - cmin))
    #             cmax_index = np.argmin(np.abs(self.cmul_species[species] - cmax))
            
    #             Kmin = self.K[np.min([cmin_index,cmax_index])]
    #             Kmax = self.K[np.max([cmin_index,cmax_index])]
                
    #             self.get_integral(integrand[species][i],self.K,xmin=Kmin,xmax=Kmax,make_plot=True,plot_title=which_convol+f', {species}, Er/v={efield}')

import math
def define_friction_coeffs(masses, charges, v_ths, Temps, dens,
                           loglambda, num_species, Smax):
    """
    Compute the classical friction coefficients lmat.

    Parameters:
      masses, charges, v_ths, Temps, dens : arrays (length num_species)
      loglambda : float
      num_species : int
      Smax : int

    Returns:
      lmat : numpy array shape ((Smax+1)*num_species, (Smax+1)*num_species)
    """

    # ---------------------------
    # Low-level helpers (local)
    # ---------------------------
    def ifactorial(n: int) -> int:
        if n < 0:
            raise ValueError("ifactorial: negative input")
        return math.factorial(n)

    def Gamma_aux(X: float) -> float:
        """
        Port of the Fortran Gamma_aux:
        - If X is positive integer -> return (X-1)! exactly.
        - If X is non-positive integer -> raise ValueError (singularity).
        - Otherwise evaluate using polynomial coefficients G (Fortran code).
        """
        PI = math.pi
        tol = 1e-12

        # Check integer
        xi = int(round(X))
        if abs(X - xi) < tol:
            # treat as integer
            if xi > 0:
                return float(ifactorial(xi - 1))
            else:
                raise ValueError("Gamma_aux: negative integer input (singularity)")

        # Non-integer case
        if abs(X) > 1.0:
            Z = abs(X)
            M = int(Z)
            R = 1.0
            for K in range(1, M + 1):
                R = R * (Z - K)
            Z = Z - M
        else:
            Z = X
            R = 1.0

        # Coefficients G (26 entries) from Fortran
        G = [1.0,
             0.5772156649015329,
             -0.6558780715202538,
             -0.0420026350340952,
             0.1665386113822915,
             -0.0421977345555443,
             -0.0096219715278770,
             0.0072189432466630,
             -0.0011651675918591,
             -0.0002152416741149,
             0.0001280502823882,
             -0.0000201348547807,
             -0.00000012504934821,
             0.00000011330272320,
             -0.000000002056338417,
             0.000000000061160950,
             0.0000000000050020075,
             -0.0000000000011812746,
             0.0000000000001043427,
             0.000000000000000077823,
             -0.000000000000000036968,
             0.00000000000000000051,
             -0.000000000000000000000206,
             -0.00000000000000000000000054,
             0.00000000000000000000000014,
             0.000000000000000000000000001]
        # Note: last two tiny entries adjusted to match Fortran's Data (approx)

        GR = G[-1]
        # Evaluate polynomial via Horner from G[24] down to G[0]
        for K in range(len(G) - 2, -1, -1):
            GR = GR * Z + G[K]

        # GR * Z might be zero if Z==0 (shouldn't happen since X not integer)
        GA = 1.0 / (GR * Z)
        if abs(X) > 1.0:
            GA = GA * R
            if X < 0.0:
                GA = -PI / (X * GA * math.sin(PI * X))
        return float(GA)

    # ---------------------------
    # eab / Eba
    # ---------------------------
    sqrtpi = math.sqrt(math.pi)

    def calc_eab(k: int, Q: float) -> float:
        Xab = Q ** (-(k + 0.5))
        Gamma_k = Gamma_aux(k + 0.5)
        eab = Gamma_k * Xab / sqrtpi
        return float(eab)

    def calc_Eba(k: int, chi: float, Q: float) -> float:
        coeff = ifactorial(k) / (2.0 * chi)
        Etmp = 0.0
        for jtest in range(0, k + 1):
            Etmp += (chi ** (1 + 2 * jtest)) * calc_eab(jtest, Q) / ifactorial(jtest)
        Eba = coeff * Etmp
        return float(Eba)

    # ---------------------------
    # Astar / Bstar
    # ---------------------------
    def calc_Astar(qtest: int, mtest: int, theta: float, mu: float, chi: float, Q: float) -> float:
        iOne = 1
        iTwo = 2
        One = 1.0
        Two = 2.0

        a_big_E0 = -Two * (One - theta) / mu
        a_big_E1 = -One + (One + Two * mtest) * (One - Two * theta) / mu
        a_big_E2 = Two * theta * (mtest ** iTwo) / mu
        a_lit_e0 = Two * (One - theta) * (One / theta + One / mu)
        a_lit_e1 = (One + Two * mtest) * (One - One / mu + Two * theta / mu)
        a_lit_e2 = -Two * theta * (mtest ** iTwo) / mu

        alpha_big_E0 = calc_Eba(iOne + qtest + mtest, chi, Q)
        alpha_big_E1 = calc_Eba(qtest + mtest, chi, Q)
        if (qtest == 0) and (mtest == 0):
            alpha_big_E2 = 0.0
        else:
            alpha_big_E2 = calc_Eba(qtest + mtest - iOne, chi, Q)

        alpha_lit_e0 = (chi ** (iOne + iTwo * (iTwo + qtest + mtest))) * calc_eab(iTwo + qtest + mtest, Q) / chi
        alpha_lit_e1 = (chi ** (iOne + iTwo * (iOne + qtest + mtest))) * calc_eab(iOne + qtest + mtest, Q) / chi
        alpha_lit_e2 = (chi ** (iOne + iTwo * (qtest + mtest))) * calc_eab(qtest + mtest, Q) / chi

        Astar = (a_big_E0 * alpha_big_E0 + a_big_E1 * alpha_big_E1 + a_big_E2 * alpha_big_E2
                 + a_lit_e0 * alpha_lit_e0 + a_lit_e1 * alpha_lit_e1 + a_lit_e2 * alpha_lit_e2)
        return float(Astar)

    def calc_Bstar(qtest: int, mtest: int, theta: float, mu: float, chi: float, Q: float) -> float:
        iOne = 1
        iTwo = 2
        iThree = 3
        iFive = 5
        One = 1.0
        Two = 2.0
        Three = 3.0
        Four = 4.0
        Five = 5.0
        Eight = 8.0

        b_e0 = Two * mu / (theta ** iTwo)
        b_e1 = - Four / Five
        b_p0 = Four + Eight / Five * mtest + Four / Three * mu / theta - Eight / Three / theta
        b_m0 = - Eight / Three * mu / theta + Four / Three / theta
        b_m1 = Eight / Five

        beta_e0 = (chi ** (iTwo * qtest + iFive)) * calc_eab(iTwo + qtest + mtest, Q)
        beta_e1 = (chi ** (iTwo * qtest + iFive)) * calc_eab(iThree + qtest + mtest, Q)

        back_p0 = 0.0
        for j in range(0, qtest + 1):
            back_p0 += ifactorial(qtest) / ifactorial(j) * (chi ** (iTwo * j)) * calc_eab(iTwo + mtest + j, Q)

        back_m0 = 0.0
        back_m1 = 0.0
        for j in range(0, mtest + 1):
            back_m0 += 1.0 / ifactorial(j) * calc_eab(iTwo + qtest + j, Q)
            back_m1 += 1.0 / ifactorial(j) * calc_eab(iThree + qtest + j, Q)

        beta_p0 = 0.5 * (chi ** iThree) * back_p0
        beta_m0 = 0.5 * (chi ** (iTwo * qtest + iFive)) * ifactorial(mtest) * back_m0
        beta_m1 = 0.5 * (chi ** (iTwo * qtest + iFive)) * ifactorial(mtest) * back_m1

        Bstar = (b_e0 * beta_e0 + b_e1 * beta_e1 + b_p0 * beta_p0
                 + b_m0 * beta_m0 + b_m1 * beta_m1)
        return float(Bstar)

    # ---------------------------
    # Mak / Nab
    # ---------------------------
    def calc_Mak_Ji(irow: int, icol: int, theta_ak: float, mu_ak: float, chi_ak: float, Q_ak: float) -> float:
        Gamma_p = Gamma_aux(irow + 2.5)
        Gamma_k = Gamma_aux(icol + 2.5)
        Apk = 0.0
        for qtest in range(0, irow + 1):
            Gamma_q = Gamma_aux(qtest + 2.5)
            cpq = ((-1.0) ** qtest) * Gamma_p / (Gamma_q * ifactorial(irow - qtest) * ifactorial(qtest))
            for mtest in range(0, icol + 1):
                Gamma_m = Gamma_aux(mtest + 2.5)
                ckm = ((-1.0) ** mtest) * Gamma_k / (Gamma_m * ifactorial(icol - mtest) * ifactorial(mtest))
                Astar = calc_Astar(qtest, mtest, theta_ak, mu_ak, chi_ak, Q_ak)
                Apk += cpq * ckm * Astar
        return 2.0 * Apk

    def calc_Nab_Ji(ptest: int, ktest: int, theta: float, mu: float, chi: float, Q: float) -> float:
        Gamma_p = Gamma_aux(ptest + 2.5)
        Gamma_k = Gamma_aux(ktest + 2.5)
        Bpk = 0.0
        for qtest in range(0, ptest + 1):
            Gamma_q = Gamma_aux(qtest + 2.5)
            cpq = ((-1.0) ** qtest) * Gamma_p / (Gamma_q * ifactorial(ptest - qtest) * ifactorial(qtest))
            for mtest in range(0, ktest + 1):
                Gamma_m = Gamma_aux(mtest + 2.5)
                ckm = ((-1.0) ** mtest) * Gamma_k / (Gamma_m * ifactorial(ktest - mtest) * ifactorial(mtest))
                Bstar = calc_Bstar(qtest, mtest, theta, mu, chi, Q)
                Bpk += cpq * ckm * Bstar
        return (2.0 / chi) * Bpk

    # ---------------------------
    # Main lmat build
    # ---------------------------
    eps0 = 8.854187817e-12
    pi_val = math.pi
    tau_coeff = 3.0 * math.sqrt(pi_val) * pi_val * eps0 ** 2

    dim = (Smax + 1) * num_species
    lmat = np.zeros((dim, dim), dtype=float)

    for ispec1 in range(num_species):
        ma = masses[ispec1]
        Ta = Temps[ispec1]
        vth_a = v_ths[ispec1]
        charge_a = charges[ispec1]
        na = dens[ispec1]

        for ispec2 in range(num_species):
            lab = np.zeros((Smax + 1, Smax + 1), dtype=float)

            mb = masses[ispec2]
            Tb = Temps[ispec2]
            vth_b = v_ths[ispec2]
            charge_b = charges[ispec2]
            nb = dens[ispec2]

            chi_ab = vth_b / vth_a
            theta_ab = Tb / Ta
            mu_ab = mb / ma
            Q_ab = 1.0 + chi_ab ** 2

            for irow in range(0, Smax + 1):
                for icol in range(0, Smax + 1):
                    val = 0.0

                    # Diagonal M-sum only if same species
                    if ispec1 == ispec2:
                        M_ovr_tau_sum = 0.0
                        for spec_k in range(num_species):
                            mk = masses[spec_k]
                            Tk = Temps[spec_k]
                            vth_k = v_ths[spec_k]
                            nk = dens[spec_k]
                            charge_k = charges[spec_k]

                            tau_ak = (tau_coeff * ma ** 2 * vth_a ** 3 /
                                      (nk * charge_a ** 2 * charge_k ** 2 * loglambda))

                            chi_ak = vth_k / vth_a
                            theta_ak = Tk / Ta
                            mu_ak = mk / ma
                            Q_ak = 1.0 + chi_ak ** 2

                            Mak = calc_Mak_Ji(irow, icol, theta_ak, mu_ak, chi_ak, Q_ak)
                            M_ovr_tau_sum += Mak / tau_ak

                        val += na * ma * M_ovr_tau_sum

                    # N term (always present)
                    tau_ab = (tau_coeff * ma ** 2 * vth_a ** 3 /
                              (nb * charge_a ** 2 * charge_b ** 2 * loglambda))
                    Nab = calc_Nab_Ji(irow, icol, theta_ab, mu_ab, chi_ab, Q_ab)
                    val += na * ma * (Nab / tau_ab)

                    lab[irow, icol] = val

            ind1 = ispec1 * (Smax + 1)
            ind2 = ispec2 * (Smax + 1)
            lmat[ind1:ind1 + (Smax + 1), ind2:ind2 + (Smax + 1)] = lab

    return lmat

  
# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)