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
    
    def __init__(self, surface, eps_rel=0.03):
        self.eps_rel = eps_rel
        
        #check surface is an integer
        if not isinstance(surface,int):
            print('ERROR: surface must be an integer')
            exit(1)
        else:
            self.surface = surface
        
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
        self.L11  = 0.5*(self.L11m+self.L11p)
        self.L31  = 0.5*(self.L31m+self.L31p)
        self.L33  = 0.5*(self.L33m+self.L33p)
        
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
                'L11': '$D_{11}^*K^{3/2}~~[m^{-1}~T^{-2}]$',
                'L31': '$D_{31}^*K$',
                'L33': '$D_{33}^*K^{1/2}~~[m~T^2]$'
                }
    
        for plot_var in ['L11', 'L31', 'L33']:
            
            print(f'Plotting {plot_var}')
            
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
            if(plot_var=='L11' or plot_var=='L33'):
                ax.set_yscale('log')
            if(plot_var=='L31'):
                ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
            ax.set_title('r/a=??')
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
            
            
    def compute_PENTA_coeffs(self,wout_filename):
        # Computes PENTA input coefficients lstar, mstar and nstar
        # size of lstar, mstar, nstar is n_cmul x n_efield
        # Check DKES/PENTA documentation to see the definition of lstar, mstar, nstar
        # In the documentation, D_ij^* corresponds to self Lvar, which are the species-independent DKES coefficientes
        print('\n#############################################################################')
        print('###################   Computing PENTA coefficients  #########################')
        print('#############################################################################')
        
        ######  WARNING: as of now this assumes an hydrogen plasma, qa=e_charge  #####
        print('\nWARNING: THIS ASSUMES A PLASMA WITH Z=1, qa=echarge')
        
        # Read Pfirsch-Schluter flow from external file
        # To do later...
        print('\nFailed to read Pfirsch-Schluter flow, <U^2>, from external file')
        print('Assuming <U^2>=0\n')
        self.Usq = 0
        
        ##########################################################
        ####### Read Bsq from VMEC wout file #####################
        # maybe there is a better way through libstell library?
        self.get_Bsq(wout_filename)
        # get r/a for the surface; this will be needed in plot_PENTA_integrands
        self.get_roa(wout_filename)
        
        #compute PENTA lstar
        aux = 1 - 1.5*self.cmul*self.L33/self.Bsq
        self.lstar = self.L11 - (2./3.)*self.cmul*self.Usq + (1.5*self.cmul*self.L31*self.L31/self.Bsq)/aux
        self.lstar = self.lstar / (EC*EC)
        
        #compute PENTA mstar
        self.mstar = self.cmul*self.cmul*self.L33 / aux
        
        #compute PENTA nstar
        self.nstar = self.cmul*self.L31 / aux
        self.nstar = self.nstar / EC
        
    def get_Bsq(self,wout_filename):
        import netCDF4 as nc
        try:
            dataset = nc.Dataset(wout_filename, 'r')
            Bsq_half = dataset.variables['bdotb'][:]
            dataset.close()
        except:
            print('\nERROR: Could not read wout file')
            sys.exit(0)
        ##########################################################
        ##########################################################
        
        #interpolate Bsq from half to full grid
        Bsq_full = self.h2f(Bsq_half)
        
        #get Bsq at the required surface
        # note the -1 because python starts in 0
        self.Bsq = Bsq_full[self.surface-1]
        print(f'Bsq={self.Bsq}')
        
    def get_roa(self,wout_filename):
        import netCDF4 as nc
        try:
            dataset = nc.Dataset(wout_filename, 'r')
            phi = dataset.variables['phi'][:]
            dataset.close()
        except:
            print('\nERROR: Could not read wout file')
            sys.exit(0)

        self.roa = np.sqrt(phi[self.surface-1]/phi[-1])
        
        print(f'r/a = {self.roa}')
     
    def plot_PENTA_coeffs(self):
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
            ax.set_title('GIGA (V120)')
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
    
    def write_PENTA_coeffs_to_files(self,where_to):
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
                    file.write(f'{value}\n')
                

                
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
                
                
                
            
            
        
        
                
                
        
    
    def plot_PENTA_integrands_energy_conv(self,intj,plasma_class):
        # This function computes the integrand of the energy convolution as in PENTA for each efield
        # and plots it as function of K
        # Integrand = sqrt(K) * exp(-K) * (K-5/2)^{intj-1} * [lstar,m,star,nstar] * K^{3/2}
        # This function requires computing collisionality nu_D
        # We also spline interpolate lstar,mstar and star as function of cmul for each efield
        
        import matplotlib.pyplot as plt
        
        print('\n ############################################################')
        print('############# COMPUTING INTEGRANDS AS IN PENTA #################')
        print('############################################################')
        
        Kmin = 1e-4      #PENTA default value
        Kmax = 10        #PENTA default value
        numKsteps = 100  #PENTA default value
        K = np.linspace(Kmin,Kmax,numKsteps)
        
        cmul = {}
        for species in plasma_class.list_of_species:
            cmul_temp = []
            for k in K:
                vth = plasma_class.get_thermal_speed(species,self.roa)
                vparticle = vth * np.sqrt(k)
                nu = plasma_class.get_collisionality(species,self.roa,vparticle)
                cmul_temp.append( nu / vparticle )
            
            cmul[species] = cmul_temp
        
        # plot here cmul vs K for all species
        cmul_min = np.min(self.cmul)
        cmul_max = np.max(self.cmul)
        fig, ax = plt.subplots(figsize=(13,11))
        for species in plasma_class.list_of_species:
            plt.rc('font', size=18)
            plt.plot(K,cmul[f'{species}'],'o-',label=f'{species}')
            plt.plot(K,np.full_like(K,cmul_min),'-r',linewidth=2)
            plt.plot(K,np.full_like(K,cmul_max),'-r')
            ax.set_yscale('log')
            ax.set_ylabel(r'$\nu_D/v~~[m^{-1}]$')
            ax.set_xlabel(f'K')
            ax.legend()      
            ax.grid()
        #plt.show()

        
        # f_j(K)*K^(3/2)
        fix_func = np.sqrt(K) * np.exp(-K) * (K-2.5)**(intj-1) * K**1.5
        
        # fig, ax = plt.subplots(figsize=(8,6))
        # ax.plot(K,fix_func,'o-')
        # #ax.set_yscale('log')
        # ax.set_ylabel(r'$\sqrt{K}e^{-K}\left(K-5/2\right)^{j-1}\,K^{3/2}$')
        # ax.set_xlabel(f'K')   
        # ax.set_title(f'j={intj}')
        # ax.grid()    
        # plt.show()
        
        # # ok, let's make spline for lstar for each efield!
        from scipy.interpolate import interp1d #, splrep, BSpline
        
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
            lstar_interp1d = interp1d(xlog,yl,kind='quadratic',bounds_error=False,fill_value=0.0)
            nstar_interp1d = interp1d(xlog,yn,kind='quadratic',bounds_error=False,fill_value=0.0)
            logmstar_interp1d = interp1d(xlog,ylogm,kind='quadratic',bounds_error=False,fill_value=0.0)
            
            xspline = np.logspace(np.log10(x[0]),np.log10(x[-1]),100)
            
            fig, ax = plt.subplots(1,3,figsize=(17,6))
            ax[0].plot(x,yl,'ob')
            ax[0].plot(xspline, lstar_interp1d(np.log(xspline)),'red',label='spline')
            ax[0].set_yscale('log')
            ax[0].set_xscale('log')
            ax[0].set_ylabel(r'lstar')
            ax[0].set_xlabel(f'cmul')   
            ax[0].set_title(f'Er/v={efield}')
            ax[0].grid()
            ax[0].legend()
            
            ax[1].plot(x,yn,'ob')
            ax[1].plot(xspline, nstar_interp1d(np.log(xspline)),'red',label='spline')
            #ax.set_yscale('log')
            ax[1].set_xscale('log')
            ax[1].set_ylabel(r'nstar')
            ax[1].set_xlabel(f'cmul')   
            ax[1].set_title(f'Er/v={efield}')
            ax[1].grid()   
            #ax[1].legend()
            
            ax[2].plot(x,ylogm,'ob')
            ax[2].plot(xspline, logmstar_interp1d(np.log(xspline)),'red',label='spline')
            #ax.set_yscale('log')
            ax[2].set_xscale('log')
            ax[2].set_ylabel('ln(mstar)')
            ax[2].set_xlabel(f'cmul')   
            ax[2].set_title(f'Er/v={efield}')
            ax[2].grid()   
            #ax[1].legend()
            
            plt.tight_layout(pad=2)
            
            # full integrand of l*
            fig,ax = plt.subplots(figsize=(8,6))
            for species in plasma_class.list_of_species:
                ax.plot(K,lstar_interp1d(np.log(cmul[species]))*fix_func,'o-',label=f'{species}')       
                ax.set_xlabel('K')
                ax.set_ylabel(fr'$f_{intj}(K)~l^*(K)~K^{{3/2}}$')
                ax.set_title(f'Er/v={efield}')
            plt.legend()
            
            # full integrand of n*
            fig,ax = plt.subplots(figsize=(10,6))
            for species in plasma_class.list_of_species:
                ax.plot(K,nstar_interp1d(np.log(cmul[species]))*fix_func,'o-',label=f'{species}')       
                ax.set_xlabel('K')
                ax.set_ylabel(f'$f_{intj}(K)~n^*(K)~K^{{3/2}}$')
                ax.set_title(f'Er/v={efield}')
            plt.legend()
            
            # full integrand of m*
            fig,ax = plt.subplots(figsize=(10,6))
            for species in plasma_class.list_of_species:
                ax.plot(K,np.exp(logmstar_interp1d(np.log(cmul[species])))*fix_func,'o-',label=f'{species}')       
                ax.set_xlabel('K')
                ax.set_ylabel(f'$f_{intj}(K)~m^*(K)~K^{{3/2}}$')
                ax.set_title(f'Er/v={efield}')
                ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
            plt.legend()
            plt.show()
            
            # fig, ax = plt.subplots(1,3,figsize=(17,6))
            # for species in plasma_class.list_of_species:
            #     ax[0].plot(K,lstar_interp1d(np.log(cmul[species]))*fix_func,'o-',label=f'{species}')       
            #     ax[0].set_xlabel('K')
            #     ax[0].set_ylabel(fr'$f_{intj}(K)~l^*(K)~K^{{3/2}}$')
            #     ax[0].set_title(f'Er/v={efield}')
            #     ax[0].grid()
            #     #ax[0].legend()
            # for species in plasma_class.list_of_species:
            #     ax[1].plot(K,nstar_interp1d(np.log(cmul[species]))*fix_func,'o-',label=f'{species}')       
            #     ax[1].set_xlabel('K')
            #     ax[1].set_ylabel(f'$f_{intj}(K)~n^*(K)~K^{{3/2}}$')
            #     ax[1].set_title(f'Er/v={efield}')
            #     ax[1].grid()
            # for species in plasma_class.list_of_species:
            #     ax[2].plot(K,np.exp(logmstar_interp1d(np.log(cmul[species])))*fix_func,'o-',label=f'{species}')       
            #     ax[2].set_xlabel('K')
            #     ax[2].set_ylabel(f'$f_{intj}(K)~m^*(K)~K^{{3/2}}$')
            #     ax[2].set_title(f'Er/v={efield}')
            #     ax[2].grid()
            #     ax[2].legend()
            # plt.tight_layout(pad=2)
            # plt.show()        
        
        
        
        
        
        
        
    
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