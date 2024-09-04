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
    
    def __init__(self, eps_rel=0.03):
        self.eps_rel = eps_rel
        
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
            
            
    def compute_PENTA_coeffs(self,wout_filename,surface):
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
        self.Bsq = Bsq_full[surface-1]
        print(f'Bsq={self.Bsq}')
        
        #compute PENTA lstar
        aux = 1 - 1.5*self.cmul*self.L33/self.Bsq
        self.lstar = self.L11 - (2./3.)*self.cmul*self.Usq + (1.5*self.cmul*self.L31*self.L31/self.Bsq)/aux
        self.lstar = self.lstar / (EC*EC)
        
        #compute PENTA mstar
        self.mstar = self.cmul*self.cmul*self.L33 / aux
        
        #compute PENTA nstar
        self.nstar = self.cmul*self.L31 / aux
        self.nstar = self.nstar / EC        
     
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
    
    #def write_PENTA_coeffs_to_file(self)
    
    #### computes integrands assuming an hydrogen plasma!!!
    def plot_PENTA_integrands_energy_conv(self,intj,plasma_class,rho):
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
                vth = plasma_class.get_thermal_speed(species,rho)
                vparticle = vth * np.sqrt(k)
                nu = plasma_class.get_collisionality(species,rho,vparticle)
                cmul_temp.append( nu / vparticle )
            
            cmul[species] = cmul_temp
        
        # plot here cmul vs K for all species
        fig, ax = plt.subplots(figsize=(8,6))
        for species in plasma_class.list_of_species:
            plt.rc('font', size=18)
            plt.plot(K,cmul[f'{species}'],'o-',label=f'{species}')
            ax.set_yscale('log')
            ax.set_ylabel(r'$\nu_D/v~~[m^{-1}]$')
            ax.set_xlabel(f'K')
            ax.legend()      
            ax.grid()
        plt.show()

        
        # plot f_j(K)*K^(3/2)
        fix_func = np.sqrt(K) * np.exp(-K) * (K-2.5)**(intj-1) * K**1.5
        
        fig, ax = plt.subplots(figsize=(8,6))
        ax.plot(K,fix_func,'o-')
        #ax.set_yscale('log')
        ax.set_ylabel(r'$\sqrt{K}e^{-K}\left(K-5/2\right)^{j-1}\,K^{3/2}$')
        ax.set_xlabel(f'K')   
        ax.set_title(f'j={intj}')
        ax.grid()    
        plt.show()
        
        # # ok, let's make spline for lstar for each efield!
        from scipy.interpolate import interp1d, splrep, BSpline
        from matplotlib.ticker import FuncFormatter
        
        for i in range(self.nefield):
            i1 = i*self.ncmul
            i2 = i1 + self.ncmul
            
            efield = self.efield[i1]
            
            x = self.cmul[i1:i2]
            y = self.lstar[i1:i2]
            
            # quadratic spline as in PENTA. Assuming log_interp = true
            xlog = np.log(x)
            lstar_interp1d = interp1d(xlog,y,kind='quadratic',bounds_error=False,fill_value=0.0)
            #lstar_Bspline = splrep(x,y,s=len(x)+100,k=2)
            
            xspline = np.logspace(-5,2,100)
            
            fig, ax = plt.subplots(figsize=(8,6))
            ax.plot(x,y,'ob')
            ax.plot(xspline, lstar_interp1d(np.log(xspline)),'red')
            #ax[0,0].plot(xspline, BSpline(*lstar_Bspline)(xspline,extrapolate=False),'green')
            ax.set_yscale('log')
            ax.set_xscale('log')
            ax.set_ylabel(r'lstar')
            ax.set_xlabel(f'cmul')   
            ax.set_title(f'spline, Er/v={efield}')
            ax.grid()    
            #plt.show()
            
            # plot lstar and full integrand as a function of K for all species
            # num_of_species = len(plasma_class.list_of_species)
            # fig,ax = plt.subplots(num_of_species,2)
            # plt.rc('font', size=16)
            # for species, iss in zip(plasma_class.list_of_species,np.arange(num_of_species)):
            #     ax[iss,0].plot(K,lstar_interp1d(np.log(cmul[species])),'ob-',label=f'Er/v={efield}')
            #     ax[iss,0].set_ylabel(r'lstar')
            #     ax[iss,0].set_xlabel(f'K')   
            #     ax[iss,0].set_title(f'{species}')
            #     ax[iss,0].grid()
                
            #     ax[iss,1].plot(K,lstar_interp1d(np.log(cmul[species]))*fix_func,'ob-') 
            #     ax[iss,1].set_xlabel('K')
            #     ax[iss,1].set_ylabel(r'$f_j(K)~l^*(K)~K^{3/2}$') 
            #     ax[iss,1].set_title(f'{species}, Er/v={efield}')
            #     ax[iss,1].grid()
            # plt.tight_layout
            # plt.show()
            
            # full integrand in the same plot
            fig,ax = plt.subplots(figsize=(8,6))
            for species in plasma_class.list_of_species:
                ax.plot(K,lstar_interp1d(np.log(cmul[species]))*fix_func,'o-',label=f'{species}')       
                ax.set_xlabel('K')
                ax.set_ylabel(f'$f_{intj}(K)~l^*(K)~K^{3/2}$')
                ax.set_title(f'Er/v={efield}')
            plt.legend()
            plt.show()
        # when it comes to mstar don't forget to take the log as in PENTA !!!!!
    
        
        
        
        
        
        
        
        
    
    def perp_coll_freq(self,n,m,Z,T,v,species_num,electron_num):
        """Calc perpendicular collision frequency of species species_num

		This routine takes a 1D arrays of n,m,Z,T,v which contains the
        density, mass, charge number, Temperature and velocity of all species
        in a plasma. Then it computes nu_D of species indicated by the 
        number species_num, which refers to the index in the arrays of 
        the species we are computing the collision frequency
        
        This is the collision frequency of the pitch-angle scatering
        operator and is defined in 
        S. P. Hirshman and D. J. Sigmar, Nucl. Fusion 21, 1079 (1981)

		Parameters
		----------
		n : list
			density in m^-3
        m : list
			mass in kg
        Z : list
			charge number; if electron, Z=-1
        T : list
			temperature in eV
        v : list
			velocity in m/s
        species num : integer
            identifies the species we want the collision freq of: nu_D_num_index
        electron_num : integer
            identifies the index of the electrons
		
        Returns
		----------
		nu_D : real
			Perpendicular collision frequency [s^-1]
		"""
  
        from scipy.special import erf
        
        EC = 1.602176634E-19 # Electron charge [C]
        EPS0 = 8.8541878188E-12 # Vacuum permittivity [F/m]
          
        Za = Z[species_num]
        ma = m[species_num]
        va = v[species_num]
        
        ne = n[electron_num]
        Te = T[electron_num]
        
        #compute loglambda as in PENTA
        if(Te>50):
            loglambda = 25.3 - 1.15*np.log10(ne/1e6) + 2.3*np.log10(Te)
        else:
            loglambda = 23.4 - 1.15*np.log10(ne/1e6) + 3.45*np.log10(Te)
        
        nu_D = 0.0
        for nb,mb,Zb,Tb,vb in zip(n,m,Z,T,v):
            
            vthb = np.sqrt(2*EC*Tb/mb)
            xb = vb/vthb
            
            Hb = (1-0.5/(xb*xb))*erf(xb) + np.exp(-xb*xb)/(xb*np.sqrt(np.pi))
            
            num = nb*(EC**4)*Za*Za*Zb*Zb*loglambda*Hb
            den = ma*ma*va*va*va*4*np.pi*EPS0*EPS0
            
            nu_D = nu_D + num/den
            
        return nu_D
            
            
        
        
        
        
        

      

     


  
# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)