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
        
        # get ncmul, nefield
        self.ncmul = np.sum(self.efield == self.efield[0])
        self.nefield = int(np.round(len(self.cmul)/self.ncmul,decimals=0))
        
        #Check convergence of coefficients
        self.check_convergence(self.L11m,self.L11p,'L11')
        self.check_convergence(self.L31m,self.L31p,'L13')
        self.check_convergence(self.L33m,self.L33p,'L33')
        
        #Take average of coefficients
        self.L11  = 0.5*(self.L11m+self.L11p)
        self.L31  = 0.5*(self.L31m+self.L31p)
        self.L33  = 0.5*(self.L33m+self.L33p)
        
    def check_convergence(self,Lm,Lp,Lvar):
        # checks if |Lm-Lp|/Lp < epsilon
        # if not, gives a warning in the command line
        # Lvar is 'L11', 'L13' or 'L33'
        
        if Lvar not in {'L11', 'L13', 'L33'}:
            print('Warning: Lvar can only be of the type L11, L13, or L33')
            sys.exit(1)
        else:
            Lvarm = Lvar+'m'
            Lvarp = Lvar+'p'
        
        if Lvar=='L11':
            self.rel_error_L11 = abs(self.L11m-self.L11p) / self.L11p
        elif Lvar=='L31':
            self.rel_error_L31 = abs(self.L31m-self.L31p) / self.L31p
        elif Lvar=='L33':
            self.rel_error_L33 = abs(self.L33m-self.L33p) / self.L33p
            
        #give indexes where rel error larger than epsilon
        idx_eps = np.where(self.rel_error_L11 > self.eps_rel)[0]
        
        #in case there is any index, print in the command line where this happens
        if idx_eps.size > 0:
            print('\nRELATIVE ERROR OF '+Lvar+' >',self.eps_rel*100,'%')
            for idx in idx_eps:
                print(f'{"cmul":<12} {"efield":<12} {Lvarm:<12} {Lvarp:<12} {"rel-error (%)":<15}')
                print(f'{self.cmul[idx]:<12.3E} {self.efield[idx]:<12.3E} {self.L11m[idx]:<12.4E} {self.L11p[idx]:<12.4E} {self.rel_error_L11[idx]*100:<15.2f}')
            print(f'Making the average anyway: {Lvar} = 0.5*({Lvarp} + {Lvarm})')
            
    def compute_PENTA_coeffs(self,wout_filename,surface):
        # Computes PENTA input coefficients lstar, mstar and nstar
        # size of lstar, mstar, nstar is n_cmul x n_efield
        # Check DKES/PENTA documentation to see the definition of lstar, mstar, nstar
        # In the documentation, D_ij^* corresponds to self.Lij, which are the species-independent DKES coefficientes
        print('\n#############################################################################')
        print('###################   Computing PENTA coefficients  #########################')
        print('#############################################################################')
        
        ######  WARNING: as of now this assumes an hydrogen plasma, qa=e_charge  #####
        print('\nWARNING: THIS ASSUMES AN HYDROGEN PLASMA, qa=echarge')
        
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
        
        from libstell.vmec import VMEC
        #interpolate Bsq from half to full grid
        Bsq_full = VMEC.h2f(self,Bsq_half)
        
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
        

        
        
        
        
        

      

     


  
# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)