"""
This library provides a python class for reading and handling DKES
results data. It also prepares the input PENTA coefficients
"""

import numpy as np
import sys

# Constants
EC = 1.602176634E-19 # Electron charge [C]

# VMEC Class
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
        
        #Check convergence of coefficients
        self.check_convergence(self.L11m,self.L11p,'L11')
        self.check_convergence(self.L31m,self.L31p,'L13')
        self.check_convergence(self.L33m,self.L33p,'L33')
        
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
            print(' ')
            print('RELATIVE ERROR OF '+Lvar+' >',self.eps_rel*100,'%')
            for idx in idx_eps:
                print(f'{"cmul":<12} {"efield":<12} {Lvarm:<12} {Lvarp:<12} {"rel-error (%)":<15}')
                print(f'{self.cmul[idx]:<12.3E} {self.efield[idx]:<12.3E} {self.L11m[idx]:<12.4E} {self.L11p[idx]:<12.4E} {self.rel_error_L11[idx]*100:<15.2f}')
            print(f'Making the average anyway: {Lvar} = 0.5*({Lvarp} + {Lvarm})')

      

     


  
# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)