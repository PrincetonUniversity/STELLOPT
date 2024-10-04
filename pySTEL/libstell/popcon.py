"""
This library provides a python class for creating stellarator POPCON plots
"""

import numpy as np
import sys

# Constants
EC = 1.602176634E-19 # Electron charge [C]

# popcon Class
class POPCON:
    
    def __init__(self, B, a, R, iota, n_exp, T_exp, n0_array, T0_array):
        # no_array given in units of m^-3
        # T0_array is units of eV
        
        from scipy import integrate
        
        # Define tauiss04
        self.tauiss04 = lambda P,n: 0.134*a**2.28*R**0.64*P**-0.64*n**0.54*B**0.84*iota**0.41
        
        # Differential volume. This works well for large aspect-ratio; not so good otherwise
        self.dVdrho = lambda rho: 4*np.pi*np.pi*R*rho*a*a
        self.Volume, _ = integrate.quad(lambda x: self.dVdrho(x),0.0,1.0)
        
        print(rf'Volume: {self.Volume}')
        
        self.n_avg = self.get_averaged_from_profile(n0_array,0.01*n0_array,n_exp)
        self.T_avg = self.get_averaged_from_profile(T0_array,0.01*T0_array,T_exp)
        
        self.N, self.T = np.meshgrid(self.n_avg, self.T_avg)
        
        # sets self.RHS, which is a 2D array w/ same shape as n,T meshgrid
        # currently, RHS has:
        # - Bremstrahlung radiation losses
        # - Conduction losses based on ISS04
        # - alpha heating
        self.set_RHS()
        
        # finds Pext for each point of self.RHS
        self.find_Pext()      

        # Plot popcon
        self.plot_popcon()
        
    def get_averaged_from_profile(self,peak_vals,edge_vals,exponent):
        
        from scipy import integrate
        
        #check that length of(peak_vals) == length of(edge_vals)
        if( len(peak_vals) != len(edge_vals) ):
            print('ERROR: array of peak vals does not have the same size as edge_vals')
            exit(0)
            
        avg = []
        
        for peak,edge in zip(peak_vals,edge_vals):
                
            integral_top, _ = integrate.quad(lambda x: self.dVdrho(x)*(edge + (peak-edge)*(1-x**exponent)),0.0,1.0)

            avg.append( integral_top / self.Volume )
        
        print(avg)
        return np.array(avg)
    
    def set_RHS(self):
        
        # Bremstrahlung
        RHS_1 = self.get_Bremsstrahlung()
        
        
        
        




# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)