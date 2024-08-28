"""
This library provides a python class for reading and handling DKES
results data. It also prepares the input PENTA coefficients
"""

import numpy as np

# Constants
EC = 1.602176634E-19 # Electron charge [C]

# VMEC Class
class DKES:
    
    def __init__(self, eps_rel=0.03):
        self.eps_rel = eps_rel
        
    def read_DKES_results(self,filename):
        
        dkes = np.loadtxt('results.GIGA_v120',skiprows=2)

       

     


  
# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)