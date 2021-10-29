
import numpy as np
import matplotlib.pyplot as plt

def get_lineColors(simulations, scan, k_fixed, k_range, plotted_modes):
    ''' The number of lines defines the colors so the colors of the lines change gradually with ky. '''
    
    # Initialize the number of colors
    dim_colors = 0
    
    # Go through the simulations to see how many modes are plotted
    for simulation in simulations:

        # Divide the modes in stable/unstable and converged/unconverged and get the modes that will actually be plotted
        vec_k = simulation.get_modesForAOneDimensionalScan(scan, k_fixed, k_range, plotted_modes) 
        dim_colors += len(vec_k) 
        
    # Define the colors based on the number of modes
    color = plt.cm.jet(np.linspace(0,1,dim_colors)) #@undefinedvariable
    return color
