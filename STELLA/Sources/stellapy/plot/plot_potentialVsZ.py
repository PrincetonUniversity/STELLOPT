 
#================================================================
# Functions to show the evolutions of omega in function of t 
#================================================================
      
# Load modules
import numpy as np
import matplotlib.pyplot as plt 
 
# Load personal modules
from stellapy.utils.decorators import verbose_wrapper 
from stellapy.plot.utils import load_plotbox2d
from stellapy.simulations.utils.get_simulation import get_simulation
from stellapy.simulations.utils.get_experiment import get_experiment
from stellapy.plot.utils import get_axisOfScan 
from stellapy.plot.utils import get_lineColors
 
@verbose_wrapper 
def plot_potentialVsZ(\
            # Specify which simulations to plot
                research=None,\
                experiment_id="All experiments",\
                simulation_id="All simulations",\
            # Specify data range
                units="normalized",\
                x_quantity="z",\
                y_quantity="phi_real",\
                x_range=None,\
                y_range=None,\
                plotted_modes="unstable",\
                kx_range=-0.0,\
                ky_range=[0,100],\
            # Labels
                x_label=None,\
                y_label=None,\
                title=None,\
            # For the GUI the figure object already exists 
                ax=None,\
                Progress=None,\
                root=None,\
            # Appearance options
                font_size=20,\
                handle_length=1):
    ''' Plot phi(z) to see if <nfield_periods> was big enough. 
 
    The z-axis can be displaced in zeta <zeta>, in z <z>, in poloidal turns <pol> and toroidal turns <tor>.
 
    Data files needed
    -----------------
    *.in
    *.out.nc or *.out.h5
    *.final_fields
     
    Parameters
    ----------
    units : {normalized; SI units}
    y_quantity : {phi2; phi_real; phi_imag}
    x_quantity : {zeta; z; pol}
 
    Notes
    -----
    Structure of the *.final_fields file
    # z   z-thet0   aky   akx   real(phi)   imag(phi)   real(apar)   imag(apar)   z_eqarc-thet0
 
    '''

    # Update the progress bar of the GUI
    if Progress: Progress.move(0,"Plot the potential versus z.")
                    
    # Decide whether to scan modes along the x-axis or y-axis
    scan, k_fixed, k_range = get_axisOfScan(kx_range, ky_range, root) 
    if scan == False: return 
    
    # Set the labels, and change the standard labels depending on the units
    label = determine_labels(x_label, y_label, title, units)
    load_plotbox2d(x_label=label[x_quantity], y_label=label[y_quantity], title=label["title"], ax=ax); 
    
    # Get a reference to the experiment and its simulations
    experiment  = get_experiment(research, experiment_id)[0]
    simulations = get_simulation(experiment, simulation_id)
    
    # The number of lines defines the colors so the colors of the lines change gradually with k
    color = get_lineColors(simulations, scan, k_fixed, k_range, plotted_modes); plot_i = 0

    # Save the axis limits
    xlims=[0,0]; ylims=[0,0]

    # Iterate over the simulations
    for simulation in simulations:
            
        # Get the modes of this simulation that need to be plotted
        vec_k = simulation.get_modesForAOneDimensionalScan(scan, k_fixed, k_range, plotted_modes) 
        
        # Get the z and potential data from the simulation
        z_data  = simulation.zeta if x_quantity=="zeta" else (simulation.z_kxky if x_quantity=="z" else simulation.z_poloidal)
        y_data  = simulation.phi_real if y_quantity=="phi_real" else (simulation.phi_imag if y_quantity=="phi_imag" else simulation.phi2)

        # Iterate over the modes
        for k in vec_k:
            
            # Update the progress bar of the GUI
            if Progress: Progress.move(vec_k.index(k)/len(vec_k)*100,"Plotting modes ("+str(vec_k.index(k))+"/"+str(len(vec_k))+")")

            # Labels for the legend
            simulation_id = simulation.marker_label.replace('_', ' ')
            if isinstance(kx_range, float): label = "$k_y\\rho_i = " + "{0:.2f}".format(k) + "$"
            if isinstance(ky_range, float): label = "$k_x\\rho_i = " + "{0:.2f}".format(k) + "$"
            if len(simulations) > 1:  label = label + ": " + simulation_id                
            if len(vec_k)==1: label = simulation_id
            
            # Get the z and potential data for this specific mode and remove the NaN and Inf values
            i_kx, i_ky = simulation.get_indicesOfMode(scan, k_fixed, k)
            vec_z = z_data[:,i_kx,i_ky][np.isfinite(y_data[:,i_kx,i_ky])]
            vec_y = y_data[:,i_kx,i_ky][np.isfinite(y_data[:,i_kx,i_ky])]
             
            # Plot the normalized potential since we are interested in the shape
            try: vec_y = vec_y/np.nanmax([np.nanmax(vec_y), abs(np.nanmin(vec_y))])
            except: print("Warning: The potential vector was empty.")
             
            # Plot omega(t) and gamma(t)
            ax.plot(vec_z, vec_y, lw=2, linestyle='-', color=color[plot_i], label=label, clip_on=False); plot_i += 1
            
            # Keep track of the axis limits
            xlims = [min(z_data.min(), xlims[0]), max(z_data.max(), xlims[1])]
            ylims = [min(y_data.min(), ylims[0]), max(y_data.max(), ylims[1])]
     
    # If there were modes to be plotted, rescale and show the figure.
    if xlims != [0,0] and ylims != [0,0]: 
        
        # Rescale the axis
        if x_range==None: x_range = xlims
        if y_range==None: y_range = [0,1] if y_quantity=="phi2" else [-1,1]
        rescale_axes(ax, x_range, y_range, units, font_size, handle_length)       

        # Show the figure if we're not in a GUI
        if not root: plt.show()  

    # If there were no modes, just clear the axis
    else:
        ax.clear()
        print("WARNING: There were no modes to be plotted.")
    
    # End
    if True: return     
        
#################################################################
#                        METHODS
#################################################################

 
def rescale_axes(ax, x_range, y_range, units, font_size, handle_length):
       
    # Change the range of the axes
    ax.set_xlim(x_range) 
    ax.set_ylim(y_range)
     
    # Change the ticker
    if units=="normalized":
        ax.ticklabel_format(style='plain', axis='x')
        ax.ticklabel_format(style='plain', axis='y')
    if units=="SI units": 
        ax.ticklabel_format(style='sci', scilimits=(0,0), axis='x')
        ax.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
         
    # Add the legend
    ax.legend(\
        loc='upper center', \
        bbox_to_anchor=(1.11, 1.03), \
        shadow=True, \
        ncol=1, \
        labelspacing=0.0, \
        prop={'size': font_size}, \
        handlelength = handle_length)
    return     

#----------------------------   
def determine_labels(x_label, y_label, title,  units):
    ''' Changes the labels if the units changed. '''
    
    label = {'z' : x_label, 'zeta' : x_label, 'pol' : x_label, 'title' : title,\
             'phi_real' : y_label, 'phi_imag' : y_label, 'phi2' : y_label}
    
    if label["z"] in [None, '$z/a$', '$z$ [m]']:   
        if units=="normalized": label["z"]  = '$z/a$' 
        if units=="SI units":   label["z"]  = '$z$ [m]'
    if label["zeta"] in [None, '$\\zeta$']:   
        if units=="normalized": label["zeta"]  = '$\\zeta$' 
        if units=="SI units":   label["zeta"]  = '$\\zeta$' 
    if label["pol"] in [None, '$Poloidal turns$' ]:   
        if units=="normalized": label["pol"] = '$Poloidal turns$' 
        if units=="SI units":   label["pol"] = '$Poloidal turns$' 
    if label["phi_real"] in [None, 'Re($\\phi)/Re($\\phi)_{max}$']:     
        if units=="normalized": label["phi_real"] = 'Re($\\phi)/Re($\\phi)_{max}$'
        if units=="SI units":   label["phi_real"] = 'Re($\\phi)/Re($\\phi)_{max}$'
    if label["phi_imag"] in [None, 'Im($\\phi)/Re($\\phi)_{max}$']:   
        if units=="normalized": label["phi_imag"] = 'Im($\\phi$)/Im($\\phi)_{max}$'
        if units=="SI units":   label["phi_imag"] = 'Im($\\phi$)/Im($\\phi)_{max}$'
    if label["phi2"] in [None, '$|\\phi|^2/|\\phi_{max}|^2$']:   
        if units=="normalized": label["phi2"] = '$|\\phi|^2/|\\phi_{max}|^2$'
        if units=="SI units":   label["phi2"] = '$|\\phi|^2/|\\phi_{max}|^2$'
    return label
         
         



