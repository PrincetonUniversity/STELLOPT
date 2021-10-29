
#================================================================
# Functions to show the evolutions of omega in function of t 
#================================================================
     
# Load modules
import numpy as np
import matplotlib.pyplot as plt

# Personal modules
from stellapy.plot.utils import load_plotbox2d
from stellapy.utils.decorators import verbose_wrapper
from stellapy.simulations.utils.get_simulation import get_simulation
from stellapy.simulations.utils.get_experiment import get_experiment
from stellapy.plot.utils import get_axisOfScan 
from stellapy.plot.utils import get_lineColors

@verbose_wrapper 
def plot_frequencyVsTime(\
            # Specify which simulations to plot
                research=None,\
                experiment_id="All experiments",\
                simulation_id="All simulations",\
            # Specify data range
                x_range=None,\
                y_range=None,\
                quantity="omega",\
                units="normalized",\
                kx_range=-0.0,\
                ky_range=[0,100],\
                plotted_modes="unstable",\
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
    ''' Return data and plot *.omega file: time  ky  kx  Re[om]  Im[om]  Re[omavg] Im[omavg] 
    
    Parameters
    ----------
    research: object
        Contains the selected experiments and simulations.
    experiment_id, simulation_id: str
        Is used to select which simulations/experiments to plot.
    plotted_modes: {"unstable", "stable", "unconverged", "converged"}
        Determines which modes to plot.
    units: {"normalized", "SI units"}
        Determines the units of the data.
    '''
    
    # Update the progress bar of the GUI
    if Progress: Progress.move(0,"Plot omega or gamma versus time.")
                    
    # Decide whether to scan modes along the x-axis or y-axis
    scan, k_fixed, k_range = get_axisOfScan(kx_range, ky_range, root) 
    if scan == False: return 
    
    # Set the labels, and change the standard labels depending on the units
    label = determine_labels(x_label, y_label, title, units)
    load_plotbox2d(x_label=label['time'], y_label=label[quantity], title=label["title"], ax=ax); 
    
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
                
        # Iterate over the modes
        for k in vec_k:
            
            # Update the progress bar of the GUI
            if Progress: Progress.move(vec_k.index(k)/len(vec_k)*100,"Plotting modes ("+str(vec_k.index(k))+"/"+str(len(vec_k))+")")
            
            # Get the indices of the mode
            i_kx, i_ky = simulation.get_indicesOfMode(scan, k_fixed, k)

            # Labels for the legend
            simulation_id = simulation.marker_label.replace('_', ' ')
            if isinstance(kx_range, float): label = "$k_y\\rho_i = " + "{0:.2f}".format(k) + "$"
            if isinstance(ky_range, float): label = "$k_x\\rho_i = " + "{0:.2f}".format(k) + "$"
            if len(simulations) > 1:  label = label + ": " + simulation_id                
            if len(vec_k)==1: label = simulation_id
            
            # Get the time axis of the mode
            time  = simulation.time_kxky[:,i_kx,i_ky]
            ydata = simulation.omega_kxky[:,i_kx,i_ky] if quantity=="omega" else simulation.gamma_kxky[:,i_kx,i_ky]
            
            # Only plot as long as omega and gamma are not NaN
            time  = time[np.isfinite(ydata)]
            ydata = ydata[np.isfinite(ydata)]
            
            # Denormalize the data if necessary
            if units=="SI units": ref = { "length" : simulation.referenceUnits["length"], "vthermal" : simulation.referenceUnits["vthermal"]}     
            if units!="SI units": ref = { "length" : 1, "vthermal" : 1}  
            time  = time*ref["length"]/ref["vthermal"]
            ydata = ydata/ref["length"]*ref["vthermal"]
            
            # Plot omega(t) and gamma(t)
            ax.plot(time, ydata, lw=2, linestyle='-', color=color[plot_i], label=label); plot_i += 1
            
            # Keep track of the axis limits
            xlims = [min(time.min(), xlims[0]), max(time.max(), xlims[1])]
            ylims = [min(abs(ydata).min(), ylims[0]), max(1.5*abs(ydata[-1]).max(), ylims[1])]           
                        
    
    # If there were modes to be plotted, rescale and show the figure.
    if xlims != [0,0] and ylims != [0,0]:
        
        # Rescale the axis
        if x_range==None: x_range = xlims
        if y_range==None: y_range = ylims
        rescale_axes(ax, x_range, y_range, units, font_size, handle_length, simulation_id)  

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

def rescale_axes(ax, x_range, y_range, units, font_size, handle_length, simulation_id):
      
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
        
    # Add the legend for the labels, make more room for the long labels
    if simulation_id!="All simulations": bbox_to_anchor=(1.15, 1.03)
    if simulation_id=="All simulations": bbox_to_anchor=(1.2, 1.03)
    ax.legend(\
        loc='upper center', \
        bbox_to_anchor=bbox_to_anchor, \
        shadow=True, \
        ncol=1, \
        labelspacing=0.0, \
        prop={'size': font_size}, \
        handlelength = handle_length)
    return

#---------------------------------------
def determine_labels(x_label, y_label, title,  units):
    ''' Changes the labels if the units changed. '''
    
    label = {'time' : x_label, 'omega' : y_label, 'gamma' : y_label, 'title' : title}
    if label["time"] in [None, '$t\, v_{\mathrm{th}}/a$', '$t$ [ms]']:   
        if units=="normalized": label["time"]  = '$t\, v_{\mathrm{th}}/a$' 
        if units=="SI units":   label["time"]  = '$t$ [ms]'
    if label["omega"] in [None, '$\\omega\, a/v_{\mathrm{th}}$', '$\\omega\, [s$^{-1}$]']:   
        if units=="normalized": label["omega"] = '$\\omega\, a/v_{\mathrm{th}}$'
        if units=="SI units":   label["omega"] = '$\\omega$ [s$^{-1}$]'
    if label["gamma"] in [None, '$\\gamma\, a/v_{\mathrm{th}}$', '$\\gamma$ [s$^{-1}$]']:     
        if units=="normalized": label["gamma"] = '$\\gamma\, a/v_{\mathrm{th}}$'
        if units=="SI units":   label["gamma"] = '$\\gamma$ [s$^{-1}$]'
    return label
        
    
    
    
    
