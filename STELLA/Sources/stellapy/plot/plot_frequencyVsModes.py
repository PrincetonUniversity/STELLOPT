
#================================================================
# Plot omega(k) or gamma(k)
#================================================================

# Load the modules
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.interpolate import interp1d

# Personal modules
from stellapy.plot.utils import load_plotbox2d 
from stellapy.utils.decorators import verbose_wrapper
from stellapy.simulations.utils.get_simulation import get_simulation
from stellapy.simulations.utils.get_experiment import get_experiment
from stellapy.plot.utils import get_axisOfScan 

#===========================
# plot omega(k) or gamma(k)
#===========================

@verbose_wrapper
def plot_frequencyVsModes(\
            # Specify which simulations to plot
                research=None,\
                experiment_id="All experiments",\
                simulation_id="All simulations",\
            # Specify data range
                y_quantity='omega',\
                x_range=None,\
                y_range=None,\
            # Details of the modes
                kx_range=-0.0, \
                ky_range=[0,100], \
                k_value=9.0, \
                k_delta=1.0, \
                plotted_modes="unstable",\
                lineardata="average",\
            # Labels
                x_label=None,\
                y_label=None,\
                title=None,\
            # For the GUI the figure object already exists 
                show_figure = False,\
                ax=None, \
                Progress=None,\
                root=None,\
            # Toggles
                units="normalized", \
                show_error=False, \
                maxima=False, \
                interpolate=False):
    ''' Plot the frequency/omega or growthrate/gamma versus the wave number. 
    
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
    lineardata: {"average", "last"}
        Either plot the last time value of omega/gamma or the average
        over the last 10% of time.
    '''


    # Update the progress bar of the GUI
    if Progress: Progress.move(0,"Plot omega or gamma versus time.")
         
    # Decide whether to scan modes along the x-axis or y-axis
    scan, k_fixed, k_range = get_axisOfScan(kx_range, ky_range, root) 
    if scan == False: return 
        
    # Set the labels, and change the standard labels depending on the units
    label = determine_labels(x_label, y_label, kx_range, ky_range, title,  units)
    load_plotbox2d(x_label=label['x'], y_label=label[y_quantity], title=label["title"], ax=ax) 

    # Keep track of the labels that are already used and save the axis limits
    labels = []; xlims=[0,0]; ylims=[0,0]; 
    
    # Get the experiments
    experiments = get_experiment(research, experiment_id)
    
    # Get the total number of varied values, to determine whether we need to use markers
    # to diferentiate between different simulations or whether the line colors are sufficient
    number_ofVariedValues = 0
    for experiment in experiments:
        number_ofVariedValues = np.max([number_ofVariedValues,len(experiment.variedValues)])
    
    # Iterate over the experiments
    for experiment in experiments:
        
        # Get the simulations
        simulations = get_simulation(experiment, simulation_id)
        
        # Iterate over the simulations
        for simulation in simulations:
            
            # Update the progress bar of the GUI
            if Progress: i = experiments.index(experiment)*len(simulations)+simulations.index(simulation); length=len(experiments)+len(simulations)
            if Progress: Progress.move(i/length*100,"Plotting spectra ("+str(i)+"/"+str(length)+")")
            
            # Get the modes of this simulation that need to be plotted and their indices in the (kx,ky) matrixes
            modes  = simulation.get_modesForAOneDimensionalScan(scan, k_fixed, k_range, plotted_modes) 
            i_kxky = simulation.get_indicesOfModes(scan, k_fixed, modes)  
            
            # Get the omega or gamma at the last time value
            if y_quantity=="omega" and lineardata=="average": y_data = simulation.omega_avg[i_kxky]
            if y_quantity=="omega" and lineardata=="last":    y_data = simulation.omega_last[i_kxky]
            if y_quantity=="gamma" and lineardata=="average": y_data = simulation.gamma_avg[i_kxky]
            if y_quantity=="gamma" and lineardata=="last":    y_data = simulation.gamma_last[i_kxky]
            if y_quantity=="gamma/ky**2" and lineardata=="average": y_data = simulation.gamma_avg[i_kxky]/[m**2 for m in modes]
            if y_quantity=="gamma/ky**2" and lineardata=="last":    y_data = simulation.gamma_last[i_kxky]/[m**2 for m in modes]
            
            # Initiate the errors on the calculation of omega/gamma at the last time step
            y_error = np.empty((2, len(modes))); y_error[:, :] = np.NaN
                
            # Get the error bars
            if y_quantity=="omega": y_error[0,:] = simulation.omega_min[i_kxky]
            if y_quantity=="omega": y_error[1,:] = simulation.omega_max[i_kxky]
            if y_quantity=="gamma": y_error[0,:] = simulation.gamma_min[i_kxky]
            if y_quantity=="gamma": y_error[1,:] = simulation.gamma_max[i_kxky]
            if y_quantity=="gamma/ky**2": y_error[0,:] = simulation.gamma_min[i_kxky]/[m**2 for m in modes]
            if y_quantity=="gamma/ky**2": y_error[1,:] = simulation.gamma_max[i_kxky]/[m**2 for m in modes]
                
            # Remove the Inf and Nan data points
            modes   = np.array(modes)[np.isfinite(y_data)]
            y_error = y_error[:, np.isfinite(y_data)]
            y_data  = y_data[np.isfinite(y_data)]
                
            # Plot the lines with the same color for each experiment if there are multiple experiments
            # or if we only varied 1 parameter since then we don't want to use any markers
            if len(modes)>0 and (len(experiments) > 1 or number_ofVariedValues==1):
                
                # Add a label to the legend for the first simulation of this experiment
                label_exp  = experiment.line_label if simulations.index(simulation) == 0 else ""

                # If we interpolate, don't plot the lines
                lw = 2 if not interpolate else 0

                # Plot the data
                ax.plot(-100, -100, lw=2, linestyle="-", color=experiment.line_color, label=label_exp)
                ax.plot(modes, y_data, lw=lw, linestyle="-", color=experiment.line_color)
                
            # Next plot the data points with the same color for each varied_value, to differentiate between simulations
            if len(modes)>0:
                
                # Add a label to the legend for the first simulation of this experiment
                label_var = simulation.marker_label
                label_var = label_var.replace("_"," ") if "$" not in label_var else label_var
                label_var = "" if (label_var in labels) else label_var; labels.append(label_var)
                marker_style = simulation.marker_style 
                marker_color = simulation.marker_color 
            
                # If there is only one experiment, add lines in the same color as the marker
                if not( (len(experiments) > 1) or (number_ofVariedValues==1) ): lw = 2
                if      (len(experiments) > 1) or (number_ofVariedValues==1)  : lw = 0
                if interpolate                                                : lw = 0
            
                # Plot the actual linear data points obtained by stella as markers (mfc = markerfacecolor; mec = markeredgecolor]   
                ax.plot(modes, y_data, lw=lw ,linestyle='-', color=marker_color, ms=5 ,\
                        marker=marker_style, mec=marker_color, mfc=marker_color, label=label_var)
    
                # Add the error bars
                if show_error==True:
                    print()
                    print("Modes", experiment.id)
                    print(y_error)
                    ax.errorbar(modes, y_data, yerr=y_error, fmt='o', color=simulation.marker_color, capsize=2)
    
                # Calculate the maximum growth rate
                if interpolate or maxima:
                    color = experiment.line_color if (len(experiments) > 1 or number_ofVariedValues==1) else marker_color 
                    plot_maximumGrowthrates(modes, y_data, y_quantity, ax, color, interpolate)
                
                # Keep track of the axis limits
                xlims = [np.nanmin([np.nanmin(modes), xlims[0]]), np.nanmax([1.1*np.nanmax(modes), xlims[1]])]
                ylims = [np.nanmin([np.nanmin(y_data), ylims[0]]), np.nanmax([1.1*np.nanmax(y_data), ylims[1]])]
                    
    # Shade ky_range in the plot to indicate where the reflectometry machine measures.
    if k_value and k_delta and units=="SI units":
        ax.axvspan(k_value-k_delta, k_value+k_delta, alpha=0.5, color='navy') 
        
    # Change appearance plot 
    ax.autoscale()
    ax.set_title(title)
    ax.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
    
    # Add the legend and sort its labels
    handles, labels = ax.get_legend_handles_labels()
    labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
    ax.legend(handles, labels, loc='best',labelspacing=0.0, prop={'size':20}, framealpha=0.9)
    
    # Rescale the axis
    if x_range==None: x_range = xlims
    if y_range==None: y_range = ylims
    ax.set_xlim(x_range) 
    ax.set_ylim(y_range) 
    
    # Show the figure
    if show_figure: plt.show()
    return 

#################################################################
#                        METHODS
#################################################################

def determine_labels(x_label, y_label, kx_range, ky_range, title,  units):
    ''' Set the labels and titels for the (omega, ki) or (gamma, ki) plot '''
    
    label = {'x' : x_label, 'omega' : y_label, 'gamma' : y_label, 'gamma/ky**2' : y_label, 'title' : title}
    
    if isinstance(ky_range, float) and label["x"] in [None, '$k_{x}\\rho_i$', '$k_{x}$ [cm$^{-1}$]']:
        if units=="normalized": label["x"] = '$k_{x}\\rho_i$'
        if units=="SI units":   label["x"] = '$k_{x}$ [cm$^{-1}$]'

    if isinstance(kx_range, float) and label["x"] in [None, '$k_{y}\\rho_i$', '$k_{y}$ [cm$^{-1}$]']:
        if units=="normalized": label["x"] = '$k_{y}\\rho_i$'
        if units=="SI units":   label["x"] = '$k_{y}$ [cm$^{-1}$]'

    if label["omega"] in [None, '$\\omega a/v_{\\mathrm{th},i}$', '$\\omega$ [s$^{-1}$]']:   
        if units=="normalized": label["omega"] = '$\\omega a/v_{\\mathrm{th},i}$'
        if units=="SI units":   label["omega"] = '$\\omega$ [s$^{-1}$]'

    if label["gamma"] in [None, '$\\gamma a/v_{\\mathrm{th},i}$', '$\\gamma$  [s$^{-1}$]']:     
        if units=="normalized": label["gamma"] = '$\\gamma a/v_{\\mathrm{th},i}$'
        if units=="SI units":   label["gamma"] = '$\\gamma$  [s$^{-1}$]'

    if label["gamma/ky**2"] in [None, '$\\gamma a/v_{\\mathrm{th},i} * (1/k_y^2)$', '$\\gamma/k_y^2$  [s$^{-1}$]']:     
        if units=="normalized": label["gamma/ky**2"] = '$\\gamma a/v_{\\mathrm{th},i} * (1/k_y^2)$'
        if units=="SI units":   label["gamma/ky**2"] = '$\\gamma/k_y^2$  [s$^{-1}$?]'
        
    return label
        
        
#---------------------------------------------
def plot_maximumGrowthrates(x_data, y_data, y_quantity, ax, color, interpolate):
    
    # Calculate the maxima after interpolation
    if interpolate:
        
        # First interpolate the data
        xnew = np.linspace(np.nanmin(x_data), np.nanmax(x_data), num=100, endpoint=True)
        function = interp1d(x_data, y_data, kind='cubic')
        
        # Plot the interpolation
        ax.plot(xnew, function(xnew), lw=2 ,linestyle='-', color=color) 
        
        # Now find and plot the maxima
        y = function(xnew)
        peaks, _ = find_peaks(y)
        ax.plot(xnew[peaks], y[peaks],lw=0,marker='*',mec='black',mfc='gold',label="",ms=12,zorder=10)

    # Calculate the maxima directly
    if not interpolate:
        x_data = list(x_data); y_data = list(y_data)
        maxima = [y for y in y_data[1:-1] if (y>y_data[y_data.index(y)-1] and y>y_data[y_data.index(y)+1]) ]
        modes =  [x_data[y_data.index(m)] for m in maxima]
        sorted_indexes = np.array(maxima).argsort()
        maxima = [maxima[i] for i in sorted_indexes]
        modes  = [modes[i]  for i in sorted_indexes]
        ax.plot(modes,maxima,lw=0,marker='*',mec='black',mfc='gold',label="",ms=12,zorder=10)
    
    

