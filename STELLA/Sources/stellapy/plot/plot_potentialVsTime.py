 
#================================================================
# Plot flux(t) for non-linear runs
#================================================================
# TODO: the potential is not normalized if units="normalized"

# Load the modules
import numpy as np
import matplotlib.pyplot as plt

# Personal modules
from stellapy.plot.utils import load_plotbox2d  
from stellapy.utils.decorators import verbose_wrapper
from stellapy.simulations.utils.get_simulation import get_simulation
from stellapy.simulations.utils.get_experiment import get_experiment 
 
@verbose_wrapper 
def plot_potentialVsTime(\
            # Specify which simulations to plot
                research=None,\
                experiment_id="All experiments",\
                simulation_id="All simulations",\
            # Specify data range
                x_range=None,\
                y_range=None,\
            # Labels
                x_label=None,\
                y_label=None,\
                title=None,\
            # For the GUI the figure object already exists 
                show_figure = True,\
                ax=None,\
                Progress=None,\
            # Toggles
                log=False,\
                units="normalized"):
    
    ''' Plot potential(t).  '''

    # Update the progress bar of the GUI
    if Progress: Progress.move(0,"Plot potential versus time.")
     
    # Set the labels, and change the standard labels depending on the units
    label = determine_labels(x_label, y_label, title, units)
    load_plotbox2d(x_label=label['x'], y_label=label["y"], title=label["title"], ax=ax) 

    # Keep track of the labels that are already used and save the axis limits
    labels = []; xlims=[0,0]; ylims=[0,0]; 
    
    # Get the total number of varied values, to determine whether we need to use markers
    # to diferentiate between different simulations or whether the line colors are sufficient
    number_ofVariedValues = 0
    experiments = get_experiment(research, experiment_id)
    for experiment in experiments:
        number_ofVariedValues = np.max([number_ofVariedValues,len(experiment.variedValues)])
     
    # Iterate over the experiments
    for experiment in experiments:
                     
        # First plot the lines with the same color for each experiment
        if (len(experiments) > 1) or (number_ofVariedValues==1):
            
            # Iterate over the simulations
            simulations = get_simulation(experiment, simulation_id)
            for simulation in simulations:
                         
                # Get the data
                vec_time = simulation.time
                vec_phi2 = simulation.phi2
                style    = '-'
                
                # Put time time axis in milliseconds if we plot in SI units
                if units!="normalized": vec_time=vec_time*1000  
                 
                # Add a label to the legend for the first simulation of this experiment
                label_exp = experiment.line_label if simulations.index(simulation) == 0 else ""

                # Plot the data
                ax.plot(vec_time, vec_phi2, lw=3, linestyle=style, color=experiment.line_color, label=label_exp) 
 
    # Plot the data points with the same color for each varied_value
    for experiment in experiments:
        
            # Iterate over the simulations
            simulations = get_simulation(experiment, simulation_id)
            for simulation in simulations:
                     
                # Get the data
                vec_time = simulation.time
                vec_phi2 = simulation.phi2
                style    = '-'
                
                # Put time time axis in milliseconds if we plot in SI units
                if units!="normalized": vec_time=vec_time*1000  
                 
                # Add a label to the legend for the first simulation of this experiment
                label_var = simulation.marker_label
                label_var = label_var.replace("_"," ") if "$" not in label_var else label_var
                label_var = "" if (label_var in labels) else label_var; labels.append(label_var)
                label_var1 = label_var if not ((len(experiments) > 1) or (number_ofVariedValues==1)) else ""
                label_var2 = label_var if (len(experiments) > 1) or (number_ofVariedValues==1) else ""
                marker_style = simulation.marker_style 
                marker_color = simulation.marker_color 
                 
                # If there is only one experiment, add lines in the same color as the marker
                if not( (len(experiments) > 1) or (number_ofVariedValues==1) ): lw = 2
                if      (len(experiments) > 1) or (number_ofVariedValues==1)  : lw = 0
 
                # First print the lines with the shot label: for ions and electrons
                ax.plot(vec_time, vec_phi2, lw=lw, linestyle=style, color=marker_color, label=label_var1) 
 
                # Keep track of the axis limits
                xlims = [np.nanmin([np.nanmin(vec_time), xlims[0]]), np.nanmax([1.1*np.nanmax(vec_time), xlims[1]])]
                ylims = [np.nanmin([np.nanmin(vec_phi2), ylims[0]]), np.nanmax([1.1*np.nanmax(vec_phi2), ylims[1]])]
                
                # Next print the symbols/markers with the rho label, print the markers differently for log and normal scales
                if (len(experiments) > 1) or (number_ofVariedValues==1):
                    vec_flux, vec_time = get_markerAtEveryXSteps(vec_phi2, vec_time, log, units)
                    ax.plot(vec_time, vec_flux, lw=0, marker=marker_style, mec=marker_color, mfc=marker_color, label=label_var2) 
 
    # Change appearance plot 
    ax.autoscale()
    ax.set_title(title)
    ax.set_xlim(xmin=0)
    ax.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
    ax.legend(labelspacing=0.0, prop={'size':20}, handlelength=1.5)
    if log == True: ax.set_yscale('log')
    
    # Rescale the axis
    if x_range==None: x_range = xlims
    if y_range==None: y_range = ylims
    ax.set_xlim(x_range) 
    ax.set_ylim(y_range) 
    
    # Show the figure
    if show_figure: plt.show()
    if True: return


#################################################################
#                        METHODS
#################################################################

def determine_labels(x_label, y_label, title, units):
    ''' Set the labels and titels for the potential versus time plot.'''
    
    # Save the labels in one dictionary
    label = {'x' : x_label, 'y' : y_label, 'title' : title}

    # Determine the label of the x-axis
    label["x"] = "$t\, v_{\mathrm{th}}/a$" if units=="normalized" else "$t$ [ms]"
    label["y"] = "$|\\phi^2|$" if units=="normalized" else "$|\\phi^2|$ ???"
    return label
 
 
#================================================================
def get_markerAtEveryXSteps(vec_flux, vec_time, log, units):
 
    # For normal scales, plot markers after the exponentional growth, for log scales, plot markers on the exponentional growth
    if log == True:
        if units=="normalized":  time_step = 15
        if units!="normalized":  time_step = 0.05 
        time_of_maximum = vec_time[np.argmax(vec_flux)]
        vec_flux = vec_flux[vec_time <= time_of_maximum]
        vec_time = vec_time[vec_time <= time_of_maximum]
    if log == False:  
        if units=="normalized":  time_step = 25
        if units!="normalized":  time_step = 0.05 
        time_of_maximum = vec_time[np.argmax(vec_flux)]
        vec_flux = vec_flux[vec_time >= time_of_maximum]
        vec_time = vec_time[vec_time >= time_of_maximum]
 
    # Plot a marker at t=time_step, 2*time_step, ...
    marker_times = np.arange(vec_time[0],vec_time[-1]+time_step, time_step)
    marker_times[-1] = marker_times[-1] + 1000
    temp_flux = []
    temp_time = []
    index = 0
    for time in vec_time:
        if time >= marker_times[index]:
            temp_time.append(time)
            temp_flux.append(vec_flux[list(vec_time).index(time)])
            index = index + 1
 
    return np.array(temp_flux), np.array(temp_time)


