 
#================================================================
# Plot flux(t) for non-linear runs
#================================================================
# TODO: Fluxes aren't normalized.
 
# Load the modules
import numpy as np
import matplotlib.pyplot as plt

# Personal modules
from stellapy.plot.utils import load_plotbox2d  
from stellapy.utils.decorators import verbose_wrapper
from stellapy.simulations.utils.get_simulation import get_simulation
from stellapy.simulations.utils.get_experiment import get_experiment 
 
@verbose_wrapper 
def plot_fluxesVsTime(\
            # Specify which simulations to plot
                research=None,\
                experiment_id="All experiments",\
                simulation_id="All simulations",\
                sat_flux=None,\
            # Specify data range
                species=[0],\
                y_quantity='q_flux',\
                x_range=None,\
                y_range=None,\
                t_range=None,\
            # Labels
                x_label=None,\
                y_label=None,\
                title=None,\
            # For the GUI the figure object already exists 
                show_figure = True,\
                ax=None,\
                Progress=None,\
            # Toggles
                show_restartTimes=False,\
                show_vspan=True,\
                normalize=False,\
                log=False,\
                units="normalized"):
    
    ''' Plot heat_flux(t) for each shot 
     
 
    Parameters
    ----------
    t_range :  dict[experiment.id][simulation.id][specie] = tuple of floats
    sat_flux : dict[experiment.id][simulation.id][specie] = float
     
    '''

    # Update the progress bar of the GUI
    if Progress: Progress.move(0,"Plot flux versus time.")
     
    # Set the labels, and change the standard labels depending on the units
    label = determine_labels(x_label, y_label, title, units, species, normalize)
    load_plotbox2d(x_label=label['x'], y_label=label[y_quantity], title=label["title"], ax=ax) 

    # Keep track of the labels that are already used and save the axis limits
    labels = []; xlims=[0,0]; ylims=[0,0]; 
    
    # Plot the electron and ion labels
    if len(species)>1:
        ax.plot(-1, -1, lw=3, linestyle='-', color='black', label='Ions') 
        ax.plot(-1, -1, lw=3, linestyle=':', color='black', label='Electrons') 
    
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
                
                # Iterate over the species
                for specie in species:
                         
                    # Get the data
                    vec_time = simulation.vec_time
                    vec_flux = getattr(simulation, y_quantity)[:,specie]
                    style    = ':' if (specie=='1' and len(species)>1) else '-'
                    
                    # Normalize the data
                    if sat_flux and normalize: 
                        vec_flux = vec_flux/sat_flux[experiment.id][simulation.id][specie]
                        try: 
                            vec_time = vec_time/t_range[experiment.id][simulation.id]["time peak"]
                        except:
                            t_range[experiment.id][simulation.id]["time peak"] = calculate_peakTime(vec_time, vec_flux) 
                            vec_time = vec_time/t_range[experiment.id][simulation.id]["time peak"]
                    
                    # Put time time axis in milliseconds if we plot in SI units
                    if units!="normalized": vec_time=vec_time*1000  
                     
                    # Add a label to the legend for the first simulation of this experiment
                    label_exp = experiment.line_label if simulations.index(simulation) == 0 else ""

                    # Plot the data
                    ax.plot(vec_time, vec_flux, lw=3, linestyle=style, color=experiment.line_color, label=label_exp) 
 
    # Plot the data points with the same color for each varied_value
    for experiment in experiments:
        
            # Iterate over the simulations
            simulations = get_simulation(experiment, simulation_id)
            for simulation in simulations:
                
                # Iterate over the species
                for specie in species:
                     
                    # Get the data
                    vec_time = simulation.vec_time
                    vec_flux = getattr(simulation, y_quantity)[:,specie]
                    style    = ':' if (specie=='1' and len(species)>1) else '-'
                    
                    # Normalize the data
                    if sat_flux and normalize: 
                        vec_flux = vec_flux/sat_flux[experiment.id][simulation.id][specie]
                        try: 
                            vec_time = vec_time/t_range[experiment.id][simulation.id]["time peak"]
                        except:
                            t_range[experiment.id][simulation.id]["time peak"] = calculate_peakTime(vec_time, vec_flux) 
                            vec_time = vec_time/t_range[experiment.id][simulation.id]["time peak"]
                        
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
                    ax.plot(vec_time, vec_flux, lw=lw, linestyle=style, color=marker_color, label=label_var1) 
     
                    # Keep track of the axis limits
                    xlims = [np.nanmin([np.nanmin(vec_time), xlims[0]]), np.nanmax([1.1*np.nanmax(vec_time), xlims[1]])]
                    ylims = [np.nanmin([np.nanmin(vec_flux), ylims[0]]), np.nanmax([1.1*np.nanmax(vec_flux), ylims[1]])]
                    
                    # Next print the symbols/markers with the rho label, print the markers differently for log and normal scales
                    if (len(experiments) > 1) or (number_ofVariedValues==1):
                        vec_flux, vec_time = get_markerAtEveryXSteps(vec_flux, vec_time, log, units)
                        ax.plot(vec_time, vec_flux, lw=0, marker=marker_style, mec=marker_color, mfc=marker_color, label=label_var2) 
 

    # Show the area where the saturated flux is calculated
    if show_restartTimes:    
        for experiment in experiments:
            simulations = get_simulation(experiment, simulation_id)
            for simulation in simulations:
                for time in simulation.restart_times:
                    print("RESTART TIME", experiment.id, simulation.id, time)
                    ax.axvline(x=time, color=experiment.line_color)                         
                         
    # Show the area where the saturated flux is calculated
    if sat_flux and show_vspan:    
        for experiment in experiments:
            simulations = get_simulation(experiment, simulation_id)
            for simulation in simulations:
                for specie in species:
                    t_start = t_range[experiment.id][simulation.id][specie][0]
                    t_stop  = t_range[experiment.id][simulation.id][specie][1]
                    if normalize:
                        t_start = t_start/t_range[experiment.id][simulation.id]["time peak"]
                        t_stop  = t_stop/t_range[experiment.id][simulation.id]["time peak"]
                    ax.axvspan(t_start, t_stop, alpha=0.5, color=simulation.marker_color)
 
    # Add a line to indicate the saturated flux mean
    if sat_flux:
        for experiment in experiments:
            simulations = get_simulation(experiment, simulation_id)
            for simulation in simulations:
                for specie in species:
                     
                    # The label can show the exact value
                    label = "" #'$Q_{sat}/Q_{gB} =$'+"{:.2f}".format(sat_flux[simulation][specie])
                     
                    # The color either matches the experiment, the varied value or the species
                    if (len(experiments) > 1) or (number_ofVariedValues==1): color = experiment.line_color
                    else : color = simulation.line_color
                    style    = ':' if (specie=='1' and len(species)>1) else '-'
                     
                    # Only plot the saturated flux over the time range where it is averaged over
                    t_start = t_range[experiment.id][simulation.id][specie][0]
                    t_stop  = t_range[experiment.id][simulation.id][specie][1]
                    vec_time = np.linspace(t_start, t_stop, 10) 
                     
                    # Manually construct the saturated flux vector
                    saturated_flux = np.ones(len(vec_time))*sat_flux[experiment.id][simulation.id][specie]
                     
                    # Plot the saturated flux
                    ax.plot(vec_time, saturated_flux, lw=2, linestyle=style, color=color, label=label) 
                    
                    
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

def calculate_peakTime(vec_time, vec_flux):
    ''' Calculate the time of the peak of the ramp up. '''
    time = vec_time[vec_flux>0.001] 
    fluxes = vec_flux[vec_flux>0.001]
    fluxes = [ fluxes[i+1]-fluxes[i] for i in range(len(fluxes)-1)]
    fluxes = [ i for i in range(len(fluxes)) if fluxes[i]<0]
    time_peak = time[fluxes[0]]
    return time_peak

def determine_labels(x_label, y_label, title, units, species, normalize):
    ''' Set the labels and titels for the (omega, ki) or (gamma, ki) plot '''
    
    # Save the labels in one dictionary
    label = {'x' : x_label, 'q_flux' : y_label, 'v_flux' : y_label, 'p_flux' : y_label, 'title' : title}

    # Determine the label of the x-axis
    if label["x"] is None and not normalize:
        label["x"] = "$t\, v_{\mathrm{th}}/a$" if units=="normalized" else "$t$ [ms]"
    if label["x"] is None and normalize:
        label["x"] = "$(t\, v_{\mathrm{th}}/a) / t_{peak}$" if units=="normalized" else "$t/t_{peak}$"

    # Get the species subscript
    if len(species)==1:
        if int(species[0]) == 0: s = "_i" 
        if int(species[0]) == 1: s = "_e" 
    else:  s = "_s" 

    # Get the correct label for the y-axis which depends on the specie
    if not normalize:
        if label["q_flux"]==None: label["q_flux"] = "$Q"+s+"/Q_{gB}$" if units=="normalized" else "$Q"+s+"$ [W/m$^2$]"
        if label["p_flux"]==None: label["p_flux"] = "$\\Gamma"+s+"/\\Gamma_{gB}$" if units=="normalized" else "$\\Gamma"+s+"$ [W/m$^2$]"
        if label["v_flux"]==None: label["v_flux"] = "$\\Pi"+s+"/\\Pi_{gB}$" if units=="normalized" else "$\\Pi"+s+"$ [W/m$^2$]"
    if normalize:
        if label["q_flux"]==None: label["q_flux"] = "$Q"+s+"/(Q_{gB} \\cdot Q_{sat})$" if units=="normalized" else "$Q"+s+"/Q_{sat}$"
        if label["p_flux"]==None: label["p_flux"] = "$\\Gamma"+s+"/(\\Gamma_{gB} \\cdot \\Gamma_{sat})$" if units=="normalized" else "$\\Gamma"+s+"/\\Gamma_{sat}$"
        if label["v_flux"]==None: label["v_flux"] = "$\\Pi"+s+"/(\\Pi_{gB} \\cdot \\Pi_{sat})$" if units=="normalized" else "$\\Pi"+s+"/\\Pi_{sat}$"
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


