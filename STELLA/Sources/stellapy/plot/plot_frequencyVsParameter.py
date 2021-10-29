
#================================================================
# Plot omega(k) or gamma(k)
#================================================================

# Load the modules
import numpy as np
import matplotlib.pyplot as plt

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
def plot_frequencyVsParameter(\
            # Specify which simulations to plot
                research=None,\
                experiment_id="All experiments",\
                simulation_id="All simulations",\
                parameter_knob="vmec_parameters",\
                parameter_key="rho",\
            # Specify data range
                y_quantity='omega',\
                x_range=None,\
                y_range=None,\
            # Details of the modes
                kx_range=-0.0, \
                ky_range=[0,100], \
                k_value=9.0, \
                lineardata="average",\
            # Labels
                x_label=None,\
                y_label=None,\
                title=None,\
            # For the GUI the figure object already exists 
                show_error = False,\
                show_figure = False,\
                ax=None, \
                Progress=None,\
                root=None):
    ''' Plot the frequency/omega or growthrate/gamma versus the wave number. 
    
    Parameters
    ----------
    research: object
        Contains the selected experiments and simulations.
    experiment_id, simulation_id: str
        Is used to select which simulations/experiments to plot.
    kvalue: {"max", integer}
        Determines which modes to plot.
    units: {"normalized", "SI units"}
        Determines the units of the data.
    lineardata: {"average", "last"}
        Either plot the last time value of omega/gamma or the average
        over the last 10% of time.
    '''


    # Update the progress bar of the GUI
    if Progress: Progress.move(0,"Plot omega or gamma versus parameter.")
         
    # Decide whether to scan modes along the x-axis or y-axis
    scan, k_fixed, k_range = get_axisOfScan(kx_range, ky_range, root) 
    if scan == False: return 
        
    # Set the labels, and change the standard labels depending on the units
    label = determine_labels(x_label, y_label, kx_range, ky_range, title,  "normalized", k_value, parameter_key)
    load_plotbox2d(x_label=label['x'], y_label=label[y_quantity], title=label["title"], ax=ax) 

    # Keep track of the labels that are already used and save the axis limits
    labels = []; xlims=[0,0]; ylims=[0,0]; 
    
    # Get the same amount of colors as experiments
    experiments = get_experiment(research, experiment_id)
    color = plt.cm.jet(np.linspace(0,1,len(experiments))) #@undefinedvariable
    
    # Iterate over the experiments
    for experiment in experiments:
        
        #==================
        # GET THE DATA
        #==================
        
        # Get the simulations, the parameters and the omega/gamma at k_value
        simulations    = get_simulation(experiment, simulation_id)
        y_dataAtKvalue = [np.NaN]*len(simulations)
        y_error        = np.empty((2, len(simulations))); y_error[:, :] = np.NaN
        parameters     = [] 
        
        # Get the requested parameter
        for simulation in simulations:
            parameters.append(float(simulation.inputParameters[parameter_knob][parameter_key]))
        
        # Sort the simulations by this parameter
        sorted_indexes = list(np.array(parameters).argsort(kind="heapsort"))
        simulations = [simulations[i] for i in sorted_indexes]
        parameters  = [parameters[i] for i in sorted_indexes]
        
        # Iterate over the simulations
        for simulation in simulations:
            
            # Get the index of the simulation
            i = simulations.index(simulation)
            
            # Update the progress bar of the GUI 
            if Progress: Progress.move(i/len(simulations)*100,"Doing analysis ("+str(i)+"/"+str(len(simulations))+")")
                
            # Get the modes of this simulation that need to be plotted and their indices in the (kx,ky) matrixes
            modes  = simulation.get_modesForAOneDimensionalScan(scan, k_fixed, k_range, plotted_modes="unstable") 
            i_kxky = simulation.get_indicesOfModes(scan, k_fixed, modes)    
            
            # Get the omega or gamma at the last time value
            if y_quantity=="omega" and lineardata=="average": y_data = simulation.omega_avg[i_kxky]
            if y_quantity=="omega" and lineardata=="last":    y_data = simulation.omega_last[i_kxky]
            if y_quantity=="gamma" and lineardata=="average": y_data = simulation.gamma_avg[i_kxky]
            if y_quantity=="gamma" and lineardata=="last":    y_data = simulation.gamma_last[i_kxky]
            if y_quantity=="gamma/ky**2" and lineardata=="average": y_data = simulation.gamma_avg[i_kxky]/[m**2 for m in modes]
            if y_quantity=="gamma/ky**2" and lineardata=="last":    y_data = simulation.gamma_last[i_kxky]/[m**2 for m in modes]
            
            # Get the error bars
            if y_quantity=="omega": y_error_min = simulation.omega_min[i_kxky]
            if y_quantity=="omega": y_error_max = simulation.omega_max[i_kxky]
            if y_quantity=="gamma": y_error_min = simulation.gamma_min[i_kxky]
            if y_quantity=="gamma": y_error_max = simulation.gamma_max[i_kxky]
            if y_quantity=="gamma/ky**2": y_error_min = simulation.gamma_min[i_kxky]/[m**2 for m in modes]
            if y_quantity=="gamma/ky**2": y_error_max = simulation.gamma_max[i_kxky]/[m**2 for m in modes]
    
            # Get the maximum of omega/gamma or omega/gamma at a specific ky
            if k_value == "max":
                try:
                    gamma_data = simulation.gamma_avg[i_kxky]
                    index = np.nanargmax(gamma_data)
                    y_dataAtKvalue[i] = y_data[index] 
                    y_error[0,i] = y_error_min[index]
                    y_error[1,i] = y_error_max[index]
                except: pass
                
            if k_value != "max":
                mask = np.isin(modes, k_value)
                if k_value in modes:  
                    y_dataAtKvalue[i] = y_data[mask]
                    y_error[0,i] = y_error_min[mask]
                    y_error[1,i] = y_error_max[mask]
                else: pass


        #==================
        # PLOT THE DATA
        #==================

        # Update the progress bar of the GUI 
        if Progress: Progress.move(i/len(simulations)*100,"Plotting analysis ("+str(i)+"/"+str(len(simulations))+")")
            
        # Plot the data points with the same color for each varied_value 
        if len(modes)>0:
            # Add a label to the legend for the first simulation of this experiment
            label_var = experiment.line_label 
            label_var = label_var.replace("_"," ") if "$" not in label_var else label_var
            label_var = "" if (label_var in labels) or (len(experiments)<=1) else label_var
            labels.append(label_var)
            marker_color = color[experiments.index(experiment)]
        
            # Plot the actual linear data points obtained by stella as markers (mfc = markerfacecolor; mec = markeredgecolor]   
            print()
            print("Parameter", experiment.id)
            print(parameters)
            print(y_dataAtKvalue)
            ax.plot(parameters, y_dataAtKvalue, lw=2 ,linestyle='-', color=marker_color, ms=5 ,\
                        marker="o", mec="black", mfc="white", label=label_var)
            
            # Add the error bars
            if show_error==True:
                ax.errorbar(parameters, y_dataAtKvalue, yerr=y_error, fmt='o', color=marker_color, capsize=2)
            
            # Keep track of the axis limits
            xlims = [np.nanmin([np.nanmin(parameters), xlims[0]]), np.nanmax([np.nanmax(parameters), xlims[1]])]
            ylims = [np.nanmin([np.nanmin(y_dataAtKvalue), ylims[0]]), np.nanmax([np.nanmax(y_dataAtKvalue), ylims[1]])]
            
    # Change appearance plot 
    ax.autoscale()
    ax.set_title(title)
    ax.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
    ax.legend(loc='best',labelspacing=0.0, prop={'size':20}, framealpha=0.9)
    
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

def determine_labels(x_label, y_label, kx_range, ky_range, title,  units, k_value, parameter_key):
    ''' Set the labels and titels for the (omega, ki) or (gamma, ki) plot '''
    
    # Save the labels in one dictionary
    label = {'x' : x_label, 'omega' : y_label, 'gamma' : y_label, 'gamma/ky**2' : y_label,  'title' : title}
    
    # Find out the label of the x-axis
    if   parameter_key=="rho":      label["x"] = "$\\rho$"
    elif parameter_key=="tprim":    label["x"] = "$a/L_{T_i}$"
    elif parameter_key=="tiprim":   label["x"] = "$a/L_{T_i}$"
    elif parameter_key=="teprim":   label["x"] = "$a/L_{T_e}$"
    elif parameter_key=="fprim":    label["x"] = "$a/L_{n}$"
    elif parameter_key=="delta t":  label["x"] = "$\\Delta t$"
    elif parameter_key=="delt":     label["x"] = "$\\Delta t$"
    elif parameter_key=="nmu":      label["x"] = "$n_{\\mu}$"
    elif parameter_key=="nvgrid":   label["x"] = "$n_{v}$"
    else:                           label["x"] = parameter_key

    # Set the label for the y-axis
    if label["omega"]==None or True:  
        if k_value=="max":  
            if units=="normalized": label["omega"] = '$\\omega_{max} a/v_{\\mathrm{th},i}$'
            if units=="SI units":   label["omega"] = '$\\omega_{max}$ [s$^{-1}$]' 
        if k_value!="max":  
            if units=="normalized": label["omega"] = '$\\omega a/v_{\\mathrm{th},i}$ at $k_y\,=\,$'+str(k_value)
            if units=="SI units":   label["omega"] = '$\\omega$ [s$^{-1}$] at $k_y\,=\,$'+str(k_value)

    if label["gamma"]==None or True:     
        if k_value=="max":  
            if units=="normalized": label["gamma"] = '$\\gamma_{max} a/v_{\\mathrm{th},i}$'
            if units=="SI units":   label["gamma"] = '$\\gamma_{max}$ [s$^{-1}$]' 
        if k_value!="max":  
            if units=="normalized": label["gamma"] = '$\\gamma a/v_{\\mathrm{th},i}$ at $k_y\,=\,$'+str(k_value)
            if units=="SI units":   label["gamma"] = '$\\gamma$ [s$^{-1}$] at $k_y\,=\,$'+str(k_value)
            
    if label["gamma/ky**2"]==None or True:     
        if k_value=="max":  
            if units=="normalized": label["gamma/ky**2"] = '$\\gamma_{max} a/v_{\\mathrm{th},i} * (1/k_y^2)$'
            if units=="SI units":   label["gamma/ky**2"] = '$\\gamma_{max}/k_y^2$ [s$^{-1}$]' 
        if k_value!="max":  
            if units=="normalized": label["gamma/ky**2"] = '$\\gamma a/v_{\\mathrm{th},i} * (1/k_y^2)$ at $k_y\,=\,$'+str(k_value)
            if units=="SI units":   label["gamma/ky**2"] = '$\\gamma/k_y^2$ [s$^{-1}$] at $k_y\,=\,$'+str(k_value)
  
    return label

