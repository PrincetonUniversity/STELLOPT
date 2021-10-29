 
#================================================================
# Plot gamma(rho) for non-linear runs
#================================================================

# Load the modules
import numpy as np
import matplotlib.pyplot as plt

# Personal modules
from stellapy.plot.utils import load_plotbox2d  
from stellapy.utils.decorators import verbose_wrapper
from stellapy.simulations.utils.get_simulation import get_simulation
from stellapy.simulations.utils.get_experiment import get_experiment 
from stellapy.utils import initiate_nesteddict
 
@verbose_wrapper 
def plot_saturatedfluxVsRho(\
            # Specify which simulations to plot
                research=None,\
                experiment_id="All experiments",\
                simulation_id="All simulations",\
                parameter_knob="vmec_parameters",\
                parameter_key="rho",\
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
                log=False,\
                units="normalized"):
    ''' Plot saturated_flux(rho) for each shot, which is mean(heat_flux) from t=t_start to t=t_last. 
    
    Parameters
    ----------
    t_range : dict[experiment.id][simulation.id][specie] = tuple of floats
    
    
    Returns
    -------
    sat_flux : dict[experiment.id][simulation.id][specie] = float
    '''
     
    # Update the progress bar of the GUI
    if Progress: Progress.move(0,"Plot saturated flux versus parameter.")

    # Set the labels, and change the standard labels depending on the units
    label = determine_labels(x_label, y_label, title, units, species, parameter_key)
    load_plotbox2d(x_label=label['x'], y_label=label[y_quantity], title=label["title"], ax=ax) 

    # Keep track of the labels that are already used and save the axis limits
    labels = []; xlims=[0,0]; ylims=[0,0]; 
    
    #======================
    # GET THE TIME RANGES
    #======================
    
    # Make sure t_range contains a time range for each simulation and experiment
    if t_range is None:
        t_range = initiate_nesteddict()
    
    # Iterate over the experiments
    experiments = get_experiment(research, experiment_id)
    for experiment in experiments:
        
        # Iterate over the simulations
        simulations = get_simulation(experiment, simulation_id)
        for simulation in simulations:
            
            # Iterate over the species
            for s in range(simulation.dim_species):
                
                # Check whether a time range already exists
                timeRangeExists = False
                if experiment.id in t_range.keys():
                    if simulation.id in t_range[experiment.id].keys():
                        if s in t_range[experiment.id][simulation.id].keys():
                            timeRangeExists = True
            
                # If no t_range is assigned, choose the last 10% of the time vector 
                if not timeRangeExists:
                    t_start = simulation.vec_time[int(len(simulation.vec_time)*9/10)] 
                    t_range[experiment.id][simulation.id][s] = [t_start, np.NaN]
                    
                # For now always make t_end = t_last
                t_range[experiment.id][simulation.id][s][1] = np.NaN
                                
    # Delete the experiments and simulations that are no longer present
    experiment_ids = [ e.id for e in research.experiments]
    experiment_iids = list(t_range.keys())
    for experiment_iid in experiment_iids:
        if experiment_iid not in experiment_ids: 
            del t_range[experiment_iid]
        if experiment_iid in experiment_ids:
            simulation_ids = [ s.id for s in research.experiments[experiment_ids.index(experiment_iid)].simulations]
            simulation_iids = list(t_range[experiment_iid].keys())
            for simulation_iid in simulation_iids:
                if simulation_iid not in simulation_ids: 
                    try: del t_range[simulation_iid]
                    except: pass 

    #==================
    # GET THE DATA
    #==================

    # Initiate the saturated flux which is plotted in function of a parameter
    sat_flux = initiate_nesteddict() 
    parameters = initiate_nesteddict() 
    simulation_ids = initiate_nesteddict() 
    
    # Iterate over the experiments
    experiments = get_experiment(research, experiment_id)
    for experiment in experiments:
        parameters[experiment.id] = []
        simulation_ids[experiment.id] = []
        
        # Iterate over the simulations
        simulations = get_simulation(experiment, simulation_id)
        for simulation in simulations:
            
            # Get the requested parameter
            parameters[experiment.id].append(float(simulation.inputParameters[parameter_knob][parameter_key]))
            simulation_ids[experiment.id].append(simulation.id)
            
            # Iterate over the species
            for s in range(simulation.dim_species):
                
                # Get the data
                vec_time = simulation.vec_time
                vec_flux = getattr(simulation, y_quantity)[:, s]
                
                # Get the time frame
                t_start = t_range[experiment.id][simulation.id][s][0]
                t_stop  = np.nanmin([vec_time[-1], t_range[experiment.id][simulation.id][s][1]])
                t_range[experiment.id][simulation.id][s][1] = t_stop
                
                # Calculate the saturated flux within this time frame
                if t_stop <= t_start: 
                    sat_flux[experiment.id][simulation.id][s] = np.nan
                if t_stop > t_start: 
                    sat_flux[experiment.id][simulation.id][s] = vec_flux[(vec_time>t_start) & (vec_time<t_stop)].mean()
 
    #==================
    # PLOT THE DATA
    #==================
    
    # Iterate over the experiments
    experiments = get_experiment(research, experiment_id)
    for experiment in experiments:
        
        # Update the progress bar of the GUI 
        i = experiments.index(experiment); length=len(experiments)
        if Progress: Progress.move(i/length*100,"Plotting saturated fluxes ("+str(i)+"/"+str(length)+")")
        
        # Iterate over the species
        for specie in species:
        
            # Add a label to the legend for the first simulation of this experiment
            label_var = experiment.line_label 
            label_var = label_var.replace("_"," ") if "$" not in label_var else label_var
            label_var = "Electrons: "+label_var if specie==0 and len(species)>1 else label_var
            label_var = "Ions: "+label_var if specie==1 and len(species)>1 else label_var
            label_var = "Impurity: "+label_var if specie==2 and len(species)>1 else label_var
            label_var = "" if (label_var in labels) else label_var; labels.append(label_var)
            marker_color = experiment.line_color  
            
            # Get the saturated flux, sorted by the parameter
            sorted_indexes = list(np.array(parameters[experiment.id]).argsort(kind="heapsort"))
            parameters[experiment.id] = [parameters[experiment.id][i] for i in sorted_indexes]
            saturated_flux = [sat_flux[experiment.id][simulation_ids[experiment.id][i]][specie] for i in sorted_indexes]
        
            # Plot the saturated flux versus the parameter
            ax.plot(parameters[experiment.id], saturated_flux, lw=2 ,linestyle='-',\
                    color=marker_color, ms=5, marker="o", mec="black", mfc="white", label=label_var)
            
            # Keep track of the axis limits
            xlims = [np.nanmin([np.nanmin(parameters[experiment.id]), xlims[0]]), np.nanmax([1.1*np.nanmax(parameters[experiment.id]), xlims[1]])]
            ylims = [np.nanmin([np.nanmin(saturated_flux), ylims[0]]), np.nanmax([1.1*np.nanmax(saturated_flux), ylims[1]])]
          
            
    # Change appearance plot     
    ax.autoscale()
    ax.set_title(title)
    ax.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
    ax.legend(loc='best',labelspacing=0.0, prop={'size': 20}, framealpha=0.9)
    
    # Rescale the axis
    if x_range==None: x_range = xlims
    if y_range==None: y_range = ylims
    ax.set_xlim(x_range) 
    ax.set_ylim(y_range) 
    
    # Show the figure
    if show_figure: plt.show()
    return sat_flux, t_range
 
 
#################################################################
#                        METHODS
#################################################################

def determine_labels(x_label, y_label, title, units, species, parameter_key):
    ''' Set the labels and titels for the (omega, ki) or (gamma, ki) plot '''
    
    # Save the labels in one dictionary
    label = {'x' : x_label, 'q_flux' : y_label, 'v_flux' : y_label, 'p_flux' : y_label, 'title' : title}

    # Determine the label of the x-axis
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

    # Get the species subscript
    if len(species)==1:
        if int(species[0]) == 0: s = "_i" 
        if int(species[0]) == 1: s = "_e" 
    else:  s = "_s" 

    # Get the correct label for the y-axis which depends on the specie
    if label["q_flux"]==None: label["q_flux"] = "$Q"+s+"/Q_{gB}$" if units=="normalized" else "$Q"+s+"$ [W/m$^2$]"
    if label["p_flux"]==None: label["p_flux"] = "$\\Gamma"+s+"/\\Gamma_{gB}$" if units=="normalized" else "$\\Gamma"+s+"$ [W/m$^2$]"
    if label["v_flux"]==None: label["v_flux"] = "$\\Pi"+s+"/\\Pi_{gB}$" if units=="normalized" else "$\\Pi"+s+"$ [W/m$^2$]"
    return label
 
 


