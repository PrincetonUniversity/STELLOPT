
import numpy as np
import os, configparser, pathlib
from itertools import takewhile
from stellapy.simulations.Simulation import create_simulations
from stellapy.simulations.utils.get_differenceInDictionaries import get_differenceInDictionaries
from stellapy.simulations.utils.load_labelsLinesMarkers import load_labelsLinesMarkers
from stellapy.simulations.utils.calculate_attributeWhenReadFirstTime import calculate_attributeWhenReadFirstTime
from stellapy.utils.decorators import printv 
from stellapy.data.read_netcdf import read_netcdf


#################################################################
#                        CREATE EXPERIMENTS
#################################################################

def create_experiments(\
    # To create the simulations we need their location and whether to ignore the resolution
    folders=None, input_files=None, ignore_resolution=True, \
    # To group them by experiment we need to know how many variables differ between the simulations
    # Or which variable is unique for each experiment
    number_variedVariables=1, experiment_knob="vmec_parameters", experiment_key="vmec_filename",\
    # Give some extra rules
    folderIsExperiment = False):
    ''' Divide the simulations in experiments, based on the number of varied values
    or based on a parameter that differs for each experiment. For example, when we 
    scan multiple radial positions of multiple shots, then each shot is defined by
    its unique magnetic field file. If they are run for the same magnetic field, just
    make sure they have a different symbolic link to differ between the experiments. '''
    
    # First create the simulations
    simulations = create_simulations(folders, input_files, ignore_resolution, number_variedVariables, folderIsExperiment)
    
    # Create an empty list to hold the experiment objects
    experiments = []

    # Iterate over the simulations objects
    for simulation in simulations:
        
        # Assume we have a new experiment
        newExperiment = True
        
        # Go through all the current experiments and see if the simulation belongs to one of these
        for experiment in experiments:
            
            # If the simulations are in a different folder its a new experiment. 
            if not (folderIsExperiment==True and simulation.parent != experiment.simulations[0].parent):
                
                # Look at the difference in the input parameters of the experiment and the simulation
                dict_difference = get_differenceInDictionaries(experiment.inputParameters, simulation.inputParameters)
                
                # Sometimes we rerun the exact same simulation with different modes: we have the same experiment
                # However if we want this simulation to be a new simulation, force this with number_variedValues=-1
                # If number_variedVariables==-2 then we want to sort by "1 folder = 1 experiment"
                # If number_variedVariables==-1 then we want to sort by "1 folder = 1 simulation"
                if len(list(dict_difference.keys())) == 0:
                    if number_variedVariables==-2:
                        newExperiment = True
                    if number_variedVariables!=-2:
                        newExperiment = False
                        experiment.simulations.append(simulation) 

                # The simulation belongs to the current experiment if
                # inputParameters[experiment_knob][experiment_key] is the same.
                # For example if the radial position is varied then different experiments are characterized
                # by having a different vmec_filename in the knob "vmec_parameters". 
                # Therefore if inputParameters[experiment_knob][experiment_key]
                # Is the same in the simulation and experiment, add the simulation to this experiment.
                # This rule is overriden if less than or equal to <number_variedVariables> are varied.
                elif experiment_knob in list(dict_difference.keys()):
                    if experiment_key not in dict_difference[experiment_knob]:
                        newExperiment = False
                        experiment.simulations.append(simulation)
                        varied_variable = experiment_knob + ": " + dict_difference[experiment_knob][0]
                        if varied_variable not in experiment.variedVariables:
                            experiment.variedVariables.append(varied_variable)
                            
                # The simulation is part of the experiment if up to <number_variedVariables> input variables are different
                elif len(list(dict_difference.keys()))<=number_variedVariables:
                    varied_values = 0
                    for knob in list(dict_difference.keys()):
                        for parameter in dict_difference[knob]:
                            varied_values += 1
                    if varied_values <= number_variedVariables:
                        newExperiment = False
                        experiment.simulations.append(simulation)
                        for knob in list(dict_difference.keys()):
                            for parameter in dict_difference[knob]:
                                varied_variable = knob + ": " + parameter
                                if varied_variable not in experiment.variedVariables:
                                    experiment.variedVariables.append(varied_variable)
        
        # If more than <number_variedVariables> variables are different we have a new experiment
        # If inputParameters[experiment_knob][experiment_key] is different we have a new experiment
        # If the experiments list is empty we have a new experiment
        # In this case create an <Experiment> object and add it to the list
        if newExperiment==True:
            experiments.append(Experiment(simulation))

    # For the variables that are varied, get their exact values for the labels in the graph
    # For example if experiment.variedVariables = ["vmec_parameter: torflux", "knobs: delt"]
    # then experiment.variedValues = ["rho=0.5; delt=0.1", "rho=0.5; delt=0.01"; "rho=0.6; delt=0.1"]
    for experiment in experiments:
        for simulation in experiment.simulations:
            varied_values = "" # This will be used as the label in the figure
            for variable in experiment.variedVariables:
                knob = variable.split(":")[0]
                variable = variable.split(": ")[-1]
                if knob != experiment_knob and variable != experiment_key:
                    value = simulation.inputParameters[knob][variable]
                    variable = variable.replace("_", " ") 
                    if variable=="torflux": value = str(round(np.sqrt(float(value)),2))
                    if variable=="torflux": variable = "$\\rho$"
                    if variable=="rho":     variable = "$\\rho$"
                    if variable=="delt":    variable = "$\\Delta t$"
                    if knob=='species_parameters_1':
                        if variable=="tprim": variable="$a/L_{Ti}$"
                        if variable=="fprim": variable="$a/L_{ne}$"
                    if knob=='species_parameters_2':
                        if variable=="tprim": variable="$a/L_{Te}$"
                        if variable=="fprim": variable="$a/L_{ni}$"
                    varied_values = varied_values + variable + "$\,=\,$" + str(value) + "; "
             
            # If a/Lne = a/Lni replace it by a/Ln
            if "$a/L_{ni}$" in varied_values and "$a/L_{ne}$" in varied_values:
                _elements = varied_values.split("; ")
                _elements = [e for e in _elements if e!='']
                _ionDensityGrad = float([e for e in _elements if "$a/L_{ni}$" in e][0].split("$\,=\,$")[-1])
                _others = [e for e in _elements if ("$a/L_{ni}$" not in e) and ("$a/L_{ne}$" not in e)] 
                if len(_others) != 0: varied_values = "; ".join(_others) + "; " + "$a/L_{n}$" + "$\,=\,$" + str(_ionDensityGrad)+"; "
                if len(_others) == 0: varied_values = "$a/L_{n}$" + "$\,=\,$" + str(_ionDensityGrad)+"; "
                          
            # If nothing was different because we only have one experiment, use the radial position as the label
            if varied_values=="":
                varied_values = "$\\rho$" + "$\,=\,$" +  str(experiment.inputParameters['vmec_parameters']['rho'])+"; "
            
            # The time step in the input_file isn't the actual time step
            if "$\\Delta t$" in varied_values and False:
                netcdf_data = read_netcdf(simulation.input_files[0])   
                try:    delt = (netcdf_data['vec_time'][5]-netcdf_data['vec_time'][4])/simulation.inputParameters['stella_diagnostics_knobs']['nwrite']
                except: delt = (netcdf_data['vec_time'][5]-netcdf_data['vec_time'][4])/simulation.inputs[simulation.input_files[0]]['stella_diagnostics_knobs']['nwrite']
                simulation.inputParameters['knobs']['delt'] = round(delt,4)
                varied_values = varied_values + "$\\delta t$$\,=\,$" + str(round(delt,4)) + "; "
            
            # Add the string of varied values to the dictionary
            experiment.variedValues.append(varied_values[0:-2])
            
    # Make sure we don't have the same label for each simulation
    for experiment in experiments:
        if len(set(experiment.variedValues))==1 and len(experiment.variedValues)!=1:
            simulation_ids = [simulation.id for simulation in experiment.simulations]
            common_prefix = ''.join(c[0] for c in takewhile(lambda x:  all(x[0] == y for y in x), zip(*simulation_ids)))
            common_prefix = common_prefix if common_prefix.endswith("_") else '_'.join(common_prefix.split('_')[0:-1])+"_"
            simulation_ids2 = [simulation.id.split("__")[0] for simulation in experiment.simulations]
            common_prefix2 = ''.join(c[0] for c in takewhile(lambda x:  all(x[0] == y for y in x), zip(*simulation_ids2)))
            common_prefix2 = common_prefix if common_prefix2.endswith("_") else '_'.join(common_prefix2.split('_')[0:-1])+"_"
            for i in range(len(experiment.variedValues)):
                if "__" in common_prefix:
                    experiment.variedValues[i] = str(experiment.simulations[i].id).split("__")[-1]
                if "__" not in common_prefix:
                    experiment.variedValues[i] = str(experiment.simulations[i].id).split("__")[0].split(common_prefix2)[1]
            # Add the delt t for now
            if False:
                for i in range(len(experiment.variedValues)):
                    simulation = experiment.simulations[i]
                    for input_file in simulation.input_files:
                        netcdf_data = read_netcdf(input_file)   
                        delt1 = (netcdf_data['vec_time'][5]-netcdf_data['vec_time'][4])/simulation.inputParameters['stella_diagnostics_knobs']['nwrite']
                        delt2 = (netcdf_data['vec_time'][-4]-netcdf_data['vec_time'][-5])/simulation.inputParameters['stella_diagnostics_knobs']['nwrite']
                        print("   Time step input:", simulation.inputParameters['knobs']['delt'])
                        print("   Time step start:", delt1)
                        print("   Time step end:", delt2)
                    experiment.variedValues[i] =  experiment.variedValues[i] + "; $\\delta t$$\,=\,$" + str(round(delt1,4))
                    #experiment.variedValues[i] =  experiment.variedValues[i] + "; $\\delta_e t$$\,=\,$" + str(round(delt2,4))
          
#     # Sort the simulations and varied values
#     numbers = [ int("".join([s for s in value if s.isdigit()])) for value in experiment.variedValues ]
#     sorted_indexes = list(np.array(numbers).argsort())
#     experiment.variedValues = [experiment.variedValues[i] for i in sorted_indexes]
#     experiment.simulations  = [experiment.simulations[i] for i in sorted_indexes]   
            
    # Now that all the simulations are added to experiments, we an add labels and markers and write config
    for experiment in experiments:
        experiment.finish_initialization()
        
    # Print some information
    for experiment in experiments:
        printv("")
        printv("Experiment "+str(experiments.index(experiment)+1)+":")
        printv("     Experiment ID:  "+experiment.id)
        printv("     Varied Values:  "+'["'+('", "').join(experiment.variedValues)+'"]')
        printv("     Simulation IDs: "+"["+(",").join([ s.id for s in experiment.simulations])+"]")
        printv("")
    
    return experiments


#################################################################
#                      CLASS EXPERIMENTS
#################################################################

class Experiment:
    ''' Make an experiment object which groups multiple simulation objects. '''
    
    def __init__(self, simulation):
        
        # Save the simulation objects and its input parameters
        # Note that the input parameters of the simulation can differ in (kx,ky)
        # And also in the resolution if ignore_resolution = True
        self.id = simulation.id
        self.simulations = [simulation]     
        self.inputParameters = simulation.inputParameters
        
        # Save whichs [knob][parameter] is different in this experiment compared to the
        # other experiments: e.g. variedVariables = ["vmec_parameter: torflux", "knobs: delt"]
        # For each simulation add the exact value of this parameter that is different:
        # e.g. variedValues = ["rho=0.5; delt=0.1", "rho=0.5; delt=0.01"; "rho=0.6; delt=0.1"]
        self.variedVariables = []
        self.variedValues = []
        
        # Print information on the GUI (set from the GUI itself)
        self.Progress = None
      
      
    #------------------------------  
    def finish_initialization(self):
        
        # Now write the configuration file of the simulation if it doesn't exist
        self.read_configurationFile()
        
        # Get some default values for the markers and lines
        self.set_labelsLinesMarkers()
        
        # Re-add the simulations to the configuration file since this might change while the id remains the same
        simulation_ids = [s.id for s in self.simulations]
        simulations_config = self.configuration_file['GENERAL']['simulations'].split("\n")
        simulations_config = [s.replace("\t", "") for s in simulations_config if (s!='' and s!="\t")]
        for simulation_config in simulations_config:
            if simulation_config not in simulation_ids:
                try: del self.configuration_file[simulation_config]
                except: pass
                try: del self.configuration_file['OVERRIDE COLORS AND FONTS FOR SIMULATION '+simulation_config]
                except: pass
        for simulation in simulation_ids:
            if simulation not in simulations_config:
                input_files = [str(i) for i in self.simulations[simulation_ids.index(simulation)].input_files]
                input_files = "\t\n" + ("\t\n").join(input_files)
                self.configuration_file[simulation] = {
                    'input files'  : input_files}
                self.configuration_file['OVERRIDE COLORS AND FONTS FOR SIMULATION '+simulation] = {
                    'line_label'   : 'USE DEFAULT',\
                    'line_style'   : 'USE DEFAULT',\
                    'line_color'   : 'USE DEFAULT',\
                    'marker_label' : 'USE DEFAULT',\
                    'marker_style' : 'USE DEFAULT',\
                    'marker_color' : 'USE DEFAULT'} 
        self.configuration_file.write(open(self.path_configurationFile, 'w'))
        
        # Dont collapse header on the next line
        if True: return

#====================================================
#  GET THE TOTAL AMOUNT OF MODES IN THE EXPERIMENT
#====================================================

    @calculate_attributeWhenReadFirstTime 
    def total_kx(self):        self.get_modes();    return self.total_kx
    @calculate_attributeWhenReadFirstTime 
    def total_ky(self):        self.get_modes();    return self.total_ky
    def get_modes(self):
        ''' Sort the modes by stable/unstable and converged/unconverged. '''
        
        # Get the total amount of modes
        self.total_kx, self.total_ky = [], []
        for simulation in self.simulations:          
            self.total_kx += list(simulation.vec_kx)
            self.total_ky += list(simulation.vec_ky)
        self.total_kx = list(set(self.total_kx)); self.total_kx.sort()
        self.total_ky = list(set(self.total_ky)); self.total_ky.sort()
        if True: return
    

#=======================================
#       STYLE OF THE EXPERIMENT
#=======================================

    def set_labelsLinesMarkers(self, line_label=None, marker_label=None, line_style=None, line_color=None, marker_style=None, marker_color=None):
        
        # If the styles are given apply them
        if line_label is not None:
            self.line_label = line_label[self.id]
            self.marker_label = marker_label[self.id]
            self.line_style = line_style[self.id]
            self.line_color = line_color[self.id]
            self.marker_style = marker_style[self.id]
            self.marker_color = marker_color[self.id]
        
        # Otherise load some basic styles
        else: 
            line_label, marker_label, line_style, line_color, marker_style, marker_color = load_labelsLinesMarkers([self.id], self.variedValues)
            self.line_label = line_label[self.id]
            self.marker_label = marker_label[self.id]
            self.line_style = line_style[self.id]
            self.line_color = line_color[self.id]
            self.marker_style = marker_style[self.id]
            self.marker_color = marker_color[self.id]
        
        # If the simulation was ran before, there is a text file, then read the manually set styles
        style = self.configuration_file['OVERRIDE COLORS AND FONTS']
        if style['line_label']   != 'USE DEFAULT':  self.line_label   = style['line_label']
        if style['line_color']   != 'USE DEFAULT':  self.line_color   = style['line_color']
        if style['line_style']   != 'USE DEFAULT':  self.line_style   = style['line_style']
        if style['marker_label'] != 'USE DEFAULT':  self.marker_label = style['marker_label']
        if style['marker_color'] != 'USE DEFAULT':  self.marker_color = style['marker_color']
        if style['marker_style'] != 'USE DEFAULT':  self.marker_style = style['marker_style']

        # Save the style of the simulation as an indicator for manual inputs
        self.configuration_file['COLORS AND FONTS'] = {
            'line_label'   : self.line_label,\
            'line_style'   : self.line_style,\
            'line_color'   : self.line_color,\
            'marker_label' : self.marker_label,\
            'marker_style' : self.marker_style,\
            'marker_color' : self.marker_color} 
        
        # Write the new data to the configuration file
        self.configuration_file.write(open(self.path_configurationFile, 'w'))
        
        # Change the lines/markers of the underlying simulation objects based on self.variedValues
        for simulation in self.simulations:
            variedValue = self.variedValues[self.simulations.index(simulation)]
            simulation.set_labelsLinesMarkers(line_label[variedValue], marker_label[variedValue], line_style[variedValue], \
                                              line_color[variedValue], marker_style[variedValue], marker_color[variedValue])
        

        # Dont collapse header on the next line
        if True: return
    
#============================================
#  SAVE THE SIMULATION TO THE CONFIG FOLDER
#============================================

    def read_configurationFile(self):
        ''' Check whether <config.ini> exists, if not write it, CONFIG will be read after calling this function.
        Reading the file will check it existence and whether the paths are correct. '''
        
        # Create a configuration object
        self.configuration_file = configparser.ConfigParser() 
        
        # Get the location of the configuration file
        path_stellaGUI = pathlib.Path(os.path.dirname(os.path.abspath(__file__)).split("stellapy")[0])
        self.path_configurationFile =  path_stellaGUI / "stellapy/config/experiments" / str("experiment_"+self.id+".ini")
                
        # Read the configuration file of the simulation
        self.configuration_file.read(self.path_configurationFile)
        
        # The configuration file didn't exist when it was read
        if "GENERAL" not in self.configuration_file:
            self.create_defaultConfigurationFile()
        return
     
    #--------------------------------------
    def create_defaultConfigurationFile(self):
        ''' Write a default configuration file and save it as "stellapy/config/experiments/id.ini". 
        The goal is to use this file to reload the data in the GUI. '''
        
        # Get the simulations of the experiment
        simulations = [s.id for s in self.simulations]
        simulations = "\t\n" + ("\t\n").join(simulations)
            
        # Remember which experiment this is   
        self.configuration_file['GENERAL'] = {
            'experiment id' : self.id, \
            'simulations': simulations}
        
        # Style the experiment
        self.configuration_file['OVERRIDE COLORS AND FONTS'] = {
            'line_label'   : 'USE DEFAULT',\
            'line_style'   : 'USE DEFAULT',\
            'line_color'   : 'USE DEFAULT',\
            'marker_label' : 'USE DEFAULT',\
            'marker_style' : 'USE DEFAULT',\
            'marker_color' : 'USE DEFAULT'} 
        
        # Create a section for each simulation
        for simulation in self.simulations:
            input_files = [str(i) for i in simulation.input_files]
            input_files = "\t\n" + ("\t\n").join(input_files)
            self.configuration_file[simulation.id] = {
                'input files' : input_files}
            
            # Style the simulation
            self.configuration_file['OVERRIDE COLORS AND FONTS FOR SIMULATION '+simulation.id] = {
                'line_label'   : 'USE DEFAULT',\
                'line_style'   : 'USE DEFAULT',\
                'line_color'   : 'USE DEFAULT',\
                'marker_label' : 'USE DEFAULT',\
                'marker_style' : 'USE DEFAULT',\
                'marker_color' : 'USE DEFAULT'} 
            
        # Write the configuration file
        self.configuration_file.write(open(self.path_configurationFile, 'w'))
        
        # Dont collapse header on the next line
        if True: return
    



