
import numpy as np
import os, configparser, pathlib
from itertools import takewhile
from stellapy.utils.decorators import printv
from stellapy.simulations.utils.load_labelsLinesMarkers import load_labelsLinesMarkers
from stellapy.simulations.Experiment import create_experiments
from stellapy.config import turnOnVerboseWrapper_configurationFile

def create_research(\
    # To create the simulations we need their location and whether to ignore the resolution
    folders=None, input_files=None, ignore_resolution=True, \
    # To group them by experiment we need to know how many variables differ between the simulations
    # Or which variable is unique for each experiment
    number_variedVariables=1, experiment_knob="vmec_parameters", experiment_key="vmec_filename",\
    # Give the research a name in order to save it as a configuration file
    research_name = "DEFAULT", folderIsExperiment=False,\
    # Save some information from the plots
    data=None):
    ''' Group multiple experiments so they can be saved as a pickly object and handed to plotting functions. '''

    # First create the experiments
    experiments = create_experiments(folders, input_files, ignore_resolution, number_variedVariables, \
                                     experiment_knob, experiment_key, folderIsExperiment)
    
    # Sort the simulations in ascending order of the corresponding string variedValues
    for experiment in experiments: 
        try:
            numbers = [ float(value.split("$")[-1]) for value in experiment.variedValues ] 
            sorted_indexes = list(np.array(numbers).argsort()) 
            experiment.simulations = [experiment.simulations[i] for i in sorted_indexes] 
            experiment.variedValues = [experiment.variedValues[i] for i in sorted_indexes] 
        except: 
            sorted_indexes = list(np.array(experiment.variedValues).argsort())  
            experiment.simulations = [experiment.simulations[i] for i in sorted_indexes] 
            experiment.variedValues = [experiment.variedValues[i] for i in sorted_indexes] 
    
    # Look for the total amount of varied values throughout the experiments
    # e.g. experiment 1 has variedValues = ["rho=0.5", "rho=0.6"] and 
    # e.g. experiment 2 has variedValues = ["rho=0.5", "rho=0.7"] then
    # we want in total 3 colors for the lines/markers based on varied values 
    all_variedValues = []
    for experiment in experiments: 
            all_variedValues = all_variedValues + experiment.variedValues

    # Get the unique values
    all_variedValues = list(set(all_variedValues)) 
    all_variedValues.sort()
    
    # Create a research object
    research = Research(experiments, research_name, all_variedValues, number_variedVariables, \
                        ignore_resolution, data, experiment_knob, experiment_key)

    # Sort the experiments based on their label
    line_label = [ e.line_label for e in research.experiments]
    sorted_indexes = np.array(line_label).argsort()
    research.experiments = [research.experiments[i] for i in sorted_indexes]  
    return research

#################################################################
#                      CLASS RESEARCH
#################################################################
    
class Research:
    ''' Make a research object which groups multiple experiment objects. '''
    
    def __init__(self, experiments, research_name, variedValues, number_variedVariables,\
                 ignore_resolution, data, experiment_knob, experiment_key):
        
        # Save the experiments and the name of the research
        self.experiments = experiments
        self.id = research_name.replace(" ", "_")
        self.variedValues = variedValues
        self.number_variedVariables = number_variedVariables
        self.ignore_resolution = ignore_resolution
        self.experiment_knob = experiment_knob
        self.experiment_key = experiment_key
        
        # Save some data from the plots
        self.data = data
    
        # Now update the experiment labels
        self.update_ids()        
        
        # Get the input_files, folders and file names for the GUI
        self.get_foldersAndNamesForGui()
        
        # Get some default values for the markers and lines
        self.set_labelsLinesMarkers()
        
        # Write this data to a configuration file
        self.create_defaultConfigurationFile()
        if True: return
      
#################################################################
#                        METHODS
#################################################################

#========
#  GUI
#========
         
    def get_foldersAndNamesForGui(self):
        # Get the input files
        self.input_files = []
        for experiment in self.experiments:
            for simulations in experiment.simulations: 
                self.input_files += simulations.input_files 
        self.input_files = sorted(self.input_files)
        # Split the input files in folder/input_file pairs
        self.files   = [ i.name   for i in self.input_files ]
        self.folders = [ i.parent for i in self.input_files ] 
        self.e_id    = [ i.name   for i in self.input_files ]
        self.s_id    = [ i.name   for i in self.input_files ]
        # Now move the "run01" folders to the file name
        for folder in self.folders:
            if "run" in folder.name and folder.name.replace("run","").isdigit():
                index = self.folders.index(folder)
                self.files[index]   = folder.name+"/"+self.files[index]
                self.folders[index] = folder.parent
        # Make sure the folders are strings and reduce folder to the path behind
        if len(self.folders) != 0:
            common_prefix = pathlib.Path(os.path.commonprefix(self.folders))
            common_prefix = str(self.folders[0].parent)+"/" if common_prefix.name in str(self.folders[0]) else common_prefix
            self.folders = [ str(f).replace(str(common_prefix), "") for f in self.folders]
            # Add the experiment and simulation id to the files
            for experiment in self.experiments:
                for simulation in experiment.simulations:
                    for input_file in simulation.input_files:  
                        index_i = self.input_files.index(input_file)
                        self.e_id[index_i] = experiment.id
                        self.s_id[index_i] = simulation.marker_label.replace("$","").replace("\,"," ")
                        
        # Now get the unique folders
        self.unique_folders = list(set(self.folders))
                  
    #--------------------
    def update_ids(self):
        
        # Find the common prefix of the experiments
        experiment_ids = [experiment.id for experiment in self.experiments]
        common_prefix = ''.join(c[0] for c in takewhile(lambda x:  all(x[0] == y for y in x), zip(*experiment_ids)))
        common_prefix = common_prefix if common_prefix.endswith("_") else '_'.join(common_prefix.split('_')[0:-1])+"_"
        experiment_ids2 = [experiment.id.split("__")[0] for experiment in self.experiments]
        common_prefix2 = ''.join(c[0] for c in takewhile(lambda x:  all(x[0] == y for y in x), zip(*experiment_ids2)))
        common_prefix2 = common_prefix if common_prefix2.endswith("_") else '_'.join(common_prefix2.split('_')[0:-1])+"_"
        
        # The ID is everything behind the common prefix.
        if len(self.experiments)==1:
            for experiment in self.experiments:
                experiment.id = str(experiment.id).split("_")[0]
        if len(self.experiments)>1:
            for experiment in self.experiments:
                if common_prefix==experiment_ids[0] or common_prefix=="":
                    experiment.id = str(experiment.id).split("__")[-1].replace("nl_","").replace("th_","").split("_")[0]
                else: 
                    if "__" in common_prefix: # Then the folder names are identical
                        experiment.id = str(experiment.id).split("__")[-1]
                    if "__" not in common_prefix: # The folder names differ
                        experiment.id = str(experiment.id).split("__")[0].split(common_prefix2)[1]

        # New experiment id system: stella knob and key
        experiment_ids = []
        for experiment in self.experiments:
            inputs = experiment.simulations[0].inputParameters
            try:
                value = inputs[self.experiment_knob][self.experiment_key]
                if self.experiment_key in ['vmec_filename']:
                    experiment_ids.append(str(value))
                else:
                    variable = self.experiment_key
                    if variable=="torflux": value = str(round(np.sqrt(float(value)),2))
                    if variable=="torflux": variable = "$\\rho$"
                    if variable=="rho":     variable = "$\\rho$"
                    if variable=="delt":    variable = "$\\Delta t$"
                    if self.experiment_knob=='species_parameters_1':
                        if variable=="tprim": variable="$a/L_{Ti}$"
                        if variable=="fprim": variable="$a/L_{ne}$"
                    if self.experiment_knob=='species_parameters_2':
                        if variable=="tprim": variable="$a/L_{Te}$"
                        if variable=="fprim": variable="$a/L_{ni}$"
                    experiment_ids.append(variable + " = " + str(value))
            except: pass
        # Only use this if it produces unique experiment id's
        print(experiment_ids)
        if len(experiment_ids) != 0 and len(experiment_ids)==len(list(set(experiment_ids))): 
            for experiment in self.experiments:
                experiment.id = experiment_ids[self.experiments.index(experiment)]
        return 
    
    #------------------------
    def print_research(self):
        
        # Print some information
        printv("")
        printv("##############################")
        printv("            RESEARCH     ")
        printv("##############################")
        
        for experiment in self.experiments:
            printv("")
            printv(" Experiment "+str(self.experiments.index(experiment)+1)+": "+experiment.id)
            printv(" --------------------------------")
            if self.number_variedVariables==-1: printv("   Each simulation is assumed to be its own experiment at the following radial position:  ")
            if self.number_variedVariables==1:  printv("   There is 1 varied parameter with the following values:  ")
            if self.number_variedVariables>1:   printv("   There are "+str(self.number_variedVariables)+" varied parameters with the following values:  ")
            for varied_value in experiment.variedValues:
                printv("          "+varied_value)
            printv("   The simulations associated to this experiment are:  ")
            for simulation in experiment.simulations:
                printv("          "+'"'+simulation.id)
            printv("")
            
        # Dont collapse header on the next line
        if True: return

#=======================================
#       STYLE OF THE RESEARCH
#=======================================

    def set_labelsLinesMarkers(self):
        
        # Get the ids of the experiments
        experiment_ids = [experiment.id for experiment in self.experiments]
               
        # Otherise load the style 
        line_label, marker_label, line_style, line_color, marker_style, marker_color = load_labelsLinesMarkers(experiment_ids, self.variedValues)
        
        # Change the lines/markers of the underlying experiment objects based on self.variedValues
        for experiment in self.experiments:
            experiment.set_labelsLinesMarkers(line_label, marker_label, line_style, \
                                              line_color, marker_style, marker_color)
        
        # Dont collapse header on the next line
        if True: return
        
#============================================
#  SAVE THE SIMULATION TO THE CONFIG FOLDER
#============================================
    
    def create_defaultConfigurationFile(self):
        ''' Write a default configuration file and save it as "stella/stellapy/simulations/id.ini". 
        The goal is to use this file to reload the data in the GUI. '''
        
        # Create a configuration object
        self.configuration_file = configparser.ConfigParser() 
        
        # Get the location of the configuration file
        path_stellaGUI = pathlib.Path(os.path.dirname(os.path.abspath(__file__)).split("stellapy")[0])
        self.path_configurationFile =  path_stellaGUI / "stellapy/config/research" / str("research_"+self.id+".ini")
        
        # Rewrite the input_files to make them look orderly
        experiments = [e.id for e in self.experiments]
        experiments = "\t\n" + ("\t\n").join(experiments)
        
        # Remember which simulation this is   
        self.configuration_file['GENERAL'] = {
            'research id' : self.id, \
            'experiments': experiments}
        
        # Create a section for each experiment
        for experiment in self.experiments:
            simulations = [s.id for s in experiment.simulations]
            simulations = "\t\n" + ("\t\n").join(simulations)
            self.configuration_file[experiment.id] = {
                'simulations' : simulations}
            
            # Style the simulation
            self.configuration_file['OVERRIDE COLORS AND FONTS FOR EXPERIMENT '+experiment.id] = {
                'line_label'   : 'USE DEFAULT',\
                'line_style'   : 'USE DEFAULT',\
                'line_color'   : 'USE DEFAULT'} 
            
            # Create a section for each simulation
            for simulation in experiment.simulations:
                input_files = [str(i) for i in simulation.input_files]
                input_files = "\t\n" + ("\t\n").join(input_files)
                self.configuration_file[simulation.id] = {
                    'input files' : input_files}
                
                # Style the simulation
                self.configuration_file['OVERRIDE COLORS AND FONTS FOR SIMULATION '+simulation.id] = {
                    'marker_label' : 'USE DEFAULT',\
                    'marker_style' : 'USE DEFAULT',\
                    'marker_color' : 'USE DEFAULT'} 
            
        # Write the configuration file
        self.configuration_file.write(open(self.path_configurationFile, 'w'))
        
        # Dont collapse header on the next line
        if True: return
        
   
    
#================
# TEST THE CLASS
#================
if False:
    turnOnVerboseWrapper_configurationFile()
    folders1 = "/home/hanne/CIEMAT/RUNS/Finished/LinearMap_fprim8tprim8_Scan1to20"
    folders2 = "/home/hanne/CIEMAT/RUNS/Finished/LinearMap_fprim8tprim8_Scan1to20_1mode1file"
    folders3 = "/home/hanne/CIEMAT/RUNS/Finished/LinearMap_fprim8tprim8_ScanDt_around20"
    folders4 = ["/home/hanne/Dropbox/stella/stellapy/examples/Shot13_B348_linear/w7xr348+252_0001_oldcode",\
                "/home/hanne/Dropbox/stella/stellapy/examples/Shot13_B348_linear/w7xr348+252_0002_oldcode",\
                "/home/hanne/Dropbox/stella/stellapy/examples/Shot13_B348_linear/w7xr348+252_0003_oldcode",\
                "/home/hanne/Dropbox/stella/stellapy/examples/Shot13_B348_linear/w7xr348+252_0004_oldcode"]
    folders5 = folders4 +\
               ["/home/hanne/Dropbox/stella/stellapy/examples/Shot17_B169_linear/w7xr169+252_0001_oldcode",\
               "/home/hanne/Dropbox/stella/stellapy/examples/Shot17_B169_linear/w7xr169+252_0002_oldcode",\
               "/home/hanne/Dropbox/stella/stellapy/examples/Shot17_B169_linear/w7xr169+252_0003_oldcode",\
               "/home/hanne/Dropbox/stella/stellapy/examples/Shot17_B169_linear/w7xr169+252_0004_oldcode"]
    folders6 = ["/home/hanne/Dropbox/stella/stellapy/examples/Nonlinear_adiabatic/nl_W7xr348+s17_0001"]

if False:
    print("\n####################################")
    print("              FOLDER 1")
    print("####################################")
    research = create_research(folders=folders1, input_files=None, ignore_resolution=True,\
                number_variedVariables=1, experiment_knob="vmec_parameters", experiment_parameter="vmec_filename",\
                research_name="linearmap  20 Modes in 1 file")
    research.print_research()
    for experiment in research.experiments:
        for s in experiment.simulations:
            a = s.time_kxky
            b = s.vec_kx
            c = s.omega_kxky  
            d = s.ref_a 
            e = s.prof_n 
            f = s.q_flux_kxky
            g = s.referenceUnits 
            h = s.p_fluxQL_kxky 
            i = s.omega_avg 
            j = s.p_flux_last
            k = s.stable_modes
            
if False:
    print("\n####################################")
    print("              FOLDER 2")
    print("####################################")
    research = create_research(folders=folders2, input_files=None, ignore_resolution=True,\
                number_variedVariables=1, experiment_knob="vmec_parameters", experiment_parameter="vmec_filename",\
                research_name="linearmap  20 modes in 20 files")
    research.print_research()

if False:           
    print("\n####################################")
    print("              FOLDER 3")
    print("#####################################")
    research = create_research(folders=folders3, input_files=None, ignore_resolution=False,\
                number_variedVariables=1, experiment_knob="vmec_parameters", experiment_parameter="vmec_filename",\
                research_name="linearmap check convergence dt")
    research.print_research()

if False:     
    print("\n####################################")
    print("              FOLDER 4")
    print("#####################################")
    research = create_research(folders=folders4, input_files=None, ignore_resolution=True,\
                number_variedVariables=1, experiment_knob="vmec_parameters", experiment_parameter="vmec_filename",\
                research_name="four rhos")
    research.print_research()

if False:     
    print("\n####################################")
    print("              FOLDER 5")
    print("#####################################")
    research = create_research(folders=folders5, input_files=None, ignore_resolution=True,\
                number_variedVariables=1, experiment_knob="vmec_parameters", experiment_parameter="vmec_filename",\
                research_name="two shots")
    research.print_research()
    
if False:     
    print("\n####################################")
    print("              FOLDER 6")
    print("#####################################")
    research = create_research(folders=folders6, input_files=None, ignore_resolution=True,\
                number_variedVariables=1, experiment_knob="vmec_parameters", experiment_parameter="vmec_filename",\
                research_name="two shots")
    research.print_research()
    
    