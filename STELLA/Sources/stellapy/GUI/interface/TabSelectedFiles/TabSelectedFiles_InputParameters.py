###################################################################
#            CLASS FOR THE FRAMES "INPUT PARAMETERS"
###################################################################
''' SUBFRAME INPUTPARAMETERS ON TAB 1

Manage the "Input parameters" LabelFrame in the middle and right column of the tab 
"Simulations" which is attached to the root. Display the input parameters corresponding
to the selected simulations. Display default values in gray, values that differ
from the default in white/black and values that differ between the selected simulations
in red to emphasize that they are different. 

In the frame "All parameters" it is possible to toggle between the knobs from 
the stella code in order to look at the values of all the possible stella inputs.


Absolute path of class InputParameters:
----------------------------------------
root.tab3.class_inputParameters

'''

import tkinter as tk
from tkinter import ttk
from stellapy.GUI.graph_tools.ToolTip import ToolTip 
from stellapy.data import load_defaultInputParameters

class InputParameters():
    '''
    Calling this class initiates frame_inputs, frame_wavenumber, frame_geometry, 
    frame_resolution, frame_species, and frame_allParam on the tab "Simulations".
    See module "TabSelectedFiles" in the TabSelectedFiles package for their creation.
    
    
    Attributes
    ----------
    defaultParameters : dict[
    
    Parent widgets
    --------------
    root:               Tk()
    tabheader:          ttk.Notebook(root)
    tab1:               ttk.Frame(tabheader)          
    frame_left:         ttk.Frame(tab1)
    frame_inputs:       ttk.LabelFrame(frame_middle)  
    frame_wavenumber:   ttk.LabelFrame(frame_middle) 
    frame_geometry:     ttk.LabelFrame(frame_middle) 
    frame_species:      ttk.LabelFrame(frame_right) 
    frame_allParam:     ttk.LabelFrame(frame_right)    
    
    
    Widgets
    -------
    '''
    

#################################################################
#                          WIDGETS
#################################################################

    def __init__(self,tab3):
        '''
        Initialize the frames that will display the input parameters of the simulations
            - 
        
        Initialize the class variables that will be used by the class methods.
            -
        '''

        #====================================================================
        # Safe info from the root so that it can be passed on or used locally
        #====================================================================
        
        self.root = tab3.root 
        self.tab  = tab3

        #===========
        # VARIABLES
        #===========

        # Read the default parameters of the stella code
        self.defaultParameters = load_defaultInputParameters()  #@UndefinedVariable
        
        # Initiate widgets list for the scrollable simulations
        self.list_knobs = {}
        self.list_parameters = {}
        self.list_tooltip = {}
        self.lbl_parameter = {} 
        self.lbl_code = {} 
        self.lbl_value = {} 
        self.values_species = {}
        self.values_parameters = {}
        self.title_species = {}
        self.frame_species = {}         

        # If no simulations are selected, show the default values
        self.data_GUI = {}
        for knob in self.defaultParameters.keys():
            self.data_GUI[knob] = {}
            for key, value in self.defaultParameters[knob].items():
                self.data_GUI[knob][key] = [value] # Turn it into lists

        #===============================
        # WIDGETS CREATION FOR SUBFRAME
        #===============================

        # Tabbed display for species 1, 2, 3, ....
        self.tab_species = ttk.Notebook(self.tab.frame_species, style='species.TNotebook')
        self.initiate_speciesTab("tab 1", "  Ions  ")
        tk.Grid.rowconfigure(self.tab_species, 0, weight=0) 
        tk.Grid.columnconfigure(self.tab_species, 0, weight=1)    
        self.tab_species.grid(row=0, column=0, sticky="nsew")

        # Add a dropdown menu to the <All parameters> frame
        frame = self.tab.frame_allParam
        self.option_knob = tk.StringVar()
        options_knob = sorted(("zgrid_parameters and vpamu_grids_parameters","geo_knobs","parameters", "vmec_parameters",\
                                        "time_advance_knobs", "kt_grids_range_parameters", "physics_flags and dist_fn_knobs", "init_g_knobs", \
                                        "knobs", "species_knobs and species_parameters_1", "stella_diagnostics_knobs", "millergeo_parameters",\
                                        "layouts_knobs", "neoclassical_input", "sfincs_input", "kt_grids_box_parameters and kt_grids_knobs",\
                                        "knobs_more", "init_g_knobs_more","millergeo_parameters_more", "sfincs_input_more"))
        self.option_knob.set(options_knob[4])
        self.popupMenu_knob = ttk.OptionMenu(frame, self.option_knob, options_knob[3], *options_knob, style='gray.TMenubutton')
        tk.Grid.columnconfigure(frame, 0, weight=1)  
        tk.Grid.columnconfigure(frame, 1, weight=1)  
        tk.Grid.rowconfigure(frame,    0, weight=0) 
        tk.Grid.rowconfigure(frame,    1, weight=1) # For the parameter list
        self.popupMenu_knob.grid(row=0, column=0, sticky="nsew", columnspan=2)
        def change_dropdown(*args): # on change dropdown value
            self.update_allParameters()
        self.option_knob.trace('w', change_dropdown) # link function to change dropdown

        # Initiate the other 5 data frames
        self.initiate_data()

        #=================
        # CONFIGURE FRAME
        #=================
        tk.Grid.rowconfigure(   self.tab.frame_species, 0, weight=1)       
        tk.Grid.columnconfigure(self.tab.frame_species, 0, weight=1) 

        # Prevent indentation
        if True: return

#================================================================
# Template to initiate the data frames
#================================================================

    def initiate_dataFrame(self, frame_id, frame):

        # Configure the frame
        tk.Grid.columnconfigure(frame, 0, weight=1)    
        tk.Grid.columnconfigure(frame, 1, weight=1)  
        
        # Initiate the variables 
        self.values_species[frame_id] = []
        self.values_parameters[frame_id] = []
        self.lbl_parameter[frame_id] = []
        self.lbl_value[frame_id] = [] 

        # Get the variables
        list_tooltip    = self.list_tooltip[frame_id]
        list_parameters = self.list_parameters[frame_id]
        list_knobs      = self.list_knobs[frame_id]

        # Fill the frame
        for i in range(len(list_parameters)): 
            # Configure row
            weight = 0 if ("tab" in frame_id) else 1
            row    = i if (frame_id != "allParam") else i+1
            tk.Grid.rowconfigure(frame, i, weight=weight)          
            # Add parameter name
            oldStyle = ttk.Style()
            oldStyle.configure('code.TLabel', font=("Courier New", 12))
            textvariable = tk.StringVar(value=list_parameters[i])
            self.values_parameters[frame_id].append(textvariable)
            self.lbl_parameter[frame_id].append(ttk.Label(frame, textvariable=textvariable, style='code.TLabel', width=20))
            ToolTip(self.root, self.lbl_parameter[frame_id][i], list_tooltip[i])
            # Add corresponding value name
            textvariable = tk.StringVar(value=self.defaultParameters[list_knobs[i]][list_parameters[i]])
            self.values_species[frame_id].append(textvariable)
            self.lbl_value[frame_id].append(ttk.Label(frame, textvariable=textvariable, style='code.TLabel'))
            self.lbl_value[frame_id][i].config(foreground="gray45")
            # Add the widgets <lbl_simulation>, <chk_plotSim> and <btn_removeSim> to the frame
            self.lbl_parameter[frame_id][-1].grid(row=row, column=0, sticky="nesw", padx=(5,2), pady=6)
            self.lbl_value[frame_id][-1].grid(row=row, column=1, sticky="ne", padx=(2,5), pady=6)


        # Prevent indentation
        if True: return
        
#================================================================
# Initiate the data frames with corresponding data
#================================================================

    def initiate_speciesTab(self, tab_id, title):
        ''' Each time this is called, a new species tab is added '''

        # Data that will be displayed
        self.list_knobs[tab_id]      = ["species_parameters_1"]*6
        self.list_parameters[tab_id] = ["z", "mass", "dens", "temp", "fprim", "tprim"]
        self.list_tooltip[tab_id]    = ["Charge state", "Mass", "Density", "Temperature", "Density gradient", "Temperature gradient"]

        # For the species we make a tabbed view
        self.frame_species[tab_id] = ttk.Frame(self.tab_species, padding=(10,10,10,10))
        self.tab_species.add(self.frame_species[tab_id], text=title) 
        frame = self.frame_species[tab_id]

        # Fill the tab with data
        self.initiate_dataFrame(tab_id, frame)


    #-------------------------
    def initiate_data(self):

        # Fill the inputs frame
        self.list_knobs["inputs"]      = [ "species_knobs", "physics_flags", "time_advance_knobs"] 
        self.list_parameters["inputs"] = [ "nspec", "nonlinear", "explicit_option"]
        self.list_tooltip["inputs"]    = [ "Number of species.", \
                            "In linear simulations each mode (k_x,k_y) evolves independently, \nin non-linear simulations the modes are coupled.", \
                            "Numerical scheme for the time advance."]
        self.initiate_dataFrame("inputs", self.tab.frame_inputs)

        # Fill the wavenumber frame
        self.list_knobs["wavenumber"]      = [ "kt_grids_range_parameters" ]*6
        self.list_parameters["wavenumber"] = [ "naky", "aky_min", "aky_max", "nakx", "akx_min", "akx_max" ] 
        self.list_tooltip["wavenumber"]    = [ "Number of modes along y.", "Minimum wave number along y.", "Maximum wave number along y.", \
                            "Number of modes along x.", "Minimum wave number along x.", "Maximum wave number along x." ]
        self.initiate_dataFrame("wavenumber", self.tab.frame_wavenumber)

        # Fill the geometry frame
        self.list_knobs["geometry"]      = [ "vmec_parameters" ]*3
        self.list_parameters["geometry"] = [ "vmec_filename", "torflux", "nfield_periods" ]
        self.list_tooltip["geometry"]    = [ "?", "?", "?" ]
        self.initiate_dataFrame("geometry", self.tab.frame_geometry)

        # Fill the resolution frame
        self.list_knobs["resolution"]      = [ "zgrid_parameters"]*2 + [ "vpamu_grids_parameters" ]*2 + [ "knobs" ]*2
        self.list_parameters["resolution"] = [ "nzed", "nperiod", "nvgrid", "nmu", "delt", "nstep" ]
        self.list_tooltip["resolution"]    = [ "?", "?", "?", "?", "?", "?" ]
        self.initiate_dataFrame("resolution", self.tab.frame_resolution)

        # Fill the all parameters frame
        self.list_knobs["allParam"]      = ["knobs"]*15
        self.list_parameters["allParam"] = [ "nstep", "delt", "fphi", "fapar", "fbpar", "delt_option", "zed_upwind", \
                                                        "vpa_upwind", "time_upwind", "avail_cpu_time",\
                                                        "cfl_cushion", "delt_adjust", "mat_gen", "mat_read", "fields_kxkyz"]
        self.list_tooltip["allParam"]    = ["?"]*15
        self.initiate_dataFrame("allParam", self.tab.frame_allParam)

        # Prevent indentation
        if True: return
        
#################################################################
#                          METHODS
#################################################################

#================================================================
# Update the frames
#================================================================

    def update_frame(self):
        self.update_simulations() 
        self.update_dataFrames(["inputs", "wavenumber", "geometry", "resolution"])
        self.update_speciesFrame()
        self.update_allParameters()

    #--------------------
    def update_simulations(self):
        ''' Collect the data from the input files, where we overwrite the default values '''
        self.data_GUI = {}
        for experiment in self.root.Research.experiments:
            for simulation in experiment.simulations:
                for input_file in simulation.input_files + ["everything"]:
                    if input_file=="everything": data_simulation = simulation.inputParameters
                    if input_file!="everything": data_simulation = simulation.inputs[input_file]
                    for knob in data_simulation.keys():
                        for key, value in data_simulation[knob].items():
                            if knob in self.data_GUI.keys():
                                if key in self.data_GUI[knob].keys():
                                    self.data_GUI[knob][key].append(data_simulation[knob][key])
                                else: 
                                    self.data_GUI[knob][key] = [data_simulation[knob][key]]
                            else:
                                self.data_GUI[knob] = {}
                                self.data_GUI[knob][key] = [data_simulation[knob][key]]
                                
        # If no simulations are selected, show the default values
        if self.data_GUI == {}: 
            for knob in self.defaultParameters.keys():
                self.data_GUI[knob] = {}
                for key, value in self.defaultParameters[knob].items():
                    self.data_GUI[knob][key] = [value] # Turn it into lists
        return 
    
    #----------------------------------
    def update_dataFrames(self, frame_ids):
        for frame_id in frame_ids:
            numberOfWidgets = len(self.values_species[frame_id])
            for i in range(numberOfWidgets):
                knob      = self.list_knobs[frame_id][i]
                parameter = self.list_parameters[frame_id][i]
                default_text = self.defaultParameters[knob][parameter]
                try: current_list = self.data_GUI[knob][parameter]
                except: current_list = [default_text]
                current_list = sorted(list(set(current_list))) # Only look at unique values
                if len(current_list) == 1 and current_list[0] == default_text: # If it is the default value, replace it
                    self.values_species[frame_id][i].set(str(current_list[0])) 
                    self.lbl_value[frame_id][i].config(foreground="gray45")
                if len(current_list) != 1 or current_list[0] != default_text: # If it is anything but the default value, change color and replace
                    if len(current_list) == 1:
                        self.values_species[frame_id][i].set(str(current_list[0])) # Print a value
                        self.lbl_value[frame_id][i].config(foreground=self.root.color['fg'])
                    if len(current_list) > 1:
                        text = '[%s]' % ', '.join(map(str, current_list))
                        if len(text) > 17:
                            splitted_text = text.split(",") # split into list of individual items
                            splitted_text = "\n".join([", ".join(splitted_text[i:i+5]) for i in range(0,len(splitted_text),5)])
                            ToolTip(self.root, self.lbl_value[frame_id][i], splitted_text)
                            text = text[0:15]+" ..."
                        self.values_species[frame_id][i].set(text) # Remove quotes when printing list
                        self.lbl_value[frame_id][i].config(foreground="red")
                
        return
    
    #--------------------------
    def update_speciesFrame(self):

        # If there are more species
        frame_ids = list(self.lbl_parameter.keys())
        tab_ids = [tab for tab in frame_ids if ("tab" in tab) ]
        nspec = self.data_GUI["species_knobs"]["nspec"]
        nspec = sorted(list(set(nspec)))            # Only look at unique values
        max_nspec = max(max(nspec), 2)              # Even with 1 specie we have adiabatic electrons set through nine and tite
        try: unique_folders = self.root.Research.unique_folders
        except: unique_folders = []
        if len(unique_folders)==0: max_nspec=1      # When all simulations are removed, remove the second species tab
        # Add more tabs if there aren't enough
        if len(tab_ids) < max_nspec: 
            for i in range(len(tab_ids), max_nspec):
                tab_id = "tab "+str(i+1)
                self.initiate_speciesTab(tab_id, "Species "+str(i+1))
                self.list_knobs[tab_id]      = ["species_parameters_"+str(i+1)]*6
                self.list_parameters[tab_id] = self.list_parameters["tab 1"]
                self.list_tooltip[tab_id]    = self.list_tooltip["tab 1"]
        # Remove tabs if there are too many
        if max_nspec < len(tab_ids): 
            for i in range(max_nspec, len(tab_ids)):
                self.tab_species.forget(self.frame_species["tab "+str(i+1)])
                self.frame_species["tab "+str(i+1)].destroy()
                del self.lbl_parameter["tab "+str(i+1)]                                      # Because we get the ids from here
        # Update the tab ids
        frame_ids = list(self.lbl_parameter.keys())
        tab_ids = [tab for tab in frame_ids if ("tab" in tab) ]

        # Change the name if the tabs or named wrongly
        if not max_nspec==1:
            if len(nspec) == 1 and nspec[0] == 1: # Assume there are ions and adiabatic electrons
                self.tab_species.tab(1, text = 'Adiabatic electrons')
                # If it was set for kinetic electrons before, over write it
                if "nine" not in self.list_parameters["tab 2"]:
                    self.list_knobs["tab 2"]      = ["species_parameters_a"]*6
                    self.list_parameters["tab 2"] = ["z", "mass", "dens", "temp", "nine", "tite"]
                    self.list_tooltip["tab 2"]    = ["Charge state", "Mass", "Density", "Temperature", "Density ratio", "Temperature ratio"]
                    self.values_parameters["tab 2"][4].set("nine") 
                    self.values_parameters["tab 2"][5].set("tite") 
            elif len(nspec) == 1 and nspec[0] == 2: # Assume second species is always kinetic electrons 
                self.tab_species.tab(self.frame_species["tab 2"], text = 'Kinetic electrons')
                # If it was set for adiabatic electrons before, over write it
                if "nine" in self.list_parameters["tab 2"]:
                    self.list_knobs["tab 2"]      = ["species_parameters_2"]*6
                    self.list_parameters["tab 2"] = ["z", "mass", "dens", "temp", "fprim", "tprim"]
                    self.list_tooltip["tab 2"]    = ["Charge state","Mass","Density","Temperature","Density gradient","Temperature gradient"]
                    self.values_parameters["tab 2"][4].set("fprim") 
                    self.values_parameters["tab 2"][5].set("tprim") 
            elif len(nspec) == 1 and nspec[0] == 3: # Assume second species is always kinetic electrons and third impurities
                self.tab_species.tab(self.frame_species["tab 2"], text = 'Kinetic electrons')
                self.tab_species.tab(self.frame_species["tab 3"], text = 'Impurities')
                # If it was set for adiabatic electrons before, over write it
                if "nine" in self.list_parameters["tab 2"]:
                    self.list_knobs["tab 2"]      = ["species_parameters_2"]*6
                    self.list_parameters["tab 2"] = ["z", "mass", "dens", "temp", "fprim", "tprim"]
                    self.list_tooltip["tab 2"]    = ["Charge state","Mass","Density","Temperature","Density gradient","Temperature gradient"]
                    self.values_parameters["tab 2"][4].set("fprim") 
                    self.values_parameters["tab 2"][5].set("tprim") 
            elif len(nspec) == 2 and (1  in nspec and 2 in nspec): # Assume we have both adiabtic and kinetic simulations
                self.tab_species.tab(self.frame_species["tab 2"], text = 'Adiabatic or kinetic electrons')
                # If it was set for adiabatic electrons before, over write it
                if "nine" in self.list_parameters["tab 2"]:
                    self.list_knobs["tab 2"]      = ["species_parameters_2"]*6
                    self.list_parameters["tab 2"] = ["z", "mass", "dens", "temp", "fprim", "tprim"]
                    self.list_tooltip["tab 2"]    = ["Charge state","Mass","Density","Temperature","Density gradient","Temperature gradient"]
                    self.values_parameters["tab 2"][4].set("fprim") 
                    self.values_parameters["tab 2"][5].set("tprim") 
            else:
                print("Warning: The current configuration of species can not be printed yet in the species screen.", len(nspec), nspec)
        # Update the data in the frames
        self.update_dataFrames(tab_ids)
    

    #------------------------------  
    def update_allParameters(self):

        # The know that needs to be shown
        knob = self.option_knob.get()

        # If "more" in know show parameters after the first 15
        if "more" in knob:
            knob = knob.split("_more")[0]; switch=30
            parameters = sorted(list(self.data_GUI[knob].keys()))[14:-1]
        elif " and " in knob:
            knob1 = knob.split(" and ")[0]; knob2 = knob.split(" and ")[1]; knob=knob1
            parameters = sorted(list(self.data_GUI[knob1].keys())) + [" "]; switch=len(parameters)
            parameters = parameters + sorted(list(self.data_GUI[knob2].keys()))
        else:
            parameters = sorted(list(self.data_GUI[knob].keys())); switch=30
        # Add data to frame
        frame_id = "allParam"
        numberOfWidgets = 15
        numberOfParameters = min(len(parameters), 15)
        # Replace current data with correct data
        for i in range(numberOfParameters):
            if switch:
                if i==switch-1:  knob=" "
                if i==switch:    knob=knob2
            default_text = self.defaultParameters[knob][parameters[i]]
            current_list = self.data_GUI[knob][parameters[i]]
            current_list = sorted(list(set(current_list))) # Only look at unique values
            if len(current_list) == 1 and current_list[0] == default_text: # If it is the default value, replace since we show another knob
                self.values_parameters[frame_id][i].set(str(parameters[i])) 
                self.values_species[frame_id][i].set(str(current_list[0])) 
                self.lbl_value[frame_id][i].config(foreground="gray45")
            if len(current_list) != 1 or current_list[0] != default_text: # If it is anything but the default value, change color and replace
                if len(current_list) == 1:
                    self.values_parameters[frame_id][i].set(str(parameters[i])) 
                    self.values_species[frame_id][i].set(str(current_list[0])) 
                    self.lbl_value[frame_id][i].config(foreground=self.root.color['fg'])
                if len(current_list) > 1:
                    text = '[%s]' % ', '.join(map(str, current_list))
                    if len(text) > 15:
                        splitted_text = text.split(",") # split into list of individual items
                        splitted_text = "\n".join([", ".join(splitted_text[i:i+5]) for i in range(0,len(splitted_text),5)])
                        ToolTip(self.root, self.lbl_value[frame_id][i], splitted_text)
                        text = text[0:15]+" ..."
                    self.values_parameters[frame_id][i].set(str(parameters[i])) 
                    self.values_species[frame_id][i].set(text) # Remove quotes when printing list
                    self.lbl_value[frame_id][i].config(foreground="red")
        # Remove text on excessive labels
        for i in range(numberOfParameters, numberOfWidgets):
            self.values_parameters[frame_id][i].set(" ")
            self.values_species[frame_id][i].set(" ") 
        return























