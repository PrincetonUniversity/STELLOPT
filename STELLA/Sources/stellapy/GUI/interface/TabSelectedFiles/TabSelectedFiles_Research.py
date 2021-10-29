####################################################################
#              CLASS FOR THE FRAME "RESEARCH"
####################################################################
''' SUBFRAME RESEARCH ON TAB 1 '''

import gc
import pickle
import tkinter as tk
from tkinter import ttk
from bisect import bisect
from stellapy.utils import initiate_nesteddict
from stellapy.config import CONFIG
from stellapy.GUI.interface import update_GUI
from stellapy.GUI.graph_tools.ScrollableFrame import ScrollableFrame
from stellapy.simulations.Research import create_research
from stellapy.config import turnOnVerboseWrapper_configurationFile, turnOffVerboseWrapper_configurationFile #@unresolvedimport
from stellapy.GUI.graph_tools import PAD_TITLE2, PAD_LABEL2, PAD_ENTRY2  #@UnusedImport 

####################################################################
# INITIALIZE THE FRAME WITH INFORMATION ON THE SELECTED SIMULATIONS
####################################################################

class Research():
    '''
    Calling this class initiates the "Research" LabelFrame on the middle left side 
    of the tab "Simulations" attached to the root.
        

    Parent widgets
    --------------
    root:               Tk()
    tabheader:          ttk.Notebook(root)
    tab1:               ttk.Frame(tabheader)          
    frame_left:         ttk.Frame(tab1)
    frame_research:     ttk.LabelFrame(frame_left)    --> frame

    '''

#################################################################
#                          WIDGETS
#################################################################

    def __init__(self,tab1):
        '''
        Initialize the widgets in the "Research" LabelFrame of the tab "Simulations".
            - A scrollable frame to show the experiments and simulations
            - A checkbox to ignore the resolution
            - A entry box to choose the name of the research (for reloading)
            - A button to load previous researches
            - A button to save the current research (will be quite automatic tho)

        '''

        #====================================================================
        # Safe info from the root so that it can be passed on or used locally
        #====================================================================
        
        self.tab1 = tab1            
        self.root = tab1.root 
        self.frame = tab1.frame_research
                    
        #===========
        # VARIABLES
        #===========
        
        # Basic variables
        self.research_name = "DEFAULT"
        self.ignore_resolution = True
        self.number_variedVariables = 1
        self.folderIsExperiment = False
        self.experiment_knob = "vmec_parameters"
        self.experiment_key = "vmec_filename"
        self.tree_experiments = []
        
        # Dropdown menus
        self._options_resolution = ["   Ignore resolution", "   Include resolution"]
        self._options_variedParam = ["  0 varied parameters",\
                                     "  1 varied parameter", "  2 varied parameters", "  3 varied parameters",
                                     "  4 varied parameters", "  5 varied parameters"]
        
        # Knob and key options
        self.options_keys = {
            "zgrid_parameters" : \
                ['nzed', 'nperiod', 'ntubes', 'boundary_option', \
                 'zed_equal_arc', 'shat_zero', 'nzgrid'],\
            "geo_knobs" : \
                ['geo_option', 'overwrite_bmag', 'overwrite_gradpar', 'overwrite_gds2', \
                 'overwrite_gds21', 'overwrite_gds22', 'overwrite_gds23', 'overwrite_gds24', \
                 'overwrite_gbdrift', 'overwrite_cvdrift', 'overwrite_gbdrift0', 'geo_file'],\
            "vmec_parameters" : \
                ['vmec_filename', 'alpha0', 'zeta_center', 'nfield_periods', 'torflux', \
                 'surface_option', 'verbose', 'zgrid_scalefac', 'zgrid_refinement_factor'],\
            "parameters" : \
                ['beta', 'vnew_ref', 'rhostar', 'zeff', 'tite', 'nine' ],\
            "vpamu_grids_parameters" : 
                ['nvgrid', 'vpa_max', 'nmu', 'vperp_max', 'equally_spaced_mu_grid'],\
            "dist_fn_knobs" : 
                ['adiabatic_option'],\
            "time_advance_knobs" : 
                ['explicit_option', 'xdriftknob', 'ydriftknob', 'wstarknob', 'flip_flo'],\
            "kt_grids_knobs" : 
                ['grid_option'],\
            "kt_grids_box_parameters" : 
                ['nx', 'ny', 'dkx', 'dky', 'jtwist', 'y0', 'naky', 'nakx'],\
            "kt_grids_range_parameters" : 
                ['nalpha', 'naky', 'nakx', 'aky_min', 'aky_max', 'akx_min', \
                 'akx_max', 'theta0_min', 'theta0_max'],\
            "physics_flags" : 
                ['full_flux_surface', 'include_mirror', 'nonlinear', \
                 'include_parallel_nonlinearity', 'include_parallel_streaming'],\
            "init_g_knobs" : 
                ['tstart', 'scale', 'ginit_option', 'width0', 'refac', 'imfac', \
                 'den0', 'par0', 'tperp0', 'den1', 'upar1', 'tpar1', 'tperp1', \
                 'den2', 'upar2', 'tpar2', 'tperp2', 'phiinit', 'zf_init', 'chop_side',\
                 'left', 'even', 'restart_file', 'restart_dir', 'read_many'],\
            "knobs" : 
                ['nstep', 'delt', 'fapar', 'fbpar', 'delt_option', 'zed_upwind', \
                 'vpa_upwind', 'time_upwind', 'avail_cpu_time', 'cfl_cushion',\
                 'delt_adjust', 'mat_gen', 'mat_read', 'fields_kxkyz', 'stream_implicit', \
                 'mirror_implicit', 'driftkinetic_implicit',  'mirror_semi_lagrange', \
                 'mirror_linear_interp', 'maxwellian_inside_zed_derivative', 'stream_matrix_inversion'],\
            "species_knobs" : 
                ['nspec', 'species_option'],\
            "species_parameters_1" : 
                ['z', 'mass', 'dens', 'tprim', 'fprim', 'd2ndr2', 'd2Tdr2', 'type' ],\
            "species_parameters_2" : 
                ['z', 'mass', 'dens', 'tprim', 'fprim', 'd2ndr2', 'd2Tdr2', 'type' ],\
            "species_parameters_3" : 
                ['z', 'mass', 'dens', 'tprim', 'fprim', 'd2ndr2', 'd2Tdr2', 'type' ],\
            "species_parameters_4" : 
                ['z', 'mass', 'dens', 'tprim', 'fprim', 'd2ndr2', 'd2Tdr2', 'type' ],\
            "stella_diagnostics_knobs" : 
                ['nwrite', 'navg', 'nmovie', 'nsave', 'save_for_restart', 'write_omega', \
                 'write_phi_vs_time', 'write_gvmus', 'write_gzvs', 'write_kspectra', \
                 'write_moments', 'flux_norm', 'write_fluxes_kxky' ],\
            "millergeo_parameters" : 
                ['rhoc', 'rmaj', 'shift', 'qinp', 'shat', 'kappa', 'kapprim', 'tri',\
                 'triprim', 'rgeo', 'betaprim', 'betadbprim', 'd2qdr2', 'd2psidr2', \
                 'nzed_local', 'read_profile_variation', 'write_profile_variation'],\
            "layouts_knobs" : 
                ['xyzs_layout', 'vms_layout'],\
            "neoclassical_input" : 
                ['include_neoclassical_terms', 'nradii', 'drho', 'neo_option'],\
            "sfincs_input" : \
                ['read_sfincs_output_from_file', 'nproc_sfincs', 'irad_min', 'irad_max', 'calculate_radial_electric_field',\
                'includeXDotTerm', 'includeElectricFieldTermInXiDot', 'magneticDriftScheme', 'includePhi1',\
                'includePhi1InKineticEquation', 'geometryScheme', 'VMECRadialOption', 'equilibriumFile',\
                'coordinateSystem', 'inputRadialCoordinate', 'inputRadialCoordinateForGradients', 'aHat',\
                'psiAHat', 'Delta', 'nu_n', 'dPhiHatdrN', 'Er_window', 'nxi', 'nx', 'Ntheta', 'Nzeta']}
        self.options_knobs = list(self.options_keys.keys())
        
        #============================
        # WIDGETS CREATION FOR FRAME
        #============================
        
        # The list of researches will be displayed inside a treeview widget
        self.tree = ttk.Treeview(self.frame, show="tree", columns=("#0"))
        style = ttk.Style()
        style.configure("Treeview.Heading", font=(None, 100))

        # "Save research" button which allows to set the name of the research.
        self.button_saveResearch = ttk.Button(self.frame, text="Save research", width=14)
        self.button_saveResearch.config(command=lambda: self.save_research())

        # "Load research" button to choose one of the previous researches
        self.button_loadResearch = ttk.Button(self.frame, text="Load research", width=14)
        self.button_loadResearch.config(command=lambda: self.load_research())

        # "Ignore resolution" and "Include resolution" option
        self.popup_resolution = tk.Menu(self.frame, tearoff=0)
        self.popup_resolution.add_command(label="Ignore resolution", command= lambda: update_btn("Ignore resolution"))
        self.popup_resolution.add_command(label="Include resolution", command= lambda: update_btn("Include resolution")) 
        self.button_resolution = ttk.Button(self.frame, text="Ignore resolution", width=16)        
        def update_btn(x):
            # When an option is chosen from the popup menu, change the text on the button
            self.button_resolution.config(text=x)
            # Get the setting for ignoring the resolution
            if x=="Ignore resolution":  self.ignore_resolution = True
            if x=="Include resolution": self.ignore_resolution = False
            # Remake the simulation object with the new setting
            self.create_researchObject()
            # Update the GUI
            update_GUI(self.root)
        def btn_popup(event):
            try:
                self.popup_resolution.tk_popup(event.x_root, event.y_root, 0)
            finally:
                self.popup_resolution.grab_release()
        self.button_resolution.bind("<Button-1>", btn_popup)
        
        # "X varied parameters" option        
        self.popup_variedParam = tk.Menu(self.frame, tearoff=0)
        self.popup_variedParam.add_command(label="0 varied parameters", command= lambda: update_btn2("0 varied parameters")) 
        self.popup_variedParam.add_command(label="1 varied parameter",  command= lambda: update_btn2("1 varied parameter")) 
        self.popup_variedParam.add_command(label="2 varied parameters", command= lambda: update_btn2("2 varied parameters")) 
        self.popup_variedParam.add_command(label="3 varied parameters", command= lambda: update_btn2("3 varied parameters")) 
        self.popup_variedParam.add_command(label="4 varied parameters", command= lambda: update_btn2("4 varied parameters")) 
        self.popup_variedParam.add_command(label="5 varied parameters", command= lambda: update_btn2("5 varied parameters"))
        self.button_variedParam = ttk.Button(self.frame, text="1 varied parameter", width=16)  
        self.button_variedParam.bind("<Button-3>", self.popup_menu)      
        def update_btn2(x):
            # When an option is chosen from the popup menu, change the text on the button
            self.button_variedParam.config(text=x)
            # Get the setting for ignoring the resolution
            if x=="0 varied parameters":        self.number_variedVariables = 0
            if x=="1 varied parameter":         self.number_variedVariables = 1
            if x=="2 varied parameters":        self.number_variedVariables = 2
            if x=="3 varied parameters":        self.number_variedVariables = 3
            if x=="4 varied parameters":        self.number_variedVariables = 4
            if x=="5 varied parameters":        self.number_variedVariables = 5
            # Remake the simulation object with the new setting
            self.create_researchObject()
            # Update the GUI
            update_GUI(self.root)
        def btn_popup2(event):
            try:
                self.popup_variedParam.tk_popup(event.x_root, event.y_root, 0)
            finally:
                self.popup_variedParam.grab_release()
        self.button_variedParam.bind("<Button-1>", btn_popup2)
    
        
        #=================
        # CONFIGURE FRAME
        #=================
        tk.Grid.rowconfigure(   self.frame, 0, weight=1) # Scrollable canvas
        tk.Grid.rowconfigure(   self.frame, 1, weight=0) # 3 buttons to edit the simulation selection
 
        tk.Grid.columnconfigure(self.frame, 0, weight=1) # Column for button_saveResearch
        tk.Grid.columnconfigure(self.frame, 1, weight=1) # Column for button_loadResearch
        tk.Grid.columnconfigure(self.frame, 2, weight=1) # Column for button_resolution
        tk.Grid.columnconfigure(self.frame, 3, weight=1) # Column for button_variedParam

        #======================
        # WIDGETS ARRANGEMENT
        #=====================
        self.tree.grid(                in_=self.frame, row=0, column=0, padx=8, pady=4, sticky='nesw', columnspan=4)
        self.button_saveResearch.grid( in_=self.frame, row=1, column=0, padx=8, pady=4, sticky='nw', ipady=7)
        self.button_loadResearch.grid( in_=self.frame, row=1, column=1, padx=8, pady=4, sticky='nw', ipady=7) 
        self.button_resolution.grid(   in_=self.frame, row=1, column=2, padx=8, pady=4, sticky='nw', ipady=7) 
        self.button_variedParam.grid(  in_=self.frame, row=1, column=3, padx=8, pady=4, sticky='nw', ipady=7) 
        
        #===========================================
        # POPUP MENU ON X VARIED PARAMETERS BUTTON
        #============================================
        self.rightclick_menu = tk.Menu(self.root, tearoff=0)
        self.rightclick_menu.add_command(label='Add "1 folder=1 experiment" option', command = lambda: self.update_menu(1))
        self.rightclick_menu.add_command(label='Remove "1 folder=1 experiment" option', command = lambda: self.update_menu(2))
        self.rightclick_menu.add_command(label='Sort experiments by knob and key', command = lambda: self.update_menu(3))

        #============
        # TREE VIEW
        #============
        
        # Add columns  
        self.tree.heading("#0", text="Simulation",anchor=tk.W) 
        self.tree["displaycolumns"] = ("#0")

        # When the tab is visible, make sure the columns have a good size
        def resize_columns(*args): 
            self.tree.column("#0", width=int(self.tree.winfo_width()), stretch=tk.YES)
        self.frame.bind("<Visibility>", resize_columns)
    
        # Prevent folding of header comment on next line
        if True: return

#==================================================
# For the popup menu on x varied parameters button
#===================================================

    def popup_menu(self, event):
        ''' Rightclick event linked to the scrollable canvas. '''
        try:        self.rightclick_menu.tk_popup(event.x_root, event.y_root)
        finally:    self.rightclick_menu.grab_release()
        return
        
    #------------------------
    def update_menu(self, i):
        if i==1: self.folderIsExperiment = True
        if i==2: self.folderIsExperiment = False
        if i==3: self.openMenuToChooseKnobAndKey()
        if len(self.root.input_files) != 0 and i in [1,2]: 
            self.create_researchObject()
            update_GUI(self.root)
        return 
    
    #---------------------
    def openMenuToChooseKnobAndKey(self):
        
        # Create a top window to ask for the parameter knob and key
        research_window = tk.Toplevel(bg=self.root.color['bg'])
        research_window.title("Sort experiments by stella knob and key")
        research_window.withdraw()
        self.root.eval(f'tk::PlaceWindow {str(research_window)} center')
        
        # When updating the stella knob, set the options list for the stella key
        def update_key(*args):
            knob = self.var_knob.get()
            self.mnu_key['menu'].delete(0, 'end')
            for knob in self.options_keys[knob]:
                self.mnu_key['menu'].add_command(label=knob, command=tk._setit(self.var_key, knob))

        # Choose the stella knob
        self.lbl_title  = ttk.Label(research_window, text="Choose the stella knob and key:")
        self.var_knob = tk.StringVar(value=self.options_knobs[2]); width=30; 
        self.lbl_knob = ttk.Label(research_window, text="    Knob: ", style='prefTitle.TLabel')
        self.mnu_knob = ttk.OptionMenu(research_window, self.var_knob, self.options_knobs[2], *self.options_knobs, style='option.TMenubutton')
        self.mnu_knob["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
        self.mnu_knob.config(width=width)
        self.var_knob.trace('w', update_key) # link function to a change of the dropdown options
        
        # Choose the stella key
        knob = self.var_knob.get()
        self.var_key = tk.StringVar(value=self.options_keys[knob][0])
        self.lbl_key = ttk.Label(research_window, text="    Key: ", style='prefTitle.TLabel') 
        self.mnu_key = ttk.OptionMenu(research_window, self.var_key, self.options_keys[knob][0], *self.options_keys[knob], style='option.TMenubutton')
        self.mnu_key["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
        self.mnu_key.config(width=width)

        # Configure the frame
        tk.Grid.rowconfigure(   research_window, 0, weight=0) 
        tk.Grid.rowconfigure(   research_window, 1, weight=0) 
        tk.Grid.rowconfigure(   research_window, 2, weight=0) 
        tk.Grid.columnconfigure(research_window, 0, weight=0)  
        tk.Grid.columnconfigure(research_window, 0, weight=0)  
        
        # Place the widgets in the frame
        self.lbl_title.grid(row=0, column=0, padx=50,     pady=(50,5), sticky='nesw', ipady=2, ipadx=5, columnspan=2)
        self.lbl_knob.grid( row=1, column=0, padx=(50,5), pady=(5,5),  sticky='nesw', ipady=2, ipadx=5)
        self.mnu_knob.grid( row=1, column=1, padx=(5,50), pady=(5,5),  sticky='nesw', ipady=2, ipadx=5)
        self.lbl_key.grid(  row=2, column=0, padx=(50,5), pady=(5,50),  sticky='nesw', ipady=2, ipadx=5)
        self.mnu_key.grid(  row=2, column=1, padx=(5,50), pady=(5,50), sticky='nesw', ipady=2, ipadx=5)

        # Closing event
        def on_closing(*args):
            # Save the stella knob and key
            self.experiment_knob = self.var_knob.get()
            self.experiment_key  = self.var_key.get() 
            # Create a new research object
            if len(self.root.input_files) != 0: 
                self.create_researchObject()
                update_GUI(self.root)
            # Close the window
            research_window.destroy()
        research_window.protocol("WM_DELETE_WINDOW", on_closing)
        research_window.deiconify()
        if True: return 
    
#====================
# Research object
#====================

    def create_researchObject(self):
        
        # First delete the previous object
        try: data = self.root.Research.data
        except: data = {}
        try: del self.root.Research; gc.collect()
        except: pass
                
        # For the selected files, create a Research object
        self.root.Research = create_research(\
            # To create the simulations we need their location and whether to ignore the resolution
            folders=None, input_files=self.root.input_files, ignore_resolution = self.ignore_resolution, \
            # To group them by experiment we need to know how many variables differ between the simulations
            # Or which variable is unique for each experiment
            number_variedVariables=self.number_variedVariables, experiment_knob=self.experiment_knob, experiment_key=self.experiment_key,\
            # Give the research a name in order to save it as a configuration file
            research_name = self.research_name, folderIsExperiment = self.folderIsExperiment,\
            # Save some information from the plots
            data=data)
        
        # Print some information to the command prompt
        turnOnVerboseWrapper_configurationFile()
        self.root.Research.print_research()
        turnOffVerboseWrapper_configurationFile()
        
        
    #---------------------------
    def save_research(self):
        
        # Create a top window to ask for the research name
        research_window = tk.Toplevel(bg=self.root.color['bg'])
        research_window.title("Save research")
        self.root.eval(f'tk::PlaceWindow {str(research_window)} center')
        
        # Label and entry for the research
        font = ("Courier New", 11); width=10; 
        self.var_research  = tk.StringVar(value=self.research_name)
        self.lbl_research  = ttk.Label(research_window, text="Choose the research name:")
        self.ent_research  = ttk.Entry(research_window, textvariable=self.var_research, font=font,  style='opt_valueCBold.TEntry', width=width)
        
        # Configure the frame
        tk.Grid.rowconfigure(   research_window, 0, weight=1) 
        tk.Grid.rowconfigure(   research_window, 1, weight=1) 
        tk.Grid.columnconfigure(research_window, 0, weight=1)  
        
        # Place the widgets in the frame
        self.lbl_research.grid(row=0, column=0, padx=50, pady=(50,5), sticky='nesw', ipady=7, ipadx=10)
        self.ent_research.grid(row=1, column=0, padx=(50,50), pady=(5,50), sticky='nesw', ipady=2, ipadx=5)
        
        # Closing event
        def on_closing(*args):
            if self.var_research.get() != "DEFAULT" and self.root.input_files != []:
                # Create a new research object
                self.research_name = self.var_research.get().replace(" ", "_")
                self.create_researchObject()   
                self.research_name = "DEFAULT"
                # Save it with pickle
                pickle_path = CONFIG['PATHS']['stellapy']+"/config/research/"+str("research_"+self.root.Research.id+".pickle")
                pickly_file = open(pickle_path, "wb")
                pickle.dump(self.root.Research, pickly_file)            
                pickly_file.close()
            # Close the window
            research_window.destroy()
        research_window.protocol("WM_DELETE_WINDOW", on_closing)
        self.ent_research.bind('<Return>', on_closing)
        
    #-----------------------
    def load_research(self):
        print("LOAD RESEARCH")
        # Choose files and start selection in the standard run folder.
        title           = "Select pickle file"
        filetypes       = (("in files","*.pickle"),("all files","*.*")) 
        initialdir      = CONFIG['PATHS']['stellapy'] + "/config/research/"
        selected_files  = tk.filedialog.askopenfilenames(initialdir=initialdir,title=title,filetypes=filetypes)

        # Only continue if a file was selected
        if len(selected_files) > 0 and selected_files!=None:
            
            # Remove the input_files that are displayed now
            self.root.tab_Simulations.class_simulations.clear_simulations()
        
            # Open the file with pickle
            pickly_file = open(selected_files[0], 'rb')
            self.root.Research = pickle.load(pickly_file)            
            pickly_file.close()
            
            # Display the current input_files
            for experiment in self.root.Research.experiments:
                for simulation in experiment.simulations:
                        self.root.input_files += simulation.input_files 
            self.root.tab_Simulations.class_simulations.update_treeView()
            
            # Update the resolution button
            self.ignore_resolution = self.root.Research.ignore_resolution
            if self.ignore_resolution == True:  x=  "Ignore resolution" 
            if self.ignore_resolution == False: x = "Include resolution"
            self.button_resolution.config(text=x)
            
            # Update the varied parameters button 
            self.number_variedVariables = self.root.Research.number_variedVariables 
            if self.number_variedVariables==-1: x = "No grouping"      
            if self.number_variedVariables==0:  x = "0 varied parameters" 
            if self.number_variedVariables==1:  x = "1 varied parameter"  
            if self.number_variedVariables==2:  x = "2 varied parameters"
            if self.number_variedVariables==3:  x = "3 varied parameters" 
            if self.number_variedVariables==4:  x = "4 varied parameters"
            if self.number_variedVariables==5:  x = "5 varied parameters"
            self.button_variedParam.config(text=x)
            
            # Remember the stella kob and key
            self.experiment_knob = self.root.Research.experiment_knob
            self.experiment_key  = self.root.Research.experiment_key
            
            # Update the GUI   
            update_GUI(self.root)
            return 
        
        # Prevent the collapsing of the header comment on the next line
        if True: return
        
#====================
# Update canvas
#====================

    def clear_simulations(self):
        ''' Remove all simulations from the research frame. '''
        
        # Remove all the items from the tree view
        self.tree.delete(*self.tree.get_children())
        self.tree_experiments = []
        
        # Remove the Research object
        self.root.Research = type('Dummy', (object,), {'content':{}})()
        self.root.Research.experiments = []
        return

    #----------------------------------
    def update_treeView(self):
        '''
        Update the scrollable canvas which displays the selected input files, 
        the selection is sorted by the parents folders. First make a frame
        and title label for each folder, next add the label, cross button and 
        checkmark widget for each input file or overwrite the existing labels
        with the correct data, overwriting reduces flickering of the GUI.
        
        Attributes
        ----------
        input_files: dict[folder][input_file] = tk.IntVar()
            Go thorough the simulations to display them on the GUI
            
        widgets : dict['folder']['title', 'label', 'button', 'check'] = list of widgets; dict['folder']['var_lbl'] = list of tk.StringVar
            Stores the newly created widgets for the simulations 
            
        frame_folders: dict['folder'] = ttk.Frame(master=self.scrollableCanvas.scrollable_frame)
            Stores the newly created widgets for the folders
        '''

        # Only update when there are input files
        if len(self.root.input_files) != 0:
            
            # Remove all the items from the tree view
            self.tree.delete(*self.tree.get_children())
            self.tree_experiments = []
            
            # Create the research object
            self.create_researchObject()
    
            # Iterate over the folders
            for experiment in self.root.Research.experiments:

                # Check whether the experiment is already in the treeview
                # If it is not, add the experiment to the treeview
                if not self.tree.exists(experiment.id):
                    contents = [self.tree.item(child)["text"] for child in self.tree.get_children("")]
                    self.tree_experiments.append(self.tree.insert(parent="", index=bisect(contents, experiment.id), iid=experiment.id, text="Experiment: "+experiment.id, values=("","")))
                    self.tree.item(experiment.id, open=True, tags=['bold'])
                    self.tree.tag_configure("bold", font=(None, 11, 'bold'))
                    
                # Add the simulations to the experiment
                for tree_experiment in self.tree_experiments:
                    self.tree.selection_set(tree_experiment)  
                    experiment_iid = self.tree.selection()[0]
                    if experiment_iid==experiment.id:
                        for simulation in experiment.simulations:
                            # If the simulation is not in the treeview, add it
                            if not self.tree.exists(simulation.id): 
                                i = experiment.simulations.index(simulation)
                                text = experiment.variedValues[i].replace("\\","").replace("$","").replace(",", " ")
                                contents = [self.tree.item(child)["text"] for child in self.tree.get_children(experiment.id)]
                                self.tree.insert(tree_experiment, bisect(contents, simulation.id), iid=simulation.id, text=text, values=(" "))           
                            # If the input file is already in the treeview, make sure it has the correct columns
                            if self.tree.exists(simulation.id):
                                i = experiment.simulations.index(simulation)
                                text = experiment.variedValues[i].replace("\\","").replace("$","").replace(",", " ")
                                self.tree.item(simulation.id, text=text, values=(" "))

            # Remove focus from items
            for item in self.tree.selection():
                self.tree.selection_remove(item)

