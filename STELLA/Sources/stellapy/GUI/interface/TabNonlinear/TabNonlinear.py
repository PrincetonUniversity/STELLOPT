
# TODO: Change "impurity x" to something that recognizes the impurity

#################################################################
#                   CLASS FOR THE FIRST TAB
#################################################################
''' 
'''

# Load modules
import tkinter as tk
from tkinter import ttk 
from bisect import bisect

# For the figures
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Personal modules
from stellapy.utils import initiate_nesteddict
from stellapy.GUI.graph_tools import Progress
from stellapy.GUI.graph_tools import CanvasForGraphs, PoppedOutWindow
from stellapy.GUI.graph_tools import PAD_TITLE2, PAD_LABEL2, PAD_ENTRY2 #@unusedimport 

# Plots
from stellapy.plot import plot_fluxesVsTime
from stellapy.plot import plot_saturatedfluxVsRho
from stellapy.plot import plot_potentialVsTime

#################################################################
#                   CLASS FOR THE FIRST TAB
#################################################################
class TabNonlinear:
 
#################################################################
#                          WIDGETS
#################################################################

    # Initate the tab for "Select simulations"
    def __init__(self, tab5): 

        #======================================================
        # Safe info from the tab and root that can be passed on
        #======================================================

        self.window = tab5                      # Direct accesses to the window of tab "Select Simulations" 
        self.root = tab5.root                   # Needed to center windows based on the root screen
        self.dict_awthemes = tab5.root.awthemes # To make the tk widgets look like the ttk widgets

        #===========
        # VARIABLES
        #===========

        # For the widgets
        self.replot = False
        self.plotted_inputFiles = None
        self.initiated_canvas = None
        self.axis_id = 0
        self.tree_experiments = {}
        self.frame_species = {} 
        self.tree = {} 
        
        # Data which needs to be transferred between plotting functions
        self.t_range = None
        self.sat_flux = None
        
        # Identifiers of the data that needs to be plotted
        self.experiment = None
        self.simulation = None
        self.experiment_id = "All experiments"
        self.simulation_id = "All simulations"
        self.knob1 = "vmec_parameters"
        self.key1 = "rho"
        self.knob2 = "-"
        self.key2 = "-"
        self.show_vspan = True
        self.normalizeBySaturation = False
        
        #===========
        # FRAMES
        #===========

        # Create the frames and add them to the tab3
        self.frame_graph1   = ttk.LabelFrame(self.window, text="   Fluxes versus time  ", **self.dict_awthemes['labelframe2'])
        self.frame_graph2   = ttk.LabelFrame(self.window, text="   Satured flux versus parameter  ", **self.dict_awthemes['labelframe2'])
        self.frame_graph3   = ttk.LabelFrame(self.window, text="   Potential squared versus time  ", **self.dict_awthemes['labelframe2'])
        self.frame_options  = ttk.Frame(self.window) 

        # Configure the frames     
        tk.Grid.rowconfigure(   self.window, 0, weight=4,  uniform="yes") 
        tk.Grid.rowconfigure(   self.window, 1, weight=10, uniform="yes") 
        tk.Grid.columnconfigure(self.window, 0, weight=1,  uniform="for sure") 
        tk.Grid.columnconfigure(self.window, 1, weight=1,  uniform="for sure") 
 
        # Add the frames to the main window
        self.frame_options.grid(in_=self.window, row=0, column=0, padx=(20,20), pady=(5,10),stick='NSEW', columnspan=2)
        self.frame_graph1.grid( in_=self.window, row=1, column=0, padx=(20,10), pady=(10,20),stick='NSEW')
        self.frame_graph2.grid( in_=self.window, row=1, column=1, padx=(10,20), pady=(10,20),stick='NSEW')
        self.frame_graph3.grid( in_=self.window, row=1, column=0, padx=(20,10), pady=(10,20),stick='NSEW')
        self.frame_graph3.grid_remove()
        
        # Initiate the widgets in the frame with the options
        self.initiate_frameOptions()
        self.root.update_idletasks()
        
        # Bind focus on every click on a widget: this makes sure the entry widgets become unfocused
        self.frame_graph1.bind("<1>", lambda event: self.frame_graph1.focus_set())
        self.frame_graph2.bind("<1>", lambda event: self.frame_graph2.focus_set())
        self.frame_graph3.bind("<1>", lambda event: self.frame_graph2.focus_set())
        self.frame_options.bind("<1>", lambda event: self.frame_options.focus_set()) 
        
        #=============
        # LOAD FIGURE
        #=============
        
        # When the tab is visible, make sure the correct figure is loaded
        self.frame_options.bind("<Visibility>", self.load_figure)
    
        # Prevent folding of header comment on next line
        if True: return
    
          
#=======================================
# Initiate the bar with extra options
#=======================================

    def initiate_frameOptions(self):
    
        # Create the sub_labelframes in options plus a button to apply any changes
        self.subframe_data  = ttk.LabelFrame(self.frame_options, text="   Simulations   ", **self.dict_awthemes['labelframe2'])
        self.subframe_plots = ttk.LabelFrame(self.frame_options, text="   Options   ",       **self.dict_awthemes['labelframe2'])
        self.subframe_time  = ttk.Frame(self.frame_options)
        self.frame_progress = ttk.Frame(self.frame_options)

        # Add the progress object
        self.class_progress = Progress(self, anchor="fill", length=300)
        self.class_progress.move(0,"Please select simulations.")       
        
        # Configure the frame
        tk.Grid.rowconfigure(   self.frame_options, 0, weight=1) 
        tk.Grid.rowconfigure(   self.frame_options, 1, weight=0)      
        tk.Grid.columnconfigure(self.frame_options, 0, weight=0)   
        tk.Grid.columnconfigure(self.frame_options, 1, weight=0)   
        tk.Grid.columnconfigure(self.frame_options, 2, weight=1) 
        
        # Add the frames to <frame_options>
        self.subframe_data.grid( row=0, column=0, padx=(0,5), pady=(20,0), stick='NSEW')
        self.subframe_plots.grid(row=0, column=1, padx=(5,5), pady=(20,0), stick='NSEW')
        self.subframe_time.grid( row=0, column=2, padx=(5,0), pady=(0,0), stick='NSEW', rowspan=2)
        self.frame_progress.grid(row=1, column=0, padx=(0,5), pady=(5,0), stick='NSEW', columnspan=2)

        # Fill the subframes with widgets
        self.initiate_frameOptions_plotsFrame()
        self.initiate_frameOptions_dataFrame()
        self.initiate_frameOptions_timeFrame()
        
    #------------------------------------------- 
    def initiate_frameOptions_dataFrame(self):
        
        # Initiate the attributes
        self.options_experiment = ["First plot data"]
        self.options_simulation = ["First plot data"]
        self.options_simulationsids = []
        self.options_fluxes  = ["Heat flux", "Momentum flux", "Particle flux", "Potential squared"]     # Choices of the flux to be plotted
        self.options_species = ["Ions", "Electrons", "Impurity 1", "Impurity 2", "All species"]   # Choices of the flux to be plotted
        self.options_par1 = ["rho", "tiprim", "teprim", "fprim", "delta t", "nmu", "nvgrid"]
        self.options_par2 = ["-", "rho", "tiprim", "teprim", "fprim", "delta t", "nmu", "nvgrid"]
        
        # Choose between the experiments
        self.var_experiment = tk.StringVar(value=self.options_experiment[0]); width=20
        self.lbl_experiment = ttk.Label(self.subframe_data, text="Experiment: ")
        self.mnu_experiment = ttk.OptionMenu(self.subframe_data, self.var_experiment, self.options_experiment[0], *self.options_experiment, style='option.TMenubutton')
        self.mnu_experiment["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
        self.mnu_experiment.config(width=width)
        self.var_experiment.trace('w', self.change_plottedExperiment) 
        
        # Choice between the simulations
        self.var_simulation = tk.StringVar(value=self.options_simulation[0])
        self.lbl_simulation = ttk.Label(self.subframe_data, text="Simulation: ")
        self.mnu_simulation = ttk.OptionMenu(self.subframe_data, self.var_simulation, self.options_simulation[0], *self.options_simulation, style='option.TMenubutton')
        self.mnu_simulation["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
        self.mnu_simulation.config(width=width)
        self.var_simulation.trace('w', self.change_plottedSimulation)
        
        # Choose the parameter(s) for the "Saturated flux versus parameter" plot
        self.var_par1 = tk.StringVar(value=self.options_par1[0]); width=10
        self.var_par2 = tk.StringVar(value=self.options_par2[0])
        self.lbl_par1 = ttk.Label(self.subframe_data, text="Parameter 1:    ")
        self.lbl_par2 = ttk.Label(self.subframe_data, text="Parameter 2:    ")
        self.mnu_par1 = ttk.OptionMenu(self.subframe_data, self.var_par1, self.options_par1[0], *self.options_par1, style='option.TMenubutton')
        self.mnu_par2 = ttk.OptionMenu(self.subframe_data, self.var_par2, self.options_par2[0], *self.options_par2, style='option.TMenubutton')
        self.mnu_par1["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
        self.mnu_par2["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
        self.mnu_par1.config(width=width)
        self.mnu_par2.config(width=width)
        def update_parameter1(*args): self.update_parameter(i=1); return 
        def update_parameter2(*args): self.update_parameter(i=2); return
        self.var_par1.trace('w', update_parameter1) 
        self.var_par2.trace('w', update_parameter2) 
        
        # Choose the fluxes: heat, momentum, particle or the potential
        def update_yQuantity(*args): 
            if self.var_norm.get()==0: self.normalizeBySaturation = False
            if self.var_norm.get()==1: self.normalizeBySaturation = True
            self.replot=True; self.plot_graph(1, None, "right")
            self.replot=True; self.plot_graph(0, None, "left")
            self.update_treeViewSaturatedFluxes()
        self.var_flux = tk.StringVar(value=self.options_fluxes[0])
        self.mnu_flux = ttk.OptionMenu(self.subframe_data, self.var_flux, self.options_fluxes[0], *self.options_fluxes, style='option.TMenubutton')
        self.mnu_flux["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
        self.mnu_flux.config(width=30)
        self.lbl_flux = ttk.Label(self.subframe_data, text="Y-quantity: ")
        self.var_flux.trace('w', update_yQuantity) 
        
        # Choose the species: ions, electrons, impurity 1, impurity 2, ...
        self.var_specie = tk.StringVar(value=self.options_species[0])
        self.mnu_specie = ttk.OptionMenu(self.subframe_data, self.var_specie, self.options_species[0], *self.options_species, style='option.TMenubutton')
        self.mnu_specie["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
        self.mnu_specie.config(width=width)
        self.lbl_specie = ttk.Label(self.subframe_data, text="Specie: ")
        
        # Configure the frame
        tk.Grid.rowconfigure(self.subframe_data, 0, weight=0)  
        tk.Grid.rowconfigure(self.subframe_data, 1, weight=0)  
        tk.Grid.rowconfigure(self.subframe_data, 2, weight=0)  
        tk.Grid.rowconfigure(self.subframe_data, 3, weight=0)  
        tk.Grid.rowconfigure(self.subframe_data, 4, weight=0)  
        tk.Grid.rowconfigure(self.subframe_data, 5, weight=0)  
        tk.Grid.columnconfigure(self.subframe_data, 0, weight=1) 
        
        # Place the widgets in the frame
        self.lbl_experiment.grid(row=0, column=0, stick='NSEW', padx=(0,0), pady=(2,2))
        self.mnu_experiment.grid(row=0, column=1, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=2, ipady=1)
        self.lbl_simulation.grid(row=1, column=0, stick='NSEW', padx=(0,0), pady=(2,2))
        self.mnu_simulation.grid(row=1, column=1, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=2, ipady=0)
        self.lbl_flux.grid(      row=2, column=0, stick='NSEW', padx=(0,0), pady=(2,2))
        self.mnu_flux.grid(      row=2, column=1, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=2, ipady=2)
        self.lbl_specie.grid(    row=3, column=0, stick='NSEW', padx=(0,0), pady=(2,2))
        self.mnu_specie.grid(    row=3, column=1, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=2, ipady=2)
        self.lbl_par1.grid(      row=4, column=0, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=1, ipady=0)
        self.mnu_par1.grid(      row=4, column=1, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=1, ipady=0)
        
    #------------------------------------------      
    def initiate_frameOptions_plotsFrame(self):
        
        # Initiate the attributes
        self.species = [0]
        self.y_quantity = "q_flux"
        
        # Set the default value of the plot out of range so we don't start with a plot
        self.var_plot  = tk.IntVar(value=0)
        
        # Add the options
        def update_vspan(): 
            if self.var_vspan.get()==0: self.show_vspan = False
            if self.var_vspan.get()==1: self.show_vspan = True
            self.replot=True; self.plot_graph(1, None, "right")
            self.replot=True; self.plot_graph(0, None, "left")
        self.var_vspan = tk.IntVar(value=1)
        self.chk_vspan = ttk.Checkbutton(self.subframe_plots, text=" Show time frames", variable=self.var_vspan, command=update_vspan)
        def update_norm(): 
            if self.var_norm.get()==0: self.normalizeBySaturation = False
            if self.var_norm.get()==1: self.normalizeBySaturation = True
            if self.var_flux.get() in ["Heat flux", "Momentum flux", "Particle flux"]:
                self.replot=True; self.plot_graph(1, None, "right")
                self.replot=True; self.plot_graph(0, None, "left")
            if self.var_flux.get() in ["Potential squared"]:
                self.replot=True; self.plot_graph(2, None, "left")
            self.update_treeViewSaturatedFluxes()
        self.var_norm  = tk.IntVar(value=0)
        self.chk_norm  = ttk.Checkbutton(self.subframe_plots, text=" Normalize by saturation", variable=self.var_norm, command=update_norm)
        
        # Configure the frame
        tk.Grid.rowconfigure(self.subframe_plots, 0, weight=0) 
        tk.Grid.rowconfigure(self.subframe_plots, 1, weight=0) 
        tk.Grid.columnconfigure(self.subframe_plots, 0, weight=0) 
        
        # Add the options to the frame
        self.chk_vspan.grid(row=0, column=0, stick='NSEW', **PAD_LABEL2)
        self.chk_norm.grid( row=1, column=0, stick='NSEW', **PAD_LABEL2)
        
    #-------------------------------------------       
    def initiate_frameOptions_timeFrame(self):
        ''' Choose the time frame to measure the saturated fluxes and
        display some information about the saturation of the fluxes. '''
        
        # Tabbed display for species 1, 2, 3, ....
        self.tab_species = ttk.Notebook(self.subframe_time, style='species.TNotebook')
        self.initiate_speciesTab(0,"  Ions  ")
        tk.Grid.rowconfigure(self.subframe_time, 0, weight=1) 
        tk.Grid.columnconfigure(self.subframe_time, 0, weight=1)  
        tk.Grid.rowconfigure(self.tab_species, 0, weight=1) 
        tk.Grid.columnconfigure(self.tab_species, 0, weight=1)    
        self.tab_species.grid(row=0, column=0, sticky="NSEW")
        return 
    
    #-------------------------------------------    
    def initiate_speciesTab(self, specie_id, title):
        
        # For each specie we make a tabbed view
        self.frame_species[specie_id] = ttk.Frame(self.tab_species, padding=(10,10,10,10))
        self.tab_species.add(self.frame_species[specie_id], text=title) 
        self.tree_experiments[specie_id] = []
        
        # The list of time frames will be displayed inside a treeview widget
        self.tree[specie_id] = ttk.Treeview(self.frame_species[specie_id], show="tree",columns=("#0","#1","#2"))
        self.tree[specie_id].grid(row=0, column=0, padx=0, pady=0, sticky='NESW')
        tk.Grid.rowconfigure(self.frame_species[specie_id], 0, weight=1)
        tk.Grid.columnconfigure(self.frame_species[specie_id], 0, weight=1)
            
        # Create a treeview widget and bind keypresses to the widget
        self.tree[specie_id].bind("<Double-1>", self.onDoubleClick)
        
        # Add columns to the tree view 
        self.tree[specie_id].heading("#0", text="Simulation",     anchor=tk.W) 
        self.tree[specie_id].heading("#1", text="Time frame",     anchor=tk.W) 
        self.tree[specie_id].heading("#2", text="Saturated flux", anchor=tk.W) 
        self.tree[specie_id]["displaycolumns"] = ("#0", "#1", "#2")

        # When the tab is visible, make sure the columns have a good size
        def resize_columns(*args): 
            self.tree[specie_id].column("#0", width=int(self.tree[specie_id].winfo_width()*4/10), stretch=tk.YES)
            self.tree[specie_id].column("#1", width=int(self.tree[specie_id].winfo_width()*3/10), stretch=tk.YES) 
            self.tree[specie_id].column("#2", width=int(self.tree[specie_id].winfo_width()*3/10), stretch=tk.YES)   
        self.frame_species[specie_id].bind("<Visibility>", resize_columns)
        if True: return


#=======================================
# Initiate site bar with extra options
#=======================================

    def initiate_canvas(self):
        ''' The toolbar requires the class to have a reset_graph and popout_window function.
        The optionswindow requires the Graph objects to have attributes: ax, x_name, x_key, y_name, y_key, range, label
        '''
        
        # Create a list for the graph classes, create the figure for the canvas and initiate the canvas.
        self.initiated_canvas = True
        self.Graph = [None, None, None]      
        self.Canvas = [None, None, None]      
        self.figure1 = plt.figure("nonlinear1")  
        self.figure2 = plt.figure("nonlinear2")  
        self.figure3 = plt.figure("nonlinear3")    
        self.figure1.set_tight_layout(False) 
        self.figure2.set_tight_layout(False)
        self.figure3.set_tight_layout(False)
        
        # Put the canvasses on the screen: left canvas has axis_id=0; right canvas has axis_id=1
        CanvasForGraphs(self.root, self.frame_graph1, self.root.tab_Nonlinear, self.figure1, axis_id=0)
        CanvasForGraphs(self.root, self.frame_graph2, self.root.tab_Nonlinear, self.figure2, axis_id=1)
        CanvasForGraphs(self.root, self.frame_graph3, self.root.tab_Nonlinear, self.figure3, axis_id=2)

        # Initiate the plotting class for the main window
        self.initiate_plottingClass()

    #-------------------------------------------
    def initiate_plottingClass(self, figure=None):
        ''' Initiate the plotting object to hold all the plotting attributes.
        This function is seperate from initiate_canvas since it is called by poppedout_window as well. '''
        
        # Change the canvas color
        self.figure1.patch.set_facecolor(self.root.color['canvas']) 
        self.figure2.patch.set_facecolor(self.root.color['canvas']) 
        self.figure3.patch.set_facecolor(self.root.color['canvas'])
        
        # Set the color of the axes, ticks and background
        plt.rcParams['text.color']       = self.root.color['fg']
        plt.rcParams['axes.edgecolor']   = self.root.color['fg']
        plt.rcParams['axes.labelcolor']  = self.root.color['fg']
        plt.rcParams['xtick.color']      = self.root.color['fg']
        plt.rcParams['ytick.color']      = self.root.color['fg']
        plt.rcParams['axes.facecolor']   = self.root.color['canvas']
        plt.rcParams['figure.facecolor'] = self.root.color['canvas']
        
        # Initiate the plotting class for figures on the GUI
        if figure is None:
            # Make the axis instance and the required attributes for the options window through the class <graph>
            self.Graph[0] = graph(self.figure1, self.var_plot, 0)
            self.Graph[1] = graph(self.figure2, self.var_plot, 1)
            self.Graph[2] = graph(self.figure3, self.var_plot, 2)
            # Put the empty figure on the GUI
            self.Canvas[0].draw_idle(); self.root.update_idletasks()
            self.Canvas[1].draw_idle(); self.root.update_idletasks()
            self.Canvas[2].draw_idle(); self.root.update_idletasks()
        
        # Initiate the plotting class for a popped out window
        if figure is not None:
            self.root.graph_poppedOut.append(graph(figure, self.var_plot, 0))
            
        # Prevent indentention of header comment on next lines
        if True: return 
        
#################################################################
#                          METHODS
#################################################################

#====================
# For the tree view
#====================

    def update_speciesFrame(self):

        # Find the maximum number of species
        dim_species = 0
        for experiment in self.root.Research.experiments:
            dim_species = max(dim_species, experiment.simulations[0].dim_species)
            
        # If there are more species, add more tabs 
        frame_ids = list(self.frame_species.keys()) 
        
        # When all simulations are removed, remove the additional species tabs
        try: unique_folders = self.root.Research.unique_folders
        except: unique_folders = []
        if len(unique_folders)==0: dim_species=1      
        
        # Add more tabs if there aren't enough
        if len(frame_ids) < dim_species: 
            for i in range(len(frame_ids), dim_species):
                self.initiate_speciesTab(i, "Species "+str(i+1))
                
        # Remove tabs if there are too many
        if dim_species < len(frame_ids): 
            for i in range(dim_species, len(frame_ids)):
                self.tab_species.forget(self.frame_species["tab "+str(i+1)])
                self.frame_species["tab "+str(i+1)].destroy()                                       
        
        # Change the name if the tabs or named wrongly
        if not dim_species==1:
            if dim_species == 2: # Assume the second species is always kinetic electrons 
                self.tab_species.tab(self.frame_species[1], text = 'Kinetic electrons')
            elif dim_species== 3: # Assume second species is always kinetic electrons and third impurities
                self.tab_species.tab(self.frame_species[1], text = 'Kinetic electrons')
                self.tab_species.tab(self.frame_species[2], text = 'Impurities')
        return 
    
    #---------------------------------------              
    def update_treeViewSaturatedFluxes(self):
        ''' Update the time frames to calculate the saturated flux. '''

        # Only update when there are input files
        if len(self.root.input_files) != 0:
            
            # Remove all the items from the tree view
            for specie_id in range(len(self.tree_experiments)):
                self.tree[specie_id].delete(*self.tree[specie_id].get_children())
                self.tree_experiments[specie_id] = []
    
            # Iterate over the experiments
            for experiment in self.root.Research.experiments:
                
                # Iterate over the species to fill each tab
                for specie_id in range(experiment.simulations[0].dim_species):
    
                    # Check whether the experiment is already in the treeview
                    if not self.tree[specie_id].exists(experiment.id):
                        contents = [self.tree[specie_id].item(child)["text"] for child in self.tree[specie_id].get_children("")]
                        index = bisect(contents, experiment.id)
                        text = experiment.line_label
                        self.tree_experiments[specie_id].append(self.tree[specie_id].insert(parent="", index=index, iid=experiment.id, text=text, values=("", "")))
                        self.tree[specie_id].item(experiment.id, open=True)
            
                    # Add the values to the experiment: marker label; time frame, Q_sat
                    for tree_experiment in self.tree_experiments[specie_id]:
                        self.tree[specie_id].selection_set(tree_experiment)  
                        experiment_iid = self.tree[specie_id].selection()[0]
                        if experiment_iid==experiment.id:
                            for simulation in experiment.simulations:
                                m = simulation.marker_label.replace("$", "").replace("\,", " ")
                                try: q = "Q = " + "{:.2e}".format(self.sat_flux[experiment.id][simulation.id][specie_id])
                                except: q = "Q = ?"
                                try: 
                                    t = self.t_range[experiment.id][simulation.id][specie_id]
                                    if self.normalizeBySaturation:
                                        time_peak = self.t_range[experiment.id][simulation.id]["time peak"]
                                        t = [round(t[0]/time_peak,2), round(t[1]/time_peak,2)]   
                                    if not self.normalizeBySaturation:
                                        t = [int(t[0]), int(t[1])]
                                    t = "t = [" + str(t[0]) + ", " + str(t[1]) + "]"
                                except: t = "[ - ]"
                                contents = [self.tree[specie_id].item(child)["text"] for child in self.tree[specie_id].get_children(experiment.id)]
                                if self.tree[specie_id].exists(simulation.id):
                                    self.tree[specie_id].item(simulation.id, text=m, values=(t, q)) 
                                if not self.tree[specie_id].exists(simulation.id):
                                    self.tree[specie_id].insert(tree_experiment, bisect(contents, m), iid=simulation.id, text=m, values=(t, q))           
                                

        # Remove focus from items
        for item in self.tree[specie_id].selection():
            self.tree[specie_id].selection_remove(item)
        return 
        
    #----------------------------------
    def update_treeViewAfterChange(self, values, simulation_iid, experiment_iid, specie_id):
        ''' After manually changing the time frames, update the figure. '''
        
        # Identify the simulation
        simulation = None; experiment = None
        for experiment_temp in self.root.Research.experiments:
            for simulation_temp in experiment_temp.simulations: 
                if simulation_iid == simulation_temp.id:
                    experiment = experiment_temp
                    simulation = simulation_temp
                    
        # Add the normalization factor
        norm = self.t_range[experiment.id][simulation.id]["time peak"] if self.normalizeBySaturation else 1
                        
        # Update the t_range object
        try: 
            time_range = [float(values[0].split("[")[-1].split(",")[0]),float(values[0].split(",")[-1].split("]")[0])]
            time_range = [time_range[0]*norm, time_range[1]*norm]
            self.t_range[experiment_iid][simulation_iid][specie_id] = time_range
            time_range = "["+str(time_range[0]*norm)+", "+str(time_range[1]*norm)+"]"
            simulation.simulation_file['TIME FRAME'][str(specie_id)] = time_range
            simulation.simulation_file.write(open(simulation.simulation_path, 'w')) 
        except: pass
        
        # Replot the figure
        self.replot = True; self.plot_graph(1, None, "right")
        self.replot = True; self.plot_graph(0, None, "left")
        return 
    
    #----------------------------------
    def onDoubleClick(self, event):
        ''' Executed, when a row is double-clicked. Opens EntryPopup above the item's column, 
        so it is possible to change the text. '''
    
        # Close previous popup if there is one
        try: self.entryPopup.on_return()
        except: pass
    
        # Get the treeview that was selected
        specie_id = self.tab_species.index(self.tab_species.select())
    
        # Select row and column that was clicked on
        rowid = self.tree[specie_id].identify_row(event.y)
        columnid = self.tree[specie_id].identify_column(event.x)
    
        # Get column position info
        x,y,width,height = self.tree[specie_id].bbox(rowid, columnid)
    
        # Place Entry popup properly         
        text = self.tree[specie_id].item(rowid, 'text')
        values = self.tree[specie_id].item(rowid, 'values')
        self.entryPopup = EntryPopup(self.root, self.tree[specie_id], specie_id, rowid, columnid, text, values, self.update_treeViewAfterChange)
        self.entryPopup.place(x=x, y=y+height//2, anchor=tk.W)
        if True: return 
     
#====================
# Functions
#====================

    def load_figure(self, *args):
        '''  When the tab is visible, make sure the correct figure is loaded. 
        If the simulations were changed, replot the graph. '''        
        # If the canvas isn't loaded, load them
        if self.initiated_canvas==None:
            self.initiate_canvas()
        
        # Load figure
        plt.figure("nonlinear")
        
        # Try to load the saved t_ranges in the simulation.ini files
        self.t_range = initiate_nesteddict()
        for experiment in self.root.Research.experiments:
            for simulation in experiment.simulations:
                if "TIME FRAME" in simulation.simulation_file:
                    for specie in simulation.simulation_file['TIME FRAME'].keys():
                        time = simulation.simulation_file['TIME FRAME'][specie]
                        if specie != "time peak": 
                            time = time.split("[")[-1].split("]")[0].split(",")
                            time = [ float(t) for t in time ]
                            self.t_range[experiment.id][simulation.id][int(specie)] = time 
                        if specie == "time peak": 
                            self.t_range[experiment.id][simulation.id][specie] = float(time)
        
        # Connect the Progress class to the research object
        self.root.Research.Progress = self.class_progress
        for experiment in self.root.Research.experiments:
            experiment.Progress = self.class_progress
            for simulation in experiment.simulations:
                simulation.Progress = self.class_progress
                
        # For the saturated flux versus paramater plot
        if len(self.root.Research.experiments)!=0:
            if len(experiment.variedVariables)!=0:
                experiment = self.root.Research.experiments[0]
                variable = experiment.variedVariables[0]
                self.knob1 = variable.split(":")[0]
                self.key1 = variable.split(": ")[-1]
                self.var_par1.set(self.key1)
                
        # Make sure the simulations/experiments/species are updated
        if self.root.input_files != []:
            self.display_plottedModesAndExperiments()
            self.update_speciesFrame()
        return
    
    #-------------------------      
    def show_progress(self, status="start", poppedout_id=None):
        # While plotting show a progress bar turn the cursor into a waiting cursor.
        if status=="start" and poppedout_id == None: 
            self.class_progress.move(0,"Start plotting.")
            self.root.config(cursor="watch")
            self.root.update_idletasks() 
        # When finished show the canvas and change the cursor back.
        if status=="finished" and poppedout_id == None: 
            self.class_progress.move(100, "No tasks are running.")
            self.root.config(cursor="")
        # When there are no simulations revert to the loading bar
        if status=="nothing" and poppedout_id == None: 
            self.class_progress.move(0,"Please select simulations.")
        return
    
    #----------------------
    def popout_window(self, axis_id):
        '''Replot the figure in a seperate window when the button on the toolbar is clicked.'''
        
        # Create a popped out window
        poppedout_window = PoppedOutWindow(self.root)
        poppedout_id = poppedout_window.poppedout_id
        
        # Add a fitting title
        if self.var_plot.get() == 1: poppedout_window.set_title("Stellapy: Fluxes versus time")
        if self.var_plot.get() == 2: poppedout_window.set_title("Stellapy: Fluxes versus (kx,ky)")

        # Initiate the plotting class
        self.initiate_plottingClass(figure=poppedout_window.figure)
        
        # Initiate the canvas        
        CanvasForGraphs(self.root, poppedout_window.frame, self.root.tab_Nonlinear, poppedout_window.figure, poppedout_id=poppedout_id)
        
        # Parse the data from the main canvas to the popped out canvas
        self.root.graph_poppedOut[poppedout_id].range = self.Graph[0].range
        self.root.graph_poppedOut[poppedout_id].label = self.Graph[0].label
        
        # Now plot the figure on the popped out canvas
        if axis_id==0: self.replot = True; self.plot_graph(None, poppedout_id, "left")
        if axis_id==1: self.replot = True; self.plot_graph(None, poppedout_id, "right")
        
        # Update screen
        self.root.canvasPoppedOut[poppedout_id].draw_idle()
        return 
    
    #----------------------------
    def change_plot(self, *args):
        
        # Set the title of the labelframe of the graph
        if self.var_plot.get() in [1]: 
            self.frame_graph.config(text="   Fluxes versus time  ")
        if self.var_plot.get() in [2]: 
            self.frame_graph.config(text="   Fluxes versus (kx,ky)  ")    
        self.replot = True; self.plot_graph(0, None, "left")
        self.replot = True; self.plot_graph(1, None, "right")
        return 
    
    #---------------------------------------
    def change_plottedExperiment(self, *args):
        
        # Change the experiment that is plotted
        self.experiment_id = self.var_experiment.get()
    
        # Get a reference to the experiment
        self.experiment = None
        for experiment in self.root.Research.experiments:
            if experiment.id == self.experiment_id:
                self.experiment = experiment
        if self.experiment==None and self.experiment_id!="All experiments":
            self.experiment = self.root.Research.experiments[0]
            self.experiment_id = self.experiment.id
        if self.experiment==None and self.experiment_id=="All experiments":
            self.experiment = self.root.Research.experiments[0]
            
        # Replot the graph
        self.Graph[0].load_defaults()
        self.Graph[1].load_defaults()
        return
        
    #------------------------------------
    def change_plottedSimulation(self, *args):
        
        # Change the simulation that is plotted
        self.simulation_id = self.var_simulation.get()
        self.simulation_id = self.options_simulationsids[self.options_simulations.index(self.simulation_id)]
        # Get a reference to the simulation
        self.simulation = None
        for simulation in self.experiment.simulations:
            if simulation.id == self.simulation_id:
                self.simulation = simulation
        if self.simulation==None and self.simulation_id!="All simulations":
            self.simulation = self.experiment.simulations[0]
            self.simulation_id = self.simulation.id
        if self.simulation==None and self.simulation_id=="All simulations":
            self.simulation = self.experiment.simulations[0] 
        self.replot=True
        return 
    
    #--------------------------
    def update_parameter(self, i):
        
        # Load the defaults when we touch this
        self.Graph[1].load_defaults()
        
        # Get the parameter
        if i==1: var_par = self.var_par1.get()
        if i==2: var_par = self.var_par2.get()
        
        # Assign the correct parameter
        if var_par=="rho":      knob = "vmec_parameters";        key="rho"
        if var_par=="tprim":    knob = "species_parameters_1";   key="tprim"
        if var_par=="tiprim":   knob = "species_parameters_1";   key="tprim"
        if var_par=="teprim":   knob = "species_parameters_2";   key="tprim"
        if var_par=="fprim":    knob = "species_parameters_1";   key="fprim"
        if var_par=="delta t":  knob = "knobs";                  key="delt"
        if var_par=="delt":     knob = "knobs";                  key="delt"
        if var_par=="nmu":      knob = "vpamu_grids_parameters"; key="nmu"
        if var_par=="nvgrid":   knob = "vpamu_grids_parameters"; key="nvgrid"
        if var_par=="-":        knob = "-";                      key="-"
        
        # Save the knobs and keys
        if i==1: self.knob1 = knob;     self.key1 = key
        if i==2: self.knob2 = knob;     self.key2 = key
        
        # Count the parameters
        if self.var_par2.get() == "-": self.dim_parameters=1
        if self.var_par2.get() != "-": self.dim_parameters=2
        print("CHANGED PARAMETER: ", self.knob1, self.key1)
        if True: return
        
#===================
# Update the graph 
#===================

    def reset_graph(self, axis_id, poppedout_id=None):
        ''' Method called by the "Reset graph" button on the canvas. '''
        
        # Perhaps the button was clicked on the popped out window
        if poppedout_id == None: Graph = self.Graph[0]
        if poppedout_id != None: Graph = self.root.graph_poppedOut[poppedout_id]
        
        # Reset the axis
        Graph.range["x"] = None
        Graph.range["y"] = None

        # Now plot the graph 
        if axis_id==0: self.replot = True; self.plot_graph(1, poppedout_id, "right")
        if axis_id==0: self.replot = True; self.plot_graph(0, poppedout_id, "left")
        if axis_id==1: self.replot = True; self.plot_graph(1, poppedout_id, "right")
        return
    
    #-----------------------------------------   
    def plot_graph(self, axis_id, poppedout_id=None, plot="left"):
        
        # Get the data: specie, flux, parameter, ...
        if self.var_specie.get()=="All species": 
            dim_species = self.root.Research.experiments[0].simulations[0].dim_species
            self.species = [ i for i in range(dim_species)]
        if self.var_specie.get()!="All species":
            self.species = [int(self.options_species.index(self.var_specie.get()))]
        if self.var_flux.get() == "Heat flux":     self.y_quantity = "q_flux"
        if self.var_flux.get() == "Momentum flux": self.y_quantity = "v_flux"
        if self.var_flux.get() == "Particle flux": self.y_quantity = "p_flux" 
        
        # Perhaps the button was clicked on the popped out window
        if poppedout_id == None: Graph = self.Graph[axis_id]
        if poppedout_id != None: Graph = self.root.graph_poppedOut[poppedout_id]
        
        # For now load the defaults
        Graph.load_defaults()

        # Only proceed if there are input_files
        if self.root.input_files != []:
            
            # Get the input_files that need to be plotted
            input_files = self.root.Research.input_files

            # If there are files plot them, if they changed plot them, if its a popout plot them
            if (input_files != [] and self.plotted_inputFiles != input_files) or poppedout_id!=None or self.replot==True:
                
                # Clear the axis and remember which input_files are plotted
                Graph.ax.clear()
                self.plotted_inputFiles = input_files.copy() 
                self.replot = False
                
                # While plotting show a progress bar turn the cursor into a waiting cursor
                self.show_progress("start", poppedout_id) 
                
                # Plot the heatflux versus time
                if plot=='left' and self.var_flux.get() in ["Heat flux", "Momentum flux", "Particle flux"]:
                    self.frame_graph1.grid()
                    self.frame_graph3.grid_remove()
                    plot_fluxesVsTime(\
                            # Specify which simulations to plot
                                research=self.root.Research,\
                                experiment_id=self.experiment_id,\
                                simulation_id=self.simulation_id,\
                                sat_flux=self.sat_flux,\
                            # Specify data range
                                species=self.species,\
                                y_quantity=self.y_quantity,\
                                x_range=Graph.range["x"],\
                                y_range=Graph.range["y"],\
                                t_range=self.t_range,\
                            # Labels
                                x_label=Graph.label["x"],\
                                y_label=Graph.label["y"],\
                                title=Graph.label["title"],\
                            # For the GUI the figure object already exists 
                                show_figure = False,\
                                ax=Graph.ax, \
                                Progress=self.class_progress,\
                            # Toggles
                                show_vspan=self.show_vspan,\
                                normalize=self.normalizeBySaturation,\
                                log=False,\
                                units="normalized")
                    
                # Plot the heatflux versus time
                if plot=='left' and self.var_flux.get() in ["Potential squared"]:
                    self.frame_graph3.grid()
                    self.frame_graph1.grid_remove()
                    plot_potentialVsTime(\
                            # Specify which simulations to plot
                                research=self.root.Research,\
                                experiment_id=self.experiment_id,\
                                simulation_id=self.simulation_id,\
                            # Specify data range
                                x_range=Graph.range["x"],\
                                y_range=Graph.range["y"],\
                            # Labels
                                x_label=Graph.label["x"],\
                                y_label=Graph.label["y"],\
                                title=Graph.label["title"],\
                            # For the GUI the figure object already exists 
                                show_figure = False,\
                                ax=Graph.ax,\
                                Progress=self.class_progress,\
                            # Toggles
                                log=False,\
                                units="normalized")
                        
                # Plot the figure of the saturated heat flux before plotting heatflux versus time!
                if plot=='right':
                    self.sat_flux, self.t_range = \
                        plot_saturatedfluxVsRho(\
                            # Specify which simulations to plot
                                research=self.root.Research,\
                                experiment_id=self.experiment_id,\
                                simulation_id=self.simulation_id,\
                                parameter_knob=self.knob1,\
                                parameter_key=self.key1,\
                            # Specify data range
                                species=self.species,\
                                y_quantity=self.y_quantity,\
                                x_range=Graph.range["x"],\
                                y_range=Graph.range["y"],\
                                t_range=self.t_range,\
                            # Labels
                                x_label=Graph.label["x"],\
                                y_label=Graph.label["y"],\
                                title=Graph.label["title"],\
                            # For the GUI the figure object already exists 
                                show_figure = False,\
                                ax=Graph.ax, \
                                Progress=self.class_progress,\
                            # Toggles
                                log=False,\
                                units="normalized")

                    # Save the t_range to the simulation.ini files
                    for experiment in self.root.Research.experiments:
                        for simulation in experiment.simulations: 
                            if not 'TIME FRAME' in simulation.simulation_file:
                                simulation.simulation_file['TIME FRAME'] = {}
                            for specie in self.t_range[experiment.id][simulation.id].keys():
                                time = self.t_range[experiment.id][simulation.id][specie]
                                if specie != "time peak": time = "["+str(time[0])+", "+str(time[1])+"]"
                                if specie == "time peak": time = str(time)
                                simulation.simulation_file['TIME FRAME'][str(specie)] = time
                            simulation.simulation_file.write(open(simulation.simulation_path, 'w'))
                            
                # Update the plotted modes
                self.display_plottedModesAndExperiments() 
                        
                # Display the saturated fluxes on the GUI
                self.update_treeViewSaturatedFluxes()
                    
                # Update the <graph> class so the option window can use them
                Graph.update() 
                
                # When finished show the canvas and change the cursor back.
                self.show_progress("finished", poppedout_id) 
                
        # If no simulations have been selected, clear the current figure
        if self.root.input_files == [] or self.root.Research.input_files == []:
            Graph.ax.clear()
            self.show_progress("nothing", poppedout_id) 
        
        # Update screen
        if poppedout_id==None: self.Canvas[axis_id].draw_idle() 
        if poppedout_id!=None: self.root.canvasPoppedOut[poppedout_id].draw_idle()
        if True: return

#==========================================
# Change what is displayed on the GUI
#==========================================

    def display_plottedModesAndExperiments(self):
        ''' When the graph is plotted, adjust the menus on the GUI of the modes and experiments. '''
        
        # Get the experiment plotted by the function
        self.options_experiments = ["All experiments"] + [ e.id for e in self.root.Research.experiments ]
        if self.experiment==None:
            self.experiment = self.root.Research.experiments[0]
            self.var_experiment.set(self.options_experiments[0])        
        if self.experiment.id not in self.options_experiments:
            self.experiment = self.root.Research.experiments[0]
            self.var_experiment.set(self.options_experiments[0])      
        if self.experiment_id not in self.options_experiments:
            self.experiment = self.root.Research.experiments[0]
            self.var_experiment.set(self.options_experiments[0])  
        if self.experiment_id == "All experiments":
            self.experiment = self.root.Research.experiments[0]
        for experiment in self.root.Research.experiments:
            if experiment.id == self.experiment_id:
                self.experiment = experiment
                self.var_experiment.set(self.experiment.id)
            
        self.mnu_experiment['menu'].delete(0, 'end')
        for experiment in self.options_experiments:
            self.mnu_experiment['menu'].add_command(label=experiment, command=tk._setit(self.var_experiment, experiment))

        # Get the simulation plotted by the function
        self.options_simulationsids = ["All simulations"] + [ s.id for s in self.experiment.simulations ]  
        self.options_simulations = ["All simulations"] + [ s.id.split("__")[-1] for s in self.experiment.simulations ]  
        self.options_simulations = ["All simulations"] + [ v for v in self.experiment.variedValues ]  
        self.options_simulations = [ s.replace('\\', '').replace(',', '').replace('$', '') for s in self.options_simulations ] 
        if self.simulation==None:
            self.simulation = self.experiment.simulations[0]
            self.var_simulation.set(self.options_simulations[0])
        if self.simulation.id not in self.options_simulationsids:
            self.simulation = self.experiment.simulations[0]
            self.var_simulation.set(self.options_simulations[0])
        if self.simulation_id not in self.options_simulationsids:
            self.simulation = self.experiment.simulations[0]
            self.var_simulation.set(self.options_simulations[0])
        for simulation in self.experiment.simulations:
            if simulation.id == self.simulation_id:
                self.simulation = simulation
                self.var_simulation.set(self.options_simulations[1+self.experiment.simulations.index(simulation)])
        self.mnu_simulation['menu'].delete(0, 'end')
        for simulation in self.options_simulations:
            self.mnu_simulation['menu'].add_command(label=simulation, command=tk._setit(self.var_simulation, simulation))
        if True: return
    
#################################################################
#                          Classes
#################################################################

# Make the axis instance and the required attributes for the options window
# Either the plotting function is a class or we manually give it some class attributes like here
class graph:
    
    def __init__(self, figure, var_plot, i):
        plt.figure(figure.number)
        grid_specifications = gridspec.GridSpec(1, 1, figure=figure)
        grid_specifications.update(top=0.95, left=0.15, right=0.85, bottom=0.15)
        self.ax = plt.subplot(grid_specifications[0])
        self.figure = figure
        self.range = {"units" : "normalized", "x_scale" : "linear", "y_scale" : "linear"}    
        self.label = {"x" : None, "y" : None, "title" : None}
        self.layout = {"fontsize" : "N.A.", 'handlelength' : "N.A."}
        self.x_key  = "x"
        self.y_key  = "y"
        self.x_name = "Heat flux"
        self.y_name = "Frequency"
        self.range["x"] = None
        self.range["y"] = None
        
        # Save the graph and id
        self.var_plot = var_plot
        self.id = i
        return 

    #---------------
    def load_defaults(self):     
        self.range["x"] = None
        self.range["y"] = None
        for key in self.label.keys(): 
            self.label[key] = None
        return
    
    #---------------
    def update(self):
        
        # Get the ranges, labels and titles
        self.range["x"]     = self.ax.get_xlim()
        self.range["y"]     = self.ax.get_ylim()
        self.label["x"]     = self.ax.get_xlabel()
        self.label["y"]     = self.ax.get_ylabel()
        self.label["title"] = self.ax.get_title()
        
        # The titles of the option window sections
        if self.var_plot.get()==1: self.y_name = "Heat flux"
        if self.var_plot.get()==2: self.y_name = "Particle flux"
        if self.var_plot.get()==2: self.y_name = "Momentum flux"
        return


class EntryPopup(ttk.Entry):

    def __init__(self, root, tree, specie_id, iid, columnid, text, values, update_treeViewAfterChange, **kw):
        super().__init__(tree, style='tree.TEntry', **kw)
        
        # Save the information about the treeview
        self.specie_id = specie_id
        self.tree = tree 
        self.update_tree = update_treeViewAfterChange
        self.iid = iid
        self.columnid = columnid
        self.text = text
        self.values = list(values)
        self['exportselection'] = False
        
        # Insert the text of the cell in the entry widget: either its text or one of the values
        if "0" in self.columnid: self.insert(0, text) 
        for i in range(1,len(values)+1):
            if str(i) in self.columnid: self.insert(0, values[i-1]) 

        # Bind the key presses to the entry widget
        self.focus_force()
        self.bind("<Return>", self.on_return)
        self.bind("<Control-a>", self.select_all)
        self.bind("<Escape>", lambda *ignore: self.destroy())
        return 
    
    #-------------------------
    def on_return(self, event=None):
        
        # Get the text/value that was changed  
        if "0" in self.columnid: self.text = self.get()
        for i in range(1,len(self.values)+1):
            if str(i) in self.columnid: 
                self.values[i-1] = self.get()
                
        # Replace the text/value in the treeview with this new text/value
        self.tree.item(self.iid, text=self.text, values = self.values) 

        # Also update the t_range
        simulation_iid = self.tree.selection()[0]
        experiment_iid = self.tree.parent(simulation_iid)
        self.update_tree(self.values, simulation_iid, experiment_iid, self.specie_id)
        self.destroy()
        return 
    
    #-------------------------
    def select_all(self, *ignore):
        ''' Set selection on the whole text '''
        self.selection_range(0, 'end')

        # returns 'break' to interrupt default key-bindings
        return 'break'

