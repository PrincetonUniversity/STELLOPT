
# TODO: kx ky tabs didn't update

#################################################################
#                   CLASS FOR THE FIRST TAB
#################################################################
''' TAB CONVERGENCE

Investigate the convergence in time and space of modes simulated with stella.


Absolute path of class TabSelectedFiles:
----------------------------------------
root.tab_Convergence
'''


# Load modules
import numpy as np
import tkinter as tk
from tkinter import ttk

# Personal modules
from stellapy.GUI.graph_tools import Progress
from stellapy.GUI.graph_tools import display_information
from stellapy.GUI.graph_tools import PAD_TITLE2, PAD_LABEL2, PAD_ENTRY2, PAD_ENTRYR #@unusedimport
from stellapy.GUI.graph_tools import CanvasForGraphs, PoppedOutWindow

# For the figures
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Plots
from stellapy.plot import plot_frequencyVsTime
from stellapy.plot import plot_potentialVsZ

#################################################################
#                   CLASS FOR THE FIRST TAB
#################################################################
class TabConvergence:
    '''
    Initiate the frames on tab3 "Convergence" to investigate the convergence 
    in time and space of modes simulated with stella.
    
    
    Attributes: objects
    -------------------
    class_inputParameters: InputParameters
    class_simulations:     Simulations
    class_progress:        Progress
    

    Parent widgets
    --------------
    root:               Tk()
    tabheader:          ttk.Notebook(root)
    tab3:               ttk.Frame(tabheader)
    '''
 
#################################################################
#                          WIDGETS
#################################################################

    def __init__(self, tab3): 

        #======================================================
        # Safe info from the tab3 and root that can be passed on
        #======================================================

        self.window = tab3                       # Direct accesses to the window of tab3 "Convergence" 
        self.root = tab3.root                    # Needed to center windows based on the root screen
        self.dict_awthemes = tab3.root.awthemes  # To make the tk widgets look like the ttk widgets

        #===========
        # VARIABLES
        #===========
        
        self.btn_popout = "Yes"            # Add the button to popout windows on the toolbar
        self.btn_reread = "Yes"            # Add the button to reread the data on the toolbar

        self.plotted_ykey = "gamma"         # Keep track of what is plotted now
        self.widgets = None
        self.plottedModes = {}
        self.initiated_canvas = None
        self.plotted_inputFiles = []        # Remember which input_files are plotted right now
        self.changed_folder = True
        self.folder = None                  # Option to change which folder is plotted
        self.commonPrefix = ""              # Remove the parent path of the simulation folder
        self.options_parameters = [" - "]
        self.options_kx = ["0.0"]
        self.options_ky = ["0.0", "100"]
        self.options_experiment = ["First plot data"]
        self.options_simulation = ["First plot data"]
        self.options_simulationsids = []
        self.experiment = None
        self.experiment_id = None
        self.simulation = None
        self.simulation_id = None
        self.plotted_modes = "unstable"
        self.axis_id = 0
        self.modes = {}
        self.reset = True
        self.replot = False
        
        #===========
        # FRAMES
        #===========

        # Create the frames and add them to the tab3
        self.frame_graph1   = ttk.LabelFrame(tab3, text="   Time evolution of the modes  ", **self.dict_awthemes['labelframe2'])
        self.frame_graph2   = ttk.LabelFrame(tab3, text="   Parallel mode structure of the modes  ", **self.dict_awthemes['labelframe2'])
        self.frame_options  = ttk.Frame(tab3)       
        
        # Configure the frames     
        tk.Grid.rowconfigure(   tab3, 0, weight=1) 
        tk.Grid.columnconfigure(tab3, 0, weight=10, uniform="columns") 
        tk.Grid.columnconfigure(tab3, 1, weight=3,  uniform="columns")
            
        # Add the two frames to the main window
        self.frame_graph1.grid( in_=tab3, row=0, column=0, padx=(20,5), pady=(20,20),stick='NSEW')
        self.frame_graph2.grid( in_=tab3, row=0, column=0, padx=(20,5), pady=(20,20),stick='NSEW')
        self.frame_options.grid(in_=tab3, row=0, column=1, padx=(5,20), pady=(20,20),stick='NSEW')        
        self.frame_graph2.grid_remove() # Hide this canvas until we plot figure 2
        
        # Initiate the widgets in the frame with the options
        self.initiate_frameOptions()                
         
        # Bind focus on every click on a widget: this makes sure the entry widgets become unfocused
        self.frame_graph1.bind("<1>", lambda event: self.frame_graph1.focus_set())
        self.frame_graph2.bind("<1>", lambda event: self.frame_graph2.focus_set())
        self.frame_options.bind("<1>", lambda event: self.frame_options.focus_set()) 
        self.subframe_simulations.bind("<1>", lambda event: self.subframe_simulations.focus_set()) 
        self.subframe_plots.bind("<1>", lambda event: self.subframe_plots.focus_set()) 
        self.subframe_modes.bind("<1>", lambda event: self.subframe_modes.focus_set()) 
        self.subframe_convergence.bind("<1>", lambda event: self.subframe_convergence.focus_set()) 
         
        #=============
        # LOAD FIGURE
        #=============
         
        # When the tab is visible, make sure the correct figure is loaded
        self.frame_options.bind("<Visibility>", self.load_figure)
        if True: return

#=======================================
# Initiate site bar with extra options
#=======================================

    def initiate_canvas(self):
        ''' Initiate the canvas. 
        
        This method gets called by stella_GUI.py after initiating the tabs to make sure 
        the attribute self.root.tab_Profiles already exists. 
        
        The toolbar requires the class to have a reset_graph and popout_window function
        The optionswindow requires the Graph objects to have attributes: ax, x_name, x_key, y_name, y_key, range, label
        '''
        
        # Create a list for the graph classes, create the figure for the canvas and initiate the canvas.
        self.initiated_canvas = True
        self.Graph = [None, None]      
        self.Canvas = [None, None]      
        self.figure1 = plt.figure("FrequencyVersusTime") 
        self.figure2 = plt.figure("PotentialVersusZ") 
        self.figure1.set_tight_layout(False) 
        self.figure2.set_tight_layout(False)
        
        # Put the canvas on the screen
        CanvasForGraphs(self.root, self.frame_graph1, self.root.tab_Convergence, self.figure1, axis_id=0)
        CanvasForGraphs(self.root, self.frame_graph2, self.root.tab_Convergence, self.figure2, axis_id=1)
        
        # Initiate the graph
        self.initiate_plottingClass()
    
    #-----------------------------------------------------  
    def initiate_plottingClass(self, figure=None):

        # Change the canvas color
        self.figure1.patch.set_facecolor(self.root.color['canvas'])  
        self.figure2.patch.set_facecolor(self.root.color['canvas']) 
        
        # Set the color of the axis
        plt.rcParams['text.color'] = self.root.color['fg']
        plt.rcParams['axes.edgecolor'] = self.root.color['fg']
        plt.rcParams['axes.labelcolor'] = self.root.color['fg']
        plt.rcParams['xtick.color'] = self.root.color['fg']
        plt.rcParams['ytick.color'] = self.root.color['fg']
        plt.rcParams['axes.facecolor'] = self.root.color['canvas']
        plt.rcParams['figure.facecolor'] = self.root.color['canvas']

        # Initiate the plotting class for the main window
        if figure is None:
            # Make the axis instance and the required attributes for the options window through the class <graph>
            self.Graph[0] = graph(self, self.figure1, self.var_plot, 0)
            self.Graph[1] = graph(self, self.figure2, self.var_plot, 1)

            # Put the empty figure on the GUI
            self.Canvas[0].draw_idle(); self.root.update_idletasks()
            self.Canvas[1].draw_idle(); self.root.update_idletasks()

        # Or initiate the plotting class for a popped out window
        if figure is not None:
            # Make the axis instance and the required attributes for the options window through the class <graph>
            self.root.graph_poppedOut.append(graph(self, figure, self.var_plot, 0))
        
        if True: return
        
#=======================================
# Initiate site bar with extra options
#=======================================

    def initiate_frameOptions(self):

        # Create the sub_labelframes in options plus a button to apply any changes
        self.subframe_simulations   = ttk.LabelFrame(self.frame_options, text="   Simulations  ", **self.dict_awthemes['labelframe2'])
        self.subframe_plots         = ttk.LabelFrame(self.frame_options, text="   Plot  ", **self.dict_awthemes['labelframe2'])
        self.subframe_modes         = ttk.LabelFrame(self.frame_options, text="   Modes   ", **self.dict_awthemes['labelframe2'])
        self.subframe_convergence   = ttk.LabelFrame(self.frame_options, text="   Convergence   ", **self.dict_awthemes['labelframe2'])
        self.frame_progress         = ttk.Frame(self.frame_options)
        
        # Add the progress object
        self.class_progress = Progress(self, anchor="fill", length=300)
        self.class_progress.move(0,"Please select simulations.")
        
        # Configure the frame
        tk.Grid.rowconfigure(   self.frame_options, 0, weight=0) 
        tk.Grid.rowconfigure(   self.frame_options, 1, weight=0) 
        tk.Grid.rowconfigure(   self.frame_options, 2, weight=0) 
        tk.Grid.rowconfigure(   self.frame_options, 3, weight=1)   
        tk.Grid.columnconfigure(self.frame_options, 0, weight=1) 
        
        # Add the three labelframes and the frame for the button to <frame_options>
        self.subframe_simulations.grid( row=0, column=0, padx=(5,20), pady=(0,10), stick='NSEW')
        self.subframe_plots.grid(       row=1, column=0, padx=(5,20), pady=(10,10),stick='NSEW')
        self.subframe_modes.grid(       row=2, column=0, padx=(5,20), pady=(10,10),stick='NSEW')
        self.subframe_convergence.grid( row=3, column=0, padx=(5,20), pady=(10,10), stick='NSEW')
        self.frame_progress.grid(       row=4, column=0, padx=(5,20), pady=(10,0),stick='NSEW')
        
        # Fill the four subframes with widgets
        self.initiate_frameOptions_simulationsFrame()
        self.initiate_frameOptions_plotFrame()
        self.initiate_frameOptions_modesFrame()
        self.initiate_frameOptions_convergenceFrame()
        
    #----------- Simulations -------------
    def initiate_frameOptions_simulationsFrame(self):
        
        # Choice between the experiments
        self.var_experiment = tk.StringVar(value=self.options_experiment[0]); width=15
        self.lbl_experiment = ttk.Label(self.subframe_simulations, text="Experiment:")
        self.mnu_experiment = ttk.OptionMenu(self.subframe_simulations, self.var_experiment, self.options_experiment[0], *self.options_experiment, style='option.TMenubutton')
        self.mnu_experiment["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
        self.mnu_experiment.config(width=width)
        self.var_experiment.trace('w', self.change_plottedExperiment) # link function to a change of the dropdown options
        
        # Choice between the simulations
        self.var_simulation = tk.StringVar(value=self.options_simulation[0])
        self.lbl_simulation = ttk.Label(self.subframe_simulations, text="Simulation:")
        self.mnu_simulation = ttk.OptionMenu(self.subframe_simulations, self.var_simulation, self.options_simulation[0], *self.options_simulation, style='option.TMenubutton')
        self.mnu_simulation["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
        self.mnu_simulation.config(width=width)
        self.var_simulation.trace('w', self.change_plottedSimulation) # link function to a change of the dropdown options
            
        # Show the range of the kx and the k values
        self.var_kxmin = tk.StringVar(value=self.options_kx[0]); width=4
        self.var_kxmax = tk.StringVar(value=self.options_kx[-1])
        self.var_kymin = tk.StringVar(value=self.options_ky[0])
        self.var_kymax = tk.StringVar(value=self.options_ky[-1])
        self.lbl_kxmin = ttk.Label(self.subframe_simulations, text="kx min:", style='opt_sign.TLabel')
        self.lbl_kxmax = ttk.Label(self.subframe_simulations, text=" kx max:", style='opt_sign.TLabel')
        self.lbl_kymin = ttk.Label(self.subframe_simulations, text="ky min:", style='opt_sign.TLabel')
        self.lbl_kymax = ttk.Label(self.subframe_simulations, text=" ky max:", style='opt_sign.TLabel')
        self.mnu_kxmin = ttk.OptionMenu(self.subframe_simulations, self.var_kxmin, self.options_kx[0], *self.options_kx, style='option.TMenubutton')
        self.mnu_kxmin["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
        self.mnu_kxmin.config(width=width)
        self.mnu_kxmax = ttk.OptionMenu(self.subframe_simulations, self.var_kxmax, self.options_kx[0], *self.options_kx, style='option.TMenubutton')
        self.mnu_kxmax["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
        self.mnu_kxmax.config(width=width)
        self.mnu_kymin = ttk.OptionMenu(self.subframe_simulations, self.var_kymin, self.options_ky[0], *self.options_ky, style='option.TMenubutton')
        self.mnu_kymin["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
        self.mnu_kymin.config(width=width)
        self.mnu_kymax = ttk.OptionMenu(self.subframe_simulations, self.var_kymax, self.options_ky[1], *self.options_ky, style='option.TMenubutton')
        self.mnu_kymax["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
        self.mnu_kymax.config(width=width)
        self.var_kxmin.trace('w', self.update_selectedModes) # link function to a change of the dropdown options
        self.var_kxmax.trace('w', self.update_selectedModes) # link function to a change of the dropdown options
        self.var_kymin.trace('w', self.update_selectedModes) # link function to a change of the dropdown options
        self.var_kymax.trace('w', self.update_selectedModes) # link function to a change of the dropdown options

        # Configure the frame
        tk.Grid.columnconfigure(self.subframe_simulations, 0, weight=0) 
        tk.Grid.columnconfigure(self.subframe_simulations, 1, weight=1, uniform="2") 
        tk.Grid.columnconfigure(self.subframe_simulations, 2, weight=0) 
        tk.Grid.columnconfigure(self.subframe_simulations, 3, weight=1, uniform="2") 
        
        # Place the widgets in the frame
        self.lbl_experiment.grid(row=0, column=0, stick='NSEW', padx=(0,0), pady=(2,2), columnspan=2)
        self.mnu_experiment.grid(row=0, column=2, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=2, ipady=1, columnspan=2)
        self.lbl_simulation.grid(row=1, column=0, stick='NSEW', padx=(0,0), pady=(2,2), columnspan=2)
        self.mnu_simulation.grid(row=1, column=2, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=2, ipady=0, columnspan=2)
        self.lbl_kxmin.grid(row=2, column=0, stick='W',    padx=(0,0), pady=(2,2), ipadx=1, ipady=0)
        self.mnu_kxmin.grid(row=2, column=1, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=1, ipady=0)
        self.lbl_kxmax.grid(row=2, column=2, stick='W',    padx=(0,0), pady=(2,2), ipadx=1, ipady=0)
        self.mnu_kxmax.grid(row=2, column=3, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=1, ipady=0)
        self.lbl_kymin.grid(row=3, column=0, stick='W',    padx=(0,0), pady=(2,2), ipadx=1, ipady=0)
        self.mnu_kymin.grid(row=3, column=1, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=1, ipady=0)
        self.lbl_kymax.grid(row=3, column=2, stick='W',    padx=(0,0), pady=(2,2), ipadx=1, ipady=0)
        self.mnu_kymax.grid(row=3, column=3, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=1, ipady=0)

    
        
    #------------- Plots -------------    
    def initiate_frameOptions_plotFrame(self):
        
        # Set the default value of the plot out of range so we don't start with a plot
        self.var_plot  = tk.IntVar(value=2) 
        
        # Add the possible plots in two sections: "Time convergence" and "Space convergence"
        self.lbl_time  = ttk.Label(self.subframe_plots, text="Convergence in time", style='prefTitle.TLabel')
        self.rbn_omega = ttk.Radiobutton(self.subframe_plots, text='  Frequency versus time')
        self.rbn_gamma = ttk.Radiobutton(self.subframe_plots, text='  Growth rate versus time')
        self.lbl_space = ttk.Label(self.subframe_plots, text="Confinement along the fluxtube", style='prefTitle.TLabel')
        self.rbn_real  = ttk.Radiobutton(self.subframe_plots, text='  Real potential versus z')
        self.rbn_imag  = ttk.Radiobutton(self.subframe_plots, text='  Imaginary potential versus z',)
        self.rbn_phi2  = ttk.Radiobutton(self.subframe_plots, text='  Potential squared versus z',)
        
        # Add the values and commands to the radiobuttons
        self.rbn_omega.config(value=1, variable=self.var_plot, command=lambda: self.change_plot())
        self.rbn_gamma.config(value=2, variable=self.var_plot, command=lambda: self.change_plot())
        self.rbn_real.config( value=3, variable=self.var_plot, command=lambda: self.change_plot())
        self.rbn_imag.config( value=4, variable=self.var_plot, command=lambda: self.change_plot())
        self.rbn_phi2.config( value=5, variable=self.var_plot, command=lambda: self.change_plot())
        
        # Configure the frame
        tk.Grid.columnconfigure(self.subframe_plots, 0, weight=1) 
        
        # Add the options to the frame
        self.lbl_time.grid(   row=0, column=0, **PAD_TITLE2)
        self.rbn_omega.grid(  row=1, column=0, **PAD_LABEL2)
        self.rbn_gamma.grid(  row=2, column=0, **PAD_LABEL2)
        # Add the options to the frame
        self.lbl_space.grid(  row=3, column=0, **PAD_TITLE2)
        self.rbn_real.grid(   row=4, column=0, **PAD_LABEL2)
        self.rbn_imag.grid(   row=5, column=0, **PAD_LABEL2)
        self.rbn_phi2.grid(   row=6, column=0, **PAD_LABEL2)
        

    #-------------- Modes -------------
    def initiate_frameOptions_modesFrame(self):
        
        # Set the default value of the plotted modes 
        self.var_stab = tk.IntVar(value=2)

        # Stability and convergence
        self.rbn_stable      = ttk.Radiobutton(self.subframe_modes, text='  Stable modes')
        self.rbn_unstable    = ttk.Radiobutton(self.subframe_modes, text='  Unstable modes')
        self.rbn_converged   = ttk.Radiobutton(self.subframe_modes, text='  Unstable modes (converged)')
        self.rbn_unconverged = ttk.Radiobutton(self.subframe_modes, text='  Unstable modes (unconverged)')
        
        # Add the values and commands to the radiobuttons
        self.rbn_stable.config(     value=1, variable=self.var_stab, command=self.change_plottedmodes)
        self.rbn_unstable.config(   value=2, variable=self.var_stab, command=self.change_plottedmodes)
        self.rbn_converged.config(  value=3, variable=self.var_stab, command=self.change_plottedmodes)
        self.rbn_unconverged.config(value=4, variable=self.var_stab, command=self.change_plottedmodes)
        
        # Configure the frame
        tk.Grid.columnconfigure(self.subframe_plots, 0, weight=1) 
        tk.Grid.columnconfigure(self.subframe_plots, 1, weight=0) 
        
        # Add the options to the frame
        self.rbn_stable.grid(     row=4, column=0, **PAD_LABEL2, columnspan=2)
        self.rbn_unstable.grid(   row=5, column=0, **PAD_LABEL2, columnspan=2)
        self.rbn_unconverged.grid(row=6, column=0, **PAD_LABEL2, columnspan=2)
        self.rbn_converged.grid(  row=7, column=0, **PAD_LABEL2, columnspan=2)


    #---------- Convergence ---------------        
    def initiate_frameOptions_convergenceFrame(self):
        
        # Some information about the convergence of the modes will be displayed:
        self.txt_convergence1 = tk.StringVar(value="")
        self.txt_convergence2 = tk.StringVar(value="No simulations are selected.") 
        self.txt_convergence3 = tk.StringVar(value="")      
        self.txt_convergence4 = tk.StringVar(value="")
        
        # Create 4 label widgets with the linked StringVar variables shown above
        self.lbl_convergence1 = ttk.Label(self.subframe_convergence, textvariable=self.txt_convergence1, style='opt_paraC.TLabel')        
        self.lbl_convergence2 = ttk.Label(self.subframe_convergence, textvariable=self.txt_convergence2, style='opt_valueC.TLabel') 
        self.lbl_convergence3 = ttk.Label(self.subframe_convergence, textvariable=self.txt_convergence3, style='opt_paraC.TLabel')  
        self.lbl_convergence4 = ttk.Label(self.subframe_convergence, textvariable=self.txt_convergence4, style='opt_valueC.TLabel')

        # Configrue the frame
        tk.Grid.rowconfigure(   self.subframe_convergence, 0, weight=0) 
        tk.Grid.rowconfigure(   self.subframe_convergence, 1, weight=0) 
        tk.Grid.rowconfigure(   self.subframe_convergence, 2, weight=0)
        tk.Grid.rowconfigure(   self.subframe_convergence, 3, weight=0)
        tk.Grid.columnconfigure(self.subframe_convergence, 0, weight=1) 
        
        # Add the label widgets to the frame
        self.lbl_convergence1.grid( in_=self.subframe_convergence, row=0, column=0, padx=(5,5), pady=(5,2),stick='NSEW')
        self.lbl_convergence2.grid( in_=self.subframe_convergence, row=1, column=0, padx=(5,5), pady=(2,2),stick='NSEW')
        self.lbl_convergence3.grid( in_=self.subframe_convergence, row=2, column=0, padx=(5,5), pady=(2,2),stick='NSEW')
        self.lbl_convergence4.grid( in_=self.subframe_convergence, row=3, column=0, padx=(5,5), pady=(2,2),stick='NSEW')
        if True: return
        
#################################################################
#                          METHODS
#################################################################

    def load_figure(self, *args):
        '''  When the tab is visible, make sure the correct figure is loaded. 
        If the simulations were changed, replot the graph. '''
        
        # If the canvas isn't loaded, load them
        if self.initiated_canvas==None:
            self.initiate_canvas()
            
        # When hitting F1 replot the graph
        def plot_graph(*args): self.replot=True; self.plot_graph(self.axis_id, None); return
        self.root.bind('<KeyPress-F1>', plot_graph)
        self.root.bind('<Control-s>', plot_graph)
        
        # Load figure
        if self.axis_id==0: plt.figure("FrequencyVersusTime")
        if self.axis_id==1: plt.figure("PotentialVersusZ")
        
        # Connect the Progress class to the research object
        self.root.Research.Progress = self.class_progress
        for experiment in self.root.Research.experiments:
            experiment.Progress = self.class_progress
            for simulation in experiment.simulations:
                simulation.Progress = self.class_progress
                
        # Make sure the plotted experiment is one of the experiments
        if self.experiment != None:
            self.options_experiments = [ e.id for e in self.root.Research.experiments ]
            if self.experiment.id not in self.options_experiments:
                self.experiment = None
        
        # Make sure we display the correct experiments
        if self.root.input_files != []:
            self.display_plottedModesAndExperiments()
        self.show_progress("finished", None) 
        return
    
    #------------------------        
    def popout_window(self, axis_id):
        '''Replot the figure in a seperate window when the button on the toolbar is clicked.'''
        
        # Create a popped out window
        poppedout_window = PoppedOutWindow(self.root)
        poppedout_id = poppedout_window.poppedout_id
        
        # Add a fitting title
        if self.var_plot.get() == 1: poppedout_window.set_title("Stellapy: Time evolution of the frequency of the modes")
        if self.var_plot.get() == 2: poppedout_window.set_title("Stellapy: Time evolution of the growth rate of the modes")
        if self.var_plot.get() == 3: poppedout_window.set_title("Stellapy: Parallel mode structure of the real part of the electrostatic fluctuations")
        if self.var_plot.get() == 4: poppedout_window.set_title("Stellapy: Parallel mode structure of the imaginary part of the electrostatic fluctuations")
        if self.var_plot.get() == 5: poppedout_window.set_title("Stellapy: Parallel mode structure of the electrostatic fluctuations")
        
        # Initiate the plotting class
        self.initiate_plottingClass(figure=poppedout_window.figure)

        # Initiate the canvas in self.frame_graph          
        CanvasForGraphs(self.root, poppedout_window.frame, self.root.tab_Convergence, poppedout_window.figure, poppedout_id=poppedout_id)
        
        # Parse the data from the main canvas to the popped out canvas
        self.root.graph_poppedOut[poppedout_id].kx        = self.Graph[0].kx
        self.root.graph_poppedOut[poppedout_id].ky        = self.Graph[0].ky
        self.root.graph_poppedOut[poppedout_id].range     = self.Graph[0].range
        self.root.graph_poppedOut[poppedout_id].label     = self.Graph[0].label
        self.root.graph_poppedOut[poppedout_id].update_quantitiesAndKeys()
        self.root.graph_poppedOut[poppedout_id].grid_specifications.update(top=0.92, left=0.1, right=0.8, bottom=0.15)
        
        # Now plot the figure on the popped out canvas
        self.plot_graph(None,poppedout_id)
        
    #--------------------------------    
    def change_plot(self):
        ''' Change the plot ("omega vs t" or "phi vs z"). '''
        
        # Switch between the canvasses and update the quantities and keys
        if self.var_plot.get() in [1,2]: 
            self.axis_id = 0
            self.frame_graph1.grid()
            self.frame_graph2.grid_remove()
            self.Graph[0].load_defaults()
            self.Graph[0].update_quantitiesAndKeys()
            self.replot = True
            self.plot_graph(self.axis_id, None)
        if self.var_plot.get() in [3,4,5]: 
            self.axis_id = 1
            self.frame_graph2.grid()
            self.frame_graph1.grid_remove()
            self.Graph[1].load_defaults()
            self.Graph[1].update_quantitiesAndKeys()
            self.replot = True
            self.plot_graph(self.axis_id, None)
        return

    #----------------------------
    def change_plottedmodes(self, *args):
        
        # Change the modes that will be plotted
        if self.var_stab.get()==1: self.plotted_modes = "stable"
        if self.var_stab.get()==2: self.plotted_modes = "unstable"
        if self.var_stab.get()==3: self.plotted_modes = "converged"
        if self.var_stab.get()==4: self.plotted_modes = "unconverged" 
        
        # Replot the graph
        self.replot = True
        self.plot_graph(self.axis_id, None)
        return
    
    #-----------------------
    def change_plottedExperiment(self, *args):
        
        # Change the experiment that is plotted
        self.experiment_id = self.var_experiment.get()
    
        # Get a reference to the experiment
        self.experiment = None
        for experiment in self.root.Research.experiments:
            if experiment.id == self.experiment_id:
                self.experiment = experiment
        if self.experiment==None:
            self.experiment = self.root.Research.experiments[0]
            self.experiment_id = self.experiment.id
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
        return 
    
    #---------------------------------------------------
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
                
#===================
# Update the graph 
#===================

    def reset_graph(self, axis_id=None, poppedout_id=None):
        
        # Perhaps the button was clicked on the popped out window
        if poppedout_id == None: Graph = self.Graph[axis_id]
        if poppedout_id != None: Graph = self.root.graph_poppedOut[poppedout_id]
            
        # Reset the axis
        Graph.load_defaults()
        self.replot="DONT"; self.var_kxmin.set(0.0)
        self.replot="DONT"; self.var_kxmax.set(0.0)
        self.replot="DONT"; self.var_kymin.set(0.0)
        self.replot="DONT"; self.var_kymax.set(100)
                    
        # Replot the graph
        self.replot = True
        self.plot_graph(axis_id, poppedout_id)

    #----------------------  
    def plot_graph(self,  axis_id=None, poppedout_id=None):
        
        # Perhaps the button was clicked on the popped out window
        if poppedout_id == None: Graph = self.Graph[axis_id]
        if poppedout_id != None: Graph = self.root.graph_poppedOut[poppedout_id]
        
        # Only proceed if there are input_files
        if self.root.input_files != []:
            
            # Get the input_files that need to be plotted
            input_files = self.root.Research.input_files
            
            # If the selection of input_files changed, then reset the axis
            if (self.plotted_inputFiles != input_files):
                Graph.load_defaults()
            
            # If there are files plot them, if they changed plot them, if its a popout plot them
            if (input_files != [] and self.plotted_inputFiles != input_files) or poppedout_id!=None or self.replot==True:
                
                # Clear the axis and remember which input_files are plotted
                Graph.ax.clear()
                self.plotted_inputFiles = input_files.copy() 
                self.replot = False
                
                # While plotting show a progress bar turn the cursor into a waiting cursor
                self.show_progress("start", poppedout_id) 
                
                # Plot the frequency versus time
                if self.var_plot.get() in [1,2]: 
                    plot_frequencyVsTime(\
                        # Specify which simulations to plot
                            research=self.root.Research,\
                            experiment_id=self.experiment_id,\
                            simulation_id=self.simulation_id,\
                        # Specify data range
                            x_range=Graph.range["x"],\
                            y_range=Graph.range["y"],\
                            quantity=Graph.y_key,\
                            units=Graph.range["units"],\
                            kx_range=Graph.kx,\
                            ky_range=Graph.ky,\
                            plotted_modes=self.plotted_modes,\
                        # Labels
                            x_label=Graph.label["x"],\
                            y_label=Graph.label["y"],\
                        # For the GUI the figure object already exists  
                            ax=Graph.ax,\
                            Progress=self.class_progress,\
                            root=self.root,\
                        # Appearance options
                            font_size=20,\
                            handle_length=1)
                if self.var_plot.get() in [3,4,5]: 
                    plot_potentialVsZ(\
                        # Specify which simulations to plot
                            research=self.root.Research,\
                            experiment_id=self.experiment_id,\
                            simulation_id=self.simulation_id,\
                        # Specify data range
                            units="normalized",\
                            x_quantity=Graph.x_key,\
                            y_quantity=Graph.y_key,\
                            x_range=Graph.range["x"],\
                            y_range=Graph.range["y"],\
                            plotted_modes=self.plotted_modes,\
                            kx_range=Graph.kx,\
                            ky_range=Graph.ky,\
                        # Labels
                            x_label=Graph.label["x"],\
                            y_label=Graph.label["y"],\
                            title=None,\
                        # For the GUI the figure object already exists 
                            ax=Graph.ax,\
                            Progress=self.class_progress,\
                            root=self.root,\
                        # Appearance options
                            font_size=20,\
                            handle_length=1)
                    
                # Update the Graph for the options window
                Graph.update_rangesAndLabels()

                # Update the plotted modes
                self.display_plottedModesAndExperiments()  
                    
                # Update the not converged simulations
                self.write_informationConvergence()
                
                # When finished show the canvas and change the cursor back.
                self.show_progress("finished", poppedout_id)                 
            
        # If no simulations have been selected, clear the current figure
        if self.root.input_files == [] or self.root.Research.input_files == []:
            Graph.ax.clear()
            Graph.load_defaults()
            self.show_progress("nothing", poppedout_id) 
            self.txt_convergence2.set("Please select simulations.")
    
        # Update screen
        if poppedout_id==None: self.Canvas[axis_id].draw_idle()
        if poppedout_id!=None: self.root.canvasPoppedOut[poppedout_id].draw_idle()
            
        return

    #----------------------------------
    def write_informationConvergence(self):
        
        # Make sure the elements are strings that are rounded to 2 digits
        vec_k            = [ str(round(k,2)) for k in self.simulation.vec_kAll ]
        vec_kUnstable    = [ str(round(k,2)) for k in self.simulation.vec_kUnstable ]  
        vec_kNotConverge = [ str(round(k,2)) for k in self.simulation.vec_kNotConverge ]
        vec_kStable      = [ str(round(k,2)) for k in self.simulation.vec_kStable ]
        
        # Print information if no modes are plotted
        if (self.var_stab.get()==1 and len(vec_kStable)==0) or (self.var_stab.get()!=1 and len(vec_kUnstable)==0):
            if self.var_stab.get()==1: # We want to see the stable modes, but there are none!
                self.txt_convergence1.set("From the "+str(len(vec_k))+" modes, all are unstable.")
            if self.var_stab.get()!=1: # We want to see the unstable modes, but there are none!
                self.txt_convergence1.set("From the "+str(len(vec_k))+" modes, all are stable.")
            self.txt_convergence2.set("There is nothing to plot.")
        
        # Overwrite the default message with information about the number of stable or unstable modes 
        else:  
            if len(vec_kStable) > 0 and self.var_stab.get()==1: # We want to see the stable modes
                txt_notconverge1 = " are stable:" if len(vec_kStable)>1 else " is stable:"
                txt_notconverge1 = "From the "+str(len(vec_k))+" modes, "+str(len(vec_kStable))+txt_notconverge1
                txt_notconverge2 = "[" + vec_kStable[0] + ", ..., " + vec_kStable[-1] + "]"
            if len(vec_kUnstable) > 0 and self.var_stab.get()!=1: # We want to see the unstable modes   
                txt_notconverge1 = " are unstable:" if len(vec_kUnstable)>1 else " is unstable:"
                txt_notconverge1 = "From the "+str(len(vec_k))+" modes, "+str(len(vec_kUnstable))+txt_notconverge1
                txt_notconverge2 = "[" + vec_kUnstable[0] + ", ..., " + vec_kUnstable[-1] + "]"
            self.txt_convergence1.set(txt_notconverge1)
            self.txt_convergence2.set(txt_notconverge2)
            
        # Print information about the number of unconverged modes
        if len(vec_kNotConverge) > 0:
            txt_notconverge3 = " modes have not converged:" if len(vec_kNotConverge)>1 else " mode has not converged:"
            txt_notconverge3 = "And " + str(len(vec_kNotConverge)) + txt_notconverge3
            txt_notconverge4 = "[" + "\n".join([", ".join(vec_kNotConverge[i:i+6]) for i in range(0,len(vec_kNotConverge),6)]) + "]"
            self.txt_convergence3.set(txt_notconverge3)
            self.txt_convergence4.set(txt_notconverge4)
        else:  
            self.txt_convergence3.set("All simulations have converged.")
            self.txt_convergence4.set("")


        return
       
    #-------------------------
    def display_plottedModesAndExperiments(self):
        ''' When the graph is plotted, adjust the menus on the GUI of the modes and experiments. '''
        
        print("--------------")
        # Get the experiment plotted by the function
        self.reset = False
        self.options_experiments = [ e.id for e in self.root.Research.experiments ]
        if self.experiment==None:
            self.experiment = self.root.Research.experiments[0]
            self.experiment_id = self.experiment.id
        if self.experiment.id not in self.options_experiments:
            self.experiment = self.root.Research.experiments[0]
            self.experiment_id = self.experiment.id
        for experiment in self.root.Research.experiments:
            if experiment.id == self.experiment_id:
                self.experiment = experiment
                self.var_experiment.set(self.experiment.id)
                
        self.mnu_experiment['menu'].delete(0, 'end')
        for experiment in self.options_experiments:
            self.mnu_experiment['menu'].add_command(label=experiment, command=tk._setit(self.var_experiment, experiment))

        # Get the simulation plotted by the function
        self.reset = False
        self.options_simulationsids = ["All simulations"] + [ s.id for s in self.experiment.simulations ]  
        self.options_simulations = ["All simulations"] + [ s.id.split("__")[-1] for s in self.experiment.simulations ]  
        self.options_simulations = ["All simulations"] + [ v for v in self.experiment.variedValues ]  
        self.options_simulations = [ s.replace('\\', '').replace(',', '').replace('$', '') for s in self.options_simulations ] 
        if self.simulation==None:
            self.simulation = self.experiment.simulations[0]
            self.simulation_id = self.simulation.id if self.simulation_id != "All simulations" else self.simulation_id
        if self.simulation.id not in self.options_simulationsids:
            self.simulation = self.experiment.simulations[0]
            self.simulation_id = self.simulation.id if self.simulation_id != "All simulations" else self.simulation_id
        for simulation in self.experiment.simulations:
            if simulation.id == self.simulation_id:
                self.simulation = simulation
                self.var_simulation.set(self.options_simulations[1+self.experiment.simulations.index(simulation)])
        self.mnu_simulation['menu'].delete(0, 'end')
        for simulation in self.options_simulations:
            self.mnu_simulation['menu'].add_command(label=simulation, command=tk._setit(self.var_simulation, simulation))
        
        # Update the options for the kx and ky ranges
        if self.simulation_id == "All simulations":
            self.options_kx = self.experiment.total_kx
            self.options_ky = self.experiment.total_ky
        if self.simulation_id != "All simulations":
            self.options_kx = self.simulation.vec_kx
            self.options_ky = self.simulation.vec_ky
        self.mnu_kxmin['menu'].delete(0, 'end')
        self.mnu_kxmax['menu'].delete(0, 'end')
        self.mnu_kymin['menu'].delete(0, 'end')
        self.mnu_kymax['menu'].delete(0, 'end')
        for kx in self.options_kx:
            self.mnu_kxmin['menu'].add_command(label=round(kx,2), command=tk._setit(self.var_kxmin, kx))
            self.mnu_kxmax['menu'].add_command(label=round(kx,2), command=tk._setit(self.var_kxmax, kx))
        for ky in self.options_ky:
            self.mnu_kymin['menu'].add_command(label=round(ky,2), command=tk._setit(self.var_kymin, ky))
            self.mnu_kymax['menu'].add_command(label=round(ky,2), command=tk._setit(self.var_kymax, ky))

        # Update the GUI variables
        if isinstance(self.Graph[self.axis_id].kx, list): 
            kx_min = round(np.max([self.Graph[self.axis_id].kx[0], self.options_kx[0]]),2)
            self.replot="DONT"; self.var_kxmin.set(kx_min)
            kx_max = round(np.min([self.Graph[self.axis_id].kx[1], self.options_kx[-1]]),2)
            self.replot="DONT"; self.var_kxmax.set(kx_max) 
        if isinstance(self.Graph[self.axis_id].ky, list): 
            ky_min = round(np.max([self.Graph[self.axis_id].ky[0], self.options_ky[0]]),2)
            self.replot="DONT"; self.var_kymin.set(ky_min)
            ky_max = round(np.min([self.Graph[self.axis_id].ky[1], self.options_ky[-1]]),2)
            self.replot="DONT"; self.var_kymax.set(ky_max) 
        
        self.reset = True
        return
    
#==========================================
# Change appearance of the graph
#==========================================
    def update_selectedModes(self, *args): 
        ''' Tell the graph class which modes to plot.'''
        
        # Get the kx and k ranges set in the GUI 
        kxmin = float(self.var_kxmin.get())
        kxmax = float(self.var_kxmax.get())
        kymin = float(self.var_kymin.get())
        kymax = float(self.var_kymax.get())
        
        # Make sure we don't scan over both kx and ky
        if kxmin != kxmax and kymin != kymax:
            message = """You can not scan over kx and ky simultaneously\n
            Please make sure one of the ranges is a fixed number."""
            display_information(self.root, "WARNING", message) 

        # Tell the graph class which modes to plot.
        self.Graph[0].kx = [ kxmin, kxmax ] 
        self.Graph[0].ky = [ kymin, kymax ] 
        self.Graph[1].kx = [ kxmin, kxmax ] 
        self.Graph[1].ky = [ kymin, kymax ] 
        if kxmin == kxmax:
            self.Graph[0].kx = kxmin
            self.Graph[1].kx = kxmin
        elif kymin == kymax:
            self.Graph[0].ky = kymin
            self.Graph[1].ky = kymin
            
        # Replot the graph if the change was made from the GUI
        if self.replot!="DONT":
            self.replot = True
            self.plot_graph(self.axis_id, None)
        if self.replot=="DONT": 
            self.replot = False
        if True: return

#################################################################
#                          Classes
#################################################################

# Make the axis instance and the required attributes for the options window
# Either the plotting function is a class or we manually give it some class attributes like here
class graph:
    
    def __init__(self, tab, figure, var_plot, i):
        
        # Attributes for every figure object
        plt.figure(figure.number)
        self.grid_specifications = gridspec.GridSpec(1, 1, figure=figure)
        self.grid_specifications.update(top=0.95, left=0.05, right=0.8, bottom=0.1)
        self.ax = plt.subplot(self.grid_specifications[0])
        self.figure = figure
        self.label = {"x" : None, "y" : None, "title" : None}
        self.range = {"units" : "normalized", "x_scale" : "linear", "y_scale" : "linear", "x" : None, "y" : None}    
        self.layout = {"fontsize" : "N.A.", 'handlelength' : "N.A."}
        self.tab = tab
        
        # The x and y-axis quantities and keys
        if i==0:
            self.x_name = "Time"
            self.x_key  = "time"
            self.y_name = "Growthrate"
            self.y_key  = "gamma"
        if i==1:     
            self.x_name = "z"
            self.x_key  = "z"              
            self.y_name = "phi_real"
            self.y_key  = "phi_real" 
            
        # Set the ranges and labels to None when the graph is reset, this will triger the defaults
        self.load_defaults()
            
        # Save the graph and id
        self.var_plot = var_plot
        self.id = i
        return 
    
    #-------------------
    def load_defaults(self):
        self.kx = 0.0
        self.ky = [0,100]
        self.range["x"] = None
        self.range["y"] = None
        for key in self.label.keys(): 
            self.label[key] = None
        return
    
    #-------------------
    def update_quantitiesAndKeys(self):
        if self.var_plot.get() == 1:     
            self.x_name = "Time"
            self.x_key  = "time"              
            self.y_name = "Omega"
            self.y_key  = "omega"    
        if self.var_plot.get() == 2:     
            self.x_name = "Time"
            self.x_key  = "time"              
            self.y_name = "Gamma"
            self.y_key  = "gamma"  
        if self.var_plot.get() == 3:     
            self.x_name = "z"
            self.x_key  = "z"              
            self.y_name = "phi_real"
            self.y_key  = "phi_real" 
        if self.var_plot.get() == 4:     
            self.x_name = "z"
            self.x_key  = "z"              
            self.y_name = "phi_imag"
            self.y_key  = "phi_imag"         
        if self.var_plot.get() == 5:     
            self.x_name = "z"
            self.x_key  = "z"              
            self.y_name = "phi2"
            self.y_key  = "phi2" 
        return
    #---------------
    def update_rangesAndLabels(self):
        
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
        
        # If we plot all simulations, make more room fot the labels
        if self.tab.simulation_id=="All simulations":
            self.grid_specifications.update(top=0.95, left=0.1, right=0.75, bottom=0.1)
        if self.tab.simulation_id!="All simulations":
            self.grid_specifications.update(top=0.95, left=0.1, right=0.8, bottom=0.1)
        return

    #----------------------------------
    def change_units(self):
        # Replot the graph
        self.load_defaults()
        self.ax.clear()
        self.tab.replot = True
        self.tab.plot_graph(self.tab.axis_id, None)
        self.update_rangesAndLabels
        
    



 





