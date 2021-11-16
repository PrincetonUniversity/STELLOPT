
# Load modules
import tkinter as tk
from tkinter import ttk 
from stellapy.GUI.graph_tools import Progress
from stellapy.GUI.interface.TabLinear.Canvasses import Canvasses
from stellapy.GUI.interface.TabLinear.Simulations import Simulations
from stellapy.GUI.interface.TabLinear.PlotLinearSpectrum import PlotLinearSpectrum
from stellapy.GUI.interface.TabLinear.PlotParameterInfluence import PlotParameterInfluence    

#################################################################
#                   CLASS FOR THE TAB
#################################################################
class TabLinear:


    def __init__(self, tab4): 

        #======================================================
        # Safe info from the tab and root that can be passed on
        #======================================================

        self.window = tab4                      # Direct accesses to the window of tab "Select Simulations" 
        self.root = tab4.root                   # Needed to center windows based on the root screen
        self.dict_awthemes = tab4.root.awthemes # To make the tk widgets look like the ttk widgets

        #===========
        # VARIABLES
        #===========
        
        # Only initiate the canvasses when the tab is used
        self.initiated_canvas = None
        
        # Store the Graph and Canvas objects in a list
        self.Graph = [None, None]      
        self.Canvas = [None, None]    

        #=======================================
        # FRAMES FOR OPTIONS, PLOT 1 and PLOT 2
        #========================================

        # Create the frames and add them to the tab3
        self.frame_graph1   = ttk.LabelFrame(self.window, text="   Linear spectrum  ", **self.dict_awthemes['labelframe2'])
        self.frame_graph2   = ttk.LabelFrame(self.window, text="   Parameter influence  ", **self.dict_awthemes['labelframe2'])
        self.frame_options  = ttk.Frame(self.window) 

        # Configure the frames     
        tk.Grid.rowconfigure(   self.window, 0, weight=0) 
        tk.Grid.rowconfigure(   self.window, 1, weight=1) 
        tk.Grid.columnconfigure(self.window, 0, weight=1,  uniform="for sure") 
        tk.Grid.columnconfigure(self.window, 1, weight=1,  uniform="for sure") 
 
        # Add the frames to the main window
        self.frame_options.grid(in_=self.window, row=0, column=0, padx=(20,20), pady=(20,10),stick='NSEW', columnspan=2)
        self.frame_graph1.grid( in_=self.window, row=1, column=0, padx=(20,10), pady=(10,20),stick='NSEW')
        self.frame_graph2.grid( in_=self.window, row=1, column=1, padx=(10,20), pady=(10,20),stick='NSEW')
        
        # Bind focus on every click on a widget: this makes sure the entry widgets become unfocused
        self.frame_graph1.bind("<1>", lambda event: self.frame_graph1.focus_set())
        self.frame_graph2.bind("<1>", lambda event: self.frame_graph2.focus_set())
        self.frame_options.bind("<1>", lambda event: self.frame_options.focus_set()) 
        
        #========================
        # FRAMES FOR THE OPTIONS
        #========================
    
        # Create the subframes in the options frame
        self.subframe_simulations = ttk.LabelFrame(self.frame_options, text="   Simulations  ", **self.dict_awthemes['labelframe2'])
        self.subframe_tabs        = ttk.Frame(self.frame_options)
        self.frame_progress       = ttk.Frame(self.frame_options)

        # Add the progress object
        self.Progress = Progress(self, anchor="fill", length=300)
        self.Progress.move(0,"Please select simulations.")
        
        # Configure the frame
        tk.Grid.rowconfigure(   self.frame_options, 0, weight=1)  
        tk.Grid.rowconfigure(   self.frame_options, 1, weight=0)     
        tk.Grid.columnconfigure(self.frame_options, 0, weight=0)   
        tk.Grid.columnconfigure(self.frame_options, 1, weight=1)  
        
        # Add the labelframes to <frame_options>
        self.subframe_simulations.grid( row=0, column=0, padx=(0,5), pady=(0,0), stick='NSEW')
        self.frame_progress.grid(       row=1, column=0, padx=(0,5), pady=(5,0), stick='NSEW')
        self.subframe_tabs.grid(        row=0, column=1, padx=(5,0), pady=(0,0), stick='NSEW', rowspan=2)
        
        # Fill the subframes with widgets
        self.Simulations = Simulations(self, self.subframe_simulations)
        
        #================================
        # FRAMES FOR THE TABS IN OPTIONS
        #=================================
        
        # Tabbed display for "Plot_linearSpectrum", "Plot_parameterInfluence" and "Plot_linearMap"
        self.tab_plots = ttk.Notebook(self.subframe_tabs, style='species.TNotebook')  
        self.tab_plots.grid(row=0, column=0, sticky="NSEW") 
        tk.Grid.rowconfigure(self.subframe_tabs, 0, weight=1)  
        tk.Grid.columnconfigure(self.subframe_tabs, 0, weight=1)  

        # Create frames for each tab
        self.frame_plotOne   = ttk.Frame(self.tab_plots, padding=(0,0,0,0)) #(10,10,10,10)
        self.frame_plotTwo   = ttk.Frame(self.tab_plots, padding=(0,0,0,0)) #(10,10,10,10)
        self.tab_plots.add(self.frame_plotOne,   text="Linear spectrum") 
        self.tab_plots.add(self.frame_plotTwo,   text="Parameter influence") 
               
        # Fill the subframes with widgets
        self.PlotLinearSpectrum = PlotLinearSpectrum(self, self.frame_plotOne)
        self.PlotParameterInfluence = PlotParameterInfluence(self, self.frame_plotTwo)

        #=============
        # LOAD FIGURE
        #=============
        
        # When the tab is visible, make sure the correct figure if loaded
        self.frame_graph1.bind("<Visibility>", self.load_figure) 
        if True: return
    
#################################################################
#                          METHODS
#################################################################

    def load_figure(self, *args):
        '''  When the tab is visible, make sure the correct figure is loaded. 
        If the simulations were changed, replot the graph. '''

        # Load the canvasses only when the tab is visible, to save time.
        if self.initiated_canvas==None:
            self.CanvasClass = Canvasses(self)
            self.initiated_canvas = True
        
        # Connect the Progress class to the research object
        self.root.Research.Progress = self.Progress
        for experiment in self.root.Research.experiments:
            experiment.Progress = self.Progress
            for simulation in experiment.simulations:
                simulation.Progress = self.Progress
            
        # Make sure we display the correct experiments
        if self.root.input_files != []:
            self.Simulations.display_plottedModesAndExperiments()
        self.show_progress("finished", None) 
        
        # Try to automatically detect the correct knob and key
        Plot = self.PlotParameterInfluence
        if Plot.knob == "-" and len(self.root.Research.experiments)!=0:
            if len(self.Simulations.experiment.variedVariables)!=0:
                variable = self.Simulations.experiment.variedVariables[0]
                Plot.knob = variable.split(":")[0]
                Plot.key = variable.split(": ")[-1]
                Plot.var_parameter.set(Plot.key)
           
        return
    
    #-------------------------------------  
    def show_progress(self, status="start", poppedout_id=None):
        # While plotting show a progress bar turn the cursor into a waiting cursor.
        if status=="start" and poppedout_id == None: 
            self.Progress.move(0,"Start plotting.")
            self.root.config(cursor="watch")
            self.root.update_idletasks() 
        # When finished show the canvas and change the cursor back.
        if status=="finished" and poppedout_id == None: 
            self.Progress.move(100, "No tasks are running.")
            self.root.config(cursor="")
        # When there are no simulations revert to the loading bar
        if status=="nothing" and poppedout_id == None: 
            self.Progress.move(0,"Please select simulations.")
        return
    
    #-------------------------------------  
    def reset_graph(self, axis_id, poppedout_id=None):
        self.CanvasClass.reset_graph(axis_id, poppedout_id)
        
    #-------------------------------------  
    def popout_window(self, axis_id):
        self.CanvasClass.popout_window(axis_id)


