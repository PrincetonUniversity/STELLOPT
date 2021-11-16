
import numpy as np
import tkinter as tk
from tkinter import ttk  
from stellapy.GUI.graph_tools.display_information import display_information
        
#################################################################
#                  CLASS TO CHOOSE SIMULATIONS
#################################################################
   
class Simulations: 
    
    def __init__(self, parent, frame):
        ''' This options frame controls which experiments/simulations are plotted. '''
        
        # Make the parents available
        self.tab = parent
        self.root = parent.root
        
        # Attributes
        self.experiment = None
        self.simulation = None
        self.experiment_id = "All experiments"
        self.simulation_id = "All simulations"
        self.options_kx = ["0.0"]
        self.options_ky = ["0.0", "100"]
        self.options_experiment = ["First plot data"]
        self.options_simulation = ["First plot data"]
        self.options_simulationsids = []
        
        # Variables
        self.var_experiment = tk.StringVar(value=self.options_experiment[0])
        self.var_simulation = tk.StringVar(value=self.options_simulation[0])
        self.var_kxmin = tk.StringVar(value=self.options_kx[0])
        self.var_kxmax = tk.StringVar(value=self.options_kx[-1])
        self.var_kymin = tk.StringVar(value=self.options_ky[0])
        self.var_kymax = tk.StringVar(value=self.options_ky[-1])
        
        # Configure the frame
        tk.Grid.columnconfigure(frame, 0, weight=1) 
        
        # Create subframes
        self.frame1 = ttk.Frame(frame)
        self.frame2 = ttk.Frame(frame)   
        self.frame1.grid(row=0, column=0, padx=(0,0), pady=(0,0), stick='NSEW')
        self.frame2.grid(row=2, column=0, padx=(0,0), pady=(0,0), stick='NSEW')
        
        # Choice between the experiments
        self.lbl_experiment = ttk.Label(self.frame1, text="Experiment: "); width=20
        self.mnu_experiment = ttk.OptionMenu(self.frame1, self.var_experiment, self.options_experiment[0], *self.options_experiment, style='option.TMenubutton')
        self.mnu_experiment["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
        self.mnu_experiment.config(width=width)
        self.var_experiment.trace('w', self.change_plottedExperiment) # link function to a change of the dropdown options
        
        # Choice between the simulations
        self.lbl_simulation = ttk.Label(self.frame1, text="Simulation: ")
        self.mnu_simulation = ttk.OptionMenu(self.frame1, self.var_simulation, self.options_simulation[0], *self.options_simulation, style='option.TMenubutton')
        self.mnu_simulation["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
        self.mnu_simulation.config(width=width)
        self.var_simulation.trace('w', self.change_plottedSimulation) # link function to a change of the dropdown options
        
        # Show the range of the kx and the k values
        self.lbl_kxmin = ttk.Label(self.frame2, text="kx min: ", style='opt_sign.TLabel'); width=4
        self.lbl_kxmax = ttk.Label(self.frame2, text=" kx max: ", style='opt_sign.TLabel')
        self.lbl_kymin = ttk.Label(self.frame2, text="ky min: ", style='opt_sign.TLabel')
        self.lbl_kymax = ttk.Label(self.frame2, text=" ky max: ", style='opt_sign.TLabel')
        self.mnu_kxmin = ttk.OptionMenu(self.frame2, self.var_kxmin, self.options_kx[0], *self.options_kx, style='option.TMenubutton')
        self.mnu_kxmin["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
        self.mnu_kxmin.config(width=width)
        self.mnu_kxmax = ttk.OptionMenu(self.frame2, self.var_kxmax, self.options_kx[0], *self.options_kx, style='option.TMenubutton')
        self.mnu_kxmax["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
        self.mnu_kxmax.config(width=width)
        self.mnu_kymin = ttk.OptionMenu(self.frame2, self.var_kymin, self.options_ky[0], *self.options_ky, style='option.TMenubutton')
        self.mnu_kymin["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
        self.mnu_kymin.config(width=width)
        self.mnu_kymax = ttk.OptionMenu(self.frame2, self.var_kymax, self.options_ky[1], *self.options_ky, style='option.TMenubutton')
        self.mnu_kymax["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
        self.mnu_kymax.config(width=width)
        self.var_kxmin.trace('w', self.change_selectedModes) # link function to a change of the dropdown options
        self.var_kxmax.trace('w', self.change_selectedModes) # link function to a change of the dropdown options
        self.var_kymin.trace('w', self.change_selectedModes) # link function to a change of the dropdown options
        self.var_kymax.trace('w', self.change_selectedModes) # link function to a change of the dropdown options

        # Configure the subframe 
        tk.Grid.columnconfigure(self.frame1, 1, weight=1) 
        tk.Grid.columnconfigure(self.frame2, 1, weight=1, uniform="2") 
        tk.Grid.columnconfigure(self.frame2, 3, weight=1, uniform="2") 
        
        # Place the widgets in the frame
        self.lbl_experiment.grid(row=0, column=0, stick='NSEW', padx=(0,0), pady=(2,2))
        self.mnu_experiment.grid(row=0, column=1, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=2, ipady=1)
        self.lbl_simulation.grid(row=1, column=0, stick='NSEW', padx=(0,0), pady=(2,2))
        self.mnu_simulation.grid(row=1, column=1, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=2, ipady=0)
        self.lbl_kxmin.grid(row=0, column=0, stick='W',    padx=(0,0), pady=(2,2), ipadx=1, ipady=0)
        self.mnu_kxmin.grid(row=0, column=1, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=1, ipady=0)
        self.lbl_kxmax.grid(row=0, column=2, stick='W',    padx=(0,0), pady=(2,2), ipadx=1, ipady=0)
        self.mnu_kxmax.grid(row=0, column=3, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=1, ipady=0)
        self.lbl_kymin.grid(row=1, column=0, stick='W',    padx=(0,0), pady=(2,2), ipadx=1, ipady=0)
        self.mnu_kymin.grid(row=1, column=1, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=1, ipady=0)
        self.lbl_kymax.grid(row=1, column=2, stick='W',    padx=(0,0), pady=(2,2), ipadx=1, ipady=0)
        self.mnu_kymax.grid(row=1, column=3, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=1, ipady=0)
        if True: return


#################################################################
#                          METHODS
#################################################################

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
            
        # Reset the graph classes
        self.tab.Graph[0].load_defaults()
        self.tab.Graph[1].load_defaults()
        self.tab.Graph[2].load_defaults()
        
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
    
    #--------------------------------------
    def change_selectedModes(self, *args): 
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
        self.tab.Graph[0].kx = [ kxmin, kxmax ] 
        self.tab.Graph[0].ky = [ kymin, kymax ] 
        self.tab.Graph[1].kx = [ kxmin, kxmax ] 
        self.tab.Graph[1].ky = [ kymin, kymax ] 
        self.tab.Graph[2].kx = [ kxmin, kxmax ] 
        self.tab.Graph[2].ky = [ kymin, kymax ] 
        
        # Make sure one of the ranges is a float
        if kxmin == kxmax:
            self.tab.Graph[0].kx = kxmin
            self.tab.Graph[1].kx = kxmin
            self.tab.Graph[2].kx = kxmin
        elif kymin == kymax:
            self.tab.Graph[0].ky = kymin
            self.tab.Graph[1].ky = kymin
            self.tab.Graph[2].ky = kymin
        return

    #-------------------------------------------
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
            
        # Reset the options in the menu
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
        
        # Reset the options in the menu
        self.mnu_simulation['menu'].delete(0, 'end')
        for simulation in self.options_simulations:
            self.mnu_simulation['menu'].add_command(label=simulation, command=tk._setit(self.var_simulation, simulation))
        
        # Update the options for the kx and ky ranges
        self.options_kx = self.experiment.total_kx
        self.options_ky = self.experiment.total_ky
        if self.experiment_id == "All experiments":
            for experiment in self.root.Research.experiments:
                self.options_kx = sorted(list(set(self.options_kx + experiment.total_kx)))
                self.options_ky = sorted(list(set(self.options_ky + experiment.total_ky)))
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

        # Update the options for the k ranges for the extraction frames!
        for Plot in [self.tab.PlotLinearMap, self.tab.PlotParameterInfluence]:
            mnu_kvalue = Plot.Extraction.mnu_kvalue
            options_k  = Plot.options_k
            var_kvalue = Plot.var_specificKvalue
            if isinstance(self.tab.Graph[0].kx, list): options_k  = self.experiment.total_ky
            if isinstance(self.tab.Graph[0].ky, list): options_k  = self.experiment.total_ky
            mnu_kvalue['menu'].delete(0, 'end')
            for k in options_k:
                mnu_kvalue['menu'].add_command(label=round(k,2), command=tk._setit(var_kvalue, k))
            
        # Update the GUI variables
        if isinstance(self.tab.Graph[0].kx, list): 
            kx_min = round(np.max([self.tab.Graph[0].kx[0], self.options_kx[0]]),2)
            self.var_kxmin.set(kx_min)
            kx_max = round(np.min([self.tab.Graph[0].kx[1], self.options_kx[-1]]),2)
            self.var_kxmax.set(kx_max) 
        if isinstance(self.tab.Graph[0].ky, list): 
            ky_min = round(np.max([self.tab.Graph[0].ky[0], self.options_ky[0]]),2)
            self.var_kymin.set(ky_min)
            ky_max = round(np.min([self.tab.Graph[0].ky[1], self.options_ky[-1]]),2) 
            self.var_kymax.set(ky_max) 

        return

