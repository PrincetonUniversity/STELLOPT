
# Load modules
import tkinter as tk
from tkinter import ttk 
from stellapy.GUI.graph_tools import PAD_TITLE2, PAD_LABEL2, PAD_ENTRY2, PAD_ENTRYR
from stellapy.plot.plot_frequencyVsParameter import plot_frequencyVsParameter

#################################################################
#          CLASS TO PLOT THE INFLUENCE OF A PARAMETER
#################################################################
   
class PlotParameterInfluence: 
    
    def __init__(self, parent, frame):
        
        # Make the parents available
        self.tab = parent
        self.root = parent.root
        self.frame = frame
        
        # Keep track of what is plotted
        self.replot = False
        self.plotted_inputFiles = None
        
        # Attributes for the plot
        self.y_quantity = "gamma"
        self.k_value = "max"
        self.show_error = False
        self.knob = "-"
        self.key  = "-"

        # Option lists
        self.options_par = ["rho", "tiprim", "teprim", "fprim", "delta t", "nmu", "nvgrid", "nzed", "nzgrid"]
        self.options_k    = ["-"]
        
        # Variables to switch the plotting options
        self.var_yQuantity = tk.IntVar(value=2)
        self.var_parameter = tk.StringVar(value=self.options_par[0])
        self.var_specificKvalue = tk.StringVar(value=self.options_k[0])
        self.var_extraction = tk.IntVar(value=1)
        self.var_plotErrorBars = tk.IntVar(value=0)
        
        # Create a subframe for each set of plotting options
        self.subframe_yQuantity = ttk.Frame(frame)
        self.subframe_parameter = ttk.Frame(frame)
        self.subframe_extraction = ttk.Frame(frame)
        self.subframe_options = ttk.Frame(frame)
        
        # Configure the main frame to hold the subframes
        tk.Grid.rowconfigure(   frame, 0, weight=1)      
        tk.Grid.rowconfigure(   frame, 1, weight=0)      
        tk.Grid.rowconfigure(   frame, 2, weight=1)      
        tk.Grid.rowconfigure(   frame, 3, weight=1)      
        tk.Grid.columnconfigure(frame, 0, weight=0) # subframe_yQuantity
        tk.Grid.columnconfigure(frame, 1, weight=0) # subframe_parameter
        tk.Grid.columnconfigure(frame, 2, weight=0) # subframe_extraction
        tk.Grid.columnconfigure(frame, 3, weight=0) # subframe_options

        # Add the subframes to the main frame
        self.subframe_yQuantity.grid(  row=1, column=0, padx=25, pady=(0,0), stick='NSEW', rowspan=2)
        self.subframe_parameter.grid(  row=1, column=1, padx=25, pady=(0,0), stick='NSEW')
        self.subframe_extraction.grid( row=1, column=2, padx=25, pady=(0,0), stick='NSEW')
        self.subframe_options.grid(    row=1, column=3, padx=25, pady=(0,0), stick='NSEW', rowspan=2)
        
        # Fill the frames with their options
        self.YQuantity  = self.Frame_YQuantity(self, self.subframe_yQuantity)
        self.Parameter  = self.Frame_Parameter(self, self.subframe_parameter)
        self.Extraction = self.Frame_Extraction(self, self.subframe_extraction)
        self.Options    = self.Frame_Options(self, self.subframe_options)
        
        # When the tab is visible, make sure its canvas and frame are loaded
        self.frame.bind("<Visibility>", self.load_figure) 
        return     
    
    #------------------- 
    def load_figure(self, event):
        
        # Try to automatically detect the correct knob and key
        if self.knob == "-" and len(self.root.Research.experiments)!=0:
            if len(self.tab.Simulations.experiment.variedVariables)!=0:
                variable = self.tab.Simulations.experiment.variedVariables[0]
                self.knob = variable.split(":")[0]
                self.key = variable.split(": ")[-1]
                self.var_parameter.set(self.key)
        
    class Frame_YQuantity:
        
        def __init__(self, parent, frame):
            
            # Make the parents available
            self.plot = parent 
            self.tab = parent.tab 
        
            # Configure the frame <self.subframe_yQuantity>
            tk.Grid.columnconfigure(frame, 0, weight=1) 
            
            # Create the widgets
            self.lbl_spectrum  = ttk.Label(frame, text="Quantity on the Y-axis", style='prefTitle.TLabel')
            self.rbn_freq      = ttk.Radiobutton(frame, text='  Frequency')
            self.rbn_growth    = ttk.Radiobutton(frame, text='  Growth rate')
            self.rbn_growthky2 = ttk.Radiobutton(frame, text='  Growth rate on ky**2')
    
            # Add the values and commands to the radiobuttons
            self.rbn_freq.config(      value=1, variable=self.plot.var_yQuantity, command=self.change_yQuantity)
            self.rbn_growth.config(    value=2, variable=self.plot.var_yQuantity, command=self.change_yQuantity)
            self.rbn_growthky2.config( value=3, variable=self.plot.var_yQuantity, command=self.change_yQuantity)
            
            # Add the options to the frame
            self.lbl_spectrum.grid( row=1,  column=0, **PAD_TITLE2)
            self.rbn_freq.grid(     row=2,  column=0, **PAD_LABEL2)
            self.rbn_growth.grid(   row=3,  column=0, **PAD_LABEL2)
            self.rbn_growthky2.grid(row=4,  column=0, **PAD_LABEL2)
            return 
        
        #---------------------
        def change_yQuantity(self, *args):
            
            # When changing the plots, reset the axis of the graph
            self.tab.Graph[1].load_defaults()
            self.tab.Graph[1].update_quantitiesAndKeys()
            
            # Extract what needs to be plotted on the y-axis
            y_quantity = self.plot.var_yQuantity.get()
            if y_quantity==1: self.plot.y_quantity = "omega";       header = "Frequency versus parameter"
            if y_quantity==2: self.plot.y_quantity = "gamma";       header = "Growthrate versus parameter" 
            if y_quantity==3: self.plot.y_quantity = "gamma/ky**2"; header = "Growthrate/ky**2 versus parameter"
            
            # Change the header of the frame and plot the graph
            self.tab.frame_graph2.config(text="   "+header+"  ")
            self.plot.replot = True; self.plot.plot_graph(None)
            return
        
    class Frame_Parameter:
        
        def __init__(self, parent, frame):

            # Make the parents available
            self.plot = parent 
            self.tab = parent.tab
            self.root = parent.tab.root
        
            # Configure the frame <self.subframe_parameter>
            tk.Grid.columnconfigure(frame, 0, weight=1) 
            
            # Choose the parameter for the x-axis 
            self.lbl_prm = ttk.Label(frame, text="Parameters", style='prefTitle.TLabel'); width=10
            self.mnu_par = ttk.OptionMenu(frame, self.plot.var_parameter, self.plot.options_par[0], *self.plot.options_par, style='option.TMenubutton')
            self.mnu_par["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
            self.mnu_par.config(width=width)
            self.plot.var_parameter.trace('w', self.update_parameter) # link function to a change of the dropdown options

            # Add the widgets to the frame
            self.lbl_prm.grid(row=0, column=0, **PAD_TITLE2)
            self.mnu_par.grid(row=1, column=0, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=1, ipady=0)

        def update_parameter(self, *args):
            
            # Load the defaults when we touch this
            self.tab.Graph[1].load_defaults()
            
            # Get the parameter
            var_par = self.plot.var_parameter.get() 
            
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
            if var_par=="nzed":     knob = "zgrid_parameters";       key="nzed"
            if var_par=="nzgrid":   knob = "zgrid_parameters";       key="nzgrid"
            if var_par=="-":        knob = "-";                      key="-"
            
            # Save the knobs and keys
            self.plot.knob = knob;     self.plot.key = key 
            return
        
    class Frame_Extraction:
        
        def __init__(self, parent, frame):

            # Make the parents available
            self.plot = parent 
            self.tab = parent.tab
            self.root = parent.tab.root
            
            # Configure the frame <self.subframe_extraction>
            tk.Grid.columnconfigure(frame, 0, weight=0) 
            tk.Grid.columnconfigure(frame, 1, weight=0) 
            tk.Grid.columnconfigure(frame, 2, weight=1) 
            
            # Choose which growth rate we extract: the maximum one or at a specific wavenumber
            self.lbl_extraction = ttk.Label(frame, text="Extraction", style='prefTitle.TLabel')
            self.rbn_gammamax = ttk.Radiobutton(frame, text='  Plot the most unstable mode')
            self.rbn_gammak   = ttk.Radiobutton(frame, text='  Plot the mode at ky =')
      
            # Add the commands to the radio buttons
            self.rbn_gammamax.config( value=1, variable=self.plot.var_extraction, command=self.change_extraction)
            self.rbn_gammak.config(   value=2, variable=self.plot.var_extraction, command=self.change_extraction)
      
            # Drop down menu to choose the k-value
            self.mnu_kvalue = ttk.OptionMenu(frame, self.plot.var_specificKvalue, self.plot.options_k[0], *self.plot.options_k, style='option.TMenubutton'); width=4
            self.mnu_kvalue["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
            self.mnu_kvalue.config(width=width)
            self.plot.var_specificKvalue.trace('w', self.update_selectedMode) # link function to a change of the dropdown options
        
            # Add the widgets to the frame
            self.lbl_extraction.grid(row=0,  column=0, **PAD_TITLE2, columnspan=3)
            self.rbn_gammamax.grid(  row=1,  column=0, **PAD_LABEL2, columnspan=3)
            self.rbn_gammak.grid(    row=2,  column=0, **PAD_LABEL2)
            self.mnu_kvalue.grid(    row=2,  column=1, **PAD_LABEL2) 
            return 
        
        #-----------------------------------------
        def update_selectedMode(self, *args):
            self.plot.k_value = float(self.plot.var_specificKvalue.get())
            self.plot.replot=True; self.plot.plot_graph(None) 
                
        #-----------------------------------------
        def change_extraction(self):
            if self.plot.var_extraction.get()==1: self.plot.k_value = "max"
            if self.plot.var_extraction.get()==2: self.plot.k_value = float(self.plot.var_specificKvalue.get())
            self.plot.replot=True; self.plot.plot_graph(None) 
            return

    class Frame_Options:
        
        def __init__(self, parent, frame):

            # Make the parents available
            self.plot = parent 
            self.tab = parent.tab
              
            # Create the widgets
            self.lbl_options = ttk.Label(frame, text="Options", style='prefTitle.TLabel')
            self.chk_error = ttk.Checkbutton(frame, text=" Include error bars")
            
            # Add the commands
            self.chk_error.config(variable=self.plot.var_plotErrorBars, command=self.update_error)
            
            # Configure the frame
            tk.Grid.columnconfigure(frame, 0, weight=1) 
            
            # Add the options to the frame
            self.lbl_options.grid(row=0, column=0, **PAD_TITLE2)
            self.chk_error.grid(  row=1, column=0, **PAD_LABEL2)

        def update_error(self, *args): 
            if self.plot.var_plotErrorBars.get()==0: self.plot.show_error = False
            if self.plot.var_plotErrorBars.get()==1: self.plot.show_error = True
            self.plot.replot=True; self.plot.plot_graph(None) 
    
    def plot_graph(self, poppedout_id=None):
        
        # Perhaps the button was clicked on the popped out window
        if poppedout_id == None: Graph = self.tab.Graph[1]
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
                self.tab.show_progress("start", poppedout_id) 
                
                # Plot the graph
                plot_frequencyVsParameter(\
                    # Specify which simulations to plot
                        research=self.root.Research,\
                        experiment_id=self.tab.Simulations.experiment_id,\
                        simulation_id=self.tab.Simulations.simulation_id,\
                        parameter_knob=self.knob,\
                        parameter_key=self.key,\
                    # Specify data range
                        y_quantity=self.y_quantity,\
                        x_range=Graph.range["x"],\
                        y_range=Graph.range["y"],\
                    # Details of the modes
                        kx_range=Graph.kx, \
                        ky_range=Graph.ky, \
                        k_value=self.k_value, \
                        lineardata="average",\
                    # Labels
                        x_label=Graph.label["x"],\
                        y_label=Graph.label["y"],\
                        title=None,\
                    # For the GUI the figure object already exists 
                        show_error = self.show_error,\
                        show_figure = False,\
                        ax=Graph.ax, \
                        Progress=self.tab.Progress,\
                        root=self.root)
                    
                # Update the <graph> class so the option window can use them
                Graph.update_rangesAndLabels()
            
                # When finished show the canvas and change the cursor back.
                self.tab.show_progress("finished", poppedout_id) 
    
        # If no simulations have been selected, clear the current figure
        if self.root.input_files == [] or self.root.Research.input_files == []:
            Graph.ax.clear()
            Graph.load_defaults()
            self.tab.show_progress("nothing", poppedout_id) 
        
        # Update screen
        if poppedout_id==None: self.tab.Canvas[1].draw_idle()
        if poppedout_id!=None: self.root.canvasPoppedOut[poppedout_id].draw_idle()
        if True: return
   