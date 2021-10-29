
import tkinter as tk
from tkinter import ttk 
from stellapy.GUI.graph_tools import PAD_TITLE2, PAD_LABEL2
from stellapy.plot.plot_frequencyVsModes import plot_frequencyVsModes

#################################################################
#            CLASS TO PLOT FREQUNECY VERSUS MODES
#################################################################
 
class PlotLinearSpectrum: 
    
    def __init__(self, parent, frame):
        ''' This class manages the <plot_frequencyVsModes> function. '''
        
        # Make the parent available
        self.tab = parent
        self.root = parent.root
        self.frame = frame
        
        # Keep track of what is plotted
        self.replot = False
        self.plotted_inputFiles = None
        
        # Attributes for the plot
        self.y_quantity = "gamma"
        self.plotted_modes="unstable"
        self.lineardata = "average"
        self.show_error = False
        self.interpolate = False
        self.maxima = False
        
        # Variables to switch the plotting options
        self.var_yQuantity = tk.IntVar(value=2)
        self.var_plotErrorBars = tk.IntVar(value=0)
        self.var_maxima = tk.IntVar(value=0)
        self.var_interpolate = tk.IntVar(value=0)
        
        # Create a subframe for each set of plotting options
        self.subframe_yQuantity = ttk.Frame(frame)
        self.subframe_options = ttk.Frame(frame)
        
        # Configure the main frame to hold the subframes
        tk.Grid.rowconfigure(   frame, 0, weight=1)      
        tk.Grid.rowconfigure(   frame, 1, weight=0)      
        tk.Grid.rowconfigure(   frame, 2, weight=1)      
        tk.Grid.rowconfigure(   frame, 3, weight=1)       
        tk.Grid.columnconfigure(frame, 0, weight=0) # subframe_yQuantity
        tk.Grid.columnconfigure(frame, 1, weight=0) # subframe_options

        # Add the subframes to the main frame
        self.subframe_yQuantity.grid(row=1, column=0, padx=25, pady=(0,0), stick='NSEW')
        self.subframe_options.grid(  row=1, column=1, padx=25, pady=(0,0), stick='NSEW')
        
        # Fill the frames with their options
        self.YQuantity = self.Frame_YQuantity(self, self.subframe_yQuantity)
        self.Options   = self.Frame_Options(  self, self.subframe_options)
        return
   
    class Frame_YQuantity:
        
        def __init__(self, parent, frame):
            
            # Make the parents available
            self.plot = parent 
            self.tab = parent.tab
        
            # Configure the frame <self.subframe_yQuantity>
            tk.Grid.columnconfigure(frame, 0, weight=1) 
            
            # Create the widgets
            self.lbl_spectrum  = ttk.Label(frame, text="Quantity on the y-axis", style='prefTitle.TLabel')
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
            self.tab.Graph[0].load_defaults()
            self.tab.Graph[0].update_quantitiesAndKeys()
            
            # Extract what needs to be plotted on the y-axis
            y_quantity = self.plot.var_yQuantity.get()
            if y_quantity==1: self.plot.y_quantity = "omega";       header = "Frequency versus modes"
            if y_quantity==2: self.plot.y_quantity = "gamma";       header = "Growthrate versus modes" 
            if y_quantity==3: self.plot.y_quantity = "gamma/ky**2"; header = "Growthrate/ky**2 versus modes"
            
            # Change the header of the frame and plot the graph
            self.tab.frame_graph1.config(text="   "+header+"  ")
            self.plot.replot = True; self.plot.plot_graph(None)
            return

    class Frame_Options:
        
        def __init__(self, parent, frame):

            # Make the parents available
            self.plot = parent 
            self.tab = parent.tab
              
            # Create the widgets
            self.lbl_options = ttk.Label(frame, text="Options", style='prefTitle.TLabel')
            self.chk_error  = ttk.Checkbutton(frame, text=" Include error bars")
            self.chk_maxima = ttk.Checkbutton(frame, text=" Find the maxima")
            self.chk_interp = ttk.Checkbutton(frame, text=" Interpolate and find maxima")
            
            # Add the commands
            self.chk_error.config(variable=self.plot.var_plotErrorBars, command=self.update_error)
            self.chk_maxima.config(variable=self.plot.var_maxima,  command=self.update_maxima)
            self.chk_interp.config(variable=self.plot.var_interpolate,  command=self.update_interp)
            
            # Configure the frame
            tk.Grid.columnconfigure(frame, 0, weight=1) 
            
            # Add the options to the frame
            self.lbl_options.grid( row=0, column=0, **PAD_TITLE2)
            self.chk_error.grid(   row=1, column=0, **PAD_LABEL2)
            self.chk_maxima.grid(  row=2, column=0, **PAD_LABEL2)
            self.chk_interp.grid(  row=3, column=0, **PAD_LABEL2)
            return 
        
        #-----------------------
        def update_error(self): 
            if self.plot.var_plotErrorBars.get()==0: self.plot.show_error = False
            if self.plot.var_plotErrorBars.get()==1: self.plot.show_error = True
            self.plot.replot = True; self.plot.plot_graph(None)
            return 
        
        #-----------------------
        def update_maxima(self): 
            if self.plot.var_maxima.get()==0: self.plot.maxima = False
            if self.plot.var_maxima.get()==1: self.plot.maxima = True
            self.plot.replot = True; self.plot.plot_graph(None)
            if True: return 
                   
        #-----------------------
        def update_interp(self): 
            if self.plot.var_interpolate.get()==0: self.plot.interpolate = False
            if self.plot.var_interpolate.get()==1: self.plot.interpolate = True
            self.plot.replot = True; self.plot.plot_graph(None)
            if True: return 

#################################################################
#                          METHODS
#################################################################

    def plot_graph(self, poppedout_id=None):
        
        # Perhaps the button was clicked on the popped out window
        if poppedout_id == None: Graph = self.tab.Graph[0]
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
                
                # Plot the figure 
                plot_frequencyVsModes(\
                    # Specify which simulations to plot
                        research=self.root.Research,\
                        experiment_id=self.tab.Simulations.experiment_id,\
                        simulation_id=self.tab.Simulations.simulation_id,\
                    # Specify data range
                        y_quantity=self.y_quantity,\
                        x_range=Graph.range["x"],\
                        y_range=Graph.range["y"],\
                    # Details of the modes
                        kx_range=Graph.kx, \
                        ky_range=Graph.ky, \
                        k_value=None, \
                        k_delta=None, \
                        plotted_modes=self.plotted_modes,\
                        lineardata=self.lineardata,\
                    # Labels
                        x_label=Graph.label["x"],\
                        y_label=Graph.label["y"],\
                    # For the GUI the figure object already exists 
                        show_figure = False,\
                        ax=Graph.ax, \
                        Progress=self.tab.Progress,\
                        root=self.root,\
                    # Toggles
                        units=Graph.range["units"], \
                        show_error = self.show_error,\
                        maxima = self.maxima,\
                        interpolate = self.interpolate)
                    
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
        if poppedout_id==None: self.tab.Canvas[0].draw_idle()
        if poppedout_id!=None: self.root.canvasPoppedOut[poppedout_id].draw_idle()
        if True: return
