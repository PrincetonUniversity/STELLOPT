#################################################################
#           OPTIONS WINDOW OPENED FROM THE TOOLBAR
#################################################################

# Load modules
import tkinter as tk
from tkinter import ttk

# Load personal modules
from .ModifyStyling import PAD_TITLE, PAD_LABEL, PAD_ENTRY  #@unresolvedimport
from curses.ascii import TAB

#================
# MENU CREATION 
#================
class OptionsWindow:
    """ Window that opens when the "options" button on the CustomToolbar is pressed.
    
    Parameters
    ----------
    root : Tk()
        Root of the tkinter application.
        
    master_class : tab_Convergence, tab_Profiles, tab_Linear, ...
        Must have an attribute called Graph
        
            Graph : class_omegavst, class_potentialvsz, ...
                The plotting class, must have attributes:
                
                    x_name, y_name : str
                    range, label : dict
                    
    axis_id : int
        Multiple canvasses can be linked to one tab, this id defines which canvas to manipulate.
    
    poppedout_id : int
        Multiple popped out windows can exist, this id defines which canvas to manipulate.
    """
    

    def __init__(self, root, master_class, axis_id=0, poppedout_id=None):
        
        # Attach the root so we can carry it into the functions, get the tab for its dimensions
        self.root = root
        
        # Attach the master class: tab_Convergence, tab_Profiles, ...
        self.tab = master_class
        
        # Get the width and height of the root window + title bars and of simply the root window
        self.height = root.winfo_height() # Height of the application minus the application title
        self.width =  root.winfo_width()  # Width of the application  

        # Save which canvas this options window is linked to
        self.axis_id = axis_id
        self.poppedout_id = poppedout_id

        # Prevent indentention of header comment on next lines
        if True: return
        
#==============================
# Open the options window
#==============================
    
    def open_optionsWindow(self):
        
        # Link the plotting class 
        if self.poppedout_id==None:  self.graph = self.tab.Graph[self.axis_id]
        if self.poppedout_id!=None:  self.graph = self.root.graph_poppedOut[self.poppedout_id]
        
        # Create the preferences window and keep it on top
        self.window_options = tk.Toplevel(self.root)
        self.window_options.title("Figure options") 
        self.window_options.attributes('-topmost', 'true')
        
        # Get the height of the screen
        winx = 200; winy = 600;
        if "twin axis" in list(self.graph.layout.keys()):
            if self.graph.layout["twin axis"] == True:
                winy = 750
                
        # Center the new window in the screen
        x = self.width/2  - winx/2
        y = self.height/2 - winy/2
        self.window_options.geometry("+%d+%d" % (x, y))
        
        # Create a tabbed view with the possible settings
        self.tab_header = ttk.Notebook(self.window_options, style='header.TNotebook')
        
        # Add frames to the tab_header which are the tab windows
        self.tab_axis   = ttk.Frame(self.tab_header) 
        self.tab_curves = ttk.Frame(self.tab_header)  
        
        # Add the tabs to the tab header
        self.tab_header.add(self.tab_axis, text='Axis')
        self.tab_header.add(self.tab_curves, text='Curves') 
        self.tab_header.pack(expand=1, fill='both')
        
        # Attach the root so the classes can acces them
        self.tab_axis.root = self.root
        self.tab_curves.root = self.root 
        self.tab_axis.window = self
        self.tab_curves.window = self 
        
        # Fill the tabs with widgets through classes
        self.tabAxes     = tabAxis(self.tab_axis, self.tab, self.graph)
        self.tabCurves   = tabCurves(self.tab_curves, self.tab, self.graph)
        
        # Prevent indentention of header comment on next lines
        if True: return

#==================================
# Show changes on the graph
#==================================
 
    def update_graph(self):
        
        # Draw the canvas and update the GUI
        if self.poppedout_id==None: self.tab.Canvas[self.axis_id].draw_idle()
        if self.poppedout_id!=None: self.root.canvasPoppedOut[self.poppedout_id].draw_idle()
        self.root.update_idletasks()
        
        # Prevent indentention of header comment on next lines
        if True: return

#==================================================================
# Add the tab controling the axis parameters to the options window
#==================================================================
   
class tabAxis:
    
    def __init__(self, window, tab, graph):
        
        # Get data from the GUI
        self.root       = window.root
        self.window     = window.window
        self.tab        = tab 
        self.graph      = graph
        self.data_range = self.graph.range
        self.labels     = self.graph.label
        
        # Variables are on the axis
        self.x_name = self.graph.x_name
        self.y_name = self.graph.y_name 
        if "twin axis" in list(self.graph.layout.keys()):
            if self.graph.layout["twin axis"] == True:
                self.ytwin_name = self.graph.ytwin_name

        # Add some variables
        self.options_scale = sorted(("Linear", "Logaritmic"))
        self.options_units = sorted(("Normalized", "SI units"))
        
        # Create the frame 
        self.tab_axis = ttk.Frame(window)
        self.tab_axis.pack(expand=1, fill=tk.BOTH)
        
        # Configure the frame
        tk.Grid.rowconfigure(   self.tab_axis, 0, weight=0) # x title
        tk.Grid.rowconfigure(   self.tab_axis, 1, weight=0) # x min
        tk.Grid.rowconfigure(   self.tab_axis, 2, weight=0) # x max
        tk.Grid.rowconfigure(   self.tab_axis, 3, weight=0) # x label
        tk.Grid.rowconfigure(   self.tab_axis, 4, weight=0) # x scale
        tk.Grid.rowconfigure(   self.tab_axis, 5, weight=0) # y title
        tk.Grid.rowconfigure(   self.tab_axis, 6, weight=0) # y min
        tk.Grid.rowconfigure(   self.tab_axis, 7, weight=0) # y max
        tk.Grid.rowconfigure(   self.tab_axis, 8, weight=0) # y label
        tk.Grid.rowconfigure(   self.tab_axis, 9, weight=0) # y scale 
        tk.Grid.columnconfigure(self.tab_axis, 0, weight=1, uniform="options") 
        tk.Grid.columnconfigure(self.tab_axis, 0, weight=1, uniform="options") 

        # Add elements to the frame
        self.init_title()
        self.init_xAxis()
        self.init_yAxis() 
        if "twin axis" in list(self.graph.layout.keys()):
            if self.graph.layout["twin axis"] == True:
                self.init_yAxis_twin()
        
    def init_title(self):
        ''' Change the title of the graph. '''
        
        def update_title(event):
            self.graph.label["title"] = self.var_title.get()
            self.graph.ax.set_title(self.graph.label["title"])
            self.window.update_graph()
        
        def update_units(*args):
            if self.var_units.get()=="Normalized":  
                if self.graph.range["units"]!="normalized":
                    self.graph.range["units"]="normalized"
                    self.graph.change_units()
                    self.var_xMin.set(round(self.graph.range["x"][0],2))
                    self.var_xMax.set(round(self.graph.range["x"][1],2))
                    self.var_yMin.set(round(self.graph.range["y"][0],2))
                    self.var_yMax.set(round(self.graph.range["y"][1],2))
                    self.var_xlabel.set(self.labels["x"])
                    self.var_ylabel.set(self.labels["y"])
            if self.var_units.get()=="SI units":    
                if self.graph.range["units"]!="SI units":
                    self.graph.range["units"]="SI units"
                    self.graph.change_units()
                    self.var_xMin.set("{:.2e}".format(self.graph.range["x"][0],2))
                    self.var_xMax.set("{:.2e}".format(self.graph.range["x"][1],2))
                    self.var_yMin.set("{:.2e}".format(self.graph.range["y"][0],2))
                    self.var_yMax.set("{:.2e}".format(self.graph.range["y"][1],2))
                    self.var_xlabel.set(self.labels["x"])
                    self.var_ylabel.set(self.labels["y"])
                    
            
        # Change the graph title
        self.lbl_graph = ttk.Label(self.tab_axis, text="Graph", style='prefTitle.TLabel')
        self.var_title = tk.StringVar(value=self.labels["title"])
        self.lbl_title = ttk.Label(self.tab_axis, text="Title")
        self.ent_title = ttk.Entry(self.tab_axis, textvariable=self.var_title, width=20, style='opt_valueR.TEntry')
        self.ent_title.bind('<Return>', update_title)
        
        # Change the units
        if self.graph.range["units"]!="N.A.":
            self.var_units = tk.StringVar(value=self.options_units[0])
            self.lbl_units = ttk.Label(self.tab_axis, text="Units")
            self.mnu_units = ttk.OptionMenu(self.tab_axis, self.var_units, self.options_units[0], *self.options_units, style='option.TMenubutton')
            self.mnu_units["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'])
            if self.graph.range["units"]=="normalized": self.var_units.set(self.options_units[0])
            if self.graph.range["units"]=="SI units":   self.var_units.set(self.options_units[1])
            self.var_units.trace('w', update_units) # link function to a change of the dropdown options
    
        # Add the widgets to the frame
        self.lbl_graph.grid( row=0, column=0, columnspan=2, **PAD_TITLE)
        self.lbl_title.grid( row=1, column=0, **PAD_LABEL)
        self.ent_title.grid( row=1, column=1, **PAD_ENTRY)
        
        if self.graph.range["units"]!="N.A.":
            self.lbl_units.grid( row=2, column=0, **PAD_LABEL)
            self.mnu_units.grid( row=2, column=1, **PAD_ENTRY)
        
    def init_xAxis(self):
        ''' Change the x-axis of the graph. '''
        
        def update_xAxis(event):
            range_ = [float(self.var_xMin.get()), float(self.var_xMax.get())]
            self.graph.range["x"] = range_
            self.graph.ax.set_xlim(range_)
            self.window.update_graph()
            
        def update_xLabel(event):
            self.graph.label["x"] = self.var_xlabel.get()
            self.graph.ax.set_xlabel(self.graph.label["x"])
            self.window.update_graph()
            
        def update_xscale(*args):
            if self.var_xscale.get()=="Linear":     scale='linear'
            if self.var_xscale.get()=="Logaritmic": scale='log'
            self.graph.range["x_scale"] = scale
            self.graph.ax.set_xscale(scale)
            # Logaritmic axis needs a positive start
            if float(self.var_xMin.get()) <= 0: 
                if float(self.var_xMax.get()) > 10:    self.var_xMin.set(1) 
                elif float(self.var_xMax.get()) > 1:   self.var_xMin.set(0.1) 
                elif float(self.var_xMax.get()) > 0.1: self.var_xMin.set(0.01) 
                range_ = [float(self.var_xMin.get()), float(self.var_xMax.get())]
                self.graph.range["x"] = range_
                self.graph.ax.set_xlim(range_)
            self.window.update_graph()
        
        # Minimum and maximum of the x-axis
        self.lbl_xTitle = ttk.Label(self.tab_axis, text=self.x_name, style='prefTitle.TLabel')
        self.lbl_xMin = ttk.Label(self.tab_axis, text="Minimum")
        self.lbl_xMax = ttk.Label(self.tab_axis, text="Maximum")
        self.var_xMin = tk.StringVar(value=round(self.graph.range["x"][0],2))
        self.var_xMax = tk.StringVar(value=round(self.graph.range["x"][1],2))
        self.ent_xMin = ttk.Entry(self.tab_axis, textvariable=self.var_xMin, width=5, style='opt_valueR.TEntry')
        self.ent_xMax = ttk.Entry(self.tab_axis, textvariable=self.var_xMax, width=5, style='opt_valueR.TEntry')
        self.ent_xMin.bind('<Return>', update_xAxis)
        self.ent_xMax.bind('<Return>', update_xAxis)
        
        # Label for the x-axis
        self.var_xlabel = tk.StringVar(value=self.labels["x"])
        self.lbl_xlabel = ttk.Label(self.tab_axis, text="Label")
        self.ent_xlabel = ttk.Entry(self.tab_axis, textvariable=self.var_xlabel, width=20, style='opt_valueR.TEntry')
        self.ent_xlabel.bind('<Return>', update_xLabel)
                
        # Choice between linear and log scales for the x-axis 
        self.var_xscale = tk.StringVar(value=self.options_scale[0])
        self.lbl_xscale = ttk.Label(self.tab_axis, text="Scale")
        self.mnu_xscale = ttk.OptionMenu(self.tab_axis, self.var_xscale, self.options_scale[0], *self.options_scale, style='option.TMenubutton')
        self.mnu_xscale["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'])
        if self.graph.range["x_scale"]=="linear": self.var_xscale.set(self.options_scale[0])
        if self.graph.range["x_scale"]=="log":    self.var_xscale.set(self.options_scale[1])
        self.var_xscale.trace('w', update_xscale) # link function to a change of the dropdown options
    
        # Add the labels to the frame
        i=3
        self.lbl_xTitle.grid( row=i+0, column=0, columnspan=2, **PAD_TITLE)
        self.lbl_xMin.grid(   row=i+1, column=0, **PAD_LABEL)
        self.ent_xMin.grid(   row=i+1, column=1, **PAD_ENTRY)
        self.lbl_xMax.grid(   row=i+2, column=0, **PAD_LABEL)
        self.ent_xMax.grid(   row=i+2, column=1, **PAD_ENTRY)
        self.lbl_xlabel.grid( row=i+3, column=0, **PAD_LABEL)
        self.ent_xlabel.grid( row=i+3, column=1, **PAD_ENTRY)
        self.lbl_xscale.grid( row=i+4, column=0, **PAD_LABEL)
        self.mnu_xscale.grid( row=i+4, column=1, **PAD_ENTRY)

    def init_yAxis(self):
        ''' Change the y-axis of the graph. '''

        def update_yAxis(event):
            range_ = [float(self.var_yMin.get()), float(self.var_yMax.get())]
            self.graph.range["y"] = range_
            self.graph.ax.set_ylim(range_)
            self.window.update_graph()
            
        def update_yLabel(event):
            self.graph.label["y"] = self.var_ylabel.get()
            self.graph.ax.set_ylabel(self.graph.label["y"])
            self.window.update_graph()
            
        def update_yscale(*args):
            if self.var_yscale.get()=="Linear":     scale='linear'
            if self.var_yscale.get()=="Logaritmic": scale='log'
            self.graph.range["y_scale"] = scale
            self.graph.ax.set_yscale(scale)
            # Logaritmic axis needs a positive start
            if float(self.var_yMin.get()) <= 0: 
                if float(self.var_yMax.get()) > 10:    self.var_yMin.set(1) 
                elif float(self.var_yMax.get()) > 1:   self.var_yMin.set(0.1) 
                elif float(self.var_yMax.get()) > 0.1: self.var_yMin.set(0.01)
                range_ = [float(self.var_yMin.get()), float(self.var_yMax.get())]
                self.graph.range["y"] = range_
                self.graph.ax.set_ylim(range_)
            self.window.update_graph()             
        
        # Minimum and maximum of the x-axis 
        self.lbl_yTitle = ttk.Label(self.tab_axis, text=self.y_name, style='prefTitle.TLabel')
        self.lbl_space = ttk.Label(self.tab_axis, text="    ", style='prefTitle.TLabel')
        self.lbl_yMin = ttk.Label(self.tab_axis, text="Minimum")
        self.lbl_yMax = ttk.Label(self.tab_axis, text="Maximum")
        self.var_yMin = tk.StringVar(value=round(self.graph.range["y"][0],2))
        self.var_yMax = tk.StringVar(value=round(self.graph.range["y"][1],2))
        self.ent_yMin = ttk.Entry(self.tab_axis, textvariable=self.var_yMin, width=5, style='opt_valueR.TEntry')
        self.ent_yMax = ttk.Entry(self.tab_axis, textvariable=self.var_yMax, width=5, style='opt_valueR.TEntry')
        self.ent_yMin.bind('<Return>', update_yAxis)
        self.ent_yMax.bind('<Return>', update_yAxis)
        
        # Label for the y-axis
        self.var_ylabel = tk.StringVar(value=self.labels["y"])
        self.lbl_ylabel = ttk.Label(self.tab_axis, text="Label")
        self.ent_ylabel = ttk.Entry(self.tab_axis, textvariable=self.var_ylabel, width=20, style='opt_valueR.TEntry')
        self.ent_ylabel.bind('<Return>', update_yLabel)
        
        # Choice between linear and log scales for the x-axis 
        self.var_yscale = tk.StringVar(value=self.options_scale[0])
        self.lbl_yscale = ttk.Label(self.tab_axis, text="Scale")
        self.mnu_yscale = ttk.OptionMenu(self.tab_axis, self.var_yscale, self.options_scale[0], *self.options_scale, style='option.TMenubutton')
        self.mnu_yscale["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'])
        if self.graph.range["y_scale"]=="linear": self.var_yscale.set(self.options_scale[0])
        if self.graph.range["y_scale"]=="log":    self.var_yscale.set(self.options_scale[1])
        self.var_yscale.trace('w', update_yscale) # link function to a change of the dropdown options
    
        # Add the labels to the frame
        i=3+5
        self.lbl_yTitle.grid( row=i+0, column=0, columnspan=2, **PAD_TITLE)
        self.lbl_yMin.grid(   row=i+1, column=0, **PAD_LABEL)
        self.ent_yMin.grid(   row=i+1, column=1, **PAD_ENTRY)
        self.lbl_yMax.grid(   row=i+2, column=0, **PAD_LABEL)
        self.ent_yMax.grid(   row=i+2, column=1, **PAD_ENTRY)
        self.lbl_ylabel.grid( row=i+3, column=0, **PAD_LABEL)
        self.ent_ylabel.grid( row=i+3, column=1, **PAD_ENTRY)
        self.lbl_yscale.grid( row=i+4, column=0, **PAD_LABEL)
        self.mnu_yscale.grid( row=i+4, column=1, **PAD_ENTRY)
        self.lbl_space.grid(  row=i+5, column=0, columnspan=2, **PAD_TITLE)
        
        # Prevent indentention of header comment on next lines
        if True: return
        
    def init_yAxis_twin(self):
        ''' Change the y-axis of the graph. '''

        def update_ytwinAxis(event):
            range_ = [float(self.var_ytwinMin.get()), float(self.var_ytwinMax.get())]
            self.graph.range["ytwin"] = range_
            self.graph.ax_twin.set_ylim(range_)
            self.window.update_graph()
            
        def update_ytwinLabel(event):
            self.graph.label["ytwin"] = self.var_ytwinlabel.get()
            self.graph.ax_twin.set_ylabel(self.graph.label["ytwin"])
            self.window.update_graph()
            
        def update_ytwinscale(*args):
            if self.var_ytwinscale.get()=="Linear":     scale='linear'
            if self.var_ytwinscale.get()=="Logaritmic": scale='log'
            self.graph.range["ytwin_scale"] = scale
            self.graph.ax_twin.set_yscale(scale)
            # Logaritmic axis needs a positive start
            if float(self.var_ytwinMin.get()) <= 0: 
                if float(self.var_ytwinMax.get()) > 10:    self.var_ytwinMin.set(1) 
                elif float(self.var_ytwinMax.get()) > 1:   self.var_ytwinMin.set(0.1) 
                elif float(self.var_ytwinMax.get()) > 0.1: self.var_ytwinMin.set(0.01)
                range_ = [float(self.var_ytwinMin.get()), float(self.var_ytwinMax.get())]
                self.graph.range["ytwin"] = range_
                self.graph.ax_twin.set_ylim(range_)
            self.window.update_graph()             
        
        # Minimum and maximum of the twinned y-axis
        self.lbl_ytwinTitle = ttk.Label(self.tab_axis, text=self.ytwin_name, style='prefTitle.TLabel')
        self.lbl_space = ttk.Label(self.tab_axis, text="    ", style='prefTitle.TLabel')
        self.lbl_ytwinMin = ttk.Label(self.tab_axis, text="Minimum")
        self.lbl_ytwinMax = ttk.Label(self.tab_axis, text="Maximum")
        self.var_ytwinMin = tk.StringVar(value=round(self.graph.range["ytwin"][0],2))
        self.var_ytwinMax = tk.StringVar(value=round(self.graph.range["ytwin"][1],2))
        self.ent_ytwinMin = ttk.Entry(self.tab_axis, textvariable=self.var_ytwinMin, width=5, style='opt_valueR.TEntry')
        self.ent_ytwinMax = ttk.Entry(self.tab_axis, textvariable=self.var_ytwinMax, width=5, style='opt_valueR.TEntry')
        self.ent_ytwinMin.bind('<Return>', update_ytwinAxis)
        self.ent_ytwinMax.bind('<Return>', update_ytwinAxis)
        
        # Label for the ytwin-axis
        self.var_ytwinlabel = tk.StringVar(value=self.labels["ytwin"])
        self.lbl_ytwinlabel = ttk.Label(self.tab_axis, text="Label")
        self.ent_ytwinlabel = ttk.Entry(self.tab_axis, textvariable=self.var_ytwinlabel, width=20, style='opt_valueR.TEntry')
        self.ent_ytwinlabel.bind('<Return>', update_ytwinLabel)
        
        # Choice between linear and log scales for the x-axis 
        self.var_ytwinscale = tk.StringVar(value=self.options_scale[0])
        self.lbl_ytwinscale = ttk.Label(self.tab_axis, text="Scale")
        self.mnu_ytwinscale = ttk.OptionMenu(self.tab_axis, self.var_ytwinscale, self.options_scale[0], *self.options_scale, style='option.TMenubutton')
        self.mnu_ytwinscale["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'])
        if self.graph.range["y_scale"]=="linear": self.var_ytwinscale.set(self.options_scale[0])
        if self.graph.range["y_scale"]=="log":    self.var_ytwinscale.set(self.options_scale[1])
        self.var_ytwinscale.trace('w', update_ytwinscale) # link function to a change of the dropdown options
    
        # Add the labels to the frame
        i=3+5+5
        self.lbl_ytwinTitle.grid( row=i+0, column=0, columnspan=2, **PAD_TITLE)
        self.lbl_ytwinMin.grid(   row=i+1, column=0, **PAD_LABEL)
        self.ent_ytwinMin.grid(   row=i+1, column=1, **PAD_ENTRY)
        self.lbl_ytwinMax.grid(   row=i+2, column=0, **PAD_LABEL)
        self.ent_ytwinMax.grid(   row=i+2, column=1, **PAD_ENTRY)
        self.lbl_ytwinlabel.grid( row=i+3, column=0, **PAD_LABEL)
        self.ent_ytwinlabel.grid( row=i+3, column=1, **PAD_ENTRY)
        self.lbl_ytwinscale.grid( row=i+4, column=0, **PAD_LABEL)
        self.mnu_ytwinscale.grid( row=i+4, column=1, **PAD_ENTRY)
        self.lbl_space.grid(      row=i+5, column=0, columnspan=2, **PAD_TITLE)
        
        # Prevent indentention of header comment on next lines
        if True: return

#==================================================================
# Add the tab controling the appearance to the options window
#==================================================================
       
class tabCurves:
    
    def __init__(self, window, tab, graph):
        
        # Get data from the GUI
        self.root       = window.root
        self.tab        = tab
        self.graph      = graph
        self.layout     = graph.layout
        
        # Create the frame 
        self.tab_curves = ttk.Frame(window)
        self.tab_curves.pack(expand=1, fill=tk.BOTH)
        
        # Configure the frame
        tk.Grid.rowconfigure(   self.tab_curves, 0, weight=0) # Appearanc title
        tk.Grid.rowconfigure(   self.tab_curves, 1, weight=0) # Background
        tk.Grid.rowconfigure(   self.tab_curves, 2, weight=0) # Font size
        tk.Grid.rowconfigure(   self.tab_curves, 3, weight=0) # Handle length
        tk.Grid.columnconfigure(self.tab_curves, 0, weight=1, uniform="options") 
        tk.Grid.columnconfigure(self.tab_curves, 0, weight=1, uniform="options") 
        
        # Add elements to the frame
        self.init_appearance()
        return
    
    def init_appearance(self):
        
        # Width of the entry widget
        width=5
                
        # Change background color, font size and handle length
        self.lbl_aTitle = ttk.Label(self.tab_curves, text="Appearance", style='prefTitle.TLabel')
        self.var_bg = tk.StringVar(value=self.root.color['bg'])
        self.var_fs = tk.StringVar(value=self.layout['fontsize'])
        self.var_hl = tk.StringVar(value=self.layout['handlelength'])
        self.lbl_bg = ttk.Label(self.tab_curves, text="Background color")
        self.lbl_fs = ttk.Label(self.tab_curves, text="Font size")
        self.lbl_hl = ttk.Label(self.tab_curves, text="Handle length")
        self.ent_bg = ttk.Entry(self.tab_curves, textvariable=self.var_bg, width=width, style='opt_valueR.TEntry')
        self.ent_fs = ttk.Entry(self.tab_curves, textvariable=self.var_fs, width=width, style='opt_valueR.TEntry')
        self.ent_hl = ttk.Entry(self.tab_curves, textvariable=self.var_hl, width=width, style='opt_valueR.TEntry')
    
        # Add the labels to the frame
        self.lbl_aTitle.grid( row=1, column=0, columnspan=2, **PAD_TITLE)
        self.lbl_bg.grid( row=2, column=0, **PAD_LABEL)
        self.ent_bg.grid( row=2, column=1, **PAD_ENTRY)
        self.lbl_fs.grid( row=3, column=0, **PAD_LABEL)
        self.ent_fs.grid( row=3, column=1, **PAD_ENTRY)
        self.lbl_hl.grid( row=4, column=0, **PAD_LABEL)
        self.ent_hl.grid( row=4, column=1, **PAD_ENTRY)
  

    
    
# print(self.mnu_xscale["menu"].keys())
# ['activebackground', 'activeborderwidth', 'activeforeground', 'background', 'bd', 'bg', 
# 'borderwidth', 'cursor', 'disabledforeground', 'fg', 'font', 'foreground', 'postcommand', 
# 'relief', 'selectcolor', 'takefocus', 'tearoff', 'tearoffcommand', 'title', 'type']

