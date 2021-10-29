


#################################################################
#                   CLASS FOR THE FIRST TAB
#################################################################
''' TAB PROFILES


'''

# Load modules
import os, pathlib
import tkinter as tk
from tkinter import ttk
from PIL import Image, ImageTk 
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec 

# Personal modules 
from stellapy.data import write_profile
from stellapy.plot import plot_profilesVsRho
from stellapy.config import CONFIG
from stellapy.GUI.utils import get_initialDirectory
from stellapy.GUI.graph_tools import PAD_TITLE2, PAD_LABEL2, PAD_ENTRY2  #@UnusedImport 
from stellapy.GUI.graph_tools import CanvasForGraphs
from stellapy.GUI.graph_tools.display_information import display_information

#################################################################
#                   CLASS FOR THE FIRST TAB
#################################################################
class TabProfiles:
 
#################################################################
#                          WIDGETS
#################################################################

    # Initate the tab for "Select simulations"
    def __init__(self, tab2): 


        #======================================================
        # Safe info from the tab and root that can be passed on
        #======================================================
 
        self.window = tab2                       # Direct accesses to the window of tab "Select Simulations" 
        self.root = tab2.root                    # Needed to center windows based on the root screen
        self.dict_awthemes = tab2.root.awthemes  # To make the tk widgets look like the ttk widgets

        #===========
        # VARIABLES
        #===========
        
        self.selected_files = []                  # Saves the selected "raw_profile" or "profile.txt" files
        self.initialdir = get_initialDirectory()  # Directory to start browsing for files
        self.options_radial = ["rho", "s", "r"]   # Options for the radial coordinate
        self.btn_popout = "N.A."                  # Prevent the button to popout windows on the toolbar
        self.btn_reread = "N.A."                  # Prevent the button to popout windows on the toolbar
        self.initiated_canvas = None
        
        # Variables for the functions
        self.x1 = None
        self.a1 = None 
        self.x_col1 = None
        self.n_col1 = None
        self.Te_col1 = None
        self.Ti_col1 = None
        self.x2 = None
        self.a2 = None 
        self.x_col2 = None
        self.n_col2 = None
        self.Te_col2 = None
        self.Ti_col2 = None
        self.plot = "profiles"
        
        #===========
        # FRAMES
        #===========

        # Create the frames and add them to the tab3
        self.frame_graph    = ttk.LabelFrame(self.window, text="   Profiles  ", **self.dict_awthemes['labelframe2'])
        self.frame_options  = ttk.Frame(self.window) 

        # Configure the frames     
        tk.Grid.rowconfigure(   self.window, 0, weight=0) 
        tk.Grid.rowconfigure(   self.window, 1, weight=1) 
        tk.Grid.columnconfigure(self.window, 0, weight=1)
 
        # Add the two frames to the main window
        self.frame_options.grid(in_=self.window, row=0, column=0, padx=(20,20), pady=(20,10),stick='NSEW')
        self.frame_graph.grid(  in_=self.window, row=1, column=0, padx=(20,20), pady=(10,20),stick='NSEW')
        
        # Initiate the widgets in the frame with the options
        self.initiate_frameOptions()      
        
        # Bind focus on every click on a widget: this makes sure the entry widgets become unfocused
        self.frame_graph.bind("<1>", lambda event: self.frame_graph.focus_set())
        self.frame_options.bind("<1>", lambda event: self.frame_options.focus_set()) 
        
        #============
        # POPUP MENU
        #============
        
        self.rightclick_menu = tk.Menu(self.root, tearoff=0)
        self.rightclick_menu.add_command(label='Open file', command = lambda: self.open_file())

        #=============
        # LOAD FIGURE
        #=============
        
        self.frame_graph.bind("<Visibility>", self.load_figure)
        
        # Prevent folding of header comment on next line
        if True: return

#=======================================
# Iniate the frame holding the graph
#=======================================

    def initiate_canvas(self):
        ''' Initiate the canvas. 
        
        This method gets called by stella_GUI.py after initiating the tabs to make sure 
        the attribute self.root.tab_Profiles already exists. 
        
        The toolbar requires the class to have a reset_graph and popout_window function
        The optionswindow requires the Graph objects to have attributes: ax, x_name, x_key, y_name, y_key, range, label
                '''
        
        # Create three frames for thee canvasses
        self.frame_graph1  = ttk.Frame(self.frame_graph) 
        self.frame_graph2  = ttk.Frame(self.frame_graph) 
        self.frame_graph3  = ttk.Frame(self.frame_graph) 
        self.frame_graph4  = ttk.Frame(self.frame_graph) 
        tk.Grid.rowconfigure(   self.frame_graph, 0, weight=1)  
        tk.Grid.columnconfigure(self.frame_graph, 0, weight=1, uniform="columns")
        tk.Grid.columnconfigure(self.frame_graph, 1, weight=1, uniform="columns")
        tk.Grid.columnconfigure(self.frame_graph, 2, weight=1, uniform="columns")
        self.frame_graph1.grid(in_=self.frame_graph, row=0, column=0, stick='NSEW')
        self.frame_graph2.grid(in_=self.frame_graph, row=0, column=1, stick='NSEW')
        self.frame_graph3.grid(in_=self.frame_graph, row=0, column=2, stick='NSEW')
        self.frame_graph4.grid(in_=self.frame_graph, row=0, column=1, stick='NSEW', columnspan=2)
        self.frame_graph4.grid_remove() # Hide this canvas until we plot figure 4
        
        # Create a list to hold the graph objects with attributes ax, x_name, x_key, y_name, y_key, range, label
        self.Graph = [None, None, None, None]   
        self.Canvas = [None, None, None, None]  
        self.initiated_canvas = True
        
        # Create the figures                   
        self.figure1 = plt.figure("profiles1") # Figure for canvas 1   
        self.figure2 = plt.figure("profiles2") # Figure for canvas 2     
        self.figure3 = plt.figure("profiles3") # Figure for canvas 3       
        self.figure4 = plt.figure("profiles4") # A figure that extends across two canvasses   
        self.figure1.set_tight_layout(False)  
        self.figure2.set_tight_layout(False)  
        self.figure3.set_tight_layout(False) 
        self.figure4.set_tight_layout(False)
        
        # Create the canvasses
        CanvasForGraphs(self.root, self.frame_graph1, self.root.tab_Profiles, self.figure1, axis_id=0); self.root.update_idletasks()
        CanvasForGraphs(self.root, self.frame_graph2, self.root.tab_Profiles, self.figure2, axis_id=1); self.root.update_idletasks()
        CanvasForGraphs(self.root, self.frame_graph3, self.root.tab_Profiles, self.figure3, axis_id=2); self.root.update_idletasks()
        CanvasForGraphs(self.root, self.frame_graph4, self.root.tab_Profiles, self.figure4, axis_id=3); self.root.update_idletasks()
        
        # Initiate the graph
        self.initiate_plottingClass()
    
    #--------------------------------
    def initiate_plottingClass(self):

        # Change the canvas color
        self.figure1.patch.set_facecolor(self.root.color['canvas']); self.root.update_idletasks()
        self.figure2.patch.set_facecolor(self.root.color['canvas']); self.root.update_idletasks() 
        self.figure3.patch.set_facecolor(self.root.color['canvas']); self.root.update_idletasks() 
        self.figure4.patch.set_facecolor(self.root.color['canvas']); self.root.update_idletasks() 
        
        # Set the color of the axis
        plt.rcParams['text.color'] = self.root.color['fg']
        plt.rcParams['axes.edgecolor'] = self.root.color['fg']
        plt.rcParams['axes.labelcolor'] = self.root.color['fg']
        plt.rcParams['xtick.color'] = self.root.color['fg']
        plt.rcParams['ytick.color'] = self.root.color['fg']
        plt.rcParams['axes.facecolor'] = self.root.color['canvas']
        plt.rcParams['figure.facecolor'] = self.root.color['canvas']
        
        # Make the axis instance and the required attributes for the options window through the class <graph>
        self.Graph[0] = graph(self.figure1, self.var_plot, self.root.color['canvas'], 0)
        self.Graph[1] = graph(self.figure2, self.var_plot, self.root.color['canvas'], 1)
        self.Graph[2] = graph(self.figure3, self.var_plot, self.root.color['canvas'], 2)
        self.Graph[3] = graph(self.figure4, self.var_plot, self.root.color['canvas'], 3)
        
        # Put the empty figures on the GUI
        self.Canvas[0].draw_idle(); self.root.update_idletasks()
        self.Canvas[1].draw_idle(); self.root.update_idletasks()
        self.Canvas[2].draw_idle(); self.root.update_idletasks()
        self.Canvas[3].draw_idle(); self.root.update_idletasks()
    
        if True: # Prevent indentention of header comment on next lines
            return

#=======================================
# Initiate site bar with extra options
#=======================================

    def initiate_frameOptions(self):
        
        def select_files(*args):
            title           = "Select input file"
            filetypes       = (("profile files","*.dat *.txt"),("all files","*.*")) 
            selected_files  = tk.filedialog.askopenfilenames(initialdir=self.initialdir["Profiles"],title=title,filetypes=filetypes)
            for selected_file in selected_files:
                self.selected_files.append(pathlib.Path(selected_file))
            if len(self.selected_files)>2:
                self.selected_files = self.selected_files[0:2] 
            self.update_filesGUI()
          
        def reset_graph(*args):
            self.reset_graph()
            return 

        # Create the sub_labelframes in options plus a button to apply any changes
        self.subframe_profiles  = ttk.LabelFrame(self.frame_options, text="   Profiles  ", **self.dict_awthemes['labelframe2'])
        self.subframe_profile1  = ttk.LabelFrame(self.frame_options, text="   Profile 1  ", **self.dict_awthemes['labelframe2'])
        self.subframe_profile2  = ttk.LabelFrame(self.frame_options, text="   Profile 2   ", **self.dict_awthemes['labelframe2'])
        self.subframe_plots     = ttk.LabelFrame(self.frame_options, text="   Plots   ", **self.dict_awthemes['labelframe2'])
        self.subframe_buttons   = ttk.Frame(self.frame_options)
        
        # "Select profiles" button
        self.button_selectProfiles = ttk.Button(self.subframe_buttons, text="Select Profiles", width=15)
        self.button_selectProfiles.config(command=select_files)
        
        # "Write profiles" button
        self.button_writeProfiles = ttk.Button(self.subframe_buttons, text="Plot Profiles", width=15)
        self.button_writeProfiles.config(command=reset_graph)
        
        # Configure the frame
        tk.Grid.rowconfigure(   self.frame_options, 0, weight=0)  
        tk.Grid.rowconfigure(   self.frame_options, 1, weight=0)     
        tk.Grid.columnconfigure(self.frame_options, 0, weight=1) 
        tk.Grid.columnconfigure(self.frame_options, 1, weight=0) 
        tk.Grid.columnconfigure(self.frame_options, 2, weight=0) 
        tk.Grid.columnconfigure(self.frame_options, 3, weight=0) 
        
        # Add the three labelframes and the frame for the button to <frame_options>
        self.subframe_profiles.grid( row=0, column=0, padx=(10,10), pady=(0,5),stick='NSEW')
        self.subframe_buttons.grid(  row=1, column=0, padx=(10,10), pady=(5,0),stick='NSEW')
        self.subframe_plots.grid(    row=0, column=1, padx=(10,10), pady=(0,0),stick='NSEW', rowspan=2)
        self.subframe_profile1.grid( row=0, column=2, padx=(10,10), pady=(0,0),stick='NSEW', rowspan=2)
        self.subframe_profile2.grid( row=0, column=3, padx=(10,10), pady=(0,0),stick='NSEW', rowspan=2)
        
        # Add the buttons 
        self.button_selectProfiles.grid(row=0, column=0, padx=(5,10), pady=(15,5),stick='NSEW', ipadx=5, ipady=5)
        self.button_writeProfiles.grid( row=0, column=1, padx=(10,5), pady=(15,5),stick='NSEW', ipadx=5, ipady=5)
        
        # Fill the four subframes with widgets
        self.initiate_frameOptions_profilesFrame()
        self.initiate_frameOptions_profile1Frame()
        self.initiate_frameOptions_profile2Frame()
        self.initiate_frameOptions_plotsFrame()

    #---------------------------------------------
    def initiate_frameOptions_profilesFrame(self):
        
        def remove_profile1(*args):
            if len(self.selected_files)>0:
                self.selected_files.pop(0)
                self.update_filesGUI()
        
        def remove_profile2(*args):
            if len(self.selected_files)>1:
                self.selected_files.pop(1)
                self.update_filesGUI()
                return 
        
        # Add the labels for the two selected profiles
        self.var_prof1 = tk.StringVar(value="Please select a file.")
        self.var_prof2 = tk.StringVar(value="")
        self.lbl_prof1 = ttk.Label(self.subframe_profiles, textvariable=self.var_prof1, padding= (5,5,5,5))
        self.lbl_prof2 = ttk.Label(self.subframe_profiles, textvariable=self.var_prof2, padding= (5,5,5,5))
        
        # Open the file with a right-click
        self.lbl_prof1.bind("<Button-3>", lambda event, i=0: self.popup_menu(event,i))
        self.lbl_prof2.bind("<Button-3>", lambda event, i=1: self.popup_menu(event,i))

        # Add cross icon to remove the profile from the selected profiles
        location_imag = CONFIG['PATHS']['stellapy'] + "/GUI/images/cross.png" 
        icon_cross = ImageTk.PhotoImage(Image.open(location_imag).resize((18, 18), Image.ANTIALIAS))
        self.btn_prof1 = tk.Button(self.subframe_profiles, image=icon_cross, relief=tk.FLAT, **self.root.awthemes['crossButton'])
        self.btn_prof2 = tk.Button(self.subframe_profiles, image=icon_cross, relief=tk.FLAT, **self.root.awthemes['crossButton'])
        self.btn_prof1.config(command=remove_profile1)
        self.btn_prof2.config(command=remove_profile2)
        self.btn_prof1.icon = icon_cross
        self.btn_prof2.icon = icon_cross
        
        # Configure the frame
        tk.Grid.columnconfigure(self.subframe_profiles, 0, weight=1) 
        tk.Grid.columnconfigure(self.subframe_profiles, 1, weight=0) 
        
        # Add the options to the frame
        self.lbl_prof1.grid(row=0, column=0, **PAD_LABEL2)
        self.lbl_prof2.grid(row=1, column=0, **PAD_LABEL2)
        self.btn_prof1.grid(row=0, column=1)
        self.btn_prof2.grid(row=1, column=1)
        return 
    
    #--------------------------------------------
    def initiate_frameOptions_profile1Frame(self):
        
        # Text options
        font = ("Courier New", 11); width=10; 
        
        # Choose the radial coordinate
        self.var_radial1 = tk.StringVar(value=self.options_radial[0])
        self.mnu_radial1 = ttk.OptionMenu(self.subframe_profile1, self.var_radial1, self.options_radial[0], *self.options_radial, style='option.TMenubutton')
        self.mnu_radial1["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
        
        # Select the columns for rho; n; Te and Ti
        self.var_rho1 = tk.StringVar(value="column x")
        self.var_n1   = tk.StringVar(value="column x")
        self.var_Te1  = tk.StringVar(value="column x")
        self.var_Ti1  = tk.StringVar(value="column x")
        self.var_a1   = tk.StringVar(value="N.A.")
        self.lbl_rho1 = ttk.Label(self.subframe_profile1, text="Radial coordinate: ")
        self.lbl_n1   = ttk.Label(self.subframe_profile1, text="Density:")
        self.lbl_Te1  = ttk.Label(self.subframe_profile1, text="Electron temperature:   ")
        self.lbl_Ti1  = ttk.Label(self.subframe_profile1, text="Ion temperature:")
        self.lbl_a1   = ttk.Label(self.subframe_profile1, text="Minor radius:")
        self.ent_rho1 = ttk.Entry(self.subframe_profile1, textvariable=self.var_rho1, font=font,  style='opt_valueCBold.TEntry', width=width)
        self.ent_n1   = ttk.Entry(self.subframe_profile1, textvariable=self.var_n1, font=font,  style='opt_valueCBold.TEntry', width=width)
        self.ent_Te1  = ttk.Entry(self.subframe_profile1, textvariable=self.var_Te1, font=font,  style='opt_valueCBold.TEntry', width=width)
        self.ent_Ti1  = ttk.Entry(self.subframe_profile1, textvariable=self.var_Ti1, font=font,  style='opt_valueCBold.TEntry', width=width)
        self.ent_a1    = ttk.Entry(self.subframe_profile1, textvariable=self.var_a1, font=font,  style='opt_valueCBold.TEntry', width=width)
        
        # Configure the frame
        tk.Grid.columnconfigure(self.subframe_profile1, 0, weight=1) 
        tk.Grid.columnconfigure(self.subframe_profile1, 1, weight=0) 
        tk.Grid.columnconfigure(self.subframe_profile1, 2, weight=0) 
        
        # Place the widgets in the frame
        self.mnu_radial1.grid(row=0, column=1, sticky="w", padx=(0,10))
        self.lbl_rho1.grid(row=0, column=0, **PAD_LABEL2)
        self.ent_rho1.grid(row=0, column=2, **PAD_ENTRY2)
        self.lbl_n1.grid(  row=1, column=0, **PAD_LABEL2, columnspan=2)
        self.ent_n1.grid(  row=1, column=2, **PAD_ENTRY2)
        self.lbl_Te1.grid( row=2, column=0, **PAD_LABEL2, columnspan=2)
        self.ent_Te1.grid( row=2, column=2, **PAD_ENTRY2)
        self.lbl_Ti1.grid( row=3, column=0, **PAD_LABEL2, columnspan=2)
        self.ent_Ti1.grid( row=3, column=2, **PAD_ENTRY2)
        self.lbl_a1.grid(  row=4, column=0, **PAD_LABEL2, columnspan=2)
        self.ent_a1.grid(  row=4, column=2, **PAD_ENTRY2)
        return 
    
    #-------------------------------------------
    def initiate_frameOptions_profile2Frame(self):
        # Text options
        font = ("Courier New", 11); width=10; 
        
        # Choose the radial coordinate
        self.var_radial2 = tk.StringVar(value=self.options_radial[0])
        self.mnu_radial2 = ttk.OptionMenu(self.subframe_profile2, self.var_radial2, self.options_radial[0], *self.options_radial, style='option.TMenubutton')
        self.mnu_radial2["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
        
        # Select the columns for rho; n; Te and Ti
        self.var_rho2 = tk.StringVar(value="column x")
        self.var_n2   = tk.StringVar(value="column x")
        self.var_Te2  = tk.StringVar(value="column x")
        self.var_Ti2  = tk.StringVar(value="column x")
        self.var_a2   = tk.StringVar(value="N.A.")
        self.lbl_rho2 = ttk.Label(self.subframe_profile2, text="Radial coordinate: ")
        self.lbl_n2   = ttk.Label(self.subframe_profile2, text="Density")
        self.lbl_Te2  = ttk.Label(self.subframe_profile2, text="Electron temperature:   ")
        self.lbl_Ti2  = ttk.Label(self.subframe_profile2, text="Ion temperature:")
        self.lbl_a2   = ttk.Label(self.subframe_profile2, text="Minor radius:")
        self.ent_rho2 = ttk.Entry(self.subframe_profile2, textvariable=self.var_rho2, font=font,  style='opt_valueCBold.TEntry', width=width)
        self.ent_n2   = ttk.Entry(self.subframe_profile2, textvariable=self.var_n2, font=font,  style='opt_valueCBold.TEntry', width=width)
        self.ent_Te2  = ttk.Entry(self.subframe_profile2, textvariable=self.var_Te2, font=font,  style='opt_valueCBold.TEntry', width=width)
        self.ent_Ti2  = ttk.Entry(self.subframe_profile2, textvariable=self.var_Ti2, font=font,  style='opt_valueCBold.TEntry', width=width)
        self.ent_a2    = ttk.Entry(self.subframe_profile2, textvariable=self.var_a2, font=font,  style='opt_valueCBold.TEntry', width=width)
        
        # Configure the frame
        tk.Grid.columnconfigure(self.subframe_profile2, 0, weight=1) 
        tk.Grid.columnconfigure(self.subframe_profile2, 0, weight=0) 
        
        # Place the widgets in the frame
        self.mnu_radial2.grid(row=0, column=1, sticky="w", padx=(0,10))
        self.lbl_rho2.grid(row=0, column=0, **PAD_LABEL2)
        self.ent_rho2.grid(row=0, column=2, **PAD_ENTRY2)
        self.lbl_n2.grid(  row=1, column=0, **PAD_LABEL2, columnspan=2)
        self.ent_n2.grid(  row=1, column=2, **PAD_ENTRY2)
        self.lbl_Te2.grid( row=2, column=0, **PAD_LABEL2, columnspan=2)
        self.ent_Te2.grid( row=2, column=2, **PAD_ENTRY2)
        self.lbl_Ti2.grid( row=3, column=0, **PAD_LABEL2, columnspan=2)
        self.ent_Ti2.grid( row=3, column=2, **PAD_ENTRY2)
        self.lbl_a2.grid(  row=4, column=0, **PAD_LABEL2, columnspan=2)
        self.ent_a2.grid(  row=4, column=2, **PAD_ENTRY2)
        return       
    
    #-------------------------------------------       
    def initiate_frameOptions_plotsFrame(self):
        
        def change_plot(*args):
            self.change_plot()
            return
        
        # Set the default value of the plot out of range so we don't start with a plot
        self.var_plot  = tk.IntVar(value=1)
        
        # Add the possible plots in two sections: "Time convergence" and "Space convergence"
        self.rbn_prof    = ttk.Radiobutton(self.subframe_plots, text='  Profiles')
        self.rbn_grad    = ttk.Radiobutton(self.subframe_plots, text='  Gradients')
        self.rbn_norm    = ttk.Radiobutton(self.subframe_plots, text='  Normalized gradients')
        self.rbn_ratio   = ttk.Radiobutton(self.subframe_plots, text='  Normalized gradients and temperature ratio')
        self.rbn_compare = ttk.Radiobutton(self.subframe_plots, text='  Compare the profiles and files')
        
        # Add the values and commands to the radiobuttons
        self.rbn_prof.config(     value=1, variable=self.var_plot, command=change_plot)
        self.rbn_grad.config(     value=2, variable=self.var_plot, command=change_plot)
        self.rbn_norm.config(     value=3, variable=self.var_plot, command=change_plot)
        self.rbn_ratio.config(    value=4, variable=self.var_plot, command=change_plot)
        self.rbn_compare.config(  value=5, variable=self.var_plot, command=change_plot)
        
        # Configure the frame
        tk.Grid.rowconfigure(self.subframe_plots, 0, weight=1) 
        tk.Grid.rowconfigure(self.subframe_plots, 1, weight=1) 
        tk.Grid.rowconfigure(self.subframe_plots, 2, weight=1) 
        tk.Grid.rowconfigure(self.subframe_plots, 3, weight=1)  
        tk.Grid.rowconfigure(self.subframe_plots, 4, weight=1)
        tk.Grid.columnconfigure(self.subframe_plots, 0, weight=1) 
        
        # Add the options to the frame
        self.rbn_prof.grid(  row=0, column=0, **PAD_LABEL2)
        self.rbn_grad.grid(  row=1, column=0, **PAD_LABEL2) 
        self.rbn_norm.grid(  row=2, column=0, **PAD_LABEL2) 
        self.rbn_ratio.grid( row=3, column=0, **PAD_LABEL2) 
        self.rbn_compare.grid(  row=4, column=0, **PAD_LABEL2) 

        if True: # Prevent folding of header comment on next line
            return

        
#################################################################
#                          METHODS
#################################################################

    def load_figure(self, *args):
        '''  When the tab is visible, make sure the correct figure is loaded. 
        If the simulations were changed, replot the graph. '''
        
        # If the canvas isn't loaded, load them
        if self.initiated_canvas==None:
            self.initiate_canvas()
        
        # Load figure
        plt.figure("profiles")
        return

    #------------------------------
    def reset_graph(self, **kwargs):
        self.read_variablesGUI()
        self.write_profiles()
        self.plot_profiles()
    
    #------------------------------   
    def popup_menu(self, event, i):
        ''' Rightclick event linked to each self._widgets[folder]['label'][i]
        
        When rightclicking on a file, it saves and opents it.
        '''

        # Bind the button information to self
        self.rightclick_menu_file = self.selected_files[i]
        
        # Pop-up the menu at the mouse location
        try:
            self.rightclick_menu.tk_popup(event.x_root, event.y_root)
        finally:
            self.rightclick_menu.grab_release()
            
        return 
    
    #-------------------
    def open_file(self):
        ''' Command executed when choosing an option of self.rightclick_menu. '''
        
        # If the file exists, open it with the texteditor written in the CONFIG file
        if os.path.isfile(self.rightclick_menu_file):           
            if os.system(CONFIG['GENERAL']['texteditor'] + ' ' + self.rightclick_menu_file) != 0: # Return code is non-zero if it failed to execute
                display_information(self.root,'Error','The text editor '+CONFIG['GENERAL']['texteditor'] + ' is not installed, chose another one.')
        return 
    
    #-------------------
    def change_plot(self):
        
        # Set the title of the labelframe of the graph
        if self.var_plot.get() in [1]: 
            self.frame_graph.config(text="   Profiles  ")
        if self.var_plot.get() in [2]: 
            self.frame_graph.config(text="   Gradients  ")
        if self.var_plot.get() in [3]: 
            self.frame_graph.config(text="   Normalized gradients  ")
        if self.var_plot.get() in [4]: 
            self.frame_graph.config(text="   Normalized gradients and temperature ratio ")
        if self.var_plot.get() in [5]: 
            self.frame_graph.config(text="   Compare the profiles and files  ")
        
        # Update which plot is plotted
        if self.var_plot.get() == 1:     
            self.plot = "profiles"   
        if self.var_plot.get() == 2:     
            self.plot = "gradients"   
        if self.var_plot.get() == 3:     
            self.plot = "normalized gradients"   
        if self.var_plot.get() == 4:     
            self.plot = "a/Ln"   
        if self.var_plot.get() == 5:     
            self.plot = "relative comparison"
        
        self.plot_profiles()
        
    #-------------------------
    def update_filesGUI(self):
        if len(self.selected_files)==0:
            self.var_prof1.set("")
            self.var_prof2.set("")
        if len(self.selected_files)==1: 
            self.var_prof1.set(self.selected_files[0].name)
            self.var_prof2.set("")
        if len(self.selected_files)==2: 
            self.var_prof1.set(self.selected_files[0].name)
            self.var_prof2.set(self.selected_files[1].name)
        return 
    
    # --------------------------
    def read_variablesGUI(self):
        
        # Read the variables on the GUI
        for source_file in self.selected_files:  
            if self.selected_files.index(source_file)==0:
                x = self.var_radial1.get()
                a = float(self.var_a1.get()) if x=='r' else None
                x_col  = int(self.var_rho1.get().split("column ")[-1]) if self.var_rho1.get() != "column x" else None
                n_col  = int(self.var_n1.get(  ).split("column ")[-1]) if self.var_n1.get(  ) != "column x" else None
                Te_col = int(self.var_Te1.get( ).split("column ")[-1]) if self.var_Te1.get( ) != "column x" else None
                Ti_col = int(self.var_Ti1.get( ).split("column ")[-1]) if self.var_Ti1.get( ) != "column x" else None

            if self.selected_files.index(source_file)==1:
                x = self.var_radial2.get()
                a = float(self.var_a2.get()) if x=='r' else None
                x_col  = int(self.var_rho2.get().split("column ")[-1]) if self.var_rho2.get() != "column x" else None
                n_col  = int(self.var_n2.get(  ).split("column ")[-1]) if self.var_n2.get(  ) != "column x" else None
                Te_col = int(self.var_Te2.get( ).split("column ")[-1]) if self.var_Te2.get( ) != "column x" else None
                Ti_col = int(self.var_Ti2.get( ).split("column ")[-1]) if self.var_Ti2.get( ) != "column x" else None
        
            # Save the variables
            if self.selected_files.index(source_file)==0:
                self.x1 = x 
                self.a1 = a 
                self.x_col1 = x_col
                self.n_col1 = n_col
                self.Te_col1 = Te_col
                self.Ti_col1 = Ti_col        
            if self.selected_files.index(source_file)==1:
                self.x2 = x 
                self.a2 = a 
                self.x_col2 = x_col
                self.n_col2 = n_col
                self.Te_col2 = Te_col
                self.Ti_col2 = Ti_col
        return
    
    #------------------------   
    def write_profiles(self):
        
        for raw_file in self.selected_files: 

            if self.selected_files.index(raw_file)==0:
                x = self.x1;           a = self.a1
                x_col = self.x_col1;   n_col = self.n_col1
                Ti_col = self.Ti_col1; Te_col = self.Te_col1
            if self.selected_files.index(raw_file)==1:
                x = self.x2;           a = self.a2
                x_col = self.x_col2;   n_col = self.n_col2
                Ti_col = self.Ti_col2; Te_col = self.Te_col2
            
            x_col, dens_col, Ti_col, Te_col = write_profile(\
                # Data of the profiles
                raw_file=raw_file, \
                # Parameter to measure distance along the cross-section of the plasma
                x=x, \
                a=a, \
                # Columns of the rho, density and temperature data
                x_col=x_col, dens_col=n_col, Ti_col=Ti_col, Te_col=Te_col, \
                # The positions throughout the cross-section to calculate the profiles
                rho_values=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], \
                # Number of digits of the calculated profiles and gradients
                digits=5)
            
            if self.selected_files.index(raw_file)==0:
                self.var_rho1.set("column "+str(x_col))
                self.var_n1.set("column "+str(dens_col))
                self.var_Te1.set("column "+str(Te_col))
                self.var_Ti1.set("column "+str(Ti_col))
            if self.selected_files.index(raw_file)==1:
                self.var_rho2.set("column "+str(x_col))
                self.var_n2.set("column "+str(dens_col))
                self.var_Te2.set("column "+str(Te_col))
                self.var_Ti2.set("column "+str(Ti_col))

        return
    
    #--------------------------
    def plot_profiles(self):
        
        # Load the figure again
        plt.figure("profiles")
        
        # We can have canvasses 1,2,3 on display 
        if self.var_plot.get() in [1,2,3,4] and self.Graph[1].plotted==False:  
            
            # Load the correct canvas
            self.Graph[0].plotted = True;  self.Graph[1].plotted = True
            self.Graph[2].plotted = True;  self.Graph[3].plotted = False
            self.frame_graph4.grid_remove()
            self.frame_graph2.grid()
            self.frame_graph3.grid()

        # We can have canvasses 1,4 on display
        if self.var_plot.get() in [5] and self.Graph[3].plotted==False: 
            
            # Load the correct canvas
            self.Graph[0].plotted = True;  self.Graph[3].plotted = True
            self.Graph[1].plotted = False; self.Graph[2].plotted = False
            self.frame_graph2.grid_remove()
            self.frame_graph3.grid_remove()
            self.frame_graph4.grid()
            
        # The first three canvasses don't have a twin axis
        if self.var_plot.get() in [1,2,3]:
            self.Graph[0].ax_twin.tick_params(axis='y', labelcolor=self.root.color['canvas'], colors=self.root.color['canvas'])
            self.Graph[1].ax_twin.tick_params(axis='y', labelcolor=self.root.color['canvas'], colors=self.root.color['canvas'])
            self.Graph[2].ax_twin.tick_params(axis='y', labelcolor=self.root.color['canvas'], colors=self.root.color['canvas'])
                    
        # Clear the axis
        self.Graph[0].ax.clear(); self.Graph[0].ax_twin.clear()
        self.Graph[1].ax.clear(); self.Graph[1].ax_twin.clear()
        self.Graph[2].ax.clear(); self.Graph[2].ax_twin.clear()
        self.Graph[3].ax.clear(); self.Graph[3].ax_twin.clear()
        
        # Load the correct axis
        if self.var_plot.get() in [1,2,3,4]:
            ax1 = self.Graph[0].ax; ax1_twin = self.Graph[0].ax_twin
            ax2 = self.Graph[1].ax; ax2_twin = self.Graph[1].ax_twin
            ax3 = self.Graph[2].ax; ax3_twin = self.Graph[2].ax_twin
        if self.var_plot.get() in [5]:
            ax1 = self.Graph[0].ax; ax1_twin = self.Graph[0].ax_twin; ax3 = None
            ax2 = self.Graph[3].ax; ax2_twin = self.Graph[3].ax_twin; ax3_twin = None
        
        # Add the plot
        plot_profilesVsRho(\
            # Get the data 
            raw_files=self.selected_files, \
            # Define the radial coordinate                
            x1=self.x1, a1=self.a1, \
            x2=self.x2, a2=self.a2, \
            # State in which col1umns the data can be found
            x_col1=self.x_col1, dens_col1=self.n_col1, Ti_col1=self.Ti_col1, Te_col1=self.Te_col1, \
            x_col2=self.x_col2, dens_col2=self.n_col2, Ti_col2=self.Ti_col2, Te_col2=self.Te_col2, \
            # If it is called from the GUI, the three axis exists
            ax1 = ax1, ax1_twin = ax1_twin,\
            ax2 = ax2, ax2_twin = ax2_twin,\
            ax3 = ax3, ax3_twin = ax3_twin,\
            plot = self.plot)
        
        # Update the <graph> classes so the option window can use them
        self.Graph[0].update()
        self.Graph[1].update()
        self.Graph[2].update()
        self.Graph[3].update()
        
        # Update screen
        self.Canvas[0].draw_idle(); self.root.update_idletasks()
        self.Canvas[1].draw_idle(); self.root.update_idletasks()
        self.Canvas[2].draw_idle(); self.root.update_idletasks()
        self.Canvas[3].draw_idle(); self.root.update_idletasks()
        
        # Prevent indentation
        if True: return
    
#################################################################
#                          Classes
#################################################################

# Make the axis instance and the required attributes for the options window
# Either the plotting function is a class or we manually give it some class attributes like here
class graph:
    
    def __init__(self, figure, var_plot, canvas_color, i):
        plt.figure(figure.number)
        grid_specifications = gridspec.GridSpec(1, 1)
        grid_specifications.update(top=0.9, left=0.15, right=0.85, bottom=0.15)
        self.ax = plt.subplot(grid_specifications[0])
        self.ax_twin = self.ax.twinx()
        self.ax_twin.tick_params(axis='y', colors=canvas_color, labelcolor=canvas_color)
        self.plotted = False
        self.range = {"units" : "N.A.", "x_scale" : "linear", "y_scale" : "linear"}    
        self.label = {}
        self.layout = {"fontsize" : "N.A.", 'handlelength' : "N.A.", 'twin axis' : False}
        self.x_key  = "x"
        self.y_key  = "y"
        self.ytwin_key  = "ytwin"
        self.x_name = "Radial Coordinate"
        self.y_name = ""
        self.ytwin_name = "Electron to ion temperature ratio"
        if i==0: self.y_name = "Density"
        if i==1: self.y_name = "Electron temperature"
        if i==2: self.y_name = "Ion temperature"
        
        # Save the graph and id
        self.var_plot = var_plot
        self.id = i
        return 
    
    #---------------
    def update(self):
        
        # Get the ranges, labels and titles
        self.range["x"]     = self.ax.get_xlim()
        self.range["y"]     = self.ax.get_ylim()
        self.label["x"]     = self.ax.get_xlabel()
        self.label["y"]     = self.ax.get_ylabel()
        self.label["title"] = self.ax.get_title()
        self.range["ytwin"] = self.ax_twin.get_ylim()
        self.label["ytwin"] = self.ax_twin.get_ylabel()
        
        # Define what is on the axis for the options window
        if self.var_plot.get()==0: 
            self.layout['twin axis'] = False
            if self.id==0: self.y_name = "Density"
            if self.id==1: self.y_name = "Electron temperature"
            if self.id==2: self.y_name = "Ion temperature"
        if self.var_plot.get()==1: 
            self.layout['twin axis'] = False
            if self.id==0: self.y_name = "Density gradient"
            if self.id==1: self.y_name = "Electron temperature gradient"
            if self.id==2: self.y_name = "Ion temperature gradient"
        if self.var_plot.get()==2: 
            self.layout['twin axis'] = False
            if self.id==0: self.y_name = "Normalized density gradient"
            if self.id==1: self.y_name = "Normalized electron temperature gradient"
            if self.id==2: self.y_name = "Normalized ion temperature gradient"
        if self.var_plot.get()==3: 
            self.layout['twin axis'] = True
            if self.id==0: self.y_name = "Normalized gradients"
            if self.id==1: self.y_name = "Normalized gradients"
            if self.id==2: self.y_name = "Normalized gradients"
        if self.var_plot.get()==4: 
            self.layout['twin axis'] = True
            if self.id==0: self.y_name = "Normalized gradients"
            if self.id==3: self.y_name = "Relative normalized gradients" 
        return






