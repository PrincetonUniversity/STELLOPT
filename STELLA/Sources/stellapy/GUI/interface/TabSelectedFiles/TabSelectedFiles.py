
#################################################################
#                   CLASS FOR THE FIRST TAB
#################################################################
''' SELECT SIMULATIONS ON TAB 1

Initiate the frames on the first tab "Simulations".
- Initiate the frames for "simulations", "progress" and "input parameters"
- Fill the frames through the classes [InputParameters, Simulations, Progress]


Absolute path of class TabSelectedFiles:
----------------------------------------
root.tab_Simulations
    
    
Attributes: widgets
-------------------
frame_simulation, frame_progress, frame_inputs, frame_wavenumber:  ttk.Frame(frame_left/frame_middle/frame_right)  
frame_geometry, frame_resolution, frame_species, frame_allParam:   ttk.Frame(frame_left/frame_middle/frame_right)  


Attributes: classes
-------------------
class_inputParameters: InputParameters
class_simulations:     Simulations
class_progress:        Progress
'''

# Load modules
import tkinter as tk
from tkinter import ttk 
from stellapy.GUI.graph_tools import Progress
from .TabSelectedFiles_Research import Research 
from .TabSelectedFiles_Simulations import Simulations 
from .TabSelectedFiles_InputParameters import InputParameters

#################################################################
#               CLASS FOR THE FIRST TAB
#################################################################
class TabSelectedFiles:
    '''
    Initiate the frames on the first tab "Simulations".
    - Initiate the frames for "simulations", "progress" and "input parameters"
    - Fill the frames through the classes [InputParameters, Simulations, Progress]

    
    Attributes: objects
    -------------------
    class_inputParameters: InputParameters
    class_simulations:     Simulations
    class_progress:        Progress
    

    Parent widgets
    --------------
    root:               Tk()
    tabheader:          ttk.Notebook(root)
    tab1:               ttk.Frame(tabheader)
    

    Local variables
    ---------------  
    frame_left: ttk.Frame(tab1)
        Holds [frame_simulation, frame_progress]
    frame_middle: ttk.Frame(tab1)
        Holds [frame_inputs, frame_wavenumber, frame_geometry, frame_resolution]
    frame_right: ttk.Frame(tab1)
        Holds [frame_species, frame_allParam]
   
   
    Attributes: widgets
    ------------------- 
    frame_simulation, frame_progress, frame_inputs, frame_wavenumber:  ttk.Frame(frame_left/frame_middle/frame_right)  
    frame_geometry, frame_resolution, frame_species, frame_allParam:   ttk.Frame(frame_left/frame_middle/frame_right)  
    '''
 
#################################################################
#                          WIDGETS
#################################################################

    # Initate the tab for "Select simulations"
    def __init__(self,tab1): 

        #==================================
        # Safe info from the tab and root 
        #==================================
           
        self.root = tab1.root                   # Needed to center windows based on the root screen
        dict_awthemes = tab1.root.awthemes      # To make the tk widgets look like the ttk widgets        

        #======================
        # Create the subframes
        #======================

        # Create three subframes to order the widgets in three columns
        frame_left   = ttk.Frame(tab1)
        frame_middle = ttk.Frame(tab1)
        frame_right  = ttk.Frame(tab1)
        
        # Configure the window
        tk.Grid.rowconfigure(   tab1, 0, weight=1) 
        tk.Grid.columnconfigure(tab1, 0, weight=5, uniform="columns") 
        tk.Grid.columnconfigure(tab1, 1, weight=3, uniform="columns")
        tk.Grid.columnconfigure(tab1, 2, weight=3, uniform="columns")
        
        # Attach the three subframes to the window
        frame_left.grid(    in_=tab1, row=0, column=0, sticky='NSEW')
        frame_middle.grid(  in_=tab1, row=0, column=1, sticky='NSEW')
        frame_right.grid(   in_=tab1, row=0, column=2, sticky='NSEW')

        # Configure the left frame
        tk.Grid.rowconfigure(   frame_left,  0, weight=1, uniform="columns") # frame_simulation
        tk.Grid.rowconfigure(   frame_left,  1, weight=1, uniform="columns") # frame_research
        tk.Grid.rowconfigure(   frame_left,  2, weight=0) # frame_progress
        tk.Grid.columnconfigure(frame_left,  0, weight=1) 

        # Configure the middle frame
        tk.Grid.rowconfigure(   frame_middle, 0, weight=1) # frame_inputs        
        tk.Grid.rowconfigure(   frame_middle, 1, weight=1) # frame_wavenumber
        tk.Grid.rowconfigure(   frame_middle, 2, weight=1) # frame_geometry
        tk.Grid.rowconfigure(   frame_middle, 3, weight=1) # frame_resolution
        tk.Grid.columnconfigure(frame_middle, 0, weight=1) 

        # Configure the right frame
        tk.Grid.rowconfigure(   frame_right, 0, weight=0) # frame_species            
        tk.Grid.rowconfigure(   frame_right, 1, weight=1) # frame_allParam
        tk.Grid.columnconfigure(frame_right, 0, weight=1)  
        
        #=========================
        # Create the subsubframes
        #=========================

        # Show the progress that is being made by the GUI when it runs update_GUI()
        self.frame_progress = ttk.Frame(frame_left)
        self.class_progress = Progress(self, anchor="fill", length=300)
        self.class_progress.move(0,"Please select simulations.")  

        # Show the input parameters of the selected files  
        text_inputs     = "    Simulation   "
        text_wavenumber = "    Wavenumbers   "
        text_geometry   = "    Geometry   "
        text_resolution = "    Resolution    "
        text_allParam   = "    All parameters    "
        self.frame_inputs            = ttk.LabelFrame(frame_middle, text=text_inputs,      **dict_awthemes['labelframe2'])
        self.frame_wavenumber        = ttk.LabelFrame(frame_middle, text=text_wavenumber,  **dict_awthemes['labelframe2'])
        self.frame_geometry          = ttk.LabelFrame(frame_middle, text=text_geometry,    **dict_awthemes['labelframe2'])
        self.frame_resolution        = ttk.LabelFrame(frame_middle, text=text_resolution,  **dict_awthemes['labelframe2'])
        self.frame_species           = ttk.Frame(frame_right)
        self.frame_allParam          = ttk.LabelFrame(frame_right,  text=text_allParam,    **dict_awthemes['labelframe2'])
        self.class_inputParameters   = InputParameters(self)  # Class to fill the <Input parameters> frame

        # Show the selected simulations and allows to add or remove simulations  
        text = "    Simulations     "
        self.frame_simulation        = ttk.LabelFrame(frame_left, text=text, **dict_awthemes['labelframe'])
        self.class_simulations       = Simulations(self)    # Class to fill the <Simulations> frame

        # Show the selected simulations and allows to add or remove simulations  
        text = "    Research     "
        self.frame_research          = ttk.LabelFrame(frame_left, text=text, **dict_awthemes['labelframe'])
        self.class_research          = Research(self)    # Class to fill the <Simulations> frame
        
        # Attach the 9 subsubframes to the 3 subframes
        self.frame_simulation.grid(  in_=frame_left,   row=0, column=0, padx=(20,5), pady=(30,5),sticky='NSEW')
        self.frame_research.grid(    in_=frame_left,   row=1, column=0, padx=(20,5), pady=(5,5), sticky='NSEW')
        self.frame_progress.grid(    in_=frame_left,   row=2, column=0, padx=(20,5), pady=(5,20),sticky='NSEW')
        self.frame_inputs.grid(      in_=frame_middle, row=0, column=0, padx=(5,5),  pady=(30,5),sticky='NSEW')
        self.frame_wavenumber.grid(  in_=frame_middle, row=1, column=0, padx=(5,5),  pady=(5,5), sticky='NSEW')
        self.frame_geometry.grid(    in_=frame_middle, row=2, column=0, padx=(5,5),  pady=(5,5), sticky='NSEW')
        self.frame_resolution.grid(  in_=frame_middle, row=3, column=0, padx=(5,5),  pady=(5,20),sticky='NSEW')
        self.frame_species.grid(     in_=frame_right,  row=0, column=0, padx=(5,20), pady=(5,5), sticky='NSEW')
        self.frame_allParam.grid(    in_=frame_right,  row=1, column=0, padx=(5,20), pady=(5,20),sticky='NSEW')


        # When the tab is visible,  Connect the Progress class to the research object
        def load_progressBar(*args):
            self.root.Research.Progress = self.class_progress
            for experiment in self.root.Research.experiments:
                experiment.Progress = self.class_progress
                for simulation in experiment.simulations:
                    simulation.Progress = self.class_progress
        frame_left.bind("<Visibility>", load_progressBar)






