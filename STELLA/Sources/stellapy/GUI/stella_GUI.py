#!/usr/bin/env python3

#================================================================
########################  STELLA GUI ############################
#================================================================
''' STELLA GUI

Script to run the STELLA GUI. 
- Generate the root window
- Fill the tabs through the classes [TabSelectedFiles, TabProfiles, TabConvergence, TabGraphs, TabVideos]
- Check whether the configuration file is correct with check_pathsToCodeAndSimulations()
- The color scheme is linked to the theme and set with ModifyStyling()
- Create the preference window "..." with PreferenceWindow()
- Run root.mainloop() to make the root window listen to any button presses or keyboard inputs.


Widgets
-------
root:                 Tk()
tab_header:           ttk.Notebook(root)
tab1, ..., tab5:      ttk.Frame(tab_header)


Classes
-------
root.tab_Simulations: TabSelectedFiles(tab1)   
root.tab_Profiles:    TabProfiles(tab2)       
root.tab_Convergence: TabConvergence(tab3)   
root.tab_Graphs:      TabGraphs(tab4)         
root.tab_Videos:      TabVideos(tab5)
'''

#----------- TO DO ---------------
# TODO: Be able to change the experiment labels and marker labels in the GUI
# TODO: Read the time and z-axis more efficiently
# TODO: When saving the research save some more things to speed up the reading (save dimensions, unique time axis)
# TODO: Add a plot of the surface averages of the fluxes
# TODO: Add a movie of the fluxes versus (kx,ky) 

# Load modules
import psutil
import os, sys
import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt

# Customize the warnings
import warnings
def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return '%s:%s: %s: %s\n' % (filename, lineno, category.__name__, message)
warnings.formatwarning = warning_on_one_line

# Tell python where to find the personal modules, and load them
path_stellaGUI = os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0]
sys.path.append(path_stellaGUI)
from stellapy.config import CONFIG, check_pathsToCodeAndSimulations, turnOnVerboseWrapper_configurationFile, turnOffVerboseWrapper_configurationFile #@unresolvedimport
from stellapy.GUI.graph_tools import ModifyStyling 
from stellapy.GUI.interface.PreferencesWindow import PreferenceWindow
from stellapy.GUI.interface.TabSelectedFiles import TabSelectedFiles 
from stellapy.GUI.interface.TabProfiles import TabProfiles 
from stellapy.GUI.interface.TabConvergence import TabConvergence
from stellapy.GUI.interface.TabLinear import TabLinear  
from stellapy.GUI.interface.TabNonlinear import TabNonlinear  
from stellapy.GUI.interface.TabVideos import TabVideos 

#======================
# ROOT WINDOW CREATION 
#======================

# Application: main window which is called root
root = tk.Tk(className="Stellapy")
root.geometry()
root.title("Stellapy: graphical environment for stella")

#============================================
# CHECK WHETHER THE PATHS ARE SET CORRECTLY
#============================================

# Will look at the paths set in the configuration file
check_pathsToCodeAndSimulations(root)

# If the files are set we can add the icon to the application.
divider = '\\' if (os.name == 'nt') else '/'
if os.getcwd().split(divider)[1] == 'home':
    root.iconphoto(False, tk.PhotoImage(file=CONFIG['PATHS']['stellapy']+"GUI/images/stellarator_long.png"))

#=========================
# DEFAULT "SAVE" LOCATION 
#=========================

import matplotlib as mpl 
mpl.rcParams["savefig.directory"] = CONFIG["PATHS"]["Figures"]

#================
# CLOSING EVENT
#================

# Closing event: make sure to close both the infinite GUI loop AND the infinite pyplot loop!
def on_closing():
    turnOnVerboseWrapper_configurationFile()        # When using the command prompt, the functions have wrappers, this also saves the current configuration
    root.destroy()                                  # Destroy the tkinter window
    plt.close()                                     # Make sure the infinite plotting loop is closed too
root.protocol("WM_DELETE_WINDOW", on_closing)

#==================
# GLOBAL VARIABLES
#==================

# Research holds a lists of experiments, which each holds a list of simulations, 
# which each consist of multiple input files. 
root.input_files = []
root.Research = type('Dummy', (object,), {'content':{}})()
root.Research.data = {}
root.Research.experiments = []

#===================================
# FILL THE ROOT WINDOW AND STYLE IT 
#===================================

# Dont flash the screen while its loading so wait with showing the GUI
root.withdraw()

# Make sure the command prompt stays clean when the GUI is running, unless you want to debug
turnOffVerboseWrapper_configurationFile()

# Styling of the windows and widgets
ModifyStyling(root, theme=CONFIG["COLORS AND FONTS"]["Theme"])

# Create the header for the tabs
tab_header = ttk.Notebook(root, style='header.TNotebook')

# Add frames to the tab_header which are the tab windows
root.tab1 = ttk.Frame(tab_header) 
root.tab2 = ttk.Frame(tab_header) 
root.tab3 = ttk.Frame(tab_header) 
root.tab4 = ttk.Frame(tab_header)  
root.tab5 = ttk.Frame(tab_header)   
root.tab6 = ttk.Frame(tab_header)   
root.tab7 = ttk.Frame(tab_header)    

# Make the root accessible for the tabs
for tab in [root.tab1, root.tab2, root.tab3, root.tab4, root.tab5, root.tab6, root.tab7]:
    tab.root = root 

# Add the tabs to the tab header
tab_header.add(root.tab1, text='Simulations')
tab_header.add(root.tab2, text='Plasma profiles')
tab_header.add(root.tab3, text='Linear time traces')
tab_header.add(root.tab4, text='Linear spectra')
tab_header.add(root.tab5, text='Nonlinear time traces')
tab_header.add(root.tab6, text='Nonlinear surface plots')
tab_header.add(root.tab7, text='Videos')
tab_header.pack(expand=1, fill='both')

# Keep track of the plotting classes for the poppedout windows
root.graph_poppedOut = [] # Keep track of the plotting classes so we can access class.ax with the Optionswindow
root.canvasPoppedOut = [] # Keep track of the canvasses of the poppedout windows so we can draw_idle() on it

# Fill the tabs with widgets
root.tab_Simulations = TabSelectedFiles(root.tab1) # Class TabSelectedFiles makes the window for tab1 
root.tab_Profiles    = TabProfiles(root.tab2)      # Class TabProfiles makes the window for tab2
root.tab_Convergence = TabConvergence(root.tab3)   # Class TabConvergence makes the window for tab3 
root.tab_Linear      = TabLinear(root.tab4)        # Class tab_Linear makes the window for tab5  
root.tab_Nonlinear   = TabNonlinear(root.tab5)     # Class tab_Nonlinear makes the window for tab6  
root.tab_Videos      = TabVideos(root.tab7)        # Class TabVideos makes the window for tab7    
root.update_idletasks()

# Now that all elements are initialized, show the GUI
root.deiconify()

# Create a "dot dot dot" button which opens the preference window
PreferenceWindow(root)
root.update_idletasks()
tab_header.select(root.tab1)

#====================
# CHECK MEMORY USAGE
#====================

def check_memory():
    process = psutil.Process(os.getpid())
    print("MEMORY USAGE:", process.memory_info().rss*10**(-6), "MB")  # in bytes
    root.after(5000, check_memory)

#====================
# INFINITE GUI LOOP
#====================
 
# Tell Python to run the Tkinter event loop. This method listens for events, such as button clicks or keypresses 
# Blocks any code that comes after it from running until the window it’s called on is closed
# check_memory()
root.mainloop()







