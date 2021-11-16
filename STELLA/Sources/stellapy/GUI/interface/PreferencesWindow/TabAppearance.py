

#################################################################
#                       APPEARANCE TAB
#################################################################
# Load modules
import tkinter as tk
from tkinter import ttk

# Load personal modules
from stellapy.config import CONFIG #@unresolvedimport
from stellapy.GUI.utils import restart_program
from stellapy.GUI.graph_tools import ModifyStyling
from stellapy.GUI.graph_tools import PAD_TITLE, PAD_LABEL #@unresolvedimport

#======================
# Appearance tab class 
#======================
class TabAppearance:
    
    def __init__(self, tab):

        # Attach the root so we can carry it into the functions
        self.root = tab.root

        # Create the frame 
        self.tab_appearance = ttk.Frame(tab)
        self.tab_appearance.pack(expand=1, fill=tk.BOTH)
        
        # Configure the frame
        tk.Grid.rowconfigure(   self.tab_appearance, 0, weight=0) 
        tk.Grid.rowconfigure(   self.tab_appearance, 1, weight=0) 
        tk.Grid.columnconfigure(self.tab_appearance, 0, weight=1) 
        
        # Change the theme of the GUI: "awdark" or "awlight"
        self.init_themeOptions()
        
#========
# Theme
#========
    def init_themeOptions(self):
    
        # Read the configuration file to know the current theme
        if   CONFIG["COLORS AND FONTS"]["Theme"] == "awlight": int_theme = 1
        elif CONFIG["COLORS AND FONTS"]["Theme"] == "awdark":  int_theme = 2
        else:                                                  int_theme = 2 # If something went wrong choose dark theme
        
        # Set the default value of the theme 
        self.var_theme  = tk.IntVar(value=int_theme)
        
        # In the section "Theme" give the options "Light" and "Dark"
        self.lbl_theme  = ttk.Label(self.tab_appearance, text="Theme", style='prefTitle.TLabel')
        self.rbn_awlight = ttk.Radiobutton(self.tab_appearance, text='  Light', variable=self.var_theme, value=1)
        self.rbn_awdark  = ttk.Radiobutton(self.tab_appearance, text='  Dark', variable=self.var_theme, value=2)
        
        # Add the options to the frame
        self.lbl_theme.grid(    row=0, column=0, **PAD_TITLE)
        self.rbn_awlight.grid(  row=1, column=0, **PAD_LABEL)
        self.rbn_awdark.grid(   row=2, column=0, **PAD_LABEL)
        
        
    def apply_changesTheme(self):
        
        # Get the chosen theme
        chosen_theme = self.var_theme.get()
        
        # Get the current theme from the configuration file
        if   CONFIG["COLORS AND FONTS"]["Theme"] == "awlight": current_theme = 1
        elif CONFIG["COLORS AND FONTS"]["Theme"] == "awdark":  current_theme = 2
        else:                                                  current_theme = 2
        
        # Change the theme if it is different now and save it to the configuration file
        if self.var_theme.get() != current_theme:
            
            if chosen_theme==1: 
                CONFIG["COLORS AND FONTS"]["Theme"] = "awlight"
                restart_program()
            if chosen_theme==2: 
                CONFIG["COLORS AND FONTS"]["Theme"] = "awdark"
                restart_program()
            

        
        
        
        
        
        
        


