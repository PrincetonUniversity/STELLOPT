

#################################################################
#                       ABOUT TAB
#################################################################

# Load modules
import tkinter as tk
from tkinter import ttk

# Personal modules
from stellapy.GUI.graph_tools import PAD_TITLE, PAD_LABEL  #@unresolvedimport

#==================
# About tab class 
#==================
class TabAbout:
    
    def __init__(self, tab):

        # Create the frame 
        self.tab_about = ttk.Frame(tab)
        self.tab_about.pack(expand=1, fill=tk.BOTH)
        
        # Configure the frame
        tk.Grid.rowconfigure(   self.tab_about, 0, weight=0)   
        tk.Grid.rowconfigure(   self.tab_about, 1, weight=0)
        tk.Grid.rowconfigure(   self.tab_about, 2, weight=0)
        tk.Grid.rowconfigure(   self.tab_about, 3, weight=0)  
        tk.Grid.columnconfigure(self.tab_about, 0, weight=1) 
        
        # Fill it with a text explaining the GUI
        self.lbl_version = ttk.Label(self.tab_about, text="Stellapy Version 1.0", style='prefTitle.TLabel')
        self.lbl_author = ttk.Label(self.tab_about,text=\
            "\n  This GUI is being created by Hanne Thienpondt, \
             \n  with the help of Jose-Manuel García-Regaña,\
             \n  Bob Davies and Antonio González Jerez.\
             \n \
             \n  Email: Hanne.Thienpondt@outlook.com")
        self.lbl_themes1 = ttk.Label(self.tab_about, text="Themes", style='prefTitle.TLabel')        
        self.lbl_themes2 = ttk.Label(self.tab_about,text=\
            "\n  Light and dark themes available at https://sourceforge.net/projects/tcl-awthemes/.\
             \n \n")
        
        
        # Add the labels to the frame
        self.lbl_version.grid( row=0, column=0, **PAD_TITLE)
        self.lbl_author.grid(  row=1, column=0, **PAD_LABEL)
        self.lbl_themes1.grid( row=2, column=0, **PAD_TITLE)
        self.lbl_themes2.grid( row=3, column=0, **PAD_LABEL)




