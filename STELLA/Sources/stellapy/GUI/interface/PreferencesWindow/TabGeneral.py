

#################################################################
#                       GENERAL TAB
#################################################################

# TODO: Change os.environ to CONGIG["PATH"]

# Load modules
import tkinter as tk
from tkinter import ttk

# Load personal modules
from stellapy.config import CONFIG #@unresolvedimport
from stellapy.GUI.graph_tools import PAD_TITLE, PAD_LABEL, PAD_ENTRY  #@unresolvedimport

#===================
# General tab class 
#===================
class TabGeneral:
    
    def __init__(self, tab):
        
        # Create the frame 
        self.tab_general = ttk.Frame(tab)
        self.tab_general.pack(expand=1, fill=tk.BOTH)
        
        # Configure the frame
        tk.Grid.rowconfigure(   self.tab_general, 0, weight=0)     # Text editor
        tk.Grid.rowconfigure(   self.tab_general, 1, weight=0)     # Paths
        tk.Grid.columnconfigure(self.tab_general, 0, weight=1) 
        
        # Change the text editor to display the text files
        self.init_chooseTextEditor(self.tab_general)
        
        # Set the paths to the code and simulation directories
        self.init_pathChoices(self.tab_general)

        
#=====================
# Choose text editor
#=====================
    def init_chooseTextEditor(self, tab):
        
        # Set self as the frame
        self.frame_textEditor = ttk.Frame(tab) 
        self.frame_textEditor.grid(row=0, column=0, sticky="WE")

        # Get the current text editor
        if   CONFIG['GENERAL']['TextEditor'] == "Emacs": int_editor = 1
        elif CONFIG['GENERAL']['TextEditor'] == "Gedit": int_editor = 2
        else                                           : int_editor = 1  # If nothing is correst, set to emacs
        
        # Set the default value of the editor 
        self.var_editor  = tk.IntVar(value=int_editor)
        
        # In the section "Text editor" give the options "Emacs" and "Gedit"
        self.lbl_editor = ttk.Label(self.frame_textEditor, text="Text editor", style='prefTitle.TLabel')
        self.rbn_emacs  = ttk.Radiobutton(self.frame_textEditor, text='  Emacs', variable=self.var_editor, value=1, style='prefRadio.TRadiobutton')
        self.rbn_gedit  = ttk.Radiobutton(self.frame_textEditor, text='  Gedit', variable=self.var_editor, value=2, style='prefRadio.TRadiobutton')
        
        # Configure the frame
        tk.Grid.rowconfigure(   self.frame_textEditor, 0, weight=1)     
        tk.Grid.rowconfigure(   self.frame_textEditor, 1, weight=1)      
        tk.Grid.rowconfigure(   self.frame_textEditor, 2, weight=1)   
        tk.Grid.columnconfigure(self.frame_textEditor, 0, weight=1) 
        
        # Add the options to the frame
        self.lbl_editor.grid( row=0, column=0, **PAD_TITLE)
        self.rbn_emacs.grid(  row=1, column=0, **PAD_LABEL)
        self.rbn_gedit.grid(  row=2, column=0, **PAD_LABEL)
        
    def apply_changesTextEdiror(self):
        if self.var_editor.get()==1: CONFIG['GENERAL']['TextEditor'] = 'emacs'
        if self.var_editor.get()==2: CONFIG['GENERAL']['TextEditor'] = 'gedit'
        return

#================
# Set thes paths
#================
    def init_pathChoices(self, tab):
        
        # Set self as the frame
        self.frame_pathChoices = ttk.Frame(tab) 
        self.frame_pathChoices.grid(row=1, column=0, sticky="WE")
        
        # Create a title for the path choices
        self.lbl_paths1 = ttk.Label(self.frame_pathChoices, text="Code directory", style='prefTitle.TLabel')
        self.lbl_paths2 = ttk.Label(self.frame_pathChoices, text="Simulations directory", style='prefTitle.TLabel')
        
        # Get the current paths
        self.var_pathStella = tk.StringVar(value=CONFIG["PATHS"]["Stella"])
        self.var_pathStellapy = tk.StringVar(value=CONFIG["PATHS"]["Stellapy"])
        self.var_pathBrowsFiles = tk.StringVar(value=CONFIG["PATHS"]["Default directory to browse files"])
        self.var_pathBrowsFolders = tk.StringVar(value=CONFIG["PATHS"]["Default directory to browse folders"])

        # Createlabels and entries for the four required paths
        self.lbl_pathStella = ttk.Label(self.frame_pathChoices, text='Path to the stella folder:  ')
        self.ent_pathStella = ttk.Entry(self.frame_pathChoices, textvariable=self.var_pathStella)
        self.lbl_pathStellapy = ttk.Label(self.frame_pathChoices, text='Path to the stellapy folder:  ')
        self.ent_pathStellapy = ttk.Entry(self.frame_pathChoices, textvariable=self.var_pathStellapy)
        self.lbl_pathBrowsFiles = ttk.Label(self.frame_pathChoices, text='Default directory to browse files:')
        self.ent_pathBrowsFiles = ttk.Entry(self.frame_pathChoices, textvariable=self.var_pathBrowsFiles)
        self.lbl_pathBrowsFolders = ttk.Label(self.frame_pathChoices, text='Default directory to browse folders:')
        self.ent_pathBrowsFolders = ttk.Entry(self.frame_pathChoices, textvariable=self.var_pathBrowsFolders)
        
        # Configure the frame
        tk.Grid.rowconfigure(   self.frame_pathChoices, 0, weight=1)     
        tk.Grid.rowconfigure(   self.frame_pathChoices, 1, weight=1)      
        tk.Grid.rowconfigure(   self.frame_pathChoices, 2, weight=1)       
        tk.Grid.rowconfigure(   self.frame_pathChoices, 3, weight=1)      
        tk.Grid.rowconfigure(   self.frame_pathChoices, 4, weight=1)  
        tk.Grid.columnconfigure(self.frame_pathChoices, 0, weight=0)   # Labels
        tk.Grid.columnconfigure(self.frame_pathChoices, 1, weight=1)   # Entries
        
        # Add the labels and entries to the frame
        self.lbl_paths1.grid(           row=0, column=0, columnspan=2, **PAD_TITLE)
        self.lbl_pathStella.grid(       row=1, column=0, **PAD_LABEL)
        self.lbl_pathStellapy.grid(     row=2, column=0, **PAD_LABEL)
        self.ent_pathStella.grid(       row=1, column=1, **PAD_ENTRY)
        self.ent_pathStellapy.grid(     row=2, column=1, **PAD_ENTRY)
        
        self.lbl_paths2.grid(           row=3, column=0, columnspan=2, **PAD_TITLE)
        self.lbl_pathBrowsFiles.grid(   row=4, column=0, **PAD_LABEL)
        self.lbl_pathBrowsFolders.grid( row=5, column=0, **PAD_LABEL)
        self.ent_pathBrowsFiles.grid(   row=4, column=1, **PAD_ENTRY)
        self.ent_pathBrowsFolders.grid( row=5, column=1, **PAD_ENTRY)

    def apply_changesPaths(self):
        # TODO: Apply changes paths
        return
        
