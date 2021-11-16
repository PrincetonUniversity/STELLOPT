#================================================================
# When browse is hit: select files
#================================================================

# Load modules
import tkinter as tk
from tkinter import ttk

class Progress():

#==================================================================
# INITIALIZE THE FRAME WITH INFORMATION ON THE SELECTED inspectFiles
#==================================================================

    # Make a frame for the inspectFiles inside the window of the first tab "Select inspectFiles"
    def __init__(self, tab, anchor=None, length=None):

        #====================================================================
        # Safe info from the root so that it can be passed on or used locally
        #====================================================================
           
        self.root = tab.root 
        frame = tab.frame_progress

        #===========
        # VARIABLES
        #===========

        # Message to show progress  
        self.txt_information=tk.StringVar()               
    
        #============================
        # WIDGETS CREATION FOR FRAME
        #============================

        # Load the current style
        self.style = ttk.Style(self.root)
        
        # create progress bar
        self.bar = ttk.Progressbar(frame, orient = tk.HORIZONTAL, style="Labeledself", length=length) 
        self.bar.config(mode = 'determinate', maximum=100, value = 0)
        self.style.configure("Labeledself", text="No tasks are running.")

        #=====================
        # WIDGETS ARRANGEMENT
        #=====================
        if anchor == None:
            tk.Grid.rowconfigure(   frame, 0, weight=0) # Progress bar, make as small as possible
            tk.Grid.columnconfigure(frame, 0, weight=1) # Elongate along y
            self.bar.grid(in_=frame, row=0, column=0, padx=2, pady=0, sticky='nesw', ipady=5)
        if anchor == "center":
            tk.Grid.rowconfigure(   frame, 0, weight=1) # Spacing above
            tk.Grid.rowconfigure(   frame, 1, weight=0) # Progress bar, make as small as possible
            tk.Grid.rowconfigure(   frame, 2, weight=1) # Spacing below
            tk.Grid.columnconfigure(frame, 0, weight=1) # Spacing left
            tk.Grid.columnconfigure(frame, 1, weight=0) # Progress bar, make as small as possible
            tk.Grid.columnconfigure(frame, 2, weight=1) # Spacing right
            self.bar.grid(in_=frame, row=1, column=1, padx=2, pady=0, sticky='nesw', ipady=5)
        if anchor == "fill":
            tk.Grid.rowconfigure(   frame, 0, weight=0) # Progress bar, make as small as possible
            tk.Grid.columnconfigure(frame, 0, weight=1) # Elongate along y
            self.bar.grid(in_=frame, row=0, column=0, padx=0, pady=0, sticky='nesw', ipady=5)


#========================
# Start the progress bar
#========================
    def start(self, message):
        self.root.update_idletasks()
        self.bar['value'] = 0
        self.style.configure("Labeledself", text=message)
        self.root.update_idletasks()

#======================
# Step progress bar up
#======================
    def move(self, position, message):
        self.root.update_idletasks()
        self.bar['value'] = position
        self.style.configure("Labeledself", text=message)
        self.root.update_idletasks()

#========================
# Finish the progress bar
#========================
    def finish(self):
        self.root.update_idletasks()
        self.bar['value'] = 0
        self.style.configure("Labeledself", text="No tasks are running.")
        self.root.update_idletasks()


















