
import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from stellapy.config import CONFIG 

class PoppedOutWindow:

    def __init__(self, root):
        
        # Save the root
        self.root = root
        
        # Make the toplevel window 
        self.poppedout_window = tk.Toplevel() 
        self.poppedout_window.geometry('1200x700')
        self.poppedout_window.iconphoto(False, tk.PhotoImage(file=CONFIG['PATHS']['stellapy']+"GUI/images/qt4_editor_options.png"))
        
        # Create a frame for the canvas
        self.frame = ttk.Frame(self.poppedout_window)
        self.frame.grid(row=0, column=0, stick='NSEW')
        tk.Grid.rowconfigure(self.poppedout_window, 0, weight=1) 
        tk.Grid.columnconfigure(self.poppedout_window, 0, weight=1) 
        tk.Grid.rowconfigure(self.frame, 0, weight=1) 
        tk.Grid.columnconfigure(self.frame, 0, weight=1) 
        
        # Create the figure for the popped out window
        self.figure = plt.figure("popout "+str(len(root.graph_poppedOut)))
        self.figure.set_tight_layout(False)
        
        # Get the ID of the popped out window (at this point the graph is not in root.graph_poppedOut yet)
        self.poppedout_id = len(root.graph_poppedOut)
        
        # Closing event: remove the extra plotting class
        def on_closing():
            
            # Remove the class for the graph to save memory but keep the indexes 
            del self.root.graph_poppedOut[self.poppedout_id]
            del self.root.canvasPoppedOut[self.poppedout_id]
            self.root.graph_poppedOut.insert(self.poppedout_id, "dummy index") 
            self.root.canvasPoppedOut.insert(self.poppedout_id, "dummy index") 
            
            # Destroy the window
            self.poppedout_window.destroy()
            
            # Destroy the figure
            plt.close("popout "+str(self.poppedout_id))
            
            # Clear the list of dummy indexes when there are no popped out window
            if all(class_ == "dummy index" for class_ in self.root.graph_poppedOut):
                self.root.graph_poppedOut = []
                self.root.canvasPoppedOut = []
        
        # Add the closing event to the window
        self.poppedout_window.protocol("WM_DELETE_WINDOW", on_closing)
        return 


    def set_title(self, title):
        self.poppedout_window.title(title)
        