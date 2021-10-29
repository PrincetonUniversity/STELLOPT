
#=======================================
# Iniate the frame holding the graph
#=======================================

import tkinter as tk
from tkinter import ttk
import matplotlib; matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvas
from stellapy.GUI.graph_tools.CustomToolbar import CustomToolbar 


class CanvasForGraphs:
    
    def __init__(self, root, frame, master_class, figure, axis_id=None, poppedout_id=None, direction="x"):
        
        # Save the variables
        self.root = root
        self.frame = frame
        self.master_class = master_class
        self.figure = figure
        self.axis_id = axis_id
        self.poppedout_id = poppedout_id
        self.direction = direction
        
        # Create the frame, canvas and toolbar
        self.initiate_frameGraph()
        
        # Make the canvas available to the master class so we can draw_idle() it
        if poppedout_id==None: self.master_class.Canvas[axis_id] = self.Canvas
        if poppedout_id!=None: self.root.canvasPoppedOut.append(self.Canvas)

    #----------- Frame ----------
    def initiate_frameGraph(self):

        # Configure frame that holds the canvas
        tk.Grid.rowconfigure(   self.frame, 0, weight=1) 
        tk.Grid.columnconfigure(self.frame, 0, weight=1)  
    
        # Create frame for graph elements
        self.frame_canvas = ttk.Frame(self.frame)
        self.frame_canvas.grid(row=0, column=0, padx=(0,0), pady=(0,0), stick='NSEW') 
        
        # Configure the frame that holds the figure and toolbar
        if self.direction=="x":
            tk.Grid.rowconfigure(   self.frame_canvas, 0, weight=1) 
            tk.Grid.rowconfigure(   self.frame_canvas, 1, weight=0) 
            tk.Grid.columnconfigure(self.frame_canvas, 0, weight=1)  
        if self.direction=="y":
            tk.Grid.columnconfigure(self.frame_canvas, 0, weight=0) 
            tk.Grid.columnconfigure(self.frame_canvas, 1, weight=1) 
            tk.Grid.rowconfigure(   self.frame_canvas, 0, weight=1)  
        
        # Add the canvas and the toolbar
        self.initiate_canvas()
        self.initiate_toolbar()
        
    #----------- Canvas -----------
    def initiate_canvas(self):

        # Add the figure to a Canvas
        self.Canvas = FigureCanvas(self.figure, self.frame_canvas) 
        self.Canvas.get_tk_widget().configure(background=self.root.color['bg'])
        self.figure.patch.set_facecolor(self.root.color['bg']) 
        self.Canvas.draw()
        self.widget = self.Canvas.get_tk_widget() 
        if self.direction=="x": self.widget.grid(row=0, column=0, stick='NSEW')
        if self.direction=="y": self.widget.grid(row=0, column=1, stick='NSEW')

    #----------- Toolbar -----------
    def initiate_toolbar(self):
        
        # Add a toolbar to the canvas
        self.toolbar = CustomToolbar(self.root, self.Canvas, self.frame_canvas, self.master_class, self.axis_id, self.poppedout_id, self.direction)
        self.toolbar.config(background=self.root.color['bg'])
        self.toolbar.update()
        
