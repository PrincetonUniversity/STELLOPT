
#################################################################
#                   CLASS FOR THE FIRST TAB
#################################################################
''' 
https://blog.tecladocode.com/tkinter-scrollable-frames/#:~:text=In%20Tkinter%2C%20only%20the%20Canvas,That's%20what%20scrolling%20really%20is!
'''

# Load modules
import tkinter as tk
from tkinter import ttk

# Class ScrollableFrame
class ScrollableFrame(ttk.Frame):
    def __init__(self, container, *args, **kwargs):

        super().__init__(container, *args, **kwargs)
        
        # Create the canvas
        self.canvas = tk.Canvas(self, **container.root.awthemes['scrollableFrame'])
        self.scrollbar = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.scrollable_frame = ttk.Frame(self.canvas) #mailbox_frame

        # The <Configure> event triggers whenever the scrollable_frame changes size—i.e. usually when we add or remove widgets from within it
        self.scrollable_frame.bind("<Configure>", self.OnFrameConfigure)
        self.canvas.bind('<Configure>', self.FrameWidth)

        # Tell the canvas to actually draw the scrollable_frame inside itself
        self.canvas_frame = self.canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")

        # Configure the canvas so that when its y-position changes, the scrollbar moves:
        self.canvas.configure(yscrollcommand=self.scrollbar.set)

        # Use Pack or Grid to add the container, canvas, and scrollbar to the application window
        self.canvas.pack(side="left", fill="both", expand=True)
        self.scrollbar.pack(side="right", fill="y")

    def FrameWidth(self, event):
        canvas_width = event.width
        self.canvas.itemconfig(self.canvas_frame, width = canvas_width)

    def OnFrameConfigure(self, event):
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))
        
    def change_colors(self, root):
        self.canvas.configure(**root.awthemes['scrollableFrame'])



