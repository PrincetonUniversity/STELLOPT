
#################################################################
#                       CUSTOM TOOLBAR
#################################################################

# Load modules
import tkinter as tk
from tkinter import ttk
import matplotlib
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk  
matplotlib.rcParams["toolbar"] = "toolbar2"

# Personal modules
from .OptionsWindow import OptionsWindow  
from stellapy.config import CONFIG  
from stellapy.GUI.graph_tools.ToolTip import ToolTip  

# Custom toolbar 
class CustomToolbar(NavigationToolbar2Tk): 
    """
    Customized Navigator object
    - Removed mouse_move event.x and event.y display
    - Changed buttons layout
    
    Parameters
    ----------
    root : Tk()
        Root of the tkinter application.
    canvas_ : FigureCanvasTkAgg
        The canvas to which the toolbar will be linked.
    parent_ : ttk.Frame
        The frame to which the canvas_ is attached.
    reset_graph : Function
        Function to replot the graph on the canvas from scratch.
    poppedout_id : int
        Number to know to which topwindow this toolbar is linked.
    """ 

    def __init__(self, root, canvas_, parent_, master_class, axis_id=None, poppedout_id=None, direction="x"):

        # Safe the root and id
        self.root = root
        self.master_class = master_class
        self.axis_id = axis_id
        self.poppedout_id = poppedout_id
        self.direction = direction
        
        # Make sure we know which buttons to plot
        if not hasattr(self.master_class, 'btn_popout'):
            self.master_class.btn_popout = True
        if not hasattr(self.master_class, 'btn_reread'):
            self.master_class.btn_reread = True
        
        # Get the width of the canvas
        canvas_frame = canvas_.get_tk_widget(); 
        canvas_frame.update(); root.update_idletasks()
        self.width = canvas_frame.winfo_width()

        # Make the option window object for each toolbar, the id is used to know which canvas it is linked to
        self.class_optionWindow = OptionsWindow(self.root, self.master_class, self.axis_id, self.poppedout_id)
        
        # Safe the custom buttons
        self.customButtons = []
        
        # Initiate the standard toolbar (which is empty) and manually add it to the canvas so we can use grid()   
        self.toolitems = ()
        NavigationToolbar2Tk.__init__(self,canvas_,parent_,pack_toolbar=False)
        if self.direction=="x": self.grid(row=1, column=0, stick='NSEW')
        if self.direction=="y": self.grid(row=0, column=0, stick='NSEW')

        # There is a tk.Label object at the end of the toolbar which has the wrong color so remove it
        # See /usr/local/lib/python3.6/dist-packages/matplotlib/backends/_backend_tk.py line 521
        list_widgets = self.pack_slaves()
        for i in range(len(list_widgets)):
            list_widgets[i].destroy()
        
        # Define our own custom buttons
        self.tool_buttons = []; num = 0
        self.extratoolitems = (
            ('Save', 'Save the figure', 'icon', 'save_figure'),
            ('Home', 'Reset original view', 'icon', 'home'),
            ('Back', 'Back to previous view', 'icon', 'back'),
            ('Forward', 'Forward to next view', 'icon', 'forward'),
            ('Pan', 'Pan axes with left mouse, zoom with right', 'icon', 'pan'),
            ('Zoom', 'Zoom to rectangle', 'icon', 'zoom'),
            (None, None, None, None),
            ('Reset', 'Reset the figure', 'icon', 'reset'),
            ('Options', 'Figure options', 'icon', 'options'),
            ('PopOut', 'Open the figure in a new window', 'icon', 'popout_window'))
        
        # If the time is dark, use white icons instead of black icons
        extension = "_inv" if self.root.theme=="awdark" else ""
        self.toolbar_icons = ["filesave","home","back","forward","move","zoom_to_rect","reload","subplots","qt4_editor_options"]
        self.toolbar_icons = [ CONFIG['PATHS']['stellapy']+"GUI/images/"+icon+extension+".png" for icon in self.toolbar_icons ]
        self.toolbar_icons.insert(6, None)
        self.toolbar_icons.append(None)
        
        # Add the icons to the toolbar
        for text, tooltip_text, image_file, callback in self.extratoolitems: #@unusedvariable
            if text is None:
                self.add_customLabelSpacing(num)
            else:
                if  not (text=='PopOut' and self.master_class.btn_popout == "N.A.")\
                and not (text=='Reread' and self.master_class.btn_reread == "N.A."):
                    try:
                        button = self.add_customButtonWithIcon(text=text, file=self.toolbar_icons[num], command=getattr(self, callback), i=num)
                        if tooltip_text is not None:
                            # Custom tooltip class with a half second delay
                            ToolTip(self.root, button, tooltip_text)
                    except IndexError:
                        pass
            num+=1
    
        # Configure the grid: put fillers before and after the buttons so its centered
        if self.direction=="x": 
            tk.Grid.rowconfigure(self, 0, weight=0) 
            tk.Grid.columnconfigure(self, 0, weight=1)
            tk.Grid.columnconfigure(self, num+1, weight=1)
        if self.direction=="y": 
            tk.Grid.columnconfigure(self, 0, weight=0) 
            tk.Grid.rowconfigure(self, 0, weight=1)
            tk.Grid.rowconfigure(self, num+1, weight=1)
            tk.Grid.rowconfigure(self, num+2, weight=1)
    
    def add_customButtonWithIcon(self, text, file, command, i):
        im = tk.PhotoImage(file=file)
        self.customButtons.append(ttk.Button(master=self, text=text, image=im, command=command, style='toolbar.TButton'))
        self.customButtons[-1]._ntimage = im # Attach image to button, otherwise it disappears
        if self.direction=="x": self.customButtons[-1].grid(row=0, column=i+1)
        if self.direction=="y": self.customButtons[-1].grid(row=i+1, column=0)
        return self.customButtons[-1]
    
    def add_customLabelSpacing(self, i):
        label = ttk.Label(master=self, text="   ")
        label.grid(row=0, column=i+1)
        self.tool_buttons.append(label)
        return label

    # When the reset button is hit, replot the figure with its default ranges and labels
    def reset(self):
        self.master_class.reset_graph(axis_id=self.axis_id, poppedout_id=self.poppedout_id)
        return

    # When the options button is hit, a window opens to change the graph appearance
    def options(self):
        self.class_optionWindow.open_optionsWindow()
        return 
        
    # When the popout button is hit, the figure gets replotted in a new top window   
    def popout_window(self):
        self.master_class.popout_window(self.axis_id) 
        return 
        
    # To manually add more buttons
    def add_customButton(self, text, command, **kwargs):
        button = ttk.Button(master=self, text=text, command=command, style='toolbar.TButton', **kwargs)
        button.pack(side=tk.LEFT,fill="y", padx=5)
        return button

