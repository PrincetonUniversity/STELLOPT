

import matplotlib.pyplot as plt
from stellapy.GUI.graph_tools import CanvasForGraphs, PoppedOutWindow
from stellapy.GUI.interface.TabLinear.Graph import Graph

#################################################################
#                    CLASS FOR THE CANVAS
#################################################################
   
class Canvasses: 
    
    def __init__(self, parent):
        ''' Initiate the canvas for the matplotlib plots.
        
        The toolbar requires the class to have a reset_graph and popout_window function
        The optionswindow requires the Graph objects to have attributes: ax, x_name, x_key, y_name, y_key, range, label
        '''

        # Make the parent available
        self.tab = parent
        self.root = parent.root
        
        # Create the figures for each canvas
        self.figure1 = plt.figure("linear1") 
        self.figure2 = plt.figure("linear2")   
        self.figure1.set_tight_layout(False) 
        self.figure2.set_tight_layout(False)
        self.figure1.patch.set_facecolor(self.root.color['canvas']) 
        self.figure2.patch.set_facecolor(self.root.color['canvas']) 
        
        # Put the canvas on the screen
        CanvasForGraphs(self.root, self.tab.frame_graph1, self.root.tab_Linear, self.figure1, axis_id=0)
        CanvasForGraphs(self.root, self.tab.frame_graph2, self.root.tab_Linear, self.figure2, axis_id=1)

        # Initiate the plotting class for the main window
        self.initiate_GraphClass()

    #--------------------------------------------
    def initiate_GraphClass(self, figure=None):
        
        # Set the color of the axis
        plt.rcParams['text.color']       = self.root.color['fg']
        plt.rcParams['axes.edgecolor']   = self.root.color['fg']
        plt.rcParams['axes.labelcolor']  = self.root.color['fg']
        plt.rcParams['xtick.color']      = self.root.color['fg']
        plt.rcParams['ytick.color']      = self.root.color['fg']
        plt.rcParams['axes.facecolor']   = self.root.color['canvas']
        plt.rcParams['figure.facecolor'] = self.root.color['canvas']
        
        # Initiate the plotting class for the main window
        if figure is None:
            # Make the axis instance and the required attributes for the options window through the class <graph>
            self.tab.Graph[0] = Graph(self.tab, self.figure1, 0)
            self.tab.Graph[1] = Graph(self.tab, self.figure2, 1)
            # Put the empty figure on the GUI
            self.tab.Canvas[0].draw_idle(); self.root.update_idletasks()
            self.tab.Canvas[1].draw_idle(); self.root.update_idletasks()
        
        # Or initiate the plotting class for a popped out window
        if figure is not None:
            # Make the axis instance and the required attributes for the options window through the class <graph>
            self.root.graph_poppedOut.append(Graph(self.tab, figure, 0))
            
        if True: return 
        
#################################################################
#                          METHODS
#################################################################

    def popout_window(self, axis_id):
        '''Replot the figure in a seperate window when the button on the toolbar is clicked.'''
        
        # Create a popped out window
        poppedout_window = PoppedOutWindow(self.root)
        poppedout_id = poppedout_window.poppedout_id
        
        # Add a fitting title to the popped out window
        if axis_id==0:
            y_quantity = self.tab.PlotLinearSpectrum.y_quantity
            if y_quantity=='omega':         poppedout_window.set_title("Stellapy: Frequency spectrum")
            if y_quantity=='gamma':         poppedout_window.set_title("Stellapy: Growth rate spectrum")
            if y_quantity=='gamma/ky**2':   poppedout_window.set_title("Stellapy: Growth rate on ky**2 spectrum")
            # ...

        # Initiate the plotting class
        self.initiate_GraphClass(figure=poppedout_window.figure)
        
        # Initiate the canvas        
        CanvasForGraphs(self.root, poppedout_window.frame, self.root.tab_Linear, poppedout_window.figure, poppedout_id=poppedout_id)
        
        # Parse the data from the main canvas to the popped out canvas
        self.root.graph_poppedOut[poppedout_id].range = self.tab.Graph[axis_id].range
        self.root.graph_poppedOut[poppedout_id].label = self.tab.Graph[axis_id].label
        
        # Now plot the figure on the popped out canvas
        if axis_id==0: Plot = self.tab.PlotLinearSpectrum 
        if axis_id==1: Plot = self.tab.PlotParameterInfluence 
        Plot.replot = True; Plot.plot_graph(poppedout_id)
        
        # Update screen
        self.root.canvasPoppedOut[poppedout_id].draw_idle()
        return 
    
    #------------------------------------------------
    def reset_graph(self, axis_id, poppedout_id=None):
        ''' Method called by the "Reset graph" button on the canvas. '''
        
        # Perhaps the button was clicked on the popped out window
        if poppedout_id == None: Graph = self.tab.Graph[axis_id]
        if poppedout_id != None: Graph = self.root.graph_poppedOut[poppedout_id]
        
        # Reset the axis
        Graph.load_defaults()
        
        # Now plot the graph 
        if axis_id==0: Plot = self.tab.PlotLinearSpectrum 
        if axis_id==1: Plot = self.tab.PlotParameterInfluence 
        Plot.replot = True; Plot.plot_graph(None)
        return 

