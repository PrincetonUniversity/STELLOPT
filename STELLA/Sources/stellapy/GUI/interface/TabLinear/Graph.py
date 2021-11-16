
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

#################################################################
#                   CLASS FOR THE AX OBJECT
#################################################################

# Make the axis instance and the required attributes for the options window
# Either the plotting function is a class or we manually give it some class attributes like here
class Graph:
    
    def __init__(self, tab, figure, i):
        
        # Select the current figure
        plt.figure(figure.number)
        
        # Change the margins
        self.grid_specifications = gridspec.GridSpec(1, 1)
        if i==0: self.grid_specifications.update(top=0.95, left=0.15, right=0.95, bottom=0.15)
        if i==1: self.grid_specifications.update(top=0.95, left=0.15, right=0.95, bottom=0.15)

        # Create the axis object
        self.ax = plt.subplot(self.grid_specifications[0])
        
        # Set the variables for the OptionsWindow
        self.range = {"units" : "normalized", "x_scale" : "linear", "y_scale" : "linear"}    
        self.label = {"x" : None, "y" : None, "title" : None}
        self.layout = {"fontsize" : "N.A.", 'handlelength' : "N.A."}
        self.kx = 0.0
        self.ky = [0,100]
        
        # Set the names for the axis variables shown in the OptionsWindow
        if i==0: self.x_name = "Wavenumber";    self.y_name = "Growth rate"
        if i==1: self.x_name = "Parameter";     self.y_name = "Growth rate"
        
        # Set the ranges and labels to None when the graph is reset, this will triger the defaults
        self.load_defaults()
        
        # Save the graph and id
        self.id = i
        self.tab = tab
        return 
    
    #---------------
    def load_defaults(self):       
        self.range["x"] = None
        self.range["y"] = None
        for key in self.label.keys(): 
            self.label[key] = None
        return

    #---------------
    def update_quantitiesAndKeys(self):
        if self.id==0:
            y_quantity = self.tab.PlotLinearSpectrum.y_quantity
            if y_quantity=="omega":         self.y_name = "Frequency" 
            if y_quantity=="gamma":         self.y_name = "Growth rate" 
            if y_quantity=="gamma/ky**2":   self.y_name = "Growthrate/ky**2 versus modes"
        if self.id==1:
            self.x_name = self.tab.PlotParameterInfluence.key 
            y_quantity  = self.tab.PlotParameterInfluence.y_quantity
            if y_quantity=="omega":         self.y_name = "Frequency" 
            if y_quantity=="gamma":         self.y_name = "Growth rate" 
            if y_quantity=="gamma/ky**2":   self.y_name = "Growthrate/ky**2 versus modes"
        return
    
    #---------------
    def update_rangesAndLabels(self):

        # Get the ranges, labels and titles
        self.range["x"]     = self.ax.get_xlim()
        self.range["y"]     = self.ax.get_ylim()
        self.label["x"]     = self.ax.get_xlabel()
        self.label["y"]     = self.ax.get_ylabel()
        self.label["title"] = self.ax.get_title()






