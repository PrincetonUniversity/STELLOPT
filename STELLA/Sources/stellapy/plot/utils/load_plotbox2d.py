
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator
from .load_styleFigures import load_styleFigures

def load_plotbox2d(x_range=None, y_range=None, x_label=None, y_label=None,\
         fig_size=(18, 9), title=None, ax=None, log=False, title_size=None):
    ''' Function to create and manipulate the appearance of a plot 

    Returns
    -------
    ax : axis object

    Parameters
    ----------
    fig_size : tuple of int, optional  
    xrange, yrange : tuple of int, optional
    xlabel, ylabel, title : str, optional
    log : {False, True}
    ax : axis object
	    If None is given, a new figure is created
    '''
    
    # Set the axes colors and width
    load_styleFigures()

    # Create a new figure if no existing axis was give 
    if not ax:
        fig = plt.figure(figsize=fig_size)
        grid_specifications = gridspec.GridSpec(1, 1, figure=fig)
        grid_specifications.update(top=0.95, left=0.09, right=0.82, bottom=0.12)
        ax = plt.subplot(grid_specifications[0]) 

    # Change the appearance of the figure
    ax.xaxis.grid(color='grey', linestyle='-', linewidth=0.3)
    ax.yaxis.grid(color='grey', linestyle='-', linewidth=0.3)  
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_xlim(x_range)
    ax.set_ylim(y_range)
    ax.yaxis.set_minor_locator(AutoMinorLocator(10))
    ax.xaxis.set_minor_locator(AutoMinorLocator(10))
    if log == True:  ax.set_yscale('log');
    if log == False: ax.ticklabel_format(useOffset=False);
    if title != None: 
        if title_size == None: ax.set_title(title); print(title)
        if title_size == 'small': ax.set_title(title, fontsize=25);
    return ax
