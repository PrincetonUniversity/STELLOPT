#================================================================
# Return line and marker colors and styles
#================================================================

# Load the modules
import numpy as np
from matplotlib.pyplot import cm
from stellapy.utils import initiate_nesteddict
from itertools import takewhile

# stellapy.plot.load_lines_colors_markers
def load_labelsLinesMarkers(experiments, varied_values, colormap=None):
    ''' Load labels, lines, color and markers for the (simulation, groups) combinations.     

    Returns
    -------
    line_style, line_color, marker_color, marker_style : dictionaries with keys [simlutation] and [groups]
    '''

    # Initiate the dictionaries
    line_label   = initiate_nesteddict()
    marker_label = initiate_nesteddict()
    line_style   = initiate_nesteddict()
    line_color   = initiate_nesteddict()
    marker_color = initiate_nesteddict()
    marker_style = initiate_nesteddict() 
    
    # Initiate colors, lines, markers to chose from
    if colormap=="rainbow": l_colors = cm.rainbow( np.linspace(0,1,len( experiments ) + 1 ) )
    if colormap=="jet":     l_colors = cm.jet( np.linspace(0,1,len( experiments ) + 1 ) )
    if colormap==None:      l_colors = ['navy', 'crimson', 'black', 'darkgreen', 'purple', 'magenta']*5 
    m_colors = cm.jet( np.linspace(0,1,len( varied_values )) )
    l_styles = ['-', '--', ':','-.','densely dotted','densely dasged', 'densely']*5
    m_styles = ['o','X','D','P','s','x','d','p']*10
    
    # Find the common prefix of the experiments
    common_prefix = ''.join(c[0] for c in takewhile(lambda x:  all(x[0] == y for y in x), zip(*experiments))) 

    # If no specific case was chosen, given random colors
    for experiment in experiments:
        if common_prefix==experiments[0] or common_prefix=="":
            line_label[experiment] = str(experiment).replace("nl_","").replace("th_","").replace("_", " ")#.split("_")[0]
        if "=" in common_prefix:
            line_label[experiment] = str(experiment)
        else:
            try: line_label[experiment] = str(experiment).replace("_", " ")
            # line_label[experiment] = str(experiment).split(common_prefix)[-1].replace("_", " ")#.split("_")[0]
            except: line_label[experiment] = str(experiment).replace("_", " ")
        marker_label[experiment] = str(experiment).replace("__", ": ").replace("_", " ")
        line_color[experiment]   = l_colors[experiments.index(experiment)]  
        line_style[experiment]   = l_styles[experiments.index(experiment)]
        marker_color[experiment] = l_colors[experiments.index(experiment)]    
        marker_style[experiment] = m_styles[experiments.index(experiment)]
    
    for variation in varied_values: 
        line_label[variation]    = str(variation)
        marker_label[variation]  = str(variation)
        line_color[variation]    = m_colors[varied_values.index(variation)]  
        line_style[variation]    = '-'
        marker_color[variation]  = m_colors[varied_values.index(variation)]   
        marker_style[variation]  = m_styles[varied_values.index(variation)] 

    if experiments==varied_values and len(experiments)==1:
        return line_label[varied_values[0]], marker_label[varied_values[0]], line_style[varied_values[0]], line_color[varied_values[0]], marker_style[varied_values[0]], marker_color[varied_values[0]]
    return line_label, marker_label, line_style, line_color, marker_style, marker_color






    
