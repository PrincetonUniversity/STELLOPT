
#################################################################
#               STYLING OF THE GUI AND WIDGETS
#################################################################
'''
Non-active background: #33393b
Active background: #474d4e
'''

# Load modules
import os
import tkinter as tk
from tkinter import ttk
from pathlib import Path

# Personal modules
from stellapy.config import CONFIG #@unresolvedimport
from stellapy.GUI.utils import print_specs_buttons

# Global variables used throughout the GUI preference windows
PAD_TITLE = {"sticky": "W",  "pady" : (25,2), "padx" : (20,2)}
PAD_LABEL = {"sticky": "W",  "pady" : (2,2),  "padx" : (40,2)}
PAD_ENTRY = {"sticky": "WE", "pady" : (2,2),  "padx" : (40,20)}

# Global variables used throughout the GUI option frames
PAD_TITLE2 = {"sticky": "W",   "pady" : (2,2),  "padx" : (2,2)}
PAD_LABEL2 = {"sticky": "W",   "pady" : (1,1),  "padx" : (15,2)}
PAD_ENTRY2 = {"sticky": "WE",  "pady" : (1,1),  "padx" : (2,2)}
PAD_DIVID2 = {"sticky": "NSEW","pady" : (1,1),  "padx" : (0,0)}
PAD_ENTRYR = {"sticky": "E",  "pady" : (1,1),  "padx" : (2,2)}

# Modify the default styling
def ModifyStyling(root, theme="awdark"):
    
    # Load the themes that I downloaded from GitHub
    _load_predefinedThemes(root)
    
    # Set the current theme: either "awdark" or "awlight"
    set_predefinedTheme(root, theme)
    
    # The theme only applies to ttk widgets, so manually set the colors of the tk widgets
    load_awthemesForTkWidgets(root, theme)
    
    # Make some more personal adjustments to the predefined themes
    adjust_predefinedThemes(root, theme)
    
    # Derive some custom styles from the standard styles
    load_customStyles(root, theme)
    
def _load_predefinedThemes(root):
    ''' Load the themes that I downloaded from GitHub '''
    
    # Marconi has tk and tcl version 8.5 which don't support the themes.
    divider = '\\' if (os.name == 'nt') else '/'
    if os.getcwd().split(divider)[1] != 'marconi_work' and os.getcwd().split(divider)[1] != 'marconi':
    
        # Path of the predefined theme
        path = Path(CONFIG['PATHS']['stellapy']+'/GUI/themes/')
        path = path.as_posix()
        
        # Load a predefined theme 
        root.tk.eval("""
            set base_theme_dir """+path+"""
    
            package ifneeded awthemes 9.3.1 \
                [list source [file join $base_theme_dir awthemes.tcl]]
            package ifneeded colorutils 4.8 \
                [list source [file join $base_theme_dir colorutils.tcl]]
            package ifneeded awdark 7.7 \
                [list source [file join $base_theme_dir awdark.tcl]]
            package ifneeded awlight 7.6 \
                [list source [file join $base_theme_dir awdark.tcl]]
            # ... (you can add the other themes from the package if you want
            """)
    
        # Load the dark and light themes in this package
        root.tk.call("package", "require", 'awdark')
        root.tk.call("package", "require", 'awlight')
    

def set_predefinedTheme(root, theme):
    ''' Set the current theme: either "awdark" or "awlight" '''
    
    # Set the dark theme as the default
    style = ttk.Style(root)
    
    # Marconi has tk and tcl version 8.5 which don't support the themes.
    divider = '\\' if (os.name == 'nt') else '/'
    if os.getcwd().split(divider)[1] != 'marconi_work' and os.getcwd().split(divider)[1] != 'marconi':
        style.theme_use(theme)  
    
    # Remember the choice
    root.theme = theme
    
    # Load the corresponding colors to make more custom widgets
    if theme == "awdark":
        root.color = {
            'bg' : '#33393b',    # Background color (windows and labels)
            'bbg': '#232829',    # Button Background color (buttons)
            'fg' : '#ffffff',    # Foreground color (white text color) 
            'abg': '#474d4e',    # Active background color (when hovering over button)
            'bc' : '#000000',    # Border color for tabs (black)
            'bar': '#000000',    # Background color of the bar behind the notebook header (black)
            'canvas': '#33393b'} # Color of the plotting canvas 
        
    if theme == "awlight": 
        root.color = {
            'bg' : '#e8e8e7',    # Background color (windows and labels)
            'bbg': '#d1d1d0',    # Button Background color (buttons)
            'fg' : '#000000',    # Foreground color (black text color) 
            'abg': '#d1d1d0',    # Active background color (when hovering over button)
            'bc' : 'gray60',     # Border color for tabs (black)
            'bar': '#ffffff',    # Background color of the bar behind the notebook header (white)
            'canvas' : '#e8e8e7'} # Color of the plotting canvas (can't be too dark in dark theme)

def load_awthemesForTkWidgets(root, theme):
    ''' Create matching styles for the tk widgets since the theme is only for ttk widgets.  ''' 
             
    root.awthemes = { 
        'toplevel' : {
            'bg': root.color['bg']},\
        'label' : {
            'bg': root.color['bg'],\
            'fg': root.color['fg'],\
            'activebackground': root.color['bg'],\
            'bd': '1'},\
        'button' : {
            'bg': root.color['bbg'],\
            'fg': root.color['fg'],\
            'activebackground': root.color['abg']},\
        'labelframe' : {
            'borderwidth' : 10,\
            'padding' : (10, 10, 10, 10),\
            'relief' : tk.SUNKEN},\
        'labelframe2' : {
            'borderwidth' : 10,\
            'padding' : (10, 2, 10, 2),\
            'relief' : tk.SUNKEN},\
        'scrollableFrame' : {
            'bg' : root.color['bg'], \
            'bd' : 0, \
            'highlightcolor' : root.color['bg'], \
            'highlightthickness' : 0},\
        'crossButton' : {
            'width' : 20, \
            'height' : 20, \
            'bg' : root.color['bg'], \
            'borderwidth' : 0, \
            'highlightthickness' : 0, \
            'activebackground': root.color['abg']}}


def adjust_predefinedThemes(root, theme):
    ''' Make some more personal adjustments to the predefined themes. '''
    
    # Get the current style
    style = ttk.Style(root)

    # Manually configure the ttk widgets that don't match the current theme or that I want to be different
    style.configure('TNotebook.Tab',\
        font = ('MS Reference Sans Serif', '11', 'italic'),\
        highlightthickness = 0,\
        focuscolor=root.color['bg'])
    
    style.map('TNotebook.Tab',\
        background=[('active', root.color['bg']), ('disabled', root.color['bg'])]) # Background of tab
    
    style.configure('TLabelframe',\
        borderwidth = 5,\
        bordercolor=root.color['bc'],\
        lightcolor=root.color['bg'],\
        darkcolor=root.color['bg'],\
        highlightthickness = 5)
    
    style.configure('TLabelframe.Label',\
        font='arial 12 bold')
    
    style.configure('TFrame',\
        highlightcolor='green',\
        bd=0, borderwidth = 0,\
        highlightthickness = 0)
    
    style.configure('TLabel',\
        borderwidth = 0,\
        highlightthickness = 0)
    
    style.configure('TCheckbutton', 
        borderwidth = 0,\
        highlightthickness = 0,\
        focuscolor=root.color['bg']) # Focus is the dashed line around when selected
    
    style.map('TCheckbutton', 
        background=[('active', root.color['bg'])]) # Dynamic options
    
    style.configure('TProgressbar',\
        bg = root.color['bg'])
    
    style.configure("Treeview", \
        highlightthickness=0, \
        borderwidth = 0,\
        bd=0)
    
    style.layout("Treeview", [('Treeview.treearea', {'sticky': 'nswe'})]) 
    
    style.configure("Treeview.Heading", \
        font=('MS Reference Sans Serif', 10, 'bold')) 


def load_customStyles(root, theme):
        
    # Get the current style
    style = ttk.Style(root)
    
    # Style the main tab control of the GUI (see stella_GUI.py)
    style.configure('header.TNotebook', \
                    bordercolor=root.color['bc'],      # Border color around notebook
                    tabmargins=[15, 15, 2, 0],         # [left margin, upper margin, right margin, margin between tab and frames]
                    background=root.color['bar'])      # color behind the notebook
    style.configure('header.TNotebook.Tab', \
                    background =root.color['bbg'],
                    bordercolor=root.color['bc'],      # Border color around tab
                    font=('MS Reference Sans Serif', '11', 'italic'), # Font of the text in the tab
                    padding=[5, 1],                    # Padding around the text in the tab
                    foreground=root.color['fg'],       # Change text font of tab\\ 
                    focuscolor=root.color['bg'])       # Dotted line around tab when selected
    style.map('header.TNotebook.Tab', 
                    background=[('selected', root.color['bg'])], 
                    expand=[("selected", [2, 3, 3, 0])])

    # Style the sub tab control which is placed on each main tab
    style.configure('subheader.TNotebook', \
                    bordercolor=root.color['bc'],       # Border color around notebook
                    tabmargins=[15, 15, 2, 0],          # [left margin, upper margin, right margin, margin between tab and frames]
                    background= root.color['bg'])       # color behind the notebook
    style.configure('subheader.TNotebook.Tab', \
                    background = root.color['bbg'],
                    bordercolor= root.color['bc'],      # Border color around tab
                    font=('MS Reference Sans Serif', '11', 'italic'), # Font of the text in the tab
                    padding=[5, 1],                     # Padding around the text in the tab
                    foreground= root.color['fg'],       # Change text font of tab\\ 
                    focuscolor= root.color['bg'])       # Dotted line around tab when selected
    style.map('subheader.TNotebook.Tab', 
                    background=[('selected', root.color['bg'])], 
                    expand=[("selected", [2, 3, 3, 0])])
    
    # Buttons: remove dashed line when selected
    style.configure('TRadiobutton', focuscolor=root.color['bg'])
    style.map('TButton', focuscolor=[('pressed', root.color['bbg']), ('active', root.color['bbg']), ('selected', root.color['bbg']), ('disabled', root.color['bbg']), ('focus', root.color['bbg'])])   
    style.map('TRadiobutton', focuscolor=[('pressed', root.color['bg']), ('active', root.color['bg']), ('selected', root.color['bg']), ('disabled', root.color['bg'])])   
    
    # For option frames we have:
    #    the parameter which we put in text font
    #    the value which we put in code font
    #    the frame around value2 to change its borders
    style.configure('title_LBold.TLabel', font=("MS Reference Sans Serif", 12, "bold"), justify=tk.LEFT, anchor=tk.LEFT)
    style.configure('title_CBold.TLabel', font=("MS Reference Sans Serif", 12, "bold"), justify=tk.CENTER, anchor=tk.CENTER)
    style.configure('opt_paraC.TLabel', font=("MS Reference Sans Serif", 10), justify=tk.CENTER, anchor=tk.CENTER)
    style.configure('opt_paraL.TLabel', font=("MS Reference Sans Serif", 10), justify=tk.LEFT, anchor=tk.LEFT)
    style.configure('opt_valueC.TLabel', font=("Courier New", 11), justify=tk.CENTER, anchor=tk.CENTER)
    style.configure('opt_valueL.TLabel', font=("Courier New", 11), justify=tk.LEFT, anchor=tk.LEFT)
    style.configure('opt_sign.TLabel',  font=("Courier New", 12), justify=tk.CENTER, anchor=tk.CENTER) # For omega/gamma/...
    style.configure('opt_sign13.TLabel',  font=("Courier New", 13), justify=tk.CENTER, anchor=tk.CENTER) # For omega/gamma/...
    style.configure('opt_frame.TFrame', borderwidth=6, relief=tk.SUNKEN, background=root.color['bbg'])
    bbg = root.color['bbg'] # lightcolor is the bordercolor
    entry_background = {"padding":(10,4,10,4), "lightcolor":bbg, "background":bbg, "fieldbackground":bbg, "bordercolor":bbg, "highlightcolor":bbg, "focuscolor":bbg}
    style.configure('opt_valueR.TEntry', font=("Courier New", 12), justify=tk.RIGHT, anchor=tk.RIGHT, **entry_background)
    style.configure('opt_valueC.TEntry', font=("Courier New", 12), justify=tk.CENTER, anchor=tk.CENTER, **entry_background)
    style.configure('opt_valueCBold.TEntry',font=("Courier New", 12, "bold"), justify=tk.CENTER, anchor=tk.CENTER, borderwidth=3, relief=tk.FLAT, **entry_background)
    entry_background = {"lightcolor":bbg, "background":bbg, "fieldbackground":bbg, "bordercolor":bbg, "highlightcolor":bbg, "focuscolor":bbg}
    style.configure('tree.TEntry', **entry_background)
    bbg = root.color['bc'] # lightcolor is the bordercolor, this is the selected border color
    style.map('opt_valueR.TEntry',lightcolor=[('pressed',bbg), ('active',bbg), ('selected',bbg), ('disabled',bbg), ('focus',bbg)])   
    style.map('opt_valueC.TEntry', lightcolor=[('pressed',bbg), ('active',bbg), ('selected',bbg), ('disabled',bbg), ('focus',bbg)])   
    style.map('opt_valueCBold.TEntry', lightcolor=[('pressed',bbg), ('active',bbg), ('selected',bbg), ('disabled',bbg), ('focus',bbg)])   
    
    # Buttons for on the toolbar
    style.configure('toolbar.TButton',\
                    bordercolor=root.color['bar'],\
                    padding=[6,6]) # only the first one along x does something, y is expanded to fit
    
    # Titles and radiobuttons in the preference menu
    style.configure('prefTitle.TLabel', font=("MS Reference Sans Serif", 10, "bold"), justify=tk.LEFT, anchor=tk.LEFT)
    style.configure('prefRadio.TRadiobutton', font=("MS Reference Sans Serif", 10), justify=tk.LEFT, anchor=tk.LEFT)
    
    # Title of folder in the simulations scrollable frame
    style.configure('title.TLabel', font='arial 12 italic bold', foreground="gray3")
    
    # Tabbed display for the species in the simulations tab
    style.configure('species.TNotebook', \
                    bordercolor=root.color['bc'],   # Border color around notebook (the window frame)
                    darkcolor=root.color['bg'],\
                    lightcolor=root.color['bg'],\
                    highlightcolor=root.color['bg'],\
                    tabmargins=[0, 0, 0, 0],        # [left margin, upper margin, right margin, margin between tab and frames] # [15, 15, 2, 0]
                    background=root.color['bg'])    # color behind the notebook
    style.configure('species.TNotebook.Tab', \
                    background=root.color['bbg'],
                    bordercolor=root.color['bc'],   # Border color around tab
                    darkcolor=root.color['bg'],\
                    lightcolor=root.color['bg'],\
                    borderwidth=0,
                    relief=tk.FLAT,\
                    font=('MS Reference Sans Serif', '11', 'italic'), # Font of the text in the tab
                    padding=[5, 1],                 # Padding around the text in the tab
                    foreground=root.color['fg'],    # Change text font of tab\\ 
                    focuscolor=root.color['bg'])    # Dotted line around tab when selected
    style.map('species.TNotebook.Tab', 
              background=[('selected', root.color['bg'])],\
              lightcolor=[('selected', root.color['bg'])])
    
    # Dropdown menu in <all parameters> in simulations tab
    style.configure('gray.TMenubutton', \
                    foreground=root.color['fg'], \
                    font=('arial', 11, 'bold'))
    
    # Dropdown menu in Options winsow
    style.configure('option.TMenubutton', \
                    background=root.color['bbg'],\
                    foreground=root.color['fg'])
                    
    # Button next to the notebook header to open the preference window
    style.configure('dot.TFrame', \
                    bordercolor=root.color['bar'], \
                    highlightcolor=root.color['bar'], \
                    background=root.color['bar'], \
                    lightcolor=root.color['bar'], \
                    darkcolor=root.color['bar'], \
                    focuscolor=root.color['bar'])
    
    # Create a button that will open the preferences window
    style.configure('dot.TButton', \
                    font=("Courier New", 13, "bold"), \
                    bordercolor=root.color['bar'], \
                    highlightcolor=root.color['bar'], \
                    foreground=root.color['fg'], \
                    background=root.color['bar'], \
                    lightcolor=root.color['bar'], \
                    darkcolor=root.color['bar'], \
                    focuscolor=root.color['bar'])
        
    
    # Add a label with information to the progressbar style
    style.layout("Labeledself",
                 [('Labeledself.trough',
                   {'children': [('Labeledself.pbar',
                                  {'side': 'left', 'sticky': 'ns'}),
                                 ("Labeledself.label",
                                  {"sticky": "w"})],
                   'sticky': 'nswe'})])
        
def _unused_code(root):

    # Print the style of the widgets to check their colors
    if False:
        print_specs_buttons(ttk.Style(root))




