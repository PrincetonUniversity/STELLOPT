''' CONFIGURATION FILE

This module is also a script which will automatically read the configuration file 
when the module is loaded. This assures that the constant CONFIG can be imported 
in other modules and that the corresponding configuration file will already exist 
and is read in.


Constants
---------
CONFIG: configparser.ConfigParser()   
    Configuration dictionary that will hold the configuration of the GUI.
'''

#################################################################
#               EDIT CONFIGURATION FILE
#################################################################

# Load modules
import tkinter as tk
import os, configparser
from pathlib import Path

# Create a configuration object
CONFIG = configparser.ConfigParser()                                                   

# Select the location of the configuration file: depends on the computer
divider = '\\' if (os.name == 'nt') else '/'
path_stellaGUI = os.path.dirname(os.path.abspath(__file__)).split("stellapy")[0]
if os.getcwd().split(divider)[1] == 'marconi_work' \
or os.getcwd().split(divider)[1] == 'marconi' \
or os.getcwd().split(divider)[1] == 'marconi_scratch':
    NAME_CONFIGURATIONFILE = Path(path_stellaGUI+"/stellapy/config/config_marconi.ini")
if os.getcwd().split(divider)[1] == 'home':
    NAME_CONFIGURATIONFILE = Path(path_stellaGUI+"/stellapy/config/config.ini")

#===============================================
# Edit configuration file
#===============================================

def read_configurationFile(root=None):
    ''' Check whether <config.ini> exists, if not write it, CONFIG will be read after calling this function.
    Reading the file will check it existence and whether the paths are correct. '''
    
    # The configuration file has not been read
    if "GENERAL" not in CONFIG:
        CONFIG.read(NAME_CONFIGURATIONFILE)
    
    # The configuration file didn't exist when it was read
    if "GENERAL" not in CONFIG:
        create_defaultConfigurationFile()
        return
    
    # The configuration file is already read
    if "GENERAL" in CONFIG:
        
        # Get the path to the current file
        path_currentFile = os.path.dirname(os.path.abspath(__file__))
        
        # Check whether the paths are set correctly
        if CONFIG['PATHS']['stella'] != path_currentFile.split("stellapy")[0] or \
           CONFIG['PATHS']['Stellapy'] != path_currentFile.split("config")[0]:
            message = "The path to the stella folder was not set correctly. \n \
                       The configuration file will be replaced with a default configuration file."
            if root: display_information(root, title="Warning", message=message) # For the GUI
            create_defaultConfigurationFile()
        return 
        

def create_defaultConfigurationFile():
    ''' Write a default configuration file and save it as "stella/stellapy/config/config.ini". '''
    
    # Get the path to the current file
    path_currentFile = os.path.dirname(os.path.abspath(__file__))
    
    # Default section
    CONFIG['DEFAULT'] = {
        'code' : 'Stellapy',\
        'use_verbosewrapper' : 'True'} # Only for the GUI it will be turned off programatically

    CONFIG['GENERAL'] = {
        'TextEditor' : 'emacs'}
    
    CONFIG['COLORS AND FONTS'] = {
        'Theme' : 'awlight'} 
    
    CONFIG['PATHS'] = {
        'Stella' : path_currentFile.split("stellapy")[0],\
        'Stellapy' : path_currentFile.split("config")[0],\
        'Default directory to browse files': path_currentFile.split("config")[0]+'examples/LinearKyScan_OneKyPerFile_AdiabaticElectrons_W7xr348',\
        'Default directory to browse folders': path_currentFile.split("config")[0]+'examples'} 

    if os.getcwd().split(divider)[1] == 'home':
        CONFIG['PATHS']['RUNS'] = '/home/hanne/CIEMAT/RUNS'
        CONFIG['PATHS']['NEWRUNS'] = '/home/hanne/CIEMAT/NEWRUNS' 
        CONFIG['PATHS']['RUNS_SUPERCOMPUTER'] = '/marconi_work/FUA34_KINCIEMA/hanne/RUNS'
        CONFIG['PATHS']['FIGURES'] = '/home/hanne/CIEMAT/RUNS/Figures'

    if os.getcwd().split(divider)[1] == 'marconi_work' or os.getcwd().split(divider)[1] == 'marconi':
        CONFIG['PATHS']['RUNS'] = '/marconi_work/FUA34_KINCIEMA/hanne/RUNS'
        CONFIG['PATHS']['NEWRUNS'] = 'DONT USE THESE COMMANDS ON THE SUPERCOMPUTER' 
        CONFIG['PATHS']['RUNS_SUPERCOMPUTER'] = '/DONT USE THESE COMMANDS ON THE SUPERCOMPUTER'
        CONFIG['PATHS']['FIGURES'] = '/marconi_work/FUA34_KINCIEMA/hanne/RUNS/Figures'
    
    # Write the configuration file
    CONFIG.write(open(NAME_CONFIGURATIONFILE, 'w'))
    
def write_configurationFile():
    ''' Write the configuration file. '''
    with open(NAME_CONFIGURATIONFILE, 'w') as configfile:
        CONFIG.write(configfile)
    
def check_pathsToCodeAndSimulations(root):
    ''' Check whether the paths in the configuration file point to the code. '''
    read_configurationFile(root)
    return

def turnOffVerboseWrapper_configurationFile(status="False"):
    '''  Turn the verbose_wrapper off in the configuration file. '''
    CONFIG['DEFAULT']['use_verboseWrapper'] = status
    write_configurationFile()
        
def turnOnVerboseWrapper_configurationFile():
    '''  Turn the verbose_wrapper on in the configuration file. '''
    turnOffVerboseWrapper_configurationFile(status="True")

def display_information(root,title=None,message=None):
    ''' Display a window with information.

    Parameters
    ----------
    self : class object
        self.root should be linked to the root window in order to center the window on the screen.
    message : str
        Message to be displayed.
    '''

    if not title:   title="Please insert a title"
    if not message: message="Please insert a message."
    information_window = tk.Toplevel(bg=root.color['bg'])
    information_window.title(title)
    root.eval(f'tk::PlaceWindow {str(information_window)} center')
    return
 
 
#===============================================
# MAIN SCRIPT
#===============================================
read_configurationFile()


        
