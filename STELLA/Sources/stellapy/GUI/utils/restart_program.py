
import os, sys
from stellapy.config import turnOnVerboseWrapper_configurationFile

def restart_program():
    """Restarts the current program.
    Note: this function does not return. Any cleanup action (like
    saving data) must be done before calling this function."""
    
    # When using the command prompt, the functions have wrappers
    # This command also saves the current configuration to "config.ini"
    turnOnVerboseWrapper_configurationFile()       
    
    # Now restart the program 
    python = sys.executable
    os.execl(python, python, *sys.argv)


