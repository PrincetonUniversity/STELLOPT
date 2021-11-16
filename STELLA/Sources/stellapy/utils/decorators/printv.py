
from stellapy.config import CONFIG 

def printv(message):
    
    # Get the indentation for the text printed to the command prompt
    from stellapy.utils.decorators.verbose_wrapper import indent
    
    # Only print if verbose==True in the configuration file.
    if CONFIG['DEFAULT']['use_verbosewrapper'] != "False":
        print(indent, message)
        return
    else:
        return

