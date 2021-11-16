
from stellapy.config import CONFIG 

def toggle_verboseOption():
    
    # Add a header to the function with the function name  
    if CONFIG['DEFAULT']['use_verbosewrapper'] != "False":
        return True
    else:
        return False
