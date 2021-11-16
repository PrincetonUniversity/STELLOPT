import warnings
from stellapy.GUI.graph_tools.display_information import display_warning

def get_axisOfScan(kx_range, ky_range, root):
    ''' Return scan={kx,ky}; k_fixed; k_range. '''  
    
    # Make sure one of the ranges is actually a float
    if isinstance(kx_range, list) and isinstance(ky_range, list):
        if kx_range[0]==kx_range[1]: 
            kx_range=kx_range[0]
        elif ky_range[0]==ky_range[1]: 
            ky_range=ky_range[0] 
        else: 
            print("TEST", kx_range, ky_range)
            warning = "Can not scan over both kx and ky, make sure kx_range or ky_range is a float."
            warnings.warn("\n"+warning, stacklevel=2); 
            display_warning(root, warning) 
            return False, False, False
    
    # Scan over ky  
    if isinstance(kx_range, float): 
        scan = 'ky'
        k_fixed = kx_range
        k_range = ky_range
    
    # Scan over kx
    elif isinstance(kx_range, float): 
        scan = 'kx'
        k_fixed = ky_range
        k_range = kx_range
    
    # Return the data      
    return scan, k_fixed, k_range
        
    