

import os, h5py
import numpy as np
from stellapy.data.read_wout import read_wout
from stellapy.data.read_vmecgeo import read_vmecgeo
from stellapy.utils.get_filesInFolder import get_filesInFolder 
from stellapy.utils.decorators import verbose_wrapper
    
@verbose_wrapper # stellapy.utils.wout
def reduce_sizeWout(folder, write_again=True):
    ''' Write a wout...h5 file based on the wout...nc file with a reduced selection of data.  '''    
    
    # Make sure we have a folder
    if folder is None:
        return
    
    # Read the "wout*" file
    wout_files = get_filesInFolder(folder, start="wout", end=".nc")   

    # Save the data of each "*.out.nc" file to a "*.out.h5" file
    if wout_files:
        for wout_file in wout_files:
            file_path = wout_file.with_suffix('.h5')
            if not os.path.isfile(file_path) or write_again:
                h5_wout = h5py.File(file_path, 'w')
    
                # ---------------------------
                # Read the relevant wout data
                wout_data = read_wout(wout_file)
    
                # Add the wout vectors to the h5 file
                for quant in ['iotas', 'phi', 'iotaf']:
                    try: h5_wout.create_dataset(quant, data=wout_data[quant])	
                    except: pass
    
                # Add the wout dimensions to the h5 file
                for quant in ['ns', 'b0', 'Aminor_p', 'iota']:
                    try: h5_wout.attrs.create(quant, data=wout_data[quant])	
                    except: pass
    
                # ------------------------------
                # Read the relevant vmec_geo data
                geo_data = read_vmecgeo(wout_file)
    
                # Add the vmecgeo dimensions to the h5 file
                for quant in ['rhotor', 'qinp', 'shat', 'aref', 'Bref', 'z_scalefac', 'rhoc', 'dxdpsi', 'dydalpha']:
                    try:
                        h5_wout.attrs.create(quant, data=geo_data[quant])	
                    except:
                        h5_wout.attrs.create(quant, data=np.string_("None"))
                # Close the file
                h5_wout.close()
    
                print("    Woutnc file is saved as "+file_path.name)

    return 


