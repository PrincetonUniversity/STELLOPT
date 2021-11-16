
import numpy as np
import os, h5py, netCDF4
from stellapy.utils import get_filesInFolder 

#===========================================
#  READ THE WOUT FILE
#===========================================

def read_wout(folder, svalue=0.5):
    ''' Read the "wout_*.nc" file and return the variables of the magnetic field at s=svalue. '''

    # If we are in the rerun folder, the wout file is in the folder above
    if "run" in folder.name:
        if folder.name.split("run")[-1].isdigit():
            folder = folder.parent

    # Get the "wout*.nc" or "wout*.h5" file in the folder
    if os.path.isfile(folder):
        wout_file_paths = [folder]
    elif get_filesInFolder(folder, start="wout", end=".nc"): 
        wout_file_paths = get_filesInFolder(folder, start="wout", end=".nc")
    elif get_filesInFolder(folder, start="wout", end=".h5"): 
        wout_file_paths = get_filesInFolder(folder, start="wout", end=".h5")
    else:
        print("WARNING No 'wout*.nc' or 'wout*.h5' file found in" + str(folder) + "."); 
        return { 'b0' : 1, 'iota' : 1}
        
    # Go through the magnetic field files
    for wout_file_path in wout_file_paths:

        # READ "WOUT...NC" file 
        if os.path.isfile(wout_file_path) and wout_file_path.suffix=='.nc':
            wout_variables = read_woutFromNcFile(wout_file_path, svalue)
            
        # READ "WOUT...h5" file 
        elif os.path.isfile(wout_file_path) and wout_file_path.suffix=='.h5'\
        or os.path.islink(wout_file_path) and os.path.isfile(wout_file_path.with_suffix(".h5")):
            wout_variables = read_woutFromH5File(wout_file_path)
            
        # ERROR IF BOTH FILES DON'T EXIST
        else: 
            print("WARNING No 'wout*.nc' or 'wout*.h5' file found in" + wout_file_path + "."); 
            return { 'b0' : 1, 'iota' : 1}
    
    return wout_variables

#-------------------------------------
def read_woutFromNcFile(wout_file_path, svalue):
        
    # Load an extra function
    from stellapy.data.get_gridDivisionsAndSize import get_iota

    # Read the netcdf file
    dataset  = netCDF4.Dataset(wout_file_path, mode='r')

    # Safe the data to a dictionary
    wout_variables              = {} 
    wout_variables['ns']        = int(dataset.variables['ns'][:])       # number of divisions along s=(r/a)^2=rho^2 during the VMEC run
    wout_variables['b0']        = int(dataset.variables['b0'][:])       # Magnetic field averaged along the axis
    wout_variables['phi']       = dataset.variables['phi'][:]           # Toroidal flux on full mesh 
    wout_variables['Aminor_p']  = dataset.variables['Aminor_p'][:]
    wout_variables['iotas']     = dataset.variables['iotas'][:]
    wout_variables['iotaf']     = dataset.variables['iotaf'][:]         # Iota on full mesh
    wout_variables['iota']      = get_iota(None, svalue=svalue, iotas=wout_variables['iotas'], ns=wout_variables['ns'])
    return wout_variables

#------------------------------------
def read_woutFromH5File(wout_file_path):

    # If there is a fake link to a .nc file, get the h5 file
    if os.path.islink(wout_file_path):
        wout_file_path = wout_file_path.with_suffix(".h5")
    
    # Load the h5 file
    h5_wout = h5py.File(wout_file_path, 'r')

    # Safe the data to a dictionary
    wout_variables = {}

    # Read the wout vectors from the h5 file
    for quant in ['iotas', 'phi', 'iotaf']:
        try:    wout_variables[quant] = h5_wout[quant][:]
        except: wout_variables[quant] = None
            
    # Add the wout dimensions to the h5 file
    for quant in ['ns', 'b0', 'Aminor_p', 'iota']:
        try:    wout_variables[quant] = h5_wout.attrs[quant]
        except: wout_variables[quant] = None
    
    # Close the file
    h5_wout.close()
    return wout_variables
    if True: return 

#==================================================
#  ATTACH THE WOUT DATA TO THE SIMULATION OBJECT
#==================================================

def get_woutData(self):
    """ Attach the sign of the magnetic field to the simulation object. """
    woutParameters = read_wout(self.input_files[0].parent)
    try: self.sign_B = np.sign(woutParameters['b0'])
    except: self.sign_B = np.sign(1)
    self.iota   = woutParameters['iota']
    del woutParameters
    




