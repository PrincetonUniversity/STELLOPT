
# TODO: To save memory, save only unique time and z axis in dictionaries and use references

import numpy as np
import h5py, copy, os
from scipy.io import netcdf as scnetcdf

#===========================================
#  READ THE NETCDF FILE
#===========================================

def read_netcdf(netcdf_file, data=["dimensions", "time", "z", "phi2", "fluxes", "density"]):
    ''' Read the "*out.nc" or "*out.h5" file of the stella simulation.
    
    Parameters
    ----------
    data : ["dimensions", "standard", "fluxes", "density"]
        Split the reading up in multiple parts since most plotting functions only 
        require a limited set of data, therefore save memory and processing time.
    '''
    
    # Check whether the netcdf_file is indeed the netcdf_file or the reduced netcdf file
    if   os.path.isfile(netcdf_file) and netcdf_file.suffix==".nc": pass
    elif os.path.isfile(netcdf_file) and netcdf_file.suffix==".h5": pass
    
    # Otherwise, turn the input file into a netcdf file
    elif os.path.isfile(netcdf_file.with_suffix('.out.h5')): netcdf_file = netcdf_file.with_suffix('.out.h5')
    elif os.path.isfile(netcdf_file.with_suffix('.out.nc')): netcdf_file = netcdf_file.with_suffix('.out.nc')
    else: return 
        
    # READ "OUT.NC" file 
    if netcdf_file.suffix =='.nc':
        return read_netcdfFromNcFile(netcdf_file, data)

    # READ "OUT.h5" file 
    elif netcdf_file.suffix =='.h5':
        return read_netcdfFromH5File(netcdf_file, data)
    
    return

#-------------------------------------
def read_netcdfFromNcFile(netcdf_file, data):
    ''' Read the "*out.nc" file of the stella simulation.'''
    
    # Open the "out.nc" file
    netcdf_file = scnetcdf.netcdf_file(netcdf_file,'r')

    # Safe the data to a dictionary, use deep copies, otherwise the netcdf file refuses to close since mmap=True
    netcdf_data = {}

    # Read the dimensions of the axes as well as the values of the modes
    if "dimensions" in data:
        netcdf_data['dim_species'] = copy.deepcopy(np.copy(netcdf_file.dimensions['species']))
        netcdf_data['dim_kx']      = copy.deepcopy(np.copy(netcdf_file.dimensions['kx']))
        netcdf_data['dim_ky']      = copy.deepcopy(np.copy(netcdf_file.dimensions['ky']))
        netcdf_data['dim_z']       = copy.deepcopy(np.copy(netcdf_file.dimensions['zed']))
        netcdf_data['vec_kx']      = copy.deepcopy(np.copy(netcdf_file.variables['kx'][:]))
        netcdf_data['vec_ky']      = copy.deepcopy(np.copy(netcdf_file.variables['ky'][:]))
        netcdf_data['dim_time']    = copy.deepcopy(np.copy(netcdf_file.dimensions['t']))
        if netcdf_data['dim_time']==None:
            vec_time = copy.deepcopy(np.copy(netcdf_file.variables['t'][:]))
            netcdf_data['dim_time'] = len(vec_time)
            #print("WARNING: There was a bug in dim_time variable in the netcdf file.")
            
    # Get the time axis
    if "time" in data:
        netcdf_data['vec_time'] = copy.deepcopy(np.copy(netcdf_file.variables['t'][:]))
        
    # Get the z axis
    if "z" in data:
        netcdf_data['vec_z'] = copy.deepcopy(np.copy(netcdf_file.variables['zed'][:])) 
        
    # Get the electrostatic potential averaged over z as function of (t,kx,ky)
    if "phi2" in data:
        netcdf_data['phi2_vs_kxky'] = copy.deepcopy(np.copy(netcdf_file.variables['phi2_vs_kxky'][:]))
    
    # Get the fluxes averaged over z as function of (ky,kx,t)
    if "fluxes" in data:
        for key in ['qflx_kxky', 'pflx_kxky', 'vflx_kxky']:
            try:    netcdf_data[key] = copy.deepcopy(np.copy(netcdf_file.variables[key][:]))                        
            except: netcdf_data[key] = np.NaN 
    
    # Get the density fluctuation (dim_species,kx,ky)?
    if "density" in data:
        netcdf_data['vec_density'] = copy.deepcopy(np.copy(netcdf_file.variables['dens'][:])*1e19)
    
    # Close the file
    netcdf_file.close()
    
    # Return the data read from the netcdf file
    return netcdf_data

#-------------------------------------
def read_netcdfFromH5File(netcdf_file, data):
    ''' Read the "*out.h5" file of the stella simulation.'''
    
    # Open the "out.h5" file
    h5_netcdf = h5py.File(netcdf_file, 'r')

    # Safe the data to a dictionary
    netcdf_data = {}
    
    # Read the dimensions of the axes as well as the values of the modes
    if "dimensions" in data:
        for quant in ['dim_species', 'dim_kx', 'dim_ky', 'dim_time', 'dim_z']:
            netcdf_data[quant] = h5_netcdf.attrs[quant]
        for quant in ['vec_kx', 'vec_ky']:
            netcdf_data[quant] = h5_netcdf[quant][:]
    
    # Get the time axis
    if "time" in data:    
        netcdf_data['vec_time'] = h5_netcdf['vec_time'][:]
        
    # Get the z axis
    if "z" in data:
        netcdf_data['vec_z'] = h5_netcdf['vec_z'][:]
        
    # Get the electrostatic potential averaged over z as function of (t,kx,ky)
    if "phi2" in data:
        netcdf_data['phi2_vs_kxky'] = h5_netcdf['phi2_vs_kxky'][:]
            
    # Get the fluxes averaged over z as function of (ky,kx,t)
    if "fluxes" in data:
        for quant in ['qflx_kxky', 'pflx_kxky', 'vflx_kxky']:
            try:    netcdf_data[quant] = h5_netcdf[quant][:]
            except: netcdf_data[quant] = np.array([np.NaN])
    
    # Get the density fluctuation (dim_species,kx,ky)?
    if "density" in data:
        try:    netcdf_data['vec_density'] = h5_netcdf['vec_density'][:]
        except: netcdf_data['vec_density'] =  np.NaN
        
    # Close the file
    h5_netcdf.close()
    
    # Return the data read from the netcdf file
    return netcdf_data
    if True: return
           
#==================================================
#  ATTACH THE NETCDF DATA TO THE SIMULATION OBJECT
#==================================================
 
def get_netcdfDimensions(self):
    """ Add the kx, ky, time and z dimensions as attributes to the simulation object."""
    
    # Initiate the attributes
    self.vec_kx, self.vec_ky = [], []
    self.vec_kxPerFile, self.vec_kyPerFile = {}, {}
    self.dim_species, self.dim_time, self.dim_z, self.dim_kx, self.dim_ky = 0, 0, 0, 0, 0
    self.dim_timePerFile, self.dim_zPerFile, self.dim_kxPerFile, self.dim_kyPerFile = {}, {}, {}, {}
     
    # Read the data for each input file that might differ between the input files
    # Each input file usually contains different modes (kx,ky)
    # But the time and spatial resolution might also be changed giving different vec_time and vec_z
    # Therefore dim_time is the maximum dimension among the input files, to construct the arrays
    # Remove the zero from the time axis since all output files start counting from delta t      
    for input_file in self.input_files:
        
        # Show the reading progress
        if self.Progress: i = self.input_files.index(input_file); length = len(self.input_files)
        if self.Progress: self.Progress.move(i/length*100,"Determining dimensions ("+str(i)+"/"+str(length)+")")
        
        # Read the netcdf data to determine the maximum dimensions of the input files
        netcdf_data       = read_netcdf(input_file, ["dimensions"])           
        self.dim_species  = np.max([netcdf_data['dim_species'], self.dim_species])    
        self.dim_time     = np.max([netcdf_data['dim_time']-1, self.dim_time])           
        self.dim_z        = np.max([netcdf_data['dim_z'], self.dim_z])     
        self.vec_kx      += list(netcdf_data['vec_kx'])
        self.vec_ky      += list(netcdf_data['vec_ky'])
        
        # For reading other output files we need the exact dimensions of each input file      
        self.dim_timePerFile[input_file] = netcdf_data['dim_time']-1    
        self.dim_zPerFile[input_file]    = netcdf_data['dim_z']  
        self.dim_kxPerFile[input_file]   = len(list(netcdf_data['vec_kx']))  
        self.dim_kyPerFile[input_file]   = len(list(netcdf_data['vec_ky']))  
        self.vec_kxPerFile[input_file]   = list(netcdf_data['vec_kx'])
        self.vec_kyPerFile[input_file]   = list(netcdf_data['vec_ky'])
        del netcdf_data
        
    # Get the total amount of modes and sort them
    self.vec_kx = list(set(self.vec_kx));   self.vec_kx.sort()
    self.vec_ky = list(set(self.vec_ky));   self.vec_ky.sort()
    
    # Get the dimensions of the modes
    self.dim_kx = len(self.vec_kx)
    self.dim_ky = len(self.vec_ky)
    return

#-----------------------------
def get_netcdfTimeVsKxKy(self):
    """ Add the time axis for each mode (kx,ky) of a linear simulation as an attribute of the simulation object. """
    
    # Initiate the attribute
    self.time_kxky = np.empty((self.dim_time, self.dim_kx, self.dim_ky)); self.time_kxky[:,:,:] = np.NaN
    
    # Read the data for each input file
    for input_file in self.input_files:
        
        # Show the reading progress
        if self.Progress: i = self.input_files.index(input_file); length = len(self.input_files)
        if self.Progress: self.Progress.move(i/length*100,"Reading the time axis ("+str(i)+"/"+str(length)+")")
        
        # Read the netcdf data to get time, z and phi2
        netcdf_data = read_netcdf(input_file, ["time"])         
        dim_time = self.dim_timePerFile[input_file]
        for kx in self.vec_kxPerFile[input_file]:
            for ky in self.vec_kyPerFile[input_file]:
                index_kx = self.vec_kx.index(kx)
                index_ky = self.vec_ky.index(ky) 
                self.time_kxky[0:dim_time, index_kx, index_ky] = netcdf_data['vec_time'][1:]   
        del netcdf_data
    return

#-----------------------------
def get_netcdfZVsKxKy(self):
    """ Add the z axis for each mode (kx,ky) as an attribute of the simulation object. """
    
    # Initiate the attributes: the time, z and phi2 can differ for each input file
    self.z_kxky = np.empty((self.dim_z, self.dim_kx, self.dim_ky));    self.z_kxky[:,:,:]    = np.NaN
    
    # Read the data for each input file
    for input_file in self.input_files:
        
        # Show the reading progress
        if self.Progress: i = self.input_files.index(input_file); length = len(self.input_files)
        if self.Progress: self.Progress.move(i/length*100,"Reading the z axis ("+str(i)+"/"+str(length)+")")
        
        # Read the netcdf data to get time, z and phi2
        netcdf_data = read_netcdf(input_file, ["z"])        
        dim_z    = self.dim_zPerFile[input_file]
        for kx in self.vec_kxPerFile[input_file]:
            for ky in self.vec_kyPerFile[input_file]:
                index_kx = self.vec_kx.index(kx)
                index_ky = self.vec_ky.index(ky) 
                self.z_kxky[0:dim_z, index_kx, index_ky] = netcdf_data['vec_z']   
        del netcdf_data 
    return

#-----------------------------
def get_netcdfPhi2VsKxKy(self):
    """ Add the potential squared for each mode (kx,ky) as an attribute of the simulation object. """
    
    # Initiate the attributes: the time, z and phi2 can differ for each input file
    self.phi2_kxky = np.empty((self.dim_time, self.dim_kx, self.dim_ky)); self.phi2_kxky[:,:,:] = np.NaN
    
    # Read the data for each input file
    for input_file in self.input_files:
        
        # Show the reading progress
        if self.Progress: i = self.input_files.index(input_file); length = len(self.input_files)
        if self.Progress: self.Progress.move(i/length*100,"Reading the potential ("+str(i)+"/"+str(length)+")")
        
        # Read the netcdf data to get time, z and phi2
        netcdf_data = read_netcdf(input_file, ["phi2"])         
        dim_time = self.dim_timePerFile[input_file]
        for kx in self.vec_kxPerFile[input_file]:
            for ky in self.vec_kyPerFile[input_file]:
                i_kx     = self.vec_kxPerFile[input_file].index(kx)
                i_ky     = self.vec_kyPerFile[input_file].index(ky)
                index_kx = self.vec_kx.index(kx)
                index_ky = self.vec_ky.index(ky) 
                self.phi2_kxky[0:dim_time, index_kx, index_ky] = netcdf_data['phi2_vs_kxky'][1:, i_kx, i_ky]
        del netcdf_data
    return

#---------------------------- 
def get_netcdfFluxesKxKy(self):
    """ Read the fluxes versus (kx,ky) from the netcdf file of the simulation. """
    
    # Initiate the attributes: Save the fluxes per (kx,ky)
    self.pflx_kxky = np.empty((self.dim_time, self.dim_kx, self.dim_ky)); self.pflx_kxky[:,:,:] = np.NaN
    self.qflx_kxky = np.empty((self.dim_time, self.dim_kx, self.dim_ky)); self.qflx_kxky[:,:,:] = np.NaN
    self.vflx_kxky = np.empty((self.dim_time, self.dim_kx, self.dim_ky)); self.vflx_kxky[:,:,:] = np.NaN
    
    # Read the data for each input file
    for input_file in self.input_files:
        
        # Show the reading progress
        if self.Progress: i = self.input_files.index(input_file); length = len(self.input_files)
        if self.Progress: self.Progress.move(i/length*100,"Reading the fluxes ("+str(i)+"/"+str(length)+")")
        
        # Read the netcdf data to get the fluxes
        netcdf_data = read_netcdf(input_file, ["fluxes"])         
        dim_time = self.dim_timePerFile[input_file]   
        for kx in self.vec_kxPerFile[input_file]:
            for ky in self.vec_kyPerFile[input_file]:
                i_kx     = self.vec_kxPerFile[input_file].index(kx)
                i_ky     = self.vec_kyPerFile[input_file].index(ky)
                index_kx = self.vec_kx.index(kx)
                index_ky = self.vec_ky.index(ky) 
                self.pflx_kxky[0:dim_time, index_kx, index_ky] = netcdf_data['pflx_kxky'][:, i_kx, i_ky]
                self.qflx_kxky[0:dim_time, index_kx, index_ky] = netcdf_data['qflx_kxky'][:, i_kx, i_ky]
                self.vflx_kxky[0:dim_time, index_kx, index_ky] = netcdf_data['vflx_kxky'][:, i_kx, i_ky] 
        del netcdf_data
    return 

#---------------------------
def get_netcdfDensity(self):
    """ Read the density versus (species, kx,ky)? from the netcdf file of the simulation. """    
    # TODO: Not sure what the dimensions are of the density variable, thus this could be wrong.
    
    # Initiate the attributes: Save the density per (kx,ky)
    self.density_kxky = np.empty((self.dim_species, self.dim_kx, self.dim_ky)); self.density_kxky[:,:,:] = np.NaN
    
    # Read the data for each input file
    for input_file in self.input_files:   
        
        # Read the netcdf data to get density(time, kx, ky)
        netcdf_data = read_netcdf(input_file, ["density"])         
        for kx in self.vec_kxPerFile[input_file]:
            for ky in self.vec_kyPerFile[input_file]:
                i_kx     = self.vec_kxPerFile[input_file].index(kx)
                i_ky     = self.vec_kyPerFile[input_file].index(ky)
                index_kx = self.vec_kx.index(kx)
                index_ky = self.vec_ky.index(ky) 
                self.density_kxky[:, index_kx, index_ky] = netcdf_data['vec_density'][:, i_kx, i_ky]
        del netcdf_data
    return
