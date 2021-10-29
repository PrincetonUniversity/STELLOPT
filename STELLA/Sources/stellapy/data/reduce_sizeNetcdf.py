

import os, h5py
from stellapy.data.read_netcdf import read_netcdf
from stellapy.utils.get_filesInFolder import get_filesInFolder
from stellapy.utils.decorators import verbose_wrapper


@verbose_wrapper 
def reduce_sizeNetcdf(folder, write_again=True):
    ''' Write a out.h5 file based on the out.nc file with a reduced selection of data.  '''    
    
    # Make sure we have a folder
    if folder is None:
        return
    
    # Get the netcdf files inside this folder, dont look at subfolders!
    netcdf_files = get_filesInFolder(folder, end="out.nc")
    if netcdf_files: netcdf_files = [f for f in netcdf_files if not "/" in str(f).split(str(folder)+"/")[-1]]

    # Save the data of each "*.out.nc" file to a "*.out.h5" file
    if netcdf_files:
        for netcdf_file in netcdf_files:
            file_path = netcdf_file.with_suffix(".h5")
            if not os.path.isfile(file_path) or write_again:
    
                h5_netcdf = h5py.File(file_path, 'w')
    
                # Read the relevant netcdf data
                netcdf_data = read_netcdf(netcdf_file)
    
                # Add the netcdf vectors to the h5 file
                for quant in ['vec_z', 'vec_time', 'vec_ky', 'vec_kx', 'vec_density', 'phi2_vs_kxky', \
                              'qflx_kxky', 'pflx_kxky', 'vflx_kxky', 'phi_vs_t', 'phi2']:
                    try:    h5_netcdf.create_dataset(quant, data=netcdf_data[quant])	
                    except: pass 
    
                # Add hte netcdf dimensions to the h5 file
                for quant in ['dim_species', 'dim_kx', 'dim_ky', 'dim_time', 'dim_z']:
                    h5_netcdf.attrs.create(quant, data=netcdf_data[quant])	
    
                # Close the file
                h5_netcdf.close()
    
                print("    Netcdf file is saved as "+file_path.name)             
    return 


