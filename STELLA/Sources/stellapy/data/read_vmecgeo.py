

import os, h5py
from stellapy.utils import get_filesInFolder
from stellapy.utils.decorators import verbose_wrapper, printv

#===========================================
#  READ THE VMEC_GEO FILE
#===========================================

@verbose_wrapper
def read_vmecgeo(input_file):
    ''' Read the "*.vmec_geo" file and return the a_ref and B_ref. '''
    
    # Find the omega file corresponding to the input file
    vmec_file = input_file.parent / "vmec_geo"
    geometry_file = input_file.with_suffix(".geometry")

    # Read the *vmec_geo file and get the first two lines (header + values)
    if os.path.isfile(vmec_file):
        geo_data = read_vmecgeoFromVmecgeoFile(vmec_file)

    # Read the *.geometry file and get the first two lines (header + values)
    elif os.path.isfile(geometry_file):
        geo_data = read_vmecgeoFromGeometryFile(geometry_file)

    # Data is saved to the reduced "wout*.h5" file
    elif get_filesInFolder(input_file.parent, start="wout", end=".h5"): 
        geo_data = read_vmecgeoFromH5File(input_file)

    # ERROR IF BOTH FILES DON'T EXIST
    else: printv("No '*.vmec_geo' or '*.geometry' or 'wout*.h5' file found."); return {}
    
    # Return the data
    return geo_data  


#-------------------------
def read_vmecgeoFromVmecgeoFile(vmec_file):
    
    # Read the file
    vmec_file   = open(vmec_file, 'r')
    header_line = vmec_file.readline() # @UnusedVariable
    first_line  = vmec_file.readline().split()
    vmec_file.close()

    # Safe the data to a dictionary:  [rhotor; qinp; shat; aref; Bref; z_scalefac]
    geo_data               = {}
    geo_data['rhotor']     = float(first_line[0])
    geo_data['qinp']       = float(first_line[1])
    geo_data['shat']       = float(first_line[2])
    geo_data['aref']       = float(first_line[3])
    geo_data['Bref']       = float(first_line[4])    # B_ref*(pi*a^2)=toroidal_flux_for_the_last_closed_flux_surface_over_two_p
    geo_data['z_scalefac'] = float(first_line[5])  
    return geo_data

#----------------------------
def read_vmecgeoFromGeometryFile(geometry_file):
    
    # Read the file
    vmec_file   = open(geometry_file, 'r')
    header_line = vmec_file.readline()              # @UnusedVariable
    first_line  = vmec_file.readline().split()
    vmec_file.close()

    # Safe the data to a dictionary  [#; rhoc; qinp; shat; rhotor; aref; bref; dxdpsi; dydalpha]
    geo_data               = {}
    geo_data['rhoc']       = float(first_line[1])
    geo_data['qinp']       = float(first_line[2])
    geo_data['shat']       = float(first_line[3])
    geo_data['rhotor']     = float(first_line[4])
    geo_data['aref']       = float(first_line[5])
    geo_data['Bref']       = float(first_line[6])    # B_ref*(pi*a^2)=toroidal_flux_for_the_last_closed_flux_surface_over_two_p
    geo_data['dxdpsi']     = float(first_line[7])  
    geo_data['dydalpha']   = float(first_line[8])
    return geo_data

#------------------------------
def read_vmecgeoFromH5File(input_file):

    # Read the file
    vmec_path = input_file.parent / get_filesInFolder(input_file.parent, start="wout", end=".h5")[0]
    h5_vmec = h5py.File(vmec_path, 'r') 

    # Safe the data to a dictionary
    geo_data = {}

    # Add the wout dimensions to the h5 file
    for quant in ['rhotor', 'qinp', 'shat', 'aref', 'Bref', 'z_scalefac', 'rhoc', 'dxdpsi', 'dydalpha']:
        try:
            geo_data[quant] = h5_vmec.attrs[quant]
        except:
            geo_data[quant] = None

    # Close the file
    h5_vmec.close()
    return geo_data
    if True: return 

#====================================================
#  ATTACH THE VMECGEO DATA TO THE SIMULATION OBJECT
#====================================================

def get_vmecgeoData(self):
    vmec = read_vmecgeo(self.input_files[0])
    self.ref_a = vmec['aref']
    self.ref_B = vmec['Bref']