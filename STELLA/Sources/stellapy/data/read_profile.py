
#===============================================================
# Read variables from the profile.txt file
#===============================================================

import sys, math
import numpy as np 
from stellapy.utils import get_filesInFolder 
from stellapy.utils.decorators import verbose_wrapper

#===========================================
#  READ THE PROFILE FILE
#===========================================

@verbose_wrapper
def read_profile(folder=None, rho=None, file_identifier=None):
    ''' Read the 'profile.txt' file and return profile_data = {s1_dens; s1_temp}
    
    Returns
    -------
    profile_data : dict[s1_dens, s1_temp]

    Notes   
    -----
    Structure of the 'Profile_Data.txt'' file
         0       1     2     3      4      5     6     7        8
        rho  torflux  n_e   fprim   Ti  Ti_prim  Te  Te_prim  Ti/Te
    '''
    
    # Read files inside folder
    profile_files = get_filesInFolder(folder, start="profile")
    profile_files = [ profile_file for profile_file in profile_files if ("krange" not in profile_file) ]
    if len(profile_files)>1 and not file_identifier:
        print('Multiple "profile" files were found and not specified which one to use.', read_profile, sys._getframe().f_lineno) 
        return
    elif len(profile_files)>1 and file_identifier:
        for profile in profile_files:
            if str(file_identifier) in profile:
                profile_file = profile
                print("Read reference units for file_identifier", file_identifier, "from", profile)
    elif len(profile_files)==0:
        profile_data = {'s1_dens' : 1, 's1_temp' : 1}
        return profile_data 
    else:
        profile_file = profile_files[0]

    # Read the 'profile.txt' file 
    try:    file_data = np.loadtxt(profile_file, dtype='float').reshape(-1, 10)
    except: file_data = np.loadtxt(profile_file, dtype='float').reshape(-1, 9)
    
    # Find the profile at the given rho value
    index_rho = np.where([ math.isclose(float(rho_vec), float(rho), rel_tol=1e-2)  for rho_vec in file_data[:,0] ])[0]
    if len(list(index_rho))==0:
        print("The profile.txt file does not have profile data for rho="+str(rho), read_profile, sys._getframe().f_lineno) 
        return
    
    # Return the profile density and temperature at rho=rho
    return {'s1_dens' : float(file_data[index_rho,2]),\
            's1_temp' : float(file_data[index_rho,4])} 
    
    # Dont collapse header
    if True: return

#====================================================
#  ATTACH THE PROFILE DATA TO THE SIMULATION OBJECT
#====================================================

def get_profileData(self):
    
    # Check if there are profile files
    profile_file = get_filesInFolder(self.input_files[0].parent, start="profile", end="txt")
    
    # If there are, read the density and temperature at rho=rho
    if profile_file: 
        profile_file=profile_file[0]
        profile = read_profile(self.input_files[0].parent, rho=self.inputParameters["vmec_parameters"]["rho"])
        self.prof_n = profile['s1_dens']
        self.prof_T = profile['s1_temp']
        
    # Otherwise use the values of the input file
    else:
        self.prof_n = self.inputParameters['species_parameters_1']['dens']
        self.prof_T = self.inputParameters['species_parameters_1']['temp']
        if self.prof_n==1 and self.prof_T==1:
            print('\nWARNING: reference units are calculated with n=1 and T=1 in stella and \
                   \n         no "profile.txt" was present to impose real values for the density and temperature.') 
    return 
