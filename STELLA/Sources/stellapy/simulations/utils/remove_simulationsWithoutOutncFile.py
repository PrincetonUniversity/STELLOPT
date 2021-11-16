
#===============================================================
# Remove cases without an *.out.nc file
#===============================================================

import os.path  
from stellapy.utils.decorators import printv, verbose_wrapper 

@verbose_wrapper 
def remove_simulationsWithoutOutncFile(input_files):
    ''' Analysis can only be performed if *out.nc file or *out.h5 file is present.

    Returns
    -------       
    input_files : list of PosixPaths
        List of input_files of which the existence of *out.nc or *out.h5 files is checked.
    '''

    # Select the cases without an *out.nc file
    printv("Checking existence of *out.nc file corresponding inputs.")
    count=0; not_valid=[]
    for input_file in input_files:
        if os.path.isfile(input_file.with_suffix('.out.nc')):     
            count = count + 1  
        elif os.path.isfile(input_file.with_suffix('.out.h5')):                 
            count = count + 1
        else:
            printv("    File "+input_file.name+" removed from input list.")
            not_valid.append(input_file)
    printv("    Number of found input files with corresponding output *.nc = "+str(count))  

    # Remove the selected cases
    for input_file in not_valid:
        i = input_files.index(input_file)
        input_files.pop(i) 
        
    return input_files
