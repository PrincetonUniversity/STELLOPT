from stellapy.config import CONFIG

def get_initialDirectory():
    '''
    Returns the initial directories to select simulations in select_files() and select_folders()
    
    
    Returns
    -------
    initialdir: dict['Files', 'Folders'] = str       
    '''

    # Start looking for files and folders in the preset directories
    example_dir1  = CONFIG['PATHS']['Default directory to browse files']
    example_dir2  = CONFIG['PATHS']['Default directory to browse folders'] 
    example_dir3  = CONFIG['PATHS']['Stellapy'] + "/examples/Profiles" 
    initialdir  = { "Files" :  example_dir1, "Folders" :  example_dir2, "Profiles" :  example_dir3}
    return initialdir
