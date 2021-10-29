
import numpy as np
import scipy.interpolate
from stellapy.utils.decorators import verbose_wrapper  

@verbose_wrapper
def write_profile(
        # Data of the profiles
        raw_file, \
        # Parameter to measure distance along the cross-section of the plasma
        x='rho', a=None, \
        # Columns of the rho, density and temperature data
        x_col=None, dens_col=None, Ti_col=None, Te_col=None, \
        # The positions throughout the cross-section to calculate the profiles
        rho_values=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], \
        # Number of digits of the calculated profiles and gradients
        digits=5):
    ''' 
    Write the profiles and gradients of the electron density (n_e), electron
    temperature (T_e) and ion temperature (T_i) at specific rho values to a txt file.
    
    If the ".dat" file has a header like this:
        "rho, reff [m], ne [1e19 m-3], Te [keV], CXRS Ti [keV]"
    That specifically mentions "rho", "ne", "Te" and "Ti", the code
    can determine which columns they are in atuomatically.

    Parameters
    ----------
    folder, raw_file : str
        Specify the folder in the RUNS directory and the raw_file with extension .dat
        that contains the data from which the profiles and gradients can be calculated.
    x : {'rho', 's', 'r'}
        Specify the units of the x-axis, this will be converted into rho.
        If [x] = 'r' also specify the minor radius [a]
    x_col, dens_col, Ti_col, Te_col : int
        Specify in which columns the quantities can be found
    rho_values : list of floats 
        Give a list of rho_values to calculate n_e, T_e, T_i, f_prim, Ti_prim, Te_prim


    Example bash command
    --------------------
    write_profiles --x=0 --n=2 --Te=3 --Ti=4


    Structure of "profile.txt"
    --------------------------
         0       1     2     3      4      5     6     7        8
        rho  torflux  n_e   fprim   Ti  Ti_prim  Te  Te_prim  Ti/Te
    '''

    # Read the data file containig the profiles
    try: 
        data = np.loadtxt(raw_file, dtype='float')
    except:
        data = np.loadtxt(raw_file, dtype=np.str, delimiter='\t', skiprows=1)
        data = np.char.replace(data, ',', '.').astype(np.float64)
    
    # If the columns are not specified, try to guess
    if (x_col is None) or (dens_col is None) or (Ti_col is None) or (Te_col is None):
        x_col, dens_col, Ti_col, Te_col  = guess_columns(raw_file)
        if (x_col is None) or (dens_col is None) or (Ti_col is None) or (Te_col is None):
            print("Please specify the columns of the variables, the code is not able to determine it automatically.")
            return
    
    # Turn rho_value into a list if its a single number
    if isinstance(rho_values, float) or isinstance(rho_values, int):
        rho_values = [rho_values]

    # Change the unit along the x-axis to rho=r/a=sqrt(s)
    if x == 'r':  data[:,x_col] = data[:,x_col]/a
    if x == 's':  data[:,x_col] = np.sqrt(data[:,x_col])

    # Calculate n_e, T_e and T_i and their gradients in function of rho, r or s
    profiles = {'rho' : [], 'torflux' : [], 'n_e' : [], 'Te' : [], 'Ti' : [], \
                'f_prim' : [], 'Te_prim' : [], 'Ti_prim' : [], 'TiTe' : [], 'TeTi' : []}
    
    for rho_val in rho_values: 
        profiles['rho'].append(float(rho_val))
        profiles['torflux'].append(float(rho_val)**2)

        if dens_col:
            Y_profile_at_rho, gradient, grad_profile_at_rho = calculate_gradients(data,rho_col=x_col,Y_col=dens_col,rho_val=rho_val)[0:3]
            profiles['n_e'].append(Y_profile_at_rho)
            profiles['f_prim'].append(-grad_profile_at_rho)

        if Ti_col: 
            Y_profile_at_rho, gradient, grad_profile_at_rho = calculate_gradients(data,rho_col=x_col,Y_col=Ti_col,rho_val=rho_val)[0:3]
            profiles['Ti'].append(Y_profile_at_rho)
            profiles['Ti_prim'].append(-grad_profile_at_rho)

        if Te_col:
            Y_profile_at_rho, gradient, grad_profile_at_rho = calculate_gradients(data,rho_col=x_col,Y_col=Te_col,rho_val=rho_val)[0:3]
            profiles['Te'].append(Y_profile_at_rho)
            profiles['Te_prim'].append(-grad_profile_at_rho) 
    
        if Ti_col and Te_col:
            profiles['TiTe'].append(profiles['Ti'][-1]/profiles['Te'][-1])
            profiles['TeTi'].append(profiles['Te'][-1]/profiles['Ti'][-1])

    # If n_e, T_e or T_i is missing, replace the data with NAN
    Nan_array = np.empty((len(rho_values),))*np.NaN
    if dens_col is None:   profiles['n_e'], profiles['f_prim'] = Nan_array, Nan_array
    if Ti_col is None:     profiles['Ti'], profiles['Ti_prim'] = Nan_array, Nan_array
    if Te_col is None:     profiles['Te'], profiles['Te_prim'] = Nan_array, Nan_array

    # Round the numbers to x digits
    for key in profiles.keys():
        profiles[key] = [round(value, digits) for value in profiles[key]]

    # Write the profile data to a text file
    extension = "" if "profile_" not in str(raw_file) else str(raw_file).split("profile_")[-1].split(".")[0]
    output_path = raw_file.parent / str('profile_'+extension+'.txt')
    output_header = '  rho      torflux        n_e        fprim         Ti        Ti_prim      Te       Te_prim      Ti/Te       Te/Ti'
    output_data = (profiles['rho'], profiles['torflux'], profiles['n_e'], profiles['f_prim'], profiles['Ti'], \
                   profiles['Ti_prim'], profiles['Te'], profiles['Te_prim'], profiles['TiTe'], profiles['TeTi'])
    output_data = np.array(output_data).transpose()
    np.savetxt(output_path, output_data, delimiter='\t', header=output_header, fmt="%8."+str(digits)+"f")
    print("-----------------------------------------------")
    print('The profile data is saved to ', output_path)
    print("-----------------------------------------------")
    return x_col, dens_col, Ti_col, Te_col


#=====================================================================
# Calculate the gradients a/L_y from the profile data at a certain rho
#=====================================================================

@verbose_wrapper
def calculate_gradients(data,rho_col,Y_col,rho_val=None,smooth=0.0,steps=100):
    '''Calculation of the gradients of the profile Y in function of rho=r/a=sqrt(s)

    Parameters
    ----------
    data : 2D array
        Contains the raw profile data
    rho_col, Y_col : int
        data[:,x_col] contains the values for the profile of x
    rho_val : float
        return the value and gradient of the profile at this specific rho

    Returns
    -------
    {profile_at_loc, a_over_Ly_at_loc, Y_profile, a_over_Ly}
    
    Notes
    -----
    a/L_Y = a*dY/dr*1/Y = a*(drho/dr)*dlog(Y)/drho = dlogY/drho = dY/drho*1/Y

    Normalize this Y_N = Y/Y_ref or Y = Y_N*Y_ref or 5 = 1*5
    
    a/L_Y = dY/drho*1/Y = dY_N Y_ref/drho * 1 / (Y_N Y_ref) = dY_N / drho * 1/Y_N 
    '''

    # Get the data along rho and Y
    rho, Y              = data[:,rho_col], data[:,Y_col]

    # Length of interpolation and step size
    x_vec               = np.linspace(0,1.,steps)

    # Build the profiles by using a B-spline representation of a 1-D curve 
    tck                 = scipy.interpolate.splrep(rho,Y,s=smooth)  
    Y_profile           = scipy.interpolate.splev(x_vec,tck,der=0)                              
    length_profile      = -scipy.interpolate.splev(x_vec,tck,der=1)/Y_profile
    grad_profile        = scipy.interpolate.splev(x_vec,tck,der=1)    

    # Get the value of the profile and the gradient at a specific rho
    if rho_val is not None:
        Y_profile_at_rho      = float(scipy.interpolate.splev(rho_val,tck,der=0))
        length_profile_at_rho = -float(scipy.interpolate.splev(rho_val,tck,der=1)/Y_profile_at_rho)
        grad_profile_at_rho   = float(scipy.interpolate.splev(rho_val,tck,der=1))
        return Y_profile_at_rho, grad_profile_at_rho, length_profile_at_rho, Y_profile, grad_profile, length_profile, x_vec

    return Y_profile, grad_profile, length_profile, x_vec


#=====================================================================
# Guess the columns of the data
#=====================================================================

def guess_columns(file_path):
    ''' If the columns of the variables are not specified, try to guess them based on 
    the order of "rho", "ne", "Te" and "Ti" in the file. Assume a line like this is present: 
    "rho, reff [m], ne [1e19 m-3], Te [keV], CXRS Ti [keV]"
    '''
    
    ordered_variables = ["rho"]
    file_text = open(file_path, "r").read() 
    if len(file_text.split("rho")) > 2:
        return None, None, None, None 
    else:
        before_rho = file_text.split("rho")[0]
        after_rho  = file_text.split("rho")[1]
        if ("ne" not in before_rho) and ("ne" not in after_rho)\
        or ("ne" in before_rho) and ("ne" in after_rho):
            return None, None, None, None 
        else:
            if "ne" in before_rho: ordered_variables.insert(0, "ne")
            if "ne" in after_rho:  ordered_variables.append("ne")
            before_ne = file_text.split("ne")[0]
            after_ne  = file_text.split("ne")[1]
            if ("Te" not in before_ne) and ("Te" not in after_ne)\
            or ("Te" in before_ne) and ("Te" in after_ne):
                return None, None, None, None 
            else:
                if "Te" in before_rho and "Te" in before_ne: ordered_variables.insert(0, "Te")
                if "Te" in before_rho and "Te" in after_ne: ordered_variables.insert(1, "Te")
                if "Te" in after_rho  and "Te" in before_ne: ordered_variables.insert(1, "Te")
                if "Te" in after_rho  and "Te" in after_ne: ordered_variables.append("Te")
                before_Te = file_text.split("Te")[0]
                after_Te  = file_text.split("Te")[1]
                if ("Ti" not in before_ne) and ("Ti" not in after_ne)\
                or ("Ti" in before_ne) and ("Ti" in after_ne):
                    return None, None, None, None  
                else:
                    if "Ti" in before_rho and "Ti" in before_ne and "Ti" in before_Te: ordered_variables.insert(0, "Ti") 
                    if "Ti" in before_rho and "Ti" in before_ne and "Ti" in after_Te:  ordered_variables.insert(1, "Ti") 
                    if "Ti" in before_rho and "Ti" in after_ne and "Ti" in before_Te:  ordered_variables.insert(1, "Ti") 
                    if "Ti" in before_rho and "Ti" in after_ne and "Ti" in after_Te:   ordered_variables.insert(2, "Ti") 
                    if "Ti" in after_rho and "Ti" in before_ne and "Ti" in before_Te:  ordered_variables.insert(1, "Ti") 
                    if "Ti" in after_rho and "Ti" in before_ne and "Ti" in after_Te:   ordered_variables.insert(2, "Ti") 
                    if "Ti" in after_rho and "Ti" in after_ne and "Ti" in before_Te:   ordered_variables.insert(2, "Ti") 
                    if "Ti" in after_rho and "Ti" in after_ne and "Ti" in after_Te:    ordered_variables.append("Ti")
                    
    # There could be more columns present between the desired columns
    # "rho, reff [m], ne [1e19 m-3], Te [keV], CXRS Ti [keV]"
    header = ordered_variables[0] + file_text.split(ordered_variables[0])[-1]
    header = header.split(ordered_variables[-1])[0] + ordered_variables[-1]
    header = header.split(",")
    
    # We found the columns
    x_col    = [header.index(text) for text in header if ("rho" in text)][0]
    dens_col = [header.index(text) for text in header if ("ne" in text)][0]
    Ti_col   = [header.index(text) for text in header if ("Ti" in text)][0]
    Te_col   = [header.index(text) for text in header if ("Te" in text)][0] 
    print()
    print("The columns of the data were found automatically, but please check if it is correct:")
    print("    rho: column", x_col, "     ne: column", dens_col, "    Ti: column", Ti_col, "    Te: column", Te_col)
    print()
    return x_col, dens_col, Ti_col, Te_col 


# #===================================================================
# # Convert (kmin, kmax) from SI units to normalized units 
# #===================================================================
# 
# @verbose_wrapper
# def add_krangeToProfile(folder, profile_file='profile.txt', lineardata_file=None, k_min=8, k_max=10):
#     ''' Convert (kmin, kmax) from SI units to normalized units used in stella.
# 
#     The results are printed to "profile_krange.txt"
#     
#     Stella has to be run already once since the normalization needs B_ref and a_ref found in *vmec_geo. 
# 
#     Structure of the 'profile_krange.txt file
#          0       1     2     3      4      5     6     7        8     9     10
#         rho  torflux  n_e   fprim   Ti  Ti_prim  Te  Te_prim  Ti/Te  kmin  kmax
#     '''
#     # Look for a linear data file
#     if not lineardata_file:
#         lineardata_file = get_filesInFolder(folder, start="lineardata")[0]
# 
#     # Read the data file containig the profiles
#     file_path = get_fullPathOfFolder(folder) + '/profile.txt'
#     file_data = np.loadtxt(file_path, dtype='float').reshape(-1, 10)
#     rho_values = file_data[:,0]
# 
#     # If no (k_min, k_max) are specified, take a standard range for reflectometry
#     if k_min is None or k_max is None:
#         k_min, k_max = 8.0, 10.0
# 
#     # Initiate the krange and add values for each rho
#     krange = {'min' : np.empty(len(rho_values)), 'max' : np.empty(len(rho_values))}
# 
#     # Transfer (kmin, kmax) in SI units to normalized units
#     ref = get_referenceUnits(folder=folder, lineardata_file=lineardata_file, profile_file=profile_file)
#     for simulation in ref.keys():
#         rho = ref[simulation]['rho']
#         index = np.where(np.isclose(rho_values, rho))[0]
#         krange['min'][index] = k_min*ref['rho_i']*100 # from cm^-1 to normalized units
#         krange['max'][index] = k_max*ref['rho_i']*100 # from cm^-1 to normalized units
#         
#     # Write the profile data to a text file
#     output_path = get_fullPathOfFolder(folder) + '/profile_krange.txt'
#     output_header = '  rho      torflux        n_e        fprim         Ti        Ti_prim      Te       Te_prim' + \
#                     '      Ti/Te        kmin        kmax'
#     output_data = np.append(file_data, convert_to_2D_array(np.array(krange['min'])), axis=1)
#     output_data = np.append(output_data, convert_to_2D_array(np.array(krange['max'])), axis=1)
#     np.savetxt(output_path, output_data, delimiter='\t', header=output_header, fmt="%8.4f")
#     print("-----------------------------------------------")
#     print('The profile data with krange is saved to ', output_path)
#     print("-----------------------------------------------")
#     return
