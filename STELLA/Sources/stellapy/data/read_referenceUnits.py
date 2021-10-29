
import numpy as np
import scipy.constants as sc
from stellapy.utils.decorators import verbose_wrapper

#===========================
# GET THE REFERENCE UNITS
#===========================

@verbose_wrapper  
def read_referenceUnits(inputParameters, ref_a, ref_B, prof_n, prof_T, verbose=False):
    ''' Return the reference values used to normalize all the stella data.
    
    Since stella runs at n=1 and T=1 the actual density and temperature need 
    to be read from the "profile.txt" file.
        
    Three possible configurations for the origin of the data
        (1) Have a 'lineardata.h5' in <folder>
        (2) Have a '*.in' and '*.vmec_geo' and perhaps a 'profile.txt' in <folder>

    Returns
    -------
    referenceUnits : dict[simulation][e, proton_m, Z, charge, temp, mass, dens, length, B, vthermal, omega, rho]
        Dictionary containing the reference units for the simulations of experiment "input_file".
    '''   
    
    # Get the input parameters that we need
    input_values = { 's1_z'     : inputParameters['species_parameters_1']['z'],\
                     's1_mass'  : inputParameters['species_parameters_1']['mass'],\
                     's1_tempI' : inputParameters['species_parameters_1']['temp'],\
                     's1_densI' : inputParameters['species_parameters_1']['dens'],\
                     's1_temp'  : prof_T,\
                     's1_dens'  : prof_n,\
                     'ref_a'    : ref_a,\
                     'ref_B'    : ref_B,\
                     'rho'      : inputParameters['vmec_parameters']['rho']}

    # Reference units that are independent of the case
    ref             = {}
    ref['e']        = sc.value('elementary charge')           # [C]
    ref['proton_m'] = sc.value('proton mass')                 # [kg]      proton mass
    
    # Save at which rho these reference units are valid
    ref['rho']      = input_values['rho']  
    
    # Units used in the input file
    unit            = {}
    unit['charge']  = ref['e']                                # [C]       when using adiabatic electrons
    unit['temp']    = 1000*ref['e']                           # [J]       1 keV  [J=kg*m²/s²]
    unit['mass']    = ref['proton_m']                         # [kg]      proton mass
    unit['dens']    = 1.E19                                   # [kg/m^3]   

    # The reference units are based on the input values
    ref['Z']        = input_values['s1_z']                      # [/]       taken from *.in (charge state of species 1)
    ref['charge']   = input_values['s1_z']*unit['charge']       # [C]       taken from *.in (charge of species 1)
    ref['temp']     = input_values['s1_temp']*unit['temp']      # [J]       taken from *.in (temp state of species 1)   
    ref['mass']     = input_values['s1_mass']*unit['mass']      # [kg]      taken from *.in (mass of species 1)
    ref['dens']     = input_values['s1_dens']*unit['dens']      # [1/m^3]   taken from *.in (density of species 1) 
    ref['length']   = input_values['ref_a']                     # [m]       taken from *.vmec_geo (effective minor radius)
    ref['B']        = input_values['ref_B']                     # [T]       taken from *.vmec_geo (volume averaged B) 

    # Calculate other reference units (cross-checked with Jose's version)
    ref['vthermal'] = np.sqrt(2*ref['temp']/ref['mass'])      # [m/s]
    ref['omega']    = ref['Z']*ref['e']*ref['B']/ref['mass']  # [1/s]     [T=kg/(C*s)]
    ref['rho_i']    = ref['vthermal']/ref['omega']            # [m]

    ref['flux_part'] = (ref['rho_i']/ref['length'])**2.0*ref['dens']*ref['vthermal']                                # [1/m^3]*[m/s] = 1/(m^2*s)
    ref['flux_mom']  = (ref['rho_i']/ref['length'])**2.0*ref['dens']*ref['length']*ref['mass']*ref['vthermal']**2.0 # [kg/m^2]*[m^2/s^2]=kg^2/(s^2)
    ref['flux_heat'] = (ref['rho_i']/ref['length'])**2.0*ref['dens']*ref['temp']*ref['vthermal']                # [1/m^3]*[m/s]*J=J/(m^2*s)=kg/s^3=W/m^2

    # Show the reference units in the command prompt
    if verbose: return ref, input_values
    
    # Return the reference units
    return ref 
    if True: return

#======================================================
#  ATTACH THE REFERENCE UNITS TO THE SIMULATION OBJECT
#======================================================
    
def get_referenceUnits(self):
    self.referenceUnits = read_referenceUnits(self.inputParameters, self.ref_a, self.ref_B, self.prof_n, self.prof_T)
    










