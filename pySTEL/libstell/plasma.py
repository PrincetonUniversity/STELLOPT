"""
This library provides a python class for multi-species plasmas
"""
import sys

# Constants
EC = 1.602176634E-19 # Electron charge [C]
DA = 1.66053906660E-27 # Dalton
ME = 9.1093837E-31 # Electron mass [kg]

import numpy as np

class PLASMA:
    
    def __init__(self,list_of_species=None):
        
        self.species_database = ['electrons','hydrogen','deuterium','tritium','helium3','helium4']
        self.mass_database = {
            'electrons' : ME,
            'hydrogen'  : 1.007276466621*DA,
            'deuterium' : 2.01410177811*DA,
            'tritium'   : 3.01604928*DA,
            'helium3'   : 3.0160293*DA,
            'helium4'   : 4.002603254*DA
        }
        self.charge_database = {
            'electrons' : -EC,
            'hydrogen'  : EC,
            'deuterium' : EC,
            'tritium'   : EC,
            'helium3'   : 2.0*EC,
            'helium4'   : 2.0*EC
        }
        self.Zcharge_database = {
            'electrons' : -1,
            'hydrogen'  : 1,
            'deuterium' : 1,
            'tritium'   : 1,
            'helium3'   : 2,
            'helium4'   : 2
        }
        
        self.mass = {}
        self.charge = {}
        self.Zcharge = {}
        self.density = {}
        self.temperature = {}
        
        if list_of_species is not None:    
            self.check_species_exist(list_of_species)
            self.list_of_species = list_of_species
            
            self.give_mass_to_species(list_of_species)
            self.give_charge_to_species(list_of_species)
            self.give_Zcharge_to_species(list_of_species)
            
            print(f'Plasma created with species: {", ".join(self.list_of_species)}')
        else:
            print('\nPlasma created without species. Please add species using the function add_species(...)')
            self.list_of_species = []       
    
    def check_species_exist(self, species_list):
        # Ensure species_list is a list of strings
        if isinstance(species_list, str):
            species_list = [species_list]  # Convert single string to list
        elif not isinstance(species_list, list) or not all(isinstance(item, str) for item in species_list):
            raise ValueError("Species_list must be a string or a list of strings.")
        
        # Check if species are in database
        for species in species_list:
            if species not in self.species_database:
                print(f"ERROR: Species {species} does not exist in the database. Leaving the program.")
                exit(1)  # Exit the program with a status code of 1 (indicating error)
                
    def give_mass_to_species(self,species_list):
        
        for species in species_list:
            self.mass[species]  = self.mass_database[f'{species}'] 
                
    def give_charge_to_species(self,species_list,charge=None):
        
        for species in species_list:
            self.charge[species]  = self.charge_database[f'{species}']
    
    def give_Zcharge_to_species(self,species_list,Zcharge=None):
        
        for species in species_list:
            self.Zcharge[species]  = self.Zcharge_database[f'{species}'] 
                
    def set_density(self,species,n0,nedge,exponent):
        
        #check if species exist in list_of_species
        if species not in self.list_of_species:
            print(f"ERROR: Species {species} is not in the plasma.")
            exit(1)
            
        # Check if species exists in the dictionary, if not, create an empty dictionary for it
        if species not in self.density:
            self.density[species] = {}
        
        # Set the value for the specific location
        profile_info = ['n0','nedge','exponent']
        profile_vals = [n0,nedge,exponent]
        for info,val in zip(profile_info,profile_vals):
            self.density[species][info] = val
            
        print(f'\nDensity profile of {species}: n[m-3] = {nedge} + {n0-nedge}*(1-rho^{exponent})')
        
    def set_temperature(self,species,T0,Tedge,exponent):
        
        #check if species exist in list_of_species
        if species not in self.list_of_species:
            print(f"ERROR: Species {species} is not in the plasma.")
            exit(1)
            
        # Check if species exists in the dictionary, if not, create an empty dictionary for it
        if species not in self.temperature:
            self.temperature[species] = {}
        
        # Set the value for the specific location
        profile_info = ['T0','Tedge','exponent']
        profile_vals = [T0,Tedge,exponent]
        for info,val in zip(profile_info,profile_vals):
            self.temperature[species][info] = val
            
        print(f'\nTemperature profile of {species}: T[eV] = {Tedge} + {T0-Tedge}*(1-rho^{exponent})')
    
    def get_density(self,species,rho):
        # rho can be a number or a list of numbers
        
        #check if species exist in list_of_species
        if species not in self.list_of_species:
            print(f"ERROR: Species {species} is not in the plasma.")
            exit(1)
        
        n0 = self.density[species]['n0']
        nedge = self.density[species]['nedge']
        exponent = self.density[species]['exponent']
        
        rho = np.array(rho)
        
        dens = nedge + (n0-nedge)*(1-rho**exponent)
        
        return dens
    
    def get_temperature(self,species,rho):
        # rho can be a number or a list of numbers
        
        #check if species exist in list_of_species
        if species not in self.list_of_species:
            print(f"ERROR: Species {species} is not in the plasma.")
            exit(1)
        
        T0 = self.temperature[species]['T0']
        Tedge = self.temperature[species]['Tedge']
        exponent = self.temperature[species]['exponent']
        
        rho = np.array(rho)
        
        temp = Tedge + (T0-Tedge)*(1-rho**exponent)
        
        return temp
    
    def get_thermal_speed(self,species,rho):
        
        temp = self.get_temperature(species,rho)
        vth = np.sqrt(2*EC*temp/self.mass[species])
        
        return vth
    
    def get_collisionality(self,species,rho,vtest=None):
        # computes collisionality of test particle of 'species' 
        # in a termal bath of all the other species (including itself)
        # The function used is: PENTA collisionality
        
        # vtest and rho must have the same size
        # if vtest is not provided, it is assumed that vtest=vth
        
        from collisions import COLLISIONS

        coll = COLLISIONS()
        
        #if vtest is not given, assume thermal speed
        if vtest is not None:
            rho_a = np.atleast_1d(rho)    
            vtest_a = np.atleast_1d(vtest) 
            
            #check vtest and rho are of the same size
            if vtest_a.shape != rho_a.shape:
                print('ERROR: vtest must have the same dimension as rho')
                exit(1)          
        else:
            rho_a = np.atleast_1d(rho) 
            vtest_a = self.get_thermal_speed(species,rho_a)
        
        # collisionfreq_PENTA does not accept arrays
        # so need to make a loop in rho
        nu_D = []
        for r,vt in zip(rho_a,vtest_a):
            
            # arrays are organized as: first element corresponds to species we want the collisionality
            # The order of the others are arbitrary
            
            m = np.array([self.mass[species]] + [self.mass[sp] for sp in self.list_of_species if sp != species])
            Z = np.array([self.Zcharge[species]] + [self.Zcharge[sp] for sp in self.list_of_species if sp != species])
            T = np.array([self.get_temperature(species,r)] + [self.get_temperature(sp,r) for sp in self.list_of_species if sp != species])
            n = np.array([self.get_density(species,r)] + [self.get_density(sp,r) for sp in self.list_of_species if sp != species])
            
            #compute loglambda as in PENTA
            Te = self.get_temperature('electrons',r)
            ne = self.get_temperature('electrons',r)
            if(Te>50):
                loglambda = 25.3 - 1.15*np.log10(ne/1e6) + 2.3*np.log10(Te)
            else:
                loglambda = 23.4 - 1.15*np.log10(ne/1e6) + 3.45*np.log10(Te)
            clog = np.full(len(m),loglambda)
            
            nu = np.sum( coll.collisionfreq_PENTA(vt,m,Z,T,n,clog) )
            
            nu_D.append( nu )
            
        # if the input rho was a scalar, then convert the result back to scalar
        # otherwise return the nu_D array
        
        if( np.isscalar(rho) ):
            return nu_D[0]
        else:
            return nu_D
    
    # # TO DO
    # def check_quasi_neutrality(self,rho):
    #     #checks if sum(n_i*Zi)=0        
           
    def print_plasma(self):
        print(f'Current species in plasma are: {", ".join(self.list_of_species)}')


# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)