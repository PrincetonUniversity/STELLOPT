"""
This library provides a python class for multi-species plasmas
"""
import sys

# Constants
EC = 1.602176634E-19 # Electron charge [C]
DA = 1.66053906660E-27 # Dalton
ME = 9.1093837E-31 # Electron mass [kg]

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
        
        self.species_mass = []
        self.species_charge = []
        self.species_Zcharge = []
        
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
                print(f"ERROR: Species {species} does not exist in the database. Please create your plasma manually. Leaving the program.")
                exit(1)  # Exit the program with a status code of 1 (indicating error)
                
    def give_mass_to_species(self,species_list,mass=None):
        
        # Ensure species_list is a list of strings
        if isinstance(species_list, str):
            species_list = [species_list]  # Convert single string to list
        elif not isinstance(species_list, list) or not all(isinstance(item, str) for item in species_list):
            raise ValueError("Species_list must be a string or a list of strings.")
        
        if mass is not None:
            self.species_mass.append( mass )
        else:
            for species in species_list:
                self.species_mass.append( self.mass_database[f'{species}'] )
                
    def give_charge_to_species(self,species_list,charge=None):
        
        # Ensure species_list is a list of strings
        if isinstance(species_list, str):
            species_list = [species_list]  # Convert single string to list
        elif not isinstance(species_list, list) or not all(isinstance(item, str) for item in species_list):
            raise ValueError("Species_list must be a string or a list of strings.")
        
        if charge is not None:
            self.species_charge.append( charge )
        else:
            for species in species_list:
                self.species_charge.append( self.charge_database[f'{species}'] )
    
    def give_Zcharge_to_species(self,species_list,Zcharge=None):
        
        # Ensure species_list is a list of strings
        if isinstance(species_list, str):
            species_list = [species_list]  # Convert single string to list
        elif not isinstance(species_list, list) or not all(isinstance(item, str) for item in species_list):
            raise ValueError("Species_list must be a string or a list of strings.")
        
        if Zcharge is not None:
            self.species_Zcharge.append( Zcharge )
        else:
            for species in species_list:
                self.species_Zcharge.append( self.Zcharge_database[f'{species}'] )        
        
        
    # TO FINISH
    # for now can only add one species at a time, ie, species_name, mass and charge
    # are not lists
    # def add_species(self,species_name,mass,charge,Zcharge):
        
    #     if not isinstance(species_name, str):
    #         print(f"ERROR: species_name must be a string. Leaving the program.")
    #         exit(1) 
    #     else:       
    #         self.list_of_species.append( species_name )
    #         self.give_mass_to_species(species_name,mass)
    #         self.give_charge_to_species(species_name,charge)
    #         self.give_Zcharge_to_species(species_name,Zcharge)
     
    def print_plasma(self):
        print(f'Current species in plasma are: {", ".join(self.list_of_species)}')


# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)