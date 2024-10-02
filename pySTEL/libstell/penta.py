"""
This library provides a python class for reading and handling PENTA3
results data.
"""

import numpy as np
import sys

# Constants
EC = 1.602176634E-19 # Electron charge [C]

# PENTA Class
class PENTA:
    
    def __init__(self, folder_path, plasma=None, Zions=None):
        #folder_path is a path to the folder containnig the following PENTA3 results files:
        # - flows_vs_Er
        # - flows_vs_roa
        # - fluxes_vs_Er
        # - fluxes_vs_roa
        # - Jprl_vs_roa
        # - contra_vs_roa
        # - sigmas_vs_roa
        
        print('\nPENTA class being created...')
        
        self.folder_path = folder_path
        
        #check if a plasma class is given. If not check if Zions is given
        if plasma is None and Zions is None:
            print('Could not create PENTA class: Need to provide a plasma class or an array with char number (Z) of the ions')
            exit(1)
        elif(plasma is not None and Zions is not None):
            print('Both plasma class and Zions provided. Considering plasma class and discarding Zions array')
            self.list_of_species = plasma.list_of_species
            print(f'List of species: {self.list_of_species}')
            self.Zions = [plasma.Zcharge[species] for species in self.list_of_species]
            #remove electrons
            self.Zions = self.Zions[1:]
            print(f'Zions={self.Zions}')
        elif(plasma is not None and Zions is None):
            self.list_of_species = plasma.list_of_species
            print(f'List of species: {self.list_of_species}')
            self.Zions = [plasma.Zcharge[species] for species in self.list_of_species]
            #remove electrons
            self.Zions = self.Zions[1:]
            print(f'Zions={self.Zions}')
        elif(plasma is None and Zions is not None):
            #check if size of Zions is compatible with number of ions in results files
            self.check_size_Zions(Zions)
            self.list_of_species = ['e'] + [f'i{k}' for k in range(1,len(Zions)+1)]
            print(f'List of species: {self.list_of_species}')
            self.Zions = Zions
            print(f'Zions={self.Zions}')
            
    def check_size_Zions(self,Zions):
        #checks if len(Zions) is the same as the number of ions in the file fluxes_vs_Er
        
        filename = self.folder_path + '/fluxes_vs_Er'
            
        penta = np.loadtxt(filename,skiprows=2)        
        num_ion_species = len(penta[0,:]) - 3
        
        if( len(Zions) is not num_ion_species):
            print(f'ERROR: There are {num_ion_species} species, but only provided Zcharge of {len(Zions)} !!')
            exit(0)




# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)