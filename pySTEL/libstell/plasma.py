"""
This library provides a python class for multi-species plasmas
"""
import sys

# Constants
EC = 1.602176634E-19 # Electron charge [C]
DA = 1.66053906660E-27 # Dalton
ME = 9.1093837E-31 # Electron mass [kg]
MP = 1.672621637E-27 # Proton mass [kg]
EPS0 = 8.8541878188E-12 # Vacuum permittivity [F/m]

import numpy as np

class PLASMA:
    
    def __init__(self,list_of_species):
        
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
       
        self.check_species_exist(list_of_species)
        self.list_of_species = list_of_species
        
        self.give_mass_to_species(list_of_species)
        self.give_charge_to_species(list_of_species)
        self.give_Zcharge_to_species(list_of_species)
        
        #create self.ion_species and self.num_ion_species
        self.set_ion_species()
        
        print(f'Plasma created with species: {", ".join(self.list_of_species)}')      
    
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
    
    def get_density_der(self,species,rho):
        # get derivative of density, dn/drho
        # rho can be a number or a list of numbers
        
        #check if species exist in list_of_species
        if species not in self.list_of_species:
            print(f"ERROR: Species {species} is not in the plasma.")
            exit(1)
        
        n0 = self.density[species]['n0']
        nedge = self.density[species]['nedge']
        exponent = self.density[species]['exponent']
        
        rho = np.array(rho)
        
        dens_der = (n0-nedge)*(-exponent*rho**(exponent-1))
        
        return dens_der
    
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
    
    def get_temperature_der(self,species,rho):
        # get derivative of temperature, dT/drho
        # rho can be a number or a list of numbers
        
        #check if species exist in list_of_species
        if species not in self.list_of_species:
            print(f"ERROR: Species {species} is not in the plasma.")
            exit(1)
        
        T0 = self.temperature[species]['T0']
        Tedge = self.temperature[species]['Tedge']
        exponent = self.temperature[species]['exponent']
        
        rho = np.array(rho)
        
        temp_der = (T0-Tedge)*(-exponent*rho**(exponent-1))
        
        return temp_der
    
    def get_thermal_speed(self,species,rho):
        
        temp = self.get_temperature(species,rho)
        vth = np.sqrt(2*EC*temp/self.mass[species])
        
        return vth
    
    def get_averaged_profile(self,species,which_profile,dVdrho):
        # which profile is either 'density' or 'temperature'
        # dVdrho is a function
        
        from scipy.integrate import trapezoid
        
        rho = np.linspace(0,1,100)
        
        if(which_profile=='density'):
            profile = self.get_density(species,rho)
        elif(which_profile=='temperature'):
            profile = self.get_temperature(species,rho)
        else:
            print(f'ERROR: profile is either "density" or "temperature". Cannot be {which_profile}')
            exit(0)
            
        volume = trapezoid(dVdrho(rho),rho)
        
        avg_profile = trapezoid(profile * dVdrho(rho), rho) / volume
        
        return avg_profile  
    
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
    
    def check_quasi_neutrality(self):
        #checks if sum(n_j*Z_j)=0
        
        rho = np.linspace(0,1,100)
        
        nZ = 0.0
        for species in self.list_of_species:
            nZ += self.get_density(species,rho)*self.Zcharge[species]
        
        #normalize
        nZ = nZ / self.get_density('electrons',rho)   
          
        print(f'\nsum(n_jZ_j)/ne (rho) = {nZ}')       
    
    def set_ion_species(self):
        
        # ion species in the order that appears in list_of_species
        self.ion_species = [species for species in self.list_of_species if species != 'electrons']
        
        print(f'Ion species: {self.ion_species}')
        
        self.num_ion_species = len(self.ion_species)
        
    
    def write_plasma_profiles_to_PENTA1(self,rho,filename=None):
        # first line: number of rhos
        # from second line: 
        # 1st column: rho=r/a
        # 2nd column: ne (in units of 10^18 m^-3)
        # 3rd column: Te (eV)
        # 4th and 5th: Ti and ni of ion1
        # 6th and 7th: Ti and ni of ion 2
        # etc
        # NOTICE it's not a mistake: for the electrons the order is n and T, while for the ions is T and n
        
        if filename is None:
            filename = 'plasma_profiles.dat'
        
        fact = 1e18
        
        ne = self.get_density('electrons',rho) / fact
        Te = self.get_temperature('electrons',rho)
        
        #Precompute ion densities and temperatures
        ion_densities = {ion_s: self.get_density(ion_s, rho)/fact for ion_s in self.ion_species}
        ion_temperatures = {ion_s: self.get_temperature(ion_s, rho) for ion_s in self.ion_species}
        
        with open(filename, 'w') as file:
            # Write the size of rho as the first line
            file.write(f'{rho.size}\n')
            
            # Write the data to the file
            for i in range(rho.size):
                # First column: electron density, second column: electron temperature
                row_data = [rho[i], ne[i], Te[i]]
                
                # Append ion densities and temperatures
                for ion_s in self.ion_species:
                    row_data.extend([ion_temperatures[ion_s][i], ion_densities[ion_s][i]])
                
                # Write the row to the file, formatted as space-separated values
                file.write(" ".join(map(str, row_data)) + '\n')
    
    def write_plasma_profiles_to_PENTA3(self,rho,filename=None):
        # first line: number of rhos
        # from second line: 
        # 1st column: rho=r/a
        # 2nd column: ne (in units of 10^18 m^-3)
        # 3rd column: Te (eV)
        # 4th and 5th: ni and Ti of ion1
        # 6th and 7th: ni and Ti of ion 2
        # etc
        
        if filename is None:
            filename = 'plasma_profiles.dat'
        
        fact = 1e18
        
        ne = self.get_density('electrons',rho) / fact
        Te = self.get_temperature('electrons',rho)
        
        #Precompute ion densities and temperatures
        ion_densities = {ion_s: self.get_density(ion_s, rho)/fact for ion_s in self.ion_species}
        ion_temperatures = {ion_s: self.get_temperature(ion_s, rho) for ion_s in self.ion_species}
        
        with open(filename, 'w') as file:
            # Write the size of rho as the first line
            file.write(f'{rho.size}\n')
            
            # Write the data to the file
            for i in range(rho.size):
                # First column: electron density, second column: electron temperature
                row_data = [rho[i], ne[i], Te[i]]
                
                # Append ion densities and temperatures
                for ion_s in self.ion_species:
                    row_data.extend([ion_densities[ion_s][i],ion_temperatures[ion_s][i]])
                
                # Write the row to the file, formatted as space-separated values
                file.write(" ".join(map(str, row_data)) + '\n')
                
                           
    def write_PENTA_namelist(self,filename=None):
        
        if filename is None:
            filename = 'ion_params'

        #mass of ion to proton mass
        miomp = {}
        for ion in self.ion_species:
            miomp[ion] = self.mass[ion] / MP
            
        with open(filename, 'w') as file:
            # write NAMELIST header
            file.write('&ion_params\n')
            # write number of species
            file.write(f'num_ion_species={self.num_ion_species}\n')
            # write Zcharge
            Zcharge = [self.Zcharge[ion] for ion in self.ion_species]
            file.write('Z_ion_init=' + ','.join(f'{z:.6f}d0' for z in Zcharge)+ ',\n')
            # write mass
            mass = [miomp[ion] for ion in self.ion_species]
            file.write('miomp_init=' + ','.join(f'{m:.6f}d0' for m in mass) + ',\n')
            # close NAMELIST
            file.write('&end\n')         
           
    def print_plasma(self):
        print(f'Current species in plasma are: {", ".join(self.list_of_species)}')
        
    def write_THRIFT_profiles(self,time_array,filename=None):
        # generates a profiles file for THRIFT
        # assumes that profiles are constant in time and writes the electron 
        # and ions density (m^-3) and temperature (eV) profiles
        # saved in the plasma class at all times given in time_array
        
        # parts of this routine copied from STELLOPT/BENCHMARKS/THRIFT_TEST/make_profiles.py
        
        import h5py
        import sys
        
        if filename is None:
            filename = 'profiles.h5'
            
        rho = np.linspace(0,1,100)
            
        nrho = len(rho)
        nt   = len(time_array)
        
        if(nt < 4):
            print('PLEASE GIVE time_array with more times (>3)')
            exit(0)
            
        nZ   = self.num_ion_species
        
        Zcharge_ions = np.array( [self.Zcharge[ion] for ion in self.ion_species], dtype=float )
        mass_ions    = [self.mass[ion] for ion in self.ion_species]
        
        # print(f'Zcharge_ions={Zcharge_ions}')
        # print(f'mass_ions={mass_ions}')
        
        # ne
        ne = np.array( self.get_density('electrons',rho) )
        ne = np.tile(ne, (nt,1)).T
        
        # Te
        Te = np.array( self.get_temperature('electrons',rho) )
        Te = np.tile(Te, (nt,1)).T
        
        # ni
        ni = np.array( [self.get_density(ion,rho) for ion in self.ion_species] )
        ni = np.repeat(ni[:,:,np.newaxis], nt, axis=2)
        ni = np.transpose(ni, axes=[1,2,0])
        
        # Ti
        Ti = np.array( [self.get_temperature(ion,rho) for ion in self.ion_species] )
        Ti = np.repeat(Ti[:,:,np.newaxis], nt, axis=2)
        Ti = np.transpose(Ti, axes=[1,2,0])

        hf = h5py.File(filename, 'w')
        
        hf.create_dataset('nrho', data=nrho)
        hf.create_dataset('nt', data=nt)
        hf.create_dataset('nion', data=nZ)
        hf.create_dataset('raxis_prof', data=rho)
        hf.create_dataset('taxis_prof', data=time_array)
        hf.create_dataset('Z_prof', data=Zcharge_ions)
        hf.create_dataset('mass_prof', data=mass_ions)
        hf.create_dataset('ne_prof', data=ne)
        hf.create_dataset('te_prof', data=Te)
        hf.create_dataset('ni_prof', data=ni)
        hf.create_dataset('ti_prof', data=Ti)

        hf.close()
        
        print(f'{filename} created with success!')
        
    def get_pressure_polynomial_coefficients(self,deg_fit=10):
        # this computes the AM coefficients and the PRES_SCALE scalar for a VMEC input
        # assuming that PMASS_TYPE = 'power_series'
        # p = sum(n_j*T_j), j=all species, and pressure comes in Pa
        
        # note that the AM coefficients correspond to the polynomial in s=roa^2
        # so we will make a polyfit with deg_fit
        
        import matplotlib.pyplot as plt
        from scipy.interpolate import CubicSpline
        
        roa = np.linspace(0,1,100)
        s = roa**2
        
        p=0
        for species in self.list_of_species:
            
            n = self.get_density(species,roa)
            T = self.get_temperature(species,roa) 
            
            # total pressure polynomial in Pascal units
            p += n*T*EC
        
        # make a polynomial fit. 
        # in the future it might be better to use splines which VMEC accepts?    
        p_pol = np.polyfit(s,p,deg_fit)
        p_pol = np.poly1d(p_pol)
            
        # normalize p_pol by p(s=0) and use this value as PRES_SCALE
        PRES_SCALE = p_pol(0.0)
        p_pol = p_pol / PRES_SCALE
        
        AM = p_pol.c[-1:0:-1]
        
        plt.plot(s,p,'.-')
        plt.plot(s,PRES_SCALE*p_pol(s),'-')
        plt.xlabel('s')
        plt.ylabel('pressure [Pa]')
        plt.show
        
        print(f'PRES_SCALE = {PRES_SCALE:.7e}')
        print(f'AM = {AM}')
        
        # make a cubic spline
        cs = CubicSpline(s,p)
        s_VMEC = np.linspace(0,1,50)
        
        pres = cs(s_VMEC)   
        pres = pres / PRES_SCALE             
        
        
        print("PMASS_TYPE = 'cubic_spline' ")
        print(f'PRES_SCALE = {PRES_SCALE:.7e}')
        print(f'AM_AUX_S = {s_VMEC}')
        print(f'AM_AUX_F = {pres}')
        
        return AM,PRES_SCALE

        
    def print_SFINCS_namelist(self,roa):
        # prints in the command line physics and species parameters namelist of for SFINCS
        
        n_species = np.array( [self.get_density(species,rho=roa) for species in self.list_of_species] )
        T_species = np.array( [self.get_temperature(species,rho=roa) for species in self.list_of_species] )
        m_species = np.array( [self.mass[species] for species in self.list_of_species] )
        Z_species = np.array( [self.Zcharge[species] for species in self.list_of_species] )
        
        nder_species = np.array( [self.get_density_der(species,rho=roa) for species in self.list_of_species] )
        Tder_species = np.array( [self.get_temperature_der(species,rho=roa) for species in self.list_of_species] )

        ## reference values
        nBar = np.max(n_species) # pick largest value. must be in m^-3
        mBar = MP # must be in kg 
        TBar = np.max(T_species) # must be in eV
        
        vBar = np.sqrt(2*TBar*EC/mBar)
        
        print(f'nBar = {nBar}')
        print(f'vBar = {vBar}')
        
        # mandatory values -- these values (BBar=1 and RBar=1) are mandatory when mag field is read from VMEC wout file (see SFINCS documentation)
        BBar = 1.0
        RBar = 1.0
        
        # print(f'jbs.B = {-1.2e-4*EC*nBar*vBar*BBar}')
        
        # compute loglambda as in PENTA
        Te = self.get_temperature('electrons',rho=roa)
        ne = self.get_density('electrons',rho=roa)
        if(Te>50):
            loglambda = 25.3 - 1.15*np.log10(ne/1e6) + 2.3*np.log10(Te)
        else:
            loglambda = 23.4 - 1.15*np.log10(ne/1e6) + 3.45*np.log10(Te)
        
        nuHat = 4*np.sqrt(2*np.pi)*nBar*EC**4*loglambda / ( 3*(4*np.pi*EPS0)**2 * np.sqrt(mBar) * (EC*TBar)**1.5 )
        
        # compute physics parameters
        Delta = mBar*vBar / (EC*BBar*RBar)
        alpha = 1.0
        nu_n = nuHat * RBar/vBar
        
        mHats = m_species/mBar
        nHats = n_species/nBar
        THats = T_species/TBar
        dNHatdrNs = nder_species/nBar
        dTHatdrNs = Tder_species/TBar
        
        # print output
        
        print('&speciesParameters')
        print(f"Zs = {' '.join(map(str, Z_species))}")
        print(f"mHats = {' '.join(map(str, mHats))}")
        print(f"nHats = {' '.join(map(str, nHats))}")
        print(f"THats = {' '.join(map(str, THats))}")
        print(f"dNHatdrNs = {' '.join(map(str, dNHatdrNs))}")
        print(f"dTHatdrNs = {' '.join(map(str, dTHatdrNs))}")
        print('/')

        print('&physicsParameters')
        print(f'Delta = {Delta}')
        print(f'alpha = {alpha}')
        print(f'nu_n = {nu_n}')

        print('dont forget the rest of the parameters...')
        
        
    def print_SFINCS_list_namelist(self,roa_list,folder_path):
        # saves input.namlist inside folder_path/surface_k
        
        import os
        
        ## reference values
        nBar = 1e20  # must be in m^-3
        mBar = MP    # must be in kg 
        TBar = 1e3   # must be in eV
        
        vBar = np.sqrt(2*TBar*EC/mBar)
        
        # mandatory values -- these values (BBar=1 and RBar=1) are mandatory when mag field is read from VMEC wout file (see SFINCS documentation)
        BBar = 1.0
        RBar = 1.0
        
        for k,roa in enumerate(roa_list):
            
            folder_name = f'surface_{k+1}'  # Folders will be surface_1, surface_2, etc.
            os.makedirs(folder_path+'/'+folder_name, exist_ok=True)
        
            n_species = np.array( [self.get_density(species,rho=roa) for species in self.list_of_species] )
            T_species = np.array( [self.get_temperature(species,rho=roa) for species in self.list_of_species] )
            m_species = np.array( [self.mass[species] for species in self.list_of_species] )
            Z_species = np.array( [self.Zcharge[species] for species in self.list_of_species] )
            
            nder_species = np.array( [self.get_density_der(species,rho=roa) for species in self.list_of_species] )
            Tder_species = np.array( [self.get_temperature_der(species,rho=roa) for species in self.list_of_species] )

            # compute loglambda as in PENTA
            Te = self.get_temperature('electrons',rho=roa)
            ne = self.get_density('electrons',rho=roa)
            if(Te>50):
                loglambda = 25.3 - 1.15*np.log10(ne/1e6) + 2.3*np.log10(Te)
            else:
                loglambda = 23.4 - 1.15*np.log10(ne/1e6) + 3.45*np.log10(Te)
            
            nuHat = 4*np.sqrt(2*np.pi)*nBar*EC**4*loglambda / ( 3*(4*np.pi*EPS0)**2 * np.sqrt(mBar) * (EC*TBar)**1.5 )
            
            # compute physics parameters
            Delta = mBar*vBar / (EC*BBar*RBar)
            alpha = 1.0
            nu_n = nuHat * RBar/vBar
            
            mHats = m_species/mBar
            nHats = n_species/nBar
            THats = T_species/TBar
            dNHatdrNs = nder_species/nBar
            dTHatdrNs = Tder_species/TBar
            
            file_content = f"""! Input file for SFINCS version 3.

&general
RHSmode = 1
ambipolarSolve = .true. !.false.
!Er_min = -10
!Er_max = 10
!NEr_ambipolarSolve = 10
/

&geometryParameters
geometryScheme = 5

inputRadialCoordinate = 3
rN_wish = {roa}
inputRadialCoordinateForGradients = 3   !the radial coordinate of the gradients given in species parameters is rN=sqrt(PHI/PHIEDGE)

VMECRadialOption = 1  !get the nearest available flux surface from VMEC HALF grid 
equilibriumFile = "wout_beta_2.nc"
min_Bmn_to_load = 1e-4
/

&speciesParameters
Zs = {' '.join(map(str, Z_species))}
mHats = {' '.join(map(str, mHats))}
nHats = {' '.join(map(str, nHats))}
THats = {' '.join(map(str, THats))}
dNHatdrNs = {' '.join(map(str, dNHatdrNs))}
dTHatdrNs = {' '.join(map(str, dTHatdrNs))}
/

&physicsParameters
Delta = {Delta}
alpha = {alpha}
nu_n = {nu_n}

dPhiHatdrN = 1.0 ! since we are looking for the ambipolar value, I guess this one is free??

collisionOperator = 0

includeXDotTerm = .true.
includeElectricFieldTermInXiDot = .true.
useDKESExBDrift = .false.

includePhi1 = .false.  ! well, at some point this should be True to be more complete...

magneticDriftScheme = 1  ! this includes poloidal and toroidal magnetic drifts
/

&resolutionParameters
Ntheta = 23 ! needs to be an odd number
Nzeta = 91 ! needs to be an odd number and at low collisionality might be needeed to be of the order 100 to converge

Nxi = 70
Nx = 6
solverTolerance = 1d-6
/

&otherNumericalParameters
/

&preconditionerOptions
/

&export_f
  export_full_f = .false.
  export_delta_f = .false.
/
"""
            
             # Define the file path
            file_path = os.path.join(folder_path+'/'+folder_name, 'input.namelist')
            
            # Write the content to the file
            with open(file_path, 'w') as f:
                f.write(file_content)

            print(f"Created {file_path}")
        

# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)