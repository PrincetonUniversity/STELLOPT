"""
This library provides a python class for solving density 
and pressure transport equations
"""

import numpy as np
import sys
import matplotlib.pyplot as plt
sys.path.insert(1,'/home/antonio/STELLOPT/pySTEL/libstell')

from plasma import PLASMA

# Constants
EC = 1.602176634E-19 # Electron charge [C]
EPS0 = 8.8541878188E-12 # Vacuum permittivity [F/m]

class PLASMA_SOLVER:
    
    def __init__(self, list_of_species, tau_fast_alphas=None, tau_thermal_alphas=None):
        
        from collections import defaultdict
        
        # later this can be changed
        valid_species = {'electrons', 'deuterium', 'tritium'}
        invalid_species = set(list_of_species) - valid_species
        if invalid_species:
            raise ValueError(f"Invalid species found: {invalid_species}")
        else:
            self.list_of_species = list_of_species
            
        # create plasma class with list_of_species
        self.plasma = PLASMA(list_of_species)

        # in case taus are provided, will solve alphas density using a simple model
        # this influences ne (and alphas Bremsstrahlung)
        if((tau_fast_alphas is None) or (tau_thermal_alphas is None)):
            self.solving_alphas_density = False
        else:
            self.solving_alphas_density = True
            self.tau_fast_alphas = tau_fast_alphas
            self.tau_thermal_alphas = tau_thermal_alphas
            
        # initialize dictionaries
        self.edge_density_BC = {}
        self.edge_pressure_BC = {}
        
        self.initial_density = {}
        self.initial_pressure = {}
        
        self.energy_sources = defaultdict(lambda: defaultdict(dict))
        self.density_sources = defaultdict(lambda: defaultdict(dict))
        
        self.heat_fluxes_info = defaultdict(lambda: defaultdict(dict))
        self.particle_fluxes_info = defaultdict(lambda: defaultdict(dict))
        
        
        print(f'Solvers for pressure of {self.list_of_species} INITIALIZED!')
        print(f'SOLVING FOR ALPHAS DENSITY: {self.solving_alphas_density}')
        
            

    
    def set_edge_boundary_condition(self,field: str,species: str,val: float):
        # set edge boundary Dirichlet boundary condition
        # field can be 'density', 'temperature' or 'pressure'
        
        # check species exist in list_of_species
        if species not in self.list_of_species:
            raise ValueError(f"ERROR: Species {species} is not in the plasma.")
        
        if field=='density':
            self.edge_density_BC[species] = val
        elif field=='pressure':
            self.edge_pressure_BC[species] = val
        elif field=='temperature':
            if (species not in self.edge_density_BC):
                raise ValueError('Need to set density BC before setting temperature BC')
            else:
                self.edge_pressure_BC[species] = self.edge_density_BC[species] * EC * val
        else: 
            raise ValueError('field not valid...')
                
        
    def set_initial_profile(self,field: str, species: str, rho_vals, profile_vals):
        
        from scipy.interpolate import CubicSpline
        
        # check species exist in list_of_species
        if species not in self.list_of_species:
            print(f"ERROR: Species {species} is not in the plasma.")
            exit(1)
        
        # check rho_vals are in the range [0,1]
        if( np.any((rho_vals<0) | (rho_vals>1))):
            print('ERROR: rho_vals must be in the domain [0,1]')
            exit(0)
                            
        # create interpolating function
        if field=='density':
            self.initial_density[species] = CubicSpline(rho_vals,profile_vals)
        elif field=='pressure':
            self.initial_pressure[species] = CubicSpline(rho_vals,profile_vals)
        elif field=='temperature':
            if (species not in self.initial_density):
                raise ValueError('Need to set density initial profile before setting temeprature initial profile')
            else:
                pressure = self.initial_density[species](rho_vals) * EC * profile_vals
                self.initial_pressure[species] = CubicSpline(rho_vals,pressure)
        else:
            raise ValueError('field is not valid...')
        
    def set_equilibrium(self,type: str,wout_path=None,aminor=None,Rmajor=None):
        
        from libstell.vmec import VMEC
        from scipy.interpolate import CubicSpline
        
        match type:
            case 'VMEC':
               if(wout_path is None):
                   print('ERROR: For a VMEC equilibrium, wout_path must be given!')
                   exit(0)
               else:
                    self.wout_path = wout_path
                    # get dVdr from file; use VMEC class
                    vmec_out = VMEC()
                    vmec_out.read_wout(wout_path)
                    
                    # 4pi^2*dVds -- vmec.py already does h2f, so vprime is in full grid
                    vp = vmec_out.vp[:].flatten()
                    
                    self.aminor = vmec_out.aminor
                    
                    roa = np.sqrt(vmec_out.phi / vmec_out.phi[-1])
                    roa = roa.flatten()
                    
                    #dVdr analytic = dVds * 2\rho / a
                    dVdr_analytic = (2*np.pi)**2 * vp * 2.*roa / self.aminor
                    
                    self.dVdr = CubicSpline(roa,dVdr_analytic)
                    
                    self.B = np.sqrt(np.squeeze(vmec_out.bdotb)[0])   
            case 'cylindrical':
                if(aminor is None or Rmajor is None):
                    print('ERROR: For a cylindrical equilibrium, Rmajor and aminor must be given')
                    exit(0)
                else:
                    self.aminor=aminor
                    self.dVdr = lambda rho: 4*np.pi*np.pi*Rmajor*aminor  * rho
                    
    def set_energy_source(self,species,source_type, total_power=None, sigma_rho=None, fraction_alpha_heating=None, cte_source=None, time_dependent_factor=None):
        # electrons: 'Bremsstrahlung', 'Coll_Heat_Exchange', 'Er', 'external', 'alpha_heating'
        # ions: 'Coll_Heat_Exchange', 'Er', 'external', 'alpha_heating'
        # 'constant' is for benchmarking
        
        from scipy.interpolate import CubicSpline
        import inspect
        
        # check species exist in list_of_species
        if species not in self.list_of_species:
            print(f"ERROR: Species {species} is not in the plasma.")
            exit(1)
         
        match source_type:
            case 'Bremsstrahlung':
                if(species != 'electrons'): 
                    print('ERROR: Bremsstrahlung is only source for electrons')
                    exit(0)
                else:
                    self.energy_sources[species][source_type] = {}
            #
            case 'external_gaussian':
                if((total_power is None) or (sigma_rho is None)):
                    print('ERROR: Need to provide total_power [W] and sigma_rho for gaussian external source')
                    exit(1) 
                else:
                    self.energy_sources[species][source_type] = {'total_power' : total_power, 'sigma_rho' : sigma_rho}
            #
            case 'time_dependent_gaussian':
                if((total_power is None) or (sigma_rho is None) or (time_dependent_factor is None)):
                    print('ERROR: Need to provide total_power [W], sigma_rho and a time depenedent factof for time-dependent gaussian')
                    exit(1) 
                else:
                    self.energy_sources[species][source_type] = {'total_power' : total_power, 'sigma_rho' : sigma_rho, 'time_factor': time_dependent_factor }
            #
            case 'Coll_Heat_Exchange':
                self.energy_sources[species][source_type] = {}
            #
            case 'Er':
                self.energy_sources[species][source_type] = {}
            #
            case 'alpha_heating':
                if(fraction_alpha_heating is None):
                    print('ERROR: fraction_alpha_heating is needed. Usually ~80% electrons, 20% ions')
                    exit(0)
                else:
                    self.energy_sources[species][source_type] = {'fraction_alpha_heating': fraction_alpha_heating}
            #
            case 'constant':
                if(cte_source is None):
                    print('ERROR: cte_source is needed in order to generate a constant source.')
                    exit(0)
                else:
                    self.energy_sources[species][source_type] = {'cte_source' : cte_source}
            #
            case _:
                print(f'ERROR: Source type {source_type} is NOT possible')
                exit(0)
                
    def set_density_source(self,species,source_type, Smax=None, rho_0=None, sigma_rho=None, cte_source=None, time_dependent_factor=None):
        
        # check species exist in list_of_species
        if species not in self.list_of_species:
            raise ValueError(f"ERROR: Species {species} is not in the plasma.")
         
        match source_type:
            case 'external_gaussian':
                if((Smax is None) or (sigma_rho is None) or (rho_0 is None)):
                    print('ERROR: Need to provide Smax [part/m^3], rho_0 and sigma_rho for gaussian external source')
                    exit(1) 
                else:
                    self.density_sources[species][source_type] = {'Smax' : Smax, 'rho_0' : rho_0, 'sigma_rho' : sigma_rho}
            #
            case 'time_dependent_gaussian':
                if((Smax is None) or (rho_0 is None) or (sigma_rho is None) or (time_dependent_factor is None)):
                    print('ERROR: Need to provide Smax [par/m^3], rho_0, sigma_rho and a time depenedent factof for time-dependent gaussian')
                    exit(1) 
                else:
                    self.density_sources[species][source_type] = {'Smax' : Smax, 'rho_0' : rho_0, 'sigma_rho' : sigma_rho, 'time_factor': time_dependent_factor }
            #
            case 'alpha_generation':
                self.density_sources[species][source_type] = {}
            #
            case 'constant':
                if(cte_source is None):
                    print('ERROR: cte_source is needed in order to generate a constant source.')
                    exit(0)
                else:
                    self.density_sources[species][source_type] = {'cte_source' : cte_source}
            #
            case _:
                print(f'ERROR: Source type {source_type} is NOT possible')
                exit(0)
                
    def set_heat_fluxes(self, type: str,dkes_folder=None,surfaces=None,theta=None,chi=None,chi_base=None,aLT_critical=None,alpha=None,stiffness=None,chi_electrons=None):
        # sets type of fluxes
        # OPTION1: type='dkespenta'; dkes_folder and surfaces(list of integers) must be provided
        # OPTION2: type='diffusive'; chi must be provided (heat diffusivity; assumes same val for all species)
            
        match type:
            case 'dkespenta':
                #checks that dkes_folder and surfaces are provided
                if((dkes_folder is None) or (surfaces is None) or (theta is None)):
                    print('ERROR: dkes_folder, surfaces and theta must be provided!!')
                    exit(0)
                    
                self.heat_fluxes_info['type'] = type
                self.heat_fluxes_info[type]['dkes_folder'] = dkes_folder
                self.heat_fluxes_info[type]['surfaces'] = surfaces
                self.heat_fluxes_info[type]['theta'] = theta
                
            case 'diffusive':
                #checks that diffusion coefficients are provided
                if(chi is None):
                    print('ERROR: chi must be provided!')
                    exit(0)
                self.heat_fluxes_info['type'] = type
                self.heat_fluxes_info[type]['chi'] = chi
                
            case 'beurskens':
                if( (chi_base is None) or (aLT_critical is None) or (alpha is None) or (stiffness is None) or (chi_electrons is None) ):
                    print('ERROR: chi_base, aLT_critical, alpha, stiffness and chi_electrons must be given!')
                    exit(0)
                self.heat_fluxes_info['type'] = type
                self.heat_fluxes_info[type]['chi_base'] = chi_base
                self.heat_fluxes_info[type]['aLT_critical'] = aLT_critical
                self.heat_fluxes_info[type]['alpha'] = alpha
                self.heat_fluxes_info[type]['stiffness'] = stiffness
                self.heat_fluxes_info[type]['chi_electrons'] = chi_electrons
                
    def set_particle_fluxes(self, type: str,dkes_folder=None,surfaces=None,theta=None,Dn=None):
        # sets type of fluxes
        # OPTION1: type='dkespenta'; dkes_folder and surfaces(list of integers) must be provided
        # OPTION2: type='diffusive'; Dn and chi must be provided (partical and heat collisional diffusion coefficients)
            
        match type:
            case 'dkespenta':
                #checks that dkes_folder and surfaces are provided
                if((dkes_folder is None) or (surfaces is None) or (theta is None)):
                    print('ERROR: dkes_folder, surfaces and theta must be provided!!')
                    exit(0)
                    
                self.particle_fluxes_info['type'] = type
                self.particle_fluxes_info[type]['dkes_folder'] = dkes_folder
                self.particle_fluxes_info[type]['surfaces'] = surfaces
                self.particle_fluxes_info[type]['theta'] = theta
                
            case 'diffusive':
                #checks that diffusion coefficient is provided
                if(Dn is None):
                    print('ERROR: Dn must be provided!')
                    exit(0)
                self.particle_fluxes_info['type'] = type
                self.particle_fluxes_info[type]['chi'] = Dn
        
    # def check_boundary_conditions(self):
    #     #Checks ALL BCs are set and that these are consistent with initial profiles
        
    #     ## THIS FROUTINE SHOULD BE CALLED INSIDE 'RUN'