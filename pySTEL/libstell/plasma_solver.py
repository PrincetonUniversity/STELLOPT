"""
This library provides a python class for solving density 
and pressure transport equations
"""

import numpy as np
import sys
from time import perf_counter
from numba import njit
from concurrent.futures import ProcessPoolExecutor

from libstell.plasma import PLASMA
from libstell.penta import PENTA
from libstell.libpenta import LIBPENTA, _init_NEO_worker, _call_PENTA_surface_worker

# Constants
EC = 1.602176634E-19 # Electron charge [C]
EPS0 = 8.8541878188E-12 # Vacuum permittivity [F/m]

class PLASMA_SOLVER:
    """Class for solving transport equations (density and pressure)

	"""
    
    def __init__(self, list_of_species, solve_fast_alphas=False, tau_fast_alphas=0.5, constrain_nT=False, add_NEO=False):
        
        from collections import defaultdict
        
        self.list_of_species = list_of_species
            
        # create plasma class with list_of_species
        self.plasma = PLASMA(list_of_species)
        
        self.solve_fast_alphas = solve_fast_alphas
        # in case solve_fast_alphas is True but deuterium&tritium are not in the plasma
        # an error is given
        if(solve_fast_alphas and ('deuterium' not in list_of_species or 'tritium' not in list_of_species)):
            raise ValueError(f'ERROR: solve_fast_alphas was set to True, but deuterium and/or tritium not in the plasma!')
        elif(solve_fast_alphas):
            self.tau_fast_alphas = tau_fast_alphas
            
        if(constrain_nT):
            # if True, nT is assumed to be equal to nD
            self.constrain_nT = True
        else:
            self.constrain_nT = False
        
        self.add_NEO = add_NEO
        if(add_NEO):
            # Create a libpenta class. Is used when computing fluxes
            self.libPenta = LIBPENTA()
            print('Solving for NEOCLASSICAL fluxes')
            
        # initialize dictionaries
        self.edge_density_BC = {}
        self.edge_pressure_BC = {}
        
        self.initial_density = {}
        self.initial_pressure = {}
        
        self.energy_sources = defaultdict(lambda: defaultdict(dict))
        self.particle_sources = defaultdict(lambda: defaultdict(dict))
        
        self.heat_fluxes_info = defaultdict(lambda: defaultdict(dict))
        self.particle_fluxes_info = defaultdict(lambda: defaultdict(dict))

        print(f'Transport plasma solver for {self.list_of_species} INITIALIZED!')
        if(self.solve_fast_alphas): print(f'SOLVING FOR FAST ALPHAS!')
    
    def set_edge_boundary_condition(self,field: str,species: str,val: float):
        """
        Sets edge Dirichlet boundary conditions.
        Field can be 'density', 'temperature' or 'pressure'
        """
        
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
        """
        Sets profiles at t=t_start
        Field can be 'density', 'temperature' or 'pressure'
        """
        
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
                raise ValueError('Need to set density initial profile before setting temperature initial profile')
            else:
                pressure = self.initial_density[species](rho_vals) * EC * profile_vals
                self.initial_pressure[species] = CubicSpline(rho_vals,pressure)
        else:
            raise ValueError('field is not valid...')
        
    def read_restart_file(self,restart_filepath,tstart):
        """
        Reads a restart .joblib file and sets initial profiles & BC's 
        according to last iteration in file.
        If tstart of current run does not match tend of restart file
        an error is given
        """
        import joblib
        from pathlib import Path
        
        # Check if extension of restart_filepath is .joblib; if not, add
        restart_filepath = str(Path(restart_filepath).with_suffix(".joblib"))
        
        restart_solver = joblib.load(restart_filepath)
        
        # Check tstart of current simulation IS EQUAL to tend in restart file
        if( not np.isclose(tstart,restart_solver.time[-1]) ):
            raise ValueError('tstart of current simulation and tend of restart file do not match!')
        
        # Check list_of_species in file are the same as those in this run
        if (restart_solver.list_of_species != self.list_of_species):
            raise ValueError('Species in restart file different from species in current solver!')
        
        # Set initial conditions and BC's of all species
        for species in self.list_of_species:
            
            self.set_initial_profile('density', species, rho_vals=restart_solver.rho_grid, profile_vals=restart_solver.N[species][-1,:])
            self.set_initial_profile('temperature', species, rho_vals=restart_solver.rho_grid, profile_vals=restart_solver.T[species][-1,:])
            
            self.set_edge_boundary_condition('density', species, restart_solver.N[species][-1,-1])
            self.set_edge_boundary_condition('temperature', species, restart_solver.T[species][-1,-1]) 
            
        # Get alphas density 
        if(self.solve_fast_alphas):
            # check restart_solver has alphas
            if 'alphas_fast' not in restart_solver.N:
                raise ValueError('restart file does not have fast alphas density! Yet you want to solve with alphas...')
            else:
                self.alphas_fast_density_restart = restart_solver.N['alphas_fast'][-1,:]
              
    def set_equilibrium(self,type: str,wout_path=None,aminor=None,Rmajor=None,B=None):
        """
        Sets magnetic equilibrium
        Type can be 'VMEC' (need to provide path to wout file) 
        or 'cylindrical (need to provide aminor, Rmajor and B)
        If 'VMEC' then Bsq(r) as defined in wout is used; if 'cylindrical' Bsq(r)=B=te
        """
        
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
                    self.Rmajor = vmec_out.rmajor
                    
                    roa = np.sqrt(vmec_out.phi / vmec_out.phi[-1])
                    roa = roa.flatten()
                    
                    #dVdr analytic = dVds * 2\rho / a
                    dVdr_analytic = (2*np.pi)**2 * vp * 2.*roa / self.aminor
                    
                    self.dVdr = CubicSpline(roa,dVdr_analytic)
                    
                    self.Baxis = np.sqrt(np.squeeze(vmec_out.bdotb)[0])   
                    self.iota23 = CubicSpline(roa,np.squeeze(vmec_out.iotaf))(2.0/3.0)
                    
                    self.Bsq_spline = CubicSpline(roa,np.squeeze(vmec_out.bdotb))
                    
                    # stella reference magnetic field
                    self.Bref = vmec_out.phi[-1]/ (np.pi*self.aminor**2)
                    
            case 'cylindrical':
                if(aminor is None or Rmajor is None or B is None):
                    print('ERROR: For a cylindrical equilibrium, Rmajor, aminor and B must be given')
                    exit(0)
                else:
                    self.aminor = aminor
                    self.Rmajor = Rmajor
                    dVdr = lambda rho: 4*np.pi*np.pi*Rmajor*aminor  * rho
                    rho = np.linspace(0,1,100)
                    self.dVdr = CubicSpline(rho,dVdr(rho))
                    self.Baxis = B
                    self.Bref = B
                    
    def set_energy_source(self,species,source_type, total_power=None, sigma_rho=None, rho_0=None, 
        fraction_alpha_heating=None, cte_source=None, time_dependent_factor=None, lambda_function_2D=None,
        max_total_power=None, time_dependent_fusion_power=None, time_dependent_DT_temp_axis=None, time_dependent_electron_temp_axis=None,
        pidK=1.0,pidI=1.0E10,pidD=0.0,noise_level=0.00):
        """
        Sets energy sources for a given species. The source_type can be:
        'external_gaussian', 'time_dependent_gaussian', 'Coll_Heat_Exchange', 'Er', 'alpha_heating', 'constant', 'lambda_2D' and 'Bremsstrahlung' (only for electrons)
        The 'Coll_Heat_Exchange' only needs to be set for electrons; once it's set it will be computed for ALL species
        """
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
                if((total_power is None) or (sigma_rho is None) or (rho_0 is None)):
                    print('ERROR: Need to provide total_power [W], sigma_rho and rho_0 for gaussian external source')
                    exit(1) 
                else:
                    self.energy_sources[species][source_type] = {'total_power' : total_power, 'sigma_rho' : sigma_rho, 'rho_0' : rho_0}
            #
            case 'time_dependent_gaussian':
                if((total_power is None) or (sigma_rho is None) or (time_dependent_factor is None) or (rho_0 is None)):
                    print('ERROR: Need to provide total_power [W], sigma_rho, rho_0 and a time depenedent factor for time-dependent gaussian')
                    exit(1) 
                else:
                    self.energy_sources[species][source_type] = {'total_power' : total_power, 'sigma_rho' : sigma_rho, 'rho_0' : rho_0, 'time_factor': time_dependent_factor }
            #
            case 'PID_pfuse_gaussian':
                if((rho_0 is None) or (sigma_rho is None) or (max_total_power is None) or (time_dependent_fusion_power is None)):
                    raise ValueError('ERROR: Need to provide time_dependent_fusion_power, max_total_power, rho_0 and sigma_rho for PID electron density gaussian')
                else:
                    self.energy_sources[species][source_type] = {'max_total_power' : max_total_power, 'rho_0' : rho_0, 
                    'sigma_rho' : sigma_rho, 'time_dependent_fusion_power': time_dependent_fusion_power,
                    'pid_K' : pidK, 'pid_Ti' : pidI, 'pid_Td' : pidD, 'pid_I' : 0.0,
                    'noise_level' : noise_level, 'previous_error' : 0.0}
                    print(f'Using PID_pfuse_gaussian with pid_K={pidK}m^3/s, pid_tauI={pidI}s, pid_tauD={pidD}s and noise_level={noise_level}')
            #
            case 'PID_itemp_gaussian':
                if((rho_0 is None) or (sigma_rho is None) or (max_total_power is None) or (time_dependent_DT_temp_axis is None)):
                    raise ValueError('ERROR: Need to provide time_dependent_DT_temp_axis, max_total_power, rho_0 and sigma_rho for PID electron density gaussian')
                else:
                    self.energy_sources[species][source_type] = {'max_total_power' : max_total_power, 'rho_0' : rho_0, 
                    'sigma_rho' : sigma_rho, 'time_dependent_DT_temp_axis': time_dependent_DT_temp_axis,
                    'pid_K' : pidK, 'pid_Ti' : pidI, 'pid_Td' : pidD, 'pid_I' : 0.0,
                    'noise_level' : noise_level, 'previous_error' : 0.0}
                    print(f'Using PID_itemp_gaussian with pid_K={pidK}m^3/s, pid_tauI={pidI}s, pid_tauD={pidD}s and noise_level={noise_level}')
            #
            case 'PID_etemp_gaussian':
                if((rho_0 is None) or (sigma_rho is None) or (max_total_power is None) or (time_dependent_electron_temp_axis is None)):
                    raise ValueError('ERROR: Need to provide time_dependent_electron_temp_axis, max_total_power, rho_0 and sigma_rho for PID electron density gaussian')
                else:
                    self.energy_sources[species][source_type] = {'max_total_power' : max_total_power, 'rho_0' : rho_0, 
                    'sigma_rho' : sigma_rho, 'time_dependent_electron_temp_axis': time_dependent_electron_temp_axis,
                    'pid_K' : pidK, 'pid_Ti' : pidI, 'pid_Td' : pidD, 'pid_I' : 0.0,
                    'noise_level' : noise_level, 'previous_error' : 0.0}
                    print(f'Using PID_etemp_gaussian with pid_K={pidK}m^3/s, pid_tauI={pidI}s, pid_tauD={pidD}s and noise_level={noise_level}')
            #
            case 'Coll_Heat_Exchange':
                self.energy_sources[species][source_type] = {}
                self.solve_coll_heat_exchange = True
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
            case 'lambda_2D':
                if(lambda_function_2D is None):
                    print('ERROR: A 2D (r,t) lambda funcion must be provided!')
                    exit(0)        
                # check it's a lambda function with two arguments
                num_args = len(inspect.signature(lambda_function_2D).parameters)
                if( not callable(lambda_function_2D) or num_args!=2 ):     
                    print('A lambda function with two arguments, (r,t), must be given')
                    exit(0)
                else:
                    self.energy_sources[species][source_type] = {'lambda_function_2D' : lambda_function_2D}
            #
            case _:
                print(f'ERROR: Source type {source_type} is NOT possible')
                exit(0)
                
    def set_particle_source(self,species,source_type, injected_particles_per_sec=None, rho_0=None, sigma_rho=None, 
        cte_source=None, time_dependent_factor=None, lambda_function_2D=None,
        max_injected_particles_per_sec=None, time_dependent_electron_dens_axis=None, time_dependent_fusion_power=None, time_dependent_DT_temp_axis=None,
        pidK=1.0,pidI=1.0E10,pidD=0.0,noise_level=0.0):
        """
        Sets particle sources for a given species. The source_type can be:
        'external_gaussian', 'time_dependent_gaussian', 'constant', 'lambda_2D',
        'fast_alphas_source' (for He-4) and 'alpha_particles_sink' (for D and T),
        'pid_edense_gaussian' (for PID electron density control with gaussian soure),
        'protium generation' (for hydrogen)
        'deuterium_sink_protium_generation' (for deuterium)
        """
        import inspect
        
        # check species exist in list_of_species
        if species not in self.list_of_species:
            raise ValueError(f"ERROR: Species {species} is not in the plasma.")
         
        match source_type:
            case 'external_gaussian':
                if((injected_particles_per_sec is None) or (sigma_rho is None) or (rho_0 is None)):
                    raise ValueError('ERROR: Need to provide injected_particles_per_sec, rho_0 and sigma_rho for gaussian external source')
                else:
                    self.particle_sources[species][source_type] = {'injected_particles_per_sec' : injected_particles_per_sec, 'rho_0' : rho_0, 'sigma_rho' : sigma_rho}
            #
            case 'time_dependent_gaussian':
                if((injected_particles_per_sec is None) or (rho_0 is None) or (sigma_rho is None) or (time_dependent_factor is None)):
                    raise ValueError('ERROR: Need to provide injected_particles_per_sec, rho_0, sigma_rho and a time depenedent factor for time-dependent gaussian')
                else:
                    self.particle_sources[species][source_type] = {'injected_particles_per_sec' : injected_particles_per_sec, 'rho_0' : rho_0, 'sigma_rho' : sigma_rho, 'time_factor': time_dependent_factor }
            #
            case 'PID_edense_gaussian':
                if((rho_0 is None) or (sigma_rho is None) or (max_injected_particles_per_sec is None) or (time_dependent_electron_dens_axis is None)):
                    raise ValueError('ERROR: Need to provide time_dependent_electron_dens_axis, max_injected_particles_per_sec, rho_0 and sigma_rho for PID electron density gaussian')
                else:
                    self.particle_sources[species][source_type] = {'max_injected_particles_per_sec' : max_injected_particles_per_sec, 'rho_0' : rho_0, 
                    'sigma_rho' : sigma_rho, 'time_dependent_electron_dens_axis': time_dependent_electron_dens_axis,
                    'pid_K' : pidK, 'pid_Ti' : pidI, 'pid_Td' : pidD, 'pid_I' : 0.0,
                    'noise_level' : noise_level, 'previous_error' : 0.0}
                    print(f'Using PID_edense_gaussian with pid_K={pidK}m^3/s, pid_tauI={pidI}s, pid_tauD={pidD}s and noise_level={noise_level}')
            #
            case 'PID_pfuse_gaussian':
                if((rho_0 is None) or (sigma_rho is None) or (max_injected_particles_per_sec is None) or (time_dependent_fusion_power is None)):
                    raise ValueError('ERROR: Need to provide time_dependent_fusion_power, max_injected_particles_per_sec, rho_0 and sigma_rho for PID electron density gaussian')
                else:
                    self.particle_sources[species][source_type] = {'max_injected_particles_per_sec' : max_injected_particles_per_sec, 'rho_0' : rho_0, 
                    'sigma_rho' : sigma_rho, 'time_dependent_fusion_power': time_dependent_fusion_power,
                    'pid_K' : pidK, 'pid_Ti' : pidI, 'pid_Td' : pidD, 'pid_I' : 0.0,
                    'noise_level' : noise_level, 'previous_error' : 0.0}
                    print(f'Using PID_pfuse_gaussian with pid_K={pidK}m^3/s, pid_tauI={pidI}s, pid_tauD={pidD}s and noise_level={noise_level}')
            #
            case 'PID_itemp_gaussian':
                if((rho_0 is None) or (sigma_rho is None) or (max_injected_particles_per_sec is None) or (time_dependent_DT_temp_axis is None)):
                    raise ValueError('ERROR: Need to provide time_dependent_DT_temp_axis, max_total_power, rho_0 and sigma_rho for PID electron density gaussian')
                else:
                    self.particle_sources[species][source_type] = {'max_injected_particles_per_sec' : max_injected_particles_per_sec , 'rho_0' : rho_0, 
                    'sigma_rho' : sigma_rho, 'time_dependent_DT_temp_axis': time_dependent_DT_temp_axis,
                    'pid_K' : pidK, 'pid_Ti' : pidI, 'pid_Td' : pidD, 'pid_I' : 0.0,
                    'noise_level' : noise_level, 'previous_error' : 0.0}
                    print(f'Using PID_itemp_gaussian with pid_K={pidK}m^3/s, pid_tauI={pidI}s, pid_tauD={pidD}s and noise_level={noise_level}')
            #
            case 'PID_pradfrac_gaussian':
                if((rho_0 is None) or (sigma_rho is None) or (time_dependent_factor is None)):
                    raise ValueError('ERROR: Need to provide injected_particles_per_sec, rho_0, sigma_rho and a time depenedent factor for PID radiated fraction gaussian')
                else:
                    self.particle_sources[species][source_type] = {'injected_particles_per_sec' : injected_particles_per_sec, 'rho_0' : rho_0, 
                    'sigma_rho' : sigma_rho, 'time_factor': time_dependent_factor,
                    'pid_K' : pidK, 'pid_Ti' : pidI, 'pid_Td' : pidD, 'pid_I' : 0.0,
                    'noise_level' : noise_level, 'previous_error' : 0.0}
            #
            case 'fast_alphas_source':
                # check we are solving fast alphas
                if(not self.solve_fast_alphas):
                    raise ValueError('solve_fast_alphas was set to false, so fast_alphas_source does not make sense...')
                if(species != 'helium4'): 
                    raise ValueError('ERROR: fast_alphas_source is only source for helium4 (thermal helium)')
                self.particle_sources[species][source_type] = {}
            #
            case 'alpha_particles_sink':
                if(species != 'deuterium' and species != 'tritium'):
                    raise ValueError('ERROR: alpha_particles_sink is only source for deuterium and tritium')
                self.particle_sources[species][source_type] = {}
            #
            case 'protium_generation':
                if(species != 'hydrogen'):
                    raise ValueError('ERROR: protium generation is only source for hydrogen')
                self.particle_sources[species][source_type] = {}
            case 'deuterium_sink_protium_generation':
                if(species != 'deuterium'):
                    raise ValueError('ERROR: deuterium sink protium generation is only source for deuterium')
                self.particle_sources[species][source_type] = {}
            #
            case 'constant':
                if(cte_source is None):
                    raise ValueError('ERROR: cte_source is needed in order to generate a constant source.')
                else:
                    self.particle_sources[species][source_type] = {'cte_source' : cte_source}
            #
            case 'lambda_2D':
                if(lambda_function_2D is None):
                    raise ValueError('ERROR: A 2D (r,t) lambda funcion must be provided!')      
                # check it's a lambda function with two arguments
                num_args = len(inspect.signature(lambda_function_2D).parameters)
                if( not callable(lambda_function_2D) or num_args!=2 ):     
                    raise ValueError('A lambda function with two arguments, (r,t), must be given')
                else:
                    self.particle_sources[species][source_type] = {'lambda_function_2D' : lambda_function_2D}
            #
            case _:
                raise ValueError(f'ERROR: Source type {source_type} is NOT possible')
                
    def set_heat_fluxes(self, type: str,surfaces=None,chi=None,chi_base=None,aLT_critical=None,alpha=None,stiffness=None,chi_electrons=None,convective_fact=None,mass_ref_species=None):
        """
        Sets heat flux for all species. In general, the heat flux for each species is:
        Q = -n chi dT/dr + convective_fact*T*Gamma_turb
        The 'type' argument will set how chi is computed:
        'diffusive' : 'chi' is constant and equal to all species
        'beurskens' : 'chi_e' is constant for electrons; chi_ions are computed according to Beurskens model
                       using aLT_critical, stiffness and alpha (can be functions of rho)
        """
            
        match type:
            case 'dkespenta':
                raise ValueError('dkespenta heat fluxes not working')
                # #checks that dkes_folder and surfaces are provided
                # if( surfaces is None):
                #     print('ERROR: surfaces must be provided!!')
                #     exit(0)
                    
                # self.heat_fluxes_info['type'] = type
                # self.heat_fluxes_info[type]['surfaces'] = surfaces
                
            case 'diffusive':
                #checks that diffusion coefficients are provided
                if((chi is None) or (convective_fact is None)):
                    raise ValueError('ERROR: chi and convective_fact must be provided!')
                self.heat_fluxes_info['type'] = type
                self.heat_fluxes_info[type]['chi'] = chi
                self.heat_fluxes_info[type]['convective_fact'] = convective_fact
                
            case 'beurskens':
                if( (chi_base is None) or (aLT_critical is None) or (alpha is None) or (stiffness is None) or (chi_electrons is None) or (convective_fact is None) or (mass_ref_species is None)):
                    raise ValueError('ERROR: chi_base, aLT_critical, alpha, stiffness, chi_electrons, convective_fact and mass_ref_species must be given!')
                self.heat_fluxes_info['type'] = type
                self.heat_fluxes_info[type]['chi_base'] = chi_base
                self.heat_fluxes_info[type]['aLT_critical'] = aLT_critical
                self.heat_fluxes_info[type]['alpha'] = alpha
                self.heat_fluxes_info[type]['stiffness'] = stiffness
                self.heat_fluxes_info[type]['chi_electrons'] = chi_electrons
                self.heat_fluxes_info[type]['convective_fact'] = convective_fact
                self.heat_fluxes_info[type]['mass_ref_species'] = mass_ref_species        
                
            case 'dkespenta_beurskens':
                raise ValueError('dkespenta_beurskens heat fluxes not working')
                # if( (surfaces is None) or (chi_base is None) or (aLT_critical is None) or (alpha is None) or (stiffness is None) or (chi_electrons is None) or (convective_fact is None)):
                #     print('ERROR: surfaces, chi_base, aLT_critical, alpha, stiffness, chi_electrons and convective_fact must be given!')
                #     exit(0)
                # self.heat_fluxes_info['type'] = type
                # self.heat_fluxes_info[type]['surfaces'] = surfaces
                # self.heat_fluxes_info[type]['chi_base'] = chi_base
                # self.heat_fluxes_info[type]['aLT_critical'] = aLT_critical
                # self.heat_fluxes_info[type]['alpha'] = alpha
                # self.heat_fluxes_info[type]['stiffness'] = stiffness
                # self.heat_fluxes_info[type]['chi_electrons'] = chi_electrons
                # self.heat_fluxes_info[type]['convective_fact'] = convective_fact
                
    def set_particle_fluxes(self, type: str, surfaces=None, Dn=None, cn=None, mass_ref_species=None):
        """
        Sets particle flux for all species. Currently, only 'diffusive' type is implemented:
        Gamma = -Dn dn/dr with Dn constant and equal to all species
        """
            
        match type:
            case 'dkespenta':
                raise ValueError('dkespenta particle flux not working')
                # #checks that dkes_folder and surfaces are provided
                # if( surfaces is None ):
                #     print('ERROR:surfaces must be provided!!')
                #     exit(0)
                    
                # self.particle_fluxes_info['type'] = type
                # self.particle_fluxes_info[type]['surfaces'] = surfaces
                
            case 'diffusive':
                #checks that diffusion coefficient is provided
                if(Dn is None):
                    raise ValueError('ERROR: Dn must be provided!')
                self.particle_fluxes_info['type'] = type
                self.particle_fluxes_info[type]['Dn'] = Dn
                
            case 'diffusive_advective':
                #checks that Dn and cn are provided. Dn and cn can be floats or functions of rho
                if(Dn is None or cn is None):
                    raise ValueError('ERROR: Dn and cn must be provided!')
                # If are functions, check they are 1D functions
                if( callable(Dn) ):
                    import inspect
                    num_args = len(inspect.signature(Dn).parameters)
                    if( num_args != 1 ):
                        raise ValueError('ERROR: Dn must be a 1D function of rho')
                #
                if( callable(cn) ):
                    import inspect
                    num_args = len(inspect.signature(cn).parameters)
                    if( num_args != 1 ):
                        raise ValueError('ERROR: cn must be a 1D function of rho')
                #  
                self.particle_fluxes_info['type'] = type
                self.particle_fluxes_info[type]['Dn'] = Dn
                self.particle_fluxes_info[type]['cn'] = cn
                
            case 'normalized_Dn_cn_rho_aLn_dependent':
                if( (Dn is None) or (cn is None) or (mass_ref_species is None)):
                    raise ValueError('ERROR: Dn, cn and mass_ref_species must be given!')
                # checks that Dn and cn are 2D function:
                if( callable(Dn) and callable(cn) ):
                    import inspect
                    num_args = len(inspect.signature(Dn).parameters)
                    if( num_args !=2 ):
                        raise ValueError('ERROR: Dn must be a function of (rho,aLn)'   )
                    num_args = len(inspect.signature(cn).parameters)
                    if( num_args !=2 ):
                        raise ValueError('ERROR: cn must be a function of (rho,aLn)'   )
                else:
                    raise ValueError('ERROR: Dn and cn must be functions of (rho,aLn)'   )
                #
                self.particle_fluxes_info['type'] = type
                self.particle_fluxes_info[type]['Dn'] = Dn
                self.particle_fluxes_info[type]['cn'] = cn             
                self.particle_fluxes_info[type]['mass_ref_species'] = mass_ref_species             
                
    def run(self,Nr,dt,tstart,tend,tolerance=1E-2,max_subiter=12,output_filename=None,restart_filename=None,dt_save=0.1):
        """
        Run the transport solver after setting initial profiles, BCs, fluxes types and sources
        If output_filename is not None, then results will be saved in a joblib file
        If restart_filename is not None, overwrites eventually already defined initial profiles
        and BCs using profiles of the last iteration in restart file 
        The restart file is simply the output file of the previous simulation
        """
        
        # Updates initial conditions and boundary conditions with profiles in restart file
        if(restart_filename is not None):
            self.read_restart_file(restart_filename,tstart)
        
        # Check everything is set and ready to proceed with the run
        self.make_checks()
        
        # Initialize rho grid
        rho = np.linspace(0,1,Nr)
        drho = rho[1]-rho[0]

        # Initialize time grid
        Nt = round( 1+(tend-tstart)/dt )
        time = np.linspace(tstart,tend,Nt)
            
        # save grids in class
        self.time = time
        self.tstart = tstart
        self.tend = tend
        self.Nt = Nt
        self.dt = dt
        self.rho_grid = rho
        self.r_grid = rho * self.aminor
        self.drho = drho
        self.dr = drho * self.aminor
        self.Nr = Nr
        
        if(self.add_NEO):
            if self.dt_NEO is None:
                self.nsteps_per_NEO = 1
            else:
                self.nsteps_per_NEO = round(self.dt_NEO/self.dt)
                
        # print grid details
        self.print_grid_details()
        
        # initialize self.## variables
        self.initialize_variables()
        
        # setup types of fluxes
        self.setup_fluxes_type()
        
        # set fields at t=tstart
        fields_old = self.set_fields_tstart()
    
        start_time = perf_counter() 
        ### LOOP IN TIME STARTING AT t=tstart+dt ###
        for it,t in enumerate(time[1:],start=1):
            self.it = it
            
            # the first subiter corresponds to the last time iteration
            for species in self.list_of_species:
                self.N[species][it,:] = self.N[species][it-1,:]
                self.T[species][it,:] = self.T[species][it-1,:]
                self.P[species][it,:] = self.P[species][it-1,:]
                if(self.solve_fast_alphas):
                    self.N['alphas_fast'][it,:] = self.N['alphas_fast'][it-1,:]
                if(self.add_NEO):
                    self.Er[it,:] = self.Er[it-1,:]   
            
            ### SUBCYCLE
            delta_p = 10*tolerance
            subiter=1
            while(delta_p > tolerance and subiter<=max_subiter):
                self.subiter = subiter  

                # Call NEO
                if(self.add_NEO):
                    if it%self.nsteps_per_NEO==0: 
                        self.call_NEO(it)
                    else:
                        self.set_NEO_coefficients_from_previous(it)
                      
                # compute fluxes and diffusion coefficients
                self.call_fluxes(it)
                
                # set explicit sources
                for species in self.list_of_species:
                    self.set_explicit_energy_sources(species,it)
                    self.set_explicit_particle_sources(species,it)
                    
                # solve
                dens = self.solve_density_equations(it) # self.N[species][it,:] are updated inside this function
                press = self.solve_pressure_equations(it) # self.P[species][it,:] and self.T[species][it,:] are updated inside this function
                
                fields = np.concatenate((dens,press))
                delta_p = np.max( np.where( fields_old>1E-10, np.abs((fields-fields_old)/fields_old), 0 ) )
                
                fields_old = fields
                
                ion_info = self.plasma.ion_species[0]
                info_str = f"  {t:<13.3f}{subiter:<10}{self.T['electrons'][it,0]/1E3:<18.3f}{self.N['electrons'][it,0]:<20.2E}{self.T[ion_info][it,0]/1E3:<18.3f}{self.N[ion_info][it,0]:<20.2E}{delta_p:<13.2E}"
                print(info_str)
                
                subiter += 1
        
        if(output_filename is not None):
            self.call_save_output(output_filename,dt_save)

        # Shut down the persistent NEO worker pool (if any), so no worker
        # processes are left running once the solve is done.
        if(getattr(self,'add_NEO',False) and self.neo_pool is not None):
            self.neo_pool.shutdown()

        end_time = perf_counter()
        print(f'Plasma Solver took {(end_time-start_time)/60:.2f}min to run.')
    
    def make_checks(self):
        """
        Before running the simulation, checks if dVdr, initial profiles, BCs 
        and fluxes types have been set
        Checks consistency between initial profiles and Dirichlet BCs
        """
        
        # check equilibrium exists
        if(not hasattr(self,'dVdr')):
            raise ValueError('ERROR: dVdr MUST BE SET!!')
        
        for species in self.list_of_species:
            
            # check boundary conditions are set
            if (species not in self.edge_density_BC or species not in self.edge_pressure_BC):
                raise KeyError(f"Missing edge boundary condition for species: {species}")
            
            # check initial profiles are set
            if (species not in self.initial_density or species not in self.initial_pressure):
                raise KeyError(f"Missing initial profile for species: {species}")
            
            # check consistency between boundary conditions and initial profiles
            if(not np.isclose(self.initial_density[species](1),self.edge_density_BC[species],rtol=1E-7, atol=1E-12)):
                raise ValueError(f'Edge density BC not consistent w/ initial density profile')
            if(not np.isclose(self.initial_pressure[species](1),self.edge_pressure_BC[species],rtol=1E-7, atol=1E-12)):
                raise ValueError(f'Edge pressure/temperature BC not consistent w/ initial temperature profile')
            
            # tol = np.abs(self.edge_density_BC[species]) * np.finfo(float).eps
            # if( np.abs(self.initial_density[species](1)-self.edge_density_BC[species]) > 5*tol ):
            #     raise ValueError(f'Edge density BC not consistent w/ initial density profile')
            # tol = np.abs(self.edge_pressure_BC[species]) * np.finfo(float).eps
            # if( np.abs(self.initial_pressure[species](1)-self.edge_pressure_BC[species]) > 10*tol ):
            #     raise ValueError(f'Edge pressure/temperature BC not consistent w/ initial temperature profile')
            
        # check fluxes info is set
        if(not hasattr(self,'heat_fluxes_info')):
            raise KeyError('ERROR: set_heat_fluxes must be called before running!!')
        if(not hasattr(self,'particle_fluxes_info')):
            raise KeyError('ERROR: set_particle_fluxes must be called before running!!')
        
        if(self.add_NEO):
            # Check if initialize_NEO has been called
            if not hasattr(self,'DKES_nuv'):
                raise ValueError('!!! add_NEO=True but initialize_NEO has not been called !!! ')
            
    def print_grid_details(self):
        """
        Prints to the command line info about the grids and the simulation header
        """
        
        print(' ')
        print( ' ***********************')
        print(f' *  tstart = {self.tstart:5.2f}s    *')
        print(f' *  tend   = {self.tend:5.2f}s    *')
        print(f' *  dt     = {self.dt:5.3f}s    *')
        print(f' *  Nt     = {self.Nt:3}       *')
        print(f' *  drho   = {self.drho:5.3f}     *')
        if(self.add_NEO):
            print(f' *  dt_NEO = {self.dt*self.nsteps_per_NEO:5.3f}s    *')
        print( ' ***********************')
        
        print(' ')
        header_str = '  TIME [s]     NSUB      TE_AXIS [keV]     NE_AXIS [m^-3]    TI_AXIS [keV]    NI_AXIS [m^-3]    MAX(dp/p_old)'  
        print(header_str)
        print('  '+'='*len(header_str))
        
    def initialize_variables(self):
        """
        Initializes all dictionaries and arrays needed to run the simulation
        Initializes the LHS sparse matrices to solve the density and pressure equations
        """
        from collections import defaultdict
        
        Nr = self.Nr
        Nt = self.Nt
        
        self.N = {}
        self.P = {}
        self.T = {}
        #
        self.Q_turb = {}
        self.Dp = {}
        self.Dp_keep = defaultdict(lambda: defaultdict(list))
        self.cp = {}
        #
        self.Gamma_turb = {}
        self.Dn = {}
        self.cn = {}
        #
        self.Q_NEO = {}
        self.Gamma_NEO = {}
        #
        self.explicit_energy_sources = {}
        self.explicit_particle_sources = {}
        #
        if(self.add_NEO):
            self.Dn_NEO = {}
            self.cn_NEO = {}
            self.Dp_NEO = {}
            self.cp_NEO = {}
            
        self.Er = np.zeros((Nt,Nr))

        for species in self.list_of_species:
            
            self.P[species] = np.zeros((Nt,Nr))   
            self.T[species] = np.zeros((Nt,Nr)) 
            self.N[species] = np.zeros((Nt,Nr))
            
            self.Dp[species] = np.zeros((Nt,Nr)) 
            self.cp[species] = np.zeros((Nt,Nr)) 
            self.Q_turb[species] = np.zeros((Nt,Nr)) 
            self.Gamma_turb[species] = np.zeros((Nt,Nr)) 
            self.Dn[species] = np.zeros((Nt,Nr)) 
            self.cn[species] = np.zeros((Nt,Nr)) 
            self.Q_NEO[species] = np.zeros((Nt,Nr)) 
            self.Gamma_NEO[species] = np.zeros((Nt,Nr))
            
            if(self.add_NEO):
                self.Dn_NEO[species] = np.zeros((Nt,Nr)) 
                self.cn_NEO[species] = np.zeros((Nt,Nr)) 
                self.Dp_NEO[species] = np.zeros((Nt,Nr)) 
                self.cp_NEO[species] = np.zeros((Nt,Nr)) 
                
            
            self.explicit_energy_sources[species] = {}
            for source_type in self.energy_sources[species].keys():
                self.explicit_energy_sources[species][source_type] = np.zeros((Nt,Nr))
            self.explicit_particle_sources[species] = {}
            for source_type in self.particle_sources[species].keys():
                self.explicit_particle_sources[species][source_type] = np.zeros((Nt,Nr))
            
        if(self.solve_fast_alphas):
            self.N['alphas_fast'] = np.zeros((Nt,Nr))
            
        # Initialize tridiagonal matrices which will be used as LHS of density and pressure equations
        self.LHS_pressure,self.pressure_lower_block, self.pressure_main_block, self.pressure_upper_block = initialize_LHS_pressure(Nr,len(self.list_of_species))
        self.LHS_density,self.density_main_diag,self.density_lower_diag,self.density_upper_diag = initialize_LHS_density(Nr)
            
    def set_fields_tstart(self):
        """
        Sets density, temperature and pressure at t=tstart according to the initial profiles,
        which were either given through set_initial_profiles or through the restart file
        Sets explicit sources and fluxes at t=tstart
        """
        
        rho_grid = self.rho_grid
        
        for species in self.list_of_species:
        
            self.P[species][0,:] = self.initial_pressure[species](rho_grid)     
            self.N[species][0,:] = self.initial_density[species](rho_grid)
               
            self.T[species][0,:] = self.P[species][0,:] / (EC*self.N[species][0,:])
        
        # N_alphas are set to 0.0 at t=tstart unless read from restart file
        if(self.solve_fast_alphas):
            try:
                self.N['alphas_fast'][0,:] = self.alphas_fast_density_restart
                print('Reading alphas density from restart file...')
            except:
                self.N['alphas_fast'][0,:] = 0.0
            
        # set sources at t=0
        for species in self.list_of_species:
                self.set_explicit_energy_sources(species,it=0)
                self.set_explicit_particle_sources(species,it=0)
                
        # set NEO coefficients
        if(self.add_NEO):
            self.call_NEO(it=0)
        
        # set fluxes at t=0
        self.call_fluxes(it=0)  
        
        # return array with (dens,press) at t=tstart
        all_fields = []
        for species in self.list_of_species:
            all_fields.append(self.N[species][0,:])
        for species in self.list_of_species:
            all_fields.append(self.P[species][0,:])
        
        return np.concatenate(all_fields)
    
    def setup_fluxes_type(self):
        """
        Selects the flux computation functions based on 
        self.particle_fluxes_info['type'] and self.heat_fluxes_info['type']
        """
        
        ptype = self.particle_fluxes_info['type']
        htype = self.heat_fluxes_info['type']

        if (ptype in ['dkespenta'] or htype in ['dkespenta', 'dkespenta_beurskens']):
            raise ValueError('The DKES+PENTA functionality is broken!!')

        # Particle flux function
        if ptype == 'diffusive':
            self.particle_flux_func = self.compute_diffusive_particle_flux
        elif ptype == 'diffusive_advective':
            self.particle_flux_func = self.compute_diffusive_advective_particle_flux
        elif ptype == 'dkespenta':
            self.particle_flux_func = self.compute_NEO_particle_flux
        elif ptype == 'normalized_Dn_cn_rho_aLn_dependent':
            self.particle_flux_func = self.compute_normalized_Dn_cn_rho_aLn_dependent_particle_flux
        else:
            raise ValueError('Unsupported particle flux type')

        # Heat flux function
        if htype == 'diffusive':
            self.heat_flux_func = self.compute_diffusive_heat_flux
        elif htype == 'beurskens':
            self.heat_flux_func = self.compute_beurskens_heat_flux
        elif htype == 'dkespenta':
            self.heat_flux_func = self.compute_NEO_heat_flux
        elif htype == 'dkespenta_beurskens':
            self.heat_flux_func = self.compute_NEO_plus_beurskens_heat_flux
        else:
            raise ValueError('Unsupported heat flux type')
             
    def call_fluxes(self,it):
        """Computes fluxes at iteration it"""
           
        # Particle Fluxes
        self.particle_flux_func(it)
        # Heat Fluxes
        self.heat_flux_func(it)
        
        if(self.add_NEO):
            # add NEO fluxes to self.Dp, self.cp, self.Dn, self.cn
            self.add_NEO_transport_coefficients(it)
            
        
    def set_explicit_energy_sources(self,species: str, it):
        """
        Computes total explicit energy source of a given species
        at iteration it using info in self.energy_sources[species]
        Returns 1D-array of same size as rho_grid
        """
        from libstell.fusion import FUSION  
        from numpy.random import rand
        fusion = FUSION()

        rho_grid = self.rho_grid
        
        for source_type in self.energy_sources[species]:
            
            aux_source = 0.0
                 
            match source_type:
                case 'Bremsstrahlung':
                    ne = self.N['electrons'][it,:]
                    Te = self.T['electrons'][it,:]
                    for ion in self.plasma.ion_species:
                        zi = self.plasma.Zcharge[ion]
                        ni = self.N[ion][it,:]

                        aux_source -= fusion.BremsstrahlungPower(zi,ni,ne,Te)
                        
                case 'external_gaussian':
                    rho_0 = self.energy_sources[species]['external_gaussian']['rho_0']
                    r0 = rho_0 * self.aminor
                    sigma_rho = self.energy_sources[species]['external_gaussian']['sigma_rho']
                    sigma_r = sigma_rho*self.aminor
                    r = self.rho_grid * self.aminor
                    P_IN = self.energy_sources[species]['external_gaussian']['total_power']
                    #
                    integrand = np.exp(-(r-r0)**2/sigma_r**2) * self.dVdr(self.rho_grid)
                    integrand = integrand.flatten()
                    #
                    cte = P_IN / np.trapz(integrand,r)
                    #
                    aux_source = cte * np.exp(-(r-r0)**2/sigma_r**2)
                    
                case 'time_dependent_gaussian':
                    rho_0 = self.energy_sources[species]['time_dependent_gaussian']['rho_0']
                    r0 = rho_0 * self.aminor
                    sigma_rho = self.energy_sources[species]['time_dependent_gaussian']['sigma_rho']
                    sigma_r = sigma_rho*self.aminor
                    r = self.rho_grid * self.aminor
                    P_IN = self.energy_sources[species]['time_dependent_gaussian']['total_power']
                    #
                    integrand = np.exp(-(r-r0)**2/sigma_r**2) * self.dVdr(self.rho_grid)
                    integrand = integrand.flatten()
                    #
                    cte = P_IN / np.trapz(integrand,r)
                    #
                    time_fact = self.energy_sources[species]['time_dependent_gaussian']['time_factor']
                    t = self.time[it]
                    aux_source = cte * np.exp(-(r-r0)**2/sigma_r**2) * time_fact(t)
                    
                case 'Coll_Heat_Exchange':
                    # this source is fully implicit and is computed inside get_LHS_pressure
                    aux_source = 0.0
                            
                case 'alpha_heating':
                    nD = self.N['deuterium'][it,:]
                    nT = self.N['tritium'][it,:]

                    fraction_alpha_heating = self.energy_sources[species]['alpha_heating']['fraction_alpha_heating']
                    
                    TD = self.T['deuterium'][it,:]
                    TT = self.T['tritium'][it,:]
                    aux_source = fraction_alpha_heating * fusion.alphaPower(nD,nT,TD,TT)
                    
                case 'lambda_2D':
                    lambda_function_2D = self.energy_sources[species][source_type]['lambda_function_2D'] #func(r,t)
                    #
                    aux_source = [lambda_function_2D(r,self.time[it]) for r in self.r_grid]
                    
                case 'Er':
                    # aux_source = self.plasma.charge[species]*self.Er_interp(rho_grid)*self.Gamma_interp[species](rho_grid)
                    raise ValueError('Er energy source only available when code coupled with dkespenta.')
                      
                case 'constant':
                    aux_source = self.energy_sources[species]['constant']['cte_source']

                case 'PID_pfuse_gaussian':
                    # These define the gaussian
                    rho_0 = self.energy_sources[species]['PID_pfuse_gaussian']['rho_0']
                    sigma_rho = self.energy_sources[species]['PID_pfuse_gaussian']['sigma_rho']
                    pid_K     = self.energy_sources[species]['PID_pfuse_gaussian']['pid_K']
                    pid_Ti    = self.energy_sources[species]['PID_pfuse_gaussian']['pid_Ti']
                    pid_Td    = self.energy_sources[species]['PID_pfuse_gaussian']['pid_Td']
                    Ival      = self.energy_sources[species]['PID_pfuse_gaussian']['pid_I']
                    power_max  = self.energy_sources[species]['PID_pfuse_gaussian']['max_total_power']
                    noise     = self.energy_sources[species]['PID_pfuse_gaussian']['noise_level']
                    error_old = self.energy_sources[species]['PID_pfuse_gaussian']['previous_error']
                    t         = self.time[it]
                    setpoint  = self.energy_sources[species]['PID_pfuse_gaussian']['time_dependent_fusion_power'](t) # Set Point
                    it1       = it #max(it - 1,2)
                    # Compute fusion power
                    nD = self.N['deuterium'][it1,:]
                    nT = self.N['tritium'][it1,:]
                    TD = self.T['deuterium'][it1,:]
                    TT = self.T['tritium'][it1,:]
                    integrand = fusion.alphaPower(nD,nT,TD,TT)*self.dVdr(rho_grid)
                    integrand = integrand.flatten()
                    p_val = max(np.trapezoid(integrand,self.r_grid),0.0)*5.0 #from alpha power to fusion power
                    if np.isnan(p_val): p_val = 0.0
                    p_val = p_val * (1.0 + (rand()-0.5)*2.0*noise)
                    # Run PID algorithm
                    control, error, Ival = self.pid_controller(setpoint, p_val, pid_K, pid_Ti, pid_Td, error_old, Ival, self.dt)
                    self.energy_sources[species]['PID_pfuse_gaussian']['pid_I'] = Ival
                    self.energy_sources[species]['PID_pfuse_gaussian']['previous_error'] = error
                    # Threshold control
                    control = np.round(control,-6) # round to nearest MW
                    control = np.clip(control,0,power_max)
                    # If clipped, then set previous integral to zero (anti wind-up)
                    if np.isclose(control,0) or np.isclose(control,power_max):
                        self.energy_sources[species]['PID_pfuse_gaussian']['pid_I'] = 0.0     
                    # Compute integrand
                    integrand = np.exp(-(rho_grid-rho_0)**2/sigma_rho**2) * self.dVdr(rho_grid)
                    integrand = integrand.flatten()
                    # 
                    cte = control / np.trapezoid(integrand,self.r_grid)
                    #
                    aux_source = cte * np.exp(-(rho_grid-rho_0)**2/sigma_rho**2)
                    
                case 'PID_etemp_gaussian':
                    # These define the gaussian
                    rho_0 = self.energy_sources[species]['PID_etemp_gaussian']['rho_0']
                    sigma_rho = self.energy_sources[species]['PID_etemp_gaussian']['sigma_rho']
                    pid_K     = self.energy_sources[species]['PID_etemp_gaussian']['pid_K']
                    pid_Ti    = self.energy_sources[species]['PID_etemp_gaussian']['pid_Ti']
                    pid_Td    = self.energy_sources[species]['PID_etemp_gaussian']['pid_Td']
                    Ival      = self.energy_sources[species]['PID_etemp_gaussian']['pid_I']
                    power_max = self.energy_sources[species]['PID_etemp_gaussian']['max_total_power']
                    noise     = self.energy_sources[species]['PID_etemp_gaussian']['noise_level']
                    error_old = self.energy_sources[species]['PID_etemp_gaussian']['previous_error']
                    t         = self.time[it]
                    setpoint  = self.energy_sources[species]['PID_etemp_gaussian']['time_dependent_electron_temp_axis'](t) # Set Point
                    it1       = it #max(it - 1,2)
                    p_val     = self.T['electrons'][it1,0] * (1.0 + (rand()-0.5)*2.0*noise)
                    # Run PID algorithm
                    control, error, Ival = self.pid_controller(setpoint, p_val, pid_K, pid_Ti, pid_Td, error_old, Ival, self.dt)
                    self.energy_sources[species]['PID_etemp_gaussian']['pid_I'] = Ival
                    self.energy_sources[species]['PID_etemp_gaussian']['previous_error'] = error
                    # Threshold control
                    control = np.round(control,-6) # round to nearest MW
                    control = np.clip(control,0,power_max)
                    # If clipped, then set previous integral to zero (anti wind-up)
                    if np.isclose(control,0) or np.isclose(control,power_max):
                        self.energy_sources[species]['PID_etemp_gaussian']['pid_I'] = 0.0  
                    # Compute integrand
                    integrand = np.exp(-(rho_grid-rho_0)**2/sigma_rho**2) * self.dVdr(rho_grid)
                    integrand = integrand.flatten()
                    # 
                    cte = control / np.trapezoid(integrand,self.r_grid)
                    #
                    aux_source = cte * np.exp(-(rho_grid-rho_0)**2/sigma_rho**2)
                    
                case 'PID_itemp_gaussian':
                    # These define the gaussian
                    rho_0     = self.energy_sources[species]['PID_itemp_gaussian']['rho_0']
                    sigma_rho = self.energy_sources[species]['PID_itemp_gaussian']['sigma_rho']
                    pid_K     = self.energy_sources[species]['PID_itemp_gaussian']['pid_K']
                    pid_Ti    = self.energy_sources[species]['PID_itemp_gaussian']['pid_Ti']
                    pid_Td    = self.energy_sources[species]['PID_itemp_gaussian']['pid_Td']
                    Ival      = self.energy_sources[species]['PID_itemp_gaussian']['pid_I']
                    power_max = self.energy_sources[species]['PID_itemp_gaussian']['max_total_power']
                    noise     = self.energy_sources[species]['PID_itemp_gaussian']['noise_level']
                    error_old = self.energy_sources[species]['PID_itemp_gaussian']['previous_error']
                    t         = self.time[it]
                    setpoint  = self.energy_sources[species]['PID_itemp_gaussian']['time_dependent_DT_temp_axis'](t) # Set Point
                    it1       = it # max(it - 1,2)
                    p_val          = 0.5*(self.T['deuterium'][it1,0]+self.T['tritium'][it1,0])* (1.0 + (rand()-0.5)*2.0*noise)
                    # Run PID algorithm
                    control, error, Ival = self.pid_controller(setpoint, p_val, pid_K, pid_Ti, pid_Td, error_old, Ival, self.dt)
                    self.energy_sources[species]['PID_itemp_gaussian']['pid_I'] = Ival
                    self.energy_sources[species]['PID_itemp_gaussian']['previous_error'] = error
                    # Threshold control
                    control = np.round(control,-6) # round to nearest MW
                    control = np.clip(control,0,power_max)
                    # If clipped, then set previous integral to zero (anti wind-up)
                    if np.isclose(control,0) or np.isclose(control,power_max):
                        self.energy_sources[species]['PID_itemp_gaussian']['pid_I'] = 0.0 
                    # Compute integrand
                    integrand = np.exp(-(rho_grid-rho_0)**2/sigma_rho**2) * self.dVdr(rho_grid)
                    integrand = integrand.flatten()
                    # 
                    cte = control / np.trapezoid(integrand,self.r_grid)
                    #
                    aux_source = cte * np.exp(-(rho_grid-rho_0)**2/sigma_rho**2)

                case _:
                    print(f'ERROR: Source type {source_type} not defined....')
                    exit(0)
                   
            self.explicit_energy_sources[species][source_type][it,:] = aux_source

    def set_explicit_particle_sources(self,species: str, it):
        """
        Computes total explicit particle source of a given species
        at iteration it using info in self.particle_sources[species]
        Returns 1D-array of same size as rho_grid
        """
        from libstell.fusion import FUSION
        from numpy.random import rand
        fusion = FUSION()
        
        rho_grid = self.rho_grid
        
        for source_type in self.particle_sources[species]:
            
            aux_source = 0.0
                     
            match source_type:
                case 'external_gaussian':
                    rho_0 = self.particle_sources[species]['external_gaussian']['rho_0']
                    sigma_rho = self.particle_sources[species]['external_gaussian']['sigma_rho']
                    injected_particles_per_sec = self.particle_sources[species]['external_gaussian']['injected_particles_per_sec']
                    #
                    integrand = np.exp(-(rho_grid-rho_0)**2/sigma_rho**2) * self.dVdr(rho_grid)
                    integrand = integrand.flatten()
                    #
                    cte = injected_particles_per_sec / np.trapz(integrand,self.r_grid)
                    #
                    aux_source = cte * np.exp(-(rho_grid-rho_0)**2/sigma_rho**2)
                    
                case 'time_dependent_gaussian':
                    rho_0 = self.particle_sources[species]['time_dependent_gaussian']['rho_0']
                    sigma_rho = self.particle_sources[species]['time_dependent_gaussian']['sigma_rho']
                    injected_particles_per_sec = self.particle_sources[species]['time_dependent_gaussian']['injected_particles_per_sec']
                    time_fact = self.particle_sources[species]['time_dependent_gaussian']['time_factor']
                    #
                    integrand = np.exp(-(rho_grid-rho_0)**2/sigma_rho**2) * self.dVdr(rho_grid)
                    integrand = integrand.flatten()
                    #
                    cte = injected_particles_per_sec / np.trapz(integrand,self.r_grid)
                    #
                    t = self.time[it]
                    aux_source = time_fact(t) * cte * np.exp(-(rho_grid-rho_0)**2/sigma_rho**2)
                    
                case 'alpha_particles_sink':
                    nD = self.N['deuterium'][it,:]
                    nT = self.N['tritium'][it,:]
                    
                    TD = self.T['deuterium'][it,:]
                    TT = self.T['tritium'][it,:]
                    
                    sigmav = fusion.sigmaBH(0.5*(TD+TT),'DT')
                    aux_source = - nD*nT*sigmav # particles/(s*m^3)
                    
                case 'fast_alphas_source':
                    aux_source = self.N['alphas_fast'][it,:] / self.tau_fast_alphas
                    
                case 'protium_generation':
                    nD = self.N['deuterium'][it,:]
                    TD = self.T['deuterium'][it,:]
                    
                    sigmav = fusion.sigmaBH(TD,'DDT')
                    aux_source = 0.5*nD*nD*sigmav # factor 1/2 due to like-particle collisions (see Freidberg for instance)
                    
                case 'deuterium_sink_protium_generation':
                    nD = self.N['deuterium'][it,:]
                    TD = self.T['deuterium'][it,:]
                    
                    sigmav = fusion.sigmaBH(TD,'DDT')
                    aux_source = - 2*0.5*nD*nD*sigmav # 2 D's disappear for each protium
                    
                case 'constant':
                    aux_source = self.particle_sources[species]['constant']['cte_source']
                    
                case 'lambda_2D':
                    lambda_function_2D = self.particle_sources[species][source_type]['lambda_function_2D'] #func(r,t)
                    #
                    aux_source = [lambda_function_2D(r,self.time[it]) for r in self.r_grid]
                    
                case 'PID_edense_gaussian':
                    rho_0 = self.particle_sources[species]['PID_edense_gaussian']['rho_0']
                    sigma_rho = self.particle_sources[species]['PID_edense_gaussian']['sigma_rho']
                    pid_K     = self.particle_sources[species]['PID_edense_gaussian']['pid_K']
                    pid_Ti    = self.particle_sources[species]['PID_edense_gaussian']['pid_Ti']
                    pid_Td    = self.particle_sources[species]['PID_edense_gaussian']['pid_Td']
                    Ival      = self.particle_sources[species]['PID_edense_gaussian']['pid_I']
                    N_IN_max  = self.particle_sources[species]['PID_edense_gaussian']['max_injected_particles_per_sec']
                    noise     = self.particle_sources[species]['PID_edense_gaussian']['noise_level']
                    error_old = self.particle_sources[species]['PID_edense_gaussian']['previous_error']
                    t         = self.time[it]
                    setpoint  = self.particle_sources[species]['PID_edense_gaussian']['time_dependent_electron_dens_axis'](t) # Set Point
                    it1       = it #max(it - 1,2)
                    p_val     = self.N['electrons'][it1,0] * (1.0 + (rand()-0.5)*2.0*noise)
                    # Run PID algorithm
                    control, error, Ival = self.pid_controller(setpoint, p_val, pid_K, pid_Ti, pid_Td, error_old, Ival, self.dt)
                    self.particle_sources[species]['PID_edense_gaussian']['pid_I'] = Ival
                    self.particle_sources[species]['PID_edense_gaussian']['previous_error'] = error
                    # Threshold control
                    control = np.clip(control,0,N_IN_max)
                    # If clipped, then set previous integral to zero (anti wind-up)
                    if np.isclose(control,0) or np.isclose(control,N_IN_max):
                        self.particle_sources[species]['PID_edense_gaussian']['pid_I'] = 0.0     
                    # Compute integrand
                    integrand = np.exp(-(rho_grid-rho_0)**2/sigma_rho**2) * self.dVdr(rho_grid)
                    integrand = integrand.flatten()
                    # 
                    cte = control / np.trapezoid(integrand,self.r_grid)
                    #
                    aux_source = cte * np.exp(-(rho_grid-rho_0)**2/sigma_rho**2)
                    
                case 'PID_pfuse_gaussian':
                    # These define the gaussian
                    rho_0 = self.particle_sources[species]['PID_pfuse_gaussian']['rho_0']
                    sigma_rho = self.particle_sources[species]['PID_pfuse_gaussian']['sigma_rho']
                    pid_K     = self.particle_sources[species]['PID_pfuse_gaussian']['pid_K']
                    pid_Ti    = self.particle_sources[species]['PID_pfuse_gaussian']['pid_Ti']
                    pid_Td    = self.particle_sources[species]['PID_pfuse_gaussian']['pid_Td']
                    Ival      = self.particle_sources[species]['PID_pfuse_gaussian']['pid_I']
                    N_IN_max  = self.particle_sources[species]['PID_pfuse_gaussian']['max_injected_particles_per_sec']
                    noise     = self.particle_sources[species]['PID_pfuse_gaussian']['noise_level']
                    error_old = self.particle_sources[species]['PID_pfuse_gaussian']['previous_error']
                    t         = self.time[it]
                    setpoint  = self.particle_sources[species]['PID_pfuse_gaussian']['time_dependent_fusion_power'](t) # Set Point
                    it1       = it #max(it - 1,2)
                    # Compute fusion power
                    nD = self.N['deuterium'][it1,:]
                    nT = self.N['tritium'][it1,:]
                    TD = self.T['deuterium'][it1,:]
                    TT = self.T['tritium'][it1,:]
                    integrand = fusion.alphaPower(nD,nT,TD,TT)*self.dVdr(rho_grid)
                    integrand = integrand.flatten()
                    p_val = max(np.trapezoid(integrand,self.r_grid),0.0)*5.0 #from alpha power to fusion power
                    if np.isnan(p_val): p_val = 0.0
                    p_val = p_val * (1.0 + (rand()-0.5)*2.0*noise)
                    # Run PID algorithm
                    control, error, Ival = self.pid_controller(setpoint, p_val, pid_K, pid_Ti, pid_Td, error_old, Ival, self.dt)
                    self.particle_sources[species]['PID_pfuse_gaussian']['pid_I'] = Ival
                    self.particle_sources[species]['PID_pfuse_gaussian']['previous_error'] = error
                    # Threshold control
                    control = np.clip(control,0,N_IN_max)
                    # If clipped, then set previous integral to zero (anti wind-up)
                    if np.isclose(control,0) or np.isclose(control,N_IN_max):
                        self.particle_sources[species]['PID_pfuse_gaussian']['pid_I'] = 0.0
                    # Compute integrand
                    integrand = np.exp(-(rho_grid-rho_0)**2/sigma_rho**2) * self.dVdr(rho_grid)
                    integrand = integrand.flatten()
                    # 
                    cte = control / np.trapezoid(integrand,self.r_grid)
                    #
                    aux_source = cte * np.exp(-(rho_grid-rho_0)**2/sigma_rho**2)
                
                case 'PID_itemp_gaussian':
                    # These define the gaussian
                    rho_0     = self.particle_sources[species]['PID_itemp_gaussian']['rho_0']
                    sigma_rho = self.particle_sources[species]['PID_itemp_gaussian']['sigma_rho']
                    pid_K     = self.particle_sources[species]['PID_itemp_gaussian']['pid_K']
                    pid_Ti    = self.particle_sources[species]['PID_itemp_gaussian']['pid_Ti']
                    pid_Td    = self.particle_sources[species]['PID_itemp_gaussian']['pid_Td']
                    Ival      = self.particle_sources[species]['PID_itemp_gaussian']['pid_I']
                    N_IN_max  = self.particle_sources[species]['PID_itemp_gaussian']['max_injected_particles_per_sec']
                    noise     = self.particle_sources[species]['PID_itemp_gaussian']['noise_level']
                    error_old = self.particle_sources[species]['PID_itemp_gaussian']['previous_error']
                    t         = self.time[it]
                    setpoint  = self.particle_sources[species]['PID_itemp_gaussian']['time_dependent_DT_temp_axis'](t) # Set Point
                    it1       = it #max(it - 1,2)
                    p_val     = 0.5*(self.T['deuterium'][it1,0]+self.T['tritium'][it1,0])* (1.0 + (rand()-0.5)*2.0*noise)
                    # Run PID algorithm
                    control, error, Ival = self.pid_controller(setpoint, p_val, pid_K, pid_Ti, pid_Td, error_old, Ival, self.dt)
                    self.particle_sources[species]['PID_itemp_gaussian']['pid_I'] = Ival
                    self.particle_sources[species]['PID_itemp_gaussian']['previous_error'] = error
                    # If clipped, then set previous integral to zero (anti wind-up)
                    if np.isclose(control,0) or np.isclose(control,N_IN_max):
                        self.particle_sources[species]['PID_edense_gaussian']['pid_I'] = 0.0                  
                    # Threshold control
                    control = np.clip(control,0,N_IN_max)
                    # Compute integrand
                    integrand = np.exp(-(rho_grid-rho_0)**2/sigma_rho**2) * self.dVdr(rho_grid)
                    integrand = integrand.flatten()
                    # 
                    cte = control / np.trapezoid(integrand,self.r_grid)
                    #
                    aux_source = cte * np.exp(-(rho_grid-rho_0)**2/sigma_rho**2)
                    
                case 'PID_pradfrac_gaussian':
                    # These define the gaussian
                    rho_0     = self.particle_sources[species]['PID_pradfrac_gaussian']['rho_0']
                    sigma_rho = self.particle_sources[species]['PID_pradfrac_gaussian']['sigma_rho']
                    time_fact = self.particle_sources[species]['PID_pradfrac_gaussian']['time_factor']
                    pid_K     = self.particle_sources[species]['PID_pradfrac_gaussian']['pid_K']
                    pid_Ti    = self.particle_sources[species]['PID_pradfrac_gaussian']['pid_Ti']
                    pid_Td    = self.particle_sources[species]['PID_pradfrac_gaussian']['pid_Td']
                    Ival      = self.particle_sources[species]['PID_pradfrac_gaussian']['pid_I']
                    N_IN      = self.particle_sources[species]['PID_pradfrac_gaussian']['injected_particles_per_sec']
                    noise     = self.particle_sources[species]['PID_pradfrac_gaussian']['noise_level']
                    error_old = self.particle_sources[species]['PID_pradfrac_gaussian']['previous_error']
                    t              = self.time[it]
                    setpoint       = self.particle_sources[species]['PID_pradfrac_gaussian']['time_factor'](t) # Set Point
                    it1            = max(it - 1,2)
                    dt             = self.dt
                    # Compute ECRH power
                    integrand = 0.0
                    for ttype in ['external_gaussian','time_dependent_gaussian','PID_etemp_gaussian']:
                        if ttype in self.explicit_energy_sources['electrons'].keys():
                            integrand += self.explicit_energy_sources['electrons'][ttype][it1,:]
                    P_ECRH = max(np.trapezoid(integrand,self.r_grid),0.0)
                    if np.isnan(P_ECRH): P_ECRH = 0.0
                    # Compute Bremstahlung
                    integrand = 0.0
                    ne = self.N['electrons'][it1,:]
                    Te = self.T['electrons'][it1,:]
                    for ion in self.plasma.ion_species:
                        zi = self.plasma.Zcharge[ion]
                        ni = self.N[ion][it1,:]
                        integrand += fusion.BremsstrahlungPower(zi,ni,ne,Te)
                    P_BREM = max(np.trapezoid(integrand,self.r_grid),0.0)
                    if np.isnan(P_BREM): P_BREM = 0.0
                    # Compute fusion power
                    P_ALPHA = 0.0
                    if False:
                        nD = self.N['deuterium'][it1,:]
                        nT = self.N['tritium'][it1,:]
                        TD = self.T['deuterium'][it1,:]
                        TT = self.T['tritium'][it1,:]
                        integrand = fusion.alphaPower(nD,nT,TD,TT)*self.dVdr(rho_grid)
                        integrand = integrand.flatten()
                        P_ALPHA = max(np.trapezoid(integrand,self.r_grid),0.0)*5.0 #alpha to neutron
                        if np.isnan(P_ALPHA): P_ALPHA = 0.0
                    # Compute Fraction
                    p_val = P_BREM / (P_ALPHA+P_ECRH+1.0)
                    print(p_val,P_BREM,P_ALPHA,P_ECRH)
                    p_val = p_val * (1.0 + (rand()-0.5)*2.0*noise)
                    control, error, Ival = self.pid_controller(setpoint, p_val, pid_K, pid_Ti, pid_Td, error_old, Ival, dt)
                    self.particle_sources[species]['PID_pradfrac_gaussian']['pid_I'] = Ival
                    self.particle_sources[species]['PID_pradfrac_gaussian']['previous_error'] = error
                    # Adjust U 
                    control = max(control,0)
                    control = min(control,N_IN)
                    # Compute integrand
                    integrand = np.exp(-(rho_grid-rho_0)**2/sigma_rho**2) * self.dVdr(rho_grid)
                    integrand = integrand.flatten()
                    # 
                    cte = control / np.trapezoid(integrand,self.r_grid)
                    #
                    aux_source = cte * np.exp(-(rho_grid-rho_0)**2/sigma_rho**2)
                    
            # bookeeping
            self.explicit_particle_sources[species][source_type][it,:] = aux_source

    def pid_controller(self, setpoint, pv, kp, tau_i, tau_d, previous_error, integral, dt):
        """
        The PID control algorithm for feedback control.
        """
        error = setpoint - pv
        integral += error * dt
        derivative = (error - previous_error) / dt
        control = kp * (error + integral/tau_i + tau_d * derivative)
        return control, error, integral
            
    def compute_diffusive_heat_flux(self,it):
        """
        Computes the coefficients Dp (diffusion) and cp (convection) of
        each species assuming a simple diffusive mode, which are then used to 
        fill the LHS_pressure matrix. The coefficients are computed such that 
        the heat flux is Q = -Dp*dp/dr + cp*p. This allows for the pressure to 
        be evolved in time according to LoDestro's method
        """
        from scipy.interpolate import CubicSpline
        
        r_grid = self.r_grid
        
        for species in self.list_of_species:
            
            chi = self.heat_fluxes_info['diffusive']['chi']
            convective_fact = self.heat_fluxes_info['diffusive']['convective_fact']
            
            # p_r = CubicSpline(r_grid,self.P[species][it,:])
            # dpdr = p_r.derivative()
            # dpdr = dpdr(r_grid)
            p_r = self.P[species][it,:]
            dpdr = akima_derivative(r_grid,p_r)
            
            # n_r = CubicSpline(r_grid,self.N[species][it,:])
            # dndr = n_r.derivative()
            n_r = self.N[species][it,:]
            dndr = akima_derivative(r_grid,n_r)
            
            self.Dp[species][it,:] = chi
            
            c = (chi/n_r)*dndr + convective_fact*self.Gamma_turb[species][it,:]/n_r
            c[0] = 0.0
            
            self.cp[species][it,:] = c
            
            # this is for bookeeping
            self.Q_turb[species][it,:] = -chi * dpdr + p_r*( (chi/n_r)*dndr + convective_fact*self.Gamma_turb[species][it,:]/n_r)
            
    def compute_diffusive_particle_flux(self,it):
        """
        Computes the coefficient Dn (diffusion) which is used to fill the LHS_density matrices. 
        The coefficient is computed such that the particle flux is
        Gamma = -Dn*dn/dr. Here we assume that advection is zero and that Dn is a constant
        """
        from scipy.interpolate import CubicSpline
        
        r_grid = self.r_grid
        
        for species in self.list_of_species:
            
            Dn = self.particle_fluxes_info['diffusive']['Dn']
            
            # n_r = CubicSpline(r_grid,self.N[species][it,:])
            n_r = self.N[species][it,:]
            # dndr = n_r.derivative()
            # dndr = dndr(r_grid)
            dndr = akima_derivative(r_grid,n_r)
            
            self.Dn[species][it,:] = Dn
            
            self.cn[species][it,:] = 0.0
                        
            # this is used in heat flux (Q=-n\chi*dT/dr + convective_fact*T*Gamma_turb)
            self.Gamma_turb[species][it,:] = -Dn * dndr
            
    def compute_diffusive_advective_particle_flux(self,it):
        """
        Computes the coefficient Dn (diffusion) and cn (convection) which are used to fill the LHS_density matrices. 
        The coefficient is computed such that the particle flux is
        Gamma = -Dn*dn/dr + cn*n. Dn and cn can be constants or functions of rho
        """
        from scipy.interpolate import CubicSpline
        
        r_grid = self.r_grid
        
        Dn = self.particle_fluxes_info['diffusive_advective']['Dn']
        cn = self.particle_fluxes_info['diffusive_advective']['cn']
        
        if callable(Dn):
            Dn = Dn(self.rho_grid)
        elif isinstance(Dn, (float, int)):
            pass
        else:
            raise ValueError('ERROR: Dn can only be a function of rho or an integer/float!')
        
        if callable(cn):
            cn = cn(self.rho_grid)
        elif isinstance(cn, (float, int)):
            pass
        else:
            raise ValueError('ERROR: cn can only be a function of rho or an integer/float!')
        
        for species in self.list_of_species:
            
            # n_r = CubicSpline(r_grid,self.N[species][it,:])
            n_r = self.N[species][it,:]
            # dndr = n_r.derivative()
            # dndr = dndr(r_grid)
            dndr = akima_derivative(r_grid,n_r)
            
            self.Dn[species][it,:] = Dn
            self.cn[species][it,:] = cn
            
            # Force cn(rho=0.0) to be 0.0
            self.cn[species][it,0] = 0.0
                        
            # this is used in heat flux (Q=-n\chi*dT/dr + convective_fact*T*Gamma_turb)
            self.Gamma_turb[species][it,:] = -Dn * dndr + cn * n_r
            
    def compute_normalized_Dn_cn_rho_aLn_dependent_particle_flux(self,it):
        """
        Computes the coefficients Dn (diffusion) and cn (convection) of
        each species. Here Dn and cn are the normalized diffusion coefficient and advection velocity
        which are functions of rho and a/Ln. The denormalization is:
        
        Dn = Dn_normalized * Gamma_gB * a/n_j
        cn = cn_normalized * Gamma_gB / n_j
        
        where Gamma_gB = 2.0*sqrt(2)*n_ref*sqrt(m_ref)*(EC*T_ref)**1.5 / (EC*Bref*aminor)**2
        
        The particle flux is then:
        Gamma = -Dn*dn/dr + cn*n. 
        This allows for the density to be evolved in time according to LoDestro's method
        """
        
        r_grid = self.r_grid
        
        Dn_normalized_func = self.particle_fluxes_info['normalized_Dn_cn_rho_aLn_dependent']['Dn'] # This is a 2D function of (rho,aLn)
        cn_normalized_func = self.particle_fluxes_info['normalized_Dn_cn_rho_aLn_dependent']['cn'] # This is a 2D function of (rho,aLn)
        mref = self.particle_fluxes_info['normalized_Dn_cn_rho_aLn_dependent']['mass_ref_species']
        
        for species in self.list_of_species:
            
            n_r = self.N[species][it,:]
            dndr = polyfit_derivative_fast(r_grid,n_r,deg=12)
            # dndr = akima_derivative(r_grid,n_r)
            
            a_Ln = - self.aminor * dndr / n_r
            
            Dn_normalized = Dn_normalized_func(self.rho_grid, a_Ln)
            cn_normalized = cn_normalized_func(self.rho_grid, a_Ln)  
            
            nref = n_r
            Tref = self.T[species][it,:]
            
            Gamma_gB = 2.0*np.sqrt(2)*nref*np.sqrt(mref)*(EC*Tref)**1.5 / (EC*self.Bref*self.aminor)**2       
            
            self.Dn[species][it,:] = Dn_normalized * Gamma_gB * self.aminor / n_r
            self.cn[species][it,:] = cn_normalized * Gamma_gB / n_r
            # Force cn(rho=0.0) to be 0.0
            self.cn[species][it,0] = 0.0
                        
            # this is used in heat flux (Q=-n\chi*dT/dr + convective_fact*T*Gamma_turb)
            self.Gamma_turb[species][it,:] = -self.Dn[species][it,:] * dndr + self.cn[species][it,:] * n_r
            
    def compute_beurskens_heat_flux(self,it):
        """
        Computes the coefficients Dp (diffusion) and cp (convection) of
        each species according to Beurskens model
        """
        from scipy.interpolate import CubicSpline
        
        chi_base = self.heat_fluxes_info['beurskens']['chi_base']
        chi_electrons = self.heat_fluxes_info['beurskens']['chi_electrons']
        aLT_critical = self.heat_fluxes_info['beurskens']['aLT_critical']
        alpha = self.heat_fluxes_info['beurskens']['alpha']
        stiffness = self.heat_fluxes_info['beurskens']['stiffness']
        convective_fact = self.heat_fluxes_info['beurskens']['convective_fact']
        m_ref_species = self.heat_fluxes_info['beurskens']['mass_ref_species']
        
        if callable(stiffness) and callable(aLT_critical):
            stiffness = stiffness(self.rho_grid)
            aLT_critical = aLT_critical(self.rho_grid)
        elif isinstance(stiffness, (float, int)) and isinstance(aLT_critical, (float, int)):
            pass
        else:
            raise ValueError('ERROR: stiffnes and aLTcritical can only be a function or integer/float!')
        
        chi = {}
        r_grid = self.r_grid
        
        ## electrons
        chi['electrons'] = chi_electrons * np.ones(self.Nr)
        T_electrons = self.T['electrons'][it,:]
         
        ## IONS
        for ion in self.plasma.ion_species:
            T_ion = self.T[ion][it,:]
            n_ion = self.N[ion][it,:]
            
            # T_polyfit = np.poly1d( np.polyfit(r_grid,T_ion,deg=12) )
            # dTdr_polyfit = np.poly1d( T_polyfit.deriv() )
            # dTdr_polyfit = dTdr_polyfit(r_grid)
            # dTdr = dTdr_polyfit
            dTdr = polyfit_derivative_fast(r_grid,T_ion,deg=12)
            
            a_LT = self.aminor * dTdr / T_ion
            
            a_LT_filtered = -a_LT
            
            X = a_LT_filtered - aLT_critical
            
            chi_turb = stiffness * X * np.heaviside(X,1) * (T_electrons/T_ion)**alpha
            
            # gyro-Bohm heat flux
            Tref = T_ion         # self.T[ref_species][it,:]
            mref = m_ref_species # self.plasma.mass[ref_species]
            # nref = self.N[ref_species][it,:] # we cannot use this one, as this poses issue to the impurities
            nref = n_ion
            
            Q_gB = 2.0*np.sqrt(2)*nref*np.sqrt(mref)*(EC*Tref)**2.5 / (EC*self.Bref*self.aminor)**2
            
            chi_gB = Q_gB * self.aminor/(n_ion*EC*T_ion)
            
            chi_turb = chi_gB * chi_turb
            
            chi[ion] = chi_base + chi_turb
        
        for species in self.list_of_species:
            
            # p_r = CubicSpline(r_grid,self.P[species][it,:])
            # dpdr = p_r.derivative()
            p_r = self.P[species][it,:]
            # dpdr = dpdr(r_grid)
            dpdr = akima_derivative(r_grid,p_r)
            # n_r = CubicSpline(r_grid,self.N[species][it,:])
            # dndr = n_r.derivative()
            
            n_r = self.N[species][it,:]
            # dndr = dndr(r_grid)
            dndr = akima_derivative(r_grid,n_r)
            
            D = chi[species]
            
            # save D of ALL subiter
            self.Dp_keep[species][it].append(np.array(D))
            
            # average to smooth-out eventual oscillations
            D_avg = np.mean(np.array(self.Dp_keep[species][it]), axis=0)

            self.Dp[species][it,:] = D_avg
            
            c = (chi[species]/n_r)*dndr + convective_fact*self.Gamma_turb[species][it,:]/n_r
            c[0] = 0.0
            
            # should we also average 'c' ??
            
            self.cp[species][it,:] = c
            
            # this is for bookeeping
            self.Q_turb[species][it,:] = -chi[species] * dpdr + p_r*( (chi[species]/n_r)*dndr + convective_fact*self.Gamma_turb[species][it,:]/n_r)
    
    def initialize_NEO(self, surfaces_k, DKES_coeffs_file, dt_NEO=None, Er_root_type='ion_root', dt_Er_ambipolar=None, n_workers_NEO=1):
        """
        n_workers_NEO controls how call_NEO evaluates the DKES surfaces (each
        surface's ambipolar root + neoclassical transport coefficients is an
        independent PENTA calculation -- see call_PENTA_surface in
        PENTA/Sources/penta_interface_mod.f90). n_workers_NEO=1 (default)
        evaluates them serially in this process. n_workers_NEO>1 spreads them
        across a persistent pool of that many worker processes instead
        (created once, here, and reused for every call_NEO call for the
        lifetime of this solver -- NOT recreated per time step).

        Note this must be OS processes, not threads: PENTA keeps its working
        state in Fortran module-level (SAVE) variables in libpenta.so, which
        are not safe to share across concurrent calls within a single process.
        """

        self.DKES_nuv, self.DKES_Erv, self.DKES_D11, self.DKES_D31, self.DKES_D33, self.dkes_k, self.roa_dkes_k = self.process_DKES_file(DKES_coeffs_file,surfaces_k)

        self.Er_root_type = Er_root_type
        self.dt_Er_ambipolar = dt_Er_ambipolar

        self.n_workers_NEO = n_workers_NEO
        if(n_workers_NEO is not None and n_workers_NEO > 1):
            self.neo_pool = ProcessPoolExecutor(max_workers=n_workers_NEO, initializer=_init_NEO_worker)
            print(f'Using {n_workers_NEO} worker processes for NEO (PENTA) surface calculations')
        else:
            self.neo_pool = None
            
        self.dt_NEO = dt_NEO

    # --- Original serial call_NEO, kept here (disabled) for reference/rollback. ---
    # def call_NEO(self):
    #
    #     start = perf_counter()
    #     t = self.time[self.it]
    #
    #     # Search for a new ambipolar Er root at t=0, and (if dt_Er_ambipolar is
    #     # set) every dt_Er_ambipolar of simulation time thereafter (measured
    #     # from t=0, not from tstart). In between, the Er root from the last
    #     # search is reused. A tolerance of dt/2 is used since t=0 (or a
    #     # multiple of dt_Er_ambipolar) may not land exactly on the time grid.
    #     at_t_zero = np.isclose(t, 0.0, atol=self.dt)
    #     if(self.dt_Er_ambipolar is None and self.subiter==1):
    #         look_for_ambipolar = at_t_zero
    #     else:
    #         # remainder = t % self.dt_Er_ambipolar
    #         # at_multiple = np.isclose(remainder, 0.0, atol=self.dt/2) or \
    #         #               np.isclose(remainder, self.dt_Er_ambipolar, atol=self.dt/2)
    #         # look_for_ambipolar = at_t_zero or at_multiple
    #         look_for_ambipolar = False
    #
    #     # Process kinetic profiles data
    #     ne, dnedrho = akima_interp(self.rho_grid, self.N['electrons'][self.it,:], self.roa_dkes_k)
    #     te, dtedrho = akima_interp(self.rho_grid, self.T['electrons'][self.it,:], self.roa_dkes_k)
    #
    #     ns_dkes = self.roa_dkes_k.size
    #     nion_prof = len(self.plasma.ion_species)
    #     ni      = np.empty((ns_dkes, nion_prof))
    #     dnidrho = np.empty((ns_dkes, nion_prof))
    #     ti      = np.empty((ns_dkes, nion_prof))
    #     dtidrho = np.empty((ns_dkes, nion_prof))
    #     for j, ion in enumerate(self.plasma.ion_species):
    #         ni[:,j], dnidrho[:,j] = akima_interp(self.rho_grid, self.N[ion][self.it,:], self.roa_dkes_k)
    #         ti[:,j], dtidrho[:,j] = akima_interp(self.rho_grid, self.T[ion][self.it,:], self.roa_dkes_k)
    #
    #     # Inputs
    #     Matom_prof = [self.plasma.mass[ion]    for ion in self.plasma.ion_species]
    #     Zatom_prof = [self.plasma.Zcharge[ion] for ion in self.plasma.ion_species]
    #     EparB = 0.0
    #     Er_min_Vcm = -250
    #     Er_max_Vcm = 250
    #     bsq = self.Bsq_spline(self.roa_dkes_k)
    #
    #     # Placeholders -- their values don't really matter
    #     btheta = np.zeros_like(self.roa_dkes_k)
    #     bzeta = np.zeros_like(self.roa_dkes_k)
    #     iota = np.zeros_like(self.roa_dkes_k)
    #     phip = np.zeros_like(self.roa_dkes_k)
    #     chip = np.zeros_like(self.roa_dkes_k)
    #     vp = np.zeros_like(self.roa_dkes_k)
    #
    #     if(look_for_ambipolar): self.Er = np.zeros_like(self.roa_dkes_k)
    #
    #     end = perf_counter()
    #     self.time_prepare += end-start
    #
    #     start = perf_counter()
    #
    #     output_Er, output_Dn, output_cn, output_Dp, output_cp = self.libPenta.call_PENTA(
    #         Matom_prof, Zatom_prof,
    #         ne, dnedrho, te, dtedrho, ni, dnidrho, ti, dtidrho,
    #         self.aminor, self.Rmajor, vp, chip, phip, iota, btheta, bzeta, bsq,
    #         self.dkes_k, self.roa_dkes_k,
    #         self.DKES_nuv, self.DKES_Erv, self.DKES_D11, self.DKES_D31, self.DKES_D33,
    #         Er_min_Vcm, Er_max_Vcm, EparB, self.Er, self.Er_root_type, look_for_ambipolar,
    #         self.rho_grid)
    #
    #     self.Er = output_Er
    #
    #     end = perf_counter()
    #     self.time_call_NEO += end-start

    def call_NEO(self,it):
        
        t = self.time[it]

        # Search for a new ambipolar Er root at t=0 and every dt_Er_ambipolar
        at_t_zero = np.isclose(t, 0.0, atol=self.dt/2)
        if(at_t_zero):
            look_for_ambipolar = True
        elif(self.dt_Er_ambipolar is None):
            look_for_ambipolar = False
        else:
            nsteps_per_Er = round(self.dt_Er_ambipolar/self.dt)
            look_for_ambipolar = (nsteps_per_Er % it == 0)

        # Process kinetic profiles data
        ne, dnedrho = akima_interp(self.rho_grid, self.N['electrons'][it,:], self.roa_dkes_k)
        te, dtedrho = akima_interp(self.rho_grid, self.T['electrons'][it,:], self.roa_dkes_k)

        ns_dkes = self.roa_dkes_k.size
        nion_prof = len(self.plasma.ion_species)
        ni      = np.empty((ns_dkes, nion_prof))
        dnidrho = np.empty((ns_dkes, nion_prof))
        ti      = np.empty((ns_dkes, nion_prof))
        dtidrho = np.empty((ns_dkes, nion_prof))
        for j, ion in enumerate(self.plasma.ion_species):
            ni[:,j], dnidrho[:,j] = akima_interp(self.rho_grid, self.N[ion][it,:], self.roa_dkes_k)
            ti[:,j], dtidrho[:,j] = akima_interp(self.rho_grid, self.T[ion][it,:], self.roa_dkes_k)

        # Inputs
        Matom_prof = [self.plasma.mass[ion]    for ion in self.plasma.ion_species]
        Zatom_prof = [self.plasma.Zcharge[ion] for ion in self.plasma.ion_species]
        EparB = np.zeros_like(self.roa_dkes_k)
        Er_min_Vcm = -250
        Er_max_Vcm = 250
        bsq = self.Bsq_spline(self.roa_dkes_k)
        Er_k, _ = akima_interp(self.rho_grid, self.Er[it,:], self.roa_dkes_k)

        # Placeholders -- their values don't really matter
        btheta = np.zeros_like(self.roa_dkes_k)
        bzeta = np.zeros_like(self.roa_dkes_k)
        iota = np.zeros_like(self.roa_dkes_k)
        phip = np.zeros_like(self.roa_dkes_k)
        chip = np.zeros_like(self.roa_dkes_k)
        vp = np.zeros_like(self.roa_dkes_k)

        if(self.neo_pool is None):
            # Serial: call_PENTA loops over all ns_dkes surfaces and
            # interpolates onto self.rho_grid in one shot.
            self.Er[it,:], Dn_out, cn_out, Dp_out, cp_out = self.libPenta.call_PENTA(
                Matom_prof, Zatom_prof,
                ne, dnedrho, te, dtedrho, ni, dnidrho, ti, dtidrho,
                self.aminor, self.Rmajor, vp, chip, phip, iota, btheta, bzeta, bsq,
                self.dkes_k, self.roa_dkes_k,
                self.DKES_nuv, self.DKES_Erv, self.DKES_D11, self.DKES_D31, self.DKES_D33,
                Er_min_Vcm, Er_max_Vcm, EparB, Er_k, self.Er_root_type, look_for_ambipolar,
                self.rho_grid)
        else:
            # Parallel: spread the ns_dkes independent surface calculations
            # across the persistent process pool created in initialize_NEO.
            tasks = [
                (Matom_prof, Zatom_prof,
                 ne[k], dnedrho[k], te[k], dtedrho[k], ni[k,:], dnidrho[k,:], ti[k,:], dtidrho[k,:],
                 self.aminor, self.Rmajor, vp[k], chip[k], phip[k], iota[k], btheta[k], bzeta[k], bsq[k],
                 int(self.dkes_k[k]), self.roa_dkes_k[k],
                 self.DKES_nuv, self.DKES_Erv, self.DKES_D11[k,:,:], self.DKES_D31[k,:,:], self.DKES_D33[k,:,:],
                 Er_min_Vcm, Er_max_Vcm, EparB[k], Er_k[k], self.Er_root_type, look_for_ambipolar)
                for k in range(ns_dkes)
                    ]
            results = list(self.neo_pool.map(_call_PENTA_surface_worker, tasks))
            Er_PENTA = np.array([r[1] for r in results])
            Dn_PENTA = np.array([r[2] for r in results]).T
            cn_PENTA = np.array([r[3] for r in results]).T
            Dp_PENTA = np.array([r[4] for r in results]).T
            cp_PENTA = np.array([r[5] for r in results]).T

            self.Er[it,:], Dn_out, cn_out, Dp_out, cp_out = self.libPenta.call_PENTA_interpolate(
                self.roa_dkes_k, Er_PENTA, Dn_PENTA, cn_PENTA, Dp_PENTA, cp_PENTA, self.rho_grid)

        # Assign transport coefficients
        for isp, species in enumerate(self.list_of_species):
            self.Dn_NEO[species][it,:] = Dn_out[isp,:]
            self.cn_NEO[species][it,:] = cn_out[isp,:]
            self.Dp_NEO[species][it,:] = Dp_out[isp,:]
            self.cp_NEO[species][it,:] = cp_out[isp,:]
    
    def process_DKES_file(self,DKES_coeffs_file,surfaces_k):
        """
        Reads a DKES coefficients file with the header
        dkes_k  Er_v           nu_v           D11            D31            D33
        Inner loop: nu_v
        Middle loop: Er_v
        Outer loop: dkes_k

        surfaces_k is a list/array of the dkes_k surfaces to extract.

        The file's first line has the format 'ns_surfaces xxx', where xxx is
        the total number of surfaces, and the second line is the column
        header (dkes_k Er_v nu_v D11 D31 D33).

        Returns nu_v, Er_v, DKES_D11, DKES_D31, DKES_D33, surfaces_k, rho_k,
        where nu_v and Er_v are the unique (non-repeated) 1D arrays of
        length Nu and Ne, DKES_D11, DKES_D31, DKES_D33 are 3D arrays of
        shape (len(surfaces_k), Nu, Ne), indexed as [surface, nu, Er], one
        slice per requested surface (in the order given in surfaces_k), and
        rho_k = sqrt((surfaces_k-1)/(ns_surfaces-1)) is the normalized
        radial coordinate of each requested surface.
        """
        with open(DKES_coeffs_file) as f:
            ns_surfaces = int(f.readline().split()[1])

        data = np.loadtxt(DKES_coeffs_file, skiprows=2)

        dkes_k_col = data[:,0].astype(int)
        Er_col     = data[:,1]
        nu_col     = data[:,2]
        D11_col    = data[:,3]
        D31_col    = data[:,4]
        D33_col    = data[:,5]

        dkes_k_v = np.unique(dkes_k_col)
        Er_v     = np.unique(Er_col)
        nu_v     = np.unique(nu_col)

        Nk = len(dkes_k_v)
        Ne = len(Er_v)
        Nu = len(nu_v)

        if len(data) != Nk*Ne*Nu:
            raise ValueError(f'ERROR: {DKES_coeffs_file} does not have a regular grid of dkes_k, Er_v and nu_v (found {len(data)} rows, expected {Nk}*{Ne}*{Nu}={Nk*Ne*Nu})')

        dkes_k_grid = dkes_k_col.reshape(Nk,Ne,Nu)
        Er_grid     = Er_col.reshape(Nk,Ne,Nu)
        nu_grid     = nu_col.reshape(Nk,Ne,Nu)

        if not np.array_equal(nu_grid, np.broadcast_to(nu_v,(Nk,Ne,Nu))):
            raise ValueError(f'ERROR: {DKES_coeffs_file} does not have nu_v as the inner loop!')
        if not np.array_equal(Er_grid, np.broadcast_to(Er_v[None,:,None],(Nk,Ne,Nu))):
            raise ValueError(f'ERROR: {DKES_coeffs_file} does not have Er_v as the middle loop!')
        if not np.array_equal(dkes_k_grid, np.broadcast_to(dkes_k_v[:,None,None],(Nk,Ne,Nu))):
            raise ValueError(f'ERROR: {DKES_coeffs_file} does not have dkes_k (first column) as the outer loop!')

        missing_k = [k for k in surfaces_k if k not in dkes_k_v]
        if len(missing_k) > 0:
            raise ValueError(f'ERROR: surfaces_k {missing_k} not found in {DKES_coeffs_file}!')

        D11_grid = D11_col.reshape(Nk,Ne,Nu)
        D31_grid = D31_col.reshape(Nk,Ne,Nu)
        D33_grid = D33_col.reshape(Nk,Ne,Nu)

        DKES_D11 = []
        DKES_D31 = []
        DKES_D33 = []
        for k in surfaces_k:
            idx = np.nonzero(dkes_k_v == k)[0][0]
            # grids are (Ne,Nu) for this surface; transpose to (Nu,Ne)
            DKES_D11.append(D11_grid[idx].T)
            DKES_D31.append(D31_grid[idx].T)
            DKES_D33.append(D33_grid[idx].T)

        DKES_D11 = np.stack(DKES_D11, axis=0)
        DKES_D31 = np.stack(DKES_D31, axis=0)
        DKES_D33 = np.stack(DKES_D33, axis=0)

        surfaces_k = np.asarray(surfaces_k)
        rho_k = np.sqrt((surfaces_k - 1) / (ns_surfaces - 1))

        return nu_v, Er_v, DKES_D11, DKES_D31, DKES_D33, surfaces_k, rho_k
    
    def set_NEO_coefficients_from_previous(self,it):
        """ Sets NEO transport coefficients at current it equal to previous it """
        
        for species in self.list_of_species:
            self.Dn_NEO[species][it,:] = self.Dn_NEO[species][it-1,:]
            self.cn_NEO[species][it,:] = self.cn_NEO[species][it-1,:]
            self.Dp_NEO[species][it,:] = self.Dp_NEO[species][it-1,:]
            self.cp_NEO[species][it,:] = self.cp_NEO[species][it-1,:]
        
    def add_NEO_transport_coefficients(self,it):
        """ Adds NEO transport coefficients to Dp,cp,Dn,cn """
        
        for species in self.list_of_species:
            self.Dp[species][it,:] += self.Dp_NEO[species][it,:]
            self.cp[species][it,:] += self.cp_NEO[species][it,:]
            self.Dn[species][it,:] += self.Dn_NEO[species][it,:]
            self.cn[species][it,:] += self.cn_NEO[species][it,:]

    # def compute_NEO_particle_flux(self,it):
        
    #     from scipy.interpolate import CubicSpline, Akima1DInterpolator
        
    #     root = 'ion_root'

    #     PENTA_class = PENTA(folder_path='.', plasma=self.plasma, lverb=False)
        
    #     for sp,species in enumerate(self.list_of_species):
            
    #         # Gamma = np.array( PENTA_class.Gamma_Maxw[species] )
    #         Gamma = np.array( PENTA_class.Gamma[species,root] )
    #         roa_PENTA = PENTA_class.roa[root]
            
    #         ## Include Gamma(r=0) = 0
    #         rho_extended = np.concatenate([[0.0],roa_PENTA])
    #         Gamma_extended = np.concatenate(([0.0],Gamma))
    #         # Gamma_interp = CubicSpline(rho_extended,Gamma_extended,extrapolate=True,bc_type='natural')
    #         Gamma_interp = Akima1DInterpolator(rho_extended,Gamma_extended,method='makima')
    #         Gamma_interp.extrapolate = True
    #         # This is used when computing the heat flux
    #         self.Gamma_NEO[species][it,:] = Gamma_interp(self.rho_grid)
            
    #         # Compute Dn and cn
    #         PENTA_class.set_plasma_solver_transport_coeffs()
    #         Dn = PENTA_class.Dn[species,root][:,sp] # the sp index picks the self diffusion coeff, Dn_aa
    #         cn = PENTA_class.cn[species,root][:]
            
    #         # extended Dn and cn towards the axis by setting them to 0.0
    #         Dn_extended = np.concatenate(([0.0],Dn))
    #         cn_extended = np.concatenate(([0.0],cn))
    #         roa_extended = np.concatenate(([0.0],roa_PENTA))
            
    #         Dn_extended_spline = Akima1DInterpolator(roa_extended,Dn_extended,method='makima')
    #         Dn_extended_spline.extrapolate = True
    #         #
    #         cn_extended_spline = Akima1DInterpolator(roa_extended,cn_extended,method='makima')
    #         cn_extended_spline.extrapolate = True
            
    #         self.Dn[species][it,:] = Dn_extended_spline(self.rho_grid)
    #         self.cn[species][it,:] = cn_extended_spline(self.rho_grid)
            
    # def compute_NEO_heat_flux(self,it):
        
    #     from scipy.interpolate import CubicSpline, Akima1DInterpolator
        
    #     root = 'ion_root'

    #     PENTA_class = PENTA(folder_path='.', plasma=self.plasma, lverb=False)
        
    #     for sp,species in enumerate(self.list_of_species):
            
    #         QoT = np.array( PENTA_class.QoT[species,root] )
    #         roa_PENTA = PENTA_class.roa[root]
            
    #         T_PENTA = CubicSpline(self.rho_grid, self.T[species][it,:])
    #         T_PENTA = T_PENTA(roa_PENTA)
            
    #         Q = QoT * EC * T_PENTA
            
    #         ## Include Q(r=0) = 0
    #         rho_extended = np.concatenate([[0.0],roa_PENTA])
    #         Q_extended = np.concatenate(([0.0],Q))
            
    #         Q_interp = Akima1DInterpolator(rho_extended,Q_extended,method='makima')
    #         Q_interp.extrapolate = True
            
    #         # This is used when computing the heat flux
    #         self.Q_NEO[species][it,:] = Q_interp(self.rho_grid)
            
    #         # Compute Dp and cp
    #         PENTA_class.set_plasma_solver_transport_coeffs()
    #         Dp = PENTA_class.Dp[species,root][:,sp] # the sp index picks the self diffusion coeff, Dn_aa
    #         cp = PENTA_class.cp[species,root][:]
            
    #         # extended Dp and cp towards the axis by setting them to 0.0
    #         Dp_extended = np.concatenate(([0.0],Dp))
    #         cp_extended = np.concatenate(([0.0],cp))
    #         roa_extended = np.concatenate(([0.0],roa_PENTA))
            
    #         Dp_extended_spline = Akima1DInterpolator(roa_extended,Dp_extended,method='makima')
    #         Dp_extended_spline.extrapolate = True
    #         #
    #         cp_extended_spline = Akima1DInterpolator(roa_extended,cp_extended,method='makima')
    #         cp_extended_spline.extrapolate = True
            
    #         self.Dp[species][it,:] = Dp_extended_spline(self.rho_grid)
    #         self.cp[species][it,:] = cp_extended_spline(self.rho_grid)
            
    # def compute_NEO_plus_beurskens_heat_flux(self,it):
        
    #     from scipy.interpolate import CubicSpline, Akima1DInterpolator
        
        
    #     ##############################################################################################
    #     ################################ NEO contribution ############################################
    #     ##############################################################################################
        
    #     root = 'ion_root'

    #     PENTA_class = PENTA(folder_path='.', plasma=self.plasma, lverb=False)
        
    #     for sp,species in enumerate(self.list_of_species):
            
    #         QoT = np.array( PENTA_class.QoT[species,root] )
    #         roa_PENTA = PENTA_class.roa[root]
            
    #         T_PENTA = CubicSpline(self.rho_grid, self.T[species][it,:])
    #         T_PENTA = T_PENTA(roa_PENTA)
            
    #         Q = QoT * EC * T_PENTA
            
    #         ## Include Q(r=0) = 0
    #         rho_extended = np.concatenate([[0.0],roa_PENTA])
    #         Q_extended = np.concatenate(([0.0],Q))
            
    #         Q_interp = Akima1DInterpolator(rho_extended,Q_extended,method='makima')
    #         Q_interp.extrapolate = True
            
    #         # This is used when computing the heat flux
    #         self.Q_NEO[species][it,:] = Q_interp(self.rho_grid)
            
    #         # Compute Dp and cp
    #         PENTA_class.set_plasma_solver_transport_coeffs()
    #         Dp = PENTA_class.Dp[species,root][:,sp] # the sp index picks the self diffusion coeff, Dn_aa
    #         cp = PENTA_class.cp[species,root][:]
            
    #         # extended Dp and cp towards the axis by setting them to 0.0
    #         Dp_extended = np.concatenate(([0.0],Dp))
    #         cp_extended = np.concatenate(([0.0],cp))
    #         roa_extended = np.concatenate(([0.0],roa_PENTA))
            
    #         Dp_extended_spline = Akima1DInterpolator(roa_extended,Dp_extended,method='makima')
    #         Dp_extended_spline.extrapolate = True
    #         #
    #         cp_extended_spline = Akima1DInterpolator(roa_extended,cp_extended,method='makima')
    #         cp_extended_spline.extrapolate = True
            
    #         self.Dp[species][it,:] = Dp_extended_spline(self.rho_grid)
    #         self.cp[species][it,:] = cp_extended_spline(self.rho_grid)
            
    #     ##############################################################################################
    #     ########################## Beurskens contribution ############################################
    #     ##############################################################################################
        
    #     chi_electrons = self.heat_fluxes_info['dkespenta_beurskens']['chi_electrons']
    #     aLT_critical = self.heat_fluxes_info['dkespenta_beurskens']['aLT_critical']
    #     alpha = self.heat_fluxes_info['dkespenta_beurskens']['alpha']
    #     stiffness = self.heat_fluxes_info['dkespenta_beurskens']['stiffness']
    #     convective_fact = self.heat_fluxes_info['dkespenta_beurskens']['convective_fact']
        
    #     chi = {}
        
    #     ## electrons
    #     chi['electrons'] = chi_electrons * np.ones(self.Nr)
        
    #     r_grid = self.r_grid
    #     Bsq = self.Bsq(self.rho_grid)
        
    #     ## IONS
    #     for ion in self.plasma.ion_species:
    #         T_ion = self.T[ion][it,:]
    #         T_electrons = self.T['electrons'][it,:]
            
    #         T_r = CubicSpline(r_grid,T_ion)
    #         # dTdr_non_filtered = T_r.derivative()
            
    #         T_polyfit = np.poly1d( np.polyfit(r_grid,T_ion,deg=12) )
    #         dTdr_polyfit = np.poly1d( T_polyfit.deriv() )
    #         dTdr_polyfit = dTdr_polyfit(r_grid)
            
    #         dTdr = dTdr_polyfit  
    #         # dTdr = dTdr_non_filtered(r_grid)
            
    #         a_LT = self.aminor * dTdr / T_ion
            
    #         a_LT_filtered = -a_LT
            
    #         X = a_LT_filtered - aLT_critical
            
    #         chi_turb = stiffness * X * np.heaviside(X,1) * (T_electrons/T_ion)**alpha
            
    #         mi = self.plasma.mass[ion]
    #         qi = self.plasma.charge[ion]
    
    #         chi_gB = (EC*T_ion/mi)**1.5 * mi*mi / (qi**2 * Bsq) / self.aminor
            
    #         chi_turb = chi_gB * chi_turb
            
    #         chi[ion] = chi_turb
        
    #     for species in self.list_of_species:
            
    #         p_r = CubicSpline(r_grid,self.P[species][it,:])
    #         dpdr = p_r.derivative()
            
    #         n_r = CubicSpline(r_grid,self.N[species][it,:])
    #         dndr = n_r.derivative()
            
    #         n_r = self.N[species][it,:]
    #         dndr = dndr(r_grid)
            
    #         D = chi[species]
            
    #         # save D of ALL subiter
    #         self.Dp_keep[species][it].append(np.array(D))
            
    #         # average to smooth-out eventual oscillations
    #         D_avg = np.mean(np.array(self.Dp_keep[species][it]), axis=0)
    #         D = D_avg

    #         # add Beurskens contribution
    #         self.Dp[species][it,:] += D
            
    #         c = (chi[species]/n_r)*dndr + convective_fact*self.Gamma_turb[species][it,:]/n_r
    #         c[0] = 0.0

    #         # add Beurskens contribution
    #         self.cp[species][it,:] += c
            
    #         # this is for bookeeping
    #         self.Q_turb[species][it,:] = -chi[species] * dpdr(r_grid) + p_r(r_grid)*( (chi[species]/n_r)*dndr + convective_fact*self.Gamma_turb[species][it,:]/n_r)
            
    def solve_density_equations(self,it):
        """Sets LHS matrices and RHS vectors of density equations and solves them"""
        from libstell.fusion import FUSION
        fusion = FUSION()
        
        for species in self.list_of_species:      
            if(species=='tritium' and self.constrain_nT):
                continue
            if(species=='electrons'):
                continue
            
            RHS_vector = self.N[species][it-1,:] + self.dt*self.get_explicit_particle_sources(species,it)
            # apply BC
            RHS_vector[-1] = self.edge_density_BC[species]
            
            # get LHS matrix
            LHS_matrix = self.get_LHS_density(species,it)
            
            # solve system
            self.N[species][it,:] = self.solve_sparse_system(LHS_matrix,RHS_vector)
            
        if(self.constrain_nT):
            self.N['tritium'][it,:] = self.N['deuterium'][it,:]   
            
        if(self.solve_fast_alphas):
            nD = self.N['deuterium'][it,:]
            nT = self.N['tritium'][it,:]
            TD = self.T['deuterium'][it,:]
            TT = self.T['tritium'][it,:]
            sigmav = fusion.sigmaBH(0.5*(TD+TT),'DT')
            #
            self.N['alphas_fast'][it,:] = (self.N['alphas_fast'][it-1,:] + self.dt*nD*nT*sigmav) / (1+self.dt/self.tau_fast_alphas)
        
        # update electron density from quasi neutrality
        self.N['electrons'][it,:] = 0.0
        for ion in self.plasma.ion_species:
            self.N['electrons'][it,:] += self.N[ion][it,:] * self.plasma.Zcharge[ion]        
        if(self.solve_fast_alphas):
            self.N['electrons'][it,:] += 2*self.N['alphas_fast'][it,:]  
        
        dens = []
        for species in self.list_of_species:
            dens.append(self.N[species][it,:])
            
        return np.concatenate(dens)
    
    def solve_pressure_equations(self,it):
        """
        Sets LHS matrices and RHS vectors of pressure equation and solves it. 
        The pressure equations of all species are in one single system because coll heat
        exchange couples all temperatures due to its implicit implementation
        """
        
        RHS_vector = []
        for species in self.list_of_species:
        
            g = self.P[species][it-1,:] + (2./3)*self.dt*self.get_explicit_energy_sources(species,it)
            # apply edge Dirichlet BC
            g[-1] = self.edge_pressure_BC[species]
            
            RHS_vector.append(g)
            
        RHS_vector = np.concatenate(RHS_vector)
        
        LHS_matrix = self.get_LHS_pressure(it)
        
        press = self.solve_sparse_system(LHS_matrix,RHS_vector)
        
        # update self.P and self.T
        Nr = self.Nr
        k=0
        for species in self.list_of_species:
            
            self.P[species][it,:] = press[k:(k+Nr)]
            self.T[species][it,:] = self.P[species][it,:] / (EC*self.N[species][it,:])
            
            k = k+Nr
            
        return press
          
    def get_explicit_particle_sources(self,species,it):
        """Returns 1D array with explicit particle sources"""
        explicit_source = 0.0
        for source_type in self.explicit_particle_sources[species]:
            explicit_source += self.explicit_particle_sources[species][source_type][it,:]
            
        return explicit_source
    
    def get_explicit_energy_sources(self,species,it):
        """Returns 1D array with explicit energy sources"""
        explicit_source = 0.0
        for source_type in self.explicit_energy_sources[species]:
            explicit_source += self.explicit_energy_sources[species][source_type][it,:]
            
        return explicit_source
    
    def get_LHS_density(self,species,it):
        """Updates LHS density sparse matrix and returns it"""
        from scipy.sparse import diags
        
        dr = self.dr
        Vp = self.dVdr
        Nr = self.Nr
        dt = self.dt
        
        vp = Vp(self.rho_grid)
        vp_inner = vp[1:-1]
        
        Dn = self.Dn[species][it,:]
        cn = self.cn[species][it,:]
        
        ############################################
        ############### COMPUTE LHS ################
        ############################################
        lower = np.zeros(Nr-1)
        main = np.zeros(Nr)
        upper = np.zeros(Nr-1)
        
        ## 0<r<a (inner grid, no boundary points)
        Dn_plus = (Dn[2:]+Dn[1:-1]) / 2
        Dn_minus = (Dn[0:-2]+Dn[1:-1]) / 2
        
        Vp_plus = (vp[2:]+vp[1:-1]) / 2
        Vp_minus = (vp[0:-2]+vp[1:-1]) / 2
        
        VDplus  = Vp_plus*Dn_plus / (vp_inner*dr**2)
        VDminus = Vp_minus*Dn_minus / (vp_inner*dr**2)
        
        cplus = cn[2:]*vp[2:] / (2*vp_inner*dr)
        cminus = cn[0:-2]*vp[0:-2] / (2*vp_inner*dr)

        main[1:-1] = 1.0 + dt*(VDplus + VDminus)
        upper[1:] = dt*(-VDplus + cplus)
        lower[0:-1] = dt*(-VDminus - cminus)
        
        ## r=0
        main[0] = 1.0 + dt*( 4*Dn[0]/dr**2 + 2*cn[1]/dr )
        upper[0] = -4*dt*Dn[0]/dr**2
                
        ## r=a
        main[-1] = 1.0
        lower[-1] = 0.0
        
        ####
        # LHS = diags([lower, main, upper], offsets=[-1, 0, 1], format="csr")    
        self.LHS_density.data[self.density_main_diag] = main
        self.LHS_density.data[self.density_upper_diag] = upper
        self.LHS_density.data[self.density_lower_diag] = lower
        
        return self.LHS_density
    
    def get_LHS_pressure(self,it):
        """Updates LHS pressure sparse matrix and returns it"""

        drho = self.drho
        dr = self.aminor * drho
        Vp = self.dVdr
        Nr = self.Nr
        num_species = len(self.list_of_species)
        dt_fact = (2./3.)*self.dt
        
        vp = Vp(self.rho_grid)
        vp_inner = vp[1:-1]
        
        Vp_plus = (vp[2:]+vp[1:-1]) / 2
        Vp_minus = (vp[0:-2]+vp[1:-1]) / 2
        
        lower = np.zeros(self.Nr-1)
        main = np.zeros(self.Nr)
        upper = np.zeros(self.Nr-1)
        
        for ispecies,species in enumerate(self.list_of_species):
            
            Dp = self.Dp[species][it,:]    
            cp = self.cp[species][it,:]
            
            ############################################
            ############### COMPUTE LHS ################
            ############################################
            
            ## 0<r<a (inner grid, no boundary points)
            Dp_plus = (Dp[2:]+Dp[1:-1]) / 2
            Dp_minus = (Dp[0:-2]+Dp[1:-1]) / 2
            
            VDplus  = Vp_plus*Dp_plus / (vp_inner*dr**2)
            VDminus = Vp_minus*Dp_minus / (vp_inner*dr**2)
            
            cplus = cp[2:]*vp[2:] / (2*vp_inner*dr)
            cminus = cp[0:-2]*vp[0:-2] / (2*vp_inner*dr)

            main[1:-1] = 1.0 + dt_fact*(VDplus + VDminus)
            upper[1:] = dt_fact*(-VDplus + cplus)
            lower[0:-1] = dt_fact*(-VDminus - cminus)
            
            ## r=0
            main[0] = 1.0 + dt_fact*( 4*Dp[0]/dr**2 + 2*cp[1]/dr )
            upper[0] = -4*dt_fact*Dp[0]/dr**2
            
            # DIFF_list.append( diags([lower, main, upper], offsets=[-1, 0, 1], format="csr") )
            self.LHS_pressure.data[self.pressure_main_block[ispecies]] = main
            self.LHS_pressure.data[self.pressure_lower_block[ispecies]] = lower
            self.LHS_pressure.data[self.pressure_upper_block[ispecies]] = upper
            
            
        LHS = self.LHS_pressure
        # DIFF_list = [DIFF[species] for species in self.list_of_species]
              
        # Add implicit terms from sources       
        if hasattr(self, "solve_coll_heat_exchange") and self.solve_coll_heat_exchange:
            LHS = LHS - dt_fact*self.get_collisionalHeatExchange(it)
        
        # impose Dirichlet boundary condition
        # indices of to seto to zero (end of each species block)
        rows_to_fix = np.arange(Nr - 1, num_species * Nr, Nr)
        # Efficiently zero out rows
        for row in rows_to_fix:
            start = LHS.indptr[row]
            end = LHS.indptr[row + 1]
            LHS.data[start:end] = 0.0  # zero existing entries
        # Then set diagonal elements to 1
        LHS[rows_to_fix, rows_to_fix] = 1.0
               
        return LHS
    
    def solve_sparse_system(self,matrix,vect):
        """Solves the linear system matrix.X=vect, where matrix is sparse"""
        from scipy.sparse.linalg import spsolve
        sol = spsolve(matrix,vect)
        
        return sol
    
    def get_collisionalHeatExchange(self,it):
        """
        Computes collisional heat exchange between all species and 
        constructs the coll heat exchange sparse matrix which is added 
        to LHS_pressure
        """
        from libstell.collisions import COLLISIONS
        from scipy import sparse
        
        coll = COLLISIONS()
        
        num_species = len(self.list_of_species)
        
        gamma = np.zeros((num_species,num_species,self.Nr))
        
        for is1,species1 in enumerate(self.list_of_species):

            m1 = self.plasma.mass[species1]
            Z1 = self.plasma.Zcharge[species1]
            n1 = self.N[species1][it,:]
            T1 = self.T[species1][it,:]
            
            for is2,species2 in enumerate(self.list_of_species[is1:], start=is1):
                
                m2 = self.plasma.mass[species2]
                Z2 = self.plasma.Zcharge[species2]
                n2 = self.N[species2][it,:]
                T2 = self.T[species2][it,:]
                
                # get Coulomb logarithm
                if(Z1>0 and Z2>0):
                    clog = coll.coullog_ii(m1,Z1,n1,T1,m2,Z2,n2,T2)
                elif(Z1>0 and Z2<0):
                    clog = coll.coullog_ei(n2,T2,m1,Z1,n1,T1)
                elif(Z1<0 and Z2>0):
                    clog = coll.coullog_ei(n1,T1,m2,Z2,n2,T2)
                else:
                    clog = 0.0

                const = (8/np.sqrt(np.pi))*(Z1*Z2*EC*EC)**2 * clog / (8*np.pi*EPS0**2)

                vth_s1_sqr = 2*EC*T1/m1
                vth_s2_sqr = 2*EC*T2/m2
                
                den = m1 * m2 * (vth_s1_sqr + vth_s2_sqr)**1.5
                
                gamma[is1,is2,:] = const / den
                # fill symmetric entry
                gamma[is2,is1,:] = gamma[is1,is2,:]
        
        N_arr = np.stack([self.N[s][it,:] for s in self.list_of_species])  # shape (Ns, Nr)
        W_s1_s2 = gamma * N_arr[:, np.newaxis, :]
        aux_B   = gamma * N_arr[np.newaxis, :, :]
        # W_s1_s2[is1,is2,:] = gamma[is1,is2,:]*n1      
        # aux_B[is1,is2,:] = gamma[is1,is2,:]*n2

        # Add aux_B matrix
        W_s1_s2[np.arange(num_species), np.arange(num_species), :] -= np.sum(aux_B, axis=1)

        # W_out matrix
        
        # W_out = np.zeros((num_species*self.Nr,num_species*self.Nr))
        # j=0
        # for is1 in range(num_species):
        #     for ir1 in range(self.Nr):
        #         p=0
        #         for is2 in range(num_species):
        #             for ir2 in range(self.Nr):
        #                 if(ir1==ir2):
        #                     W_out[j,p] = W_s1_s2[is1,is2,ir2]  
        #                 p=p+1
        #         j = j+1
        # W_out = sparse.csr_matrix(W_out)
        
        Nr = self.Nr
        is1, is2, ir = np.meshgrid(
            np.arange(num_species),
            np.arange(num_species),
            np.arange(Nr),
            indexing="ij"
        )

        # Flatten
        is1 = is1.ravel()
        is2 = is2.ravel()
        ir = ir.ravel()

        # Map to global indices
        rows = is1 * Nr + ir
        cols = is2 * Nr + ir
        data = W_s1_s2[is1, is2, ir]

        W_out = sparse.csr_matrix((data, (rows, cols)), shape=(num_species * Nr, num_species * Nr))
        
        return W_out
    
    # def call_PENTA3(self,it):
    #     import subprocess
    #     from concurrent.futures import ProcessPoolExecutor, as_completed
    #     import functools
        
    #     # create PLASMA class in order to write PENTA inputs       
    #     plasma_PENTA = PLASMA(self.list_of_species)
    #     for species in self.list_of_species:
    #         plasma_PENTA.set_density(species,'interp',rho_vals=self.rho_grid,n_vals=self.N[species][it,:])
    #         plasma_PENTA.set_temperature(species,'interp',rho_vals=self.rho_grid,T_vals=self.T[species][it,:])
        
    #     plasma_profiles_extension = 'transp_solver'
    #     plasma_PENTA.write_plasma_profiles_to_PENTA3(filename='plasma_profiles_'+plasma_profiles_extension+'.dat')
    #     plasma_PENTA.write_PENTA_namelist()

    #     try:
    #         surfaces = self.particle_fluxes_info['dkespenta']['surfaces']
    #     except:
    #         try:
    #             surfaces = self.heat_fluxes_info['dkespenta']['surfaces']
    #         except:
    #             try:
    #                 surfaces = self.particle_fluxes_info['dkespenta_beurskens']['surfaces']
    #             except:
    #                 surfaces = self.heat_fluxes_info['dkespenta_beurskens']['surfaces']
        
    #     time_sec = []
    #     with ProcessPoolExecutor() as executor:
    #         futures = [executor.submit(process_surfaces, surface, self.wout_path) for surface in surfaces]

    #         for future in as_completed(futures):
    #             elapsed_seconds = future.result()
    #             time_sec.append(elapsed_seconds)
    #             # print(f'Surface processed in {elapsed_seconds:.2f} seconds')
            
    #     # delete files not needed
    #     remove = 'rm ucontra* sigmas* flows_vs_Er*'
    #     subprocess.run(remove, shell=True, check=True, text=True, capture_output=True)
        
    #     # merge _surface_# files into single file
    #     merge_and_delete('fluxes_vs_roa_surface*','fluxes_vs_roa')
    #     merge_and_delete('fluxes_vs_Er_surface*','fluxes_vs_Er')
    #     merge_and_delete('flows_vs_roa_surface*','flows_vs_roa')
    #     merge_and_delete('Jprl_vs_roa_surface*','Jprl_vs_roa')
    #     merge_and_delete('particleTransportCoeffs_vs_roa_surface*','particleTransportCoeffs_vs_roa')
    #     merge_and_delete('heatTransportCoeffs_vs_roa_surface*','heatTransportCoeffs_vs_roa')
    #     merge_and_delete('plasma_profiles_check_surface*','plasma_profiles_check')
        
    def call_save_output(self,output_filename,dt_save):
        """Saves simulation in output joblib file"""
        from types import SimpleNamespace
        from pathlib import Path
        import joblib
        
        # check if extension of output_filename is .joblib; if not, add
        output_filename = str(Path(output_filename).with_suffix(".joblib"))
        
        # save the class (cannot save solver directly cause it contains lambda functions...)
        saved_class = SimpleNamespace()
        saved_class.rho_grid = self.rho_grid
        saved_class.r_grid = self.r_grid
        saved_class.dVdr = self.dVdr(self.rho_grid)
        saved_class.aminor = self.aminor
        saved_class.Rmajor = self.Rmajor
        saved_class.Baxis = self.Baxis
        saved_class.Bref  = self.Bref
        saved_class.list_of_species = self.list_of_species
        # In case of a VMEC equilibrium
        if hasattr(self,'iota23'):
            saved_class.iota2o3 = self.iota23
        
        # only save at minimum every dt=dt_save
        freq = max(1, round(dt_save / self.dt))
        sl = slice(0, None, freq)  # defines the slice once

        saved_class.time = self.time[sl]
        saved_class.Nt = len(self.time[sl])

        for attr in ('N','T','Dn','cn','Dp','cp','Q_NEO','Q_turb','Gamma_NEO','Gamma_turb'):
            setattr(saved_class, attr, {})
            for species in self.list_of_species:
                getattr(saved_class, attr)[species] = getattr(self, attr)[species][sl, :]
                
        if(self.add_NEO):
            for attr in ('Dn_NEO','cn_NEO','Dp_NEO','cp_NEO'):
                setattr(saved_class, attr, {})
                for species in self.list_of_species:
                    getattr(saved_class, attr)[species] = getattr(self, attr)[species][sl, :]
        
        if 'alphas_fast' in self.N:
            saved_class.N['alphas_fast'] = self.N['alphas_fast'][sl, :]
        
        # nested dict attributes
        nested_attrs = ['explicit_energy_sources','explicit_particle_sources']
        for species in self.list_of_species:
            for attr in nested_attrs:
                saved_class.__dict__.setdefault(attr, {})
                saved_class.__dict__[attr].setdefault(species, {})
                for type_string, arr in getattr(self, attr)[species].items():
                    saved_class.__dict__[attr][species][type_string] = arr[sl, :]
                    
            # Save particle PID's info (if they exist). This is useful for restarts
            pid_keys = [k for k in self.particle_sources[species] if k.startswith("PID_")]
            if pid_keys:
                if not hasattr(saved_class, "particle_sources"):
                    saved_class.particle_sources = {}
                if species not in saved_class.particle_sources:
                    saved_class.particle_sources[species] = {}
            for key in pid_keys:
                saved_class.particle_sources[species][key] = {
                    "previous_error": self.particle_sources[species][key]["previous_error"],
                    "pid_I": self.particle_sources[species][key]["pid_I"]}
                
            # Save energy PID's info (if they exist). This is useful for restarts
            pid_keys = [k for k in self.energy_sources[species] if k.startswith("PID_")]
            if pid_keys:
                if not hasattr(saved_class, "energy_sources"):
                    saved_class.energy_sources = {}
                if species not in saved_class.energy_sources:
                    saved_class.energy_sources[species] = {}
            for key in pid_keys:
                saved_class.energy_sources[species][key] = {
                    "previous_error": self.energy_sources[species][key]["previous_error"],
                    "pid_I": self.energy_sources[species][key]["pid_I"]}
                
        saved_class.Er = self.Er

        joblib.dump(saved_class, output_filename)
        
def merge_output_files(*output_files,concatenated_file=None):
    """ Merges sequential joblib output files into a single one. 
    Returns concatenated class and only saves concatenated joblib file if concatenated_file is not None"""
    from types import SimpleNamespace
    import joblib
    import warnings
    from copy import deepcopy
    
    ATTR_TIME_DEP = ('N','T','Dn','cn','Dp','cp','Q_NEO','Q_turb','Gamma_NEO','Gamma_turb')
    
    SOURCE_ATTRS = ('explicit_energy_sources', 'explicit_particle_sources')

    GRID_ATTRS = ('rho_grid', 'r_grid', 'dVdr')

    SCALAR_ATTRS = ('aminor', 'Rmajor', 'Baxis', 'Bref')
    
    OPTIONAL_ATTRS = ('iota23')
    
    ########################## AUX FUNCT ##################################
    def check_same(name, ref, val):
        if isinstance(ref, np.ndarray):
            if not np.allclose(ref, val):
                warnings.warn(f"Invariant mismatch in {name}")
        else:
            if ref != val:
                warnings.warn(f"Invariant mismatch in {name}")

   
    ################ CHECK TIME IS SEQUENTIAL ##################################
    time_all = []
    for file in output_files:
        solver = joblib.load(file)
        time_all.append(solver.time)
    #
    time_all = np.concatenate(time_all)
    is_sequential = np.all(time_all[1:] >= time_all[:-1])
    if(not is_sequential):
        raise ValueError('ERROR: Time is not sequential in the given files...')
    
    ######################## CONCATENATE DATA ##################################
    # ------------------------------------------------------------
    # Load first solver as reference
    # ------------------------------------------------------------
    ref_solver = joblib.load(output_files[0])
    concatenated_class = deepcopy(ref_solver)
    # ------------------------------------------------------------
    # Loop over remaining solvers
    # ------------------------------------------------------------
    for file in output_files[1:]:
        solver = joblib.load(file)

        # -------------------------------
        # (1) Check invariant attributes
        # -------------------------------
        for attr in GRID_ATTRS:
            check_same(attr, getattr(ref_solver, attr), getattr(solver, attr))

        for attr in SCALAR_ATTRS:
            check_same(attr, getattr(ref_solver, attr), getattr(solver, attr))

        if solver.list_of_species != ref_solver.list_of_species:
            warnings.warn("list_of_species mismatch")

        # -------------------------------
        # (2) Concatenate time
        # -------------------------------
        t_old = concatenated_class.time
        t_new = solver.time

        if np.isclose(t_old[-1], t_new[0]):
            time_slice = slice(1, None)
        else:
            time_slice = slice(None)

        concatenated_class.time = np.concatenate(
            [t_old, t_new[time_slice]]
        )

        # -------------------------------
        # (3) Concatenate time-dependent attributes
        # -------------------------------
        for attr in ATTR_TIME_DEP:
            ref_attr = getattr(concatenated_class, attr)
            new_attr = getattr(solver, attr)

            # species in list_of_species
            for species in ref_solver.list_of_species:
                if species in new_attr:
                    ref_attr[species] = np.concatenate(
                    [ref_attr[species],
                     new_attr[species][time_slice, :]],
                    axis=0
                )
                else:
                    warnings.warn(
                        f"{attr}: species '{species}' missing in one solver"
                    )

            # special species: alphas_fast
            if attr == 'N' and 'alphas_fast' in new_attr:
                if 'alphas_fast' not in ref_attr:
                    ref_attr['alphas_fast'] = new_attr['alphas_fast'][time_slice, :]
                else:
                    ref_attr['alphas_fast'] = np.concatenate(
                        [ref_attr['alphas_fast'],
                        new_attr['alphas_fast'][time_slice, :]],
                        axis=0
                    )
                    
        # ------------------------------------------------------------
        # (4) Concatenate explicit source terms
        #     Structure: sources[species][key][time, space]
        # ------------------------------------------------------------
        for attr in SOURCE_ATTRS:

            ref_sources = getattr(concatenated_class, attr)
            new_sources = getattr(solver, attr)

            # Loop over species in the reference solver
            for species in ref_sources.keys():

                if species not in new_sources:
                    warnings.warn(
                        f"{attr}: species '{species}' missing in one solver"
                    )
                    continue

                for key in ref_sources[species]:

                    if key not in new_sources[species]:
                        warnings.warn(
                            f"{attr}[{species}]: key '{key}' missing in one solver"
                        )
                        continue

                    ref_sources[species][key] = np.concatenate(
                        [
                            ref_sources[species][key],
                            new_sources[species][key][time_slice, :]
                        ],
                        axis=0
                    )
        
    concatenated_class.Nt = len(concatenated_class.time)
    
    if(concatenated_file is not None):
        joblib.dump(concatenated_class, concatenated_file)
    
    return concatenated_class
                    
        
# def process_surfaces(surface,wout_path):
#     import time
#     import subprocess
    
#     start_time = time.time()
    
#     type_of_write = 0
#     Er_min_V_cm = -100
#     Er_max_V_cm = 200

#     # wout_path = solver_class.wout_path
#     EparB = 0.0

#     #Sonine (Laguerre) polynomials
#     Smax = 1
    
#     plasma_profiles_extension = 'transp_solver'
    
#     extension_output_files = f'_surface_{surface}'
        
#     extension_star_files = f'surface_{surface}'

#     call_penta3 = f'~/bin/xpenta {extension_star_files} {Er_min_V_cm} {Er_max_V_cm} {surface} {type_of_write} {wout_path} {plasma_profiles_extension} {EparB} {Smax} {extension_output_files}'
    
#     result = subprocess.run(call_penta3, shell=True, check=True, text=True, capture_output=True)
#     # print(result.stdout)
#     if(result.stderr):
#         print(result.stderr)
        
#     end_time = time.time()
#     elapsed_time = (end_time - start_time)
#     return elapsed_time
                
# def merge_and_delete(pattern, output_filename):
#     """
#     Merges files matching the given pattern into a single file and deletes the originals.
    
#     Parameters:
#     pattern (str): The pattern to match files (e.g., 'fluxes_vs_roa_surface_*').
#     output_filename (str): The name of the output file.
#     """
#     import os
#     import re
#     import glob
    
#     def extract_number(filename):
#         match = re.search(r'_(\d+)$', filename)  # Extract number at the end
#         return int(match.group(1)) if match else float('inf')

#     # Find and sort matching files
#     file_list = glob.glob(pattern)
#     file_list.sort(key=extract_number)

#     if not file_list:
#         print(f"No files found matching pattern: {pattern}")
#         return

#     header_written = False
#     with open(output_filename, 'w') as outfile:
#         for filename in file_list:
#             with open(filename, 'r') as infile:
#                 lines = infile.readlines()
#                 if not header_written:
#                     outfile.write(lines[0])  # Write header
#                     outfile.write(lines[1])
#                     header_written = True
#                 outfile.writelines(lines[2:])  # Write data

#     # Delete original files
#     for filename in file_list:
#         os.remove(filename)
        
def initialize_LHS_density(Nr):
    """
    Initializes LHS density sparse matrix and returns it. This avoids having to define
    a sparse matrix at each iteration. The indexes of the sparse matrix to be
    updated are also returned
    """
    from scipy.sparse import diags
    lower = np.ones(Nr-1)
    main  = np.ones(Nr)
    upper = np.ones(Nr-1)

    # Build identical dummy blocks
    A = diags([lower, main, upper], offsets=[-1, 0, 1], format="csr")
    
    rows, cols = A.nonzero()
    mask_main  = (rows == cols)
    mask_lower = (rows == cols + 1)
    mask_upper = (rows + 1 == cols)
    
    return A,mask_main,mask_lower,mask_upper

def initialize_LHS_pressure(Nr, num_species):
    """
    Initializes LHS pressure sparse matrix and returns it. This avoids having to define
    a sparse matrix at each iteration. The indexes of the sparse matrix to be
    updated are also returned
    """
    from scipy.sparse import diags, block_diag
    # Use ones so SciPy allocates all expected data entries
    lower = np.ones(Nr-1)
    main  = np.ones(Nr)
    upper = np.ones(Nr-1)

    # Build identical dummy blocks
    blocks = [diags([lower, main, upper], offsets=[-1, 0, 1], format="csr")
              for _ in range(num_species)]

    # Combine them into one big block-diagonal sparse matrix
    A = block_diag(blocks, format="csr")
    block_size = Nr
    
    nblocks = A.shape[0] // block_size
    lower_indices = []
    main_indices  = []
    upper_indices = []

    # Precompute row pointers for CSR
    for b in range(nblocks):
        row_start = b * block_size
        row_end   = (b+1) * block_size
        rows = np.arange(row_start, row_end)
        for r in rows:
            start = A.indptr[r]
            end   = A.indptr[r+1]
            cols = A.indices[start:end]
            for i, c in zip(range(start, end), cols):
                if c == r:
                    main_indices.append(i)
                elif c == r-1:
                    lower_indices.append(i)
                elif c == r+1:
                    upper_indices.append(i)

    # Split indices by block
    def split_blocks(indices):
        per_block = []
        counts = [0] + [block_size if block_size <= len(indices) else len(indices) for _ in range(nblocks)]
        idx = 0
        for b in range(nblocks):
            block_idx = []
            while idx < len(indices) and A.indices[indices[idx]] < (b+1)*block_size:
                block_idx.append(indices[idx])
                idx += 1
            per_block.append(np.array(block_idx, dtype=int))
        return per_block

    lower_blocks = split_blocks(lower_indices)
    main_blocks  = split_blocks(main_indices)
    upper_blocks = split_blocks(upper_indices)

    return A, lower_blocks, main_blocks, upper_blocks

@njit
def akima_derivative(x, y):
    """
    Compute Akima spline derivatives (Hermite form) at points x,
    reproducing Fortran r8akherm1(ipx=0) behavior.

    Parameters
    ----------
    x : 1D array of shape (N,)
        Strictly increasing coordinate values
    y : 1D array of shape (N,)
        Function values at x

    Returns
    -------
    dy : 1D array of shape (N,)
        Akima numerical derivatives at x
    """
    n = x.size
    dy = np.zeros_like(y)
    if n < 2:
        raise ValueError("Need at least 2 points")

    # First divided differences
    m = np.empty(n - 1)
    for i in range(n - 1):
        m[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i])

    if n == 2:
        dy[0] = m[0]
        dy[1] = m[0]
        return dy

    # --- Boundary slopes (consider the case ipx=0, as in the EZsplines fortran code) ---
    cxp = m[0]
    cxpp = m[1]
    cxm = m[-1]
    cxmm = m[-2]

    dy[0] = 1.5 * cxp - 0.5 * cxpp
    dy[-1] = 1.5 * cxm - 0.5 * cxmm

    # Ghost slopes for extrapolation
    cxtrap0 = 2.0 * dy[0] - cxp
    cxtrap1 = 2.0 * dy[-1] - cxm

    # --- Interior points ---
    for i in range(1, n - 1):
        # Left slopes
        if i == 1:
            cxmm = cxtrap0
        else:
            cxmm = (y[i - 1] - y[i - 2]) / (x[i - 1] - x[i - 2])

        cxm = (y[i] - y[i - 1]) / (x[i] - x[i - 1])
        cxp = (y[i + 1] - y[i]) / (x[i + 1] - x[i])

        # Right slopes
        if i == n - 2:
            cxpp = cxtrap1
        else:
            cxpp = (y[i + 2] - y[i + 1]) / (x[i + 2] - x[i + 1])

        # Akima weights
        w1 = abs(cxp - cxpp)
        w2 = abs(cxm - cxmm)

        if (w1 + w2) == 0.0:
            dy[i] = 0.5 * (cxm + cxp)
        else:
            dy[i] = (w1 * cxm + w2 * cxp) / (w1 + w2)

    return dy

@njit
def akima_interp(x, y, xnew):
    """
    Akima interpolation of y(x) and its derivative, evaluated at arbitrary
    points xnew (unlike akima_derivative, which only returns the derivative
    at the original nodes x).

    Computes the Akima node slopes via akima_derivative, then evaluates the
    piecewise cubic Hermite interpolant (and its analytic derivative) built
    from those slopes -- equivalent to PSPLINE's r8herm1ev on top of
    r8akherm1(ipx=0) slopes. Points outside [x[0], x[-1]] are extrapolated
    using the boundary cubic segment.

    Parameters
    ----------
    x : 1D array of shape (N,)
        Strictly increasing source coordinate values
    y : 1D array of shape (N,)
        Function values at x
    xnew : 1D array of shape (M,)
        Points at which to evaluate the interpolant (need not be sorted)

    Returns
    -------
    ynew, dynew : 1D arrays of shape (M,)
        Interpolated values and derivatives at xnew
    """
    dy = akima_derivative(x, y)
    n = x.size
    m = xnew.size
    ynew = np.empty(m)
    dynew = np.empty(m)
    for k in range(m):
        xk = xnew[k]
        i = np.searchsorted(x, xk) - 1
        if i < 0:
            i = 0
        elif i > n - 2:
            i = n - 2

        h = x[i + 1] - x[i]
        t = (xk - x[i]) / h
        t2 = t * t
        t3 = t2 * t

        h00 = 2.0 * t3 - 3.0 * t2 + 1.0
        h10 = t3 - 2.0 * t2 + t
        h01 = -2.0 * t3 + 3.0 * t2
        h11 = t3 - t2
        ynew[k] = h00 * y[i] + h10 * h * dy[i] + h01 * y[i + 1] + h11 * h * dy[i + 1]

        dh00 = 6.0 * t2 - 6.0 * t
        dh10 = 3.0 * t2 - 4.0 * t + 1.0
        dh01 = -6.0 * t2 + 6.0 * t
        dh11 = 3.0 * t2 - 2.0 * t
        dynew[k] = (dh00 * y[i] + dh10 * h * dy[i] + dh01 * y[i + 1] + dh11 * h * dy[i + 1]) / h

    return ynew, dynew

@njit
def polyfit_derivative_fast(x, y, deg):
    """
    Makes a polyfit of degree deg to the points (xi,yi), computes its
    derivative at xi and returns it. This function is compiled JIT and
    is a much faster alternative to np.polyfit(...).derivative()
    when called many times (for instance inside a loop)
    """
    N = len(x)
    V = np.zeros((N, deg + 1))
    for i in range(N):
        p = 1.0
        for j in range(deg, -1, -1):
            V[i, j] = p
            p *= x[i]

    # Use least-squares solution (like np.linalg.lstsq)
    # Numba doesn't support np.linalg.lstsq, but we can emulate it via SVD
    U, s, VT = np.linalg.svd(V, full_matrices=False)
    c = np.zeros(deg + 1)
    for i in range(len(s)):
        c += (U[:, i] @ y) / s[i] * VT[i, :]

    # Derivative coefficients (decreasing powers)
    dcoeff = np.zeros(deg)
    for i in range(deg):
        dcoeff[i] = (deg - i) * c[i]

    # Evaluate derivative (Horner)
    dy = np.zeros(N)
    for k in range(N):
        val = 0.0
        for i in range(deg):
            val = val * x[k] + dcoeff[i]
        dy[k] = val

    return dy  
            
# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)      