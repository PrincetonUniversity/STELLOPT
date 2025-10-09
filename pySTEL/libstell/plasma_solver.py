"""
This library provides a python class for solving density 
and pressure transport equations
"""

import numpy as np
import sys
from time import perf_counter

from libstell.plasma import PLASMA
from libstell.penta import PENTA

# Constants
EC = 1.602176634E-19 # Electron charge [C]
EPS0 = 8.8541878188E-12 # Vacuum permittivity [F/m]

class PLASMA_SOLVER:
    
    def __init__(self, list_of_species, solve_fast_alphas=False, tau_fast_alphas=0.5, constrain_nT=False):
        
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
        
    def read_restart_file(self,restart_filepath):
        # reads a restart .joblib file and sets initial profiles & BC's according to last itertion in file
        import joblib
        from pathlib import Path
        
        # Convert to Path object
        file_path = Path(restart_filepath)

        # Check the extension
        if file_path.suffix != ".joblib":
            raise ValueError(f"Error: The file '{restart_filepath}' does not have a .joblib extension.")
        
        restart_solver = joblib.load(restart_filepath)
        
        # check list_of_species in file are the same as those in this run
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
            if 'alphas_fast' not in restart_solver.N or 'alphas_thermal' not in restart_solver.N:
                raise ValueError('restart file does not have alphas density! Yet you want to solve with alphas...')
            else:
                self.alphas_fast_density_restart    = restart_solver.N['alphas_fast'][-1,:]
              
    def set_equilibrium(self,type: str,wout_path=None,aminor=None,Rmajor=None,B=None):
        
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
                    
                    self.B0 = np.sqrt(np.squeeze(vmec_out.bdotb)[0])   
                    self.Bsq = CubicSpline(roa,np.squeeze(vmec_out.bdotb))
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
                    self.B0 = B
                    self.Bsq = lambda rho: B*B
                    
    def set_energy_source(self,species,source_type, total_power=None, sigma_rho=None, rho_0=None, fraction_alpha_heating=None, cte_source=None, time_dependent_factor=None, lambda_function_2D=None):
        # electrons: 'Bremsstrahlung', 'Coll_Heat_Exchange', 'Er', 'external', 'alpha_heating'
        # ions: 'Coll_Heat_Exchange', 'Er', 'external', 'alpha_heating'
        # 'constant' is for benchmarking

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
                
    def set_particle_source(self,species,source_type, injected_particles_per_sec=None, rho_0=None, sigma_rho=None, cte_source=None, time_dependent_factor=None, lambda_function_2D=None):
        
        import inspect
        
        # check species exist in list_of_species
        if species not in self.list_of_species:
            raise ValueError(f"ERROR: Species {species} is not in the plasma.")
         
        match source_type:
            case 'external_gaussian':
                if((injected_particles_per_sec is None) or (sigma_rho is None) or (rho_0 is None)):
                    print('ERROR: Need to provide Smax [part/(m^3*s)], rho_0 and sigma_rho for gaussian external source')
                    exit(1) 
                else:
                    self.particle_sources[species][source_type] = {'injected_particles_per_sec' : injected_particles_per_sec, 'rho_0' : rho_0, 'sigma_rho' : sigma_rho}
            #
            case 'time_dependent_gaussian':
                if((injected_particles_per_sec is None) or (rho_0 is None) or (sigma_rho is None) or (time_dependent_factor is None)):
                    print('ERROR: Need to provide Smax [par/(m^3*s)], rho_0, sigma_rho and a time depenedent factof for time-dependent gaussian')
                    exit(1) 
                else:
                    self.particle_sources[species][source_type] = {'injected_particles_per_sec' : injected_particles_per_sec, 'rho_0' : rho_0, 'sigma_rho' : sigma_rho, 'time_factor': time_dependent_factor }
            #
            case 'fast_alphas_source':
                # check we are solving fast alphas
                if(not self.solve_fast_alphas):
                    raise ValueError('solve_fast_alphas was set to false, so fast_alphas_source does not make sense...')
                self.particle_sources[species][source_type] = {}
            #
            case 'alpha_particles_sink':
                self.particle_sources[species][source_type] = {}
            #
            case 'constant':
                if(cte_source is None):
                    print('ERROR: cte_source is needed in order to generate a constant source.')
                    exit(0)
                else:
                    self.particle_sources[species][source_type] = {'cte_source' : cte_source}
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
                    self.particle_sources[species][source_type] = {'lambda_function_2D' : lambda_function_2D}
            #
            case _:
                print(f'ERROR: Source type {source_type} is NOT possible')
                exit(0)
                
    def set_heat_fluxes(self, type: str,surfaces=None,chi=None,chi_base=None,aLT_critical=None,alpha=None,stiffness=None,chi_electrons=None,convective_fact=None):
        # sets type of fluxes
        # OPTION1: type='dkespenta'; dkes_folder and surfaces(list of integers) must be provided
        # OPTION2: type='diffusive'; chi must be provided (heat diffusivity; assumes same val for all species)
            
        match type:
            case 'dkespenta':
                #checks that dkes_folder and surfaces are provided
                if( surfaces is None):
                    print('ERROR: surfaces must be provided!!')
                    exit(0)
                    
                self.heat_fluxes_info['type'] = type
                self.heat_fluxes_info[type]['surfaces'] = surfaces
                
            case 'diffusive':
                #checks that diffusion coefficients are provided
                if((chi is None) or (convective_fact is None)):
                    print('ERROR: chi and convective_fact must be provided!')
                    exit(0)
                self.heat_fluxes_info['type'] = type
                self.heat_fluxes_info[type]['chi'] = chi
                self.heat_fluxes_info[type]['convective_fact'] = convective_fact
                
            case 'beurskens':
                if( (chi_base is None) or (aLT_critical is None) or (alpha is None) or (stiffness is None) or (chi_electrons is None) or (convective_fact is None)):
                    print('ERROR: chi_base, aLT_critical, alpha, stiffness, chi_electrons and convective_fact must be given!')
                    exit(0)
                self.heat_fluxes_info['type'] = type
                self.heat_fluxes_info[type]['chi_base'] = chi_base
                self.heat_fluxes_info[type]['aLT_critical'] = aLT_critical
                self.heat_fluxes_info[type]['alpha'] = alpha
                self.heat_fluxes_info[type]['stiffness'] = stiffness
                self.heat_fluxes_info[type]['chi_electrons'] = chi_electrons
                self.heat_fluxes_info[type]['convective_fact'] = convective_fact
                
            case 'dkespenta_beurskens':
                if( (surfaces is None) or (chi_base is None) or (aLT_critical is None) or (alpha is None) or (stiffness is None) or (chi_electrons is None) or (convective_fact is None)):
                    print('ERROR: surfaces, chi_base, aLT_critical, alpha, stiffness, chi_electrons and convective_fact must be given!')
                    exit(0)
                self.heat_fluxes_info['type'] = type
                self.heat_fluxes_info[type]['surfaces'] = surfaces
                self.heat_fluxes_info[type]['chi_base'] = chi_base
                self.heat_fluxes_info[type]['aLT_critical'] = aLT_critical
                self.heat_fluxes_info[type]['alpha'] = alpha
                self.heat_fluxes_info[type]['stiffness'] = stiffness
                self.heat_fluxes_info[type]['chi_electrons'] = chi_electrons
                self.heat_fluxes_info[type]['convective_fact'] = convective_fact
                
    def set_particle_fluxes(self, type: str, surfaces=None,Dn=None):
        # sets type of fluxes
        # OPTION1: type='dkespenta'; dkes_folder and surfaces(list of integers) must be provided
        # OPTION2: type='diffusive'; Dn and chi must be provided (partical and heat collisional diffusion coefficients)
            
        match type:
            case 'dkespenta':
                #checks that dkes_folder and surfaces are provided
                if( surfaces is None ):
                    print('ERROR:surfaces must be provided!!')
                    exit(0)
                    
                self.particle_fluxes_info['type'] = type
                self.particle_fluxes_info[type]['surfaces'] = surfaces
                
            case 'diffusive':
                #checks that diffusion coefficient is provided
                if(Dn is None):
                    print('ERROR: Dn must be provided!')
                    exit(0)
                self.particle_fluxes_info['type'] = type
                self.particle_fluxes_info[type]['Dn'] = Dn
                
    def run(self,Nr,dt,tstart,tend,tolerance=1E-2,max_subiter=12,output_filename=None):
        
        from collections import defaultdict
        
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
        self.print_grid_details()
        
        # initialize self.## variables
        self.initialize_variables()
        
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
            
            ### SUBCYCLE
            delta_p = 10*tolerance
            subiter=1
            while(delta_p > tolerance and subiter<=max_subiter):
                self.subiter = subiter  

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
            self.call_save_output(output_filename)  
            
        end_time = perf_counter()   
        print(f'Plasma Solver took {(end_time-start_time)/60:.2f}min to run.')  
        
    def make_checks(self):
        
        # check equilibrium exists
        if(not hasattr(self,'dVdr')):
            print('ERROR: dVdr MUST BE SET!!')
            exit(1)
        
        for species in self.list_of_species:
            
            # check boundary conditions are set
            if (species not in self.edge_density_BC or species not in self.edge_pressure_BC):
                raise KeyError(f"Missing edge boundary condition for species: {species}")
            
            # check initial profiles are set
            if (species not in self.initial_density or species not in self.initial_pressure):
                raise KeyError(f"Missing initial profile for species: {species}")
            
            # check consistency between boundary conditions and initial profiles
            tol = np.abs(self.edge_density_BC[species]) * np.finfo(float).eps
            if( np.abs(self.initial_density[species](1)-self.edge_density_BC[species]) > 5*tol ):
                raise ValueError(f'Edge density BC not consistent w/ initial density profile')
            tol = np.abs(self.edge_pressure_BC[species]) * np.finfo(float).eps
            if( np.abs(self.initial_pressure[species](1)-self.edge_pressure_BC[species]) > 10*tol ):
                raise ValueError(f'Edge pressure/temperature BC not consistent w/ initial temperature profile')
            
        # check fluxes info is set
        if(not hasattr(self,'heat_fluxes_info')):
            print('ERROR: set_heat_fluxes must be called before running!!')
            exit(1)
        if(not hasattr(self,'particle_fluxes_info')):
            print('ERROR: set_particle_fluxes must be called before running!!')
            exit(1)
            
    def print_grid_details(self):
        
        print(' ')
        print( ' ***********************')
        print(f' *  tstart = {self.tstart:5.2f}s    *')
        print(f' *  tend   = {self.tend:5.2f}s    *')
        print(f' *  dt     = {self.dt:5.2f}s    *')
        print(f' *  Nt     = {self.Nt:3}       *')
        print(f' *  drho   = {self.drho:5.2f}     *')
        print( ' ***********************')
        
        print(' ')
        header_str = '  TIME [s]     NSUB      TE_AXIS [keV]     NE_AXIS [m^-3]    TI_AXIS [keV]    NI_AXIS [m^-3]    MAX(dp/p_old)'  
        print(header_str)
        print('  '+'='*len(header_str))
        
    def initialize_variables(self):
        
        from collections import defaultdict
        
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
        
        Nr = self.Nr
        Nt = self.Nt

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
            
            self.explicit_energy_sources[species] = {}
            for source_type in self.energy_sources[species].keys():
                self.explicit_energy_sources[species][source_type] = np.zeros((Nt,Nr))
            self.explicit_particle_sources[species] = {}
            for source_type in self.particle_sources[species].keys():
                self.explicit_particle_sources[species][source_type] = np.zeros((Nt,Nr))
            
        if(self.solve_fast_alphas):
            self.N['alphas_fast'] = np.zeros((Nt,Nr))
            
    def set_fields_tstart(self):
        
        rho_grid = self.rho_grid
        
        for species in self.list_of_species:
        
            self.P[species][0,:] = self.initial_pressure[species](rho_grid)     
            self.N[species][0,:] = self.initial_density[species](rho_grid)
               
            self.T[species][0,:] = self.P[species][0,:] / (EC*self.N[species][0,:])
        
        # N_alphas are set to ZERO at t=tstart
        # UNLESS read from restart file
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
        
        # set fluxes at t=0
        self.call_fluxes(it=0)  
        
        # return array with (dens,press) at t=tstart
        all_fields = []
        for species in self.list_of_species:
            all_fields.append(self.N[species][0,:])
        for species in self.list_of_species:
            all_fields.append(self.P[species][0,:])
        
        return np.concatenate(all_fields)
                
    def call_fluxes(self,it):
        
        if(self.particle_fluxes_info['type']=='dkespenta' or self.heat_fluxes_info['type']=='dkespenta' or self.heat_fluxes_info['type']=='dkespenta_beurskens'):
            self.call_PENTA3(it)
        
        # compute Dn_interp
        if(self.particle_fluxes_info['type']=='diffusive'):
            self.compute_diffusive_particle_flux(it)
        elif(self.particle_fluxes_info['type']=='dkespenta'):
            self.compute_NEO_particle_flux(it)
        else:
            raise ValueError('ERROR: Not available other type of particle flux...')
        
        # compute Dp_interp and cp_interp
        if(self.heat_fluxes_info['type']=='diffusive'):
            self.compute_diffusive_heat_flux(it)
        elif(self.heat_fluxes_info['type']=='beurskens'):
            self.compute_beurskens_heat_flux(it)
        elif(self.heat_fluxes_info['type']=='dkespenta'):
            self.compute_NEO_heat_flux(it)
        elif(self.heat_fluxes_info['type']=='dkespenta_beurskens'):
            self.compute_NEO_plus_beurskens_heat_flux(it)
        else:
            raise ValueError('ERROR: Not available other type of heat flux...')
        
    def set_explicit_energy_sources(self,species: str, it):
        # returns 1D-array of same size as rho_grid
        # computes sources using info in self.heat_sources[species]
        
        from libstell.fusion import FUSION
        
        fusion = FUSION()
          
        # total_explicit_energy_source = np.zeros(len(rho_grid))
        
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
                    # this source is fully implicit
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
                    raise ValueError('This is still not fully available... need to check its implementation.')
                      
                case 'constant':
                    aux_source = self.energy_sources[species]['constant']['cte_source']

                case _:
                    print(f'ERROR: Source type {source_type} not defined....')
                    exit(0)
                   
            self.explicit_energy_sources[species][source_type][it,:] = aux_source

    def set_explicit_particle_sources(self,species: str, it):
        # returns 1D-array of same size as rho_grid
        # computes sources using info in self.particle_sources[species]
        
        from libstell.fusion import FUSION
        
        fusion = FUSION()
        
        rho_grid = self.rho_grid
          
        # total_explicit_particle_source = np.zeros(len(rho_grid))
        
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
                    
                case 'constant':
                    aux_source = self.particle_sources[species]['constant']['cte_source']
                    
                case 'lambda_2D':
                    lambda_function_2D = self.particle_sources[species][source_type]['lambda_function_2D'] #func(r,t)
                    #
                    aux_source = [lambda_function_2D(r,self.time[it]) for r in self.r_grid]
                    
            # bookeeping
            self.explicit_particle_sources[species][source_type][it,:] = aux_source
            
    def compute_diffusive_heat_flux(self,it):
        # computes an interpolating function for D and c
        
        from scipy.interpolate import CubicSpline
        
        r_grid = self.r_grid
        
        for species in self.list_of_species:
            
            chi = self.heat_fluxes_info['diffusive']['chi']
            convective_fact = self.heat_fluxes_info['diffusive']['convective_fact']
            
            p_r = CubicSpline(r_grid,self.P[species][it,:])
            dpdr = p_r.derivative()
            
            n_r = CubicSpline(r_grid,self.N[species][it,:])
            dndr = n_r.derivative()
            
            n_r = self.N[species][it,:]
            dndr = dndr(r_grid)
            
            self.Dp[species][it,:] = chi
            
            c = (chi/n_r)*dndr + convective_fact*self.Gamma_turb[species][it,:]/n_r
            c[0] = 0.0
            
            self.cp[species][it,:] = c
            
            # this is for bookeeping
            self.Q_turb[species][it,:] = -chi * dpdr(r_grid) + p_r(r_grid)*( (chi/n_r)*dndr + convective_fact*self.Gamma_turb[species][it,:]/n_r)
            
    def compute_diffusive_particle_flux(self,it):
        # computes an interpolating function for Dn
        
        from scipy.interpolate import CubicSpline
        
        r_grid = self.r_grid
        
        for species in self.list_of_species:
            
            Dn = self.particle_fluxes_info['diffusive']['Dn']
            
            n_r = CubicSpline(r_grid,self.N[species][it,:])
            dndr = n_r.derivative()
            
            self.Dn[species][it,:] = Dn
            
            self.cn[species][it,:] = 0.0
                        
            # this is used in heat flux
            self.Gamma_turb[species][it,:] = -Dn * dndr(r_grid)
            
    def compute_beurskens_heat_flux(self,it):
        # uses model in [ref...]
        from scipy.interpolate import CubicSpline
        
        chi_base = self.heat_fluxes_info['beurskens']['chi_base']
        chi_electrons = self.heat_fluxes_info['beurskens']['chi_electrons']
        aLT_critical = self.heat_fluxes_info['beurskens']['aLT_critical']
        alpha = self.heat_fluxes_info['beurskens']['alpha']
        stiffness = self.heat_fluxes_info['beurskens']['stiffness']
        convective_fact = self.heat_fluxes_info['beurskens']['convective_fact']
        
        if callable(stiffness) and callable(aLT_critical):
            stiffness = stiffness(self.rho_grid)
            aLT_critical = aLT_critical(self.rho_grid)
        elif isinstance(stiffness, (float, int)) and isinstance(aLT_critical, (float, int)):
            pass
        else:
            raise ValueError('ERROR: stiffnes and aLTcritical can only be a function or integer/float!')
        
        chi = {}
        
        ## electrons
        chi['electrons'] = chi_electrons * np.ones(self.Nr)
        
        r_grid = self.r_grid
        
        ## IONS
        for ion in self.plasma.ion_species:
            T_ion = self.T[ion][it,:]
            T_electrons = self.T['electrons'][it,:]
            
            T_r = CubicSpline(r_grid,T_ion)
            dTdr_non_filtered = T_r.derivative()
            
            T_polyfit = np.poly1d( np.polyfit(r_grid,T_ion,deg=12) )
            dTdr_polyfit = np.poly1d( T_polyfit.deriv() )
            dTdr_polyfit = dTdr_polyfit(r_grid)
            
            dTdr = dTdr_polyfit  
            # dTdr = dTdr_non_filtered(r_grid)
            
            a_LT = self.aminor * dTdr / T_ion
            
            a_LT_filtered = -a_LT
            
            X = a_LT_filtered - aLT_critical
            
            chi_turb = stiffness * X * np.heaviside(X,1) * (T_electrons/T_ion)**alpha
            
            Bsq = self.Bsq(self.rho_grid)
            
            mi = self.plasma.mass[ion]
            qi = self.plasma.charge[ion]
    
            chi_gB = (EC*T_ion/mi)**1.5 * mi*mi / (qi**2 * Bsq) / self.aminor
            
            chi_turb = chi_gB * chi_turb
            
            chi[ion] = chi_base + chi_turb
        
        for species in self.list_of_species:
            
            p_r = CubicSpline(r_grid,self.P[species][it,:])
            dpdr = p_r.derivative()
            
            n_r = CubicSpline(r_grid,self.N[species][it,:])
            dndr = n_r.derivative()
            
            n_r = self.N[species][it,:]
            dndr = dndr(r_grid)
            
            D = chi[species]
            
            # save D of ALL subiter
            self.Dp_keep[species][it].append(np.array(D))
            
            # average to smooth-out eventual oscillations
            D_avg = np.mean(np.array(self.Dp_keep[species][it]), axis=0)
            D = D_avg

            self.Dp[species][it,:] = D
            
            c = (chi[species]/n_r)*dndr + convective_fact*self.Gamma_turb[species][it,:]/n_r
            c[0] = 0.0
            
            # should we also average 'c' ??
            
            self.cp[species][it,:] = c
            
            # this is for bookeeping
            self.Q_turb[species][it,:] = -chi[species] * dpdr(r_grid) + p_r(r_grid)*( (chi[species]/n_r)*dndr + convective_fact*self.Gamma_turb[species][it,:]/n_r)
            
    def compute_NEO_particle_flux(self,it):
        
        from scipy.interpolate import CubicSpline, Akima1DInterpolator
        
        root = 'ion_root'

        PENTA_class = PENTA(folder_path='.', plasma=self.plasma, lverb=False)
        
        for sp,species in enumerate(self.list_of_species):
            
            # Gamma = np.array( PENTA_class.Gamma_Maxw[species] )
            Gamma = np.array( PENTA_class.Gamma[species,root] )
            roa_PENTA = PENTA_class.roa[root]
            
            ## Include Gamma(r=0) = 0
            rho_extended = np.concatenate([[0.0],roa_PENTA])
            Gamma_extended = np.concatenate(([0.0],Gamma))
            # Gamma_interp = CubicSpline(rho_extended,Gamma_extended,extrapolate=True,bc_type='natural')
            Gamma_interp = Akima1DInterpolator(rho_extended,Gamma_extended,method='makima')
            Gamma_interp.extrapolate = True
            # This is used when computing the heat flux
            self.Gamma_NEO[species][it,:] = Gamma_interp(self.rho_grid)
            
            # Compute Dn and cn
            PENTA_class.set_plasma_solver_transport_coeffs()
            Dn = PENTA_class.Dn[species,root][:,sp] # the sp index picks the self diffusion coeff, Dn_aa
            cn = PENTA_class.cn[species,root][:]
            
            # extended Dn and cn towards the axis by setting them to 0.0
            Dn_extended = np.concatenate(([0.0],Dn))
            cn_extended = np.concatenate(([0.0],cn))
            roa_extended = np.concatenate(([0.0],roa_PENTA))
            
            Dn_extended_spline = Akima1DInterpolator(roa_extended,Dn_extended,method='makima')
            Dn_extended_spline.extrapolate = True
            #
            cn_extended_spline = Akima1DInterpolator(roa_extended,cn_extended,method='makima')
            cn_extended_spline.extrapolate = True
            
            self.Dn[species][it,:] = Dn_extended_spline(self.rho_grid)
            self.cn[species][it,:] = cn_extended_spline(self.rho_grid)
            
    def compute_NEO_heat_flux(self,it):
        
        from scipy.interpolate import CubicSpline, Akima1DInterpolator
        
        root = 'ion_root'

        PENTA_class = PENTA(folder_path='.', plasma=self.plasma, lverb=False)
        
        for sp,species in enumerate(self.list_of_species):
            
            QoT = np.array( PENTA_class.QoT[species,root] )
            roa_PENTA = PENTA_class.roa[root]
            
            T_PENTA = CubicSpline(self.rho_grid, self.T[species][it,:])
            T_PENTA = T_PENTA(roa_PENTA)
            
            Q = QoT * EC * T_PENTA
            
            ## Include Q(r=0) = 0
            rho_extended = np.concatenate([[0.0],roa_PENTA])
            Q_extended = np.concatenate(([0.0],Q))
            
            Q_interp = Akima1DInterpolator(rho_extended,Q_extended,method='makima')
            Q_interp.extrapolate = True
            
            # This is used when computing the heat flux
            self.Q_NEO[species][it,:] = Q_interp(self.rho_grid)
            
            # Compute Dp and cp
            PENTA_class.set_plasma_solver_transport_coeffs()
            Dp = PENTA_class.Dp[species,root][:,sp] # the sp index picks the self diffusion coeff, Dn_aa
            cp = PENTA_class.cp[species,root][:]
            
            # extended Dp and cp towards the axis by setting them to 0.0
            Dp_extended = np.concatenate(([0.0],Dp))
            cp_extended = np.concatenate(([0.0],cp))
            roa_extended = np.concatenate(([0.0],roa_PENTA))
            
            Dp_extended_spline = Akima1DInterpolator(roa_extended,Dp_extended,method='makima')
            Dp_extended_spline.extrapolate = True
            #
            cp_extended_spline = Akima1DInterpolator(roa_extended,cp_extended,method='makima')
            cp_extended_spline.extrapolate = True
            
            self.Dp[species][it,:] = Dp_extended_spline(self.rho_grid)
            self.cp[species][it,:] = cp_extended_spline(self.rho_grid)
            
    def compute_NEO_plus_beurskens_heat_flux(self,it):
        
        from scipy.interpolate import CubicSpline, Akima1DInterpolator
        
        
        ##############################################################################################
        ################################ NEO contribution ############################################
        ##############################################################################################
        
        root = 'ion_root'

        PENTA_class = PENTA(folder_path='.', plasma=self.plasma, lverb=False)
        
        for sp,species in enumerate(self.list_of_species):
            
            QoT = np.array( PENTA_class.QoT[species,root] )
            roa_PENTA = PENTA_class.roa[root]
            
            T_PENTA = CubicSpline(self.rho_grid, self.T[species][it,:])
            T_PENTA = T_PENTA(roa_PENTA)
            
            Q = QoT * EC * T_PENTA
            
            ## Include Q(r=0) = 0
            rho_extended = np.concatenate([[0.0],roa_PENTA])
            Q_extended = np.concatenate(([0.0],Q))
            
            Q_interp = Akima1DInterpolator(rho_extended,Q_extended,method='makima')
            Q_interp.extrapolate = True
            
            # This is used when computing the heat flux
            self.Q_NEO[species][it,:] = Q_interp(self.rho_grid)
            
            # Compute Dp and cp
            PENTA_class.set_plasma_solver_transport_coeffs()
            Dp = PENTA_class.Dp[species,root][:,sp] # the sp index picks the self diffusion coeff, Dn_aa
            cp = PENTA_class.cp[species,root][:]
            
            # extended Dp and cp towards the axis by setting them to 0.0
            Dp_extended = np.concatenate(([0.0],Dp))
            cp_extended = np.concatenate(([0.0],cp))
            roa_extended = np.concatenate(([0.0],roa_PENTA))
            
            Dp_extended_spline = Akima1DInterpolator(roa_extended,Dp_extended,method='makima')
            Dp_extended_spline.extrapolate = True
            #
            cp_extended_spline = Akima1DInterpolator(roa_extended,cp_extended,method='makima')
            cp_extended_spline.extrapolate = True
            
            self.Dp[species][it,:] = Dp_extended_spline(self.rho_grid)
            self.cp[species][it,:] = cp_extended_spline(self.rho_grid)
            
        ##############################################################################################
        ########################## Beurskens contribution ############################################
        ##############################################################################################
        
        chi_electrons = self.heat_fluxes_info['dkespenta_beurskens']['chi_electrons']
        aLT_critical = self.heat_fluxes_info['dkespenta_beurskens']['aLT_critical']
        alpha = self.heat_fluxes_info['dkespenta_beurskens']['alpha']
        stiffness = self.heat_fluxes_info['dkespenta_beurskens']['stiffness']
        convective_fact = self.heat_fluxes_info['dkespenta_beurskens']['convective_fact']
        
        chi = {}
        
        ## electrons
        chi['electrons'] = chi_electrons * np.ones(self.Nr)
        
        r_grid = self.r_grid
        
        ## IONS
        for ion in self.plasma.ion_species:
            T_ion = self.T[ion][it,:]
            T_electrons = self.T['electrons'][it,:]
            
            T_r = CubicSpline(r_grid,T_ion)
            # dTdr_non_filtered = T_r.derivative()
            
            T_polyfit = np.poly1d( np.polyfit(r_grid,T_ion,deg=12) )
            dTdr_polyfit = np.poly1d( T_polyfit.deriv() )
            dTdr_polyfit = dTdr_polyfit(r_grid)
            
            dTdr = dTdr_polyfit  
            # dTdr = dTdr_non_filtered(r_grid)
            
            a_LT = self.aminor * dTdr / T_ion
            
            a_LT_filtered = -a_LT
            
            X = a_LT_filtered - aLT_critical
            
            chi_turb = stiffness * X * np.heaviside(X,1) * (T_electrons/T_ion)**alpha
            
            Bsq = self.Bsq(self.rho_grid)
            
            mi = self.plasma.mass[ion]
            qi = self.plasma.charge[ion]
    
            chi_gB = (EC*T_ion/mi)**1.5 * mi*mi / (qi**2 * Bsq) / self.aminor
            
            chi_turb = chi_gB * chi_turb
            
            chi[ion] = chi_turb
        
        for species in self.list_of_species:
            
            p_r = CubicSpline(r_grid,self.P[species][it,:])
            dpdr = p_r.derivative()
            
            n_r = CubicSpline(r_grid,self.N[species][it,:])
            dndr = n_r.derivative()
            
            n_r = self.N[species][it,:]
            dndr = dndr(r_grid)
            
            D = chi[species]
            
            # save D of ALL subiter
            self.Dp_keep[species][it].append(np.array(D))
            
            # average to smooth-out eventual oscillations
            D_avg = np.mean(np.array(self.Dp_keep[species][it]), axis=0)
            D = D_avg

            # add Beurskens contribution
            self.Dp[species][it,:] += D
            
            c = (chi[species]/n_r)*dndr + convective_fact*self.Gamma_turb[species][it,:]/n_r
            c[0] = 0.0

            # add Beurskens contribution
            self.cp[species][it,:] += c
            
            # this is for bookeeping
            self.Q_turb[species][it,:] = -chi[species] * dpdr(r_grid) + p_r(r_grid)*( (chi[species]/n_r)*dndr + convective_fact*self.Gamma_turb[species][it,:]/n_r)
            
    def solve_density_equations(self,it):
        
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
        
        explicit_source = 0.0
        for source_type in self.explicit_particle_sources[species]:
            explicit_source += self.explicit_particle_sources[species][source_type][it,:]
            
        return explicit_source
    
    def get_explicit_energy_sources(self,species,it):
        
        explicit_source = 0.0
        for source_type in self.explicit_energy_sources[species]:
            explicit_source += self.explicit_energy_sources[species][source_type][it,:]
            
        return explicit_source
    
    def get_LHS_density(self,species,it):
        from scipy.sparse import diags
        
        drho = self.drho
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
        LHS = diags([lower, main, upper], offsets=[-1, 0, 1], format="csr")    
       
        return LHS
    
    def get_LHS_pressure(self,it):
        from scipy.sparse import diags, block_diag, csr_matrix
        
        drho = self.drho
        dr = self.aminor * drho
        Vp = self.dVdr
        Nr = self.Nr
        num_species = len(self.list_of_species)
        
        vp = Vp(self.rho_grid)
        vp_inner = vp[1:-1]
        
        DIFF = {}
        
        for species in self.list_of_species:
            
            Dp = self.Dp[species][it,:]    
            cp = self.cp[species][it,:]
 
            dt_fact = (2./3.)*self.dt
            
            ############################################
            ############### COMPUTE LHS ################
            ############################################
            lower = np.zeros(self.Nr-1)
            main = np.zeros(self.Nr)
            upper = np.zeros(self.Nr-1)
            
            ## 0<r<a (inner grid, no boundary points)
            Dp_plus = (Dp[2:]+Dp[1:-1]) / 2
            Dp_minus = (Dp[0:-2]+Dp[1:-1]) / 2
            
            Vp_plus = (vp[2:]+vp[1:-1]) / 2
            Vp_minus = (vp[0:-2]+vp[1:-1]) / 2
            
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
            
            DIFF[species] = diags([lower, main, upper], offsets=[-1, 0, 1], format="csr")    
        
        DIFF_list = [DIFF[species] for species in self.list_of_species]

        # Construct the block diagonal sparse matrix
        LHS = block_diag(DIFF_list, format="csr")
            
        # Add implicit terms from sources
        sources_implicit = np.zeros((Nr*num_species,Nr*num_species))
        
        if 'Coll_Heat_Exchange' in self.energy_sources['electrons']:
            sources_implicit -= dt_fact*self.get_collisionalHeatExchange(it)
            
        sources_implicit = csr_matrix(sources_implicit)
        LHS = LHS + sources_implicit
        
        # impose Dirichlet boundary condition
        LHS = LHS.tolil()
        for s in range(1, num_species + 1):  # s starts at 1, up to num_species
            row_idx = s * Nr - 1  # Compute the correct row index

            # Set the entire row to zero
            LHS.rows[row_idx] = []  # Clear all column indices in that row
            LHS.data[row_idx] = []  # Clear all values in that row

            # Set the diagonal element to 1
            LHS[row_idx, row_idx] = 1
        
        LHS = LHS.tocsr()
                
        return LHS
    
    def solve_sparse_system(self,matrix,vect):
        
        from scipy.sparse.linalg import spsolve
        
        #import matplotlib.pyplot as plt
        # # plot matrix
        # plt.figure(figsize=(6, 6))
        # plt.spy(matrix, markersize=5, color="black")
        # plt.show()
        
        sol = spsolve(matrix,vect)
        
        return sol
    
    def get_collisionalHeatExchange(self,it):
        # returns collisional heat exchange to use as implicit operator
        
        from libstell.collisions import COLLISIONS
        from scipy import sparse
        
        coll = COLLISIONS()
        
        num_species = len(self.list_of_species)
        
        clog = np.zeros(self.Nr)
        
        W_s1_s2 = np.zeros((num_species,num_species,self.Nr))
        aux_B = np.zeros((num_species,num_species,self.Nr))
        
        for is1,species1 in enumerate(self.list_of_species):

            m1 = self.plasma.mass[species1]
            Z1 = self.plasma.Zcharge[species1]
            n1 = self.N[species1][it,:]
            T1 = self.T[species1][it,:]
            
            for is2,species2 in enumerate(self.list_of_species):
                
                m2 = self.plasma.mass[species2]
                Z2 = self.plasma.Zcharge[species2]
                n2 = self.N[species2][it,:]
                T2 = self.T[species2][it,:]
                
                # get Coulomb logarithm
                # for ir in range(self.Nr):
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
                
                gamma = const / den

                W_s1_s2[is1,is2,:] = gamma*n1  
                
                aux_B[is1,is2,:] = gamma*n2

        # Add aux_B matrix
        W_s1_s2[np.arange(num_species), np.arange(num_species), :] -= np.sum(aux_B, axis=1)

        # W_out matrix
        W_out = np.zeros((num_species*self.Nr,num_species*self.Nr))
        
        j=0
        for is1 in range(num_species):
            for ir1 in range(self.Nr):
                p=0
                for is2 in range(num_species):
                    for ir2 in range(self.Nr):
                        if(ir1==ir2):
                            W_out[j,p] = W_s1_s2[is1,is2,ir2]  
                        p=p+1
                j = j+1
                
        W_out = sparse.csr_matrix(W_out)

        return W_out
    
    def call_PENTA3(self,it):
        import subprocess
        from concurrent.futures import ProcessPoolExecutor, as_completed
        import functools
        
        # create PLASMA class in order to write PENTA inputs       
        plasma_PENTA = PLASMA(self.list_of_species)
        for species in self.list_of_species:
            plasma_PENTA.set_density(species,'interp',rho_vals=self.rho_grid,n_vals=self.N[species][it,:])
            plasma_PENTA.set_temperature(species,'interp',rho_vals=self.rho_grid,T_vals=self.T[species][it,:])
        
        plasma_profiles_extension = 'transp_solver'
        plasma_PENTA.write_plasma_profiles_to_PENTA3(filename='plasma_profiles_'+plasma_profiles_extension+'.dat')
        plasma_PENTA.write_PENTA_namelist()

        try:
            surfaces = self.particle_fluxes_info['dkespenta']['surfaces']
        except:
            try:
                surfaces = self.heat_fluxes_info['dkespenta']['surfaces']
            except:
                try:
                    surfaces = self.particle_fluxes_info['dkespenta_beurskens']['surfaces']
                except:
                    surfaces = self.heat_fluxes_info['dkespenta_beurskens']['surfaces']
        
        time_sec = []
        with ProcessPoolExecutor() as executor:
            futures = [executor.submit(process_surfaces, surface, self.wout_path) for surface in surfaces]

            for future in as_completed(futures):
                elapsed_seconds = future.result()
                time_sec.append(elapsed_seconds)
                # print(f'Surface processed in {elapsed_seconds:.2f} seconds')
            
        # delete files not needed
        remove = 'rm ucontra* sigmas* flows_vs_Er*'
        subprocess.run(remove, shell=True, check=True, text=True, capture_output=True)
        
        # merge _surface_# files into single file
        merge_and_delete('fluxes_vs_roa_surface*','fluxes_vs_roa')
        merge_and_delete('fluxes_vs_Er_surface*','fluxes_vs_Er')
        merge_and_delete('flows_vs_roa_surface*','flows_vs_roa')
        merge_and_delete('Jprl_vs_roa_surface*','Jprl_vs_roa')
        merge_and_delete('particleTransportCoeffs_vs_roa_surface*','particleTransportCoeffs_vs_roa')
        merge_and_delete('heatTransportCoeffs_vs_roa_surface*','heatTransportCoeffs_vs_roa')
        merge_and_delete('plasma_profiles_check_surface*','plasma_profiles_check')
        
    def call_save_output(self,output_filename):
        # saves in joblib file
        from types import SimpleNamespace
        from pathlib import Path
        import joblib
        
        # check if extension of output_filename is .joblib; if not, add
        output_filename = str(Path(output_filename).with_suffix(".joblib"))
        
        # save the class (cannot save solver directly cause it contains lambda functions...)
        saved_class = SimpleNamespace()
        saved_class.rho_grid = self.rho_grid
        saved_class.r_grid = self.r_grid
        saved_class.dVdr = self.dVdr
        saved_class.aminor = self.aminor
        saved_class.Rmajor = self.Rmajor
        saved_class.B = self.B0
        saved_class.list_of_species = self.list_of_species
        
        # only save at minimum every dt=0.1s 
        freq = max(1, round(0.1 / self.dt))
        sl = slice(0, -1, freq)  # defines the slice once

        saved_class.time = self.time[sl]
        saved_class.Nt = len(self.time[sl])

        for attr in ('N','T','Dn','cn','Dp','cp','Q_NEO','Q_turb','Gamma_NEO','Gamma_turb'):
            setattr(saved_class, attr, {})
            for species in self.list_of_species:
                getattr(saved_class, attr)[species] = getattr(self, attr)[species][sl, :]
                
        # nested dict attributes
        nested_attrs = ['explicit_energy_sources','explicit_particle_sources']
        for species in self.list_of_species:
            for attr in nested_attrs:
                saved_class.__dict__.setdefault(attr, {})
                saved_class.__dict__[attr].setdefault(species, {})
                for type_string, arr in getattr(self, attr)[species].items():
                    saved_class.__dict__[attr][species][type_string] = arr[sl, :]

        joblib.dump(saved_class, output_filename)
        
def process_surfaces(surface,wout_path):
    import time
    import subprocess
    
    start_time = time.time()
    
    type_of_write = 0
    Er_min_V_cm = -100
    Er_max_V_cm = 200

    # wout_path = solver_class.wout_path
    EparB = 0.0

    #Sonine (Laguerre) polynomials
    Smax = 1
    
    plasma_profiles_extension = 'transp_solver'
    
    extension_output_files = f'_surface_{surface}'
        
    extension_star_files = f'surface_{surface}'

    call_penta3 = f'~/bin/xpenta {extension_star_files} {Er_min_V_cm} {Er_max_V_cm} {surface} {type_of_write} {wout_path} {plasma_profiles_extension} {EparB} {Smax} {extension_output_files}'
    
    result = subprocess.run(call_penta3, shell=True, check=True, text=True, capture_output=True)
    # print(result.stdout)
    if(result.stderr):
        print(result.stderr)
        
    end_time = time.time()
    elapsed_time = (end_time - start_time)
    return elapsed_time
                
def merge_and_delete(pattern, output_filename):
    """
    Merges files matching the given pattern into a single file and deletes the originals.
    
    Parameters:
    pattern (str): The pattern to match files (e.g., 'fluxes_vs_roa_surface_*').
    output_filename (str): The name of the output file.
    """
    import os
    import re
    import glob
    
    def extract_number(filename):
        match = re.search(r'_(\d+)$', filename)  # Extract number at the end
        return int(match.group(1)) if match else float('inf')

    # Find and sort matching files
    file_list = glob.glob(pattern)
    file_list.sort(key=extract_number)

    if not file_list:
        print(f"No files found matching pattern: {pattern}")
        return

    header_written = False
    with open(output_filename, 'w') as outfile:
        for filename in file_list:
            with open(filename, 'r') as infile:
                lines = infile.readlines()
                if not header_written:
                    outfile.write(lines[0])  # Write header
                    outfile.write(lines[1])
                    header_written = True
                outfile.writelines(lines[2:])  # Write data

    # Delete original files
    for filename in file_list:
        os.remove(filename)
        
            
# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)      