"""
This library provides a python class for solving pressure transport
equations
"""

import numpy as np
import sys
from time import perf_counter

# Constants
EC = 1.602176634E-19 # Electron charge [C]
EPS0 = 8.8541878188E-12 # Vacuum permittivity [F/m]

# PENTA Class
class PRESSURE_SOLVER_FULL_MATRIX:
    
    def __init__(self, plasma_class, density_functions_dict=None):
        # density_function must be as a function of (t,rho)
        # if density_functions_dict is given, it means n is not fixed in time
        # density_functions_dict is a dictionary of n(t,rho) functions for each
        
        from collections import defaultdict
        
        self.list_of_species = plasma_class.list_of_species
        
        if(density_functions_dict is not None):
            # make check on the dictionary
            self.check_density_function(density_functions_dict)
            
            # if density_function is given, density should not be set in plasma class
            for species in self.list_of_species:
                if(species in plasma_class.density):
                    print(f'ERROR" density of {species} SHOULD NOT BE SET! OTHERWISE, DO NOT GIVE density_function')
                    exit(0)
            
            self.external_density_given = True
            self.density_funcs = density_functions_dict
            
        else:
            # Check density profiles exist for all species
            for species in self.list_of_species:
                if(species not in plasma_class.density):
                    print(f'ERROR: density of {species} has not been set yet')
                    exit(0)
            
            self.external_density_given = False
                
        self.plasma = plasma_class
                
        # Initialize dicitionaries
        self.edge_bnd_cnd = {}
        self.initial_profile = {}
        self.sources = defaultdict(lambda: defaultdict(dict))
        self.fluxes_info = defaultdict(lambda: defaultdict(dict))
                
        print(f'Solvers for pressure of {self.list_of_species} INITIALIZED!')
        
    def check_density_function(self,density_functions_dict):
        from inspect import signature
        
        # Check if it's a dictionary
        if not isinstance(density_functions_dict, dict):
            raise TypeError("density_functions_dict must be a dictionary")
    
        # Check that all keys are in species_list
        invalid_keys = set(density_functions_dict.keys()) - set(self.list_of_species)
        if invalid_keys:
            raise ValueError(f"Invalid species in dictionary: {invalid_keys}")
        
        # Check that all species are present in the dictionary
        missing_keys = set(self.list_of_species) - set(density_functions_dict.keys())
        if missing_keys:
            raise ValueError(f"Missing species in dictionary: {missing_keys}")

        # Check that all values are functions that accept exactly 2 arguments
        for species, func in density_functions_dict.items():
            if not callable(func):
                raise TypeError(f"density_functions_dict[{species}] is not a function")

        # sig = signature(func)
        # if len(sig.parameters) != 2:
        #     raise TypeError(f"density_functions_dict[{species}] must accept exactly 2 arguments")
                
    def set_edge_boundary_condition(self,species: str,val: float):
        
        # check species exist in list_of_species
        if species not in self.list_of_species:
            print(f"ERROR: Species {species} is not in the plasma.")
            exit(1)
            
        self.edge_bnd_cnd[species] = val
            
    def set_initial_profile(self,species: str, rho_vals, profile_vals):
        
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
        self.initial_profile[species] = CubicSpline(rho_vals,profile_vals)
        
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

    def set_source(self,species,source_type, total_power=None, sigma_rho=None, fraction_alpha_heating=None, cte_source=None, interpolant_2D=None, time_dependent_factor=None, tau_alphas=None, tau_palphas=None, Tion_threshold=None, lambda_function_2D=None):
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
                    self.sources[species][source_type] = {}
            case 'Bremsstrahlung_alphas':
                if(species != 'electrons'): 
                    print('ERROR: Bremsstrahlung is only source for electrons')
                    exit(0)
                if( (tau_alphas is None) or (tau_palphas is None)):
                    print('ERROR: tau_alphas and tau_palphas must be given (in seconds)')
                    exit(0)
                else:
                    self.sources[species][source_type] = {'tau_alphas': tau_alphas, 'tau_palphas' : tau_palphas}
            case 'external_gaussian':
                if((total_power is None) or (sigma_rho is None)):
                    print('ERROR: Need to provide total_power [W] and sigma_rho for gaussian external source')
                    exit(1) 
                else:
                    self.sources[species][source_type] = {'total_power' : total_power, 'sigma_rho' : sigma_rho, 'Tion_threshold' : Tion_threshold }
            case 'time_dependent_gaussian':
                if((total_power is None) or (sigma_rho is None) or (time_dependent_factor is None)):
                    print('ERROR: Need to provide total_power [W], sigma_rho and a time depenedent factof for time-dependent gaussian')
                    exit(1) 
                else:
                    self.sources[species][source_type] = {'total_power' : total_power, 'sigma_rho' : sigma_rho, 'time_factor': time_dependent_factor }
            case 'Coll_Heat_Exchange':
                self.sources[species][source_type] = {}
            case 'Er':
                self.sources[species][source_type] = {}
            case 'alpha_heating':
                if(fraction_alpha_heating is None):
                    print('ERROR: fraction_alpha_heating is needed. Usually ~80% electrons, 20% ions')
                    exit(0)
                else:
                    self.sources[species][source_type] = {'fraction_alpha_heating': fraction_alpha_heating}
            case 'constant':
                if(cte_source is None):
                    print('ERROR: cte_source is needed in order to generate a constant source.')
                    exit(0)
                else:
                    self.sources[species][source_type] = {'cte_source' : cte_source}
            case 'interpolant_2D':
                if(interpolant_2D is None):
                    print('ERROR: 2D interpolanting function (interpolant) must be provided')
                    exit(0)
                num_args = len(inspect.signature(interpolant_2D).parameters)
                if(num_args != 2):
                    print('ERROR" interoplating function must have 2 args: time and space')
                    exit(0)
                else:
                    self.sources[species][source_type] = {'interpolant_2D' : interpolant_2D}
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
                    self.sources[species][source_type] = {'lambda_function_2D' : lambda_function_2D}
            case _:
                print(f'ERROR: Source type {source_type} is NOT possible')
                exit(0)
                
    def set_fluxes(self, type: str,dkes_folder=None,surfaces=None,theta=None,chi=None,chi_base=None,aLT_critical=None,alpha=None,stiffness=None,chi_electrons=None):
        # sets type of fluxes
        # OPTION1: type='dkespenta'; dkes_folder and surfaces(list of integers) must be provided
        # OPTION2: type='diffusive'; Dn and chi must be provided (partical and heat collisional diffusion coefficients)
            
        match type:
            case 'dkespenta':
                #checks that dkes_folder and surfaces are provided
                if((dkes_folder is None) or (surfaces is None) or (theta is None)):
                    print('ERROR: dkes_folder, surfaces and theta must be provided!!')
                    exit(0)
                    
                self.fluxes_info['type'] = type
                self.fluxes_info[type]['dkes_folder'] = dkes_folder
                self.fluxes_info[type]['surfaces'] = surfaces
                self.fluxes_info[type]['theta'] = theta
                
            case 'diffusive':
                #checks that diffusion coefficients are provided
                if(chi is None):
                    print('ERROR: chi must be provided!')
                    exit(0)
                self.fluxes_info['type'] = type
                self.fluxes_info[type]['chi'] = chi
                
            case 'beurskens':
                if( (chi_base is None) or (aLT_critical is None) or (alpha is None) or (stiffness is None) or (chi_electrons is None) ):
                    print('ERROR: chi_base, aLT_critical, alpha, stiffness and chi_electrons must be given!')
                    exit(0)
                self.fluxes_info['type'] = type
                self.fluxes_info[type]['chi_base'] = chi_base
                self.fluxes_info[type]['aLT_critical'] = aLT_critical
                self.fluxes_info[type]['alpha'] = alpha
                self.fluxes_info[type]['stiffness'] = stiffness
                self.fluxes_info[type]['chi_electrons'] = chi_electrons
                    
    def set_density_and_temperature(self,species,pressure,it):
        # from pressure (Pa), sets temperature (eV) in plasma class
        # sets density (m^-3) in plasma class as well
        
        rho = self.rho_grid
        
        if(self.external_density_given):
            dens_func = self.density_funcs[species]
            density = dens_func(t=self.time[it],rho=rho)
            # print(density)
            self.plasma.set_density(species,'interp',rho_vals=rho,n_vals=density)    
        else:
            density = self.plasma.get_density(species,rho)
            
        temperature =  pressure / (EC*density) # eV
        
        self.plasma.set_temperature(species,'interp',rho_vals=rho,T_vals=temperature) 
                 
    def run(self,Nr,dt,tstart,tend,tolerance=1E-2,max_subiter=12):
        
        from collections import defaultdict
        
        # Check everything is set and ready to proceed with the run
        self.make_checks()
        
        # Initialize rho grid
        rho = np.linspace(0,1,Nr)
        drho = rho[1]-rho[0]
        
        # Initialize time
        # check Nt=1+(tend-tstart)/dt is an integer
        # if( ((tend-tstart)/dt)%1 != 0 ):
        #     print(f'ERROR: dt not compatible with tstart and tend -- {((tend-tstart)/dt)%1}')
        #     exit(0)
        # else:
        Nt = round( 1+(tend-tstart)/dt )
        time = np.linspace(tstart,tend,Nt)
            
        # save grids in class
        self.time = time
        self.Nt = Nt
        self.dt = dt
        self.rho_grid = rho
        self.r_grid = rho * self.aminor
        self.drho = drho
        self.dr = drho * self.aminor
        self.Nr = Nr
            
        print(' ')
        print( ' ***********************')
        print(f' *  tstart = {tstart:5.2f}s    *')
        print(f' *  tend   = {tend:5.2f}s    *')
        print(f' *  dt     = {dt:5.2f}s    *')
        print(f' *  Nt     = {Nt:3}       *')
        print(f' *  drho   = {drho:5.2f}     *')
        print( ' ***********************')
            
        ################## INITIALIZE VARIABLES ########################################
        press = []
        self.nsubiter = np.zeros(Nt)
        self.N = {}
        self.P = {}
        self.T = {}
        self.Q_interp = {}
        self.D_interp = {}
        self.D_keep = defaultdict(lambda: defaultdict(list))
        self.c_interp = {}
        self.total_sources_explicit = {}
        self.all_sources = {}
        self.ECRH_off = False
        for species in self.list_of_species:
            p_init = self.initial_profile[species](rho)
            
            self.set_density_and_temperature(species,p_init,it=0)
            
            self.P[species] = np.zeros((Nt,Nr))
            self.P[species][0,:] = p_init
            
            press.append(p_init)
            
            self.T[species] = np.zeros((Nt,Nr))
            self.T[species][0,:] = self.plasma.get_temperature(species,rho)
            
            self.N[species] = np.zeros((Nt,Nr))
            self.N[species][0,:] = self.plasma.get_density(species,rho)
            
            self.D_interp[species] = [None]*Nt
            self.c_interp[species] = [None]*Nt
            self.Q_interp[species] = [None]*Nt
                        
            self.total_sources_explicit[species] = np.zeros((Nt,Nr))

            # Initialize arrays for each source_type in the specified species
            self.all_sources[species] = {}
            for source_type in self.sources[species].keys():
                self.all_sources[species][source_type] = np.zeros((Nt, Nr))
                
        if('Bremsstrahlung_alphas' in self.sources['electrons']):
            if(tstart>0):
                print('ERROR: Can only use Bremsstrahlung_alphas when tstart=0 because Nalphas from previous iter does not exist...')
                exit(0)
            self.Nalphas_fast = np.zeros((Nt,Nr))
            self.Nalphas_thermal = np.zeros((Nt,Nr))
            
        press = np.concatenate(press)
        p_old = press # 1E3*np.ones(Nr*len(self.list_of_species)) # so on loop 1 we don't divide by zero in delta_p
        
        ####################################################################################
                
        print(' ')
        header_str = '  TIME [s]     NSUB      TE_AXIS [keV]     SE_AXIS [MW/m^3]    TI1_AXIS [keV]    SI1_AXIS [MW/m^3]    MAX(dp/p_old)'  
        print(header_str)
        print('  '+'='*len(header_str))
        
        # t=t_start
        self.nsubiter[0] = 1
        self.call_fluxes(it=0)
        for species in self.list_of_species:
                self.total_sources_explicit[species][0,:] = self.get_sources_explicit(species,rho,it=0)
        info_str = f'  {tstart:<13.2f}{1:<10}{self.T['electrons'][0,0]/1E3:<18.3f}{'---':<20}{self.T['deuterium'][0,0]/1E3:<18.3f}{'---':<21}{0.0:<13.2E}'
        print(info_str)
        
        time_update_pressure_temperature = 0.0
        time_call_fluxes = 0.0
        time_get_sources_explicit = 0.0
        time_get_RHS_vector = 0.0
        time_get_LHS_matrix = 0.0
        self.time_get_collisionalHeatExchange = 0.0
        self.time_diffusion_matrix_build = 0.0
        self.time_set_LHS_boundary_conditions = 0.0
        time_solve_sparse_system = 0.0  
        
        ### LOOP IN TIME STARTING AT t=tstart+dt ###
        for it,t in enumerate(time[1:],start=1):
            self.it = it
            
            # the first subiter corresponds to the last time iteration
            start_time = perf_counter()
            self.update_pressure_temperature(it,press)
            end_time = perf_counter()
            time_update_pressure_temperature += end_time-start_time
            
            ### SUBCYCLE
            delta_p = 10*tolerance
            subiter=1
            while(delta_p > tolerance and subiter<=max_subiter):
                self.nsubiter[it] += 1  

                # compute self.D_interp and self.c_interp
                start_time = perf_counter()
                self.call_fluxes(it)
                end_time = perf_counter()
                time_call_fluxes += end_time-start_time
                
                # get explicit sources
                start_time = perf_counter()
                for species in self.list_of_species:
                    self.total_sources_explicit[species][it,:] = self.get_sources_explicit(species,rho,it)  # W/m^3
                end_time = perf_counter()
                time_get_sources_explicit += end_time-start_time
                
                start_time = perf_counter()                      
                RHS_vector = self.get_RHS_vector(it)
                end_time = perf_counter()
                time_get_RHS_vector += end_time-start_time
                #
                start_time = perf_counter()                      
                LHS_matrix = self.get_LHS_matrix(it)
                end_time = perf_counter()
                time_get_LHS_matrix += end_time-start_time
                
                # solve system
                start_time = perf_counter()
                press = self.solve_sparse_system(LHS_matrix,RHS_vector)
                end_time = perf_counter()
                time_solve_sparse_system += end_time-start_time
                
                #PICARD FACTOR
                # fpicard = 0.75
                # press = press*fpicard + (1-fpicard)*p_old
                
                start_time = perf_counter()
                self.update_pressure_temperature(it,press) 
                end_time = perf_counter()
                time_update_pressure_temperature += end_time-start_time
            
                delta_p = np.max( np.where( p_old>1E-10, np.abs((press-p_old)/p_old), 0 ) )
                
                p_old = press
                
                info_str = f'  {t:<13.3f}{subiter:<10}{self.T['electrons'][it,0]/1E3:<18.3f}{self.total_sources_explicit['electrons'][it,0]/1E6:<20.2E}{self.T['deuterium'][it,0]/1E3:<18.3f}{self.total_sources_explicit['deuterium'][it,0]/1E6:<21.2E}{delta_p:<13.2E}'
                print(info_str)
                
                subiter += 1
                
        ### print TIMINGS ####
        print(f' ')
        print(f'***** TIMINGS *****')
        total_time = time_update_pressure_temperature+time_call_fluxes+time_get_sources_explicit+time_get_RHS_vector+time_get_LHS_matrix+time_solve_sparse_system
        print(f'time_update_pressure_temperature = {time_update_pressure_temperature:.1f}s [{time_update_pressure_temperature/total_time*100:.1f}%]')
        print(f'time_call_fluxes = {time_call_fluxes:.1f}s [{time_call_fluxes/total_time*100:.1f}%]')
        print(f'time_get_sources_explicit = {time_get_sources_explicit:.1f}s [{time_get_sources_explicit/total_time*100:.1f}%]')
        print(f'time_get_RHS_vector = {time_get_RHS_vector:.1f}s [{time_get_RHS_vector/total_time*100:.1f}%]')
        print(f'time_get_LHS_matrix = {time_get_LHS_matrix:.1f}s [{time_get_LHS_matrix/total_time*100:.1f}%]')
        print(f'    time_get_collisionalHeatExchange = {self.time_get_collisionalHeatExchange:.1f}s [{self.time_get_collisionalHeatExchange/total_time*100:.1f}%]')
        print(f'    time_diffusion_matrix_build = {self.time_diffusion_matrix_build:.1f}s [{self.time_diffusion_matrix_build/total_time*100:.1f}%]')
        print(f'    time_set_LHS_boundary_conditions = {self.time_set_LHS_boundary_conditions:.1f}s [{self.time_set_LHS_boundary_conditions/total_time*100:.1f}%]')
        print(f'time_solve_sparse_system = {time_solve_sparse_system:.1f}s [{time_solve_sparse_system/total_time*100:.1f}%]')
        print(f'TOTAL TIME = {total_time/60:.2f}min')
        print(f' ')

    def update_pressure_temperature(self,it,press):
        # updates self.P and self.T
        
        Nr = self.Nr
        
        k=0
        for species in self.list_of_species:
            
            self.P[species][it,:] = press[k:(k+Nr)]
            
            self.set_density_and_temperature(species,self.P[species][it,:],it)
            self.T[species][it,:] = self.plasma.get_temperature(species,self.rho_grid)
            self.N[species][it,:] = self.plasma.get_density(species,self.rho_grid)
            
            k = k+Nr     
                        
    def get_sources_explicit(self,species: str,rho_grid, it):
        # returns 1D-array of same size as rho_grid
        # computes sources using info in self.sources[species]
        # 'it' necessary to save in self.all_sources
        
        from libstell.fusion import FUSION
        from collisions import COLLISIONS
        
        fusion = FUSION()
        
        # checks rho_grid is a 1D array
        rho_grid = np.asarray(rho_grid)
        if( rho_grid.ndim != 1):
            raise ValueError("rho_grid must be 1D")
        
        total_source_explicit = np.zeros(len(rho_grid))
        
        for source_type in self.sources[species]:
            
            aux_source = 0.0
                     
            match source_type:
                case 'Bremsstrahlung':
                    for ion in self.plasma.ion_species:
                        zi = self.plasma.Zcharge[ion]
                        ni = self.plasma.get_density(ion,rho_grid)
                        ne = self.plasma.get_density('electrons',rho_grid)
                        Te = self.plasma.get_temperature('electrons',rho_grid)
                        
                        aux_source -= fusion.BremsstrahlungPower(zi,ni,ne,Te)
                            
                    #save in dictionary for bookeeping
                    self.all_sources[species][source_type][it,:] = aux_source
                    
                case 'Bremsstrahlung_alphas':
                    ne = self.plasma.get_density('electrons',rho_grid)
                    Te = self.plasma.get_temperature('electrons',rho_grid)
                    Z_alpha = 2
                    tau_alpha = self.sources['electrons']['Bremsstrahlung_alphas']['tau_alphas']
                    tau_palpha = self.sources['electrons']['Bremsstrahlung_alphas']['tau_palphas']
                    nD = self.plasma.get_density('deuterium',rho_grid)
                    nT = self.plasma.get_density('tritium',rho_grid)
                    Ti = 0.5* ( self.plasma.get_temperature('deuterium', rho_grid) + self.plasma.get_temperature('tritium', rho_grid) )
                    sigmav = fusion.sigmaBH(Ti,'DT') # m^3/s
                    
                    self.Nalphas_fast[it,:] = (self.Nalphas_fast[it-1,:] + self.dt*nD*nT*sigmav) / (1+self.dt/tau_alpha)
                    self.Nalphas_thermal[it,:] = (self.Nalphas_thermal[it,:] + (self.dt/tau_alpha)*self.Nalphas_fast[it,:]) / (1+self.dt/tau_palpha)
                    
                    aux_source -= fusion.BremsstrahlungPower(Z_alpha,self.Nalphas_thermal[it,:],ne,Te)
                        
                    #save in dictionary for bookeeping
                    self.all_sources[species][source_type][it,:] = aux_source 
                
                case 'external_gaussian':
                    r0 = 0.0
                    sigma_rho = self.sources[species]['external_gaussian']['sigma_rho']
                    sigma_r = sigma_rho*self.aminor
                    r = self.rho_grid * self.aminor
                    P_IN = self.sources[species]['external_gaussian']['total_power']
                    #
                    integrand = np.exp(-(r-r0)**2/sigma_r**2) * self.dVdr(self.rho_grid)
                    integrand = integrand.flatten()
                    #
                    cte = P_IN / np.trapz(integrand,r)
                    #
                    aux_source = cte * np.exp(-(r-r0)**2/sigma_r**2)
                    
                    Tion_threshold = self.sources[species]['external_gaussian']['Tion_threshold']
                    if(Tion_threshold is not None):
                        Tion = (self.T['deuterium'][it-1,0]+self.T['tritium'][it-1,0])/2
                        if(Tion > Tion_threshold or self.ECRH_off):
                            aux_source = 0.0
                            # self.ECRH_off = True
                    
                    #save in dictionary for bookeeping
                    self.all_sources[species][source_type][it,:] = aux_source
                    
                case 'time_dependent_gaussian':
                    r0 = 0.0
                    sigma_rho = self.sources[species]['time_dependent_gaussian']['sigma_rho']
                    sigma_r = sigma_rho*self.aminor
                    r = self.rho_grid * self.aminor
                    P_IN = self.sources[species]['time_dependent_gaussian']['total_power']
                    #
                    integrand = np.exp(-(r-r0)**2/sigma_r**2) * self.dVdr(self.rho_grid)
                    integrand = integrand.flatten()
                    #
                    cte = P_IN / np.trapz(integrand,r)
                    #
                    time_fact = self.sources[species]['time_dependent_gaussian']['time_factor']
                    t = self.time[it]
                    aux_source = cte * np.exp(-(r-r0)**2/sigma_r**2) * time_fact(t)
                    
                    # print(f'INTEGRAL SOURCE = {np.trapz(cte * np.exp(-(r-r0)**2/sigma_r**2)*self.dVdr(self.rho_grid),r):.1E}')
                    # print(f'VOLUME = {np.trapz(self.dVdr(self.rho_grid),r)}')
                    
                    #save in dictionary for bookeeping
                    self.all_sources[species][source_type][it,:] = aux_source
                    
                case 'Coll_Heat_Exchange':
                    
                    # W_explicit,W_implicit = self.get_collisionalHeatExchange(species,it)
                    
                    # W_bookeeping = np.sum(W_explicit,axis=0) + np.sum(W_implicit,axis=0)*self.P[species][it,:]

                    # sum_W_explicit = np.sum(W_explicit,axis=0)
                            
                    aux_source = 0.0 #sum_W_explicit
                    
                    #save in dictionary for bookeeping
                    self.all_sources[species][source_type][it,:] = 0.0 #W_bookeeping
                            
                case 'alpha_heating':
                    nD = self.plasma.get_density('deuterium', rho_grid)
                    nT = self.plasma.get_density('tritium', rho_grid)
                    
                    # Ti = 0.5* ( self.plasma.get_temperature('deuterium', rho_grid) + self.plasma.get_temperature('tritium', rho_grid) )
                    
                    # sigmav = [fusion.sigmaBH(ti,'DT') for ti in Ti]
                    
                    # S_alpha = nD * nT * sigmav *  fusion.E_DT_He # W/m^3
                    
                    fraction_alpha_heating = self.sources[species]['alpha_heating']['fraction_alpha_heating']
                    
                    TD = self.plasma.get_temperature('deuterium', rho_grid)
                    TT = self.plasma.get_temperature('tritium', rho_grid)
                    aux_source = fraction_alpha_heating * fusion.alphaPower(nD,nT,TD,TT)
                    
                    # aux_source = S_alpha * fraction_alpha_heating
                    
                    #save in dictionary for bookeeping
                    self.all_sources[species][source_type][it,:] = aux_source
                    
                case 'Er':
                    aux_source = self.plasma.charge[species]*self.Er_interp(rho_grid)*self.Gamma_interp[species](rho_grid)
                    
                    #save in dictionary for bookeeping
                    self.all_sources[species][source_type][it,:] = aux_source
                    
                case 'constant':
                    aux_source = self.sources[species]['constant']['cte_source']
                    
                    #save in dictionary for bookeeping
                    self.all_sources[species][source_type][it,:] = aux_source
                    
                case 'interpolant_2D':
                    interpolating_func = self.sources[species]['interpolant_2D']
                    t = self.time[it]
                    aux_source = interpolating_func(t,self.rho_grid)
                    
                    #save in dictionary for bookeeping
                    self.all_sources[species][source_type][it,:] = aux_source
                    
                case 'lambda_2D':
                    lambda_function_2D = self.sources[species][source_type]['lambda_function_2D'] #func(r,t)
                    
                    aux_source = np.zeros(len(rho_grid))
                    for ir,rho in enumerate(self.rho_grid):
                        aux_source[ir] = lambda_function_2D(rho*self.aminor,self.time[it])   

                case _:
                    print(f'ERROR: Source type {source_type} not defined....')
                    exit(0)
                    
            total_source_explicit += aux_source
                    
                    
        return total_source_explicit
    
    def get_LHS_matrix_non_optimized(self,it):
        # returns
        
        from scipy.interpolate import CubicSpline
        from scipy.sparse import diags, block_diag, csr_matrix
        # from scipy.sparse.linalg import eigs
        
        drho = self.rho_grid[1]-self.rho_grid[0]
        dr = self.aminor * drho
        Vp = self.dVdr
        a = self.aminor
        Nr = self.Nr
        num_species = len(self.list_of_species)
        
        DIFF = {}
        
        start_time = perf_counter()
        for species in self.list_of_species:
            
            D_interp = self.D_interp[species][it]
            c_interp = self.c_interp[species][it]
            
            dt_fact = (2./3.)*self.dt
            
            ############################################
            ############### COMPUTE LHS ################
            ############################################
            lower = np.zeros(self.Nr-1)
            main = np.zeros(self.Nr)
            upper = np.zeros(self.Nr-1)
            
            ## r=0
            main[0] = 1.0 + dt_fact*( 4*D_interp(0)/dr**2 + 2*c_interp(drho)/dr )
            upper[0] = -4*dt_fact*D_interp(0)/dr**2
            
            ## 0<r<a
            for ir,rho in enumerate(self.rho_grid[1:-1],start=1):
                
                rplus = rho + drho/2
                rminus = rho - drho/2
                
                VDplus = Vp(rplus)*D_interp(rplus) / (Vp(rho)*dr**2)
                VDminus = Vp(rminus)*D_interp(rminus) / (Vp(rho)*dr**2)
                
                cplus  = c_interp(rho+drho)*Vp(rho+drho) / (2*Vp(rho)*dr)
                cminus = c_interp(rho-drho)*Vp(rho-drho) / (2*Vp(rho)*dr)
                
                main[ir] = 1.0 + dt_fact*(VDplus+VDminus)
                upper[ir] = dt_fact*(-VDplus+cplus)
                lower[ir-1] = dt_fact*(-VDminus-cminus)
                
                ### UPWIND SCHEME FOR ADVECTION ###
                # cplus  = c_interp(rho+drho)*Vp(rho+drho) / (Vp(rho)*dr)
                # cminus = c_interp(rho-drho)*Vp(rho-drho) / (Vp(rho)*dr)
                
                # main[ir] = 1.0 + dt_fact*(VDplus+VDminus+cplus)
                # upper[ir] = dt_fact*(-VDplus)
                # lower[ir-1] = dt_fact*(-VDminus-cminus)

            ## r=1  -- this can now be removed since Dirichlet BC/s are imposed at the end of this routine
            main[-1] = 1.0
            lower[-1] = 0.0
            
            ####
            DIFF[species] = diags([lower, main, upper], offsets=[-1, 0, 1], format="csr")    
        
        DIFF_list = [DIFF[species] for species in self.list_of_species]

        # Construct the block diagonal sparse matrix
        LHS = block_diag(DIFF_list, format="csr")
        
        end_time = perf_counter()
        self.time_diffusion_matrix_build += end_time-start_time
            
        # Add implicit terms from sources
        sources_implicit = np.zeros((Nr*num_species,Nr*num_species))
        
        start_time = perf_counter()
        if 'Coll_Heat_Exchange' in self.sources['electrons']:
            sources_implicit -= dt_fact*self.get_collisionalHeatExchange()
        end_time = perf_counter()
        self.time_get_collisionalHeatExchange += end_time-start_time
            
        sources_implicit = csr_matrix(sources_implicit)
        LHS = LHS + sources_implicit
        
        # impose Dirichlet boundary condition
        start_time = perf_counter()
        LHS = LHS.tolil()
        for s in range(1, num_species + 1):  # s starts at 1, up to num_species
            row_idx = s * Nr - 1  # Compute the correct row index

            # Set the entire row to zero
            LHS.rows[row_idx] = []  # Clear all column indices in that row
            LHS.data[row_idx] = []  # Clear all values in that row

            # Set the diagonal element to 1
            LHS[row_idx, row_idx] = 1
        
        LHS = LHS.tocsr()
        
        end_time = perf_counter()
        self.time_set_LHS_boundary_conditions += end_time-start_time
                
        return LHS
    
    def get_LHS_matrix(self,it):
        # returns
        
        from scipy.interpolate import CubicSpline
        from scipy.sparse import diags, block_diag, csr_matrix
        # from scipy.sparse.linalg import eigs
        
        drho = self.drho
        dr = self.aminor * drho
        Vp = self.dVdr
        Nr = self.Nr
        num_species = len(self.list_of_species)
        
        DIFF = {}
        
        start_time = perf_counter()
        for species in self.list_of_species:
            
            D_interp = self.D_interp[species][it]
            c_interp = self.c_interp[species][it]
            
            dt_fact = (2./3.)*self.dt
            
            ############################################
            ############### COMPUTE LHS ################
            ############################################
            lower = np.zeros(self.Nr-1)
            main = np.zeros(self.Nr)
            upper = np.zeros(self.Nr-1)
            
            ## 0<r<a
            rhos = self.rho_grid
            rplus = rhos + drho/2
            rminus = rhos - drho/2
            
            VDplus = Vp(rplus)*D_interp(rplus) / (Vp(rhos)*dr**2)
            VDminus = Vp(rminus)*D_interp(rminus) / (Vp(rhos)*dr**2)
                
            cplus  = c_interp(rhos+drho)*Vp(rhos+drho) / (2*Vp(rhos)*dr)
            cminus = c_interp(rhos-drho)*Vp(rhos-drho) / (2*Vp(rhos)*dr)
            
            main[1:] = 1.0 + dt_fact*(VDplus[1:] + VDminus[1:])
            upper = dt_fact*(-VDplus[:-1] + cplus[:-1])
            lower = dt_fact*(-VDminus[1:] - cminus[1:])
            
            ## r=0
            main[0] = 1.0 + dt_fact*( 4*D_interp(0)/dr**2 + 2*c_interp(drho)/dr )
            upper[0] = -4*dt_fact*D_interp(0)/dr**2
            
            DIFF[species] = diags([lower, main, upper], offsets=[-1, 0, 1], format="csr")    
        
        DIFF_list = [DIFF[species] for species in self.list_of_species]

        # Construct the block diagonal sparse matrix
        LHS = block_diag(DIFF_list, format="csr")
        
        end_time = perf_counter()
        self.time_diffusion_matrix_build += end_time-start_time
            
        # Add implicit terms from sources
        sources_implicit = np.zeros((Nr*num_species,Nr*num_species))
        
        start_time = perf_counter()
        if 'Coll_Heat_Exchange' in self.sources['electrons']:
            sources_implicit -= dt_fact*self.get_collisionalHeatExchange()
        end_time = perf_counter()
        self.time_get_collisionalHeatExchange += end_time-start_time
            
        sources_implicit = csr_matrix(sources_implicit)
        LHS = LHS + sources_implicit
        
        # impose Dirichlet boundary condition
        start_time = perf_counter()
        LHS = LHS.tolil()
        for s in range(1, num_species + 1):  # s starts at 1, up to num_species
            row_idx = s * Nr - 1  # Compute the correct row index

            # Set the entire row to zero
            LHS.rows[row_idx] = []  # Clear all column indices in that row
            LHS.data[row_idx] = []  # Clear all values in that row

            # Set the diagonal element to 1
            LHS[row_idx, row_idx] = 1
        
        LHS = LHS.tocsr()
        
        end_time = perf_counter()
        self.time_set_LHS_boundary_conditions += end_time-start_time
                
        return LHS
        
        
    def get_RHS_vector(self,it):
        
        g_all = []
        for species in self.list_of_species:
        
            g = self.P[species][it-1,:] + (2./3)*self.dt*self.total_sources_explicit[species][it,:]
        
            # apply edge Dirichlet boundary condition
            g[-1] = self.edge_bnd_cnd[species]
            
            g_all.append(g)
            
        g_all = np.concatenate(g_all)
        
        return g_all
    
    def solve_sparse_system(self,matrix,vect):
        
        from scipy.sparse.linalg import spsolve
        import matplotlib.pyplot as plt
        
        # # plot matrix
        # plt.figure(figsize=(6, 6))
        # plt.spy(matrix, markersize=5, color="black")
        # plt.show()
        
        sol = spsolve(matrix,vect)
        
        return sol
                                              
    def make_checks(self):
        
        # check VMEC equilibrium exists
        if(not hasattr(self,'dVdr')):
            print('ERROR: dVdr MUST BE SET!!')
            exit(1)
        
        # check boundary conditions exist for all species
        for species in self.list_of_species:
            if species not in self.edge_bnd_cnd:
                raise KeyError(f"Missing edge boundary condition for species: {species}")
            
        # check initial profiles are set for all species
        for species in self.list_of_species:
            if species not in self.initial_profile:
                raise KeyError(f"Missing initial profile for species: {species}")
            
        # check fluxes info is set
        if(not hasattr(self,'fluxes_info')):
            print('ERROR: set_fluxes must be called before running!!')
            exit(1)
            
    def call_fluxes(self,it):
        
        # compute D_interp and c_interp
        if(self.fluxes_info['type']=='dkespenta'):
            self.call_PENTA3(it)
        elif(self.fluxes_info['type']=='diffusive'):
            self.compute_diffusive_flux(it)
        elif(self.fluxes_info['type']=='beurskens'):
            self.compute_beurskens_flux(it)
        else:
            print('ERROR: Not available other type of flux...')
            exit(0)
            
    def call_PENTA3(self,it):
        import subprocess
        from concurrent.futures import ProcessPoolExecutor, as_completed
        import functools
        plasma_profiles_extension = 'transp_solver'
        
        # create plasma input
        self.plasma.write_plasma_profiles_to_PENTA3(filename='plasma_profiles_'+plasma_profiles_extension+'.dat')
        self.plasma.write_PENTA_namelist()

        surfaces = self.fluxes_info['dkespenta']['surfaces']
        
        dkes_coeffs_path = self.fluxes_info['dkespenta']['dkes_folder']
        
        time_sec = []
        with ProcessPoolExecutor() as executor:
            futures = [executor.submit(process_surfaces, surface, self.wout_path, dkes_coeffs_path) for surface in surfaces]

            for future in as_completed(futures):
                elapsed_seconds = future.result()
                time_sec.append(elapsed_seconds)
                # print(f'Surface processed in {elapsed_seconds:.2f} seconds')
            
        # delete files not needed
        remove = 'rm ucontra* sigmas* flows_vs_Er* plasma_profiles*'
        subprocess.run(remove, shell=True, check=True, text=True, capture_output=True)
        
        # merge _surface_# files into single file
        merge_and_delete('fluxes_vs_roa_surface*','fluxes_vs_roa')
        merge_and_delete('fluxes_vs_Er_surface*','fluxes_vs_Er')
        merge_and_delete('flows_vs_roa_surface*','flows_vs_roa')
        merge_and_delete('Jprl_vs_roa_surface*','Jprl_vs_roa')
            
        self.get_PENTA3_results(it)
        
    def get_PENTA3_results(self,it):
        # creates interpolating functions for Er, Gamma_r[species], Q[species], D_interp[species] and c_interp[species]
        
        import sys
        import subprocess
        from scipy.interpolate import CubicSpline
        sys.path.insert(1,'/home/antonio/STELLOPT/pySTEL/libstell')
        from penta import PENTA
        
        theta = self.fluxes_info['dkespenta']['theta']

        PENTA_class = PENTA(folder_path='.', plasma=self.plasma, lverb=False)
        
        # Create interpolating functions for Er, Gamma and QoT
        self.Er_interp = CubicSpline(PENTA_class.roa_unique,PENTA_class.Er_Maxw)
        
        self.Gamma_interp = {}
        
        for species in self.list_of_species:
            
            rho_extended = np.concatenate([[0.0],PENTA_class.roa_unique])
            
            ## Include Gamma(r=0) = 0
            Gamma_extended = np.concatenate(([0.0],PENTA_class.Gamma_Maxw[species]))
            self.Gamma_interp[species] = CubicSpline(rho_extended,Gamma_extended,extrapolate=True,bc_type='natural')
            
            Q = PENTA_class.QoT_Maxw[species] * self.plasma.get_temperature(species,PENTA_class.roa_unique) * EC # [Q] = J/(m^2*s)
            
            dpdr = self.plasma.get_pressure_der(species,PENTA_class.roa_unique) / self.aminor
            press = self.plasma.get_pressure(species,PENTA_class.roa_unique)
            press_axis = self.plasma.get_pressure(species,0.0)
            dr_large = PENTA_class.roa_unique[0]*self.aminor
            press_dr_large = self.plasma.get_pressure(species,PENTA_class.roa_unique[0])
            
            D_penta = np.where(dpdr!=0, 
                            -theta * Q / dpdr,
                            0.0)
            
            # COMPUTES D_axis assuming Q and dp/dr are zero on the axis (this is a formula resulting from Cauchy rule!)
            if(np.abs(press_dr_large-press_axis) > 1E-14):
                D_axis = -theta*0.5* (Q[0]/dr_large) / ((press_dr_large-press_axis)/dr_large**2)
            else:
                # linear interpolation
                print(' !!! ENTERING in linear interpolation')
                D_axis =  D_penta[0] - PENTA_class.roa_unique[0]*(D_penta[1]-D_penta[0])/(PENTA_class.roa_unique[1]-PENTA_class.roa_unique[0])
            
            D_extended = np.concatenate([[D_axis],D_penta])
            
            self.D_interp[species][it] = CubicSpline(rho_extended,D_extended,extrapolate=True,bc_type='natural')
            
            c_penta = (1-theta) * Q / press
            c_axis = 0.0
            c_extended = np.concatenate([[c_axis],c_penta])
            
            self.c_interp[species][it] = CubicSpline(rho_extended,c_extended)
            
            # This is for bookeeping (it's not used in the calculations)
            Q_extended = np.concatenate(([0.0],Q)) # Include Q(r=0) = 0
            self.Q_interp[species][it] = CubicSpline(rho_extended,Q_extended,extrapolate=True,bc_type='natural')
        
        ## rename file names for bookeeping
        it_subiter = f'{it:03}'
        rename = 'mv fluxes_vs_roa fluxes_vs_roa_'+it_subiter
        subprocess.run(rename, shell=True, check=True, text=True, capture_output=True)
        rename = 'mv fluxes_vs_Er fluxes_vs_Er_'+it_subiter
        subprocess.run(rename, shell=True, check=True, text=True, capture_output=True)
        rename = 'mv flows_vs_roa flows_vs_roa_'+it_subiter
        subprocess.run(rename, shell=True, check=True, text=True, capture_output=True)
        rename = 'mv Jprl_vs_roa Jprl_vs_roa_'+it_subiter
        subprocess.run(rename, shell=True, check=True, text=True, capture_output=True)
            
    def compute_diffusive_flux(self,it):
        # computes an interpolating function for Q=-n*chi*dT/dr -T*Dn*dn/dr
        
        from scipy.interpolate import CubicSpline
        import matplotlib.pyplot as plt
        
        r_grid = self.r_grid
        
        for species in self.list_of_species:
            
            chi = self.fluxes_info['diffusive']['chi']
            
            p_r = CubicSpline(r_grid,self.P[species][it,:])
            dpdr = p_r.derivative()
            
            n_r = self.plasma.get_density(species,self.rho_grid)
            dndr = self.plasma.get_density_der(species,self.rho_grid) / self.aminor
            
            D = chi * np.ones(self.Nr)

            self.D_interp[species][it] = CubicSpline(self.rho_grid,D,bc_type='natural',extrapolate=True)
            
            c = (chi/n_r)*dndr
            c[0] = 0.0
            
            self.c_interp[species][it] = CubicSpline(self.rho_grid,c,bc_type='natural',extrapolate=True)
            
            # this is for bookeeping
            Q = -chi * dpdr(r_grid) + p_r(r_grid)*(chi/n_r)*dndr
            self.Q_interp[species][it] = CubicSpline(self.rho_grid,Q)
            
    def compute_beurskens_flux(self,it):
        # uses model in [ref...]
        from scipy.interpolate import CubicSpline, UnivariateSpline, PchipInterpolator
        import matplotlib.pyplot as plt
        from scipy.signal import savgol_filter
        
        chi_base = self.fluxes_info['beurskens']['chi_base']
        chi_electrons = self.fluxes_info['beurskens']['chi_electrons']
        aLT_critical = self.fluxes_info['beurskens']['aLT_critical']
        alpha = self.fluxes_info['beurskens']['alpha']
        stiffness = self.fluxes_info['beurskens']['stiffness']
        
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
            
            B = self.B ## currently, this is only defined when using a VMEC equilibrium
            chi_gB = self.plasma.get_gyroBohm_diffusivity(ion,B,self.aminor,self.rho_grid)
            
            chi_turb = chi_gB * chi_turb
            
            chi[ion] = chi_base + chi_turb
        
        for species in self.list_of_species:
            
            p_r = CubicSpline(r_grid,self.P[species][it,:])
            dpdr = p_r.derivative()
            
            n_r = self.plasma.get_density(species,self.rho_grid)
            dndr = self.plasma.get_density_der(species,self.rho_grid) / self.aminor
            
            D = chi[species]
            
            # save D of ALL subiter
            self.D_keep[species][it].append(np.array(D))
            
            # average to smooth-out eventual oscillations
            D_avg = np.mean(np.array(self.D_keep[species][it]), axis=0)
            D = D_avg

            self.D_interp[species][it] = CubicSpline(self.rho_grid,D,bc_type='natural',extrapolate=True)
            
            c = (chi[species]/n_r)*dndr
            c[0] = 0.0
            
            # should we average 'c' ??
            
            self.c_interp[species][it] = CubicSpline(self.rho_grid,c,bc_type='natural',extrapolate=True)
            
            # this is for bookeeping
            Q = -chi[species] * dpdr(r_grid) + p_r(r_grid)*(chi[species]/n_r)*dndr
            self.Q_interp[species][it] = CubicSpline(self.rho_grid,Q)
              
    def get_collisionalHeatExchange_non_optimized(self):
        # returns collisional heat exchange to use as implicit operator
        
        from collisions import COLLISIONS
        from scipy import sparse
        
        coll = COLLISIONS()
        
        num_species = len(self.list_of_species)
        clog = np.zeros(self.Nr)
        
        W_s1_s2 = np.zeros((num_species,self.Nr,num_species,self.Nr))
        aux_B = np.zeros((num_species,self.Nr,num_species,self.Nr))

        for ir1,r1 in enumerate(self.rho_grid):
            for ir2,r2 in enumerate(self.rho_grid):
                if(ir1 != ir2):
                    W_s1_s2[is1,ir1,is2,ir2] = 0.0
                else:
                    for is1,species1 in enumerate(self.list_of_species):

                        m1 = self.plasma.mass[species1]
                        Z1 = self.plasma.Zcharge[species1]
                        n1 = self.plasma.get_density(species1,r1)
                        T1 = self.plasma.get_temperature(species1,r1)
                        
                        for is2,species2 in enumerate(self.list_of_species):
                            
                            m2 = self.plasma.mass[species2]
                            Z2 = self.plasma.Zcharge[species2]
                            n2 = self.plasma.get_density(species2,r2)
                            T2 = self.plasma.get_temperature(species2,r2)
                            
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

                            W_s1_s2[is1,ir1,is2,ir2] = gamma*n1  
                            # W_s1_s2[is1,ir1,is2,ir2] = 1.29
                            
                            aux_B[is1,ir1,is2,ir2] = gamma*n2
                            # aux_B[is1,ir1,is2,ir2] = 1.29

        # Add aux_B matrix
        for is1,_ in enumerate(self.list_of_species):
            for ir in range(self.Nr):
                W_s1_s2[is1,ir,is1,ir] -= np.sum(aux_B[is1,ir,:,ir])                 
        
        W_out = np.zeros((num_species*self.Nr,num_species*self.Nr))
        
        j=0
        for is1 in range(num_species):
            for ir1 in range(self.Nr):
                p=0
                for is2 in range(num_species):
                    for ir2 in range(self.Nr):
                        W_out[j,p] = W_s1_s2[is1,ir1,is2,ir2]
                        p=p+1
                j = j+1
                
        W_out = sparse.csr_matrix(W_out)

        return W_out
    
    def get_collisionalHeatExchange(self):
        # returns collisional heat exchange to use as implicit operator
        
        from collisions import COLLISIONS
        from scipy import sparse
        
        coll = COLLISIONS()
        
        num_species = len(self.list_of_species)
        Nr = self.Nr
        
        clog = np.zeros(self.Nr)
        
        W_s1_s2 = np.zeros((num_species,num_species,self.Nr))
        aux_B = np.zeros((num_species,num_species,self.Nr))

        rho_grid = self.rho_grid
        
    #     N1 = np.zeros((num_species,num_species,Nr))
    #     T1 = np.zeros((num_species,num_species,Nr))
    #     M1 = np.zeros((num_species,num_species,Nr))
    #     Z1 = np.zeros((num_species,num_species,Nr))
        
    #     for is1,species1 in enumerate(self.list_of_species):
    #         N1[is1,:,:] = self.plasma.get_density(species1,rho_grid)
    #         T1[is1,:,:] = self.plasma.get_temperature(species1,rho_grid)
    #         M1[is1,:,:] = self.plasma.mass[species1]
    #         Z1[is1,:,:] = self.plasma.Zcharge[species1]
            
    #     N2 = np.transpose(N1,axes=(1,0,2))
    #     T2 = np.transpose(T1,axes=(1,0,2))
    #     M2 = np.transpose(M1,axes=(1,0,2))
    #     Z2 = np.transpose(Z1,axes=(1,0,2))
        
    #    # Create a 3D matrix for 'clog' with the same shape as Z1 and Z2
    #     clog = np.zeros_like(Z1, dtype=float)  # This creates an array with the same shape as Z1, initialized to zero.

    #     # Define conditions
    #     condition_1 = (Z1 > 0) & (Z2 > 0)  # Z1 > 0 and Z2 > 0
    #     condition_2 = (Z1 > 0) & (Z2 < 0)  # Z1 > 0 and Z2 < 0
    #     condition_3 = (Z1 < 0) & (Z2 > 0)  # Z1 < 0 and Z2 > 0

    #     # Apply conditions using np.where and vectorize the operations

    #     clog = np.where(
    #         condition_1,
    #         coll.coullog_ii(M1, Z1, Z1, T1, M2, Z2, N2, T2),  # When Z1 > 0 and Z2 > 0
    #         np.where(
    #             condition_2,
    #             coll.coullog_ei(N2, T2, M1, Z1, N1, T1),  # When Z1 > 0 and Z2 < 0
    #             np.where(
    #                 condition_3,
    #                 coll.coullog_ei(N1, T1, N2, Z2, N2, T2),  # When Z1 < 0 and Z2 > 0
    #                 0.0  # Else case
    #             )
    #         )
    #     )
        
    #     const = (8/np.sqrt(np.pi))*(Z1*Z2*EC*EC)**2 * clog / (8*np.pi*EPS0**2)

    #     vth_s1_sqr = 2*EC*T1/M1
    #     vth_s2_sqr = 2*EC*T2/M2
                
    #     den = M1 * M2 * (vth_s1_sqr + vth_s2_sqr)**1.5
                
    #     gamma = const / den

    #     W_s1_s2 = gamma*N1  
    #     aux_B = gamma*N2
        
        for is1,species1 in enumerate(self.list_of_species):

            m1 = self.plasma.mass[species1]
            Z1 = self.plasma.Zcharge[species1]
            n1 = self.plasma.get_density(species1,rho_grid)
            T1 = self.plasma.get_temperature(species1,rho_grid)
            
            for is2,species2 in enumerate(self.list_of_species):
                
                m2 = self.plasma.mass[species2]
                Z2 = self.plasma.Zcharge[species2]
                n2 = self.plasma.get_density(species2,rho_grid)
                T2 = self.plasma.get_temperature(species2,rho_grid)
                
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
                # W_s1_s2[is1,ir1,is2,ir2] = 1.29
                
                aux_B[is1,is2,:] = gamma*n2
                # aux_B[is1,ir1,is2,ir2] = 1.29

        # Add aux_B matrix
        # for is1,_ in enumerate(self.list_of_species):
        #     for ir in range(self.Nr):
        #         W_s1_s2[is1,is1,ir] -= np.sum(aux_B[is1,:,ir])
        W_s1_s2[np.arange(num_species), np.arange(num_species), :] -= np.sum(aux_B, axis=1)

        
        # W_out = W_s1_s2.reshape(num_species * self.Nr, num_species * self.Nr)         
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
    
    
    def get_alpha_heating(self):
        
        from libstell.fusion import FUSION
        import scipy.sparse as sp
        
        fusion = FUSION()
        E_alpha = fusion.E_DT_He
        
        # make sure deuterium and tritium are the 2nd and 3rd species, otherwise returned matrix is wrong
        if( self.list_of_species[1] != 'deuterium' or self.list_of_species[2] != 'tritium'):
            print('ERROR" deuterium and tritium MUST BE the 2nd and 3rd species, respectively')
            exit(0)
            
        Ns = len(self.list_of_species)
        
        D_diagonal = np.zeros((Ns,self.Nr))
        T_diagonal = np.zeros((Ns,self.Nr))
        
        for iss,species in enumerate(self.list_of_species):
            
            fraction_alpha_heating = self.sources[species]['alpha_heating']['fraction_alpha_heating']
            
            for ir,rho in enumerate(self.rho_grid):

                nD = self.plasma.get_density('deuterium', rho)
                nT = self.plasma.get_density('tritium', rho)
            
                Ti = 0.5* ( self.plasma.get_temperature('deuterium', rho) + self.plasma.get_temperature('tritium', rho) )
            
                sigmav_prime = self.my_sigmaBH(Ti)
                
                D_diagonal[iss,ir] = fraction_alpha_heating * (E_alpha/EC) * sigmav_prime * nT * 0.5
                T_diagonal[iss,ir] = fraction_alpha_heating * (E_alpha/EC) * sigmav_prime * nD * 0.5
        
        # assemble everything in a sparse matrix
        blocks = blocks = [[sp.csr_matrix((self.Nr, self.Nr)) for _ in range(Ns)] for _ in range(Ns)]
        
        for iss in range(Ns):
            diag_D = sp.diags(D_diagonal[iss, :])  # Create diagonal matrix from D_diagonal[i,:]
            diag_T = sp.diags(T_diagonal[iss, :])  # Create diagonal matrix from T_diagonal[i,:]
            
            blocks[iss][1] = diag_D  # Place D_diagonal matrix in 2nd block column
            blocks[iss][2] = diag_T  # Place T_diagonal matrix in 3rd block column
      
        S_alpha = sp.bmat(blocks, format="csr")
        
        #check S_alpha has the expected dimensions of (NsxNr,NsxNr)
        expected_shape = (Ns * self.Nr, Ns * self.Nr)
        assert S_alpha.shape == expected_shape, f"Error: Matrix shape is {S_alpha.shape}, while expected is {expected_shape}"
        
        # import matplotlib.pyplot as plt
        
        # # plot matrix
        # plt.figure(figsize=(6, 6))
        # plt.spy(S_alpha, markersize=5, color="black")
        # plt.show()
        
        return S_alpha
        
    def my_sigmaBH(self, ti_eV):
        # adappted from sigmaBH in fusion class
        
        from libstell.fusion import FUSION
        
        fusion = FUSION()
        
        reaction='DT'
        C = fusion.C_DICT[reaction]
        BG = fusion.BG_DICT[reaction]
        MRC2 = fusion.MRC2_DICT[reaction]
        
        ti_kev = ti_eV * 1E-3
        zeta   = ( ( ( C[5] * ti_kev ) + C[3] ) * ti_kev + C[1] ) * ti_kev
        zeta   = zeta / ( ( ( ( C[6] * ti_kev ) + C[4] ) * ti_kev + C[2] ) * ti_kev + 1.0 )
        zeta   = 1.0 - zeta
        theta  = ti_kev / zeta
        eta    = ( 0.25 * BG * BG / theta ) ** (1.0/3.0)
        
        result = 1.0E-6 * C[0] * theta * np.sqrt( eta / ( MRC2 * ti_kev * ti_kev * ti_kev ) ) * np.exp( -3 * eta )
        
        return 1E-3*result / ti_kev
    
def process_surfaces(surface,wout_path,dkes_coeffs_path):
    import time
    import subprocess
    
    start_time = time.time()
    
    type_of_write = 0
    Er_min_V_cm = -300
    Er_max_V_cm = 300

    # wout_path = solver_class.wout_path
    EparB = 0.0

    #Sonine (Laguerre) polynomials
    Smax = 1
    
    plasma_profiles_extension = 'transp_solver'
    
    # copy coeffs files
    subprocess.run(f'cp {dkes_coeffs_path}/*surface_{surface} .',shell=True, check=True, text=True, capture_output=True)
    
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
        
        
        
            
         