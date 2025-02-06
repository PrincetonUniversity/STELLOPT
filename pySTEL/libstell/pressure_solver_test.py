"""
This library provides a python class for solving pressure transport
equations
"""

import numpy as np
import sys

# Constants
EC = 1.602176634E-19 # Electron charge [C]
EPS0 = 8.8541878188E-12 # Vacuum permittivity [F/m]

# PENTA Class
class PRESSURE_SOLVER_test:
    
    def __init__(self, plasma_class):
        
        from collections import defaultdict
        
        self.list_of_species = plasma_class.list_of_species
        
        # Check density profiles exist for all species
        for species in self.list_of_species:
            if(species not in plasma_class.density):
                print('ERROR" density of {species} has not been set yet')
                exit(0)
                
        self.plasma = plasma_class
                
        # Initialize dicitionaries
        self.edge_bnd_cnd = {}
        self.initial_profile = {}
        self.sources = defaultdict(lambda: defaultdict(dict))
        self.fluxes_info = defaultdict(lambda: defaultdict(dict))
                
        print(f'Solvers for pressure of {self.list_of_species} INITIALIZED!')
                
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
            case 'cylindrical':
                if(aminor is None or Rmajor is None):
                    print('ERROR: For a cylindrical equilibrium, Rmajor and aminor must be given')
                    exit(0)
                else:
                    self.aminor=aminor
                    self.dVdr = lambda rho: 4*np.pi*np.pi*Rmajor*aminor  * rho
            
    def set_source(self,species,source_type, total_power=None, sigma_rho=None, fraction_alpha_heating=None, cte_source=None, interpolant_2D=None, time_dependent_factor=None, tau_alphas=None, tau_palphas=None):
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
                    print('ERROR: tau_alphas and tau_alphas must be given (in seconds)')
                    exit(0)
                else:
                    self.sources[species][source_type] = {'tau_alphas': tau_alphas, 'tau_palphas' : tau_palphas}
            case 'external_gaussian':
                if((total_power is None) or (sigma_rho is None)):
                    print('ERROR: Need to provide total_power [W] and sigma_rho for gaussian external source')
                    exit(1) 
                else:
                    self.sources[species][source_type] = {'total_power' : total_power, 'sigma_rho' : sigma_rho }
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
            case _:
                print(f'ERROR: Source type {source_type} is NOT possible')
                exit(0)
                
    def set_fluxes(self,type: str,dkes_folder=None,surfaces=None,D_coeff=None):
        # sets type of fluxes
        # OPTION1: type='dkespenta'; dkes_folder and surfaces(list of integers) must be provided
        # OPTION2: type='diffusive'; D_coeff must be provided (diffusion coefficient)
            
        match type:
            case 'dkespenta':
                #checks that dkes_folder and surfaces are provided
                if((dkes_folder is None) or (surfaces is None)):
                    print('ERROR: dkes_folder and surfaces must be provided!!')
                    exit(0)
                    
                self.fluxes_info['type'] = type
                self.fluxes_info[type]['dkes_folder'] = dkes_folder
                self.fluxes_info[type]['surfaces'] = surfaces
                
            case 'diffusive':
                #checks that D_coeff is provided
                if(D_coeff is None):
                    print('ERROR: D_coeff must be provided!')
                    exit(0)
                self.fluxes_info['type'] = type
                self.fluxes_info[type]['D_coeff'] = D_coeff
                
    def set_temperature(self,species,rho,density,pressure):
        # from density (m^-3) and pressure (Pa), sets temperature (eV) in plasma class
        
        temperature =  pressure / (EC*density) # eV
        
        self.plasma.set_temperature(species,'interp',rho_vals=rho,T_vals=temperature) 
                 
    def run(self,Nr,dt,tstart,tend,theta=1.0,tolerance=1E-2,max_subiter=12):
        
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
        self.Nr = Nr
        self.theta = theta
            
        print(' ')
        print( ' ***********************')
        print(f' *  tstart = {tstart:5.2f}s    *')
        print(f' *  tend   = {tend:5.2f}s    *')
        print(f' *  dt     = {dt:5.2f}s    *')
        print(f' *  Nt     = {Nt:3}       *')
        print(f' *  drho   = {drho:5.2f}     *')
        print( ' ***********************')
            
        ################## INITIALIZE DICTIONARIES ########################################
        press = {}
        dens = {}
        self.P = {}
        self.T = {}
        self.Q = {}
        delta_p = {}
        p_old = {}
        self.D_interp = {}
        self.total_sources_explicit = {}
        self.all_sources = {}
        for species in self.list_of_species:
            press[species] = self.initial_profile[species](rho)
            dens[species] = self.plasma.get_density(species,rho)
            self.set_temperature(species,rho,dens[species],press[species])
            
            self.P[species] = np.zeros((Nt,Nr))
            self.P[species][0,:] = press[species]
            
            self.T[species] = np.zeros((Nt,Nr))
            self.T[species][0,:] = self.plasma.get_temperature(species,rho)
            
            self.Q[species] = np.zeros((Nt,Nr))
            
            self.D_interp[species] = [None]*Nt
                        
            self.total_sources_explicit[species] = np.zeros((Nt,Nr))
            
            p_old[species] = 1E3*np.ones(Nr) # so on loop 1 we don't divide by zero in delta_p
            
            # Initialize arrays for each source_type in the specified species
            self.all_sources[species] = {}
            for source_type in self.sources[species].keys():
                self.all_sources[species][source_type] = np.zeros((Nt, Nr))
                
        if('Bremsstrahlung_alphas' in self.sources['electrons']):
            self.Nalphas_fast = np.zeros((Nt,Nr))
            self.Nalphas_thermal = np.zeros((Nt,Nr))
            
        ####################################################################################
                
        print(' ')
        header_str = '  TIME [s]     NSUB      TE_AXIS [keV]     SE_AXIS [MW/m^3]    TI1_AXIS [keV]    SI1_AXIS [MW/m^3]    MAX(dp/p_old)'  
        print(header_str)
        print('  '+'='*len(header_str))
        
        # Compute sources at t=0.0 and print info
        self.update_at_start()
        
        ### LOOP IN TIME STARTING AT t=tstart+dt ###
        for it,t in enumerate(time[1:],start=1):

            ### SUBCYCLE
            delta_p_all = 10*tolerance
            subiter=1
            while(delta_p_all > tolerance and subiter<max_subiter):
            
                # sets temperature in plasma class from density and pressure for ALL species
                for species in self.list_of_species:
                    self.set_temperature(species,rho,dens[species],press[species])
                    
                    # l-1
                    self.P[species][it,:] = press[species]
                    self.T[species][it,:] = self.plasma.get_temperature(species,rho) # this is only for bookeeping
                    
                # if(t>15 and subiter==1):
                #     return
                    
                # call PENTA3 and compute interpolating functions self.Er_interp, self.Gamma_interp[species] and self.Q_interp[species]
                if(self.fluxes_info['type']=='dkespenta'):
                    self.call_PENTA3()
                elif(self.fluxes_info['type']=='diffusive'):
                    self.compute_diffusive_flux(it)
                else:
                    print('ERROR: Not available other type of flux...')
                    exit(0)
                
                #solver for each species
                for species in self.list_of_species:
                    
                    # # l-1
                    # self.P[species][it,:] = press[species] 
                    
                    # compute sources on grid (1D-array)
                    self.total_sources_explicit[species][it,:] = self.get_sources_explicit(species,rho,it)  # W/m^3
                    
                    RHS_vector = self.get_RHS_vector(species,it)
                    
                    LHS_matrix = self.get_LHS_tridig_matrix(species,it)
                    
                    # solve system
                    press[species] = self.solve_tridiagonal_system(LHS_matrix,RHS_vector)
                    
                    # self.check_NaNs_and_neg_values(species,it)
                    # TEMP: TREAT NEG VALUES
                    # eps = 1E-10
                    # press[species] = np.where(press[species] < 0, eps, press[species])

                    delta_p[species] = np.max( np.where( p_old[species]>1E-10, np.abs((press[species]-p_old[species])/p_old[species]), 0 ) )
                    
                    p_old[species] = press[species]
                    # self.P[species][it,:] = press[species] ---> I moved this to the beginning of the loop
                    
                delta_p_all = np.max([np.max(value) for value in delta_p.values()])
                
                info_str = f'  {t:<13.2f}{subiter:<10}{self.T['electrons'][it,0]/1E3:<18.3f}{self.total_sources_explicit['electrons'][it,0]/1E6:<20.2E}{self.T['deuterium'][it,0]/1E3:<18.3f}{self.total_sources_explicit['deuterium'][it,0]/1E6:<21.2E}{delta_p_all:<13.2E}'
                print(info_str)
                
                subiter += 1
                        
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
                     
            match source_type:
                case 'Bremsstrahlung':
                    aux_source = np.zeros(len(rho_grid))
                    for ir,rho in enumerate(rho_grid):
                        for ion in self.plasma.ion_species:
                            zi = self.plasma.Zcharge[ion]
                            ni = self.plasma.get_density(ion,rho)
                            ne = self.plasma.get_density('electrons',rho)
                            Te = self.plasma.get_temperature('electrons',rho)
                            
                            aux_source[ir] -= fusion.BremsstrahlungPower(zi,ni,ne,Te)
                            
                    #save in dictionary for bookeeping
                    self.all_sources[species][source_type][it,:] = aux_source
                    
                case 'Bremsstrahlung_alphas':
                    aux_source = np.zeros(len(rho_grid))
                    for ir,rho in enumerate(rho_grid):
                        ne = self.plasma.get_density('electrons',rho)
                        Te = self.plasma.get_temperature('electrons',rho)
                        Z_alpha = 2
                        tau_alpha = self.sources['electrons']['Bremsstrahlung_alphas']['tau_alphas']
                        tau_palpha = self.sources['electrons']['Bremsstrahlung_alphas']['tau_palphas']
                        nD = self.plasma.get_density('deuterium',rho)
                        nT = self.plasma.get_density('tritium',rho)
                        Ti = 0.5* ( self.plasma.get_temperature('deuterium', rho) + self.plasma.get_temperature('tritium', rho) )
                        sigmav = fusion.sigmaBH(Ti,'DT') # m^3/s
                        
                        # ideally here it should be of current time iteration, previous subiteration... TO DO LATER...
                        # ESSENTIALLY I need to update it w/ previous iter whe going to 1st subiter...
                        self.Nalphas_fast[it,ir] = (self.Nalphas_fast[it-1,ir] + self.dt*nD*nT*sigmav) / (1+self.dt/tau_alpha)
                        self.Nalphas_thermal[it,ir] = (self.Nalphas_thermal[it,ir] + (self.dt/tau_alpha)*self.Nalphas_fast[it,ir]) / (1+self.dt/tau_palpha)
                        
                        aux_source[ir] -= fusion.BremsstrahlungPower(Z_alpha,self.Nalphas_thermal[it,ir],ne,Te)
                        
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
                    
                    # if(self.T['deuterium'][it-1,0] > 12E3): aux_source = 0.0
                    
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
                    
                    W_explicit,W_implicit = self.get_collisionalHeatExchange(species,it)
                    
                    W_bookeeping = np.sum(W_explicit,axis=0) + np.sum(W_implicit,axis=0)*self.P[species][it,:]

                    sum_W_explicit = np.sum(W_explicit,axis=0)
                            
                    aux_source = sum_W_explicit
                    
                    #save in dictionary for bookeeping
                    self.all_sources[species][source_type][it,:] = W_bookeeping
                            
                case 'alpha_heating':
                    nD = self.plasma.get_density('deuterium', rho_grid)
                    nT = self.plasma.get_density('tritium', rho_grid)
                    
                    Ti = 0.5* ( self.plasma.get_temperature('deuterium', rho_grid) + self.plasma.get_temperature('tritium', rho_grid) )
                    
                    sigmav = [fusion.sigmaBH(ti,'DT') for ti in Ti]
                    
                    S_alpha = nD * nT * sigmav *  fusion.E_DT_He # W/m^3
                    
                    fraction_alpha_heating = self.sources[species]['alpha_heating']['fraction_alpha_heating']
                    
                    aux_source = S_alpha * fraction_alpha_heating
                    
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

                case _:
                    print(f'ERROR: Source type {source_type} not defined....')
                    exit(0)
                    
            total_source_explicit += aux_source
                    
                    
        return total_source_explicit
    
    def get_LHS_tridig_matrix(self,species,it):
        # returns a 3xNr array containing diagonal,lower and upper arrays of LHS matrix
        
        from scipy.interpolate import CubicSpline
        from scipy.linalg import eigvals
        
        drho = self.rho_grid[1]-self.rho_grid[0]
        dr = self.aminor * drho
        Vp = self.dVdr
        a = self.aminor
        
        Q = self.Q_interp[species]
        
        self.Q[species][it,:] = Q(self.rho_grid)  # this is for bookeeping
        
        # get interpolating function for dpdr (at l-1)
        p_interp_over_a = CubicSpline(self.rho_grid,self.P[species][it,:]/a) #divides by aminor to go from rho to r
        dpdr = p_interp_over_a.derivative()
        
        # compute D_interp (interpolating function)
        D = np.zeros(self.Nr)
        
        # dP = self.P[species][it,1]-self.P[species][it,0]
        # if(np.abs(dP)<1E-15):
        #     D[0] = 0.0
        # else:
        #     D[0] = -self.theta * Q(drho) * dr * 0.5  / dP
        
        D[1:] = np.where(dpdr(self.rho_grid[1:])!=0, 
                        -self.theta * Q(self.rho_grid[1:]) / dpdr(self.rho_grid[1:]),
                        0.0)
        D[0] = 2*D[1] - D[2]
        D_interp = CubicSpline(self.rho_grid,D)
        # D_interp = lambda rho: 1.5 #CubicSpline(self.rho_grid,D)
        # D[0] = 1.5
        
        #bookeping
        self.D_interp[species][it] = D_interp
        # print(self.D_interp[it])
        
        # compute c_interp (interpolating function)
        c = np.zeros(self.Nr)
        c[0] = 0
        c[1:] = (1-self.theta) * Q(self.rho_grid[1:]) / (p_interp_over_a(self.rho_grid[1:])*a)
        c_interp = CubicSpline(self.rho_grid,c)
        
        dt_fact = (2./3.)*self.dt
        
        # CHECK IF COLL HEAT EXCHANGE IS SET; IF YES, THEN SHOULD ADD THE IMPICIT TERM HERE
        sources_implicit_facts = np.zeros(self.Nr)
        if 'Coll_Heat_Exchange' in self.sources[species]:
            _, implicit_heat_exchange = self.get_collisionalHeatExchange(species,it)
            sources_implicit_facts -= dt_fact*np.sum(implicit_heat_exchange,axis=0)
        
        ############################################
        ############### COMPUTE LHS ################
        ############################################
        lower = np.zeros(self.Nr-1)
        main = np.zeros(self.Nr)
        upper = np.zeros(self.Nr-1)
        
        ## r=0
        main[0] = 1.0 + dt_fact*( 4*D[0]/dr**2 + 2*c[1]/dr ) + sources_implicit_facts[0]
        upper[0] = -4*dt_fact*D[0]/dr**2
        
        ## 0<r<a
        for ir,rho in enumerate(self.rho_grid[1:-1],start=1):
            
            rplus = rho + drho/2
            rminus = rho - drho/2
            
            VDplus = Vp(rplus)*D_interp(rplus) / (Vp(rho)*dr**2)
            VDminus = Vp(rminus)*D_interp(rminus) / (Vp(rho)*dr**2)
            
            cplus  = c_interp(rho+drho)*Vp(rho+drho) / (2*Vp(rho)*dr)
            cminus = c_interp(rho-drho)*Vp(rho-drho) / (2*Vp(rho)*dr)
            
            main[ir] = 1.0 + dt_fact*(VDplus+VDminus) + sources_implicit_facts[ir]
            upper[ir] = dt_fact*(-VDplus+cplus)
            lower[ir-1] = dt_fact*(-VDminus-cminus)

        ## r=1
        main[-1] = 1.0
        lower[-1] = 0.0
        
        ####
        LHS = np.zeros((3,self.Nr))
        LHS[0,1:] = upper
        LHS[1,:] = main
        LHS[2,:-1] = lower
        
        
        
        # tridiag_matrix = np.zeros((self.Nr, self.Nr))
        # np.fill_diagonal(tridiag_matrix, main)  # Main diagonal
        # np.fill_diagonal(tridiag_matrix[1:], lower)  # Lower diagonal
        # np.fill_diagonal(tridiag_matrix[:, 1:], upper)  # Upper diagonal

        # # Compute eigenvalues of the full matrix
        # eigenvalues = eigvals(tridiag_matrix)

        # # Compute eigenvalues of the inverse
        # eigenvalues_inverse = 1 / eigenvalues
        
        # print(f'lambdas = {eigenvalues_inverse}')


        
        return LHS
        
    def get_RHS_vector(self,species,it):
        
        g = self.P[species][it-1,:] + (2./3)*self.dt*self.total_sources_explicit[species][it,:]
        
        # apply edge Dirichlet boundary condition
        g[-1] = self.edge_bnd_cnd[species]
        
        return g
    
    def solve_tridiagonal_system(self,matrix,vect):
        
        from scipy.linalg import solve_banded
        
        sol = solve_banded((1,1),matrix,vect)
        
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
            
    def call_PENTA3(self):
        
        import subprocess
        
        plasma_profiles_extension = 'transp_solver'
        
        # create plasma input
        self.plasma.write_plasma_profiles_to_PENTA3(filename='plasma_profiles_'+plasma_profiles_extension+'.dat')
        self.plasma.write_PENTA_namelist()
        
        # Er_min in V/cm
        Er_min_V_cm = -300
        # Er_max in V/cm
        Er_max_V_cm = 300

        wout_path = self.wout_path
        EparB = 0.0

        #Sonine (Laguerre) polynomials
        Smax = 1

        surfaces = self.fluxes_info['dkespenta']['surfaces']

        #loop in surfaces
        for isurface in surfaces:
            
            if(isurface==surfaces[0]):
                type_of_write = 0
            else:
                type_of_write = 1
                
            # copy coeffs files
            dkes_coeffs_path = self.fluxes_info['dkespenta']['dkes_folder']
            subprocess.run(f'cp {dkes_coeffs_path}/*surface_{isurface} .',shell=True, check=True, text=True, capture_output=True)
                
            extension_star_files = f'surface_{isurface}'

            call_penta3 = f'~/bin/xpenta {extension_star_files} {Er_min_V_cm} {Er_max_V_cm} {isurface} {type_of_write} {wout_path} {plasma_profiles_extension} {EparB} {Smax}'
            
            result = subprocess.run(call_penta3, shell=True, check=True, text=True, capture_output=True)
            # print(result.stdout)
            # if(result.stderr):
            #     print(result.stderr)
            
        self.get_PENTA3_results()
        
    def get_PENTA3_results(self):
        # creates interpolating functions for Er, Gamma_r[species] and QoT[species]
        
        import sys
        from scipy.interpolate import CubicSpline
        sys.path.insert(1,'/home/antonio/STELLOPT/pySTEL/libstell')
        from penta import PENTA

        PENTA_class = PENTA(folder_path='.', plasma=self.plasma, lverb=False)
        
        # Create interpolating functions for Er, Gamma and QoT
        self.Er_interp = CubicSpline(PENTA_class.roa_unique,PENTA_class.Er_Maxw)
        
        self.Gamma_interp = {}
        self.Q_interp = {}
        
        for species in self.list_of_species:
            self.Gamma_interp[species] = CubicSpline(PENTA_class.roa_unique,PENTA_class.Gamma_Maxw[species])
            
            # [Q] = J/(m^2*s)
            Q = PENTA_class.QoT_Maxw[species] * self.plasma.get_temperature(species,PENTA_class.roa_unique) * EC
            self.Q_interp[species] = CubicSpline(PENTA_class.roa_unique,Q) # [self.Q] = J/(m^2*s)
            
    def compute_diffusive_flux(self,it):
        # computes an interpolating function for Q=-D*dp/dr
        
        from scipy.interpolate import CubicSpline
        
        D = self.fluxes_info['diffusive']['D_coeff']
        
        self.Q_interp = {}
        
        for species in self.list_of_species:
            DP = CubicSpline(self.rho_grid,-D*self.P[species][it,:]/self.aminor) #divides by aminor to go from rho to r
            self.Q_interp[species] = DP.derivative()
            
    def update_at_start(self):
           
        if(self.fluxes_info['type']=='dkespenta'):
            self.call_PENTA3()
            
        if(self.fluxes_info['type']=='diffusive'):
            self.compute_diffusive_flux(it=0)
            
        for species in self.list_of_species:
            self.total_sources_explicit[species][0,:] = self.get_sources_explicit(species,self.rho_grid,0)  # W/m^3
            self.Q[species][0,:] = self.Q_interp[species](self.rho_grid)
        
        # print info for t=t_starts
        info_str = f'  {self.time[0]:<13.2f}{1:<10}{self.T['electrons'][0,0]/1E3:<18.3f}{self.total_sources_explicit['electrons'][0,0]/1E6:<20.2E}{self.T['deuterium'][0,0]/1E3:<18.3f}{self.total_sources_explicit['deuterium'][0,0]/1E6:<21.2E}{0.0:<13.2E}'
        print(info_str)
        
    def check_NaNs_and_neg_values(self,species,it):
        # checks if self.P[species][:,:] has NaNs or negative values
        # if a NaN is found, program is aborted
        # if negative value is found, a warning is yield and the neg value substituted by a very small number
        
        # look for NaNs
        if np.any(np.isnan(self.P[species][it,:])):
            print('ERROR: Pressure has NaN values !! ')
            exit(1)
                
        # look for neg values
        eps = 1E-10
        self.P[species][it, :] = np.where(self.P[species][it, :] < 0, eps, self.P[species][it, :])
        
    def get_collisionalHeatExchange(self,species,it):
        # returns 2 arrays: the explicit part of W_s1_s2 and the implicit fact of W_s1_s2

        from collisions import COLLISIONS
        
        coll = COLLISIONS()
        
        num_species = len(self.list_of_species)
        clog = np.zeros(self.Nr)
        
        W_s1_s2_explicit = np.zeros((num_species,self.Nr))
        W_s1_s2_implicit = np.zeros((num_species,self.Nr))

        # for is1,species1 in enumerate(self.list_of_species):
        
        species1 = species
        m1 = self.plasma.mass[species1]
        Z1 = self.plasma.Zcharge[species1]
        n1 = self.plasma.get_density(species1,self.rho_grid)
        T1 = self.plasma.get_temperature(species1,self.rho_grid)
        
        for is2,species2 in enumerate(self.list_of_species):
            
            m2 = self.plasma.mass[species2]
            Z2 = self.plasma.Zcharge[species2]
            n2 = self.plasma.get_density(species2,self.rho_grid)
            T2 = self.plasma.get_temperature(species2,self.rho_grid)
            
            # get Coulomb logarithm
            for ir in range(self.Nr):
                if(Z1>0 and Z2>0):
                    clog[ir] = coll.coullog_ii(m1,Z1,n1[ir],T1[ir],m2,Z2,n2[ir],T2[ir])
                elif(Z1>0 and Z2<0):
                    clog[ir] = coll.coullog_ei(n2[ir],T2[ir],m1,Z1,n1[ir],T1[ir])
                elif(Z1<0 and Z2>0):
                    clog[ir] = coll.coullog_ei(n1[ir],T1[ir],m2,Z2,n2[ir],T2[ir])
                else:
                    clog[ir] = 0.0

            gamma = (Z1*Z2*EC*EC)**2 * clog / (8*np.pi*EPS0**2)

            vth_s1_sqr = 2*EC*T1/m1
            vth_s2_sqr = 2*EC*T2/m2

            num_explicit_part = gamma * n1 * n2 * T2
            num_implicit_fact = -gamma * n2
            den = m1 * m2 * (vth_s1_sqr + vth_s2_sqr)**1.5
    
            explicit_term_W_s1_s2 = (8/np.sqrt(np.pi)) * num_explicit_part / den  # eV / (s.m^3)
            explicit_term_W_s1_s2 = explicit_term_W_s1_s2 * EC  # W/m^3
            
            implicit_fact_W_s1_s2 = (8/np.sqrt(np.pi)) * num_implicit_fact / den
            implicit_fact_W_s1_s2 = implicit_fact_W_s1_s2 # *EC -- NEED TO REMOVE *EC from implicit
            
            W_s1_s2_explicit[is2,:] = 1.29*EC*(n2*self.T[species2][it-1,:])   #explicit_term_W_s1_s2
            W_s1_s2_implicit[is2,:] =   -1.29 #implicit_fact_W_s1_s2
            
        # print(f'coll_HEAT={np.sum(W_s1_s2_explicit,axis=0)}')
  
        return W_s1_s2_explicit,W_s1_s2_implicit
        
            
# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)
        
        
        
            
        