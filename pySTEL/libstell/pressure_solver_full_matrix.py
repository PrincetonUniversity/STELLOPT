"""
This library provides a python class for solving pressure transport
equations
"""

import numpy as np
import sys
import matplotlib.pyplot as plt

# Constants
EC = 1.602176634E-19 # Electron charge [C]
EPS0 = 8.8541878188E-12 # Vacuum permittivity [F/m]

# PENTA Class
class PRESSURE_SOLVER_FULL_MATRIX:
    
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
                
    def set_fluxes(self,type: str,dkes_folder=None,surfaces=None,chi=None):
        # sets type of fluxes
        # OPTION1: type='dkespenta'; dkes_folder and surfaces(list of integers) must be provided
        # OPTION2: type='diffusive'; Dn and chi must be provided (partical and heat collisional diffusion coefficients)
            
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
                #checks that diffusion coefficients are provided
                if(chi is None):
                    print('ERROR: chi must be provided!')
                    exit(0)
                self.fluxes_info['type'] = type
                self.fluxes_info[type]['chi'] = chi
                
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
        self.drho = drho
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
            
        ################## INITIALIZE VARIABLES ########################################
        press = []
        dens = {}
        self.P = {}
        self.T = {}
        self.Q = {}
        self.D_interp = {}
        self.c_interp = {}
        self.total_sources_explicit = {}
        self.all_sources = {}
        self.ECRH_off = False
        for species in self.list_of_species:
            p_init = self.initial_profile[species](rho)
            dens[species] = self.plasma.get_density(species,rho)
            self.set_temperature(species,rho,dens[species],p_init)
            
            self.P[species] = np.zeros((Nt,Nr))
            self.P[species][0,:] = p_init
            
            press.append(p_init)
            
            self.T[species] = np.zeros((Nt,Nr))
            self.T[species][0,:] = self.plasma.get_temperature(species,rho)
            
            self.Q[species] = np.zeros((Nt,Nr))
            
            self.D_interp[species] = [None]*Nt
            self.c_interp[species] = [None]*Nt
                        
            self.total_sources_explicit[species] = np.zeros((Nt,Nr))

            # Initialize arrays for each source_type in the specified species
            self.all_sources[species] = {}
            for source_type in self.sources[species].keys():
                self.all_sources[species][source_type] = np.zeros((Nt, Nr))
                
        if('Bremsstrahlung_alphas' in self.sources['electrons']):
            self.Nalphas_fast = np.zeros((Nt,Nr))
            self.Nalphas_thermal = np.zeros((Nt,Nr))
            
        p_old = 1E3*np.ones(Nr*len(self.list_of_species)) # so on loop 1 we don't divide by zero in delta_p
        press = np.concatenate(press)
            
        ####################################################################################
                
        print(' ')
        header_str = '  TIME [s]     NSUB      TE_AXIS [keV]     SE_AXIS [MW/m^3]    TI1_AXIS [keV]    SI1_AXIS [MW/m^3]    MAX(dp/p_old)'  
        print(header_str)
        print('  '+'='*len(header_str))
        
        # Compute sources at t=0.0 and print info
        # self.update_at_start()
        
        ### LOOP IN TIME STARTING AT t=tstart+dt ###
        for it,t in enumerate(time[1:],start=1):
            self.it = it
            
            ### SUBCYCLE
            delta_p = 10*tolerance
            subiter=1
            while(delta_p > tolerance and subiter<max_subiter):
                self.subiter = subiter
                # sets temperature in plasma class from density and pressure for ALL species
                k=0
                for species in self.list_of_species:
                    
                    self.P[species][it,:] = press[k:(k+Nr)]
                    
                    self.set_temperature(species,rho,dens[species],self.P[species][it,:])
                    self.T[species][it,:] = self.plasma.get_temperature(species,rho) # this is only for bookeeping
                    
                    k = k+Nr   

                # compute D_interp and c_interp
                if(self.fluxes_info['type']=='dkespenta'):
                    self.call_PENTA3(it)
                elif(self.fluxes_info['type']=='diffusive'):
                    self.compute_diffusive_flux(it)
                else:
                    print('ERROR: Not available other type of flux...')
                    exit(0)
                
                #solver for each species
                for species in self.list_of_species:
                    self.total_sources_explicit[species][it,:] = self.get_sources_explicit(species,rho,it)  # W/m^3
                                        
                RHS_vector = self.get_RHS_vector(it)
                LHS_matrix = self.get_LHS_matrix(it)

                # solve system
                press = self.solve_sparse_system(LHS_matrix,RHS_vector)
                
                #PICARD FACTOR
                # fpicard = 0.75
                # press = press*fpicard + (1-fpicard)*p_old
            
                delta_p = np.max( np.where( p_old>1E-10, np.abs((press-p_old)/p_old), 0 ) )
                
                p_old = press
                
                info_str = f'  {t:<13.2f}{subiter:<10}{self.T['electrons'][it,0]/1E3:<18.3f}{self.total_sources_explicit['electrons'][it,0]/1E6:<20.2E}{self.T['deuterium'][it,0]/1E3:<18.3f}{self.total_sources_explicit['deuterium'][it,0]/1E6:<21.2E}{delta_p:<13.2E}'
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
    
    def get_LHS_matrix(self,it):
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

            ## r=1
            main[-1] = 1.0
            lower[-1] = 0.0
            
            ####
            DIFF[species] = diags([lower, main, upper], offsets=[-1, 0, 1], format="csr")    
        
        DIFF_list = [DIFF[species] for species in self.list_of_species]

        # Construct the block diagonal sparse matrix
        LHS = block_diag(DIFF_list, format="csr")
            
        # Add implicit terms from sources
        sources_implicit = np.zeros((Nr*num_species,Nr*num_species))
        
        if 'Coll_Heat_Exchange' in self.sources['electrons']:
            sources_implicit -= dt_fact*self.get_collisionalHeatExchange()
            
        sources_implicit = csr_matrix(sources_implicit)
        LHS = LHS + sources_implicit
                
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

        print(f'Total time: {np.max(time_sec):.1f}s')
            
        # delete files not needed
        remove = 'rm ucontra* sigmas* plasma_profiles_check*'
        result = subprocess.run(remove, shell=True, check=True, text=True, capture_output=True)
        
        #
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

        PENTA_class = PENTA(folder_path='.', plasma=self.plasma, lverb=False)
        
        # Create interpolating functions for Er, Gamma and QoT
        self.Er_interp = CubicSpline(PENTA_class.roa_unique,PENTA_class.Er_Maxw)
        
        self.Gamma_interp = {}
        self.Q_interp = {}
        
        for species in self.list_of_species:
            
            rho_extended = np.concatenate([[0.0],PENTA_class.roa_unique])
            
            ## Include Gamma(r=0) = 0
            Gamma_extended = np.concatenate(([0.0],PENTA_class.Gamma_Maxw[species]))
            self.Gamma_interp[species] = CubicSpline(rho_extended,Gamma_extended,extrapolate=True,bc_type='natural')
            
            Q = PENTA_class.QoT_Maxw[species] * self.plasma.get_temperature(species,PENTA_class.roa_unique) * EC # [Q] = J/(m^2*s)
            
            # dndr_penta = self.plasma.get_density_der(species,PENTA_class.roa_unique) / self.aminor
            # dTdr_penta = self.plasma.get_temperature_der(species,PENTA_class.roa_unique) / self.aminor
            # n_penta = self.plasma.get_density(species,PENTA_class.roa_unique)
            # T_penta = self.plasma.get_temperature(species,PENTA_class.roa_unique)
            # p_penta_axis = EC * self.plasma.get_density(species,0.0)*self.plasma.get_temperature(species,0.0)
            # p_penta_dr =   EC * self.plasma.get_density(species,PENTA_class.roa_unique[0])*self.plasma.get_temperature(species,PENTA_class.roa_unique[0])
            # dpdr_penta = EC * (n_penta*dTdr_penta + T_penta*dndr_penta)
            
            dpdr = self.plasma.get_pressure_der(species,PENTA_class.roa_unique) / self.aminor
            press = self.plasma.get_pressure(species,PENTA_class.roa_unique)
            press_axis = self.plasma.get_pressure(species,0.0)
            dr_large = PENTA_class.roa_unique[0]*self.aminor
            press_dr_large = self.plasma.get_pressure(species,PENTA_class.roa_unique[0])
            
            D_penta = np.where(dpdr!=0, 
                            -self.theta * Q / dpdr,
                            0.0)
            
            # COMPUTES D_axis assuming Q and dp/dr are zero on the axis (this is a formula resulting from Cauchy rule!)
            if(np.abs(press_dr_large-press_axis) > 1E-14):
                D_axis = -self.theta*0.5* (Q[0]/dr_large) / ((press_dr_large-press_axis)/dr_large**2)
            else:
                # linear interpolation
                print(' !!! ENTERING in linear interpolation')
                D_axis =  D_penta[0] - PENTA_class.roa_unique[0]*(D_penta[1]-D_penta[0])/(PENTA_class.roa_unique[1]-PENTA_class.roa_unique[0])
            
            D_extended = np.concatenate([[D_axis],D_penta])
            
            self.D_interp[species][it] = CubicSpline(rho_extended,D_extended,extrapolate=True,bc_type='natural')
            
            c_penta = (1-self.theta) * Q / press
            c_axis = 0.0
            c_extended = np.concatenate([[c_axis],c_penta])
            
            self.c_interp[species][it] = CubicSpline(rho_extended,c_extended)
            
            # This is for bookeeping (it's not used in the calculations)
            Q_extended = np.concatenate(([0.0],Q)) # Include Q(r=0) = 0
            self.Q_interp[species] = CubicSpline(rho_extended,Q_extended,extrapolate=True,bc_type='natural')
            
            ## rename file names for bookeeping
            # it_subiter = f'{it:03}_{self.subiter:03}'
            # rename = 'mv fluxes_vs_roa fluxes_vs_roa_'+it_subiter
            # result = subprocess.run(rename, shell=True, check=True, text=True, capture_output=True)
            # rename = 'mv fluxes_vs_Er fluxes_vs_Er_'+it_subiter
            # result = subprocess.run(rename, shell=True, check=True, text=True, capture_output=True)
            # rename = 'mv flows_vs_roa flows_vs_roa_'+it_subiter
            # result = subprocess.run(rename, shell=True, check=True, text=True, capture_output=True)
            # rename = 'mv Jprl_vs_roa Jprl_vs_roa_'+it_subiter
            # result = subprocess.run(rename, shell=True, check=True, text=True, capture_output=True)
            
    def compute_diffusive_flux(self,it):
        # computes an interpolating function for Q=-n*chi*dT/dr -T*Dn*dn/dr
        
        from scipy.interpolate import CubicSpline
        import matplotlib.pyplot as plt
        
        chi = self.fluxes_info['diffusive']['chi']
        
        self.Q_interp = {}
        
        r_grid = self.rho_grid * self.aminor
        
        for species in self.list_of_species:
            
            p_r = CubicSpline(r_grid,self.P[species][it,:])
            dpdr = p_r.derivative()
            
            Q = -chi*dpdr(r_grid)
            
            D = np.zeros(self.Nr)
            
            D[1:] = np.where(dpdr(r_grid[1:])!=0, 
                            -self.theta * Q[1:] / dpdr(r_grid[1:]),
                            0.0)
            D[0] = 2*D[1] - D[2]

            self.D_interp[species][it] = CubicSpline(self.rho_grid,D,bc_type='natural',extrapolate=True)
            
            c = np.zeros(self.Nr)
            c[0] = 0
            c[1:] = (1-self.theta) * Q[1:] / p_r(r_grid[1:])
            
            self.c_interp[species][it] = CubicSpline(self.rho_grid,c,bc_type='natural',extrapolate=True)
            
            # this is for bookeeping
            self.Q_interp[species] = CubicSpline(self.rho_grid,Q)
            
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
        
    def get_collisionalHeatExchange(self):
        # returns 2 arrays: the explicit part of W_s1_s2 and the implicit fact of W_s1_s2
        
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
        
        # print(aux_B[:,0,:,0])      
        # print(W_s1_s2[:,0,:,0])                  
        
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

        # print(W_out)
                

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
        
        
        
            
         