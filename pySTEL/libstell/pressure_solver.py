"""
This library provides a python class for solving pressure transport
equations
"""

import numpy as np
import sys

# Constants
EC = 1.602176634E-19 # Electron charge [C]

# PENTA Class
class PRESSURE_SOLVER:
    
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
        self.bnd_cnds = defaultdict(dict)
        self.initial_profile = {}
        self.sources = defaultdict(lambda: defaultdict(dict))
        self.fluxes_info = defaultdict(lambda: defaultdict(dict))
                
        print(f'Solvers for pressure of {self.list_of_species} INITIALIZED!')
                
    def set_boundary_condition(self,species: str,where: str,type: str,val: float):
        # where is 'axis' or 'edge'
        # type is 'Dir' or 'Neu'
        # val is the value
        
        # check species exist in list_of_species
        if species not in self.list_of_species:
            print(f"ERROR: Species {species} is not in the plasma.")
            exit(1)
            
        # save boundary condition information in dictionary
        match where, type:
            case "axis" | "edge", "Dir" | "Neu":
                self.bnd_cnds[species][where] = {'type': type, 'val' : val}
                # self.bnd_cnds[species][where]['val'] = val
            case _:
                raise ValueError(f"Invalid where (axis/edge) or type (Dir/Neu)")
            
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
        
    def set_VMEC_equilibrium(self,wout_path=None):
        
        from libstell.vmec import VMEC
        from scipy.interpolate import CubicSpline
        
        if(wout_path is None):
            print('ASSUMING CYLINDRICAL COORDINATES: dVdr=r (EXCLUDES ALREADY 2pi*R0*2pi)')
            self.dVdr = lambda r: r 
        else:
            self.wout_path = wout_path
            # get dVdr from file; use VMEC class
            vmec_out = VMEC()
            vmec_out.read_wout(wout_path)
            
            # 4pi^2*dVds -- vmec.py already does h2f, so vprime is in full grid
            vp = vmec_out.vp[:].flatten()
            
            self.a_VMEC = vmec_out.aminor
            
            roa = np.sqrt(vmec_out.phi / vmec_out.phi[-1])
            roa = roa.flatten()
            
            #dVdr analytic = dVds * 2\rho / a
            dVdr_analytic = (2*np.pi)**2 * vp * 2.*roa / self.a_VMEC
            
            self.dVdr = CubicSpline(roa,dVdr_analytic)
            
    def set_source(self,species,source_type, rho_vals=None, source_vals=None, fraction_alpha_heating=None):
        # electrons: 'Bremsstrahlung', 'Coll_Heat_Exchange', 'Er', 'external', 'alpha_heating'
        # ionts: 'Coll_Heat_Exchange', 'Er', 'external', 'alpha_heating'
        
        from scipy.interpolate import CubicSpline
        
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
            case 'external':
                if((rho_vals is None) or (source_vals is None)):
                    print('ERROR: Need to provide rho_vals and source_vals for external source')
                    exit(1) 
                else:
                    self.sources[species][source_type] = {'interp_func' : CubicSpline(rho_vals,source_vals) }
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
            case _:
                print(f'ERROR: Source type {source_type} is NOT possible')
                exit(0)
                
    def set_fluxes(self,type: str,dkes_folder=None,surfaces=None,D_coeff=None):
        # sets type of fluxes
        # OPTION1: type='dkespenta'; dkes_folder and surfaces(list of integers) must be provided
        # OPTION2: type='diffusive'; D_coeff must be provided (diffusion coefficient of pressure in ?? units )
            
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
                
    def set_temperature(self,species,rho,density,pressure):
        # from density (m^-3) and pressure (Pa), sets temperature (eV) in plasma class
        
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
        ## check Nt=1+(tend-tstart)/dt is an integer
        if( (tend-tstart)%dt != 0 ):
            print('ERROR: (tend-tstart)%dt != 0')
            exit(0)
        else:
            Nt = int( 1+(tend-tstart)/dt )
            time = np.linspace(tstart,tend,Nt)
            
        print(' ')
        print( ' ***********************')
        print(f' *  tstart = {tstart:5.2f}s    *')
        print(f' *  tend = {tend:5.2f}s      *')
        print(f' *  dt = {dt:5.2f}s        *')
        print(f' *  drho = {drho:5.2f}       *')
        print( ' ***********************')
            
        # Initialize density and pressure using initial_profile
        press = {}
        dens = {}
        self.P = {}
        delta_p = {}
        self.all_sources = defaultdict(lambda: defaultdict(list))
        for species in self.list_of_species:
            press[species] = self.initial_profile[species](rho)
            dens[species] = self.plasma.get_density(species,rho)
            
            self.P[species] = np.zeros((Nt,Nr))
            self.P[species][0,:] = press[species]
        
        print(' ')
        header_str = '  TIME [s]     NSUB      PE_AXIS [MPa]     SE_AXIS [MW/m^3]    PI1_AXIS [MPa]    SI1_AXIS [MW/m^3]    MAX(dp/p_old)'  
        print(header_str)
        print('  '+'='*len(header_str))
        # print info from t=0
        # print(f'  {time[0]:<13.2f} {'1':<10} {self.P['electrons'][0,0]/1E6:<18.2E} {self.P['deuterium'][0,0]/1E6:<18.2E} {0.0:<13.2E}')
        
        ### LOOP IN TIME STARTING AT t=tstart+dt ###
        p_old = 1E3*np.ones(Nr) # so on loop 1 we don't divide by zero
        for it,t in enumerate(time[1:],start=1):

            ### SUBCYCLE
            delta_p_all = 10*tolerance
            subiter=1
            while(delta_p_all > tolerance or subiter<max_subiter):
            
                # sets temperature in plasma class from density and pressure for ALL species
                for species in self.list_of_species:
                    self.set_temperature(species,rho,dens[species],press[species])
                    
                # call PENTA3 and compute interpolating functions self.Er_interp, self.Gamma_interp[species] and self.Q_interp[species]
                self.call_PENTA3()
                
                #solver for each species (this can be parallelized... numba??)   
                for species in self.list_of_species:
                    
                    # compute sources on grid (1D-array)
                    total_sources = self.get_sources(species,rho)  # W/m^3
                    
                    # compute divergence-of-heat-flux term
                    div_flux = self.get_div_flux(species,rho) # J/ (m^3.s)
                    
                    # time evolution (don't forget 3/2 term)
                    press[species] = self.P[species][it-1,:] + (2./3)*dt*(-div_flux+total_sources)
                    
                    ### NEED TO THINK : BOUNDARY CONDITIONS before OR after TIME EVOLUTION?
                    ### OR EVEN: before AND after ???!!!
                    
                    # impose Dirichlet Boundary Conditions
                    if (self.bnd_cnds[species]['axis']['type'] == 'Dir'):
                        press[species][0] = self.bnd_cnds[species]['axis']['val']
                    if (self.bnd_cnds[species]['edge']['type'] == 'Dir'):
                        press[species][-1] = self.bnd_cnds[species]['edge']['val']
                        
                    # impose Neumann Boundary conditions (on axis it comes naturally that dp/dr=0)
                    if (self.bnd_cnds[species]['edge']['type'] == 'Neu'):
                        press[species][-1] = (2./3)*(2*self.P[species][it,-2] - 0.5*self.P[species][it,-3] + self.bnd_cnds[species]['edge']['val'])
                    if (self.bnd_cnds[species]['axis']['type'] == 'Neu' and self.bnd_cnds[species]['axis']['val'] != 0):
                        print('Neumann value != 0 on axis is NOT possible !!')
                        exit(0)

                    delta_p[species] = np.max( np.where( p_old>1E-10, np.abs((press[species]-p_old)/p_old), 0 ) )
                    
                    p_old = press[species]
                    self.P[species][it,:] = press[species]
                    
                delta_p_all = np.max([np.max(value) for value in delta_p.values()])
                
                info_str = f'  {t:<13.2f}{subiter:<10}{self.P['electrons'][it,0]/1E6:<18.2E}{self.all_sources['electrons'][-1,0]/1E6:<20.2E}{self.P['deuterium'][it,0]/1E6:<18.2E}{self.all_sources['deuterium'][-1,0]/1E6:<21.2E}{delta_p_all:<13.2E}'
                print(info_str)
                
                subiter += 1
                
                        

    def get_sources(self,species: str,rho_grid):
        # returns 1D-array of same size as rho_grid
        # computes sources using info in self.sources[species]
        
        from libstell.fusion import FUSION
        from collisions import COLLISIONS
        
        fusion = FUSION()
        
        # checks rho_grid is a 1D array
        rho_grid = np.asarray(rho_grid)
        if( rho_grid.ndim != 1):
            raise ValueError("rho_grid must be 1D")
        
        total_source = np.zeros(len(rho_grid))
        
        for source_type in self.sources[species]:
            
            aux_source = np.zeros(len(rho_grid))         
            match source_type:
                case 'Bremsstrahlung':
                    for ir,rho in enumerate(rho_grid):
                        for ion in self.plasma.ion_species:
                            zi = self.plasma.Zcharge[ion]
                            ni = self.plasma.get_density(ion,rho)
                            ne = self.plasma.get_density('electrons',rho)
                            Te = self.plasma.get_temperature('electrons',rho)
                            
                            aux_source[ir] -= fusion.BremsstrahlungPower(zi,ni,ne,Te)
                
                case 'external':
                    aux_source = self.sources[species]['external']['interp_func'](rho_grid)
                    
                case 'Coll_Heat_Exchange':
                    collisions = COLLISIONS()
                    for ir,rho in enumerate(rho_grid):
                        n1 = self.plasma.get_density(species,rho)
                        T1 = self.plasma.get_temperature(species,rho)
                        m1 = self.plasma.mass[species]
                        Z1 = self.plasma.Zcharge[species]
                        for species2 in self.list_of_species:
                            n2 = self.plasma.get_density(species2,rho)
                            T2 = self.plasma.get_temperature(species2,rho)
                            m2 = self.plasma.mass[species2]
                            Z2 = self.plasma.Zcharge[species2]
                            
                            aux_source[ir] += collisions.collisionalHeatExchange(n1,T1,m1,Z1,n2,T2,m2,Z2)
                            
                case 'alpha_heating':
                    nD = self.plasma.get_density('deuterium', rho_grid)
                    nT = self.plasma.get_density('tritium', rho_grid)
                    
                    Ti = 0.5* ( self.plasma.get_temperature('deuterium', rho_grid) + self.plasma.get_temperature('tritium', rho_grid) )
                    
                    sigmav = [fusion.sigmaBH(ti,'DT') for ti in Ti]
                    
                    S_alpha = nD * nT * sigmav *  fusion.E_DT_He # W/m^3
                    
                    fraction_alpha_heating = self.sources[species]['alpha_heating']['fraction_alpha_heating']
                    
                    aux_source = S_alpha * fraction_alpha_heating
                    
                case 'Er':
                    aux_source = self.plasma.charge[species]*self.Er_interp(rho_grid)*self.Gamma_interp[species](rho_grid)

                case _:
                    print(f'ERROR: Source type {source_type} not defined....')
                    exit(0)
                    
            total_source += aux_source
            
            #save in dictionary for bookeeping
            self.all_sources[species][source_type].append(aux_source)
                    
                    
        return total_source
                    
    def get_div_flux(self,species: str,rho_grid):
        # returns 1D-array of same size as rho_grid
        # computes div of heat flux term using self.QoT_interp[species] and self.dVdr
        # uses a 2nd order central finite difference scheme (where fluxes are conserved)
        
        drho = rho_grid[1]-rho_grid[0]
        dr = self.a_VMEC * drho
        
        Q = self.Q_interp[species]
        Vp = self.dVdr
        a = self.a_VMEC
        
        div_flux = np.zeros(len(rho_grid))
        
        # r=0
        div_flux[0] = 2*Q(drho) / dr
        
        for ir,rho in enumerate(rho_grid[1:-1],start=1):
            
            rp = rho + drho/2
            rm = rho - drho/2
            
            div_flux[ir] = ( Vp(rp)*Q(rp) - Vp(rm)*Q(rm) ) / (Vp(rho)*dr)
        
        # r=a (backward finite difference)
        div_flux[-1] = 1.5*Vp(a)*Q(a) - 2*Vp(a-drho)*Q(a-drho) + 0.5*Vp(a-2*drho)*Q(a-2*drho)
        div_flux[-1] = div_flux[-1] / (Vp(a)*dr)
            
        return div_flux
                                              
    def make_checks(self):
        
        # check VMEC equilibrium exists
        if(not hasattr(self,'dVdr')):
            print('dVdr MUST BE SET!!')
            exit(1)
        
        # check boundary conditions exist for all species
        for species in self.list_of_species:
            if species not in self.bnd_cnds:
                raise KeyError(f"Missing boundary conditions for species: {species}")
            
        # check initial profiles are set for all species
        for species in self.list_of_species:
            if species not in self.initial_profile:
                raise KeyError(f"Missing initial profile for species: {species}")
            
        # check fluxes info is set
        if(not hasattr(self,'fluxes_info')):
            print('set_fluxes must be called before running!!')
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
            # print(result.stderr)
            
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
            
# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)
        
        
        
            
        