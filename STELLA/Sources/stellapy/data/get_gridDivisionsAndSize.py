''' MODULE GET GRID DIVISIONS 

Obtain the divisions in kx and ky for non-linear runs.

For non-linear runs we simulate a box of size (x0, 2*pi*y0) around the magnetic field line.
The input parameters are (y0, nzed, ny and nx); then (x0, nakx, naky) are calculated by stella.

The resolution of the simulation is:
    - nzed along the magnetic field line
    - ny along the y-axis with length 2*pi*y0  
    - nx along the x-axis with length x0 ~ 1/hat(s)  with  hat(s) = -2s/iota  with  s = (r/a)^2
    
Theory
------
The size of the grid and it's divisions are based on the standard parallel boundary condition.
Their values are calculated according to the paper: "The parallel boundary condition 
for turbulence simulations in low magnetic shear devices." by M. F. Martin [1]
Or according to "Definitions for GS2 full-flux-surface stellarator geometry" by
M. Landreman which can be found in "stella/geo/vmec_interface/doc" [2]
'''

import numpy as np
import math as math

#################################################################
#                       GRID DIVISIONS
#################################################################

def get_naky(ny):
    ''' Returns the number of modes <ky> along the y-direction, when <ny> is given. 
    Note that only positive ky values are simulated due to the mirror condition. 
    The real amount of ky modes is (2*naky-1).'''
    naky = (ny-1)/3 + 1
    return int(naky)


def get_nakx(nx):
    ''' Returns the number of modes <kx> along the x-direction, when <nx> is given. '''
    nakx = 2*((nx-1)/3) + 1
    return int(nakx)


def get_dky(y0):
    ''' Returns the step size between the <ky> modes. '''
    return 1.0/y0


def get_dkx(woutfile, svalue, nfield_periods, y0, jtwist=None, shat=None, iota=None):
    ''' Returns the step size between the <kx> modes. '''
    
    # Non-quantized boundary condition assumed to be periodic instead of linked boundary conditions if zero magnetic shear
    if abs(get_shat(woutfile, svalue)) <= shat_zero():
        jtwist = get_jtwist(woutfile, svalue, nfield_periods)
        dky    = get_dky(y0)
        dkx    = dky/jtwist
        
    # The division along x is based on the twist-and-shift boundary condition
    else:
        dky        = get_dky(y0)
        jtwist     = get_jtwist(woutfile, svalue, nfield_periods) if jtwist is None else jtwist
        geo_factor = twist_and_shift_geo_fac(woutfile, svalue, nfield_periods, shat=shat, iota=iota)
        dkx        = dky*geo_factor/jtwist
    return abs(dkx)


def get_aky(ny, y0):
    ''' Returns a list of the <ky> values simulated by stella, when <ny> and <y0> are given. '''

    # Get the step size between the <ky> values and the number of <ky> values
    dky  = get_dky(y0)
    naky = get_naky(ny)

    # Create a list with <naky> values of <ky> 
    aky = [ round(iky*dky,2) for iky in np.arange(0, naky) ] 
    return aky


def get_akx(woutfile, svalue, nfield_periods, y0, nx, jtwist=-1):
    ''' Returns a list of the <kx> values simulated by stella, when <nx>, <y0> and <s> are given. '''
    
    # Get the step size between the <kx> values and the number of <kx> values
    dkx  = get_dkx(woutfile, svalue, nfield_periods, y0, jtwist)
    nakx = get_nakx(nx)
    
    # get the ikx index corresponding to kx_max
    ikx_max = int(get_nakx(nx)/2+1)
    
    # kx goes from zero to kx_max down to zero and then from -kx_max to -|kx_min|    
    akx = [ round(ikx*dkx,2) for ikx in np.arange(0, ikx_max) ]       
    akx = akx + [ round((ikx-nakx)*dkx,2) for ikx in np.arange(ikx_max, nakx) ]
    return akx

#################################################################
#                      GRID SIZE
#################################################################

def get_Lx(woutfile, svalue, nfield_periods, y0, jtwist=-1, nperiod=1, nfp=5):
    ''' Returns the length of the box along x.
    
    Definitions
    -----------
    Lx = 2*pi/dkx = 2*pi*jtwist/(factor*dky) = -jtwist/(shat*iota*field_period_ratio*dky)
    Lx = L/(P*shat*dky) since Lx/Ly = jtwist/(2*np.pi*P*shat)
    '''
    
    # Get the factors
    dkx  = get_dkx(woutfile, svalue, nfield_periods, y0, jtwist=None, shat=None, iota=None) 
    dky  = get_dky(y0)
    iota = get_iota(woutfile, svalue)
    shat = get_shat(woutfile, svalue=svalue)
    jtwist = get_jtwist(woutfile, svalue, nfield_periods)
    field_period_ratio = nfield_periods/nfp

    # Calculation through the divisiona long kx
    Lx = 2*np.pi/dkx
    
    # Direct calculation through Lx/Ly = jtwist/(2*np.pi*P*shat) so Lx = L/(P*shat*dky)
    P = -iota*field_period_ratio
    Lx = jtwist/(shat*P*dky)
    return Lx


def get_Lx_over_Ly(woutfile, svalue, nfield_periods, y0=10, nfp=5):
    ''' Returns the perpendicular aspect ratio of the domain (see eq.11 in [1]).  '''

    # Official equation
    field_period_ratio = nfield_periods/nfp
    jtwist = get_jtwist(woutfile, svalue, nfield_periods)
    shat   = get_shat(woutfile, svalue=svalue)
    iota = get_iota(woutfile, svalue)
    factor = twist_and_shift_geo_fac(woutfile, svalue, nfield_periods)
                                     
    # The end of the zeta domain is defined as zeta = 2*pi*P                           
    P = factor/(2*np.pi*shat)
    P = -iota*field_period_ratio
    
    # Official equation  
    LxLy = jtwist/(2*np.pi*P*shat)
    
    # Direct calculation
    LxLy = get_Lx(woutfile, svalue, nfield_periods, y0=100)/(2*np.pi*100)
    return LxLy


#################################################################
#                   MATHEMATIC FACTORS
#################################################################

def shat_zero():
    ''' Shat_zero is minimum shat value below which the periodic boundary condition is enforced. '''
    return 1.E-2

def get_shat(woutfile, svalue):
    ''' Returns hat(s) = x/q*(dq/dx) =  -2*s/iota*diota/ds. 
    
    Definitions
    -----------
    q = 1/iota 
    x = a*sqrt(s)
    s = psi_toroidal/psi_edge
    
    Derivation
    ----------
    shat = x/q*(dq/dx) 
         = (a*sqrt(s))/(1/iota) * (d(1/iota)/d(a*sqrt(s))
         = (a*sqrt(s)*iota) * (-1/iota^2) * (2*sqrt(s)/a) * diota/ds
         = -2 * s/iota * diota/ds
    
    '''
    return (-2 * svalue / get_iota(woutfile, svalue)) * get_diotaOverds(woutfile, svalue)

def get_iota(woutfile, svalue, iotas=None, ns=None):
    ''' Calculate the rotational transform iota at s=svalue if <iotas> and <ns> are given, 
    or read it through read_wout(woutfile, s=svalue).
    
    Definition
    ----------
    iota = 1/q = psi_poloidal/psi_toroidal = R*B_theta/(r*B_phi)
    '''
    
    # Read the iotas and ns from the "wout" file which calculates iota 
    if woutfile is not None:
        from .read_wout import read_wout
        wout_data = read_wout(woutfile, svalue=svalue)
        return wout_data['iota'] 
    
    # Calculate iota
    if iotas is not None and ns is not None:
        vmec_radial_i_half, vmec_radial_w_half = get_vmecRadialWeightHalf(ns, svalue)
        return iotas[vmec_radial_i_half[0]] * vmec_radial_w_half[0] + iotas[vmec_radial_i_half[1]] * vmec_radial_w_half[1]        
           
def get_diotaOverds(woutfile, svalue):
    ''' Return diota/ds at s=svalue. '''
    
    # Read the wout data
    from .read_wout import read_wout
    wout_data = read_wout(woutfile, svalue=svalue)
    ns = wout_data['ns']
    iotaf = wout_data['iotaf']
    
    # Get the weighing factors
    vmec_radial_i_half, vmec_radial_w_half = get_vmecRadialWeightHalf(ns, svalue)
    
    # Step size of the divisions along s=(r/a)^2=rho^2 during the VMEC run
    ds = 1.0/float(ns-1)
    
    # Calculate diota/ds
    diotaOverds_on_half_grid         = np.empty(ns)*0.0 
    diotaOverds_on_half_grid[1:ns-1] = (iotaf[1:ns-1] - iotaf[0:ns-2]) / ds
    diotaOverds = diotaOverds_on_half_grid[vmec_radial_i_half[0]] * vmec_radial_w_half[0]+\
                  diotaOverds_on_half_grid[vmec_radial_i_half[0]] * vmec_radial_w_half[1]
    
    return diotaOverds


def get_jtwist(woutfile, svalue, nfield_periods, shat=None, iota=None, nfp=5):
    ''' J is a non-zero integer that can be set in the code to 
    potentially achieve a more desirable aspect ratio. '''
    factor = twist_and_shift_geo_fac(woutfile, svalue, nfield_periods, shat=shat, iota=iota, nfp=nfp)
    jtwist = max(1,int(abs(factor)+0.5))
    return jtwist

def twist_and_shift_geo_fac(woutfile, svalue, nfield_periods, nfp=5, shat=None, iota=None):
    ''' A factor used directly in the stella code: 
    
    Signs
    -----
    -sign(psi_tor) = sign(psi) = dxdpsi_sign*sign_torflux
    dxdpsi_sign = -1
    dydalpha_sign = 1
    
    Definitions
    -----------
    psi = -psi_tor
    rho = rhotor = sqrt(s)  = r/a = sqrt(psitor/psitor_LCFS) = sqrt(torflux)
    Bref = 2*abs(psi_tor_LCFS)/a^2
    a*Bref*dx/dpsi_tor = sign(psi_tor)/rhotor
    dxdpsi_N = a*Bref*dx/dpsi = -a*Bref*dx/dpsi_tor = -sign(psi_tor)/rhotor = dxdpsi_sign*sign_torflux/rhotor; 
    dydalpha_N = (dy/dalpha) / a = sign(dydalpha)*rhotor
    psi_N = -psitor/(a**2*Bref)
    drhodpsi_N = -drho/d(rho**2)*(aref**2*Bref/psitor_lcfs) = -1.0/rho = dxdpsi_sign*sign_torflux/rhotor
    field_period_ratio = nfield_periods/nfp  with nfp = 5 for Wendelstein 7-X
    
    Derivation
    -----------
    abs(twist_and_shift_geo_fac) = dkx/dky*jtwist
    twist_and_shift_geo_fac = -2.*pi*shat*drhodpsi*dydalpha/(qinp*dxdpsi*rhotor)*field_period_ratio
    
    '''

    # Define the signs of the derivatives
    dxdpsi_sign   = -1.
    dydalpha_sign =  1.
    
    # Get the factors
    dydalpha_N = dydalpha_sign*np.sqrt(svalue)                      
    drhodpsi_N = dxdpsi_sign*sign_torflux(woutfile)/np.sqrt(svalue)  # Opposite sign than in marconi 
    dxdpsi_N   = dxdpsi_sign*sign_torflux(woutfile)/np.sqrt(svalue)  # Opposite sign than in marconi 
    field_period_ratio = nfield_periods/nfp
    
    # If shat and iota are not read from vmecgeo, try to calculate them:
    if not shat:
        shat = get_shat(woutfile, svalue)
    if not iota:
        iota = get_iota(woutfile, svalue)

    # The factor calculated by marconi is: 
    factor = -2.*np.pi*shat*iota*drhodpsi_N*dydalpha_N/(dxdpsi_N*np.sqrt(svalue))*field_period_ratio
    
    # However this can be simplied too:
    factor = -2.*np.pi*shat*iota*field_period_ratio
    return factor
#################################################################
#                         SIGNS
#################################################################

def sign_torflux(woutfile):
    ''' Returns the sign of the toroidal flux. '''
    return np.sign(edge_toroidal_flux_over_2pi(woutfile))   

def edge_toroidal_flux_over_2pi(woutfile):
    ''' Returns psi_edge/(2*pi). '''
    from .read_wout import read_wout # Circular dependency 
    wout_data = read_wout(woutfile)
    toroidalFlux = wout_data['phi'] # Toroidal flux on full mesh
    dimension = np.size(toroidalFlux)
    return toroidalFlux[dimension-1] 

#################################################################
#                   OTHERS
#################################################################

def get_vmecRadialWeightHalf(ns, svalue):
    ''' Interpolation for an <s> value between two flux surfaces, 
    since this value is not in the output file. '''

    # Initiate the arrays
    vmec_radial_index_half  = np.empty(2, dtype='int')
    vmec_radial_weight_half = np.empty(2, dtype='float')
    normalizedToroidalFlux_half_grid = np.empty(ns-1, dtype='float')
    normalizedToroidalFlux_full_grid = np.linspace(0, 1, ns)

    # *_index_* referes to the position of the two neighbouring flux surfaces of the half mesh for a given svalue
    for i in np.arange(0, ns-1):
        normalizedToroidalFlux_half_grid[i]=(normalizedToroidalFlux_full_grid[i]+normalizedToroidalFlux_full_grid[i+1])*0.5

    if svalue < normalizedToroidalFlux_half_grid[0]:
        print("Warning: extrapolating beyond the end of VMEC's half grid.")
        print("Extrapolating towards the magnetic axis. Results are likely to be inaccurate.")
        # We start at element 2 since element 1 is always 0 for quantities on the half grid.
        vmec_radial_index_half[0]  = 1
        vmec_radial_index_half[1]  = 2
        vmec_radial_weight_half[0] = (normalizedToroidalFlux_half_grid[1] - svalue)/\
                                     (normalizedToroidalFlux_half_grid[1] - normalizedToroidalFlux_half_grid[0])
                                     
    elif svalue > normalizedToroidalFlux_half_grid[ns-2]:
        # We are beyond the last point of the half grid
        print("Warning: extrapolating beyond the end of VMEC's half grid.")
        print("(Extrapolating towards the last closed flux surface.) Results may be inaccurate.")
        vmec_radial_index_half[0]  = ns-2
        vmec_radial_index_half[1]  = ns-1
        vmec_radial_weight_half[0] = (normalizedToroidalFlux_half_grid[ns-2] - svalue)/\
                                     (normalizedToroidalFlux_half_grid[ns-2] -\
                                      normalizedToroidalFlux_half_grid[ns-3])
    elif svalue == normalizedToroidalFlux_half_grid[ns-2]:
        # We are exactly at the last point of the half grid
        vmec_radial_index_half[0]  = ns-2
        vmec_radial_index_half[1]  = ns-1
        vmec_radial_weight_half[0] = 0.0
    else:
        # normalizedToroidalFlux_used is inside the half grid.
        # This is the most common case.
        vmec_radial_index_half[0] = math.floor(svalue*(ns-1) + 0.5) # Different from Matt!!
        if vmec_radial_index_half[0] < 1:
            # This can occur sometimes due to roundoff error.
            vmec_radial_index_half[0] = 1
        vmec_radial_index_half[1]  = vmec_radial_index_half[0] + 1
        vmec_radial_weight_half[0] = vmec_radial_index_half[0] - svalue*(ns-1) + 1.5
    vmec_radial_weight_half[1] = 1.0 - vmec_radial_weight_half[0]
    return vmec_radial_index_half, vmec_radial_weight_half
