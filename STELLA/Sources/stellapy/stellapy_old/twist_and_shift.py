from read_wout import *
import numpy as np

def kGrid(equil=None, psitor=0.5, y0=14.14, jtwist=-1, nx=100, ny=64):
    
    xv, yv = np.meshgrid(get_akx(equil, psitor, y0, nx, jtwist),\
                         get_aky(ny, y0))
    # plot ccw path
    fig, ax = plt.subplots()
    ax.set_xlim((xv.min(),xv.max()))
    ax.set_ylim((yv.min(),yv.max()))
    ax.set_xlabel('$k_y\\rho_{r}$')
    ax.set_ylabel('$k_x\\rho_{r}$')
    ax.scatter(xv.flatten(),yv.flatten(), color='red', marker='.', cmap="RdYlGn")
    plt.show()


def twist_and_shift_geo_fac(equil=None, psitor=0.5):

    dxdpsi_sign        = -1. # Ni idea por qué. Debería ser igual a sign_toroidal_flux
    dydalpha_sign      =  1.
    dydalpha = dydalpha_sign*rhotor(psitor)
    drhodpsi = dxdpsi_sign*sign_torflux(equil)/rhotor(psitor)
    dxdpsi   = dxdpsi_sign*sign_torflux(equil)/rhotor(psitor)

    return -2.*pi*get_shat(equil, psitor)*\
           (1/get_iota(equil, psitor))*\
           drhodpsi*dydalpha/(dxdpsi*rhotor(psitor))

def get_naky(ny=100):
    return (ny-1)/3 + 1

def get_nakx(nx=64):
    return 2*((nx-1)/3) + 1

def get_aky(ny=100, y0=10.0):
    aky = empty(get_naky(ny),dtype='float')
    dky = get_dky(y0)
    for iky in arange(0, get_naky(ny)):
        aky[iky] = float(iky)*dky
    return aky

def get_akx(equil=None, psitor=0.49, y0=10.0, nx=64, jtwist=-1):
    # get the ikx index corresponding to kx_max
    ikx_max = int(get_nakx(nx)/2+1)

    # get the total number of ky values, including negative ky
    #naky_all = 2*get_naky(ny)-1

    # kx goes from zero to kx_max down to zero...
    akx = empty(get_nakx(nx), dtype='float')
    dkx = get_dkx(equil, psitor, y0, nx)
    
    for ikx in arange(0, ikx_max):
        akx[ikx] = float(ikx)*get_dkx(equil, psitor, y0, nx, jtwist)
        
    # and then from -kx_max to -|kx_min|
    for ikx in arange(ikx_max, get_nakx(nx)):
        akx[ikx] = float(ikx-get_nakx(nx))*get_dkx(equil, psitor, y0, nx, jtwist)

    return akx

def get_dky(y0=None):
    return 1.0/y0

def get_dkx(equil=None, psitor=0.49, y0=10.0, nx=64, jtwist=-1):
    if abs(get_shat(equil, psitor)) <= shat_zero():
        dkx = get_dky(y0) / float(get_jtwist())
    else:
        dkx = get_dky(y0) * abs(twist_and_shift_geo_fac(equil, psitor))/\
              float(get_jtwist(equil, psitor, jtwist))
    return dkx

def shat_zero():
    return 1E-2

def get_jtwist(equil=None, psitor=None, jtwist=-1):
    if (jtwist < 1):
        jtwist = max(1,int(abs(twist_and_shift_geo_fac(equil, psitor))+0.5))
    return jtwist



