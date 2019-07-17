# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 12:28:08 2019

@author: gawe
"""

# ===================================================================== #
# ===================================================================== #

from __future__ import absolute_import, with_statement, absolute_import, \
                       division, print_function, unicode_literals
import numpy as _np

__metaclass__ = type

mu0=4.0e-7*_np.pi

# ===================================================================== #
# ===================================================================== #


def calc_jll(vmec_data, theta, zeta, use_nyq=True):
    """
         CALC_JLL(vmec_data,theta,zeta) Calculates the parallel current density.
         This funciton takes a VMEC data structure (as read by read_vmec) and
         theta/zeta arrays as input and outputs the parallel current density.

         Example usage (Matlab)
              theta=0:2*pi/359:2*pi;
              zeta=0:2*pi/63:2*pi;
              data=read_vmec('wout.test');        % Reads VMEC wout file
              jll=calc_jll(vmec_data,theta,zeta); % Calculate the current

         Example usage (Python)
              theta=_np.linspace(0, 2*_np.pi, 360)
              zeta=_np.linspace(0, 2*_np.pi, 64)
              vmec_data=read_vmec('wout.nc')
              jll=calc_jll(vmec_data, theta, zeta)

         Maintained by: Samuel Lazerson (lazerson@pppl.gov)
         Ported to python by Jonathan Schilling (jonathan.schilling@ipp.mpg.de)
         Version:       1.01

        July 16th, 2019. GMW -  Merged differences between libstell /
                    matlabVMEC / VMECplot repo's while maintaining compatability
    """
    mkey = 'xm'    # matlabVMEC version
    nkey = 'xn'
    if use_nyq:
        mkey += '_nyq'    # libstell version
        nkey += '_nyq'
    # end if

    b  = cfunct(theta,zeta,vmec_data['bmnc'],    vmec_data[mkey],vmec_data[nkey])
    g  = cfunct(theta,zeta,vmec_data['gmnc'],    vmec_data[mkey],vmec_data[nkey])
    bu = cfunct(theta,zeta,vmec_data['bsubumnc'],vmec_data[mkey],vmec_data[nkey])
    bv = cfunct(theta,zeta,vmec_data['bsubvmnc'],vmec_data[mkey],vmec_data[nkey])
    ju = cfunct(theta,zeta,vmec_data['currumnc'],vmec_data[mkey],vmec_data[nkey])
    jv = cfunct(theta,zeta,vmec_data['currvmnc'],vmec_data[mkey],vmec_data[nkey])

    if (vmec_data['iasym']):
        b  += sfunct(theta,zeta,vmec_data['bmns'],    vmec_data[mkey],vmec_data[nkey])
        g  += sfunct(theta,zeta,vmec_data['gmns'],    vmec_data[mkey],vmec_data[nkey])
        bu += sfunct(theta,zeta,vmec_data['bsubumns'],vmec_data[mkey],vmec_data[nkey])
        bv += sfunct(theta,zeta,vmec_data['bsubvmns'],vmec_data[mkey],vmec_data[nkey])
        ju += sfunct(theta,zeta,vmec_data['currumns'],vmec_data[mkey],vmec_data[nkey])
        jv += sfunct(theta,zeta,vmec_data['currvmns'],vmec_data[mkey],vmec_data[nkey])
    # end if

    jll = (bu*ju+bv*jv)/(g*b)
    return jll

# ===================================================================== #


# Calculate Currents
def calc_curr(f):
    f['currumnc']=_np.zeros([f['ns'], f['mnmax_nyq']])
    f['currvmnc']=_np.zeros([f['ns'], f['mnmax_nyq']])

    ohs = f['ns']-1.0
    hs  = 1.0/_np.double(ohs)
    ns = f['ns']

    shalf=_np.zeros(ns)
    sfull=_np.zeros(ns)
    for i in _np.arange(1,ns):
        shalf[i] = _np.sqrt(hs*(i-0.5))   # check this...
        sfull[i] = _np.sqrt(hs*(i-0.0))   # check this

    js1 = _np.arange(2,ns)
    js  = _np.arange(1,ns-1)

    # calculation of curl B
    for mn in _np.arange(f['mnmax_nyq']):
        if (_np.mod(f['xm_nyq'][mn],2) == 1):
            t1  = 0.5*(shalf[js1] * f['bsubsmns'][js1,mn] +
                       shalf[js]  * f['bsubsmns'][js,mn]   ) / sfull[js]
            bu0 = f['bsubumnc'][js,mn]/shalf[js]
            bu1 = f['bsubumnc'][js1,mn]/shalf[js1]
            t2  = ohs*(bu1-bu0)*sfull[js]+0.25*(bu0+bu1)/sfull[js]
            bv0 = f['bsubvmnc'][js,mn]/shalf[js]
            bv1 = f['bsubvmnc'][js1,mn]/shalf[js1]
            t3  = ohs*(bv1-bv0)*sfull[js]+0.25*(bv0+bv1)/sfull[js]
        else:
            t1  = 0.5*(f['bsubsmns'][js1,mn]+f['bsubsmns'][js,mn])
            t2  = ohs*(f['bsubumnc'][js1,mn]-f['bsubumnc'][js,mn])
            t3  = ohs*(f['bsubvmnc'][js1,mn]-f['bsubvmnc'][js,mn])
        f['currumnc'][js,mn] = -_np.double(f['xn_nyq'][mn])*t1 - t3
        f['currvmnc'][js,mn] = -_np.double(f['xm_nyq'][mn])*t1 + t2

    # linear extrapolation to axis
    for i in _np.arange(f['mnmax_nyq']):
        if (f['xm_nyq'][i]<=1):
            f['currumnc'][0,i]=2.0*f['currumnc'][1,i]-f['currumnc'][2,i]
            f['currvmnc'][0,i]=2.0*f['currvmnc'][1,i]-f['currvmnc'][2,i]
        else:
            f['currumnc'][0,i]=0.0
            f['currvmnc'][0,i]=0.0

    # linear extrapolation to LCFS
    f['currumnc'][f['ns']-1,:]=2.0*f['currumnc'][f['ns']-2,:]-f['currumnc'][f['ns']-3,:]
    f['currvmnc'][f['ns']-1,:]=2.0*f['currvmnc'][f['ns']-2,:]-f['currvmnc'][f['ns']-3,:]

    # scaling to SI units
    f['currumnc']=f['currumnc']/mu0;
    f['currvmnc']=f['currvmnc']/mu0;

    if f['iasym']:
        f['currumns']=_np.zeros([f['ns'], f['mnmax_nyq']])
        f['currvmns']=_np.zeros([f['ns'], f['mnmax_nyq']])

        for mn in _np.arange(f['mnmax_nyq']):
            if (_np.mod(f['xm_nyq'][mn],2) == 1):
                t1  = 0.5*(shalf[js1] * f['bsubsmnc'][js1,mn] +
                           shalf[js]  * f['bsubsmnc'][js,mn]   ) / sfull[js]
                bu0 = f['bsubumns'][js,mn]/shalf[js1]
                bu1 = f['bsubumns'][js1,mn]/shalf[js1]
                t2  = ohs*(bu1-bu0)*sfull[js]+0.25*(bu0+bu1)/sfull[js]
                bv0 = f['bsubvmns'][js,mn]/shalf[js]
                bv1 = f['bsubvmns'][js1,mn]/shalf[js1]
                t3  = ohs*(bv1-bv0)*sfull[js]+0.25*(bv0+bv1)/sfull[js]
            else:
                t1  = 0.5*(f['bsubsmnc'][js1,mn]+f['bsubsmnc'][js,mn])
                t2  = ohs*(f['bsubumns'][js1,mn]-f['bsubumns'][js,mn])
                t3  = ohs*(f['bsubvmns'][js1,mn]-f['bsubvmns'][js,mn])
            f['currumns'][js,mn] = _np.double(f['xn_nyq'][mn])*t1 - t3
            f['currvmns'][js,mn] = _np.double(f['xm_nyq'][mn])*t1 + t2

        for i in _np.arange(f['mnmax_nyq']):
            if (f['xm_nyq'][i]<=1):
                f['currumns'][0,i]=2.0*f['currumns'][1,i]-f['currumns'][2,i]
                f['currvmns'][0,i]=2.0*f['currvmns'][1,i]-f['currvmns'][2,i]
            else:
                f['currumns'][0,i]=0.0
                f['currvmns'][0,i]=0.0

        f['currumns'][f['ns']-1,:]=2.0*f['currumns'][f['ns']-2,:]-f['currumns'][f['ns']-3,:]
        f['currvmns'][f['ns']-1,:]=2.0*f['currvmns'][f['ns']-2,:]-f['currvmns'][f['ns']-3,:]

        f['currumns']=f['currumns']/mu0;
        f['currvmns']=f['currvmns']/mu0;
    return f

# ======================================================================== #
# ======================================================================== #


# call this with the result of read_vmec:
# vmec_data=read_vmec("wout_w7x_ref_60.nc")
# finish_import(vmec_data)
# to transform the half-mesh quantities onto the full grid
# This has to be done with some care as p.ex. booz_xform rely
# on having the half-mesh quantities after the import.
# Most likely this is only useful for plotting...
def finish_import(vmec_data):
    ns=vmec_data['ns']

    # new definition => adapt to NESCOIL ???
    vmec_data['xn'] = -vmec_data['xn']

    # Put on full grid
    vmec_data['buco'] = h2f(vmec_data['buco'], ns)
    vmec_data['bvco'] = h2f(vmec_data['bvco'], ns)
    vmec_data['vp'] = h2f(vmec_data['vp'], ns)
    vmec_data['overr'] = h2f(vmec_data['overr'], ns)
    vmec_data['specw'] = h2f(vmec_data['specw'], ns)

    # Put matrix quantities on full grid
    for key in ['bmnc','gmnc','lmns','bsupumnc','bsupvmnc','bsubsmns','bsubumnc','bsubvmnc']:
        vmec_data[key][0,:] = 1.5 * vmec_data[key][1,:] - 0.5 * vmec_data[key][2,:]
        vmec_data[key][1:ns-2,:] = 0.5 * (vmec_data[key][1:ns-2,:] + vmec_data[key][2:ns-1,:])
        vmec_data[key][ns-1,:] = 2.0 * vmec_data[key][ns-2,:] - vmec_data[key][ns-3,:]
    # end for key in dict

    if vmec_data['iasym']:
        for key in ['bmns','gmns','lmnc','bsupumns','bsupvmns','bsubsmnc','bsubumns','bsubvmns']:
            vmec_data[key][0,:] = 1.5 * vmec_data[key][1,:] - 0.5 * vmec_data[key][2,:]
            vmec_data[key][1:ns-2,:] = 0.5 * (vmec_data[key][1:ns-2,:] + vmec_data[key][2:ns-1,:])
            vmec_data[key][ns-1,:] = 2.0 * vmec_data[key][ns-2,:] - vmec_data[key][ns-3,:]
        # end for key in dict
    # end if is asymmetric
# end def finish_import

# ===================================================================== #

# copied from libstell.py by SAL to circumvent libstell import here
def h2f(var, ns):
    # explicit typing here is a good touch, but it will throw errors with lists
    # ... our data is expected as numpy arrays so this isn't bad
    #     at the output this should be a numpy array
    var = _np.copy(var)
    if len(_np.shape(var)) == 1:
        var = _np.atleast_2d(var).T
    # end if
    temp = _np.zeros((ns,1), dtype=var.dtype)
    temp[0] = 1.5 * var[0] - 0.5 * var[1]
    temp[1:ns-1] = 0.5* (var[0:ns-2] + var[1:ns-1])
    temp[ns-1] = 1.5 * var[ns-2] - 0.5 * var[ns-3]
    return temp

def h2f_special(var, ns):
    # explicit typing here is a good touch, but it will throw errors with lists
    # ... our data is expected as numpy arrays so this isn't bad
    #     at the output this should be a numpy array
    var = _np.copy(var)
    # Map quantitiy from half to full grid
    temp = _np.zeros( (ns,1), dtype=var.dtype)
    # temp[0] = 1.5*var[0] - 0.5*var[1]
    temp[1:ns-1] = 0.5*(var[1:ns-1] + var[2:ns])
    temp[ns-1] = 1.5*var[ns-1] - 0.5*var[ns-2]
#        temp[1:-1] = 0.5*(var[1:ns-2] + var[2:ns-1])
#        temp[-1] = 1.5*var[ns-1] - 0.5*var[ns-2]
    return temp
#end def h2f()

def cfunct(theta,zeta,fmnc,xm,xn):
    f=0
    (ns,mn)=fmnc.shape
    lt = len(theta)
    lz = len(zeta)
    mt=_np.matmul(xm,theta.T)
    nz=_np.matmul(xn,zeta.T)
    cosmt=_np.cos(mt)
    sinmt=_np.sin(mt)
    cosnz=_np.cos(nz)
    sinnz=_np.sin(nz)
    f = _np.zeros((ns,lt,lz))

    fmn = _np.ndarray((mn,lt))
    for k in range(ns):
        fmn = _np.broadcast_to(fmnc[k,:],(lt,mn)).T
        fmncosmt=(fmn*cosmt).T
        fmnsinmt=(fmn*sinmt).T
        f[k,:,:]=_np.matmul(fmncosmt, cosnz)-_np.matmul(fmnsinmt, sinnz)
    return f

def sfunct(theta,zeta,fmnc,xm,xn):
    f=0
    (ns,mn)=fmnc.shape
    lt = len(theta)
    lz = len(zeta)
    mt=_np.matmul(xm,theta.T)
    nz=_np.matmul(xn,zeta.T)
    cosmt=_np.cos(mt)
    sinmt=_np.sin(mt)
    cosnz=_np.cos(nz)
    sinnz=_np.sin(nz)
    f = _np.zeros((ns,lt,lz))
    fmn = _np.ndarray((mn,lt))
    for k in range(ns):
        fmn = _np.broadcast_to(fmnc[k,:],(lt,mn)).T
        f[k,:,:]=_np.matmul((fmn*sinmt).T,cosnz)+_np.matmul((fmn*cosmt).T,sinnz)
    return f

# ===================================================================== #
# ===================================================================== #
