# -*- coding: utf-8 -*-
# 
# This module reads a VMEC wout file in netCDF format and calculates the currents.
# The code is heavily base on Sam Lazerson's libstell.py and the famous
# read_wout_mod.f by Steve Hirshman.
# 
# History:
#  2017-06-18/jons created file
#  2017-06-30/jons added stellarator-asymmetric current calculation
#
# corresponding author: Jonathan Schilling (jonathan.schilling@ipp.mpg.de)

from netCDF4 import Dataset
import numpy as np

mu0=4.0e-7*np.pi

# Calculate Currents
def calc_curr(f):
    f['currumnc']=np.zeros([f['ns'], f['mnmax_nyq']])
    f['currvmnc']=np.zeros([f['ns'], f['mnmax_nyq']])
    
    ohs = f['ns']-1.0
    hs  = 1.0/np.double(ohs)
    ns = f['ns']
    
    shalf=np.zeros(ns)
    sfull=np.zeros(ns)
    for i in np.arange(1,ns):
        shalf[i] = np.sqrt(hs*(i-0.5))
        sfull[i] = np.sqrt(hs*(i-0.0))
        
    js1 = np.arange(2,ns)
    js  = np.arange(1,ns-1)
    
    # calculation of curl B
    for mn in np.arange(f['mnmax_nyq']):
        if (np.mod(f['xm_nyq'][mn],2) == 1):
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
        f['currumnc'][js,mn] = -np.double(f['xn_nyq'][mn])*t1 - t3
        f['currvmnc'][js,mn] = -np.double(f['xm_nyq'][mn])*t1 + t2

    # linear extrapolation to axis
    for i in np.arange(f['mnmax_nyq']):
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
        f['currumns']=np.zeros([f['ns'], f['mnmax_nyq']])
        f['currvmns']=np.zeros([f['ns'], f['mnmax_nyq']])

        for mn in np.arange(f['mnmax_nyq']):
            if (np.mod(f['xm_nyq'][mn],2) == 1):
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
            f['currumns'][js,mn] = np.double(f['xn_nyq'][mn])*t1 - t3
            f['currvmns'][js,mn] = np.double(f['xm_nyq'][mn])*t1 + t2
            
        for i in np.arange(f['mnmax_nyq']):
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
    
# copied from libstell.py by SAL to circumvent libstell import here
def h2f(var,ns):
    import numpy as np
    temp = np.zeros((ns,1))
    temp[0] = 1.5 * var[0] - 0.5 * var[1]
    temp[1:ns-1] = 0.5* (var[0:ns-2] + var[1:ns-1])
    temp[ns-1] = 1.5 * var[ns-2] - 0.5 * var[ns-3]
    return temp


def cfunct(theta,zeta,fmnc,xm,xn):
    import numpy as np
    f=0
    (ns,mn)=fmnc.shape
    lt = len(theta)
    lz = len(zeta)
    mt=np.matmul(xm,theta.T)
    nz=np.matmul(xn,zeta.T)
    cosmt=np.cos(mt)
    sinmt=np.sin(mt)
    cosnz=np.cos(nz)
    sinnz=np.sin(nz)
    f = np.zeros((ns,lt,lz))
    
    fmn = np.ndarray((mn,lt))
    for k in range(ns):
        fmn = np.broadcast_to(fmnc[k,:],(lt,mn)).T
        fmncosmt=(fmn*cosmt).T
        fmnsinmt=(fmn*sinmt).T
        f[k,:,:]=np.matmul(fmncosmt, cosnz)-np.matmul(fmnsinmt, sinnz)
    return f
    
def sfunct(theta,zeta,fmnc,xm,xn):
    import numpy as np
    f=0
    (ns,mn)=fmnc.shape
    lt = len(theta)
    lz = len(zeta)
    mt=np.matmul(xm,theta.T)
    nz=np.matmul(xn,zeta.T)
    cosmt=np.cos(mt)
    sinmt=np.sin(mt)
    cosnz=np.cos(nz)
    sinnz=np.sin(nz)
    f = np.zeros((ns,lt,lz))
    fmn = np.ndarray((mn,lt))
    for k in range(ns):
        fmn = np.broadcast_to(fmnc[k,:],(lt,mn)).T
        f[k,:,:]=np.matmul((fmn*sinmt).T,cosnz)+np.matmul((fmn*cosmt).T,sinnz)
    return f
    

def calc_jll(vmec_data, theta, zeta ):
    # CALC_JLL(vmec_data,theta,zeta) Calculates the parallel current density.
    # This funciton takes a VMEC data structure (as read by read_vmec) and
    # theta/zeta arrays as input and outputs the parallel current density.
    
    # Example usage (Matlab)
    #      theta=0:2*pi/359:2*pi;
    #      zeta=0:2*pi/63:2*pi;
    #      data=read_vmec('wout.test');        % Reads VMEC wout file
    #      jll=calc_jll(vmec_data,theta,zeta); % Calculate the current
    #
    # Example usage (Python)
    #      theta=np.linspace(0, 2*np.pi, 360)
    #      zeta=np.linspace(0, 2*np.pi, 64)
    #      vmec_data=read_vmec('wout.nc')
    #      jll=calc_jll(vmec_data, theta, zeta)
    
    
    # Maintained by: Samuel Lazerson (lazerson@pppl.gov)
    # Ported to python by Jonathan Schilling (jonathan.schilling@ipp.mpg.de)
    # Version:       1.00
    
    b =cfunct(theta,zeta,vmec_data['bmnc'],    vmec_data['xm'],vmec_data['xn'])
    g =cfunct(theta,zeta,vmec_data['gmnc'],    vmec_data['xm'],vmec_data['xn'])
    bu=cfunct(theta,zeta,vmec_data['bsubumnc'],vmec_data['xm'],vmec_data['xn'])
    bv=cfunct(theta,zeta,vmec_data['bsubvmnc'],vmec_data['xm'],vmec_data['xn'])
    ju=cfunct(theta,zeta,vmec_data['currumnc'],vmec_data['xm'],vmec_data['xn'])
    jv=cfunct(theta,zeta,vmec_data['currvmnc'],vmec_data['xm'],vmec_data['xn'])
    
    if (vmec_data['iasym']):
        b =b +sfunct(theta,zeta,vmec_data['bmns'],    vmec_data['xm'],vmec_data['xn'])
        g =g +sfunct(theta,zeta,vmec_data['gmns'],    vmec_data['xm'],vmec_data['xn'])
        bu=bu+sfunct(theta,zeta,vmec_data['bsubumns'],vmec_data['xm'],vmec_data['xn'])
        bv=bv+sfunct(theta,zeta,vmec_data['bsubvmns'],vmec_data['xm'],vmec_data['xn'])
        ju=ju+sfunct(theta,zeta,vmec_data['currumns'],vmec_data['xm'],vmec_data['xn'])
        jv=jv+sfunct(theta,zeta,vmec_data['currvmns'],vmec_data['xm'],vmec_data['xn'])
    
    
    jll = (bu*ju+bv*jv)/(g*b)
    return jll


    
def read_vmec(filename):
    rootgrp = Dataset(filename, 'r')
    vmec_data=dict()

    vmec_data['ns']=rootgrp['/ns'][0]
    vmec_data['nfp']=rootgrp['/nfp'][0]
    vmec_data['mpol']=rootgrp['/mpol'][0]
    vmec_data['ntor']=rootgrp['/ntor'][0]
    vmec_data['mnmax']=rootgrp['/mnmax'][0]    
    vmec_data['mnmax_nyq']=rootgrp['/mnmax_nyq'][0]    
    vmec_data['iasym']=rootgrp['/lasym__logical__'][0]
    vmec_data['ierr_vmec']=rootgrp['/ier_flag'][0]
    vmec_data['wb']=rootgrp['/wb'][0]
    vmec_data['wp']=rootgrp['/wp'][0]
    vmec_data['gamma']=rootgrp['/gamma'][0]
    vmec_data['rmax_surf']=rootgrp['/rmax_surf'][0]
    vmec_data['rmin_surf']=rootgrp['/rmin_surf'][0]
    vmec_data['zmax_surf']=rootgrp['/zmax_surf'][0]
    vmec_data['aspect']=rootgrp['/aspect'][0]
    vmec_data['betatot']=rootgrp['/betatotal'][0]
    vmec_data['betapol']=rootgrp['/betapol'][0]
    vmec_data['betator']=rootgrp['/betator'][0]
    vmec_data['betaxis']=rootgrp['/betaxis'][0]
    vmec_data['b0']=rootgrp['/b0'][0]
    vmec_data['version']=rootgrp['/version_'][0]
    vmec_data['IonLarmor']=rootgrp['/IonLarmor'][0]
    vmec_data['VolAvgB']=rootgrp['/volavgB'][0]
    #vmec_data['fsql']=rootgrp['/fsql'][0]
    #vmec_data['fsqr']=rootgrp['/fsqr'][0]
    #vmec_data['fsqz']=rootgrp['/fsqz'][0]
    #vmec_data['ftolv']=rootgrp['/ftolv'][0]
    vmec_data['Aminor']=rootgrp['/Aminor_p'][0]
    vmec_data['Rmajor']=rootgrp['/Rmajor_p'][0]
    vmec_data['Volume']=rootgrp['/volume_p'][0]
    vmec_data['RBtor']=rootgrp['/rbtor'][0]
    vmec_data['RBtor0']=rootgrp['/rbtor0'][0]
    vmec_data['Itor']=rootgrp['/ctor'][0]
    vmec_data['isigng']=rootgrp['/signgs'][0] # sign of sqrt(g)

    # since machsq and pfac is declared in read_wout_mod.f, it is assigned a value
    # at the fortran-to-python interface in libstell.
    # Therefore, create this variable here as well even if no data is present.
    for key in ['machsq', 'pfac']:
        if key in rootgrp.variables.keys():
            vmec_data[key]=rootgrp['/'+key][0]
        else:
            vmec_data[key]=0.0

    # Array values 1D
    ns = vmec_data['ns']
    mnmax = vmec_data['mnmax']
    mnmax_nyq = vmec_data['mnmax_nyq']
    ns_size = (ns,1)
    mn_size = (mnmax,1)
    mnnyq_size = (mnmax_nyq,1)
    vmec_data['iotas']  =np.reshape(rootgrp['/iotas'][:], ns_size)
    vmec_data['iotaf']  =np.reshape(rootgrp['/iotaf'][:], ns_size)
    vmec_data['presf']  =np.reshape(rootgrp['/presf'][:], ns_size)
    vmec_data['phipf']  =np.reshape(rootgrp['/phipf'][:], ns_size)
    #vmec_data['qfact'] =np.reshape(rootgrp['/qfact'][:], ns_size) ## was commented out in libstell.py as well
    #vmec_data['chipf'] =np.reshape(rootgrp['/chipf'][:], ns_size) ## random data ???
    #vmec_data['chi']   =np.reshape(rootgrp['/chi'][:], ns_size) ## random data ???
    vmec_data['phi']    =np.reshape(rootgrp['/phi'][:], ns_size)
    vmec_data['mass']   =np.reshape(rootgrp['/mass'][:], ns_size)
    vmec_data['pres']   =np.reshape(rootgrp['/pres'][:], ns_size)
    vmec_data['betavol']=np.reshape(rootgrp['/beta_vol'][:], ns_size)
    vmec_data['xm']     =np.reshape(rootgrp['/xm'][:], mn_size)
    vmec_data['xn']     =np.reshape(rootgrp['/xn'][:], mn_size)
    vmec_data['xm_nyq'] =np.reshape(rootgrp['/xm_nyq'][:], mnnyq_size)
    vmec_data['xn_nyq'] =np.reshape(rootgrp['/xn_nyq'][:], mnnyq_size)
    vmec_data['phip']   =np.reshape(rootgrp['/phips'][:], ns_size)
    vmec_data['buco']   =np.reshape(rootgrp['/buco'][:], ns_size)
    vmec_data['bvco']   =np.reshape(rootgrp['/bvco'][:], ns_size)
    vmec_data['vp']     =np.reshape(rootgrp['/vp'][:], ns_size)
    #vmec_data['overr'] =np.reshape(rootgrp['/over_r'][:], ns_size) ## random data ???
    vmec_data['jcuru']  =np.reshape(rootgrp['/jcuru'][:], ns_size)
    vmec_data['jcurv']  =np.reshape(rootgrp['/jcurv'][:], ns_size)
    vmec_data['specw']  =np.reshape(rootgrp['/specw'][:], ns_size)
    vmec_data['jdotb']  =np.reshape(rootgrp['/jdotb'][:], ns_size)
    vmec_data['Dmerc']  =np.reshape(rootgrp['/DMerc'][:], ns_size)
    vmec_data['Dshear'] =np.reshape(rootgrp['/DShear'][:], ns_size)
    vmec_data['Dwell']  =np.reshape(rootgrp['/DWell'][:], ns_size)
    vmec_data['Dcurr']  =np.reshape(rootgrp['/DCurr'][:], ns_size)
    vmec_data['Dgeod']  =np.reshape(rootgrp['/DGeod'][:], ns_size)
    vmec_data['equif']  =np.reshape(rootgrp['/equif'][:], ns_size)

    for key in ['qfact', 'chipf', 'chi', 'overr']:
        if key in rootgrp.variables.keys():
            vmec_data[key]=np.reshape(rootgrp['/'+key][:], ns_size)
        else:
            vmec_data[key]=np.zeros(ns_size)
    
    vmec_data['rmnc']=rootgrp['/rmnc'][:,:]
    vmec_data['zmns']=rootgrp['/zmns'][:,:]
        
    vmec_data['lmns']=rootgrp['/lmns'][:,:]
    vmec_data['bmnc']=rootgrp['/bmnc'][:,:]
    vmec_data['gmnc']=rootgrp['/gmnc'][:,:] 
    
    vmec_data['bsupumnc']=rootgrp['/bsupumnc'][:,:]
    vmec_data['bsupvmnc']=rootgrp['/bsupvmnc'][:,:]
    vmec_data['bsubsmns']=rootgrp['/bsubsmns'][:,:]
    vmec_data['bsubumnc']=rootgrp['/bsubumnc'][:,:]
    vmec_data['bsubvmnc']=rootgrp['/bsubvmnc'][:,:]
    
    # not tested yet...    
    if (vmec_data['iasym']):
        vmec_data['rmns']=rootgrp['/rmns'][:,:]
        vmec_data['zmnc']=rootgrp['/zmnc'][:,:]
        
        vmec_data['lmnc']=rootgrp['/lmnc'][:,:]
        vmec_data['bmns']=rootgrp['/bmns'][:,:]
        vmec_data['gmns']=rootgrp['/gmns'][:,:]
        
        vmec_data['bsupumns']=rootgrp['/bsupumns'][:,:]
        vmec_data['bsupvmns']=rootgrp['/bsupvmns'][:,:]
        vmec_data['bsubsmnc']=rootgrp['/bsubsmnc'][:,:]
        vmec_data['bsubumns']=rootgrp['/bsubumns'][:,:]
        vmec_data['bsubvmns']=rootgrp['/bsubvmns'][:,:]
    else:
        vmec_data['rmns']=None
        vmec_data['zmnc']=None
        
        vmec_data['lmnc']=None
        vmec_data['bmns']=None
        vmec_data['gmns']=None
        
        vmec_data['bsupumns']=None
        vmec_data['bsupvmns']=None
        vmec_data['bsubsmnc']=None
        vmec_data['bsubumns']=None
        vmec_data['bsubvmns']=None
    
    # we are finished reading the file
    rootgrp.close()
    
    # calc currumnc, currvmnc
    # important: do this _BEFORE_ you transform the half-mesh quantities onto
    # the full mesh !!!!!!!!!
    currents=calc_curr(vmec_data)
    vmec_data['currumnc']=currents['currumnc']
    vmec_data['currvmnc']=currents['currvmnc']
    if (vmec_data['iasym']):
        vmec_data['currumns']=currents['currumns']
        vmec_data['currvmns']=currents['currvmns']
    else:
        vmec_data['currumns']=None
        vmec_data['currvmns']=None
        

    # new definition => adapt to NESCOIL ???
#    vmec_data['xn'] = -vmec_data['xn']

#    # Put on full grid
#    vmec_data['buco'] = h2f(vmec_data['buco'],ns)
#    vmec_data['bvco'] = h2f(vmec_data['bvco'],ns)
#    vmec_data['vp'] = h2f(vmec_data['vp'],ns)
#    vmec_data['overr'] = h2f(vmec_data['overr'],ns)
#    vmec_data['specw'] = h2f(vmec_data['specw'],ns)
    
#    # Put matrix quantities on full grid
#    for key in ['bmnc','gmnc','lmns','bsupumnc','bsupvmnc','bsubsmns','bsubumnc','bsubvmnc']:
#        vmec_data[key][0,:] = 1.5 * vmec_data[key][1,:] - 0.5 * vmec_data[key][2,:]
#        vmec_data[key][1:ns-2,:] = 0.5 * (vmec_data[key][1:ns-2,:] + vmec_data[key][2:ns-1,:])
#        vmec_data[key][ns-1,:] = 2.0 * vmec_data[key][ns-2,:] - vmec_data[key][ns-3,:]    
#    if vmec_data['iasym']:
#        for key in ['bmns','gmns','lmnc','bsupumns','bsupvmns','bsubsmnc','bsubumns','bsubvmns']:
#            vmec_data[key][0,:] = 1.5 * vmec_data[key][1,:] - 0.5 * vmec_data[key][2,:]
#            vmec_data[key][1:ns-2,:] = 0.5 * (vmec_data[key][1:ns-2,:] + vmec_data[key][2:ns-1,:])
#            vmec_data[key][ns-1,:] = 2.0 * vmec_data[key][ns-2,:] - vmec_data[key][ns-3,:]

    
    return vmec_data

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
    vmec_data['buco'] = h2f(vmec_data['buco'],ns)
    vmec_data['bvco'] = h2f(vmec_data['bvco'],ns)
    vmec_data['vp'] = h2f(vmec_data['vp'],ns)
    vmec_data['overr'] = h2f(vmec_data['overr'],ns)
    vmec_data['specw'] = h2f(vmec_data['specw'],ns)
    
    # Put matrix quantities on full grid
    for key in ['bmnc','gmnc','lmns','bsupumnc','bsupvmnc','bsubsmns','bsubumnc','bsubvmnc']:
        vmec_data[key][0,:] = 1.5 * vmec_data[key][1,:] - 0.5 * vmec_data[key][2,:]
        vmec_data[key][1:ns-2,:] = 0.5 * (vmec_data[key][1:ns-2,:] + vmec_data[key][2:ns-1,:])
        vmec_data[key][ns-1,:] = 2.0 * vmec_data[key][ns-2,:] - vmec_data[key][ns-3,:]    
    if vmec_data['iasym']:
        for key in ['bmns','gmns','lmnc','bsupumns','bsupvmns','bsubsmnc','bsubumns','bsubvmns']:
            vmec_data[key][0,:] = 1.5 * vmec_data[key][1,:] - 0.5 * vmec_data[key][2,:]
            vmec_data[key][1:ns-2,:] = 0.5 * (vmec_data[key][1:ns-2,:] + vmec_data[key][2:ns-1,:])
            vmec_data[key][ns-1,:] = 2.0 * vmec_data[key][ns-2,:] - vmec_data[key][ns-3,:]
            
         
            
            
            
            
            
            
            
            
            
