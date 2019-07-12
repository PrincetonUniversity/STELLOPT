# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
mu0=4*np.pi*1e-7


# copied from libstell.py by SAL to circumvent libstell import here
def h2f(var,ns):
    import numpy as np
    temp = np.zeros((ns,1))
    temp[0] = 1.5 * var[0] - 0.5 * var[1]
    temp[1:ns-1] = 0.5* (var[0:ns-2] + var[1:ns-1])
    temp[ns-1] = 1.5 * var[ns-2] - 0.5 * var[ns-3]
    return temp



filename='/home/IPP-HGW/jons/wout_w7x_ref_60.nc'

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
vmec_data['fsql']=rootgrp['/fsql'][0]
vmec_data['fsqr']=rootgrp['/fsqr'][0]
vmec_data['fsqz']=rootgrp['/fsqz'][0]
vmec_data['ftolv']=rootgrp['/ftolv'][0]
vmec_data['Aminor']=rootgrp['/Aminor_p'][0]
vmec_data['Rmajor']=rootgrp['/Rmajor_p'][0]
vmec_data['Volume']=rootgrp['/volume_p'][0]
vmec_data['RBtor']=rootgrp['/rbtor'][0]
vmec_data['RBtor0']=rootgrp['/rbtor0'][0]
vmec_data['Itor']=rootgrp['/ctor'][0]

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

# we are finished reading the file
rootgrp.close()

# calc currumnc, currvmnc
# important: do this _BEFORE_ you transform the half-mesh quantities onto
# the full mesh !!!!!!!!!
# Calculate Currents
vmec_data['currumnc']=np.zeros([vmec_data['ns'], vmec_data['mnmax_nyq']])
vmec_data['currvmnc']=np.zeros([vmec_data['ns'], vmec_data['mnmax_nyq']])

ohs = vmec_data['ns']-1.0
hs  = 1.0/np.double(ohs)
ns = vmec_data['ns']

shalf=np.zeros(ns)
sfull=np.zeros(ns)
for i in np.arange(1,ns):
    shalf[i] = np.sqrt(hs*(i-0.5))
    sfull[i] = np.sqrt(hs*(i-0.0))
    
js1 = np.arange(2,ns)
js  = np.arange(1,ns-1)

for mn in np.arange(vmec_data['mnmax_nyq']):
    if (np.mod(vmec_data['xm_nyq'][mn],2) == 1):
        #print("eq 1 at "+np.str(mn))
        t1  = 0.5*(shalf[js1] * vmec_data['bsubsmns'][js1,mn] + 
                   shalf[js]  * vmec_data['bsubsmns'][js,mn]   ) / sfull[js]
        bu0 = vmec_data['bsubumnc'][js,mn]/shalf[js]
        bu1 = vmec_data['bsubumnc'][js1,mn]/shalf[js1]
        t2  = ohs*(bu1-bu0)*sfull[js]+0.25*(bu0+bu1)/sfull[js]
        
        bv0 = vmec_data['bsubvmnc'][js,mn]/shalf[js]
        bv1 = vmec_data['bsubvmnc'][js1,mn]/shalf[js1]
        t3  = ohs*(bv1-bv0)*sfull[js]+0.25*(bv0+bv1)/sfull[js]
    else:
        t1  = 0.5*(vmec_data['bsubsmns'][js1,mn]+vmec_data['bsubsmns'][js,mn])
        t2  = ohs*(vmec_data['bsubumnc'][js1,mn]-vmec_data['bsubumnc'][js,mn])
        t3  = ohs*(vmec_data['bsubvmnc'][js1,mn]-vmec_data['bsubvmnc'][js,mn])
    vmec_data['currumnc'][js,mn] = -np.double(vmec_data['xn_nyq'][mn])*t1 - t3
    vmec_data['currvmnc'][js,mn] = -np.double(vmec_data['xm_nyq'][mn])*t1 + t2
    
    
    #print("currumnc[1]["+np.str(mn)+"] = "+np.str(vmec_data['currumnc'][1,mn]/mu0))
    
for i in np.arange(vmec_data['mnmax_nyq']):
    if (vmec_data['xm_nyq'][i]<=1):
        vmec_data['currumnc'][0,i]=2.0*vmec_data['currumnc'][1,i]-vmec_data['currumnc'][2,i]
        vmec_data['currvmnc'][0,i]=2.0*vmec_data['currvmnc'][1,i]-vmec_data['currvmnc'][2,i]
    else:
        vmec_data['currumnc'][0,i]=0.0
        vmec_data['currvmnc'][0,i]=0.0

vmec_data['currumnc'][vmec_data['ns']-1,:]=2.0*vmec_data['currumnc'][vmec_data['ns']-2,:]-vmec_data['currumnc'][vmec_data['ns']-3,:]
vmec_data['currvmnc'][vmec_data['ns']-1,:]=2.0*vmec_data['currvmnc'][vmec_data['ns']-2,:]-vmec_data['currvmnc'][vmec_data['ns']-3,:]

vmec_data['currumnc']=vmec_data['currumnc']/mu0;
vmec_data['currvmnc']=vmec_data['currvmnc']/mu0;


# new definition
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
        
np.savetxt("/mnt/jons/03_Masterarbeit/currumnc_libstell_ref60.dat", vmec_data['currumnc'], fmt="%20.20f")