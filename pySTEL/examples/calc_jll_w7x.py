# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from read_vmec import read_vmec, calc_jll, cfunct, sfunct
import tempfile as tf
from osa import Client

vmec = Client('http://esb.ipp-hgw.mpg.de:8280/services/vmec_v7?wsdl')

#from pybaseutils.utils import cart2pol

#%% local file
#filename='/mnt/jons/03_Masterarbeit/data/w7x_ref_113/wout_w7x_ref_113.nc'

#%% read directly from webservice
#vmec_id="w7x_ref_82"
beta=[0.00, 0.16, 0.32, 0.48, 0.64, 0.80]
i=0
#for vmec_id in ['w7x_ref_82', 'w7x_ref_67', 'w7x_ref_68', 'w7x_ref_69', 'w7x_ref_83', 'w7x_ref_84']:
for vmec_id in ['w7x_ref_163']:
#    filename = "test.nc"
#    wout_temp = open(filename, 'w')
#    wout_bytes=vmec.service.getVmecOutputNetcdf(vmec_id)
#    newFileByteArray = bytearray(wout_bytes)
#    wout_temp.write(newFileByteArray)
#    filename=wout_temp.name
    
    filename="Z:\\03_Masterarbeit\\data\\"+vmec_id+"\\wout_"+vmec_id+".nc"
    
    #%% import data
    vmec_data=read_vmec(filename)
    
    nu = vmec_data['mpol']*4
    nv = vmec_data['ntor']*4*vmec_data['nfp']
    
    theta = np.ndarray((nu,1))
    zeta  = np.ndarray((nv,1))
    
    for j in range(nu):
        theta[j]=2*np.pi*j/(nu-1)
    for j in range(nv):
        zeta[j] =2*np.pi*j/(nv-1)
    
    r=cfunct(theta,zeta,vmec_data['rmnc'],vmec_data['xm'],vmec_data['xn'])
    z=sfunct(theta,zeta,vmec_data['zmns'],vmec_data['xm'],vmec_data['xn'])    
    
    jll=calc_jll(vmec_data, theta, zeta)
    
    
    #%% plotting
    # toroidal angle --> idx in zeta
    # 24 is here the triangular-shaped plane
    v=24 # triangular-shaped plane
    #v=0 # bean-shaped plane
    
    plt.figure()
    #plt.pcolor(r[:,:,v]-r[0,0,24],z[:,:,v]-z[0,0,24],jll[:,:,v],cmap='jet')
    plt.pcolor(r[:,:,v],z[:,:,v],jll[:,:,v],cmap='jet', vmin=-1e5, vmax=1e5)
#    plt.pcolor(r[:,:,v],z[:,:,v],jll[:,:,v],cmap='jet')
    #plt.imshow(jll[:,:,v], interpolation='none')
    plt.xlabel('r / m')
    plt.ylabel('z / m')
    plt.axes().set_aspect('equal')
    plt.title(vmec_id+", beta="+np.str(beta[i])+"%")
    plt.colorbar()
    plt.tight_layout()
    i=i+1


#    #radial  cross section at z~~0 plane
#    plt.figure()
#    plt.plot(r[:,0,v], jll[:,0,v], '.-',r[:,23,v], jll[:,23,v], '.-')
#    plt.grid()
#
#plt.figure()
#for num_fs in np.arange(1,np.shape(r)[0]):
#    [theta, rho]=cart2pol(r[num_fs,:,v]-r[0,0,24],z[num_fs,:,v]-z[0,0,24])
#    cont_theta=np.hstack([theta[24:], theta[1:24]])
#    cont_jll=np.hstack([jll[num_fs,24:,24], jll[num_fs,1:24,24]])
#    plt.plot(cont_theta*180.0/np.pi, cont_jll)
