#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 21:55:22 2017

@author: jonathan
"""

from libstell.libstell import read_vmec, calc_jll, cfunct, sfunct
import numpy as np

import matplotlib.pyplot as plt

vmec_data=read_vmec('wout_w7x_ref_60.nc')

nu = vmec_data['mpol']*4
nv = vmec_data['ntor']*4*vmec_data['nfp']

theta = np.ndarray((nu,1))
zeta  = np.ndarray((nv,1))

for j in range(nu):
    theta[j]=2*np.pi*j/(nu-1)
for j in range(nv):
    zeta[j]=2*np.pi*j/((nv-1))
    
mu0 = 4*np.pi*1e-7

jll = calc_jll(vmec_data, theta, zeta)*mu0

r=cfunct(theta,zeta,vmec_data['rmnc'],vmec_data['xm'],vmec_data['xn'])
z=sfunct(theta,zeta,vmec_data['zmns'],vmec_data['xm'],vmec_data['xn'])

# toroidal angle --> idx in zeta
# 24 is here the triangular-shaped plane
v=24

plt.pcolormesh(r[:,:,v],z[:,:,v],jll[:,:,v],cmap='jet',shading='gouraud')
plt.xlabel('r / m')
plt.ylabel('z / m')
plt.axes().set_aspect('equal')
plt.colorbar()