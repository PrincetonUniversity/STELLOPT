#!/usr/bin/env python3
import sys, os
sys.path.insert(0, '../../pySTEL/')
import numpy as np                    #For Arrays
from math import pi
from libstell.beams3d import BEAMS3D

lfail = 0
failtol = 5.0
filename='beams3d_ORBITS_depo_adas.h5'
b3d = BEAMS3D()
b3d.read_beams3d(filename)

# Calc values
data = {}
[rho,depo]= b3d.calcDepo(ns=16)
data['Deposition_B1E1'] = depo[0,:]
data['Deposition_B1E2'] = depo[1,:]
data['Deposition_B1E3'] = depo[2,:]
data['Deposition_B2E1'] = depo[3,:]
data['Deposition_B2E2'] = depo[4,:]
data['Deposition_B2E3'] = depo[5,:]

print(f'BEAMS3D VERSION: {b3d.VERSION:4.2f}')
print('==== Vectors ====')
varlist={}
#varlist['Shinethrough']=np.array([34, 18, 10, 1, 0, 0])
varlist['Deposition_B1E1']=np.array([6.632940e+18,2.594600e+18,1.431116e+18,1.012614e+18,7.867336e+17, \
 6.228894e+17,4.855791e+17,4.210176e+17,3.390240e+17,3.075314e+17, \
 2.658876e+17,2.119868e+17,1.523728e+17,1.426544e+17,9.585534e+16, \
 5.917356e+16])
varlist['Deposition_B1E2']=np.array([1.671878e+19,6.162154e+18,3.473067e+18,2.389185e+18,1.718217e+18, \
 1.470851e+18,1.243820e+18,1.094551e+18,8.874908e+17,7.915353e+17, \
 6.467071e+17,5.569967e+17,4.268296e+17,3.553849e+17,2.420271e+17, \
 1.791050e+17])
varlist['Deposition_B1E3']=np.array([2.547626e+19,8.600564e+18,5.248852e+18,3.758316e+18,2.831966e+18, \
 2.494977e+18,1.845696e+18,1.730327e+18,1.528507e+18,1.287941e+18, \
 1.081973e+18,9.562760e+17,7.542674e+17,6.753299e+17,4.477863e+17, \
 2.704863e+17])
varlist['Deposition_B2E1']=np.array([1.502596e+19,3.208317e+18,1.637006e+18,1.215329e+18,9.440003e+17, \
 8.338557e+17,7.134530e+17,5.970941e+17,5.603523e+17,5.046326e+17, \
 4.198645e+17,3.591806e+17,2.852193e+17,2.213884e+17,1.769476e+17, \
 8.392739e+16])
varlist['Deposition_B2E2']=np.array([1.046748e+19,3.297453e+18,2.037332e+18,1.656553e+18,1.659661e+18, \
 1.577892e+18,1.382742e+18,1.425945e+18,1.295802e+18,1.216265e+18, \
 1.116363e+18,1.000230e+18,8.233047e+17,7.038711e+17,5.307589e+17, \
 2.932251e+17])
varlist['Deposition_B2E3']=np.array([5.945498e+18,2.521645e+18,1.877535e+18,1.797673e+18,1.923841e+18, \
 2.078492e+18,1.921390e+18,2.017907e+18,2.056794e+18,1.922589e+18, \
 1.900198e+18,1.785054e+18,1.552839e+18,1.315236e+18,1.034035e+18, \
 5.391188e+17])
for temp in varlist:
    act = varlist[temp]
    cal = data[temp]
    #print(np.array2string(cal,precision=6, separator=','))
    cal = np.where(act==0,0,cal)
    div = np.where(act==0,1,act)
    perct = 100*sum(abs(act-cal)/div)
    print(f'  Quantity: {temp} -- CODE -- REF. -- %')
    for i in range(len(act)):
        perct = 100*abs(act[i]-cal[i])/div[i]
        print(f'  {i} {cal[i]:7.6f} {act[i]:7.6f} {round(perct)}')
        if perct > failtol:
            lfail = 1
    print('=================')
if (lfail):
    print('  STATUS: FAIL!!!!!')
    sys.exit(-1)
else:
    print('  STATUS: PASS')
    sys.exit(0)
