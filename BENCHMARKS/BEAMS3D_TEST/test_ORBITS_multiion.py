#!/usr/bin/env python3
import sys, os
sys.path.insert(0, '../../pySTEL/')
import numpy as np                    #For Arrays
from math import pi
from libstell.beams3d import BEAMS3D

lfail = 0
failtol = 20.0
filename='beams3d_ORBITS_multiion.h5'
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
#data['Deposition_B2E3'] = depo[5,:]

print(f'BEAMS3D VERSION: {b3d.VERSION:4.2f}')
print('==== Vectors ====')
varlist={}
varlist['Deposition_B1E1']=np.array([3.628550e+18,1.291153e+18,7.927005e+17,5.199417e+17,3.713347e+17, \
 3.248485e+17,2.580924e+17,2.271131e+17,1.955029e+17,1.562149e+17, \
 1.436643e+17,1.152797e+17,8.797583e+16,8.113917e+16,5.733855e+16, \
 3.808243e+16])
varlist['Deposition_B1E2']=np.array([7.995488e+18,2.699743e+18,1.678032e+18,1.212555e+18,9.676878e+17, \
 7.032958e+17,5.964166e+17,5.138963e+17,4.508982e+17,3.738844e+17, \
 3.310177e+17,2.764420e+17,2.098081e+17,2.022517e+17,1.479625e+17, \
 9.705334e+16])
varlist['Deposition_B1E3']=np.array([1.189420e+19,4.261046e+18,2.594906e+18,1.906508e+18,1.366930e+18, \
 1.136669e+18,9.516113e+17,8.312717e+17,7.041941e+17,6.424587e+17, \
 5.336555e+17,4.899983e+17,3.621583e+17,3.108720e+17,2.561323e+17, \
 1.811121e+17])
varlist['Deposition_B2E1']=np.array([6.381469e+18,1.467143e+18,7.435923e+17,5.444115e+17,4.284223e+17, \
 3.993574e+17,3.312299e+17,3.139192e+17,2.810265e+17,2.628457e+17, \
 2.255120e+17,1.821179e+17,1.655318e+17,1.430314e+17,1.088912e+17, \
 5.881516e+16])
varlist['Deposition_B2E2']=np.array([5.168567e+18,1.582075e+18,1.002512e+18,8.891264e+17,7.349174e+17, \
 7.461785e+17,7.214803e+17,6.348985e+17,6.180024e+17,5.751017e+17, \
 5.243465e+17,4.874632e+17,4.348598e+17,3.956598e+17,3.168090e+17, \
 1.766461e+17])
#varlist['Deposition_B2E3']=np.array([3.947099e+18,1.434560e+18,1.107841e+18,9.501735e+17,9.426179e+17, \
# 9.863313e+17,9.603728e+17,9.918007e+17,9.765742e+17,9.423943e+17, \
# 8.477654e+17,8.170676e+17,7.705660e+17,6.813209e+17,5.645839e+17, \
# 3.219246e+17])
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

