#!/usr/bin/env python3
import sys, os
sys.path.insert(0, '../../pySTEL/')
import numpy as np                    #For Arrays
from math import pi
from libstell.beams3d import BEAMS3D

lfail = 0
failtol = 20.0
filename='beams3d_ORBITS_eqdsk.h5'
b3d = BEAMS3D()
b3d.read_beams3d(filename)

# Calc values
data = {}
data['Shinethrough'] = b3d.Shinethrough
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
varlist['Deposition_B1E1']=np.array([4.466062e+19,3.112300e+19,2.348636e+19,1.487319e+19,1.103912e+19, \
 8.094944e+18,6.257979e+18,5.400623e+18,4.874593e+18,4.197121e+18, \
 3.579773e+18,3.633305e+18,3.020388e+18,2.750214e+18,2.418788e+18, \
 9.718464e+17])
varlist['Deposition_B1E2']=np.array([7.714851e+19,5.602495e+19,3.913470e+19,2.883455e+19,2.080546e+19, \
 1.675960e+19,1.312089e+19,1.243215e+19,1.114408e+19,1.041841e+19, \
 9.067863e+18,7.988774e+18,7.539497e+18,7.576122e+18,6.276814e+18, \
 2.410528e+18])
varlist['Deposition_B1E3']=np.array([9.937372e+19,6.797368e+19,4.925968e+19,3.829194e+19,2.930141e+19, \
 2.275644e+19,2.072090e+19,1.771813e+19,1.748743e+19,1.591874e+19, \
 1.457257e+19,1.360660e+19,1.366905e+19,1.426168e+19,1.130261e+19, \
 5.104438e+18])
varlist['Deposition_B2E1']=np.array([4.587507e+19,3.445716e+19,2.316060e+19,1.485832e+19,1.041520e+19, \
 8.142330e+18,6.454796e+18,5.737334e+18,4.873494e+18,4.302033e+18, \
 3.765762e+18,3.134169e+18,3.052567e+18,2.869484e+18,2.405204e+18, \
 9.035431e+17])
varlist['Deposition_B2E2']=np.array([7.954348e+19,5.803747e+19,3.988500e+19,2.859940e+19,1.935425e+19, \
 1.736060e+19,1.342601e+19,1.171052e+19,1.034543e+19,1.011741e+19, \
 8.893924e+18,8.495731e+18,7.935270e+18,7.512586e+18,6.675309e+18, \
 2.552085e+18])
varlist['Deposition_B2E3']=np.array([9.830415e+19,6.944995e+19,5.248501e+19,3.670912e+19,2.849720e+19, \
 2.286985e+19,2.069946e+19,1.809512e+19,1.662198e+19,1.546932e+19, \
 1.463727e+19,1.396286e+19,1.425235e+19,1.363103e+19,1.200816e+19, \
 4.618210e+18])
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
