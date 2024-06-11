#!/usr/bin/env python3
import sys, os
sys.path.insert(0, '../../pySTEL/')
import numpy as np                    #For Arrays
from math import pi
from libstell.beams3d import BEAMS3D

lfail = 0
failtol = 25.0
filename='beams3d_ORBITS_depo_adas.h5'
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
#data['Deposition_B2E3'] = depo[5,:]

print(f'BEAMS3D VERSION: {b3d.VERSION:4.2f}')
print('==== Vectors ====')
varlist={}
#varlist['Shinethrough']=np.array([34, 18, 10, 1, 0, 0])
varlist['Deposition_B1E1']=np.array([7.601978e+18,2.676854e+18,1.370818e+18,1.029026e+18,8.044945e+17, \
 6.280154e+17,5.260080e+17,4.681687e+17,3.531504e+17,3.015005e+17, \
 2.513114e+17,1.924037e+17,1.360968e+17,1.220854e+17,7.583060e+16, \
 3.850092e+16])
varlist['Deposition_B1E2']=np.array([1.656654e+19,5.704751e+18,3.479683e+18,2.569270e+18,1.940054e+18, \
 1.580344e+18,1.264489e+18,1.027386e+18,8.501639e+17,6.994295e+17, \
 5.771098e+17,4.904445e+17,3.204218e+17,2.913479e+17,1.888035e+17, \
 9.843338e+16])
varlist['Deposition_B1E3']=np.array([2.532364e+19,9.063671e+18,5.456691e+18,3.941292e+18,3.136900e+18, \
 2.295321e+18,1.877845e+18,1.718707e+18,1.445436e+18,1.231316e+18, \
 9.485466e+17,7.910622e+17,5.996606e+17,5.189052e+17,3.057828e+17, \
 1.927655e+17])
varlist['Deposition_B2E1']=np.array([1.482073e+19,3.492539e+18,1.603240e+18,1.346279e+18,1.010841e+18, \
 8.411326e+17,7.188941e+17,6.464879e+17,5.858530e+17,5.138584e+17, \
 4.321975e+17,3.360893e+17,2.452664e+17,2.068648e+17,1.439081e+17, \
 5.480686e+16])
varlist['Deposition_B2E2']=np.array([1.485353e+19,3.945685e+18,2.543646e+18,2.066734e+18,1.906469e+18, \
 1.786775e+18,1.622447e+18,1.436244e+18,1.270898e+18,1.181512e+18, \
 1.083766e+18,8.805631e+17,6.934838e+17,5.633447e+17,4.008583e+17, \
 1.511660e+17])
# varlist['Deposition_B2E3']=np.array([5.945498e+18,2.521645e+18,1.877535e+18,1.797673e+18,1.923841e+18, \
#  2.078492e+18,1.921390e+18,2.017907e+18,2.056794e+18,1.922589e+18, \
#  1.900198e+18,1.785054e+18,1.552839e+18,1.315236e+18,1.034035e+18, \
#  5.391188e+17])
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
