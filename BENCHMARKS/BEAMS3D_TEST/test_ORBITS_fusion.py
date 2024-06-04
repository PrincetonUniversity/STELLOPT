#!/usr/bin/env python3
import sys, os
sys.path.insert(0, '../../pySTEL/')
import numpy as np                    #For Arrays
from math import pi
from libstell.beams3d import BEAMS3D

lfail = 0
failtol = 5.0
filename='beams3d_ORBITS_fusion.h5'
b3d = BEAMS3D()
b3d.read_beams3d(filename)

# Calc values
data = {}
[rho,depo]= b3d.calcDepo(ns=16)
data['Birth_He4'] = depo[0,:]
data['Birth_T'] = depo[1,:]
data['Birth_p'] = depo[2,:]
data['Birth_He3'] = depo[3,:]

print(f'BEAMS3D VERSION: {b3d.VERSION:4.2f}')
print('==== Vectors ====')
varlist={}
#varlist['Shinethrough']=np.array([34, 18, 10, 1, 0, 0])
varlist['Birth_He4']=np.array([2.812140e+17,2.751405e+17,2.590578e+17,2.410255e+17,2.184702e+17, \
 1.950605e+17,1.718746e+17,1.458198e+17,1.168888e+17,8.445555e+16, \
 5.107926e+16,2.305566e+16,6.613366e+15,8.638698e+14,2.686111e+13, \
 3.171751e+10])
varlist['Birth_T']=np.array([7.196400e+14,7.022415e+14,6.655436e+14,6.156879e+14,5.659616e+14, \
 5.110290e+14,4.560634e+14,3.966889e+14,3.283568e+14,2.468759e+14, \
 1.577814e+14,7.708903e+13,2.434126e+13,3.671389e+12,1.364800e+11, \
 2.129394e+08])
varlist['Birth_p']=np.array([7.207028e+14,7.003268e+14,6.652009e+14,6.173576e+14,5.658323e+14, \
 5.104955e+14,4.566769e+14,3.965443e+14,3.276335e+14,2.474383e+14, \
 1.575486e+14,7.711104e+13,2.442342e+13,3.666051e+12,1.374268e+11, \
 2.163493e+08])
varlist['Birth_He3']=np.array([7.512938e+14,7.290716e+14,6.912372e+14,6.414480e+14,5.857319e+14, \
 5.310180e+14,4.701472e+14,4.080294e+14,3.368531e+14,2.525191e+14, \
 1.600640e+14,7.802688e+13,2.456785e+13,3.648097e+12,1.358381e+11, \
 2.130947e+08])
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




