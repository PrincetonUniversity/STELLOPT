#!/usr/bin/env python3
import sys, os
sys.path.insert(0, '../../pySTEL/')
import numpy as np                    #For Arrays
from math import pi
from libstell.beams3d import BEAMS3D

try:
	qtCreatorPath=os.environ["STELLOPT_PATH"]
except KeyError:
	print("Please set environment variable STELLOPT_PATH")
	sys.exit(1)

lfail = 0
failtol = 1.0
filename='beams3d_ORBITS.h5'
b3d = BEAMS3D()
b3d.read_beams3d(filename)

# Calc values
rho = np.sqrt(b3d.S_lines)
rho_max = np.max(rho,axis=1)
rho_min = np.min(rho,axis=1)
data = {}
data['delta'] = rho_max-rho_min
x = b3d.R_lines - 10.0
y = b3d.Z_lines
theta = np.arctan2(y,x)
theta = np.where(theta > np.pi,theta-pi,theta)
data['turning'] = np.max(theta,axis=1)

print(f'BEAMS3D VERSION: {b3d.VERSION:4.2f}')
print('==== Vectors ====')
varlist={}
varlist['turning']=np.array([0.122773,0.178807,0.234767,0.290948,0.347227,0.403768,0.460658,0.517991, \
 0.575882,0.634417,0.69372 ,0.753918,0.815149,0.877555,0.941314,1.006603, \
 1.07365 ,1.142698,1.214049,1.288031,1.365077,1.445688,1.530542,1.620162, \
 1.716398,1.819985,1.933357,2.05977 ,2.204824,2.379497,2.613095,3.140659, \
 3.140704,3.141489,3.129374,3.140686,3.130009,3.136217,3.140605,3.140605])
varlist['delta']= np.array([0.009462,0.01375 ,0.018008,0.022235,0.026432,0.030597,0.034731,0.038832, \
 0.042898,0.046932,0.05093 ,0.054895,0.058827,0.062721,0.066581,0.070407, \
 0.074195,0.07795 ,0.081667,0.085351,0.089   ,0.092616,0.096198,0.099746, \
 0.103264,0.106749,0.110207,0.113636,0.117039,0.120418,0.123759,0.058065, \
 0.049486,0.045529,0.04275 ,0.040595,0.038827,0.037335,0.036049,0.034923])
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




