#!/usr/bin/env python3
import sys, os
sys.path.insert(0, '../../pySTEL/')
import numpy as np
from math import pi
from libstell.beams3d import BEAMS3D

lfail = 0
failtol = 1.0
filename='beams3d_ORBITS_er.h5'
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
varlist['turning']=np.array([0.731858,0.790553,0.850824,0.912077,0.974712,1.038719,1.104833,1.172821, \
 1.243094,1.316038,1.392032,1.471591,1.555411,1.644194,1.739137,1.841611, \
 1.953852,2.07912 ,2.223097,2.396782,2.630352,3.141335,3.137608,3.139102, \
 3.14006 ,3.131666,3.140788,3.139625,3.139039,3.12906 ,3.13544 ,3.131839, \
 3.140812,3.141167,3.141409,3.135179,3.141481,3.140862,3.140919,3.07051 ])
varlist['delta']= np.array([0.054275,0.058177,0.062093,0.066008,0.069844,0.073672,0.077491,0.081266, \
 0.085027,0.088759,0.092453,0.096132,0.099778,0.103395,0.106995,0.110572, \
 0.114112,0.11766 ,0.121189,0.124668,0.128161,0.054815,0.046612,0.042663, \
 0.039882,0.037724,0.035944,0.034453,0.033163,0.032031,0.031029,0.030138, \
 0.029342,0.0286  ,0.027936,0.02734 ,0.026778,0.02624 ,0.025781,0.025343])
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
