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
filename='beams3d_ORBITS_loss.h5'
b3d = BEAMS3D()
b3d.read_beams3d(filename)

# Calc values
data={}
data['end_state'] = b3d.end_state

print(f'BEAMS3D VERSION: {b3d.VERSION:4.2f}')
print('==== Vectors ====')
varlist={}
varlist['end_state']=np.array([0, 0, 0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2])
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
