#!/usr/bin/env python3
import sys, os
sys.path.insert(0, '../../pySTEL/')
import numpy as np                    #For Arrays
from math import pi
from libstell.beams3d import read_beams3d

try:
	qtCreatorPath=os.environ["STELLOPT_PATH"]
except KeyError:
	print("Please set environment variable STELLOPT_PATH")
	sys.exit(1)

lfail = 0
failtol = 5.0
filename='beams3d_ORBITS_depo.h5'
data=read_beams3d(filename)
if not data:
    print('ERROR Opening File: '+filename)
    sys.exit(0)

# Calc values

print('BEAMS3D VERSION: ' + str(round(data['VERSION'],2)))
print('==== Vectors ====')
varlist={}
varlist['Shinethrough']=np.array([34, 18, 10, 15, 4, 2])
#print(data['Shinethrough'])
for temp in varlist:
    act = varlist[temp]
    cal = data[temp]
    cal = np.where(act==0,0,cal)
    div = np.where(act==0,1,act)
    perct = max(abs(act-cal))
    print('  '+temp+': '+str(cal[0])+'   '+str(act[0])+'   '+str(int(perct))+'%')
    if perct > failtol:
        lfail = 1
print('=================')

if (lfail):
    print('  STATUS: FAIL!!!!!')
else:
    print('  STATUS: PASS')

sys.exit(0)




