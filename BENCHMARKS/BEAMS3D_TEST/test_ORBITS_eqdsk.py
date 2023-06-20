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
failtol = 10.0
filename='beams3d_ORBITS_eqdsk.h5'
data=read_beams3d(filename)
if not data:
    print('ERROR Opening File: '+filename)
    sys.exit(0)

# Calc values

print('BEAMS3D VERSION: ' + str(round(data['VERSION'],2)))
print('==== Vectors ====')
varlist={}
varlist['Shinethrough']=np.array([18.543322,13.09413 ,11.835307,19.003081,13.060472,12.109276])
#print(data['Shinethrough'])
for temp in varlist:
    act = varlist[temp]
    cal = data[temp]
    #print(np.array2string(cal,precision=6, separator=','))
    cal = np.where(act==0,0,cal)
    div = np.where(act==0,1,act)
    perct = max(abs(act-cal))
    print('  '+temp+': '+str(cal[0])+'   '+str(act[0])+'   '+str(round(perct))+'%')
    if perct > failtol:
        lfail = 1
print('=================')

if (lfail):
    print('  STATUS: FAIL!!!!!')
else:
    print('  STATUS: PASS')

sys.exit(0)




