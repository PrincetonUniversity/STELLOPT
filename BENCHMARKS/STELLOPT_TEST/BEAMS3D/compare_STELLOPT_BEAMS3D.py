#!/usr/bin/env python3
import sys, os
sys.path.insert(0, '../../../pySTEL/')
import numpy as np                    #For Arrays
from math import pi
from libstell.stellopt import read_stellopt

try:
	qtCreatorPath=os.environ["STELLOPT_PATH"]
except KeyError:
	print("Please set environment variable STELLOPT_PATH")
	sys.exit(1)

failtol = 5.0
filename='stellopt.BEAMS3D'
data=read_stellopt(filename)
if not data:
    print('ERROR Opening File: '+filename)
    sys.exit(0)

print('STELLOPT_VERSION: ' + str(data['VERSION']))
print('==== Scalars ====')
varlist={}
varlist['ORBIT_equil']=0.80675
lfail = 0;
for temp in varlist:
    act = varlist[temp]
    cal = data[temp]
    perct = 100*(abs(act-cal)/act)
    print('  '+temp+': '+str(cal)+'   '+str(act)+'   '+str(int(perct))+'%')
    if perct > failtol:
        lfail = 1
print('=================')

if (lfail):
    print('  STATUS: FAIL!!!!!')
else:
    print('  STATUS: PASS')

sys.exit(0)




