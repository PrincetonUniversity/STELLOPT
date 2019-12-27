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

failtol = 1.0
filename='stellopt.TOK_R0_RHO'
data=read_stellopt(filename)
if not data:
    print('ERROR Opening File: '+filename)
    sys.exit(0)
else:
    print('EXTENSION: '+filename)
lfail = 0;
print('==== Scalars ====')
varlist={}
varlist['R0_equil']=10.0
varlist['ASPECT_equil']=10.0
n = data['R0_equil'].shape
for temp in varlist:
    act = varlist[temp]
    cal = data[temp][n[0]-1]
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




