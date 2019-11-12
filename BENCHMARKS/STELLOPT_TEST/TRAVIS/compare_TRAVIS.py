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
filename='stellopt.TRAVIS'
data=read_stellopt(filename)
if not data:
    print('ERROR Opening File: '+filename)
    sys.exit(0)

print('STELLOPT_VERSION: ' + str(data['VERSION']))


print('==== Vectors ====')
varlist={}
lfail = 0;

#print(repr(data['ECEREFLECT_equil']))

varlist['ECEREFLECT_equil']=np.array([ 149.62942249,  416.42189189, 1218.30181253, 2379.32028834,
       3043.94881112, 2604.26633864, 2115.48316902, 1801.01946807,
       1423.31732607,  951.10885614,  969.86915244, 1659.29686333,
       2277.7807671 , 2858.81640496, 3350.37209842, 4693.89114358,
       4974.94892212, 5206.42380758, 5419.70800538, 5378.70954407,
       5127.51852269, 4841.43230567, 4468.05249166, 4090.85368817,
       3747.40805432, 3321.13850704, 2900.84645625, 2488.74347823,
       1953.71379908, 1384.1969784 ,  747.37057215])
for temp in varlist:
    act = varlist[temp]
    cal = data[temp]
    cal = np.where(act==0,0,cal)
    div = np.where(act==0,1,act)
    perct = 100*sum(abs(act-cal)/div)
    print('  '+temp+': '+str(cal[0])+'   '+str(act[0])+'   '+str(int(perct))+'%')
    if perct > failtol:
        lfail = 1
print('=================')

if (lfail):
    print('  STATUS: FAIL!!!!!')
else:
    print('  STATUS: PASS')

sys.exit(0)




