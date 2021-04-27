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

#print(data.keys())
#print(repr(data['ECEREFLECT_tradx']))
#print(repr(data['ECEREFLECT_trado']))

varlist['ECEREFLECT_tradx']=np.array([ 100.14191544,  285.21383012,  884.16879877, 1920.42636854,
       2889.45412484, 2698.37612494, 1932.0825585 , 1493.59875332,
       1080.03304214,  785.06156753,  593.39493069,  517.27333818,
        763.16261791, 1782.79894957, 2227.78326621, 3712.56188131,
       4091.79497205, 4430.40275729, 4771.69904288, 4956.69614992,
       4861.64003972, 4630.95327536, 4341.08614299, 4025.11231281,
       3727.90334093, 3344.62537682, 2959.30840487, 2572.07145188,
       2057.97744265, 1498.441035  ,  924.89983331])
varlist['ECEREFLECT_trado']=np.array([1.73981845e+00, 4.74143731e+00, 1.45895739e+01, 3.59207001e+01,
       8.03194322e+01, 1.70228520e+02, 3.83893541e+02, 6.31224848e+02,
       8.06319901e+02, 1.00839363e+03, 1.26284981e+03, 1.61271120e+03,
       1.91459104e+03, 2.19753637e+03, 2.40011486e+03, 2.59166798e+03,
       2.48922903e+03, 2.34412139e+03, 2.15584462e+03, 1.95728822e+03,
       1.78325983e+03, 1.59141312e+03, 1.34353208e+03, 1.10743056e+03,
       9.00896255e+02, 6.62340739e+02, 4.55075596e+02, 2.86963791e+02,
       1.30606126e+02, 4.23401338e+01, 1.05443887e+01])
for temp in varlist:
    act = varlist[temp]
    cal = data[temp]
    cal = np.where(act==0,0,cal)
    div = np.where(act==0,1,act)
    perct = 100*sum(abs(act-cal)/div)
    print('  '+temp+': '+str(cal[14])+'   '+str(act[14])+'   '+str(int(perct))+'%')
    if perct > failtol:
        lfail = 1
print('=================')

if (lfail):
    print('  STATUS: FAIL!!!!!')
else:
    print('  STATUS: PASS')

sys.exit(0)




