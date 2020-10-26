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
failtol = 1.0
filename='beams3d_ORBITS_fusion.h5'
data=read_beams3d(filename)
if not data:
    print('ERROR Opening File: '+filename)
    sys.exit(0)

# Calc values
rho = np.sqrt(data['S_lines'])
print(data.keys())
weight = data['Weight']
beam   = data['Beam']

wtemp = np.where(beam==1,weight,0)
data['birth_rate_He4'], temp = np.histogram(rho[:,0],weights=wtemp,bins=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
wtemp = np.where(beam==2,weight,0)
data['birth_rate_T'], temp = np.histogram(rho[:,0],weights=wtemp,bins=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
wtemp = np.where(beam==3,weight,0)
data['birth_rate_p'], temp = np.histogram(rho[:,0],weights=wtemp,bins=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
wtemp = np.where(beam==4,weight,0)
data['birth_rate_He3'], temp = np.histogram(rho[:,0],weights=wtemp,bins=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])


print('BEAMS3D VERSION: ' + str(data['VERSION']))
print('==== Vectors ====')
varlist={}
varlist['birth_rate_He4']=np.array([1.690304e+15,4.291248e+15,6.679591e+15,7.726770e+15,7.798617e+15,
 6.931776e+15,3.878610e+15,1.212546e+15,8.749673e+13,2.337517e+11])
varlist['birth_rate_T']=np.array([8.607383e+12,2.261953e+13,3.653183e+13,3.648088e+13,4.231431e+13,
 4.001636e+13,2.462698e+13,7.971497e+12,7.621906e+11,2.424376e+09])
varlist['birth_rate_p']=np.array([8.607383e+12,2.197273e+13,3.606611e+13,3.850850e+13,4.292426e+13,
 3.755047e+13,2.536364e+13,8.196269e+12,7.417118e+11,2.320129e+09])
varlist['birth_rate_He3']=np.array([8.964892e+12,2.218304e+13,3.585884e+13,4.311138e+13,4.208365e+13,
 4.054499e+13,2.459525e+13,8.636335e+12,6.508350e+11,2.605140e+09])
for temp in varlist:
    act = varlist[temp]
    cal = data[temp]
    #print(np.array2string(cal,precision=6, separator=','))
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




