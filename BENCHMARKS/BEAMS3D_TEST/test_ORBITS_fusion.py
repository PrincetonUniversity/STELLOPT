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
filename='beams3d_ORBITS_fusion.h5'
data=read_beams3d(filename)
if not data:
    print('ERROR Opening File: '+filename)
    sys.exit(0)

# Calc values
rho = np.sqrt(data['S_lines'])
#print(data.keys())
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
varlist['birth_rate_He4']=np.array([5.531974e+17,1.548293e+18,2.260494e+18,2.671512e+18,2.731105e+18,
 2.334851e+18,1.394123e+18,4.132231e+17,3.011013e+16,7.444151e+13])
varlist['birth_rate_T']=np.array([1.405567e+15,3.948818e+15,5.843987e+15,7.007070e+15,7.375928e+15,
 6.614761e+15,4.276418e+15,1.427369e+15,1.248056e+14,4.142369e+11])
varlist['birth_rate_p']=np.array([1.405789e+15,3.945669e+15,5.851363e+15,6.997016e+15,7.383014e+15,
 6.607894e+15,4.281632e+15,1.427379e+15,1.249589e+14,4.210743e+11])
varlist['birth_rate_He3']=np.array([1.469183e+15,4.105819e+15,6.047775e+15,7.267702e+15,7.602216e+15,
 6.774583e+15,4.349915e+15,1.442406e+15,1.245455e+14,4.117774e+11])
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




