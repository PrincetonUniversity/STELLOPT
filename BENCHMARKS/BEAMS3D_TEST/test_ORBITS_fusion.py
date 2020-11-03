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
varlist['birth_rate_He4']=np.array([1.390846e+15,5.372319e+15,5.512664e+15,6.764489e+15,9.216209e+15,
 6.320371e+15,3.591960e+15,9.335272e+14,7.593575e+13,2.104742e+11])
varlist['birth_rate_T']=np.array([5.696724e+12,2.693095e+13,3.148294e+13,3.392625e+13,4.941253e+13,
 3.611799e+13,2.249715e+13,6.620694e+12,6.573694e+11,2.118624e+09])
varlist['birth_rate_p']=np.array([6.387716e+12,2.688948e+13,3.083769e+13,3.574146e+13,4.799576e+13,
 3.642672e+13,2.201255e+13,6.386955e+12,6.645299e+11,1.849002e+09])
varlist['birth_rate_He3']=np.array([5.932616e+12,2.869994e+13,3.143493e+13,3.566952e+13,5.134181e+13,
 3.734044e+13,2.228550e+13,6.568891e+12,6.217918e+11,2.142508e+09])
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




