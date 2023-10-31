#!/usr/bin/env python3
import sys, os
sys.path.insert(0, '../../pySTEL/')
import numpy as np                    #For Arrays
from math import pi
from libstell.fidasim import read_fidasim

try:
	qtCreatorPath=os.environ["STELLOPT_PATH"]
except KeyError:
	print("Please set environment variable STELLOPT_PATH")
	sys.exit(1)

lfail = 0
failtol = 15.0
filename='fidasim_ORBITS_FIDASIM_distribution.h5'
data=read_fidasim(filename)
if not data:
    print('ERROR Opening File: '+filename)
    sys.exit(0)

# Calc values
data['n']=np.array([np.squeeze(data['denf'])*1e6*10*2*3.141593*3*2]) #cm^-3, R*2*pi*(rmax-rmin)*(zmax-zmin)
data['pitch']=np.squeeze(data['f'])*1e6*10*2*3.141593*3*2*2*data['energy']#cm^-3, R*2*pi*(rmax-rmin)*(zmax-zmin)

print('==== Vectors ====')
varlist={}
varlist['n']=np.array([0.0024])
varlist['pitch']=np.array([0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0002,0.0002,0.0002,0.0002,0.0002,0.0003,0.0003,0.0003,0.0003,0.0003,0.0003,0.0003,0.0004,0.0006,0.0007,0.0009,0.0011,0.0012,0.0014,0.0017,0.0018,0.0020,0.0026,0.0030,0.0040,0.0051,0.0082,0.0167])
for temp in varlist:
    act = varlist[temp]
    cal = data[temp]
    cal = np.where(act==0,0,cal)  
    perct = 100*np.sqrt(sum((act-cal)**2)/len(act))/np.mean(act[np.nonzero(act)]) #normalized root mean squared deviation
    #perct = max(abs(act-cal))
    print(f'  {temp}: {cal[-1]:.4f}   {act[-1]:.4f}   {perct:.2f}%')
    #print(cal)
    if perct > failtol:
        lfail = 1
print('=================')

if (lfail):
    print('  STATUS: FAIL!!!!!')
else:
    print('  STATUS: PASS')

sys.exit(0)




