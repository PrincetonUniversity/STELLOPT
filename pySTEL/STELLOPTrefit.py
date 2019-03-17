#!/usr/bin/env python3
import sys, os, globsys.path.insert(0, '../../../pySTEL/')
import numpy as np                    #For Arrays
from math import pi
from libstell.stellopt import read_stellopt

try:
	qtCreatorPath=os.environ["STELLOPT_PATH"]
except KeyError:
	print("Please set environment variable STELLOPT_PATH")
	sys.exit(1)

npoly=5
print('----- STELLOPT REFIT -----')
files = sorted(glob.glob('./stellopt.*'))
print('  Opening File:'+files[0])
data=read_stellopt(files[0])
ndex=data['ITER'].size
print('  Iterations Detected: '+str(ndex))
print('-----  Fitting ne -----')
s=data['NE_s'][ndex-1,:]
f=data['NE_target'][ndex-1,:]
s_fit=[0,0.04,0.16,0.36,0.64,1.0]
f_fit=np.polyfit(s,f,npoly)
print(f_fit)
