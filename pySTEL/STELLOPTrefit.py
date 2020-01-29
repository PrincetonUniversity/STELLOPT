#!/usr/bin/env python3
import sys, os, glob
sys.path.insert(0, '../../../pySTEL/') 
import ctypes as ct
import numpy as np                    #For Arrays
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from math import pi
from libstell.stellopt import read_stellopt


def fit_poly10(x,am0,am1,am2,am3,am4,am5,am6,am7,am8,am9,am10):
	f2=np.zeros(x.shape)
	for i,xx in enumerate(x):
		f2[i] = (f2[i]+am10)*xx
		f2[i] = (f2[i]+am9)*xx
		f2[i] = (f2[i]+am8)*xx
		f2[i] = (f2[i]+am7)*xx
		f2[i] = (f2[i]+am6)*xx
		f2[i] = (f2[i]+am5)*xx
		f2[i] = (f2[i]+am4)*xx
		f2[i] = (f2[i]+am3)*xx
		f2[i] = (f2[i]+am2)*xx
		f2[i] = (f2[i]+am1)*xx
		f2[i] = f2[i]+am0
		if (xx>1):
			f[i]=0
	return f


def fit_poly5(x,am0,am1,am2,am3,am4,am5):
	f2=np.zeros(x.shape)
	for i,xx in enumerate(x):
		f2[i] = (f2[i]+am5)*xx
		f2[i] = (f2[i]+am4)*xx
		f2[i] = (f2[i]+am3)*xx
		f2[i] = (f2[i]+am2)*xx
		f2[i] = (f2[i]+am1)*xx
		f2[i] = f2[i]+am0
		if (xx>1):
			f2[i]=0
	return f2

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
print('-----  Fitting Ne -----')
s=data['NE_s'][ndex-1,:]
f=data['NE_target'][ndex-1,:]
s_fit=np.array([0,0.04,0.16,0.36,0.64,1.0])
f_fit = curve_fit(fit_poly5, s, f)
f_err = np.sqrt(np.diag(f_fit[1]))
f_fit = f_fit[0]
plt.plot(s,f,'o')
plt.plot(s_fit,fit_poly5(s_fit,f_fit[0],f_fit[1],f_fit[2],f_fit[3],f_fit[4],f_fit[5]),'+')
plt.show()
