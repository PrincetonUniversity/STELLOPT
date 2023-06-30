#!/usr/bin/env python3
import sys, os
sys.path.insert(0, '../../pySTEL/')
import numpy as np                    #For Arrays
from math import pi
from libstell.thrift import read_thrift
from scipy.optimize import curve_fit

if False:
    import matplotlib.pyplot as plt
    lplot = True
else:
    lplot = False


try:
	qtCreatorPath=os.environ["STELLOPT_PATH"]
except KeyError:
	print("Please set environment variable STELLOPT_PATH")
	sys.exit(1)

lfail = 0
failtol = 1.0
filename='thrift_analytical.h5'
data=read_thrift(filename)
if not data:
    print('ERROR Opening File: '+filename)
    sys.exit(0)

fit_data = data['THRIFT_I'][:,-1]
fit_time = data['THRIFT_T']
fit_func = lambda x,I0,tau,x0 : I0*(1-np.exp(-(x-x0)/tau))
bounds = ([-np.inf,1,-10],[np.inf,1000,np.inf])
p0 = np.array([fit_data.max(),10.0,0])
popt, pcov = curve_fit(fit_func, fit_time, fit_data, p0=p0, bounds = bounds)
data['I_FIT'] = popt[0]
data['TAU'] = popt[1]
if lplot:
    plt.plot(fit_time, fit_data, 'ko', label='Data')
    x2 = np.linspace(0,max(fit_time)*1.2,1000)
    plt.plot(x2,fit_func(x2,*popt), 'r-', label='Fit: I0=%5.3f, tau=%5.3f, xo=%5.3f' % tuple(popt))
    plt.xlabel('T [s]')
    plt.ylabel('Current [A]')
    plt.legend()
    plt.show()

print('THRIFT VERSION: ' + str(data['VERSION']))
print('==== Vectors ====')
varlist={}
varlist['I_FIT'] = np.array([-21100.0])
varlist['TAU'] = np.array([6.85])
for temp in varlist:
    act = varlist[temp]
    cal = data[temp]
    #print(np.array2string(cal,precision=6, separator=','))
    cal = np.where(act==0,0,cal)
    div = np.where(act==0,1,act)
    perct = 100*sum(abs(act-cal)/div)
    print('  '+temp+': '+str(cal[0])+'   '+str(act[0])+'   '+str(round(perct))+'%')
    if perct > failtol:
        lfail = 1
print('=================')

if (lfail):
    print('  STATUS: FAIL!!!!!')
else:
    print('  STATUS: PASS')

sys.exit(0)




