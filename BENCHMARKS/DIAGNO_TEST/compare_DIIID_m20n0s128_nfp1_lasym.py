#!/usr/bin/env python3
import numpy as np                    #For Arrays
from math import sqrt

lfail = 0
failtol = 3.0
fileext='DIIID_m20n0s128_nfp1_lasym'
filename='diagno_bench.'+fileext
data = np.loadtxt(filename)

print('===== B-Field ======')
for i in range(0, 75):
    cal = sqrt(data[i][6]*data[i][6] + data[i][7]*data[i][7] + data[i][8]*data[i][8])
    act = sqrt(data[i][12]*data[i][12] + data[i][13]*data[i][13] + data[i][14]*data[i][14])
    perct = 100*(abs(act-cal)/abs(act))
    print('   '+str(cal)+'   '+str(act)+'   '+str(int(perct))+'%')
    if perct > failtol:
        lfail = 1
print('==== Flux Loops ====')
filename='diagno_flux.'+fileext+'_j'
filename2='diagno_flux.'+fileext+'_b'
file_handle = open(filename,'r')
file_handle2 = open(filename2,'r')
line1 = file_handle.readline()
line2 = file_handle2.readline()
nlines = int(line1)
for i in range(0,nlines):
	line1 = file_handle.readline()
	line2 = file_handle2.readline()
	cal   = float(line1)
	act   = float(line2)
	if abs(act) < 1E-3:
		continue
	if abs(cal) < 1E-3:
		continue
	perct = 100*(abs(act-cal)/abs(act))
	print('   '+str(cal)+'   '+str(act)+'   '+str(int(perct))+'%')
	if perct > failtol:
		lfail = 1
file_handle.close()
file_handle2.close()
#with open(filename) as myfile:
#    head = myfile.next()
#print(head)
print('=================')

if (lfail):
    print('  STATUS: FAIL!!!!!')
else:
    print('  STATUS: PASS')

quit()




