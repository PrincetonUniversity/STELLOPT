#!/usr/bin/env python3
import numpy as np                    #For Arrays
from math import sqrt

lfail = 0
failtol = 1.0
fileext='bigtok'
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
#with open(filename) as myfile:
#    head = myfile.next()
#print(head)
filename='diagno_flux.'+fileext+'_b'
print('=================')

if (lfail):
    print('  STATUS: FAIL!!!!!')
else:
    print('  STATUS: PASS')

quit()




