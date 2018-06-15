#!/usr/bin/env python3
import numpy as np                    #For Arrays
from math import sqrt

lfail = 0
failtol = 1.0
filename='diagno_bench.bigtok'
data = np.loadtxt(filename)

print('==== B-Field ====')
for i in range(0, 75):
    cal = sqrt(data[i][6]*data[i][6] + data[i][7]*data[i][7] + data[i][8]*data[i][8])
    act = sqrt(data[i][12]*data[i][12] + data[i][13]*data[i][13] + data[i][14]*data[i][14])
    perct = 100*(abs(act-cal)/abs(act))
    print('   '+str(cal)+'   '+str(act)+'   '+str(int(perct))+'%')
    if perct > failtol:
        lfail = 1
print('=================')

if (lfail):
    print('  STATUS: FAIL!!!!!')
else:
    print('  STATUS: PASS')

quit()




