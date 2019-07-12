# -*- coding: utf-8 -*-

# mir_tools.f90 in LIBSTELL contains some useful methods for converting
# coordinates between real space and flux coordinates and
# calculate useful quantitiesfrom the flux surface output of VMEC.


import ctypes as ct
import numpy.ctypeslib as npct
import numpy as np

#file = "/mnt/jons/03_Masterarbeit/data/KJM001/wout_w7x_ref_164.nc"

# http://svvmec1.ipp-hgw.mpg.de:8080/vmecrest/v1/run/w7x.14000.0_14000.0_13160.0_12950.0_12390.0_-9660.0_-9650.0
file = "/mnt/jons/03_Masterarbeit/data/tmp/wout_tmp2.nc"

# Load Libraries
#libstell = ct.cdll.LoadLibrary("/u/slazerso/src/STELLOPT_GCC/LIBSTEL/Release/libstell.so")
#libstell = ct.cdll.LoadLibrary("/home/IPP-HGW/jons/bin/libstell.so")
libstell = ct.cdll.LoadLibrary("/home/jonathan/bin/libstell.so")
#
#initialize_from_VMEC = getattr(libstell,'__mir_tools_MOD_initialize_from_vmec')
#initialize_from_VMEC.argparse=[ct.c_int, ct.c_char_p, ct.c_int, ct.c_int, ct.c_int]
#initialize_from_VMEC.restype=None
error_status = ct.c_int(0)
#nu = ct.c_int(51)
#nv = ct.c_int(101)
#initialize_from_VMEC(ct.byref(error_status), file.encode('UTF-8'), ct.byref(nu), ct.byref(nv), len(file))
#print("error_status after init: " + np.str(error_status))
#
#
#b_car_from_cyl = getattr(libstell, '__mir_tools_MOD_b_car_from_cyl')
#b_car_from_cyl.argparse=[ct.ARRAY(ct.c_double, 3), ct.ARRAY(ct.c_double, 3), ct.c_int]
#b_car_from_cyl.restype=None
#
##point_cyl=np.array([5.8, 0.0, 0.0])
##b_car=np.array([5.8, 0.0, 0.0])
#
#point_cyl=(ct.c_double*3)()
#point_cyl[0]=5.8
#b_car=(ct.c_double*3)()
#b_car_from_cyl(ct.byref(point_cyl), ct.byref(b_car), ct.byref(error_status))
#
#print("error_status after b_car_from_cyl: " + np.str(error_status))
#
#print(b_car[:])
#
#minerva_comparison=[-1.167710184123799E-16, 2.565499687444305, 0.8917631776620901]
#diff=np.linalg.norm(np.multiply(-1.0, minerva_comparison)-b_car[:])
#print("b_car_from_cyl diff to Minerva: "+np.str(diff))
#
#


#
#
#####################################################
## the following uses the methods from vmec_utils.f #
#####################################################
#
## cylindrical components of magnetic flux density B
#
#getbcyl_wout = getattr(libstell, '__vmec_utils_MOD_getbcyl_wout')
#getbcyl_wout.argparse=[ct.c_double, ct.c_double, ct.c_double,
#                       ct.c_double, ct.c_double, ct.c_double,
#                       ct.c_double, ct.c_double, ct.c_int]
#getbcyl_wout.restype=None
#
#R1=(ct.c_double)(5.8)
#Phi=(ct.c_double)(0.0)
#Z1=(ct.c_double)(0.0)
#
#Br=(ct.c_double)()
#Bphi=(ct.c_double)()
#Bz=(ct.c_double)()
#
#sflx=(ct.c_double)()
#uflx=(ct.c_double)()
#
#getbcyl_wout(ct.byref(R1), ct.byref(Phi), ct.byref(Z1),
#             ct.byref(Br), ct.byref(Bphi), ct.byref(Bz),
#             ct.byref(sflx), ct.byref(uflx), ct.byref(error_status))
#
#print("error_status after getbcyl_wout: " + np.str(error_status))
#
#print([Br, Bphi, Bz])
#
#minerva_comparison=[-1.167710184123799E-16, 2.565499687444305, 0.8917631776620901]
#diff=np.linalg.norm(np.multiply(-1.0, minerva_comparison)-[Br.value, Bphi.value, Bz.value])
#print("getbcyl_wout diff to Minerva: "+np.str(diff))

# cylindrical components of magnetic vector potential A

read_wout = getattr(libstell,'__read_wout_mod_MOD_readw_and_open')
read_wout.argparse=[ct.c_char_p, ct.c_int, ct.c_int, ct.c_int]
read_wout.restype=None
ierr = ct.c_int(0)
iopen = ct.c_int(0)
read_wout(file.encode('UTF-8'), ct.byref(ierr), iopen, len(file))


getacyl_wout = getattr(libstell, '__vmec_utils_MOD_getacyl_wout')
getacyl_wout.argparse=[ct.c_double, ct.c_double, ct.c_double,
                       ct.c_double, ct.c_double, ct.c_double,
                       ct.c_double, ct.c_double, ct.c_int]
getacyl_wout.restype=None

def getVectorPotential(pointsXYZ):

    allX = pointsXYZ[0][:]
    allY = pointsXYZ[1][:]
    allZ = pointsXYZ[2][:]
    
    A = np.zeros([3, len(allX)])
    
    for i in range(len(allX)):
        
        toroidalR = np.sqrt(allX[i]*allX[i]+allY[i]*allY[i])
        toroidalPhi = np.arctan2(allY[i], allX[i])
    
        R1=(ct.c_double)(toroidalR)
        Phi=(ct.c_double)(toroidalPhi)
        Z1=(ct.c_double)(allZ[i])
        
        Ar=(ct.c_double)()
        Aphi=(ct.c_double)()
        Az=(ct.c_double)()
        
        sflx=(ct.c_double)()
        uflx=(ct.c_double)()
        
        getacyl_wout(ct.byref(R1), ct.byref(Phi), ct.byref(Z1),
                     ct.byref(Ar), ct.byref(Aphi), ct.byref(Az),
                     ct.byref(sflx), ct.byref(uflx), ct.byref(error_status))
        
#        print("error_status after getacyl_wout: " + np.str(error_status))
        
#        print([Ar.value, Aphi.value, Az.value])

        A[0][i] = Ar.value * np.cos(toroidalPhi) - Aphi.value * np.sin(toroidalPhi);
        A[1][i] = Ar.value * np.sin(toroidalPhi) + Aphi.value * np.cos(toroidalPhi);
        A[2][i] = Az.value
    
    return A

    
testX = np.linspace(4, 8, np.int((8-4)/0.01))
pointsXYZ=np.zeros([3,len(testX)])
pointsXYZ[0][:] = testX

myA = getVectorPotential(pointsXYZ)

import matplotlib.pyplot as plt

plt.figure(1)
plt.plot(testX, myA[0][:])

plt.figure(2)
plt.plot(testX, myA[1][:])

plt.figure(3)
plt.plot(testX, myA[2][:])


#minerva_comparison=[-1.167710184123799E-16, 2.565499687444305, 0.8917631776620901]
#diff=np.linalg.norm(np.multiply(-1.0, minerva_comparison)-[Br.value, Bphi.value, Bz.value])
#print("getacyl_wout diff to Minerva: "+np.str(diff))




