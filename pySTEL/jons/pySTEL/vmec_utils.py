# -*- coding: utf-8 -*-

import ctypes as ct
import numpy as np

libstell = ct.cdll.LoadLibrary("/home/IPP-HGW/jons/bin/libstell.so")

ierr = ct.c_int(0)

# import method to read wout file
read_wout = getattr(libstell,'__read_wout_mod_MOD_readw_and_open')
read_wout.argparse=[ct.c_char_p, ct.c_int, ct.c_int, ct.c_int]
read_wout.restype=None

# import method to calculate magnetic flux density B
getbcyl_wout = getattr(libstell, '__vmec_utils_MOD_getbcyl_wout')
getbcyl_wout.argparse=[ct.c_double, ct.c_double, ct.c_double,
                       ct.c_double, ct.c_double, ct.c_double,
                       ct.c_double, ct.c_double, ct.c_int]
getbcyl_wout.restype=None

# import method to calculate magnetic vector potential A
getacyl_wout = getattr(libstell, '__vmec_utils_MOD_getacyl_wout')
getacyl_wout.argparse=[ct.c_double, ct.c_double, ct.c_double,
                       ct.c_double, ct.c_double, ct.c_double,
                       ct.c_double, ct.c_double, ct.c_int]
getacyl_wout.restype=None

# import method to calculate current density J
getjcyl_wout = getattr(libstell, '__vmec_utils_MOD_getjcyl_wout')
getjcyl_wout.argparse=[ct.c_double, ct.c_double, ct.c_double,
                       ct.c_double, ct.c_double, ct.c_double,
                       ct.c_double, ct.c_double, ct.c_int]
getjcyl_wout.restype=None

# init libstell by reading a wout file
def init(wout_filename):
    iopen = ct.c_int(0)
    read_wout(wout_filename.encode('UTF-8'), ct.byref(ierr), iopen, len(wout_filename))
    return ierr

# calculate the magnetic flux density B
# pointsXYZ[xyz][numPoints] => 3xN matrix of cartesian positions where B should be calculated
# returns B[xyz][numPoints] => 3xN matrix of cartesian components of magnetic flux density
def getMagneticField(pointsXYZ):
    
    allX = pointsXYZ[0][:]
    allY = pointsXYZ[1][:]
    allZ = pointsXYZ[2][:]
    
    B = np.zeros([3, len(allX)])
    
    for i in range(len(allX)):
        
        toroidalR = np.sqrt(allX[i]*allX[i]+allY[i]*allY[i])
        toroidalPhi = np.arctan2(allY[i], allX[i])
    
        R1=(ct.c_double)(toroidalR)
        Phi=(ct.c_double)(toroidalPhi)
        Z1=(ct.c_double)(allZ[i])
        
        Br=(ct.c_double)()
        Bphi=(ct.c_double)()
        Bz=(ct.c_double)()
        
        sflx=(ct.c_double)()
        uflx=(ct.c_double)()
        
        getbcyl_wout(ct.byref(R1), ct.byref(Phi), ct.byref(Z1),
                     ct.byref(Br), ct.byref(Bphi), ct.byref(Bz),
                     ct.byref(sflx), ct.byref(uflx), ct.byref(ierr))
        
#        print("error_status after getbcyl_wout: " + np.str(ierr))
#        print([Br.value, Bphi.value, Bz.value])

        B[0][i] = Br.value * np.cos(toroidalPhi) - Bphi.value * np.sin(toroidalPhi);
        B[1][i] = Br.value * np.sin(toroidalPhi) + Bphi.value * np.cos(toroidalPhi);
        B[2][i] = Bz.value
    
    return B

# calculate the magnetic vector potential A
# pointsXYZ[xyz][numPoints] => 3xN matrix of cartesian positions where A should be calculated
# returns A[xyz][numPoints] => 3xN matrix of cartesian components of magnetic vector potential
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
                     ct.byref(sflx), ct.byref(uflx), ct.byref(ierr))
        
#        print("error_status after getacyl_wout: " + np.str(ierr))
#        print([Ar.value, Aphi.value, Az.value])

        A[0][i] = Ar.value * np.cos(toroidalPhi) - Aphi.value * np.sin(toroidalPhi);
        A[1][i] = Ar.value * np.sin(toroidalPhi) + Aphi.value * np.cos(toroidalPhi);
        A[2][i] = Az.value
    
    return A

# calculate the current density J
# pointsXYZ[xyz][numPoints] => 3xN matrix of cartesian positions where J should be calculated
# returns J[xyz][numPoints] => 3xN matrix of cartesian components of current density
def getCurrentDensity(pointsXYZ):

    allX = pointsXYZ[0][:]
    allY = pointsXYZ[1][:]
    allZ = pointsXYZ[2][:]
    
    J = np.zeros([3, len(allX)])
    
    for i in range(len(allX)):
        
        toroidalR = np.sqrt(allX[i]*allX[i]+allY[i]*allY[i])
        toroidalPhi = np.arctan2(allY[i], allX[i])
    
        R1=(ct.c_double)(toroidalR)
        Phi=(ct.c_double)(toroidalPhi)
        Z1=(ct.c_double)(allZ[i])
        
        Jr=(ct.c_double)()
        Jphi=(ct.c_double)()
        Jz=(ct.c_double)()
        
        sflx=(ct.c_double)()
        uflx=(ct.c_double)()
        
        getjcyl_wout(ct.byref(R1), ct.byref(Phi), ct.byref(Z1),
                     ct.byref(Jr), ct.byref(Jphi), ct.byref(Jz),
                     ct.byref(sflx), ct.byref(uflx), ct.byref(ierr))
        
#        print("error_status after getjcyl_wout: " + np.str(ierr))
#        print([Jr.value, Jphi.value, Jz.value])

        J[0][i] = Jr.value * np.cos(toroidalPhi) - Jphi.value * np.sin(toroidalPhi);
        J[1][i] = Jr.value * np.sin(toroidalPhi) + Jphi.value * np.cos(toroidalPhi);
        J[2][i] = Jz.value

        # internally: sflux == jsupu1, uflx = jsupv1
#        J[0][i] = sflx.value
#        J[1][i] = uflx.value
    
    return J


vmec_id = "w7x.14000.0_14000.0_13160.0_12950.0_12390.0_-9660.0_-9650.0"
#wout_filename = "/mnt/jons/03_Masterarbeit/data/tmp/vmec_lbsubs.false/wout_w7x.14000.0_14000.0_13160.0_12950.0_12390.0_-9660.0_-9650.0.txt.nc"
#wout_filename = "/mnt/jons/03_Masterarbeit/data/tmp/vmec_lbsubs.false/wout.nc"
wout_filename = "/mnt/jons/03_Masterarbeit/data/tmp/vmec_lbsubs.true/wout.nc"
init(wout_filename)

gridExtentR = 0.8
gridExtentZ = 0.5
numGridCellsR = 250
numGridCellsZ = 170
toroidalPhi = 180*(np.pi/180.0) # 0: bean-shaped plane, 180: triangular-shaped plane
#
#grid = np.zeros([3,numGridCellsR*numGridCellsZ])
#
## get magnetic axis from VMEC webservice
#from osa import Client
#vmec = Client("http://esb.ipp-hgw.mpg.de:8280/services/vmec_v7?wsdl")
#axis = vmec.service.getMagneticAxis(vmec_id, toroidalPhi)
#
#idx=0
#for i in range(numGridCellsR):
#	for j in range(numGridCellsZ):
#
#		r = axis.x1[0] + gridExtentR*(-1.+(i+0.5)*2.0/(numGridCellsR));
#
#		grid[0][idx] = r*np.cos(toroidalPhi)
#		grid[1][idx] = r*np.sin(toroidalPhi)
#		grid[2][idx] = axis.x3[0] + gridExtentZ*(-1.+(j+0.5)*2.0/(numGridCellsZ))
#         
#		idx+=1
#        
#currentDensityXYZ = getCurrentDensity(grid)
#
#
#currentDensityOnGrid = np.zeros([numGridCellsZ, numGridCellsR])
#

#%%

import matplotlib.pyplot as plt

#
#idx=0
#for i in range(numGridCellsR):
#	for j in range(numGridCellsZ):
#
#		modJ = np.sqrt( np.power(currentDensityXYZ[0][idx],2)
#				      +np.power(currentDensityXYZ[1][idx],2)
#				      +np.power(currentDensityXYZ[2][idx],2));
#
#		currentDensityOnGrid[j][i] = modJ;
#
#		idx+=1
#
#
#plt.figure()
#plt.imshow(currentDensityOnGrid, interpolation='nearest', cmap=plt.get_cmap('jet'),
#           extent=[-1.0*gridExtentR+axis.x1[0], gridExtentR+axis.x1[0],
#		          -1.0*gridExtentZ+axis.x3[0], gridExtentZ+axis.x3[0]],
#                   vmin=0, vmax=6e5)
#plt.colorbar()
#
#plt.xlim(-1.0*gridExtentR+axis.x1[0], gridExtentR+axis.x1[0])
#plt.ylim(-1.0*gridExtentZ+axis.x3[0], gridExtentZ+axis.x3[0])
#plt.xlabel("r / m")
#plt.ylabel("z / m")
#plt.title("|J| in A m^-2")
#plt.show()
#
##%%
#idx=0
#for i in range(numGridCellsR):
#	for j in range(numGridCellsZ):
#
#		currentDensityOnGrid[j][i] = currentDensityXYZ[0][idx];
#
#		idx+=1
#
#plt.figure()
#plt.imshow(currentDensityOnGrid, interpolation='nearest', cmap=plt.get_cmap('jet'),
#           extent=[-1.0*gridExtentR+axis.x1[0], gridExtentR+axis.x1[0],
#		          -1.0*gridExtentZ+axis.x3[0], gridExtentZ+axis.x3[0]],vmin=-5e5, vmax=5e5)
#plt.colorbar()
#
#plt.xlim(-1.0*gridExtentR+axis.x1[0], gridExtentR+axis.x1[0])
#plt.ylim(-1.0*gridExtentZ+axis.x3[0], gridExtentZ+axis.x3[0])
#plt.xlabel("r / m")
#plt.ylabel("z / m")
#plt.title("J_x in A m^-2")
#plt.show()
#
##np.savetxt("/mnt/jons/03_Masterarbeit/data/tmp/jsupu_tmp2.dat", currentDensityOnGrid)
#
##%%
#idx=0
#for i in range(numGridCellsR):
#	for j in range(numGridCellsZ):
#
#		currentDensityOnGrid[j][i] = currentDensityXYZ[1][idx];
#
#		idx+=1
#
#plt.figure()
#plt.imshow(currentDensityOnGrid, interpolation='nearest', cmap=plt.get_cmap('jet'),
#           extent=[-1.0*gridExtentR+axis.x1[0], gridExtentR+axis.x1[0],
#		          -1.0*gridExtentZ+axis.x3[0], gridExtentZ+axis.x3[0]],
#                   vmin=-5e5, vmax=5e5)
#plt.colorbar()
#
#plt.xlim(-1.0*gridExtentR+axis.x1[0], gridExtentR+axis.x1[0])
#plt.ylim(-1.0*gridExtentZ+axis.x3[0], gridExtentZ+axis.x3[0])
#plt.xlabel("r / m")
#plt.ylabel("z / m")
#plt.title("J_y in A m^-2")
#plt.show()
#
##np.savetxt("/mnt/jons/03_Masterarbeit/data/tmp/jsupv_tmp2.dat", currentDensityOnGrid)
#
##%%
#idx=0
#for i in range(numGridCellsR):
#	for j in range(numGridCellsZ):
#
#		currentDensityOnGrid[j][i] = currentDensityXYZ[2][idx];
#
#		idx+=1
#
#plt.figure()
#plt.imshow(currentDensityOnGrid, interpolation='nearest', cmap=plt.get_cmap('jet'),
#           extent=[-1.0*gridExtentR+axis.x1[0], gridExtentR+axis.x1[0],
#		          -1.0*gridExtentZ+axis.x3[0], gridExtentZ+axis.x3[0]],
#                   vmin=-5e5, vmax=5e5)
#plt.colorbar()
#
#plt.xlim(-1.0*gridExtentR+axis.x1[0], gridExtentR+axis.x1[0])
#plt.ylim(-1.0*gridExtentZ+axis.x3[0], gridExtentZ+axis.x3[0])
#plt.xlabel("r / m")
#plt.ylabel("z / m")
#plt.title("J_z in A m^-2")
#plt.show()


#%%
testXYZ=np.loadtxt("/mnt/jons/03_Masterarbeit/data/tmp/testXYZ.dat")
A_xyz = getVectorPotential(testXYZ)
B_xyz = getMagneticField(testXYZ)