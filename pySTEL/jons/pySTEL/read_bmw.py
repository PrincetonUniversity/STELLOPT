# -*- coding: utf-8 -*-

from netCDF4 import Dataset
import numpy as np
import os

folder = r"Z:\03_Masterarbeit\data\tmp"
#filename = "bmw_w7x.14000.0_14000.0_13160.0_12950.0_12390.0_-9660.0_-9650.0.nc"
filename="bmw_w7x.14000.0_currents.nc"

rootgrp = Dataset(os.path.join(folder, filename), 'r')

# dimensions
r = rootgrp.dimensions['r'].size
phi = rootgrp.dimensions['phi'].size
z = rootgrp.dimensions['z'].size

# variables
version = rootgrp['/series'][0]
nfp = rootgrp['/nfp'][0]
rmin = rootgrp['/rmin'][0]
rmax = rootgrp['/rmax'][0]
zmin = rootgrp['/zmin'][0]
zmax = rootgrp['/zmax'][0]

# current density on grid
# [phi][z][r]
j_x = rootgrp['/j_x'][:]
j_y = rootgrp['/j_y'][:]
j_z = rootgrp['/j_z'][:]

## vector potential on grid
## [phi][z][r]
#ar_grid = rootgrp['/ar_grid'][:]
#ap_grid = rootgrp['/ap_grid'][:]
#az_grid = rootgrp['/az_grid'][:]
#
## magnetic flux density on grid
## [phi][z][r]
#br_grid = rootgrp['/br_grid'][:]
#bp_grid = rootgrp['/bp_grid'][:]
#bz_grid = rootgrp['/bz_grid'][:]

# toroidal index at which to plot
toroidalPhiIdx = 18


# VMEC webservice for LCFS shape
from osa import Client
vmec = Client("http://esb.ipp-hgw.mpg.de:8280/services/vmec_v7?wsdl")
lcfs = vmec.service.getFluxSurfaces('w7x.14000.0_14000.0_13160.0_12950.0_12390.0_-9660.0_-9650.0', toroidalPhiIdx/(phi-1)*2*np.pi, [1.0], 100)

lcfs_r = lcfs[0].x1
lcfs_z = lcfs[0].x3


import matplotlib.pyplot as plt


#plt.figure(1)
#plt.pcolormesh(np.linspace(rmin, rmax, r), np.linspace(zmin, zmax, z), ar_grid[toroidalPhiIdx])
#plt.plot(lcfs_r, lcfs_z, 'k')
#plt.colorbar()
#plt.title("A_r")
#plt.xlabel("r / m")
#plt.ylabel("z / m")
#
#plt.figure(2)
#plt.pcolormesh(np.linspace(rmin, rmax, r), np.linspace(zmin, zmax, z), ap_grid[toroidalPhiIdx])
#plt.plot(lcfs_r, lcfs_z, 'k')
#plt.colorbar()
#plt.title("A_phi")
#plt.xlabel("r / m")
#plt.ylabel("z / m")
#
#plt.figure(3)
#plt.pcolormesh(np.linspace(rmin, rmax, r), np.linspace(zmin, zmax, z), az_grid[toroidalPhiIdx])
#plt.plot(lcfs_r, lcfs_z, 'k')
#plt.colorbar()
#plt.title("A_z")
#plt.xlabel("r / m")
#plt.ylabel("z / m")
#
#plt.figure(4)
#plt.pcolormesh(np.linspace(rmin, rmax, r), np.linspace(zmin, zmax, z), br_grid[toroidalPhiIdx])
#plt.plot(lcfs_r, lcfs_z, 'k')
#plt.colorbar()
#plt.title("B_r")
#plt.xlabel("r / m")
#plt.ylabel("z / m")
#
#plt.figure(5)
#plt.pcolormesh(np.linspace(rmin, rmax, r), np.linspace(zmin, zmax, z), bp_grid[toroidalPhiIdx])
#plt.plot(lcfs_r, lcfs_z, 'k')
#plt.colorbar()
#plt.title("B_phi")
#plt.xlabel("r / m")
#plt.ylabel("z / m")
#
#plt.figure(6)
#plt.pcolormesh(np.linspace(rmin, rmax, r), np.linspace(zmin, zmax, z), bz_grid[toroidalPhiIdx])
#plt.plot(lcfs_r, lcfs_z, 'k')
#plt.colorbar()
#plt.title("B_z")
#plt.xlabel("r / m")
#plt.ylabel("z / m")

plt.figure(1)
plt.pcolormesh(np.linspace(rmin, rmax, r), np.linspace(zmin, zmax, z), j_x[toroidalPhiIdx])
plt.plot(lcfs_r, lcfs_z, 'k')
plt.colorbar()
plt.title("j_x")
plt.xlabel("r / m")
plt.ylabel("z / m")

plt.figure(2)
plt.pcolormesh(np.linspace(rmin, rmax, r), np.linspace(zmin, zmax, z), j_y[toroidalPhiIdx])
plt.plot(lcfs_r, lcfs_z, 'k')
plt.colorbar()
plt.title("j_y")
plt.xlabel("r / m")
plt.ylabel("z / m")

plt.figure(3)
plt.pcolormesh(np.linspace(rmin, rmax, r), np.linspace(zmin, zmax, z), j_z[toroidalPhiIdx])
plt.plot(lcfs_r, lcfs_z, 'k')
plt.colorbar()
plt.title("j_z")
plt.xlabel("r / m")
plt.ylabel("z / m")
#
#
## radial cut for comparison with Minerva
## first, find index where z value is closest to search_z
#search_z = 0.0
#z_idx = np.argmin(np.abs(np.linspace(zmin, zmax, z)-search_z))
#
#plt.figure(7)
#plt.plot(np.linspace(rmin, rmax, r), ar_grid[toroidalPhiIdx][z_idx])
#plt.xlabel("r / m")
#plt.ylabel("A / Tm")
#plt.title("A_r")
#
#plt.figure(8)
#plt.plot(np.linspace(rmin, rmax, r), ap_grid[toroidalPhiIdx][z_idx])
#plt.xlabel("r / m")
#plt.ylabel("A / Tm")
#plt.title("A_phi")
#
#
#plt.figure(9)
#plt.plot(np.linspace(rmin, rmax, r), az_grid[toroidalPhiIdx][z_idx])
#plt.xlabel("r / m")
#plt.ylabel("A / Tm")
#plt.title("A_z")
#
#plt.figure(10)
#plt.plot(np.linspace(rmin, rmax, r), br_grid[toroidalPhiIdx][z_idx])
#plt.xlabel("r / m")
#plt.ylabel("B / T")
#plt.title("B_r")
#
#plt.figure(11)
#plt.plot(np.linspace(rmin, rmax, r), bp_grid[toroidalPhiIdx][z_idx])
#plt.xlabel("r / m")
#plt.ylabel("B / T")
#plt.title("B_phi")
#
#plt.figure(12)
#plt.plot(np.linspace(rmin, rmax, r), bz_grid[toroidalPhiIdx][z_idx])
#plt.xlabel("r / m")
#plt.ylabel("B / T")
#plt.title("B_z")
