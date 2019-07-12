# -*- coding: utf-8 -*-

from netCDF4 import Dataset
import numpy as np
import os

folder = r"Z:\03_Masterarbeit\data\tmp"
filename = "bmw_w7x.14000.0_14000.0_13160.0_12950.0_12390.0_-9660.0_-9650.0.nc"

rootgrp = Dataset(os.path.join(folder, filename), 'r')

# dimensions
num_r   = rootgrp.dimensions['r'].size
num_phi = rootgrp.dimensions['phi'].size
num_z   = rootgrp.dimensions['z'].size

# variables
version = rootgrp['/series'][0]
nfp = rootgrp['/nfp'][0]
rmin = rootgrp['/rmin'][0]
rmax = rootgrp['/rmax'][0]
zmin = rootgrp['/zmin'][0]
zmax = rootgrp['/zmax'][0]

#filename = r"Z:\03_Masterarbeit\data\tmp\fort.42"
#
#if fort42 == None or filename != loaded_filename:
#    fort42 = np.loadtxt(filename).T
#    loaded_filename = filename

# from console debug output
#num_v num_u ns:          180         101          99
num_v = 180
num_u = 101
ns = 99


r = fort42[0][:].reshape([ns, num_u, num_v], order='F')
v = fort42[1][:].reshape([ns, num_u, num_v], order='F')
z = fort42[2][:].reshape([ns, num_u, num_v], order='F')

drdu = fort42[3][:].reshape([ns, num_u, num_v], order='F')
drdv = fort42[4][:].reshape([ns, num_u, num_v], order='F')

dzdu = fort42[5][:].reshape([ns, num_u, num_v], order='F')
dzdv = fort42[6][:].reshape([ns, num_u, num_v], order='F')

ju = fort42[7][:].reshape([ns, num_u, num_v], order='F')
jv = fort42[8][:].reshape([ns, num_u, num_v], order='F')

jx = fort42[9][:].reshape([ns, num_u, num_v], order='F')
jy = fort42[10][:].reshape([ns, num_u, num_v], order='F')
jz = fort42[11][:].reshape([ns, num_u, num_v], order='F')

x = np.multiply(r, np.cos(v))
y = np.multiply(r, np.sin(v))

toroidalPlotIndex = 18
phi = v[0,0,toroidalPlotIndex]/np.pi*180.0


# pcolormesh wants to have the first value in each dimension to be repeated at the end
# in case the coordinates are periodic. So we need to append them to the end of each
# dimension in the arrays to plot.
a=r[:,:,toroidalPlotIndex]
extended_r=np.vstack([a.T, a[:,0]]).T

b=z[:,:,toroidalPlotIndex]
extended_z=np.vstack([b.T, b[:,0]]).T

#%%
import matplotlib.pyplot as plt

plt.figure()
plt.pcolormesh(extended_r, extended_z, jx[:,:,toroidalPlotIndex])
cbar = plt.colorbar()
cbar.set_label("current density / A/m^2")
plt.axis("equal")
plt.xlabel("r / m")
plt.ylabel("z / m")
plt.title("j_x only plasma at phi=" + np.str(phi) + "°")
plt.tight_layout();

plt.figure()
plt.pcolormesh(extended_r, extended_z, jy[:,:,toroidalPlotIndex])
cbar = plt.colorbar()
cbar.set_label("current density / A/m^2")
plt.axis("equal")
plt.xlabel("r / m")
plt.ylabel("z / m")
plt.title("j_y only plasma at phi=" + np.str(phi) + "°")
plt.tight_layout();

plt.figure()
plt.pcolormesh(extended_r, extended_z, jz[:,:,toroidalPlotIndex])
cbar = plt.colorbar()
cbar.set_label("current density / A/m^2")
plt.axis("equal")
plt.xlabel("r / m")
plt.ylabel("z / m")
plt.title("j_z only plasma at phi=" + np.str(phi) + "°")
plt.tight_layout();

#%%

from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.gca(projection='3d')

def axisEqual3D(ax):
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)

radialIndex = 66
eachPoloidal = 2

#
for t in np.arange(0, 180, 1):
#t=toroidalPlotIndex

    
    xq = x[radialIndex,::eachPoloidal,t]
    yq = y[radialIndex,::eachPoloidal,t]
    zq = z[radialIndex,::eachPoloidal,t]
    
    jxq = jx[radialIndex,::eachPoloidal,t]
    jyq = jy[radialIndex,::eachPoloidal,t]
    jzq = jz[radialIndex,::eachPoloidal,t]
    
    modJ = np.sqrt(jxq*jxq + jyq*jyq + jzq*jzq)
    
    #arrows=ax.quiver3D(xq, yq, zq, jxq, jyq, jzq, length=0.2, pivot="tail", cmap=plt.cm.jet, linewidth=2)
    arrows=ax.quiver3D(xq, yq, zq, jxq, jyq, jzq, length=0.3, pivot="tail", cmap=plt.cm.jet, linewidth=2, arrow_length_ratio = 0.0)
    arrows.set_array(modJ)


cb=fig.colorbar(arrows, cmap=plt.cm.jet)
    

            
axisEqual3D(ax)
ax.set_xlabel("x / m")
ax.set_ylabel("y / m")
ax.set_zlabel("z / m")          

