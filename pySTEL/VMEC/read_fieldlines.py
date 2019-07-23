# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 21:46:58 2019

@author: weir
"""

import numpy as _np

try:
    from VMEC.read_hdf5 import read_hdf5
except:
    from .read_hdf5 import read_hdf5
# end try

def read_fieldlines(filename, verbose=True):
    """
    READ_FIELDLINES Reads the HDF5 file created by FIELDLINES
     This funciton reads the fieldlines file and returns the data from the
     file in a structure.

     Example usage
          data=read_fieldlines('fieldline_test.h5');  % Reads FIELDLINE HDF5 file

     Maintained by: Samuel Lazerson (lazerson@pppl.gov)
     Version:       1.2
     """
    data = read_hdf5(filename)
    data.datatype='FIELDLINES'
    data.X_lines=data.R_lines*_np.cos(data.PHI_lines)
    data.Y_lines=data.R_lines*_np.sin(data.PHI_lines)
    data.phiend = data.PHI_lines[:,data.nsteps-1]
    data.dphi = float(data.phiend)/float(data.nsteps-1)

    # Do this for VMECplot
    data.ns = data.nlines
    data.mpol = 0
    data.nu = 2
    data.ntor = 0
    data.nzeta_grid = _np.copy(data.npoinc)
    if hasattr(data,'phiaxis'):
        data.nfp= round(2*_np.pi/max(data.phiaxis))
    elif hasattr(data,'phi_grid'):
        data.nfp= _np.round(2*_np.pi/max(data.phi_grid))
    # end if
    data.input_extension=filename

    # If we have B_PHI make B_R and B_Z true B_R and B_Z
    if hasattr(data,'B_PHI'):
        # r_temp   = repmat(repmat(data.raxis,[1 data.nphi]),[1 1 data.nz])
        r_temp = _np.tile( _np.tile(data.raxis, (1, data.nphi)), (1,1,data.nz))
        data.B_R = data.B_R*data.B_PHI/r_temp
        data.B_Z = data.B_Z*data.B_PHI/r_temp
    # end if

    # Calculate rotational transform
    x = _np.zeros((data.R_lines.shape[0],data.nsteps), dtype=float) # _np.zeros_like(data.R_lines)
    y = _np.zeros((data.Z_lines.shape[0],data.nsteps), dtype=float) # _np.zeros_like(data.Z_lines)
    for ii in range(data.nsteps):  #i=1:data.nsteps
        x[:,ii] = data.R_lines[:,ii]-data.R_lines[0,ii]
        y[:,ii] = data.Z_lines[:,ii]-data.Z_lines[0,ii]
    # end for
    theta = _np.arctan2(y,x)
    dtheta = _np.diff(theta.T).T
    dtheta[dtheta<-_np.pi] = dtheta[dtheta<-_np.pi]+2*_np.pi
    dtheta[dtheta>_np.pi] = dtheta[dtheta>_np.pi]-2*_np.pi
    theta = _np.cumsum(dtheta.T).T
    # theta=[zeros(data.nlines,1) theta]
    theta = _np.hstack((_np.zeros((data.nlines,1), dtype=float), theta))
    for jj in range(1, data.nlines):  # j=2:data.nlines
        # f0 = fit(data.PHI_lines[jj,:].T,theta[jj,:].T,'poly1')
        f0 = _np.poly1d(data.PHI_lines[jj,:].T, theta[jj,:].T)
        ci = confint(f0)-f0.p1
        ci = ci[:,0]
        data.iota[jj] = f0.p1
        data.iota_err[jj] = max(abs(ci))
    # end for
    data.rho=_np.sqrt(x[:,0]**2.0+y[:,0]**2.0)
    return data
# end def fieldlines
