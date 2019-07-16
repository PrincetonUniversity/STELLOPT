# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 17:10:40 2019

@author: gawe
"""
# ======================================================================== #
# ======================================================================== #

from __future__ import absolute_import, with_statement, absolute_import, \
                       division, print_function, unicode_literals

import numpy as _np
from utils import Struct, fgetl, sscanf

# ====================================================================== #
# ====================================================================== #


def read_vmec_jxbout(filname):
    f = Struct()
    try:
        fid = open(filname, 'r')

        fgetl(fid)           # Blank Line
        line = fgetl(fid)    # Radial Poloidal Toroidal points
        val  = sscanf(line,'%*20c %d %*24c %d %*24c %d')

        f.nrad   = val[0]
        f.ntheta = val[1]
        f.nzeta  = val[2]

        line = fgetl(fid)     # mpol ntor
        val  = sscanf(line,'%*18c %d %*18c %d')
        f.mpol = val[0]
        f.ntor = val[1]
        for ii in range( 13 ): # i=1:13      %Header stuff
            fgetl(fid)
        # end

        line = fgetl(fid)    # Values
        val  = sscanf(line,'%*18c %e %*17c %e %*11c %e %*17c %e')
        f.tor_flux = val[0]
        f.fnorm    = val[1]
        f.jdotb    = val[2]
        f.bdotv    = val[3]

        line = fgetl(fid)    # percentages
        fgetl(fid)
        fgetl(fid)
        fgetl(fid)

        f.data = _np.zeros( (f.nrad,f.nzeta,f.ntheta,13), dtype=_np.float64)
        for ii in range(f.nrad): # i=1:nrad:
            for jj in range(f.nzeta): # j=1:f.nzeta
                line=fgetl(fid)    # Angle Information
                for kk in range(f.ntheta): # k=1:f.ntheta
                    line = fgetl(fid)    # Data
                    f.data[ii, jj, kk, :] = sscanf(line,'%d %12e')
                # end
            # end
        # end
    # implicitly close the file using the with open (otherwise explicitly)
    except:
        raise
    finally:
        fid.close()
    return f
# end  def read_vmec_jxbout

# ========================================================================= #
# ========================================================================= #

