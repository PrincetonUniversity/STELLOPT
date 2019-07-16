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
from utils import Struct, fgetl, fgets, sscanf

# ====================================================================== #
# ====================================================================== #


def read_vmec_mercier(filname):
    f = Struct()
    try:
#    with open(filname, 'r') as fid: # fid=fopen(filname,'r')
        fid = open(filname, 'r')
        fgetl(fid) # First Header
        fgetl(fid) # Second Header ----

        # Read first line
        line = fgetl(fid)
        val  = sscanf(line,'%e')
        f.s     = val[0]
        f.phi   = val[1]
        f.iota  = val[2]
        f.shear = val[3]
        f.vp    = val[4]
        f.well  = val[5]
        f.itor  = val[6]
        f.ditor = val[7]
        f.pres  = val[8]
        f.dpres = val[9]

        line = fgets(fid)   # includes end of line symbol \n
        while line!='':  # !strcmp(line, ''):
            val = sscanf(line.replace('\n',''),'%e')
            f.s     = _np.vstack( (f.s,     val[0]) )
            f.phi   = _np.vstack( (f.phi,   val[1]) )
            f.iota  = _np.vstack( (f.iota,  val[2]) )
            f.shear = _np.vstack( (f.shear, val[3]) )
            f.vp    = _np.vstack( (f.vp,    val[4]) )
            f.well  = _np.vstack( (f.well,  val[5]) )
            f.itor  = _np.vstack( (f.itor,  val[6]) )
            f.ditor = _np.vstack( (f.ditor, val[7]) )
            f.pres  = _np.vstack( (f.pres,  val[8]) )
            f.dpres = _np.vstack( (f.dpres, val[9]) )
            line=fgetl(fid)
        # end
        fgetl(fid)
        fgetl(fid)

        # Read first line
        line = fgetl(fid)
        val  = sscanf(line,'%e')
        # first element unused?
        f.dmerc  = val[1]
        f.dshear = val[2]
        f.dcurr  = val[3]
        f.dwell  = val[4]
        f.dgeod  = val[5]

        while True:  # not feof(fid):
            try:
                line = fgetl(fid)
            except EOFError:
                break
            except:
                raise
            # end try
            val  = sscanf(line,'%e')
            f.dmerc  = _np.vstack( (f.dmerc, val[0]) )
            f.dshear = _np.vstack( (f.dshear, val[1]) )
            f.dcurr  = _np.vstack( (f.dcurr, val[2]) )
            f.dwell  = _np.vstack( (f.dwell, val[3]) )
            f.dgeod  = _np.vstack( (f.dgeod, val[4]) )
        # end
    # implicitly close the file using the with open (otherwise explicitly)
    except:
        raise
    finally:
        fid.close()
    return f
# end

# ====================================================================== #
# ====================================================================== #




