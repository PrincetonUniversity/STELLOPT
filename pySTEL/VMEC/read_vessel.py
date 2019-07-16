# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 17:17:57 2019

@author: gawe
"""
# ===================================================================== #
# ===================================================================== #

from __future__ import absolute_import, with_statement, absolute_import, \
                       division, print_function, unicode_literals

import numpy as _np
try:
    from utils import Struct, fgetl, sscanf
except:
    from .utils import Struct, fgetl, sscanf
# end try

# ===================================================================== #
# ===================================================================== #


def read_vessel(filename=None):
    """
    READ_VESSEL(filename) Reads the vessel data into a structure
       READ_VESSEL reads a vessel data file.  The vessel data file is a file
       with the format:

        MACHINE:  LHD
        DATE:  08-12-11

        0.5094E+01  0.1016E+01  0.0000E+00
        0.5065E+01  0.1006E+01  0.0000E+00
        0.5037E+01  0.9960E+00  0.0000E+00

       Here the string after MACHINE is the name of the machine.  The DATE
       string contains information about the date of the file.  The array that
       follows contains a series of datapoints (R, Z, phi) [m, m, rad].  These
       datapoints define a piecewise defined limiter surface at each phi.

       data:
           datatype:   Identifies the structure type. 'vessel'
           machine:    Machine string from file.
           date:       Date string from file.
           coords:     Coordinate data (R,Z,Phi,Face)

       Example:
           ves_data=read_vessel('vessel.dat')

       See also plot_vessel.

       Original version (from matlabVMEC)
           Written by:     S.Lazerson (lazerson@pppl.gov)
           Version:        1.0
           Date:           1/11/11

       Ported to python on July 16th, 2019 by Gavin Weir <gavin.weir@ipp.mpg.de>
    """
    data = Struct()
    # Check arguments
    if filename is None:
        print('ERROR: read_vessel requires filename')
        return data
    # end if

    # Define datatype
    data.datatype='vessel'

    # Open File
    try:
        fid = open(filename,'r')

        # Read Header
        header_line1 = fgetl(fid)
        header_line2 = fgetl(fid)

        # Parse the header
        inds = _np.asarray(range(1, len(header_line1)), dtype=int)
        data.machine = header_line1[header_line1.find(':')+inds].strip()

        inds = _np.asarray(range(1, len(header_line2)), dtype=int)
        data.date = header_line2[header_line2.find(':')+inds].strip()

        # Scan through file
        first=1;
        group=1;
        for line in fid:
            line = fgetl(fid)
            if len(line) == 0:
                temp = sscanf(line, '%e %e %e', [1, 3])
                r = temp[0]
                z = temp[1]
                phi = temp[2]
                if first:
                    data.coords = _np.asarray([r, z, phi, group])
                    first = 0
                    groupphi = phi
                    rfirst = r
                    zfirst = z
                    phifirst = phi
                else:
                    if phi != groupphi:
                        # First add first point as last
                        temp = _np.hstack((rfirst, zfirst, phifirst, group))
                        data.coords = _np.vstack((data.coords, temp))

                        # Now set current point as first
                        rfirst = r
                        zfirst = z
                        phifirst = phi

                        # Update group number and groupphi
                        group = group + 1
                        groupphi = phi
                    # end if
                    temp = _np.hstack((r, z, phi, group))
                    data.coords = _np.vstack((data.coords, temp))
                # end if
            # end if
        # end while

        # Add last point
        temp = _np.hstack((rfirst, zfirst, phifirst, group))
        data.coords = _np.vstack((data.coords, temp))

        fid.close()
    except:
        raise
    finally:
        try:    fid.close()
        except: pass
    # end try
    return data
# end def read_vessel


# ===================================================================== #
# ===================================================================== #

