# -*- coding: utf-8 -*-
"""
    READ_VMEC(filname) This class reads the VMEC wout file.
    This function reads the wout file and returns the data from the file
    in a structure.  Based on the 'readw_only_priv.f' subroutine.
    In addtion to the raw variables, the structure also contains the
    reconstituted Fourier arrays in their full form.
    Currently this function can read VMEC files up to version 8+ (netCDF and
    text).

    Example usage
    data = read_vmec('wout.test')      # Reads VMEC wout file
    mdata = read_vmec('mercier.test')  # Reads VMEC mercier file

    Original version written in MATLAB
        Maintained by: Samuel Lazerson (lazerson@pppl.gov)
        Version:       1.8
        available through the MathWorks file exchange under matlabVMEC:
        https://www.mathworks.com/matlabcentral/fileexchange/29031-matlabvmec

    Ported to python and updated
        Maintained by: Gavin Weir (gavin.weir@ipp.mpg.de)
        Version:       0.99
        available through the pySTELL library in STELLOPT
        https://bitbucket.org/lazerson_princeton/stellopt/wiki/Home

    NOTES:
    01/05/2011      Modified output variables to use cos sine (c/s) notation
                    instead of the e notation.  Also added modifications to
                    support non-axisymmetric VMEC runs.
    01/13/2011     Overloaded to read vmec mercier file.
    01/31/2011     All quantities now mapped to full mesh.
    02/01/2011     Updated for version 8.47
    02/28/2011     Properly Handles half grid quantities (see wrfcn in pgplot)
    05/31/2011     Added support for +8.0 text files
    03/21/2012     Modified so opening netCDF via path is possible.
                   Fixed issue with mu constant when reading netCDF files.

    Original branch point for porting to python:
    03/03/2016     G.M. Weir (IPP-Greifswald) Ported to Python (old version)

    1/30/13     Added calculation of chip
    1/10/14     Modified to read 8.51 files
                Uses new methods for calculating J taking into acount the
                odd modes properly.
    3/31/14     Fixed calculation of non-stellarator symmetric J terms.
    3/1/16      Corrected variable names so text files use the new method
                calculation of current densities.

    Updated merge in port to python:
    07/12/2019     G.M. Weir updated for general use in python
"""
# ======================================================================== #
# ======================================================================== #

from __future__ import absolute_import, with_statement, absolute_import, \
                       division, print_function, unicode_literals

import numpy as _np
from utils import Struct, fscanf, fgetl

# ======================================================================== #
# ======================================================================== #




# ======================================================================== #
# ======================================================================== #


def return_fmts(fmt):
    delim = ''
    if fmt.find(',')>-1:
        delim = ','
    # end if
    fmtg = ''.join(['%g'+delim for _ in range(1)]) # '%g'
    fmtd = ''.join(['%d'+delim for _ in range(1)]) # '%d'
    fmts = ''.join(['%s'+delim for _ in range(1)]) # '%s'
    fmt2 = ''.join(['%g'+delim for _ in range(2)]) # '%g%g'
    fmt3 = ''.join(['%g'+delim for _ in range(3)]) # '%g%g%g'
    fmt6 = ''.join(['%g'+delim for _ in range(6)]) # '%g%g%g%g%g%g'
    fmt7 = ''.join(['%g'+delim for _ in range(7)]) # '%g%g%g%g%g%g%g'
    fmt10 = ''.join(['%g'+delim for _ in range(10)]) # '%g%g%g%g%g%g%g%g%g%g'
    fmt11 = ''.join(['%g'+delim for _ in range(11)]) # '%g%g%g%g%g%g%g%g%g%g%g'
    fmt12 = ''.join(['%g'+delim for _ in range(12)]) # '%g%g%g%g%g%g%g%g%g%g%g%g'
    fmt13 = ''.join(['%g'+delim for _ in range(13)]) # '%g%g%g%g%g%g%g%g%g%g%g%g%g'
    fmt14 = ''.join(['%g'+delim for _ in range(14)]) # '%g%g%g%g%g%g%g%g%g%g%g%g%g%g'
    fmt20 = ''.join(['%g'+delim for _ in range(20)]) # '%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g'

    fmt2_3 = (''.join(['%d'+delim for _ in range(2)]) # '%d%d%g%g%g'
           + ''.join(['%g'+delim for _ in range(3)]) )
    fmt2_6 = (''.join(['%d'+delim for _ in range(2)]) # '%d%d%g%g%g%g%g%g'
           + ''.join(['%g'+delim for _ in range(6)]) )
    fmt2_7 = (''.join(['%d'+delim for _ in range(2)]) # '%d%d%g%g%g%g%g%g%g'
            + ''.join(['%g'+delim for _ in range(7)]) )
    fmt2_11 = (''.join(['%d'+delim for _ in range(2)]) # '%d%d%g%g%g%g%g%g%g%g%g%g%g'
           + ''.join(['%g'+delim for _ in range(11)]) )
    fmt2_14 = (''.join(['%d'+delim for _ in range(2)]) # '%d%d%g%g%g%g%g%g%g%g%g%g%g%g%g%g'
             + ''.join(['%g'+delim for _ in range(14)]) )

    fmt2_3_2_7 = (''.join(['%d'+delim for _ in range(2)]) # '%d%d%g%g%g%d%d%g%g%g%g%g%g%g'
                + ''.join(['%g'+delim for _ in range(3)])
                + ''.join(['%d'+delim for _ in range(2)])
                + ''.join(['%g'+delim for _ in range(7)]) )
    fmt2_6_2_14 = (''.join(['%d'+delim for _ in range(2)]) # '%d%d%g%g%g%g%g%g%d%d%g%g%g%g%g%g%g%g%g%g%g%g%g%g'
                + ''.join(['%g'+delim for _ in range(6)])
                + ''.join(['%d'+delim for _ in range(2)])
                + ''.join(['%g'+delim for _ in range(14)]) )

    return (fmtg, fmtd, fmts, fmt2, fmt3, fmt6, fmt7, fmt10, fmt11, fmt12,
            fmt13, fmt14, fmt20, fmt2_3, fmt2_6, fmt2_7, fmt2_11, fmt2_14,
            fmt2_3_2_7, fmt2_6_2_14)

    # =============================================== #

def read_vmec_orig(fid, fmt):
    # For 5.10 and earlier
    f = Struct()

    # Unit Conversions
#    mu0 = 4.0*_np.pi*1e-7
    dmu0=(2.0e-7)/_np.pi

    fmtg, fmtd, fmts, fmt2, fmt3, fmt6, fmt7, fmt10, fmt11, fmt12, \
        fmt13, fmt14, fmt20, fmt2_3, fmt2_6, fmt2_7, fmt2_11, fmt2_14, \
        fmt2_3_2_7, fmt2_6_2_14 = return_fmts(fmt)

    # =============================================== #

    # Read Data
    f.wb, f.wp, f.gamma, f.pfac, f.nfp, f.ns, f.mpol, f.ntor, f.mnmax, f.itfsq, \
        f.niter, f.iasym, f.ireconstruct = tuple(fscanf(fid, fmtg, 13).tolist())

    f.imse, f.itse, f.nbsets, f.nobd, f.nextcur = tuple(fscanf(fid, fmtd, 5).tolist())

    f.nstore_seq = 100
    # Error Check
    if f.ierr_vmec and (f.ierr_vmec != 4):
        print('ierr_vmec >0')
        return None
    # end if

    # Read nbfld
    if (f.nbsets > 0):
        f.nbfld = fscanf(fid, fmtg, f.nbsets)
    # end if

    # Read mgrid filename and setup other format statements
    f.mgrid_file = fscanf(fid, fmts, 1)

    # Read Arrays
    if f.iasym > 0:
        data1 = fscanf(fid, fmt2_14, [16, f.mnmax])
        data = fscanf(fid, fmt14, [14, f.mnmax*(f.ns-1)])
    else:
        data1 = fscanf(fid, fmt2_11, [13, f.mnmax])
        data = fscanf(fid, fmt11, [11, f.mnmax*(f.ns-1)])
    # end if

    #Extract Data from Arrays
    f.xm = _np.asarray(data1[0][:])
    f.xn = _np.asarray(data1[1][:])

    # Extract the data and reshape
#    data1 = _np.asarray(data1[2:][:])
    f.rmnc, f.zmns, f.lmns, f.bmn, f.gmn, f.bsubumn, f.bsubvmn, f.bsubsmn, \
        f.bsupumn, f.bsupvmn, f.currvmn \
        = tuple([_np.vstack((_np.copy(data1[ii+2][:].T), _np.reshape(data[ii][:], f.mnmax, f.ns-1, order='F').copy()))
            for ii in range(10)])
#        = tuple([_np.vstack((_np.copy(data1[ii][:].T), _np.reshape(data[ii][:], f.mnmax, f.ns-1, order='F').copy()))
#            for ii in range(10)])

    #Read the half-mesh quantities
#    data = fscanf(fid, fmt12, [12, f.ns/2])

    f.iotas, f.mass, f.pres, f.phip, f.buco, f.bvco, f.phi, f.vp, f.overr, \
        f.juru, f.jcurv, f.specw = tuple(fscanf(fid, fmt12, [12, f.ns/2]).tolist())
#        f.juru, f.jcurv, f.specw = tuple(map(_np.asarray, zip(*data)))

    data = fscanf(fid, fmt6, 6)
    f.aspect  = data[0]
    f.betatot = data[1]
    f.betapol = data[2]
    f.betator = data[3]
    f.betaxis = data[4]  # betaaxis?
    f.b0      = data[5]

    #Mercier Criterion
    data = fscanf(fid, fmt6, [6, f.ns-2])
    f.Dmerc  = data[0][:]
    f.Dshear = data[1][:]
    f.Dwell  = data[2][:]
    f.Dcurr  = data[3][:]
    f.Dgeod  = data[4][:]
    f.equif  = data[5][:]
    if (f.nextcur > 0):
        f.extcur = fscanf(fid, fmt, f.nextcur)
        f.curlabel = fscanf(fid, fmt, f.nextcur)
    #end if

    data = fscanf(fid, fmt2, [2, f.nstore_seq])
    f.sqt  = data[0][:]
    f.wdot = data[1][:]

    #Convert from Internal Units to Physical Units
    f.mass  /= dmu0
    f.pres  /= dmu0
    f.jcuru /= dmu0
    f.jcurv /= dmu0
    f.jdotb /= dmu0
    f.phi   *= -1.0       # Data and MSE Fits

    if (f.ireconstruct > 0):
        if (f.imse >= 2) or (f.itse >0):
            f.twsgt   = fscanf(fid, fmt, 1)
            f.msewgt  = fscanf(fid, fmt, 1)
            f.isnodes = fscanf(fid, fmt, 1)

            data = fscanf(fid, fmt3, [3, f.isnodes])
            f.sknots  = data[0][:]
            f.ystark  = data[1][:]
            f.y2stark = data[2][:]

            f.ipnodes = fscanf(fid, fmt, 1)

            data = fscanf(fid, fmt3, [3, f.ipnodes])
            f.pknots = data[0][:]
            f.ythom  = data[1][:]
            f.y2thom = data[2][:]

            data = fscanf(fid, fmt7, [7, (2*f.ns)-1])
            f.anglemse = data[0][:]
            f.rmid     = data[1][:]
            f.qmid     = data[2][:]
            f.shear    = data[3][:]
            f.presmid  = data[4][:]
            f.alfa     = data[5][:]
            f.curmid   = data[6][:]

            data = fscanf(fid, fmt3, [3, f.imse])
            f.rstark    = data[0][:]
            f.datastark = data[1][:]
            f.qmeas     = data[2][:]

            data = fscanf(fid, fmt2, [2, f.itse])
            f.rthom    = data[0][:]
            f.datathom = data[1][:]
        # end

        if (f.nobd > 0):
            data = fscanf(fid, fmt3, [3, f.nobd])
            f.dsiext = data[0][:]
            f.plflux = data[1][:]
            f.dsiobt = data[2][:]

            f.flmwgt = fscanf(fid, fmt, 1)
        # end

        nbfldn = _np.sum( f.nbfld[:f.nbsets] )
#        nbfldn = _np.sum( nbldf[:f.nbsets] )
        if (nbfldn > 0):
            for nn in range(f.nbsets): # n=1:nbsets
                data = fscanf(fid, fmt3, [3, f.nbfld[nn]])
                f.bcoil[:, nn]  = data[0][:]
                f.plbfld[:, nn] = data[1][:]
                f.bbc[:, nn]    = data[2][:]
            # end
            f.bcwgt = fscanf(fid, fmt, 1)
        # end
        f.phidiam = fscanf(fid, fmt, 1)
        f.delphid = fscanf(fid, fmt, 1)

        #Read Limiter and Prout Plotting Specs
        f.nsets   = fscanf(fid, fmt, 1)
        f.nparts  = fscanf(fid, fmt, 1)
        f.nlim    = fscanf(fid, fmt, 1)
        f.nsetsn  = fscanf(fid, fmt, f.nsets)
        f.pfcspec = _np.zeros( (f.nparts,_np.max(f.nxsetsn),f.nsets), dtype=_np.float64)

        for kk in range(f.nsets): # k=1:f.nsets
            for jj in range(f.nsetsn[kk]): # j=1:f.nsetsn(k)
                for ii in range(f.nparts): # i=1:f.nparts
                    f.pfcspec[ii, jj, kk] = fscanf(fid, fmt, 1)
                # end
            # end
        # end

        f.limitr = fscanf(fid, fmt, f.nlim)
        f.rlim   = _np.zeros( (_np.max(f.limitr),f.nlim), dtype=_np.float64)
        f.zlim   = _np.zeros_like(f.rlim)

        for jj in range(f.nlim): # j=1:f.nlim
            for ii in range(f.limitr[jj]): # i=1:f.limitr(j)
                data = fscanf(fid,fmt2,2)
                f.rlim[ii, jj] = data[0]
                f.zlim[ii, jj] = data[1]
            # end
        # end

        f.nrgrid = fscanf(fid, fmt, 1)
        f.nzgrid = fscanf(fid, fmt, 1)

        f.tokid = fscanf(fid, fmt, 1)
        f.rx1   = fscanf(fid, fmt, 1)
        f.rx2   = fscanf(fid, fmt, 1)
        f.zy1   = fscanf(fid, fmt, 1)
        f.zy2   = fscanf(fid, fmt, 1)
        f.conif = fscanf(fid, fmt, 1)
        f.imatch_phiedge = fscanf(fid, fmt, 1)
    # end
    return f
# end  def read_vmec_orig    # checked vs matlabVMEC on July15, 2019

# ========================================================================= %

def read_vmec_605(fid,fmt):
    #For 6.05
    f = Struct()

    # Unit Conversions
    # mu0 = 4.0*_np.pi*1e-7
    dmu0=(2.0e-7)/_np.pi

    # =============================================== #

    fmtg, fmtd, fmts, fmt2, fmt3, fmt6, fmt7, fmt10, fmt11, fmt12, \
        fmt13, fmt14, fmt20, fmt2_3, fmt2_6, fmt2_7, fmt2_11, fmt2_14, \
        fmt2_3_2_7, fmt2_6_2_14 = return_fmts(fmt)

    # =============================================== #

    #Read Data
    data = fscanf(fid,fmtg,6)
    f.wb        = data[0]
    f.wp        = data[1]
    f.gamma     = data[2]
    f.pfac      = data[3]
    f.rmax_surf = data[4]
    f.rmin_surf = data[5]

    data = fscanf(fid,fmtd,10)
    f.nfp          = data[0]
    f.ns           = data[1]
    f.mpol         = data[2]
    f.ntor         = data[3]
    f.mnmax        = data[4]
    f.itfsq        = data[5]
    f.niter        = data[6]
    f.iasym        = data[7]
    f.ireconstruct = data[8]
    f.ierr_vmec    = data[9]

    data = fscanf(fid,fmtd,5)
    f.imse       = data[0]
    f.itse       = data[1]
    f.nbsets     = data[2]
    f.nobd       = data[3]
    f.nextcur    = data[4]
    f.nstore_seq = 100

    #Error Check
    if f.ierr_vmec and (f.ierr_vmec != 4):
        print('ierr_vmec >0')
        return
    # end

    #Read nbfld
    if (f.nbsets > 0):
        f.nbfld = fscanf(fid, fmtg, f.nbsets)
    # end

    #Read mgrid filename
    f.mgrid_file = fscanf(fid, fmts, 1)

    #Read Arrays
    if f.iasym > 0:
        data1 = fscanf(fid,fmt2_14,[16, f.mnmax])
        data  = fscanf(fid,fmt14,[14, f.mnmax*(f.ns-1)])
    else:
        data1 = fscanf(fid,fmt2_11,[13, f.mnmax])
        data  = fscanf(fid,fmt11,[11, f.mnmax*(f.ns-1)])
    # end

    #Extract Data from Arrays
    f.xm = _np.asarray(data1[0][:])
    f.xn = _np.asarray(data1[1][:])

    # Extract the data and reshape
#    data1 = _np.asarray(data1[2:][:])
    f.rmnc, f.zmns, f.lmns, f.bmnc, f.gmnc, f.bsubumnc, f.bsubvmnc, f.bsubsmnc, \
        f.bsupumnc, f.bsupvmnc, f.currvmnc \
        = tuple([_np.vstack((_np.copy(data1[ii+2][:].T), _np.reshape(data[ii][:], f.mnmax, f.ns-1, order='F').copy()))
            for ii in range(10)])
#        = tuple([_np.vstack((_np.copy(data1[ii][:].T), _np.reshape(data[ii][:], f.mnmax, f.ns-1, order='F').copy()))
#            for ii in range(10)]):

    #Read the half-mesh quantities
    data = fscanf(fid, fmt12, [12, f.ns/2])
    f.iotas = data[0][:]
    f.mass  = data[1][:]
    f.pres  = data[2][:]
    f.phip  = data[3][:]
    f.buco  = data[4][:]
    f.bvco  = data[5][:]
    f.phi   = data[6][:]
    f.vp    = data[7][:]
    f.overr = data[8][:]
    f.jcuru = data[9][:]
    f.jcurv = data[10][:]
    f.specw = data[11][:]

    data = fscanf(fid, fmt6, 6)
    f.aspect  = data[0]
    f.betatot = data[1]
    f.betapol = data[2]
    f.betator = data[3]
    f.betaxis = data[4]
    f.b0      = data[5]

    #Mercier Criterion
    data = fscanf(fid, fmt6, [6, f.ns-2])
    f.Dmerc  = data[0][:]
    f.Dshear = data[1][:]
    f.Dwell  = data[2][:]
    f.Dcurr  = data[3][:]
    f.Dgeod  = data[4][:]
    f.equif  = data[5][:]

    if (f.nextcur > 0):
        f.extcur = fscanf(fid, fmtg, f.nextcur)
        f.curlabel = fscanf(fid, fmtg, f.nextcur)
    # end

    data = fscanf(fid, fmt2, [2, f.nstore_seq])
    f.sqt = data[0][:]
    f.wdot = data[1][:]

    #Convert from Internal Units to Physical Units
    f.mass /= dmu0
    f.pres /= dmu0
    f.jcuru /= dmu0
    f.jcurv /= dmu0
    f.jdotb /= dmu0
    f.phi   *= -1.0    #Data and MSE Fits

    if (f.ireconstruct > 0):
        if (f.imse >= 2) or (f.itse >0):
            f.twsgt = fscanf(fid, fmtg, 1)
            f.msewgt = fscanf(fid, fmtg, 1)
            f.isnodes = fscanf(fid, fmtd, 1)

            data = fscanf(fid, fmt3, [3, f.isnodes])
            f.sknots = data[0][:]
            f.ystark = data[1][:]
            f.y2stark = data[2][:]
            f.ipnodes = fscanf(fid, fmtd, 1)

            data = fscanf(fid, fmt3, [3, f.ipnodes])
            f.pknots = data[0][:]
            f.ythom = data[1][:]
            f.y2thom = data[2][:]

            data = fscanf(fid, fmt7, [7, (2*f.ns)-1])
            f.anglemse = data[0][:]
            f.rmid = data[1][:]
            f.qmid = data[2][:]
            f.shear = data[3][:]
            f.presmid = data[4][:]
            f.alfa = data[5][:]
            f.curmid = data[6][:]

            data = fscanf(fid, fmt3, [3, f.imse])
            f.rstark = data[0][:]
            f.datastark = data[1][:]
            f.qmeas = data[2][:]

            data = fscanf(fid, fmt2, [2, f.itse])
            f.rthom = data[0][:]
            f.datathom = data[1][:]
        # end

        if (f.nobd > 0):
            data = fscanf(fid, fmt3, [3, f.nobd])
            f.dsiext = data[0][:]
            f.plflux = data[1][:]
            f.dsiobt = data[2][:]
            f.flmwgt = fscanf(fid, fmtg, 1)
        # end
        nbfldn=_np.sum(f.nbfld[:f.nbsets])

        if (nbfldn > 0):
            for nn in range(f.nbsets): # n=1:nbsets
                data = fscanf(fid, fmt3, [3, f.nbfld[nn]])
                f.bcoil[:, nn] = data[0][:]
                f.plbfld[:, nn] = data[1][:]
                f.bbc[:, nn] = data[2][:]
            # end
            f.bcwgt = fscanf(fid, fmtg, 1)
        # end
        f.phidiam = fscanf(fid, fmtg, 1)
        f.delphid = fscanf(fid, fmtg, 1)

        #Read Limiter and Prout Plotting Specs
        f.nsets = fscanf(fid, fmtg, 1)
        f.nparts = fscanf(fid, fmtg, 1)
        f.nlim = fscanf(fid, fmtg, 1)
        f.nsetsn = fscanf(fid, fmtg, f.nsets)
        f.pfcspec=_np.zeros( (f.nparts, _np.max(f.nxsetsn), f.nsets), dtype=_np.float64)

        for kk in range(f.nsets): # k=1:f.nsets
            for jj in range(f.nsetsn[kk]): # j=1:f.nsetsn(k)
                for ii in range(f.nparts): # i=1:f.nparts
                    f.pfcspec[ii, jj, kk] = fscanf(fid, fmtg, 1)
                # end
            # end
        # end

        f.limitr = fscanf(fid, fmtg, f.nlim)
        f.rlim = _np.zeros( (_np.max(f.limitr), f.nlim), dtype=_np.float64)
        f.zlim = _np.zeros( (_np.max(f.limitr), f.nlim), dtype=_np.float64)

        for jj in range(f.nlim): # j=1:f.nlim
            for ii in range(f.limitr[jj]): #i=1:f.limitr(j)
                data = fscanf(fid, fmt2, 2)
                f.rlim[ii, jj] = data[1]
                f.zlim[ii, jj] = data[2]
            # end
        # end
        f.nrgrid = fscanf(fid, fmtg, 1)
        f.nzgrid = fscanf(fid, fmtg, 1)
        f.tokid = fscanf(fid, fmtg, 1)
        f.rx1 = fscanf(fid, fmtg, 1)
        f.rx2 = fscanf(fid, fmtg, 1)
        f.zy1 = fscanf(fid, fmtg, 1)
        f.zy2 = fscanf(fid, fmtg, 1)
        f.conif = fscanf(fid, fmtg, 1)
        f.imatch_phiedge = fscanf(fid, fmtg, 1)
    # end
    return f
# end  def read_vmec_605(fid,fmt)   # checked against matlabVMEC on July15,2019

# ========================================================================= #


def read_vmec_620(fid,fmt):
    #For 6.20
    f = Struct()

    #Unit Conversions
#    dmu0=(2.0e-7)/_np.pi

    # =============================================== #

    fmtg, fmtd, fmts, fmt2, fmt3, fmt6, fmt7, fmt10, fmt11, fmt12, \
        fmt13, fmt14, fmt20, fmt2_3, fmt2_6, fmt2_7, fmt2_11, fmt2_14, \
        fmt2_3_2_7, fmt2_6_2_14 = return_fmts(fmt)

    # =============================================== #

    #Read Data
    data = fscanf(fid, fmtg, 6)
    f.wb = data[0]
    f.wp = data[1]
    f.gamma = data[2]
    f.pfac = data[3]
    f.rmax_surf = data[4]
    f.rmin_surf = data[5]

    data = fscanf(fid, fmtd, 10)
    f.nfp = data[0]
    f.ns = data[1]
    f.mpol = data[2]
    f.ntor = data[3]
    f.mnmax = data[4]
    f.itfsq = data[5]
    f.niter = data[6]
    f.iasym = data[7]
    f.ireconstruct = data[8]
    f.ierr_vmec = data[9]

    data = fscanf(fid, fmtd, 5)
    f.imse = data[0]
    f.itse = data[1]
    f.nbsets = data[2]
    f.nobd = data[3]
    f.nextcur = data[4]
    f.nstore_seq=100

    #Error Check
    if (f.ierr_vmec) and (f.ierr_vmec != 4):
        print('ierr_vmec >0')
        return f
    # end

    #Read nbfld
    if (f.nbsets > 0):
        f.nbfld = fscanf(fid, fmtg, f.nbsets)
    # end

    #Read mgrid filename
    f.mgrid_file = fscanf(fid, fmts, 1)

    #Read Arrays
    if (f.iasym > 0):
        data1 = fscanf(fid, fmt2_14, [16, f.mnmax])
        data = fscanf(fid, fmt14, [14, f.mnmax*( f.ns-1 )])
    else:
        data1 = fscanf(fid, fmt2_11, [13, f.mnmax])
        data = fscanf(fid, fmt11, [11, f.mnmax*( f.ns-1 )])
    # end

    #Extract Data from Arrays
    f.xm = _np.asarray(data1[0][:])
    f.xn = _np.asarray(data1[1][:])

    # Extract the data and reshape
#    data1 = _np.asarray(data1[2:][:])
    f.rmnc, f.zmns, f.lmns, f.bmn, f.gmn, f.bsubumnc, f.bsubvmnc, f.bsubsmns, \
        f.bsupumnc, f.bsupvmnc, f.currvmnc \
        = tuple([_np.vstack((_np.copy(data1[ii+2][:].T), _np.reshape(data[ii][:], f.mnmax, f.ns-1, order='F').copy()))
            for ii in range(10)])
#        = tuple([_np.vstack((_np.copy(data1[ii][:].T), _np.reshape(data[ii][:], f.mnmax, f.ns-1, order='F').copy()))
#            for ii in range(10)]):

    #Read the half-mesh quantities
    data = fscanf(fid, fmt13, [13, f.ns/2])
    f.iotas = data[0][:]
    f.mass = data[1][:]
    f.pres = data[2][:]
    f.beta_vol[3][:]
    f.phip = data[4][:]
    f.buco = data[5][:]
    f.bvco = data[6][:]
    f.phi = data[7][:]
    f.vp = data[8][:]
    f.overr = data[9][:]
    f.jcuru = data[10][:]
    f.jcurv = data[11][:]
    f.specw = data[12][:]

    data = fscanf(fid, fmt6, 6)
    f.aspect = data[0]
    f.betatot = data[1]
    f.betapol = data[2]
    f.betator = data[3]
    f.betaxis = data[4]
    f.b0 = data[5]

#    f.isigna = fscanf(fid, '%d\n', 1)
    f.isigna = fscanf(fid, fmtd, 1)
    f.input_extension = fgetl(fid).strip()

    data = fscanf(fid, fmtg, 8)
    f.IonLarmor = data[0]
    f.VolAvgB = data[1]
    f.RBtor0 = data[2]
    f.RBtor = data[3]
    f.Itor = data[4]
    f.Aminor = data[5]
    f.Rmajor = data[6]
    f.Volume = data[7]

    #Mercier Criterion
    data = fscanf(fid, fmt6, [6, f.ns-2])
    f.Dmerc = data[0][:]
    f.Dshear = data[1][:]
    f.Dwell = data[2][:]
    f.Dcurr = data[3][:]
    f.Dgeod = data[4][:]
    f.equif = data[5][:]

    if (f.nextcur > 0):
        f.extcur = fscanf(fid, fmtg, f.nextcur)
        f.curlabel = fscanf(fid, fmtg, f.nextcur)
    # end

    data = fscanf(fid, fmt2, [2, f.nstore_seq])
    f.sqt = data[0][:]
    f.wdot = data[1][:]

    data = fscanf(fid, fmt2, [2, f.nstore_seq])
    f.jdotb = data[0][:]
    f.bdotgradv = data[1][:]

    #Data and MSE Fits
    if (f.ireconstruct > 0):
        if ((f.imse >= 2) or (f.itse >0)):
            f.twsgt = fscanf(fid, fmtg, 1)
            f.msewgt = fscanf(fid, fmtg, 1)
            f.isnodes = fscanf(fid, fmtd, 1)

            data = fscanf(fid, fmt3, [3, f.isnodes])
            f.sknots = data[0][:]
            f.ystark = data[1][:]
            f.y2stark = data[2][:]
            f.ipnodes = fscanf(fid, fmtd, 1)

            data = fscanf(fid, fmt3, [3, f.ipnodes])
            f.pknots = data[0][:]
            f.ythom = data[1][:]
            f.y2thom = data[2][:]

            data = fscanf(fid, fmt7, [7, (2*f.ns)-1])
            f.anglemse = data[0][:]
            f.rmid = data[1][:]
            f.qmid = data[2][:]
            f.shear = data[3][:]
            f.presmid = data[4][:]
            f.alfa = data[5][:]
            f.curmid = data[6][:]

            data = fscanf(fid, fmt3, [3, f.imse])
            f.rstark = data[0][:]
            f.datastark = data[1][:]
            f.qmeas = data[2][:]

            data = fscanf(fid, fmt2, [2, f.itse])
            f.rthom = data[0][:]
            f.datathom = data[1][:]
        # end

        if (f.nobd > 0):
            data = fscanf(fid, fmt3, [3,f.nobd])
            f.dsiext = data[0][:]
            f.plflux = data[1][:]
            f.dsiobt = data[2][:]
            f.flmwgt = fscanf(fid, fmtg, 1)
        # end
        nbfldn=_np.sum(f.nbfld[:f.nbsets])

        if (nbfldn > 0):
            for nn in range(f.nbsets): # n=1:nbsets
                data = fscanf(fid, fmt3, [3, f.nbfld[nn]])
                f.bcoil[:, nn] = data[0][:]
                f.plbfld[:, nn] = data[1][:]
                f.bbc[:, nn] = data[2][:]
            # end
            f.bcwgt = fscanf(fid, fmtg, 1)
        # end
        f.phidiam = fscanf(fid, fmtg, 1)
        f.delphid = fscanf(fid, fmtg, 1)

        #Read Limiter and Prout Plotting Specs
        f.nsets = fscanf(fid, fmtg, 1)
        f.nparts = fscanf(fid, fmtg, 1)
        f.nlim = fscanf(fid, fmtg, 1)
        f.nsetsn = fscanf(fid, fmtg, f.nsets)
        f.pfcspec=_np.zeros( (f.nparts, _np.max(f.nxsetsn), f.nsets), dtype=_np.float64)
        for kk in range(f.nsets): # k=1:f.nsets
            for jj in range(f.nsetsn[kk]): # j=1:f.nsetsn(k)
                for ii in range(f.nparts): # i=1:f.nparts
                    f.pfcspec[ii, jj, kk] = fscanf(fid, fmtg, 1)
                # end
            # end
        # end
        f.limitr = fscanf(fid, fmtg, f.nlim)

        f.rlim=_np.zeros( (_np.max(f.limitr),f.nlim), dtype=_np.float64)
        f.zlim=_np.zeros_like(f.rlim)
        for jj in range(f.nlim): # j=1:f.nlim
            for ii in range(f.limitr[jj]): # i=1:f.limitr(j)
                data = fscanf(fid, fmt2, 2)
                f.rlim[ii, jj] = data[0]
                f.zlim[ii, jj] = data[1]
            # end
        # end
        f.nrgrid = fscanf(fid, fmtg, 1)
        f.nzgrid = fscanf(fid, fmtg, 1)
        f.tokid = fscanf(fid, fmtg, 1)
        f.rx1 = fscanf(fid, fmtg, 1)
        f.rx2 = fscanf(fid, fmtg, 1)
        f.zy1 = fscanf(fid, fmtg, 1)
        f.zy2 = fscanf(fid, fmtg, 1)
        f.conif = fscanf(fid, fmtg, 1)
        f.imatch_phiedge = fscanf(fid, fmtg, 1)
    # end
    return f
# end  def read_vmec_620()    # checked against matlabVMEC on July15th,2019

# ========================================================================= #


def read_vmec_650(fid, fmt):
    #function f=read_vmec_650(fid, fmt)
    #For 6.50
    f = Struct()

    #Unit Conversions
#    dmu0=(2.0e-7)/_np.pi

    # =============================================== #

    fmtg, fmtd, fmts, fmt2, fmt3, fmt6, fmt7, fmt10, fmt11, fmt12, \
        fmt13, fmt14, fmt20, fmt2_3, fmt2_6, fmt2_7, fmt2_11, fmt2_14, \
        fmt2_3_2_7, fmt2_6_2_14 = return_fmts(fmt)

    # =============================================== #

    #Read Data
    data = fscanf(fid, fmtg, 6)
    f.wb = data[0]
    f.wp = data[1]
    f.gamma = data[2]
    f.pfac = data[3]
    f.rmax_surf = data[4]
    f.rmin_surf = data[5]

    data = fscanf(fid, fmtd, 10)
    f.nfp = data[0]
    f.ns = data[1]
    f.mpol = data[2]
    f.ntor = data[3]
    f.mnmax = data[4]
    f.itfsq = data[5]
    f.niter = data[6]
    f.iasym = data[7]
    f.ireconstruct = data[8]
    f.ierr_vmec = data[9]

    data = fscanf(fid, fmtd, 6)
    f.imse = data[0]
    f.itse = data[1]
    f.nbsets = data[2]
    f.nobd = data[3]
    f.nextcur = data[4]
    f.nstore_seq = data[5]

    #Error Check
    if (f.ierr_vmec and (f.ierr_vmec != 4)):
        print('ierr_vmec >0')
        return
    # end

    #Read nbfld
    if (f.nbsets > 0):
        f.nbfld = fscanf(fid, fmtg, f.nbsets)
    # end

    #Read mgrid filename
    f.mgrid_file = fscanf(fid, fmts, 1)

    #Read Arrays
    if f.iasym > 0:
        data1 = fscanf(fid, fmt2_14, [16, f.mnmax])
        data = fscanf(fid, fmt14, [14, f.mnmax*(f.ns-1)])
    else:
        data1 = fscanf(fid, fmt2_11, [13, f.mnmax])
        data = fscanf(fid, fmt11, [11, f.mnmax*(f.ns-1)])
    # end

    #Extract Data from Arrays
    f.xm = _np.asarray(data1[0][:])
    f.xn = _np.asarray(data1[1][:])

    # Extract the data and reshape
#    data1 = _np.asarray(data1[2:][:])
    f.rmnc, f.zmns, f.lmns, f.bmnc, f.gmnc, f.bsubumnc, f.bsubvmnc, f.bsubsmns, \
        f.bsupumnc, f.bsupvmnc, f.currvmnc \
        = tuple([_np.vstack((_np.copy(data1[ii+2][:].T), _np.reshape(data[ii][:], f.mnmax, f.ns-1, order='F').copy()))
            for ii in range(10)])
#        = tuple([_np.vstack((_np.copy(data1[ii][:].T), _np.reshape(data[ii][:], f.mnmax, f.ns-1, order='F').copy()))
#            for ii in range(10)]):

    #Read the half-mesh quantities
    data = fscanf(fid, fmt13, [13, f.ns/2])
    f.iotas = data[0][:]
    f.mass = data[1][:]
    f.pres = data[2][:]
    f.beta_vol[3][:]
    f.phip = data[4][:]
    f.buco = data[5][:]
    f.bvco = data[6][:]
    f.phi = data[7][:]
    f.vp = data[8][:]
    f.overr = data[9][:]
    f.jcuru = data[10][:]
    f.jcurv = data[11][:]
    f.specw = data[12][:]

    data = fscanf(fid, fmt6, 6)
    f.aspect = data[0]
    f.betatot = data[1]
    f.betapol = data[2]
    f.betator = data[3]
    f.betaxis = data[4]
    f.b0 = data[5]

#    f.isigna = fscanf(fid, '%d\n', 1)
    f.isigna = fscanf(fid, fmtd, 1)
    f.input_extension = fgetl(fid).strip()

    data = fscanf(fid, fmtg, 8)
    f.IonLarmor = data[0]
    f.VolAvgB = data[1]
    f.RBtor0 = data[2]
    f.RBtor = data[3]
    f.Itor = data[4]
    f.Aminor = data[5]
    f.Rmajor = data[6]
    f.Volume = data[7]

    #Mercier Criterion
    data = fscanf(fid, fmt6, [6, f.ns-2])
    f.Dmerc = data[0][:]
    f.Dshear = data[1][:]
    f.Dwell = data[2][:]
    f.Dcurr = data[3][:]
    f.Dgeod = data[4][:]
    f.equif = data[5][:]

    if (f.nextcur > 0):
        f.extcur = fscanf(fid, fmtg, f.nextcur)
        f.curlabel = fscanf(fid, fmtg, f.nextcur)
    # end

    data = fscanf(fid, fmt2, [2, f.nstore_seq])
    f.sqt = data[0][:]
    f.wdot = data[1][:]

    data = fscanf(fid, fmt2, [2, f.nstore_seq])
    f.jdotb = data[0][:]
    f.bdotgradv = data[1][:]

    #Data and MSE Fits
    if (f.ireconstruct > 0):
        if ((f.imse >= 2) or (f.itse > 0)):
            f.twsgt = fscanf(fid, fmtg, 1)
            f.msewgt = fscanf(fid, fmtg, 1)
            f.isnodes = fscanf(fid, fmtd, 1)

            data = fscanf(fid, fmt3, [3, f.isnodes])
            f.sknots = data[0][:]
            f.ystark = data[1][:]
            f.y2stark = data[2][:]
            f.ipnodes = fscanf(fid, fmtd, 1)

            data = fscanf(fid, fmt3, [3, f.ipnodes])
            f.pknots = data[0][:]
            f.ythom = data[1][:]
            f.y2thom = data[2][:]

            data = fscanf(fid, fmt7, [7, (2*f.ns)-1])
            f.anglemse = data[0][:]
            f.rmid = data[1][:]
            f.qmid = data[2][:]
            f.shear = data[3][:]
            f.presmid = data[4][:]
            f.alfa = data[5][:]
            f.curmid = data[6][:]

            data = fscanf(fid, fmt3, [3, f.imse])
            f.rstark = data[0][:]
            f.datastark = data[1][:]
            f.qmeas = data[2][:]

            data = fscanf(fid, fmt2, [2, f.itse])
            f.rthom = data[0][:]
            f.datathom = data[1][:]
        # end

        if (f.nobd > 0):
            data = fscanf(fid, fmt3, [3,f.nobd])
            f.dsiext = data[0][:]
            f.plflux = data[1][:]
            f.dsiobt = data[2][:]
            f.flmwgt = fscanf(fid, fmtg, 1)
        # end

        nbfldn = _np.sum(f.nbfld[:f.nbsets])
        if (nbfldn > 0):
            for nn in range(f.nbsets):  #for n=1:nbsets
                data = fscanf(fid, fmt3, [3, f.nbfld[nn]])
                f.bcoil[:, nn] = data[0][:]
                f.plbfld[:, nn] = data[1][:]
                f.bbc[:, nn] = data[2][:]
            # end
            f.bcwgt = fscanf(fid, fmtg, 1)
        # end
        f.phidiam = fscanf(fid, fmtg, 1)
        f.delphid = fscanf(fid, fmtg, 1)

        #Read Limiter and Prout Plotting Specs
        f.nsets = fscanf(fid, fmtg, 1)
        f.nparts = fscanf(fid, fmtg, 1)
        f.nlim = fscanf(fid, fmtg, 1)
        f.nsetsn = fscanf(fid, fmtg, f.nsets)
        f.pfcspec=_np.zeros( (f.nparts,_np.max(f.nxsetsn),f.nsets), dtype=_np.float64)
        for kk in range(f.nsets):  #k=1:f.nsets
            for jj in range(f.nsetsn[kk]):  #for j=1:f.nsetsn(k)
                for ii in range(f.nparts):  #for i=1:f.nparts
                    f.pfcspec[ii, jj, kk] = fscanf(fid, fmtg, 1)
                # end
            # end
        # end

        f.limitr = fscanf(fid, fmtg, f.nlim)

        f.rlim=_np.zeros( (_np.max(f.limitr), f.nlim), dtype=_np.float64)
        f.zlim=_np.zeros_like(f.zlim)
        for jj in range(f.nlim):  #j=1:f.nlim
            for ii in range(f.limitr[jj]):  #i=1:f.limitr(j)
                data = fscanf(fid, fmt2, 2)
                f.rlim[ii, jj] = data[0]
                f.zlim[ii, jj] = data[1]
            # end
        # end
        f.nrgrid = fscanf(fid, fmtg, 1)
        f.nzgrid = fscanf(fid, fmtg, 1)
        f.tokid = fscanf(fid, fmtg, 1)
        f.rx1 = fscanf(fid, fmtg, 1)
        f.rx2 = fscanf(fid, fmtg, 1)
        f.zy1 = fscanf(fid, fmtg, 1)
        f.zy2 = fscanf(fid, fmtg, 1)
        f.conif = fscanf(fid, fmtg, 1)
        f.imatch_phiedge = fscanf(fid, fmtg, 1)
    # end
    return f
# end  def read_vmec_650()      # checked against matlabVMEC on July15,2019

# ========================================================================= #


def read_vmec_695(fid, fmt):
    #function f=read_vmec_695(fid,fmt)
    f = Struct()

    # =============================================== #

    fmtg, fmtd, fmts, fmt2, fmt3, fmt6, fmt7, fmt10, fmt11, fmt12, \
        fmt13, fmt14, fmt20, fmt2_3, fmt2_6, fmt2_7, fmt2_11, fmt2_14, \
        fmt2_3_2_7, fmt2_6_2_14 = return_fmts(fmt)

    # =============================================== #

    f.lfreeb = 0

    data = fscanf(fid, fmt, 7)
    f.wb = data[0]
    f.wp = data[1]
    f.gamma = data[2]
    f.pfac = data[3]
    f.rmax_surf = data[4]
    f.rmin_surf = data[5]
    f.zmax_surf = data[6]

    data = fscanf(fid, fmt, 10)   # TODO:! should be fmtd
    f.nfp = data[0]
    f.ns = data[1]
    f.mpol = data[2]
    f.ntor = data[3]
    f.mnmax = data[4]
    f.itfsq = data[5]
    f.niter = data[6]
    f.iasym = data[7]
    f.ireconstruct = data[8]
    f.ierr_vmec = data[9]

    data = fscanf(fid, fmt, 6)    # TODO:! should be fmtd
    f.imse = data[0]
    f.itse = data[1]
    f.nbsets = data[2]
    f.nobd = data[3]
    f.nextcur = data[4]
    f.nstore_seq = data[5]

    #Error Check
    if (f.ierr_vmec and (f.ierr_vmec != 4)):
        print('ierr_vmec:'+ str(f.ierr_vmec))
        return
    # end

    #Read nbfld
    if (f.nbsets > 0):
        f.nbfld = fscanf(fid, fmt, f.nbsets)
    # end

    # Commented out in matlabVMEC
    if fmt.find(',')<0:
        temp = fgetl(fid)
        temp = fgetl(fid)
        f.mgrid_file = temp.strip()
    else:
        f.mgrid_file = fscanf(fid, fmts, 1)
    # end if
    # Read Arrays
    if f.iasym > 0:
        data1 = fscanf(fid, fmt2_14, [16, f.mnmax])
        data = fscanf(fid, fmt14, [14, f.mnmax*(f.ns-1)])
    else:
        data1 = fscanf(fid, fmt2_11, [13, f.mnmax])
        data = fscanf(fid, fmt11, [11, f.mnmax*(f.ns-1)])
    # end

    # Extract Data from Arrays
    f.xm = _np.asarray(data1[0][:])
    f.xn = _np.asarray(data1[1][:])

    # Extract the data and reshape
#    data1 = _np.asarray(data1[2:][:])
    if f.iasym > 0:
        # On half grid
        f.rmnc, f.zmns, f.lmns, f.bmnc, f.gmnc, f.bsubumnc, f.bsubvmnc, f.bsubsmns, \
            f.bsupumnc, f.bsupvmnc, f.currvmnc, f.rmns, f.zmnc, f.lmnc \
            = tuple([_np.vstack((_np.copy(data1[ii+2][:].T), _np.reshape(data[ii][:], f.mnmax, f.ns-1, order='F').copy()))
                for ii in range(10)])
#            = tuple([_np.vstack((_np.copy(data1[ii][:].T), _np.reshape(data[ii][:], f.mnmax, f.ns-1, order='F').copy()))
#                for ii in range(10)])
    else:
        f.rmnc, f.zmns, f.lmns, f.bmnc, f.gmnc, f.bsubumnc, f.bsubvmnc, f.bsubsmns, \
            f.bsupumnc, f.bsupvmnc, f.currvmnc \
            = tuple([_np.vstack((_np.copy(data1[ii+2][:].T), _np.reshape(data[ii][:], f.mnmax, f.ns-1, order='F').copy()))
                for ii in range(10)])
#            = tuple([_np.vstack((_np.copy(data1[ii][:].T), _np.reshape(data[ii][:], f.mnmax, f.ns-1, order='F').copy()))
#                for ii in range(10)])
    # end

    # Read the half-mesh quantities
    data = fscanf(fid, fmt13, [13, f.ns-1])
    f.iotas = data[0][:]
    f.mass = data[1][:]
    f.pres = data[2][:]
    f.beta_vol = data[3][:]
    f.phip = data[4][:]
    f.buco = data[5][:]
    f.bvco = data[6][:]
    f.phi = data[7][:]
    f.vp = data[8][:]
    f.overr = data[9][:]
    f.jcuru = data[10][:]
    f.jcurv = data[11][:]
    f.specw = data[12][:]

    data = fscanf(fid, fmt, 6)
    f.aspect = data[0]
    f.betatot = data[1]
    f.betapol = data[2]
    f.betator = data[3]
    f.betaxis = data[4]
    f.b0 = data[5]

#    f.isigna = fscanf(fid, fmt+'\n', 1)
    f.isigna = fscanf(fid, fmt, 1)
    f.input_extension = fgetl(fid).strip()

    data = fscanf(fid, fmt, 8)
    f.IonLarmor = data[0]
    f.VolAvgB = data[1]
    f.RBtor0 = data[2]
    f.RBtor = data[3]
    f.Itor = data[4]
    f.Aminor = data[5]
    f.Rmajor = data[6]
    f.Volume = data[7]

    # Mercier Criterion
    data = fscanf(fid, fmt6, [6, f.ns-2])
    f.Dmerc = data[0][:]
    f.Dshear = data[1][:]
    f.Dwell = data[2][:]
    f.Dcurr = data[3][:]
    f.Dgeod = data[4][:]
    f.equif = data[5][:]

#    f.curlabel=cell(f.nextcur, 1)
#    f.curlabel = [f.nextcur, 1]
    f.curlabel = []
    delim = None
    if fmt.find(',')>-1:
        delim = ','
    # end if
    if (f.nextcur > 0):
        f.lfreeb = 1
        f.extcur = fscanf(fid, fmt, f.nextcur)
#        fscanf(fid, '\n')
        rem = f.nextcur
        while rem > 0:
            line = fgetl(fid)
            line = line.replace('\'','').strip() # remove the spaces and apostrophes at end-points
            line = line.split(delim)          # split into a list of Coils (Coil_3, Coil_4)
#            fscanf(fid, '\n')            # skip a line seems wrong here
            index = _np.asarray(range(len(line)))
            for ii in index:
                f.curlabel.append(line[ii].strip())
            # end for
            rem -= len(index)
        # end
    # end

    data = fscanf(fid, fmt2, [2, f.nstore_seq])
    f.sqt = data[0][:]
    f.wdot = data[1][:]

    data = fscanf(fid, fmt2, [2, f.nstore_seq])
    f.jdotb = data[0][:]
    f.bdotgradv = data[1][:]

    # No Unit Conversion Necessary
    # Data and MSE Fits
    if (f.ireconstruct > 0):
        if ((f.imse >= 2) or (f.itse > 0)):
            f.twsgt = fscanf(fid, fmt, 1)
            f.msewgt = fscanf(fid, fmt, 1)
            f.isnodes = fscanf(fid, fmt, 1)

            data = fscanf(fid, fmt3, [3, f.isnodes])
            f.sknots = data[0][:]
            f.ystark = data[1][:]
            f.y2stark = data[2][:]
            f.ipnodes = fscanf(fid, fmt, 1)

            data = fscanf(fid, fmt3, [3, f.ipnodes])
            f.pknots = data[0][:]
            f.ythom = data[1][:]
            f.y2thom = data[2][:]

            data = fscanf(fid, fmt7, [7, (2*f.ns)-1])
            f.anglemse = data[0][:]
            f.rmid = data[1][:]
            f.qmid = data[2][:]
            f.shear = data[3][:]
            f.presmid = data[4][:]
            f.alfa = data[5][:]
            f.curmid = data[6][:]

            data = fscanf(fid, fmt3, [3, f.imse])
            f.rstark = data[0][:]
            f.datastark = data[1][:]
            f.qmeas = data[2][:]

            data = fscanf(fid, fmt2, [2, f.itse])
            f.rthom = data[0][:]
            f.datathom = data[1][:]
        # end

        if (f.nobd > 0):
            data = fscanf(fid, fmt3, [3, f.nobd])
            f.dsiext = data[0][:]
            f.plflux = data[1][:]
            f.dsiobt = data[2][:]
            f.flmwgt = fscanf(fid, fmt, 1)
        # end

        nbfldn = _np.sum(f.nbfld[:f.nbsets])
        if (nbfldn > 0):
            for nn in range(f.nbsets):  # for n=1:nbsets
                data = fscanf(fid, fmt3, [3, f.nbfld[nn]])
                f.bcoil[:, nn] = data[0][:]
                f.plbfld[:, nn] = data[1][:]
                f.bbc[:, nn] = data[2][:]
            # end
            f.bcwgt = fscanf(fid, fmt, 1)
        # end
        f.phidiam = fscanf(fid, fmt, 1)
        f.delphid = fscanf(fid, fmt, 1)

        # Read Limiter and Prout Plotting Specs
        f.nsets = fscanf(fid, fmt, 1)
        f.nparts = fscanf(fid, fmt, 1)
        f.nlim = fscanf(fid, fmt, 1)
        f.nsetsn = fscanf(fid, fmt, f.nsets)
        f.pfcspec=_np.zeros( (f.nparts, _np.max(f.nxsetsn), f.nsets), dtype=_np.float64)
        for kk in range(f.nsets):  # k=1:f.nsets
            for jj in range(f.nsetsn[kk]):  # j=1:f.nsetsn(k)
                for ii in range(f.nparts):  # i=1:f.nparts
                    f.pfcspec[ii, jj, kk] = fscanf(fid, fmt, 1)
                # end
            # end
        # end

        f.limitr = fscanf(fid, fmt, f.nlim)
        f.rlim = _np.zeros( (_np.max(f.limitr), f.nlim), dtype=_np.float64)
        f.zlim = _np.zeros_like(f.zlim)
        for jj in range(f.nlim):  # j=1:f.nlim
            for ii in range(f.limitr[jj]):  # i=1:f.limitr(j)
                data = fscanf(fid, fmt2, 2)
                f.rlim[ii, jj] = data[0]
                f.zlim[ii, jj] = data[1]
            # end
        # end

        f.nrgrid = fscanf(fid, fmt, 1)
        f.nzgrid = fscanf(fid, fmt, 1)
        f.tokid = fscanf(fid, fmt, 1)
        f.rx1 = fscanf(fid, fmt, 1)
        f.rx2 = fscanf(fid, fmt, 1)
        f.zy1 = fscanf(fid, fmt, 1)
        f.zy2 = fscanf(fid, fmt, 1)
        f.conif = fscanf(fid, fmt, 1)
        f.imatch_phiedge = fscanf(fid, fmt, 1)
    # end
    return f
# end  def read_vmec_695

# ========================================================================= #


def read_vmec_800(fid, fmt):
    # function f=read_vmec_800(fid,fmt)
    f = Struct()
    f.lfreeb = 0

    fmtg, fmtd, fmts, fmt2, fmt3, fmt6, fmt7, fmt10, fmt11, fmt12, \
        fmt13, fmt14, fmt20, fmt2_3, fmt2_6, fmt2_7, fmt2_11, fmt2_14, \
        fmt2_3_2_7, fmt2_6_2_14 = return_fmts(fmt)

    # =============================================== #

    data = fscanf(fid, fmtg, 7)
    f.wb = data[0]
    f.wp = data[1]
    f.gamma = data[2]
    f.pfac = data[3]
    f.rmax_surf = data[4]
    f.rmin_surf = data[5]
    f.zmax_surf = data[6]

    data = fscanf(fid, fmtd, 10)
    f.nfp = data[0]
    f.ns = data[1]
    f.mpol = data[2]
    f.ntor = data[3]
    f.mnmax = data[4]
    f.itfsq = data[5]
    f.niter = data[6]
    f.iasym = data[7]
    f.ireconstruct = data[8]
    f.ierr_vmec = data[9]
    f.mnmax_nyq = f.mnmax

    data = fscanf(fid, fmtd, 6)
    f.imse = data[0]
    f.itse = data[1]
    f.nbsets = data[2]
    f.nobd = data[3]
    f.nextcur = data[4]
    f.nstore_seq = data[5]

    # Error Check
    if (f.ierr_vmec and (f.ierr_vmec != 4)):
        print('ierr_vmec: '+str(f.ierr_vmec))
        return f
    # end

    # Read nbfld
    if (f.nbsets > 0):
        f.nbfld = fscanf(fid, fmt, f.nbsets)
    # end

    # Read mgrid filename and setup other format statements
    f.mgrid_file = fscanf(fid, fmts, 1)

    # Read Arrays
    if f.iasym > 0:
        data1 = fscanf(fid, fmt2_14, [16, f.mnmax])
        data = fscanf(fid, fmt14, [14, f.mnmax*(f.ns-1)])
    else:
        data1 = fscanf(fid, fmt2_11, [13, f.mnmax])
        data = fscanf(fid, fmt11, [11, f.mnmax*(f.ns-1)])
    # end

    # Extract Data from Arrays
    f.xm = _np.asarray(data1[0][:])
    f.xn = _np.asarray(data1[1][:])

    # Extract the data and reshape
#    data1 = _np.asarray(data1[2:][:])
    if f.iasym > 0:
        # On half grid
        f.rmnc, f.zmns, f.lmns, f.bmnc, f.gmnc, f.bsubumnc, f.bsubvmnc, f.bsubsmns, \
            f.bsupumnc, f.bsupvmnc, f.currvmnc, f.rmns, f.zmnc, f.lmnc \
            = tuple([_np.vstack((_np.copy(data1[ii+2][:].T), _np.reshape(data[ii][:], f.mnmax, f.ns-1, order='F').copy()))
                for ii in range(10)])
#            = tuple([_np.vstack((_np.copy(data1[ii][:].T), _np.reshape(data[ii][:], f.mnmax, f.ns-1, order='F').copy()))
#                for ii in range(10)])
    else:
        f.rmnc, f.zmns, f.lmns, f.bmnc, f.gmnc, f.bsubumnc, f.bsubvmnc, f.bsubsmns, \
            f.bsupumnc, f.bsupvmnc, f.currvmnc \
            = tuple([_np.vstack((_np.copy(data1[ii+2][:].T), _np.reshape(data[ii][:], f.mnmax, f.ns-1, order='F').copy()))
                for ii in range(10)])
#            = tuple([_np.vstack((_np.copy(data1[ii][:].T), _np.reshape(data[ii][:], f.mnmax, f.ns-1, order='F').copy()))
#                for ii in range(10)])
    # end

    # Read the full-mesh quantities
    data = fscanf(fid, fmt6, [6, f.ns])
    f.iotaf = data[0][:]
    f.presf = data[1][:]
    f.phipf = data[2][:]
    f.phi = data[3][:]
    f.jcuru = data[4][:]
    f.jcurv = data[5][:]

    # Read the half-mesh quantities
    data = fscanf(fid, fmt10, [10, f.ns-1])
    f.iotas = data[0][:]
    f.mass = data[1][:]
    f.pres = data[2][:]
    f.beta_vol = data[3][:]
    f.phip = data[4][:]
    f.buco = data[5][:]
    f.bvco = data[6][:]
    f.vp = data[7][:]
    f.overr = data[8][:]
    f.specw = data[9][:]

    data = fscanf(fid, fmt6, 6)
    f.aspect = data[0]
    f.betatot = data[1]
    f.betapol = data[2]
    f.betator = data[3]
    f.betaxis = data[4]
    f.b0 = data[5]

#    f.isigna = fscanf(fid, fmt+'\n', 1)  # does this mean read to format matching (fmt\n)?
#    f.input_extension = strtrim(fgetl(fid))
    f.isigna = fscanf(fid, fmt, 1)
    f.input_extension = fgetl(fid).strip()

    data = fscanf(fid, fmtg, 8)
    f.IonLarmor = data[0]
    f.VolAvgB = data[1]
    f.RBtor0 = data[2]
    f.RBtor = data[3]
    f.Itor = data[4]
    f.Aminor = data[5]
    f.Rmajor = data[6]
    f.Volume = data[7]

    # Mercier Criterion
    data = fscanf(fid, fmt6, [6, f.ns-2])
    f.Dmerc = data[0][:]
    f.Dshear = data[1][:]
    f.Dwell = data[2][:]
    f.Dcurr = data[3][:]
    f.Dgeod = data[4][:]
    f.equif = data[5][:]

#    f.curlabel = cell(f.nextcur, 1)
    f.curlabel = []
    delim = None
    if fmt.find(',')>-1:
        delim = ','
    # end if
    if (f.nextcur > 0):
        f.lfreeb = 1
        f.extcur = fscanf(fid, fmt, f.nextcur)
#        fscanf(fid, '\n')   # TODO!: need to test this
        rem = f.nextcur
        while rem > 0:
            line = fgetl(fid)
            line = line.replace('\'','').strip() # remove the spaces and apostrophes at end-points
            line = line.split(delim)          # split into a list of Coils (Coil_3, Coil_4)
#            fscanf(fid, '\n')            # skip a line seems wrong here
            index = _np.asarray(range(len(line)))
            for ii in index:
                f.curlabel.append(line[ii].strip())
            # end for
            rem -= len(index)
        # end
    # end

    data = fscanf(fid, fmt2, [2, f.nstore_seq])
    f.sqt = data[0][:]
    f.wdot = data[1][:]

    data = fscanf(fid, fmt2, [2, f.nstore_seq])
    f.jdotb = data[0][:]
    f.bdotgradv = data[1][:]

    # No Unit Conversion Necessary
    # Data and MSE Fits
    if (f.ireconstruct > 0):
        if ((f.imse >= 2) or (f.itse > 0)):
            f.twsgt = fscanf(fid, fmt, 1)
            f.msewgt = fscanf(fid, fmt, 1)
            f.isnodes = fscanf(fid, fmt, 1)

            data = fscanf(fid, fmt3, [3, f.isnodes])
            f.sknots = data[0][:]
            f.ystark = data[1][:]
            f.y2stark = data[2][:]
            f.ipnodes = fscanf(fid, fmt, 1)

            data = fscanf(fid, fmt3, [3, f.ipnodes])
            f.pknots = data[0][:]
            f.ythom = data[1][:]
            f.y2thom = data[2][:]

            data = fscanf(fid,fmt7, [7, (2*f.ns)-1])
            f.anglemse = data[0][:]
            f.rmid = data[1][:]
            f.qmid = data[2][:]
            f.shear = data[3][:]
            f.presmid = data[4][:]
            f.alfa = data[5][:]
            f.curmid = data[6][:]

            data = fscanf(fid, fmt3, [3, f.imse])
            f.rstark = data[0][:]
            f.datastark = data[1][:]
            f.qmeas = data[2][:]

            data = fscanf(fid, fmt2, [2, f.itse])
            f.rthom = data[0][:]
            f.datathom = data[1][:]
        # end

        if (f.nobd > 0):
            data = fscanf(fid, fmt3, [3,f.nobd])
            f.dsiext = data[0][:]
            f.plflux = data[1][:]
            f.dsiobt = data[2][:]
            f.flmwgt = fscanf(fid, fmt, 1)
        # end

        nbfldn = _np.sum(f.nbfld[:f.nbsets])
        if (nbfldn > 0):
            for nn in range(f.nbsets):  # for n=1:nbsets
                data = fscanf(fid, fmt3, [3, f.nbfld[nn]])
                f.bcoil[:, nn] = data[0][:]
                f.plbfld[:, nn] = data[1][:]
                f.bbc[:, nn] = data[2][:]
            # end
            f.bcwgt = fscanf(fid, fmt, 1)
        # end
        f.phidiam = fscanf(fid, fmt, 1)
        f.delphid = fscanf(fid, fmt, 1)

        # Read Limiter and Prout Plotting Specs
        f.nsets = fscanf(fid, fmt, 1)
        f.nparts = fscanf(fid, fmt, 1)
        f.nlim = fscanf(fid, fmt, 1)
        f.nsetsn = fscanf(fid, fmt, f.nsets)
        f.pfcspec = _np.zeros( (f.nparts, _np.max(f.nxsetsn), f.nsets), dtype=_np.float64)
        for kk in range(f.nsets):  # k=1:f.nsets
            for jj in range(f.nsetns[kk]):  # j=1:f.nsetsn[kk]
                for ii in range(f.nparts):  # i=1:f.nparts
                    f.pfcspec[ii, jj, kk] = fscanf(fid, fmt, 1)
                # end
            # end
        # end

        f.limitr = fscanf(fid, fmt, f.nlim)
        f.rlim = _np.zeros( (_np.max(f.limitr), f.nlim), dtype=_np.float64)
        f.zlim = _np.zeros_like(f.rlim)
        for jj in range(f.nlim):  # j=1:f.nlim
            for ii in range(f.limitr[jj]):  # i=1:f.limitr(j)
                data = fscanf(fid, fmt2, 2)
                f.rlim[ii, jj] = data[0]
                f.zlim[ii, jj] = data[1]
            # end
        # end

        f.nrgrid = fscanf(fid, fmt, 1)
        f.nzgrid = fscanf(fid, fmt, 1)
        f.tokid = fscanf(fid, fmt, 1)
        f.rx1 = fscanf(fid, fmt, 1)
        f.rx2 = fscanf(fid, fmt, 1)
        f.zy1 = fscanf(fid, fmt, 1)
        f.zy2 = fscanf(fid, fmt, 1)
        f.conif = fscanf(fid, fmt, 1)
        f.imatch_phiedge = fscanf(fid, fmt, 1)
    # end

    return f
# end  def read_vmec_800()   # checked against matlabVMEC on July15,2019

# ========================================================================= #


def read_vmec_847(fid,fmt):
    f = Struct()

    # Unit Conversions
    mu0 = 4.0*_np.pi*1e-7
#    dmu0=(2.0e-7)/_np.pi

    # =============================================== #

    fmtg, fmtd, fmts, fmt2, fmt3, fmt6, fmt7, fmt10, fmt11, fmt12, \
        fmt13, fmt14, fmt20, fmt2_3, fmt2_6, fmt2_7, fmt2_11, fmt2_14, \
        fmt2_3_2_7, fmt2_6_2_14 = return_fmts(fmt)

    # =============================================== #

    f.lfreeb = 0

    data = fscanf(fid, fmtg, 7)
    f.wb = data[0]
    f.wp = data[1]
    f.gamma = data[2]
    f.pfac = data[3]
    f.rmax_surf = data[4]
    f.rmin_surf = data[5]
    f.zmax_surf = data[6]

#    data = fscanf(fid, fmt, 11)
    data = fscanf(fid, fmtd, 11)
    f.nfp = data[0]
    f.ns = data[1]
    f.mpol = data[2]
    f.ntor = data[3]
    f.mnmax = data[4]
    f.mnmax_nyq = data[5]
    f.itfsq = data[6]
    f.niter = data[7]
    f.iasym = data[8]
    f.ireconstruct = data[9]
    f.ierr_vmec = data[10]

#    data = fscanf(fid, fmt, 6)
    data = fscanf(fid, fmtd, 6)
    f.imse = data[0]
    f.itse = data[1]
    f.nbsets = data[2]
    f.nobd = data[3]
    f.nextcur = data[4]
    f.nstore_seq = data[5]

    # Error Check
    if f.ierr_vmec and (f.ierr_vmec != 4):
        print('ierr_vmec: '+str(f.ierr_vmec))
        return None
    # end

    # Read nbfld
    if (f.nbsets > 0):
        f.nbfld = fscanf(fid, fmt, f.nbsets)
    # end

    # Read mgrid filename
    f.mgrid_file = fscanf(fid, fmts, 1)

    # Read Arrays
#    f.xm = _np.zeros( (1,f.mnmax), dtype=_np.int64)
#    f.xn = _np.zeros_like(f.xm)
    f.xm = _np.zeros( (f.mnmax,), dtype=_np.int64)
    f.xn = _np.zeros_like(f.xm)

    f.rmnc = _np.zeros( (f.mnmax, f.ns), dtype=_np.float64)
    f.zmns = _np.zeros_like(f.rmnc)
    f.lmns = _np.zeros_like(f.rmnc)

#    f.xm_nyq = _np.zeros( (1, f.mnmax_nyq), dtype=_np.int64)
    f.xm_nyq = _np.zeros( (f.mnmax_nyq,), dtype=_np.int64)
    f.xn_nyq = _np.zeros_like(f.xm_nyq)

    f.bmnc = _np.zeros( (f.mnmax_nyq, f.ns), dtype=_np.float64)
    f.gmnc = _np.zeros_like(f.bmnc)
    f.bsubumnc = _np.zeros_like(f.bmnc)
    f.bsubvmnc = _np.zeros_like(f.bmnc)
    f.bsubsmns = _np.zeros_like(f.bmnc)
    f.bsupumnc = _np.zeros_like(f.bmnc)
    f.bsupvmnc = _np.zeros_like(f.bmnc)

    if f.iasym >0:
        f.rmns = _np.zeros( (f.mnmax, f.ns), dtype=_np.float64)
        f.zmnc = _np.zeros_like(f.rmns)
        f.lmnc = _np.zeros_like(f.rmns)

        f.bmns = _np.zeros( (f.mnmax_nyq, f.ns), dtype=_np.float64)
        f.gmns = _np.zeros_like(f.bmns)
        f.bsubumns = _np.zeros_like(f.bmns)
        f.bsubvmns = _np.zeros_like(f.bmns)
        f.bsubsmnc = _np.zeros_like(f.bmns)
        f.bsupumns = _np.zeros_like(f.bmns)
        f.bsupvmns = _np.zeros_like(f.bmns)
    # end

    for ii in range(f.ns): # i=1:f.ns
        for jj in range(f.mnmax): # j=1:f.mnmax
            if ii==0:
                f.xm[jj] = fscanf(fid, fmtd, 1)
                f.xn[jj] = fscanf(fid, fmtd, 1)
            # end
            f.rmnc[jj, ii] = fscanf(fid, fmtg, 1)
            f.zmns[jj, ii] = fscanf(fid, fmtg, 1)
            f.lmns[jj, ii] = fscanf(fid, fmtg, 1)

            if f.iasym > 0:
                f.rmns[jj, ii] = fscanf(fid, fmtg, 1)
                f.zmnc[jj, ii] = fscanf(fid, fmtg, 1)
                f.lmnc[jj, ii] = fscanf(fid, fmtg, 1)
            # end
        # end

        for jj in range(f.mnmax_nyq): # j=1:f.mnmax_nyq
            if ii==0:
                f.xm_nyq[jj] = fscanf(fid, fmtd, 1)
                f.xn_nyq[jj] = fscanf(fid, fmtd, 1)
            # end

            f.bmnc[jj, ii] = fscanf(fid, fmtg, 1)
            f.gmnc[jj, ii] = fscanf(fid, fmtg, 1)
            f.bsubumnc[jj, ii] = fscanf(fid, fmtg, 1)
            f.bsubvmnc[jj, ii] = fscanf(fid, fmtg, 1)
            f.bsubsmns[jj, ii] = fscanf(fid, fmtg, 1)
            f.bsupumnc[jj, ii] = fscanf(fid, fmtg, 1)
            f.bsupvmnc[jj, ii] = fscanf(fid, fmtg, 1)

            if f.iasym > 0:
                f.bmns[jj, ii] = fscanf(fid, fmtg, 1)
                f.gmns[jj, ii] = fscanf(fid, fmtg, 1)
                f.bsubumns[jj, ii] = fscanf(fid, fmtg, 1)
                f.bsubvmns[jj, ii] = fscanf(fid, fmtg, 1)
                f.bsubsmnc[jj, ii] = fscanf(fid, fmtg, 1)
                f.bsupumns[jj, ii] = fscanf(fid, fmtg, 1)
                f.bsupvmns[jj, ii] = fscanf(fid, fmtg, 1)
            # end
        # end
    # end

    f.mnyq = _np.max(f.xm_nyq)
    f.nnyq = _np.max(f.xn_nyq)/f.nfp

    # Calculate the Currents
    f.currumnc = _np.zeros( (f.mnmax_nyq, f.ns), dtype=_np.float64)
    f.currvmnc = _np.zeros_like(f.currumnc)
    ohs = f.ns - 1
    hs = 1.0/ohs
    ns = f.ns

    shalf = _np.zeros((ns,), dtype=_np.float64)
    sfull = _np.zeros_like(shalf)
    for ii in range(1,ns):  # 2:ns
        shalf[ii] = _np.sqrt(hs*(ii-1.5))
        sfull[ii] = _np.sqrt(hs*(ii-1.0))
    # end for
    js1 = range(2,ns)    # 3:ns
    js = range(1, ns-1)  # 2:(ns-1)
    for mn in range(f.mnmax_nyq):  # 1:f.mnmax_nyq:
        if f.xm_nyq[mn]%2 == 0:     # mod(f.xm_nyq, 2) == 1
            t1 = 0.5*(shalf[js1]*f.bsubsmns[mn, js1] +
                      shalf[js]*f.bsubsmns[mn, js])/sfull[js]
            bu0 = f.bsubumnc[mn, js]/ shalf[js]
            bu1 = f.bsubumnc[mn, js1]/ shalf[js1]

            t2 = ohs*(bu1-bu0)*sfull[js]+0.25*(bu0+bu1)/sfull[js]
            bv0 = f.bsubvmnc[mn, js] / shalf[js]
            bv1 = f.bsubvmnc[mn, js1] / shalf[js1]

            t3 = ohs*(bv1-bv0)*sfull[js]+0.25*(bv0+bv1)/sfull[js]
        else:
            t1 = 0.5*(f.bsubsmns[mn, js1] + f.bsubsmns[mn, js])
            t2 = ohs*(f.bsubumnc[mn, js1] + f.bsubumnc[mn, js])
            t3 = ohs*(f.bsubvmnc[mn, js1] + f.bsubvmnc[mn, js])
        # end if
        f.currumnc[mn, js] = -_np.float64(f.xn_nyq[mn])*t1-t3
        f.currvmnc[mn, js] = -_np.float64(f.xn_nyq[mn])*t1+t2
    # end for
    # OLD way
    # for ii in range(1,f.ns-1): # i=2:f.ns-1
    #     f.currumnc[:, ii] = -f.xn_nyq.flatten()*f.bsubsmns[:, ii] - (f.bsubvmnc[:, ii+1] - f.bsubvmnc[:, ii-1])/(f.ns-1)
    #     f.currvmnc[:, ii] = -f.xm_nyq.flatten()*f.bsubsmns[:, ii] - (f.bsubumnc[:, ii+1] - f.bsubumnc[:, ii-1])/(f.ns-1)
    # # end

    f.currumnc[:, 0] = 0.0
    f.currvmnc[:, 0] = 0.0
    for ii in range(f.mnmax_nyq): # 1:f.mnmax_nyq
        if f.xm_nyq[ii] == 0:
            f.currumnc[ii, 0] = 2*f.currumnc[ii, 1] - f.currumnc[ii, 2]
            f.currvmnc[ii, 0] = 2*(f.ns-1)*f.bsubumnc[ii, 1]
        # end if
    # end for
    f.currumnc[:, -1] = 2*f.currumnc[:, -2] - f.currumnc[:, -3]
    f.currvmnc[:, -1] = 2*f.currvmnc[:, -2] - f.currvmnc[:, -3]
    f.currumnc /= mu0
    f.currvmnc /= mu0

    if f.iasym>0:
        f.currumns = _np.zeros( (f.mnmax_nyq, f.ns), dtype=_np.float64)
        f.currvmns = _np.zeros_like(f.currumns)
        for mn in range(f.mnmax_nyq):  # 1:f.mnmax_nyq:
#            if f.xm_nyq%2 == 0:     # mod(f.xm_nyq, 2) == 1
            if _np.mod(f.xm_nyq, 2) == 1:     # mod(f.xm_nyq, 2) == 1
                t1 = 0.5*(shalf[js1]*f.bsubsmnc[mn, js1] +
                          shalf[js]*f.bsubsmnc[mn, js])/sfull[js]
                bu0 = f.bsubumns[mn, js]/ shalf[js]
                bu1 = f.bsubumns[mn, js1]/ shalf[js1]

                t2 = ohs*(bu1-bu0)*sfull[js]+0.25*(bu0+bu1)/sfull[js]
                bv0 = f.bsubvmns[mn, js] / shalf[js]
                bv1 = f.bsubvmns[mn, js1] / shalf[js1]

                t3 = ohs*(bv1-bv0)*sfull[js]+0.25*(bv0+bv1)/sfull[js]
            else:
                t1 = 0.5*(f.bsubsmnc[mn, js1] + f.bsubsmnc[mn, js])
                t2 = ohs*(f.bsubumns[mn, js1] + f.bsubumns[mn, js])
                t3 = ohs*(f.bsubvmns[mn, js1] + f.bsubvmns[mn, js])
            # end if
            f.currumns[mn, js] = -_np.float64(f.xn_nyq[mn])*t1-t3
            f.currvmns[mn, js] = -_np.float64(f.xn_nyq[mn])*t1+t2
        # end for
        # OLD way
        # for ii in range(1,f.ns-1): # i=2:f.ns-1
        #     f.currumns[:, ii] = -f.xnnyq.flatten()*f.bsubsmnc[:, ii] - (f.ns-1)*(f.bsubvmns[:, ii+1] - f.bsubvmns[:, ii])
        #     f.currvmns[:, ii] = -f.xmnyq.flatten()*f.bsubsmnc[:, ii] + (f.ns-1)*(f.bsubumns[:, ii+1] - f.bsubumns[:, ii])
        # # end

        f.currumns[:, 0] = 0.0
        f.currvmns[:, 0] = 0.0
        for ii in range(f.mnmax_nyq): # 1:f.mnmax_nyq
            if f.xm_nyq[ii] == 0:
                f.currumns[ii, 0] = 2*f.currumns[ii, 1] - f.currumns[ii, 2]
                f.currvmns[ii, 0] = 2*(f.ns-1)*f.bsubumns[ii, 1]
            # end if
        # end for
        f.currumns[:, -1] = 2*f.currumns[:, -2] - f.currumns[:, -3]
        f.currvmns[:, -1] = 2*f.currvmns[:, -2] - f.currvmns[:, -3]
        f.currumns /= mu0
        f.currvmns /= mu0
    # end for

#    f.currumnc[f.xm_nyq==0, 0] = 2*f.currumnc[ f.xm_nyq == 0, 1] - f.currumnc[ f.xm_nyq == 0, 2]
#    f.currvmnc[f.xm_nyq==0, 0] = 2*f.bsubumnc[ f.xm_nyq == 0, 1]/(f.ns-1)
#
#    f.currvmnc[f.xm_nyq!=0, 0] = 0.0
#    f.currumnc[f.xm_nyq!=0, 0] = 0.0
#    f.currumnc[:, f.ns] = 2*f.currumnc[:, f.ns-1] - f.currumnc[:, f.ns-2]
#    f.currvmnc[:, f.ns] = 2*f.currvmnc[:, f.ns-1] - f.currvmnc[:, f.ns-2]
#    f.currumnc = f.currumnc/(4*_np.pi*1e-7)
#    f.currvmnc = f.currvmnc/(4*_np.pi*1e-7)
#
#    if f.iasym > 0:
#        f.currumns = _np.zeros( (f.mnmax_nyq, f.ns), dtype=_np.float64)
#        f.currvmns = _np.zeros_like(f.currumns)
#        for ii in range(1, f.ns-1):  # i=2:f.ns-1
#            f.currumns[:, ii] = -f.xn_nyq.flatten()*f.bsubsmnc[:, ii]-(f.bsubvmns[:, ii+1]-f.bsubvmns[:, ii-1])/(f.ns-1)
#            f.currvmns[:, ii] = -f.xm_nyq.flatten()*f.bsubsmnc[:, ii]-(f.bsubumns[:, ii+1]-f.bsubumns[:, ii-1])/(f.ns-1)
#        # end
#
#        f.currvmns[f.xm_nyq == 0, 0] = 2*f.bsubumns[f.xm_nyq == 0, 1]/(f.ns-1)
#        f.currumns[f.xm_nyq == 0, 0] = 2*f.currumns[f.xm_nyq == 0, 1]-f.currumns[f.xm_nyq == 0, 2]
#        f.currvmns[f.xm_nyq != 0, 0] = 0.0
#        f.currumns[f.xm_nyq != 0, 0] = 0.0
#        f.currumns[:, f.ns] = 2*f.currumns[:, f.ns-1]-f.currumns[:, f.ns-2]
#        f.currvmns[:, f.ns] = 2*f.currvmns[:, f.ns-1]-f.currvmns[:, f.ns-2]
#        f.currumns = f.currumns/(4*_np.pi*1e-7)
#        f.currvmns = f.currvmns/(4*_np.pi*1e-7)
#    # end

    # Read the full-mesh quantities
    data = fscanf(fid, fmt6, [6, f.ns])
    f.iotaf = data[0][:]
    f.presf = data[1][:]
    f.phipf = data[2][:]
    f.phi = data[3][:]
    f.jcuru = data[4][:]
    f.jcurv = data[5][:]

    # Read the half-mesh quantities
    data = fscanf(fid, fmt10, [10, f.ns-1])
    f.iotas = data[0][:]
    f.mass = data[1][:]
    f.pres = data[2][:]
    f.beta_vol = data[3][:]
    f.phip = data[4][:]
    f.buco = data[5][:]
    f.bvco = data[6][:]
    f.vp = data[7][:]
    f.overr = data[8][:]
    f.specw = data[9][:]

    data = fscanf(fid, fmt, 6)
    f.aspect = data[0]
    f.betatot = data[1]
    f.betapol = data[2]
    f.betator = data[3]
    f.betaxis = data[4]
    f.b0 = data[5]

#    f.isigna = fscanf(fid, fmt+'\n', 1)
#    f.input_extension = strtrim(fgetl(fid))
    f.isigna = fscanf(fid, fmt, 1)
    f.input_extension = fgetl(fid).strip()

    data = fscanf(fid, fmtg, 8)
    f.IonLarmor = data[0]
    f.VolAvgB = data[1]
    f.RBtor0 = data[2]
    f.RBtor = data[3]
    f.Itor = data[4]
    f.Aminor = data[5]
    f.Rmajor = data[6]
    f.Volume = data[7]

    # Mercier Criterion
    data = fscanf(fid, fmt6, [6, f.ns-2])
    f.Dmerc = data[0][:]
    f.Dshear = data[1][:]
    f.Dwell = data[2][:]
    f.Dcurr = data[3][:]
    f.Dgeod = data[4][:]
    f.equif = data[5][:]

#    f.curlabel=cell(f.nextcur, 1)
    f.curlabel = []
    delim = None
    if fmt.find(',')>-1:
        delim = ','
    # end if
    if (f.nextcur > 0):
        f.lfreeb = 1
        f.extcur = fscanf(fid, fmt, f.nextcur)
        lcurr = fscanf(fid, fmts, 1).strip()
        if lcurr.upper() == 'T':  # strcmpi(lcurr,'T'):
#            fscanf(fid, '\n')  # skipping seems wrong here?
            rem = f.nextcur
#            jj = 0
            while rem > 0:
                line = fgetl(fid)
                line = line.replace('\'','').strip() # remove the spaces and apostrophes at end-points
                line = line.split(delim)          # split into a list of Coils (Coil_3, Coil_4)
#                fscanf(fid, '\n')            # skip a line seems wrong here
                index = _np.asarray(range(len(line)))
                for ii in index:
                    f.curlabel.append(line[ii].strip())
                # end for
                rem -= len(index)
            # end
        # end
    # end

    data = fscanf(fid, fmt2, [2, f.nstore_seq])
    f.sqt = data[0][:]
    f.wdot = data[1][:]

    data = fscanf(fid, fmt2, [2, f.nstore_seq])   # running into eof here
    f.jdotb = data[0][:]
    f.bdotgradv = data[1][:]

    # No Unit Conversion Necessary
    # Data and MSE Fits
    if (f.ireconstruct > 0):
        if ((f.imse >= 2) or (f.itse >0)):
            f.twsgt = fscanf(fid, fmt, 1)
            f.msewgt = fscanf(fid, fmt, 1)
            f.isnodes = fscanf(fid, fmt, 1)

            data = fscanf(fid, fmt3, [3, f.isnodes])
            f.sknots = data[0][:]
            f.ystark = data[1][:]
            f.y2stark = data[2][:]

            f.ipnodes = fscanf(fid, fmt, 1)

            data = fscanf(fid, fmt3, [3, f.ipnodes])
            f.pknots = data[0][:]
            f.ythom = data[1][:]
            f.y2thom = data[2][:]

            data = fscanf(fid,fmt7, [7, (2*f.ns)-1])
            f.anglemse = data[0][:]
            f.rmid = data[1][:]
            f.qmid = data[2][:]
            f.shear = data[3][:]
            f.presmid = data[4][:]
            f.alfa = data[5][:]
            f.curmid = data[6][:]

            data = fscanf(fid, fmt3, [3, f.imse])
            f.rstark = data[0][:]
            f.datastark = data[1][:]
            f.qmeas = data[2][:]

            data = fscanf(fid, fmt2, [2, f.itse])
            f.rthom = data[0][:]
            f.datathom = data[1][:]
        # end

        if (f.nobd > 0):
            data = fscanf(fid, fmt3, [3,f.nobd])
            f.dsiext = data[0][:]
            f.plflux = data[1][:]
            f.dsiobt = data[2][:]

            f.flmwgt = fscanf(fid, fmt, 1)
        # end

        nbfldn=_np.sum(f.nbfld[:f.nbsets])
        if (nbfldn > 0):
            for nn in range(f.nbsets):  # n=1:nbsets
                data = fscanf(fid, fmt3, [3, f.nbfld[nn]])
                f.bcoil[:, nn] = data[0][:]
                f.plbfld[:, nn] = data[1][:]
                f.bbc[:, nn] = data[2][:]
            # end
            f.bcwgt = fscanf(fid, fmt, 1)
        # end
        f.phidiam = fscanf(fid, fmt, 1)
        f.delphid = fscanf(fid, fmt, 1)

        # Read Limiter and Prout Plotting Specs
        f.nsets = fscanf(fid, fmt, 1)
        f.nparts = fscanf(fid, fmt, 1)
        f.nlim = fscanf(fid, fmt, 1)
        f.nsetsn = fscanf(fid, fmt, f.nsets)
        f.pfcspec = _np.zeros( (f.nparts, _np.max(f.nxsetsn), f.nsets), dtype=_np.float64)
        for kk in range( f.nsets ):  # k=1:f.nsets
            for jj in range( f.nsetsn[kk]):  # j=1:f.nsetsn(k)
                for ii in range( f.nparts):  # i=1:f.nparts
                    f.pfcspec[ii, jj, kk] = fscanf(fid, fmt, 1)
                # end
            # end
        # end

        f.limitr = fscanf(fid, fmt, f.nlim)
        f.rlim = _np.zeros( (_np.max(f.limitr), f.nlim), dtype=_np.float64)
        f.zlim = _np.zeros_like(f.zlim)
        for jj in range(f.nlim):  # j=1:f.nlim
            for ii in range(f.limitr[jj]):  # ii=1:f.limitr(jj)
                data = fscanf(fid, fmt2, 2)
                f.rlim[ii, jj] = data[0]
                f.zlim[ii, jj] = data[1]
            # end
        # end

        f.nrgrid = fscanf(fid, fmt, 1)
        f.nzgrid = fscanf(fid, fmt, 1)
        f.tokid = fscanf(fid, fmt, 1)
        f.rx1 = fscanf(fid, fmt, 1)
        f.rx2 = fscanf(fid, fmt, 1)
        f.zy1 = fscanf(fid, fmt, 1)
        f.zy2 = fscanf(fid, fmt, 1)
        f.conif = fscanf(fid, fmt, 1)
        f.imatch_phiedge = fscanf(fid, fmt, 1)
    # end
    f.mgrid_mode=fscanf(fid, fmts).strip()    # if we already hit the end of file, then this will return None?
    return f
# end  def read_vmec_847    checked against matlabVMEC on July15,2019

# ========================================================================= #
# ========================================================================= #








