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

    03/03/2016     G.M. Weir (IPP-Greifswald) Ported to Python (old version)
        To do: Still needs updates!
    07/12/2019     G.M. Weir updated for general use
"""
# ======================================================================== #
# ======================================================================== #

from __future__ import absolute_import, with_statement, absolute_import, \
                       division, print_function, unicode_literals

import numpy as _np
from utils import Struct, Spline, fscanf, fgetl

# ======================================================================== #
# ======================================================================== #


def read_vmec_orig(fid, fmt):
    # For 5.10 and earlier
    f = Struct()

    # Unit Conversions
    dmu0=(2.0e-7)/_np.pi

    # Read Data
    data = fscanf(fid,'%g', 13)
    f.wb    = data[0]
    f.wp    = data[1]
    f.gamma = data[2]
    f.pfac  = data[3]
    f.nfp   = data[4]
    f.ns    = data[5]
    f.mpol  = data[6]
    f.ntor  = data[7]
    f.mnmax = data[8]
    f.itfsq = data[9]
    f.niter = data[10]
    f.iasym = data[11]
    f.ireconstruct = data[12]

    data = fscanf(fid,'%d', 5)
    f.imse    = data[0]
    f.itse    = data[1]
    f.nbsets  = data[2]
    f.nobd    = data[3]
    f.nextcur = data[4]

    f.nstore_seq = 100
    # Error Check
    if f.ierr_vmec and (f.ierr_vmec != 4):
        print('ierr_vmec >0')
        return None
    # end if

    # Read nbfld
    if (f.nbsets > 0):
        f.nbfld = fscanf(fid,'%g',f.nbsets)
    # end if

    # Read mgrid filename
    if fmt.find(',')<0:  # isempty(strfind(fmt,',')):
        f.mgrid_file = fscanf(fid, '%s', 1)
        fmt2 = '%g%g'
        fmt3 = '%g%g%g'
        fmt6 = '%g%g%g%g%g%g'
        fmt7 = '%g%g%g%g%g%g%g'
        fmt11 = '%g%g%g%g%g%g%g%g%g%g%g'
        fmt12 = '%g%g%g%g%g%g%g%g%g%g%g%g'
        fmt14 = '%g%g%g%g%g%g%g%g%g%g%g%g%g%g'
        fmt2_11 = '%d%d%g%g%g%g%g%g%g%g%g%g%g'
        fmt2_14 = '%d%d%g%g%g%g%g%g%g%g%g%g%g%g%g%g'
    else:
        f.mgrid_file = fscanf(fid, '%s,', 1)
        fmt2 = '%g,%g'
        fmt3 = '%g,%g,%g'
        fmt6 = '%g,%g,%g,%g,%g,%g'
        fmt7 = '%g,%g,%g,%g,%g,%g,%g'
        fmt11 = '%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g'
        fmt12 = '%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g'
        fmt14 = '%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g'
        fmt2_11 = '%d,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g'
        fmt2_14 = '%d,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g'
    # end

    #Read Arrays
    if f.iasym > 0:
        data = fscanf(fid, fmt2_14, [16, f.mnmax])
        data = fscanf(fid, fmt14, [14, f.mnmax*(f.ns-1)])
    else:
        data = fscanf(fid,fmt2_11,[13, f.mnmax])
        data = fscanf(fid,fmt11,[11, f.mnmax*(f.ns-1)])
    # end

    #Extract Data from Arrays
    f.xm = data1[0, :]
    f.xn = data1[1, :]

    #First reshape data
    f.rmnc = [data1[2, :].T, reshape(data[0, :], f.mnmax, f.ns-1)]
    f.zmns = [data1[3, :].T, reshape(data[1, :], f.mnmax, f.ns-1)]
    f.lmns = [data1[4, :].T, reshape(data[2, :], f.mnmax, f.ns-1)]
    f.bmn = [data1[5, :].T, reshape(data[3, :], f.mnmax, f.ns-1)]
    f.gmn = [data1[6, :].T, reshape(data[4, :], f.mnmax, f.ns-1)]
    f.bsubumn = [data1[7, :].T, reshape(data[5, :], f.mnmax, f.ns-1)]
    f.bsubvmn = [data1[8, :].T, reshape(data[6, :], f.mnmax, f.ns-1)]
    f.bsubsmn = [data1[9, :].T, reshape(data[7, :], f.mnmax, f.ns-1)]
    f.bsupumn = [data1[10, :].T, reshape(data[8, :], f.mnmax, f.ns-1)]
    f.bsupvmn = [data1[11, :].T, reshape(data[9, :], f.mnmax, f.ns-1)]
    f.currvmn = [data1[12, :].T, reshape(data[10, :], f.mnmax, f.ns-1)]

    #Read the half-mesh quantities
    data = fscanf(fid, fmt12, [12, f.ns/2])
    f.iotas = data[0, :]
    f.mass  = data[1, :]
    f.pres  = data[2, :]
    f.phip  = data[3, :]
    f.buco  = data[4, :]
    f.bvco  = data[5, :]
    f.phi   = data[6, :]
    f.vp    = data[7, :]
    f.overr = data[8, :]
    f.jcuru = data[9, :]
    f.jcurv = data[10, :]
    f.specw = data[11, :]

    data = fscanf(fid, fmt6, 6)
    f.aspect  = data[0]
    f.betatot = data[1]
    f.betapol = data[2]
    f.betator = data[3]
    f.betaxis = data[4]
    f.b0      = data[5]

    #Mercier Criterion
    data = fscanf(fid, fmt6, [6, f.ns-2])
    f.Dmerc  = data[0, :]
    f.Dshear = data[1, :]
    f.Dwell  = data[2, :]
    f.Dcurr  = data[3, :]
    f.Dgeod  = data[4, :]
    f.equif  = data[5, :]
    if (f.nextcur > 0):
        f.extcur = fscanf(fid, fmt, f.nextcur)
        f.curlabel = fscanf(fid, fmt, f.nextcur)
    #end if

    data = fscanf(fid, fmt2, [2, f.nstore_seq])
    f.sqt  = data[0,:]
    f.wdot = data[1,:]

    #Convert from Internal Units to Physical Units
    f.mass  = f.mass/dmu0
    f.pres  = f.pres/dmu0
    f.jcuru = f.jcuru/dmu0
    f.jcurv = f.jcurv/dmu0
    f.jdotb = f.jdotb/dmu0
    f.phi   = -f.phi       # Data and MSE Fits

    if (f.ireconstruct > 0):
        if (f.imse >= 2) or (f.itse >0):
            f.twsgt   = fscanf(fid, fmt, 1)
            f.msewgt  = fscanf(fid, fmt, 1)
            f.isnodes = fscanf(fid, fmt, 1)

            data = fscanf(fid, fmt3, [3, f.isnodes])
            f.sknots  = data[0, :]
            f.ystark  = data[1, :]
            f.y2stark = data[2, :]
            f.ipnodes = fscanf(fid, fmt, 1)

            data = fscanf(fid, fmt3, [3, f.ipnodes])
            f.pknots = data[0, :]
            f.ythom  = data[1, :]
            f.y2thom = data[2, :]

            data = fscanf(fid, fmt7, [7, (2*f.ns)-1])
            f.anglemse = data[0, :]
            f.rmid     = data[1, :]
            f.qmid     = data[2, :]
            f.shear    = data[3, :]
            f.presmid  = data[4, :]
            f.alfa     = data[5, :]
            f.curmid   = data[6, :]

            data = fscanf(fid, fmt3, [3, f.imse])
            f.rstark    = data[0, :]
            f.datastark = data[1, :]
            f.qmeas     = data[2, :]

            data = fscanf(fid, fmt2, [2, f.itse])
            f.rthom    = data[0, :]
            f.datathom = data[1, :]
        # end

        if (f.nobd > 0):
            data = fscanf(fid, fmt3, [3, f.nobd])
            f.dsiext = data[0, :]
            f.plflux = data[1, :]
            f.dsiobt = data[2, :]
            f.flmwgt = fscanf(fid, fmt, 1)
        # end

        nbfldn = _np.sum( nbldf[:f.nbsets] )
        if (nbfldn > 0):
            for nn in range(nbsets): # n=1:nbsets
                data = fscanf(fid, fmt3, [3, f.nbfld[nn]])
                f.bcoil[:, nn]  = data[0, :]
                f.plbfld[:, nn] = data[1, :]
                f.bbc[:, nn]    = data[2, :]
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
# end  def read_vmec_orig

#---------------------------------------------------------------- %

def read_vmec_605(fid,fmt):
    #For 6.05
    f = Struct()

    #Unit Conversions
    dmu0=(2.0e-7)/_np.pi

    #Read Data
    data = fscanf(fid,'%g',6)
    f.wb        = data[0]
    f.wp        = data[1]
    f.gamma     = data[2]
    f.pfac      = data[3]
    f.rmax_surf = data[4]
    f.rmin_surf = data[5]

    data = fscanf(fid,'%d',10)
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

    data = fscanf(fid,'%d',5)
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
        f.nbfld = fscanf(fid, '%g', f.nbsets)
    # end

    #Read mgrid filename
    f.mgrid_file = fscanf(fid, '%s', 1)

    #Read Arrays
    if f.iasym > 0:
        data = fscanf(fid,'%d%d%g%g%g%g%g%g%g%g%g%g%g%g%g%g',[16, f.mnmax])
        data  = fscanf(fid,'%g%g%g%g%g%g%g%g%g%g%g%g%g%g',[14, f.mnmax*(f.ns-1)])
    else:
        data = fscanf(fid,'%d%d%g%g%g%g%g%g%g%g%g%g%g',[13, f.mnmax])
        data  = fscanf(fid,'%g%g%g%g%g%g%g%g%g%g%g',[11, f.mnmax*(f.ns-1)])
    # end

    #Extract Data from Arrays
    f.xm = data1[0, :]
    f.xn = data1[1, :]

    #First reshape data
    f.rmnc = [data1[2, :].T, reshape(data[0, :], f.mnmax, f.ns-1)]
    f.zmns = [data1[3, :].T, reshape(data[1, :], f.mnmax, f.ns-1)]
    f.lmns = [data1[4, :].T, reshape(data[2, :], f.mnmax, f.ns-1)]
    f.bmnc = [data1[5, :].T, reshape(data[3, :], f.mnmax, f.ns-1)]
    f.gmnc = [data1[6, :].T, reshape(data[4, :], f.mnmax, f.ns-1)]
    f.bsubumnc = [data1[7, :].T, reshape(data[5, :], f.mnmax, f.ns-1)]
    f.bsubvmnc = [data1[8, :].T, reshape(data[6, :], f.mnmax, f.ns-1)]
    f.bsubsmns = [data1[9, :].T, reshape(data[7, :], f.mnmax, f.ns-1)]
    f.bsupumnc = [data1[10, :].T, reshape(data[8, :], f.mnmax, f.ns-1)]
    f.bsupvmnc = [data1[11, :].T, reshape(data[9, :], f.mnmax, f.ns-1)]
    f.currvmnc = [data1[12, :].T, reshape(data[10, :], f.mnmax, f.ns-1)]

    #Read the half-mesh quantities
    data = fscanf(fid,'%g%g%g%g%g%g%g%g%g%g%g%g',[12, f.ns/2])
    f.iotas = data[0, :]
    f.mass  = data[1, :]
    f.pres  = data[1, :]
    f.phip  = data[2, :]
    f.buco  = data[3, :]
    f.bvco  = data[4, :]
    f.phi   = data[5, :]
    f.vp    = data[6, :]
    f.overr = data[7, :]
    f.jcuru = data[8, :]
    f.jcurv = data[9, :]
    f.specw = data[10, :]

    data = fscanf(fid,'%g%g%g%g%g%g',6)
    f.aspect  = data[0]
    f.betatot = data[1]
    f.betapol = data[2]
    f.betator = data[3]
    f.betaxis = data[4]
    f.b0      = data[5]

    #Mercier Criterion
    data = fscanf(fid,'%g%g%g%g%g%g',[6, f.ns-2])
    f.Dmerc  = data[0, :]
    f.Dshear = data[1, :]
    f.Dwell  = data[2, :]
    f.Dcurr  = data[3, :]
    f.Dgeod  = data[4, :]
    f.equif  = data[5, :]

    if (f.nextcur > 0):
        f.extcur = fscanf(fid,'%g',f.nextcur)
        f.curlabel = fscanf(fid,'%g',f.nextcur)
    # end

    data = fscanf(fid,'%g%g',[2, f.nstore_seq])
    f.sqt = data[0,:]
    f.wdot = data[1,:]

    #Convert from Internal Units to Physical Units
    f.mass =  f.mass / dmu0
    f.pres =  f.pres / dmu0
    f.jcuru =  f.jcuru / dmu0
    f.jcurv =  f.jcurv / dmu0
    f.jdotb =  f.jdotb / dmu0
    f.phi   = -f.phi    #Data and MSE Fits

    if (f.ireconstruct > 0):
        if (f.imse >= 2) or (f.itse >0):
            f.twsgt = fscanf(fid, '%g', 1)
            f.msewgt = fscanf(fid, '%g', 1)
            f.isnodes = fscanf(fid, '%d', 1)

            data = fscanf(fid,'%g%g%g',[3, f.isnodes])
            f.sknots = data[0, :]
            f.ystark = data[1, :]
            f.y2stark = data[2, :]
            f.ipnodes = fscanf(fid, '%d', 1)

            data = fscanf(fid,'%g%g%g',[3, f.ipnodes])
            f.pknots = data[0, :]
            f.ythom = data[1, :]
            f.y2thom = data[2, :]

            data = fscanf(fid,'%g%g%g%g%g%g%g', [7, (2*f.ns)-1])
            f.anglemse = data[0, :]
            f.rmid = data[1, :]
            f.qmid = data[2, :]
            f.shear = data[3, :]
            f.presmid = data[4, :]
            f.alfa = data[5, :]
            f.curmid = data[6, :]

            data = fscanf(fid,'%g%g%g',[3, f.imse])
            f.rstark = data[0, :]
            f.datastark = data[1, :]
            f.qmeas = data[2, :]

            data = fscanf(fid,'%g%g',[2, f.itse])
            f.rthom = data[0, :]
            f.datathom = data[1, :]
        # end

        if (f.nobd > 0):
            data = fscanf(fid,'%g%g%g',[3, f.nobd])
            f.dsiext = data[0, :]
            f.plflux = data[1, :]
            f.dsiobt = data[2, :]
            f.flmwgt = fscanf(fid, '%g', 1)
        # end
        nbfldn=_np.sum(nbldf[:f.nbsets])

        if (nbfldn > 0):
            for nn in range(nbsets): # n=1:nbsets
                data = fscanf(fid,'%g%g%g', [3, f.nbfld[nn]])
                f.bcoil[:, nn] = data[0, :]
                f.plbfld[:, nn] = data[1, :]
                f.bbc[:, nn] = data[2, :]
            # end
            f.bcwgt = fscanf(fid, '%g', 1)
        # end
        f.phidiam = fscanf(fid, '%g', 1)
        f.delphid = fscanf(fid, '%g', 1)

        #Read Limiter and Prout Plotting Specs
        f.nsets = fscanf(fid, '%g', 1)
        f.nparts = fscanf(fid, '%g', 1)
        f.nlim = fscanf(fid, '%g', 1)
        f.nsetsn = fscanf(fid, '%g', f.nsets)
        f.pfcspec=_np.zeros( (f.nparts, _np.max(f.nxsetsn), f.nsets), dtype=_np.float64)

        for kk in range(f.nsets): # k=1:f.nsets
            for jj in range(f.nsetsn[kk]): # j=1:f.nsetsn(k)
                for ii in range(f.nparts): # i=1:f.nparts
                    f.pfcspec[ii, jj, kk] = fscanf(fid, '%g', 1)
                # end
            # end
        # end

        f.limitr = fscanf(fid,'%g',f.nlim)
        f.rlim = _np.zeros( (_np.max(f.limitr), f.nlim), dtype=_np.float64)
        f.zlim = _np.zeros( (_np.max(f.limitr), f.nlim), dtype=_np.float64)

        for jj in range(f.nlim): # j=1:f.nlim
            for ii in range(f.limitr[jj]): #i=1:f.limitr(j)
                data = fscanf(fid,'%g%g',2)
                f.rlim[ii, jj] = data[1]
                f.zlim[ii, jj] = data[2]
            # end
        # end
        f.nrgrid = fscanf(fid, '%g', 1)
        f.nzgrid = fscanf(fid, '%g', 1)
        f.tokid = fscanf(fid, '%g', 1)
        f.rx1 = fscanf(fid, '%g', 1)
        f.rx2 = fscanf(fid, '%g', 1)
        f.zy1 = fscanf(fid, '%g', 1)
        f.zy2 = fscanf(fid, '%g', 1)
        f.conif = fscanf(fid, '%g', 1)
        f.imatch_phiedge = fscanf(fid, '%g', 1)
    # end
    return f
# end  def read_vmec_605(fid,fmt)

#---------------------------------------------------------------- %

def read_vmec_620(fid,fmt):
    #For 6.20
    f = Struct()

    #Unit Conversions
    dmu0=(2.0e-7)/_np.pi

    #Read Data
    data = fscanf(fid,'%g',6)
    f.wb = data[0]
    f.wp = data[1]
    f.gamma = data[2]
    f.pfac = data[3]
    f.rmax_surf = data[4]
    f.rmin_surf = data[5]

    data = fscanf(fid,'%d',10)
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

    data = fscanf(fid,'%d',5)
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
        f.nbfld = fscanf(fid, '%g', f.nbsets)
    # end

    #Read mgrid filename
    f.mgrid_file = fscanf(fid, '%s', 1)

    #Read Arrays
    if (f.iasym > 0):
        data = fscanf(fid, '%d%d%g%g%g%g%g%g%g%g%g%g%g%g%g%g', [16, f.mnmax])
        data = fscanf(fid,'%g%g%g%g%g%g%g%g%g%g%g%g%g%g',[14, f.mnmax*( f.ns-1 )])
    else:
        data = fscanf(fid, '%d%d%g%g%g%g%g%g%g%g%g%g%g', [13, f.mnmax])
        data = fscanf(fid,'%g%g%g%g%g%g%g%g%g%g%g',[11, f.mnmax*( f.ns-1 )])
    # end

    #Extract Data from Arrays
    f.xm = data1[0, :]
    f.xn = data1[1, :]

    #First reshape data
    f.rmnc=[data1[2, :].T, reshape(data[0, :], f.mnmax, f.ns-1)]
    f.zmns=[data1[3, :].T, reshape(data[1, :], f.mnmax, f.ns-1)]
    f.lmns=[data1[4, :].T, reshape(data[2, :], f.mnmax, f.ns-1)]
    f.bmn=[data1[5, :].T, reshape(data[3, :], f.mnmax, f.ns-1)]
    f.gmn=[data1[6, :].T, reshape(data[4, :], f.mnmax, f.ns-1)]
    f.bsubumnc=[data1[7, :].T, reshape(data[5, :], f.mnmax, f.ns-1)]
    f.bsubvmnc=[data1[8, :].T, reshape(data[6, :], f.mnmax, f.ns-1)]
    f.bsubsmns=[data1[9, :].T, reshape(data[7, :], f.mnmax, f.ns-1)]
    f.bsupumnc=[data1[10, :].T, reshape(data[8, :], f.mnmax, f.ns-1)]
    f.bsupvmnc=[data1[11, :].T, reshape(data[9, :], f.mnmax, f.ns-1)]
    f.currvmnc=[data1[12, :].T, reshape(data[10, :], f.mnmax, f.ns-1)]

    #Read the half-mesh quantities
    data = fscanf(fid,'%g%g%g%g%g%g%g%g%g%g%g%g%g', [13, f.ns/2])
    f.iotas = data[0, :]
    f.mass = data[1, :]
    f.pres = data[2, :]
    f.beta_vol[3, :]
    f.phip = data[4, :]
    f.buco = data[5, :]
    f.bvco = data[6, :]
    f.phi = data[7, :]
    f.vp = data[8, :]
    f.overr = data[9, :]
    f.jcuru = data[10, :]
    f.jcurv = data[11, :]
    f.specw = data[12, :]

    data = fscanf(fid,'%g%g%g%g%g%g',6)
    f.aspect = data[0]
    f.betatot = data[1]
    f.betapol = data[2]
    f.betator = data[3]
    f.betaxis = data[4]
    f.b0 = data[5]

    f.isigna = fscanf(fid, '%d\n', 1)
    f.input_extension=strtrim(fgetl(fid))

    data = fscanf(fid,'%g',8)
    f.IonLarmor = data[0]
    f.VolAvgB = data[1]
    f.RBtor0 = data[2]
    f.RBtor = data[3]
    f.Itor = data[4]
    f.Aminor = data[5]
    f.Rmajor = data[6]
    f.Volume = data[7]

    #Mercier Criterion
    data = fscanf(fid,'%g%g%g%g%g%g', [6, f.ns-2])
    f.Dmerc = data[0, :]
    f.Dshear = data[1, :]
    f.Dwell = data[2, :]
    f.Dcurr = data[3, :]
    f.Dgeod = data[4, :]
    f.equif = data[5, :]

    if (f.nextcur > 0):
        f.extcur = fscanf(fid, '%g', f.nextcur)
        f.curlabel = fscanf(fid, '%g', f.nextcur)
    # end

    data = fscanf(fid,'%g%g', [2, f.nstore_seq])
    f.sqt = data[0, :]
    f.wdot = data[1, :]

    data = fscanf(fid,'%g%g', [2, f.nstore_seq])
    f.jdotb = data[0, :] #changed indices
    f.bdotgradv = data[1, :]

    #Data and MSE Fits
    if (f.ireconstruct > 0):
        if ((f.imse >= 2) or (f.itse >0)):
            f.twsgt = fscanf(fid, '%g', 1)
            f.msewgt = fscanf(fid, '%g', 1)
            f.isnodes = fscanf(fid, '%d', 1)

            data = fscanf(fid,'%g%g%g', [3, f.isnodes])
            f.sknots = data[0, :] #changed indices
            f.ystark = data[1, :]
            f.y2stark = data[2, :]
            f.ipnodes = fscanf(fid, '%d', 1)

            data = fscanf(fid,'%g%g%g', [3, f.ipnodes])
            f.pknots = data[0, :] #changed indices
            f.ythom = data[1, :]
            f.y2thom = data[2, :]

            data = fscanf(fid,'%g%g%g%g%g%g%g', [7, (2*f.ns)-1])
            f.anglemse = data[0, :] #changed indices
            f.rmid = data[1, :]
            f.qmid = data[2, :]
            f.shear = data[3, :]
            f.presmid = data[4, :]
            f.alfa = data[5, :]
            f.curmid = data[6, :]

            data = fscanf(fid,'%g%g%g', [3, f.imse])
            f.rstark = data[0, :] #changed indices
            f.datastark = data[1, :]
            f.qmeas = data[2, :]

            data = fscanf(fid,'%g%g', [2, f.itse])
            f.rthom = data[0, :] #changed indices
            f.datathom = data[1, :]
        # end

        if (f.nobd > 0):
            data = fscanf(fid,'%g%g%g',[3,f.nobd])
            f.dsiext = data[0, :] #changed indices
            f.plflux = data[1, :]
            f.dsiobt = data[2, :]
            f.flmwgt = fscanf(fid, '%g', 1)
        # end
        nbfldn=_np.sum(nbldf[:f.nbsets])

        if (nbfldn > 0):
            for nn in range(nbsets): # n=1:nbsets
                data = fscanf(fid,'%g%g%g', [3, f.nbfld[nn]])
                f.bcoil[:, nn] = data[0, :] #changed indices
                f.plbfld[:, nn] = data[1, :]
                f.bbc[:, nn] = data[2, :]
            # end
            f.bcwgt = fscanf(fid, '%g', 1)
        # end
        f.phidiam = fscanf(fid, '%g', 1)
        f.delphid = fscanf(fid, '%g', 1)

        #Read Limiter and Prout Plotting Specs
        f.nsets = fscanf(fid, '%g', 1)
        f.nparts = fscanf(fid, '%g', 1)
        f.nlim = fscanf(fid, '%g', 1)
        f.nsetsn = fscanf(fid,'%g',f.nsets)
        f.pfcspec=_np.zeros( (f.nparts, _np.max(f.nxsetsn), f.nsets), dtype=_np.float64)
        for kk in range(f.nsets): # k=1:f.nsets
            for jj in range(f.nsetsn[kk]): # j=1:f.nsetsn(k)
                for ii in range(f.nparts): # i=1:f.nparts
                    f.pfcspec[ii, jj, kk] = fscanf(fid, '%g', 1)
                # end
            # end
        # end
        f.limitr = fscanf(fid, '%g', f.nlim)
        f.rlim=_np.zeros( (_np.max(f.limitr),f.nlim), dtype=_np.float64)
        f.zlim=_np.zeros( (_np.max(f.limitr),f.nlim), dtype=_np.float64)

        for jj in range(f.nlim): # j=1:f.nlim
            for ii in range(f.limitr[jj]): # i=1:f.limitr(j)
                data = fscanf(fid,'%g%g',2)
                f.rlim[ii, jj] = data[0] #changed indices
                f.zlim[ii, jj] = data[1]
            # end
        # end
        f.nrgrid = fscanf(fid, '%g', 1)
        f.nzgrid = fscanf(fid, '%g', 1)
        f.tokid = fscanf(fid, '%g', 1)
        f.rx1 = fscanf(fid, '%g', 1)
        f.rx2 = fscanf(fid, '%g', 1)
        f.zy1 = fscanf(fid, '%g', 1)
        f.zy2 = fscanf(fid, '%g', 1)
        f.conif = fscanf(fid, '%g', 1)
        f.imatch_phiedge = fscanf(fid, '%g', 1)
    # end
    return f
# end  def read_vmec_620()

#---------------------------------------------------------------- %

def read_vmec_650(fid, fmt):
    #function f=read_vmec_650(fid, fmt)
    #For 6.50
    f = Struct()

    #Unit Conversions
    dmu0=(2.0e-7)/_np.pi

    #Read Data
    data = fscanf(fid, '%g', 6)
    f.wb = data[0]
    f.wp = data[1]
    f.gamma = data[2]
    f.pfac = data[3]
    f.rmax_surf = data[4]
    f.rmin_surf = data[5]

    data = fscanf(fid, '%d', 10)
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

    data = fscanf(fid, '%d', 6)
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
        f.nbfld = fscanf(fid, '%g', f.nbsets)
    # end

    #Read mgrid filename
    f.mgrid_file = fscanf(fid, '%s', 1)

    #Read Arrays
    if f.iasym > 0:
        data = fscanf(fid, '%d%d%g%g%g%g%g%g%g%g%g%g%g%g%g%g', [16, f.mnmax])
        data = fscanf(fid, '%g%g%g%g%g%g%g%g%g%g%g%g%g%g', [14, f.mnmax*(f.ns-1)])
    else:
        data = fscanf(fid, '%d%d%g%g%g%g%g%g%g%g%g%g%g', [13, f.mnmax])
        data = fscanf(fid, '%g%g%g%g%g%g%g%g%g%g%g', [11, f.mnmax*(f.ns-1)])
    # end

    #Extract Data from Arrays
    f.xm = data1[0, :]
    f.xn = data1[1, :]

    #First reshape data
    f.rmnc=[data1[2, :].T, reshape(data[0, :],f.mnmax,f.ns-1)]
    f.zmns=[data1[3, :].T, reshape(data[1, :],f.mnmax,f.ns-1)]
    f.lmns=[data1[4, :].T, reshape(data[2, :],f.mnmax,f.ns-1)]
    f.bmnc=[data1[5, :].T, reshape(data[3, :],f.mnmax,f.ns-1)]
    f.gmnc=[data1[6, :].T, reshape(data[4, :],f.mnmax,f.ns-1)]
    f.bsubumnc=[data1[7, :].T, reshape(data[5, :],f.mnmax,f.ns-1)]
    f.bsubvmnc=[data1[8, :].T, reshape(data[6, :],f.mnmax,f.ns-1)]
    f.bsubsmns=[data1[9, :].T, reshape(data[7, :],f.mnmax,f.ns-1)]
    f.bsupumnc=[data1[10, :].T, reshape(data[8, :],f.mnmax,f.ns-1)]
    f.bsupvmnc=[data1[11, :].T, reshape(data[9, :],f.mnmax,f.ns-1)]
    f.currvmnc=[data1[12, :].T, reshape(data[10, :],f.mnmax,f.ns-1)]

    #Read the half-mesh quantities
    data = fscanf(fid, '%g%g%g%g%g%g%g%g%g%g%g%g%g', [13, f.ns/2])
    f.iotas = data[0, :]
    f.mass = data[1, :]
    f.pres = data[2, :]
    f.beta_vol[3, :]
    f.phip = data[4, :]
    f.buco = data[5, :]
    f.bvco = data[6, :]
    f.phi = data[7, :]
    f.vp = data[8, :]
    f.overr = data[9, :]
    f.jcuru = data[10, :]
    f.jcurv = data[11, :]
    f.specw = data[12, :]

    data = fscanf(fid, '%g%g%g%g%g%g', 6)
    f.aspect = data[0]
    f.betatot = data[1]
    f.betapol = data[2]
    f.betator = data[3]
    f.betaxis = data[4]
    f.b0 = data[5]

    f.isigna = fscanf(fid, '%d\n', 1)
    f.input_extension=strtrim(fgetl(fid))

    data = fscanf(fid, '%g', 8)
    f.IonLarmor = data[0]
    f.VolAvgB = data[1]
    f.RBtor0 = data[2]
    f.RBtor = data[3]
    f.Itor = data[4]
    f.Aminor = data[5]
    f.Rmajor = data[6]
    f.Volume = data[7]

    #Mercier Criterion
    data = fscanf(fid, '%g%g%g%g%g%g', [6, f.ns-2])
    f.Dmerc = data[0, :]
    f.Dshear = data[1, :]
    f.Dwell = data[2, :]
    f.Dcurr = data[3, :]
    f.Dgeod = data[4, :]
    f.equif = data[5, :]

    if (f.nextcur > 0):
        f.extcur = fscanf(fid, '%g', f.nextcur)
        f.curlabel = fscanf(fid, '%g', f.nextcur)
    # end

    data = fscanf(fid, '%g%g', [2, f.nstore_seq])
    f.sqt = data[0, :]
    f.wdot = data[1, :]

    data = fscanf(fid, '%g%g', [2, f.nstore_seq])
    f.jdotb = data[0, :]
    f.bdotgradv = data[1, :]

    #Data and MSE Fits
    if (f.ireconstruct > 0):
        if ((f.imse >= 2) or (f.itse > 0)):
            f.twsgt = fscanf(fid, '%g', 1)
            f.msewgt = fscanf(fid, '%g', 1)
            f.isnodes = fscanf(fid, '%d', 1)

            data = fscanf(fid, '%g%g%g', [3, f.isnodes])
            f.sknots = data[0, :]
            f.ystark = data[1, :]
            f.y2stark = data[2, :]
            f.ipnodes = fscanf(fid, '%d', 1)

            data = fscanf(fid, '%g%g%g', [3, f.ipnodes])
            f.pknots = data[0, :]
            f.ythom = data[1, :]
            f.y2thom = data[2, :]

            data = fscanf(fid, '%g%g%g%g%g%g%g', [7, (2*f.ns)-1])
            f.anglemse = data[0, :]
            f.rmid = data[1, :]
            f.qmid = data[2, :]
            f.shear = data[3, :]
            f.presmid = data[4, :]
            f.alfa = data[5, :]
            f.curmid = data[6, :]

            data = fscanf(fid, '%g%g%g', [3, f.imse])
            f.rstark = data[0, :]
            f.datastark = data[1, :]
            f.qmeas = data[2, :]

            data = fscanf(fid, '%g%g', [2, f.itse])
            f.rthom = data[0, :]
            f.datathom = data[1, :]
        # end

        if (f.nobd > 0):
            data = fscanf(fid, '%g%g%g', [3,f.nobd])
            f.dsiext = data[0, :]
            f.plflux = data[1, :]
            f.dsiobt = data[2, :]
            f.flmwgt = fscanf(fid, '%g', 1)
        # end

        nbfldn = _np.sum(nbldf[:f.nbsets])
        if (nbfldn > 0):
            for nn in range(nbsets):  #for n=1:nbsets
                data = fscanf(fid, '%g%g%g', [3, f.nbfld[nn]])
                f.bcoil[:, nn] = data[0, :]
                f.plbfld[:, nn] = data[1, :]
                f.bbc[:, nn] = data[2, :]
            # end
            f.bcwgt = fscanf(fid, '%g', 1)
        # end
        f.phidiam = fscanf(fid, '%g', 1)
        f.delphid = fscanf(fid, '%g', 1)

        #Read Limiter and Prout Plotting Specs
        f.nsets = fscanf(fid, '%g', 1)
        f.nparts = fscanf(fid, '%g', 1)
        f.nlim = fscanf(fid, '%g', 1)
        f.nsetsn = fscanf(fid, '%g', f.nsets)
        f.pfcspec=_np.zeros( (f.nparts,_np.max(f.nxsetsn),f.nsets), dtype=_np.float64)
        for kk in range(f.nsets):  #k=1:f.nsets
            for jj in range(f.nsetsn[kk]):  #for j=1:f.nsetsn(k)
                for ii in range(f.nparts):  #for i=1:f.nparts
                    f.pfcspec[ii, jj, kk] = fscanf(fid, '%g', 1)
                # end
            # end
        # end

        f.limitr = fscanf(fid, '%g', f.nlim)
        f.rlim=_np.zeros( (_np.max(f.limitr), f.nlim), dtype=_np.float64)
        f.zlim=_np.zeros( (_np.max(f.limitr), f.nlim), dtype=_np.float64)
        for jj in range(f.nlim):  #j=1:f.nlim
            for ii in range(f.limitr[jj]):  #i=1:f.limitr(j)
                data = fscanf(fid, '%g%g', 2)
                f.rlim[ii, jj] = data[0]
                f.zlim[ii, jj] = data[1]
            # end
        # end
        f.nrgrid = fscanf(fid, '%g', 1)
        f.nzgrid = fscanf(fid, '%g', 1)
        f.tokid = fscanf(fid, '%g', 1)
        f.rx1 = fscanf(fid, '%g', 1)
        f.rx2 = fscanf(fid, '%g', 1)
        f.zy1 = fscanf(fid, '%g', 1)
        f.zy2 = fscanf(fid, '%g', 1)
        f.conif = fscanf(fid, '%g', 1)
        f.imatch_phiedge = fscanf(fid, '%g', 1)
    # end
    return f
# end  def read_vmec_650()

#---------------------------------------------------------------- %

def read_vmec_695(fid, fmt):
    #function f=read_vmec_695(fid,fmt)
    f = Struct()
    f.lfreeb = 0

    data = fscanf(fid, fmt, 7)
    f.wb = data[0]
    f.wp = data[1]
    f.gamma = data[2]
    f.pfac = data[3]
    f.rmax_surf = data[4]
    f.rmin_surf = data[5]
    f.zmax_surf = data[6]

    data = fscanf(fid, fmt, 10)
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

    data = fscanf(fid, fmt, 6)
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

    #Read mgrid filename and setup other format statements
    if fmt.find(',')<0:
        #  f.mgrid_file = fscanf(fid, '%s', 1)
        fmt2 = '%d%d%g%g%g%g%g%g%g%g%g%g%g%g%g%g'
        fmt3 = '%g%g%g%g%g%g%g%g%g%g%g%g%g%g'
        fmt4 = '%d%d%g%g%g%g%g%g%g%g%g%g%g'
        fmt5 = '%g%g%g%g%g%g%g%g%g%g%g'
        fmt6 = '%g%g%g%g%g%g%g%g%g%g%g%g%g'
        fmt7 = '%g%g%g%g%g%g'
        fmt8 = '%g%g'
        fmt9 = '%g%g%g'
        fmt10 = '%g%g%g%g%g%g%g'
    else:
        #   f.mgrid_file = fscanf(fid, '%s,', 1)
        fmt2 = '%d,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g'
        fmt3 = '%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g'
        fmt4 = '%d,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g'
        fmt5 = '%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g'
        fmt6 = '%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g'
        fmt7 = '%g,%g,%g,%g,%g,%g'
        fmt8 = '%g,%g'
        fmt9 = '%g,%g,%g'
        fmt10 = '%g,%g,%g,%g,%g,%g,%g'
    # end

    if fmt.find(',')<0:  #isempty(strfind(fmt,','))
        f.mgrid_file = fscanf(fid,'%s',1)
        fmt2 = '%g%g'
        fmt3 = '%g%g%g'
        fmt6 = '%g%g%g%g%g%g'
        fmt7 = '%g%g%g%g%g%g%g'
        fmt11 = '%g%g%g%g%g%g%g%g%g%g%g'
        fmt12 = '%g%g%g%g%g%g%g%g%g%g%g%g'
        fmt13 = '%g%g%g%g%g%g%g%g%g%g%g%g%g'
        fmt14 = '%g%g%g%g%g%g%g%g%g%g%g%g%g%g'
        fmt2_11 = '%d%d%g%g%g%g%g%g%g%g%g%g%g'
        fmt2_14 = '%d%d%g%g%g%g%g%g%g%g%g%g%g%g%g%g'
    else:
        f.mgrid_file = fscanf(fid, '%s,', 1)
        fmt2 = '%g,%g'
        fmt3 = '%g,%g,%g'
        fmt6 = '%g,%g,%g,%g,%g,%g'
        fmt7 = '%g,%g,%g,%g,%g,%g,%g'
        fmt11 = '%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g'
        fmt12 = '%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g'
        fmt13 = '%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g'
        fmt14 = '%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g'
        fmt2_11 = '%d,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g'
        fmt2_14 = '%d,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g'
    # end

    # Read Arrays
    if f.iasym > 0:
        data = fscanf(fid, fmt2_14, [16, f.mnmax])
        data = fscanf(fid, fmt14, [14, f.mnmax*(f.ns-1)])
    else:
        data = fscanf(fid, fmt2_11, [13, f.mnmax])
        data = fscanf(fid, fmt11, [11, f.mnmax*(f.ns-1)])
    # end

    # Extract Data from Arrays
    f.xm = data1[0, :]
    f.xn = data1[1, :]

    # First reshape data
    f.rmnc=[data1[2, :].T, reshape(data[0, :],f.mnmax,f.ns-1)]
    f.zmns=[data1[3, :].T, reshape(data[1, :],f.mnmax,f.ns-1)]
    f.lmns=[data1[4, :].T, reshape(data[2, :],f.mnmax,f.ns-1)] # On half grid
    f.bmnc=[data1[5, :].T, reshape(data[3, :],f.mnmax,f.ns-1)]
    f.gmnc=[data1[6, :].T, reshape(data[4, :],f.mnmax,f.ns-1)]
    f.bsubumnc=[data1[7, :].T, reshape(data[5, :],f.mnmax,f.ns-1)]
    f.bsubvmnc=[data1[8, :].T, reshape(data[6, :],f.mnmax,f.ns-1)]
    f.bsubsmns=[data1[9, :].T, reshape(data[7, :],f.mnmax,f.ns-1)]
    f.bsupumnc=[data1[10, :].T, reshape(data[8, :],f.mnmax,f.ns-1)]
    f.bsupvmnc=[data1[11, :].T, reshape(data[9, :],f.mnmax,f.ns-1)]
    f.currvmnc=[data1[12, :].T, reshape(data[10, :],f.mnmax,f.ns-1)]
    if f.iasym > 0:
        f.rmns=[data1[13, :].T, reshape(data[11, :],f.mnmax,f.ns-1)]
        f.zmnc=[data1[14, :].T, reshape(data[12, :],f.mnmax,f.ns-1)]
        f.lmnc=[data1[15, :].T, reshape(data[13, :],f.mnmax,f.ns-1)] # On half grid
    # end

    # Read the half-mesh quantities
    data = fscanf(fid,fmt13, [13, f.ns-1])
    f.iotas = data[0, :]
    f.mass = data[1, :]
    f.pres = data[2, :]
    f.beta_vol = data[3, :]
    f.phip = data[4, :]
    f.buco = data[5, :]
    f.bvco = data[6, :]
    f.phi = data[7, :]
    f.vp = data[8, :]
    f.overr = data[9, :]
    f.jcuru = data[10, :]
    f.jcurv = data[11, :]
    f.specw = data[12, :]

    data = fscanf(fid, fmt, 6)
    f.aspect = data[0]
    f.betatot = data[1]
    f.betapol = data[2]
    f.betator = data[3]
    f.betaxis = data[4]
    f.b0 = data[5]

    f.isigna = fscanf(fid, fmt+'\n', 1)
    f.input_extension = strtrim(fgetl(fid))

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
    f.Dmerc = data[0, :]
    f.Dshear = data[1, :]
    f.Dwell = data[2, :]
    f.Dcurr = data[3, :]
    f.Dgeod = data[4, :]
    f.equif = data[5, :]

    f.curlabel=cell(f.nextcur, 1)
    if (f.nextcur > 0):
        f.lfreeb = 1
        f.extcur = fscanf(fid, fmt, f.nextcur)
        fscanf(fid, '\n')
        rem = f.nextcur
        jj = 0
        while rem > 0:
            line = fgetl(fid)
            fscanf(fid,'\n')
            test = line[0]
            index = line.find(test)  # findstr(line,test)
            for ii in range(_np.size(index, axis=1)/2):  # i=1:size(index,2)/2
                f.curlabel[ii+jj] = strtrim(line[index[2*ii-1]+range(index[2*ii])-1])
            # end
            jj += _np.size(index, axis=1)/2
            rem -= _np.size(index, axis=1)/2
        # end
    # end

    data = fscanf(fid, fmt2, [2, f.nstore_seq])
    f.sqt = data[0, :]
    f.wdot = data[1, :]

    data = fscanf(fid, fmt2, [2, f.nstore_seq])
    f.jdotb = data[0, :]
    f.bdotgradv = data[1, :]

    # No Unit Conversion Necessary
    # Data and MSE Fits
    if (f.ireconstruct > 0):
        if ((f.imse >= 2) or (f.itse > 0)):
            f.twsgt = fscanf(fid, fmt, 1)
            f.msewgt = fscanf(fid, fmt, 1)
            f.isnodes = fscanf(fid, fmt, 1)

            data = fscanf(fid, fmt3, [3, f.isnodes])
            f.sknots = data[0, :]
            f.ystark = data[1, :]
            f.y2stark = data[2, :]
            f.ipnodes = fscanf(fid, fmt, 1)

            data = fscanf(fid, fmt3, [3, f.ipnodes])
            f.pknots = data[0, :]
            f.ythom = data[1, :]
            f.y2thom = data[2, :]

            data = fscanf(fid,fmt7, [7, (2*f.ns)-1])
            f.anglemse = data[0, :]
            f.rmid = data[1, :]
            f.qmid = data[2, :]
            f.shear = data[3, :]
            f.presmid = data[4, :]
            f.alfa = data[5, :]
            f.curmid = data[6, :]

            data = fscanf(fid, fmt3, [3, f.imse])
            f.rstark = data[0, :]
            f.datastark = data[1, :]
            f.qmeas = data[2, :]

            data = fscanf(fid, fmt2, [2, f.itse])
            f.rthom = data[0, :]
            f.datathom = data[1, :]
        # end

        if (f.nobd > 0):
            data = fscanf(fid, fmt3, [3, f.nobd])
            f.dsiext = data[0, :]
            f.plflux = data[1, :]
            f.dsiobt = data[2, :]
            f.flmwgt = fscanf(fid, fmt, 1)
        # end

        nbfldn = _np.sum(nbldf[:f.nbsets])
        if (nbfldn > 0):
            for nn in range(nbsets):  # for n=1:nbsets
                data = fscanf(fid, fmt3, [3, f.nbfld[nn]])
                f.bcoil[:, nn] = data[0, :]
                f.plbfld[:, nn] = data[1, :]
                f.bbc[:, nn] = data[2, :]
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
        f.zlim = _np.zeros( (_np.max(f.limitr), f.nlim), dtype=_np.float64)
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

# ---------------------------------------------------------------- %

def read_vmec_800(fid, fmt):
    # function f=read_vmec_800(fid,fmt)
    f = Struct()
    f.lfreeb = 0
    data = fscanf(fid, fmt, 7)
    f.wb = data[0]
    f.wp = data[1]
    f.gamma = data[2]
    f.pfac = data[3]
    f.rmax_surf = data[4]
    f.rmin_surf = data[5]
    f.zmax_surf = data[6]

    data = fscanf(fid, fmt, 10)
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

    data = fscanf(fid, fmt, 6)
    f.imse = data[0]
    f.itse = data[1]
    f.nbsets = data[2]
    f.nobd = data[3]
    f.nextcur = data[4]
    f.nstore_seq = data[5]

    # Error Check
    if (f.ierr_vmec and (f.ierr_vmec != 4)):
        print(strcat('ierr_vmec:',num2str(f.ierr_vmec)))
        return
    # end

    # Read nbfld
    if (f.nbsets > 0):
        f.nbfld = fscanf(fid, fmt, f.nbsets)
    # end

    # Read mgrid filename and setup other format statements
    if fmt.find(',')<0:  # isempty(strfind(fmt,','))
        f.mgrid_file = fscanf(fid,'%s',1)
        fmt2 = '%g%g'
        fmt3 = '%g%g%g'
        fmt6 = '%g%g%g%g%g%g'
        fmt7 = '%g%g%g%g%g%g%g'
        fmt10 = '%g%g%g%g%g%g%g%g%g%g'
        fmt11 = '%g%g%g%g%g%g%g%g%g%g%g'
        fmt12 = '%g%g%g%g%g%g%g%g%g%g%g%g'
        fmt13 = '%g%g%g%g%g%g%g%g%g%g%g%g%g'
        fmt14 = '%g%g%g%g%g%g%g%g%g%g%g%g%g%g'
        fmt2_11 = '%d%d%g%g%g%g%g%g%g%g%g%g%g'
        fmt2_14 = '%d%d%g%g%g%g%g%g%g%g%g%g%g%g%g%g'
    else:
        f.mgrid_file = fscanf(fid, '%s,', 1)
        fmt2 = '%g,%g'
        fmt3 = '%g,%g,%g'
        fmt6 = '%g,%g,%g,%g,%g,%g'
        fmt7 = '%g,%g,%g,%g,%g,%g,%g'
        fmt11 = '%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g'
        fmt12 = '%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g'
        fmt13 = '%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g'
        fmt14 = '%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g'
        fmt2_11 = '%d,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g'
        fmt2_14 = '%d,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g'
    # end

    # Read Arrays
    if f.iasym > 0:
        data = fscanf(fid, fmt2_14, [16, f.mnmax])
        data = fscanf(fid, fmt14, [14, f.mnmax*(f.ns-1)])
    else:
        data = fscanf(fid, fmt2_11, [13, f.mnmax])
        data = fscanf(fid, fmt11, [11, f.mnmax*(f.ns-1)])
    # end
    # Extract Data from Arrays
    f.xm = data1[0, :]
    f.xn = data1[1, :]

    # First reshape data
    f.rmnc=[data1[2, :].T, reshape(data[0, :], f.mnmax,f.ns-1)]
    f.zmns=[data1[3, :].T, reshape(data[1, :], f.mnmax,f.ns-1)]
    f.lmns=[data1[4, :].T, reshape(data[2, :], f.mnmax,f.ns-1)] # On half grid
    f.bmnc=[data1[5, :].T, reshape(data[3, :], f.mnmax,f.ns-1)]
    f.gmnc=[data1[6, :].T, reshape(data[4, :], f.mnmax,f.ns-1)]
    f.bsubumnc=[data1[7, :].T, reshape(data[5, :], f.mnmax,f.ns-1)]
    f.bsubvmnc=[data1[8, :].T, reshape(data[6, :], f.mnmax,f.ns-1)]
    f.bsubsmns=[data1[9, :].T, reshape(data[7, :], f.mnmax,f.ns-1)]
    f.bsupumnc=[data1[10, :].T, reshape(data[8, :], f.mnmax,f.ns-1)]
    f.bsupvmnc=[data1[11, :].T, reshape(data[9, :], f.mnmax,f.ns-1)]
    f.currvmnc=[data1[12, :].T, reshape(data[10, :], f.mnmax,f.ns-1)]
    if f.iasym >0:
        f.rmns = [data1[13, :].T, reshape(data[11, :], f.mnmax,f.ns-1)]
        f.zmnc = [data1[14, :].T, reshape(data[12, :], f.mnmax,f.ns-1)]
        f.lmnc = [data1[15, :].T, reshape(data[13, :], f.mnmax,f.ns-1)] # On half grid
    # end

    # Read the full-mesh quantities
    data = fscanf(fid, fmt6, [6, f.ns])
    f.iotaf = data[0, :]
    f.presf = data[1, :]
    f.phipf = data[2, :]
    f.phi = data[3, :]
    f.jcuru = data[4, :]
    f.jcurv = data[5, :]

    # Read the half-mesh quantities
    data = fscanf(fid, fmt10, [10, f.ns-1])
    f.iotas = data[0,:]
    f.mass = data[1,:]
    f.pres = data[2,:]
    f.beta_vol = data[3, :]
    f.phip = data[4, :]
    f.buco = data[5, :]
    f.bvco = data[6, :]
    f.vp = data[7, :]
    f.overr = data[8, :]
    f.specw = data[9, :]

    data = fscanf(fid, fmt, 6)
    f.aspect = data[0]
    f.betatot = data[1]
    f.betapol = data[2]
    f.betator = data[3]
    f.betaxis = data[4]
    f.b0 = data[5]
    f.isigna = fscanf(fid, fmt+'\n', 1)
    f.input_extension = strtrim(fgetl(fid))

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
    f.Dmerc = data[0, :]
    f.Dshear = data[1, :]
    f.Dwell = data[2, :]
    f.Dcurr = data[3, :]
    f.Dgeod = data[4, :]
    f.equif = data[5, :]

    f.curlabel = cell(f.nextcur, 1)
    if (f.nextcur > 0):
        f.lfreeb = 1
        f.extcur = fscanf(fid, fmt, f.nextcur)
        fscanf(fid, '\n')
        rem = f.nextcur
        jj = 0
        while rem > 0:
            line = fgetl(fid)
            fscanf(fid, '\n')
            test = line[0]
            index = findstr(line, test)
            for ii in range(_np.size(index, axis=1)/2):  # i=1:size(index, 2)/2
                f.curlabel[ii+jj] = strtrim(line[index[2*ii-1]+range(index[2*ii])-1])
            # end
            jj += _np.size(index, axis=1)/2
            rem -= _np.size(index, axis=1)/2
        # end
    # end

    data = fscanf(fid, fmt2, [2, f.nstore_seq])
    f.sqt = data[0, :]
    f.wdot = data[1, :]

    data = fscanf(fid, fmt2, [2, f.nstore_seq])
    f.jdotb = data[0, :]
    f.bdotgradv = data[1, :]

    # No Unit Conversion Necessary
    # Data and MSE Fits
    if (f.ireconstruct > 0):
        if ((f.imse >= 2) or (f.itse > 0)):
            f.twsgt = fscanf(fid, fmt, 1)
            f.msewgt = fscanf(fid, fmt, 1)
            f.isnodes = fscanf(fid, fmt, 1)

            data = fscanf(fid, fmt3, [3, f.isnodes])
            f.sknots = data[0, :]
            f.ystark = data[1, :]
            f.y2stark = data[2, :]
            f.ipnodes = fscanf(fid, fmt, 1)

            data = fscanf(fid, fmt3, [3, f.ipnodes])
            f.pknots = data[0, :]
            f.ythom = data[1, :]
            f.y2thom = data[2, :]

            data = fscanf(fid,fmt7, [7, (2*f.ns)-1])
            f.anglemse = data[0, :]
            f.rmid = data[1, :]
            f.qmid = data[2, :]
            f.shear = data[3, :]
            f.presmid = data[4, :]
            f.alfa = data[5, :]
            f.curmid = data[6, :]

            data = fscanf(fid, fmt3, [3, f.imse])
            f.rstark = data[0, :]
            f.datastark = data[1, :]
            f.qmeas = data[2, :]

            data = fscanf(fid, fmt2, [2, f.itse])
            f.rthom = data[0, :]
            f.datathom = data[1, :]
        # end

        if (f.nobd > 0):
            data = fscanf(fid, fmt3, [3,f.nobd])
            f.dsiext = data[0, :]
            f.plflux = data[1, :]
            f.dsiobt = data[2, :]
            f.flmwgt = fscanf(fid, fmt, 1)
        # end

        nbfldn = _np.sum(nbldf[:f.nbsets])
        if (nbfldn > 0):
            for nn in range(nbsets):  # for n=1:nbsets
                data = fscanf(fid, fmt3, [3, f.nbfld[nn]])
                f.bcoil[:, nn] = data[0, :]
                f.plbfld[:, nn] = data[1, :]
                f.bbc[:, nn] = data[2, :]
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
        f.zlim = _np.zeros( (_np.max(f.limitr), f.nlim), dtype=_np.float64)
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
# end  def read_vmec_800()

# ---------------------------------------------------------------- %

def read_vmec_847(fid,fmt):
    f = Struct()
    f.lfreeb = 0

    data = fscanf(fid, fmt, 7)
    f.wb = data[0]
    f.wp = data[1]
    f.gamma = data[2]
    f.pfac = data[3]
    f.rmax_surf = data[4]
    f.rmin_surf = data[5]
    f.zmax_surf = data[6]

    data = fscanf(fid, fmt, 11)
    f.nfp = data[0]
    f.ns = data[1]
    f.mpol = data[2]
    f.ntor = data[3]
    f.mnmax = data[4]
    f.mnmax_nyq = data[5]
    f.itfsq = data[5]
    f.niter = data[6]
    f.iasym = data[8]
    f.ireconstruct = data[9]
    f.ierr_vmec = data[10]

    data = fscanf(fid, fmt, 6)
    f.imse = data[0]
    f.itse = data[1]
    f.nbsets = data[2]
    f.nobd = data[3]
    f.nextcur = data[4]
    f.nstore_seq = data[5]

    # Error Check
    if f.ierr_vmec and (f.ierr_vmec != 4):
        print(strcat('ierr_vmec:',num2str(f.ierr_vmec)))
        return
    # end

    # Read nbfld
    if (f.nbsets > 0):
        f.nbfld = fscanf(fid, fmt, f.nbsets)
    # end

    # Read mgrid filename and setup other format statements
    if fmt.find(',')<0:  # isempty(strfind(fmt,',')):
        f.mgrid_file = fscanf(fid,'%s',1)
        fmt2 = '%g%g'
        fmt3 = '%g%g%g'
        fmt6 = '%g%g%g%g%g%g'
        fmt7 = '%g%g%g%g%g%g%g'
        fmt10 = '%g%g%g%g%g%g%g%g%g%g'
        fmt13 = '%g%g%g%g%g%g%g%g%g%g%g%g%g'
        fmt14 = '%g%g%g%g%g%g%g%g%g%g%g%g%g%g'
        fmt20 = '%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g'
        fmt2_3_2_7 = '%d%d%g%g%g%d%d%g%g%g%g%g%g%g'
        fmt2_6_2_14 = '%d%d%g%g%g%g%g%g%d%d%g%g%g%g%g%g%g%g%g%g%g%g%g%g'
        fmt2_7 = '%d%d%g%g%g%g%g%g%g'
        fmt2_14 = '%d%d%g%g%g%g%g%g%g%g%g%g%g%g%g%g'
    else:
        f.mgrid_file = fscanf(fid,'%s,',1)
        fmt2 = '%g,%g'
        fmt3 = '%g,%g,%g'
        fmt6 = '%g,%g,%g,%g,%g,%g'
        fmt7 = '%g,%g,%g,%g,%g,%g,%g'
        fmt13 = '%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g'
        fmt2_3 = '%d,%d,%g,%g,%g'
        fmt2_6 = '%d,%d,%g,%g,%g,%g,%g,%g'
    # end

    # Read Arrays
    f.xm = _np.zeros( (1,f.mnmax), dtype=_np.int64)
    f.xn = _np.zeros_like(f.xn)

    f.rmnc = _np.zeros( (f.mnmax,f.ns), dtype=_np.float64)
    f.zmns = _np.zeros_like(f.rmnc)
    f.lmns = _np.zeros_like(f.rmnc)

    f.xm_nyq = _np.zeros( (1,f.mnmax_nyq), dtype_np.int64)
    f.xn_nyq = _np.zeros_like(f.xm_nyq)

    f.bmnc = _np.zeros( (f.mnmax_nyq,f.ns), dtype=_np.float64)
    f.gmnc = _np.zeros_like(f.bmnc)
    f.bsubumnc = _np.zeros_like(f.bmnc)
    f.bsubvmnc = _np.zeros_like(f.bmnc)
    f.bsubsmns = _np.zeros_like(f.bmnc)
    f.bsupumnc = _np.zeros_like(f.bmnc)
    f.bsupvmnc = _np.zeros_like(f.bmnc)

    if f.iasym >0:
        f.rmns = _np.zeros( (f.mnmax,f.ns), dtype=_np.float64)
        f.zmnc = _np.zeros_like(f.rmns)
        f.lmnc = _np.zeros_like(f.rmns)

        f.bmns = _np.zeros( (f.mnmax_nyq,f.ns), dtype=_np.float64)
        f.gmns = _np.zeros_like(f.bmns)
        f.bsubumns = _np.zeros_like(f.bmns)
        f.bsubvmns = _np.zeros_like(f.bmns)
        f.bsubsmnc = _np.zeros_like(f.bmns)
        f.bsupumns = _np.zeros_like(f.bmns)
        f.bsupvmns = _np.zeros_like(f.bmns)
    # end

    for ii in range(f.ns): # i=1:f.ns
        for jj in range(f.mnmax): # j=1:f.mnmax
            if ii==1:
                f.xm[jj] = fscanf(fid, '%d', 1)
                f.xn[jj] = fscanf(fid, '%d', 1)
            # end
            f.rmnc[jj, ii] = fscanf(fid, '%g', 1)
            f.zmns[jj, ii] = fscanf(fid, '%g', 1)
            f.lmns[jj, ii] = fscanf(fid, '%g', 1)

            if f.iasym > 0:
                f.rmns[jj, ii] = fscanf(fid, '%g', 1)
                f.zmnc[jj, ii] = fscanf(fid, '%g', 1)
                f.lmnc[jj, ii] = fscanf(fid, '%g', 1)
            # end
        # end

        for jj in range(f.mnmax_nyq): # j=1:f.mnmax_nyq
            if ii==1:
                f.xm_nyq[jj] = fscanf(fid, '%d', 1)
                f.xn_nyq[jj] = fscanf(fid, '%d', 1)
            # end

            f.bmnc[jj, ii] = fscanf(fid, '%g', 1)
            f.gmnc[jj, ii] = fscanf(fid, '%g', 1)
            f.bsubumnc[jj, ii] = fscanf(fid, '%g', 1)
            f.bsubvmnc[jj, ii] = fscanf(fid, '%g', 1)
            f.bsubsmns[jj, ii] = fscanf(fid, '%g', 1)
            f.bsupumnc[jj, ii] = fscanf(fid, '%g', 1)
            f.bsupvmnc[jj, ii] = fscanf(fid, '%g', 1)

            if f.iasym > 0:
                f.bmns[jj, ii] = fscanf(fid, '%g', 1)
                f.gmns[jj, ii] = fscanf(fid, '%g', 1)
                f.bsubumns[jj, ii] = fscanf(fid, '%g', 1)
                f.bsubvmns[jj, ii] = fscanf(fid, '%g', 1)
                f.bsubsmnc[jj, ii] = fscanf(fid, '%g', 1)
                f.bsupumns[jj, ii] = fscanf(fid, '%g', 1)
                f.bsupvmns[jj, ii] = fscanf(fid, '%g', 1)
            # end
        # end
    # end

    f.mnyq = _np.max(f.xm_nyq)
    f.nnyq = _np.max(f.xn_nyq)/f.nfp

    # Calculate the Currents
    f.currumnc = _np.zeros( (f.mnmax_nyq,f.ns), dtype=_np.float64)
    f.currvmnc = _np.zeros_like(f.currumnc)

    for ii in range(1,f.ns-1): # i=2:f.ns-1
        f.currumnc[:, ii] = -f.xn_nyq.flatten()*f.bsubsmns[:, ii] - (f.bsubvmnc[:, ii+1] - f.bsubvmnc[:, ii-1])/(f.ns-1)
        f.currvmnc[:, ii] = -f.xm_nyq.flatten()*f.bsubsmns[:, ii] - (f.bsubumnc[:, ii+1] - f.bsubumnc[:, ii-1])/(f.ns-1)
    # end

    f.currvmnc[f.xm_nyq==0, 0] = 2*f.bsubumnc[ f.xm_nyq == 0, 1]/(f.ns-1)
    f.currumnc[f.xm_nyq==0, 0] = 2*f.currumnc[ f.xm_nyq == 0, 1] - f.currumnc[ f.xm_nyq == 0, 2]
    f.currvmnc[f.xm_nyq!=0, 0] = 0.0
    f.currumnc[f.xm_nyq!=0, 0] = 0.0
    f.currumnc[:, f.ns] = 2*f.currumnc[:, f.ns-1] - f.currumnc[:, f.ns-2]
    f.currvmnc[:, f.ns] = 2*f.currvmnc[:, f.ns-1] - f.currvmnc[:, f.ns-2]
    f.currumnc = f.currumnc/(4*_np.pi*1e-7)
    f.currvmnc = f.currvmnc/(4*_np.pi*1e-7)

    if f.iasym > 0:
        f.currumns = _np.zeros( (f.mnmax_nyq, f.ns), dtype=_np.float64)
        f.currvmns = zeros_like(f.currumns)
        for ii in range(1, f.ns-1):  # i=2:f.ns-1
            f.currumns[:, ii] = -f.xn_nyq.flatten()*f.bsubsmnc[:, ii]-(f.bsubvmns[:, ii+1]-f.bsubvmns[:, ii-1])/(f.ns-1)
            f.currvmns[:, ii] = -f.xm_nyq.flatten()*f.bsubsmnc[:, ii]-(f.bsubumns[:, ii+1]-f.bsubumns[:, ii-1])/(f.ns-1)
        # end

        f.currvmns[f.xm_nyq == 0, 0] = 2*f.bsubumns[f.xm_nyq == 0, 1]/(f.ns-1)
        f.currumns[f.xm_nyq == 0, 0] = 2*f.currumns[f.xm_nyq == 0, 1]-f.currumns[f.xm_nyq == 0, 2]
        f.currvmns[f.xm_nyq != 0, 0] = 0.0
        f.currumns[f.xm_nyq != 0, 0] = 0.0
        f.currumns[:, f.ns] = 2*f.currumns[:, f.ns-1]-f.currumns[:, f.ns-2]
        f.currvmns[:, f.ns] = 2*f.currvmns[:, f.ns-1]-f.currvmns[:, f.ns-2]
        f.currumns = f.currumns/(4*_np.pi*1e-7)
        f.currvmns = f.currvmns/(4*_np.pi*1e-7)
    # end

    # Read the full-mesh quantities
    data = fscanf(fid,fmt6, [6, f.ns])
    f.iotaf = data[0, :]
    f.presf = data[1, :]
    f.phipf = data[2, :]
    f.phi = data[3, :]
    f.jcuru = data[4, :]
    f.jcurv = data[5, :]

    # Read the half-mesh quantities
    data = fscanf(fid,fmt10, [10, f.ns-1])
    f.iotas = data[0, :]
    f.mass = data[1, :]
    f.pres = data[2, :]
    f.beta_vol = data[3, :]
    f.phip = data[4, :]
    f.buco = data[5, :]
    f.bvco = data[6, :]
    f.vp = data[7, :]
    f.overr = data[8, :]
    f.specw = data[9, :]

    data = fscanf(fid, fmt, 6)
    f.aspect = data[0]
    f.betatot = data[1]
    f.betapol = data[2]
    f.betator = data[3]
    f.betaxis = data[4]
    f.b0 = data[5]

    f.isigna = fscanf(fid, fmt+'\n', 1)
    f.input_extension = strtrim(fgetl(fid))

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
    data = fscanf(fid,fmt6, [6, f.ns-2])
    f.Dmerc = data[0, :]
    f.Dshear = data[1, :]
    f.Dwell = data[2, :]
    f.Dcurr = data[3, :]
    f.Dgeod = data[4, :]
    f.equif = data[5, :]

    f.curlabel=cell(f.nextcur, 1)
    if (f.nextcur > 0):
        f.lfreeb = 1
        f.extcur = fscanf(fid, fmt, f.nextcur)
        lcurr = strtrim(fscanf(fid, '%s', 1))
        if lcurr == 'T':  # strcmpi(lcurr,'T'):
            fscanf(fid, '\n')
            rem = f.nextcur
            jj = 0
            while rem > 0:
                line = fgetl(fid)
                fscanf(fid, '\n')
                test = line[0]
                index = line.find(test)  # findstr(line,test)
                for ii in range(_np.size(index, axis=1)/2):  # i=1:size(index,2)/2
                    f.curlabel[ii+jj] = strtrim(line[index[2*ii-1]+range(index[2*ii])-1])
                # end
                jj += _np.size(index, axis=1)/2
                rem -= _np.size(index, axis=1)/2
            # end
        # end
    # end

    data = fscanf(fid, fmt2, [2, f.nstore_seq])
    f.sqt = data[0, :]
    f.wdot = data[1, :]

    data = fscanf(fid, fmt2, [2, f.nstore_seq])
    f.jdotb = data[0, :]
    f.bdotgradv = data[1, :]

    # No Unit Conversion Necessary
    # Data and MSE Fits
    if (f.ireconstruct > 0):
        if ((f.imse >= 2) or (f.itse >0)):
            f.twsgt = fscanf(fid, fmt, 1)
            f.msewgt = fscanf(fid, fmt, 1)
            f.isnodes = fscanf(fid, fmt, 1)

            data = fscanf(fid, fmt3, [3, f.isnodes])
            f.sknots = data[0, :]
            f.ystark = data[1, :]
            f.y2stark = data[2, :]
            f.ipnodes = fscanf(fid, fmt, 1)

            data = fscanf(fid, fmt3, [3, f.ipnodes])
            f.pknots = data[0, :]
            f.ythom = data[1, :]
            f.y2thom = data[2, :]

            data = fscanf(fid,fmt7, [7, (2*f.ns)-1])
            f.anglemse = data[0, :]
            f.rmid = data[1, :]
            f.qmid = data[2, :]
            f.shear = data[3, :]
            f.presmid = data[4, :]
            f.alfa = data[5, :]
            f.curmid = data[6, :]

            data = fscanf(fid, fmt3, [3, f.imse])
            f.rstark = data[0, :]
            f.datastark = data[1, :]
            f.qmeas = data[2, :]

            data = fscanf(fid, fmt2, [2, f.itse])
            f.rthom = data[0, :]
            f.datathom = data[1, :]
        # end

        if (f.nobd > 0):
            data = fscanf(fid, fmt3, [3,f.nobd])
            f.dsiext = data[0, :]
            f.plflux = data[1, :]
            f.dsiobt = data[2, :]
            f.flmwgt = fscanf(fid, fmt, 1)
        # end

        nbfldn=_np.sum(nbldf[:f.nbsets])
        if (nbfldn > 0):
            for nn in range(nbsets):  # n=1:nbsets
                data = fscanf(fid, fmt3, [3, f.nbfld[nn]])
                f.bcoil[:, nn] = data[0, :]
                f.plbfld[:, nn] = data[1, :]
                f.bbc[:, nn] = data[2, :]
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

        f.limitr = fscanf(fid, fmt,  f.nlim)
        f.rlim = _np.zeros( (_np.max(f.limitr), f.nlim), dtype=_np.float64)
        f.zlim = _np.zeros( (_np.max(f.limitr), f.nlim), dtype=_np.float64)
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
    f.mgrid_mode=strtrim(fscanf(fid, '%s'))
    return f
# end  def read_vmec_847

# ---------------------------------------------------------------- %
# ---------------------------------------------------------------- %


def read_vmec_mercier(filname):

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
        line = fgetl(fid)

        while line!='':  # !strcmp(line, ''):
            val = sscanf(line,'%e')
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
        f.dmerc  = val[0]  # - TWICE?
        f.dshear = val[1]
        f.dcurr  = val[2]
        f.dwell  = val[3]
        f.dgeod  = val[4]

        while not feof(fid):
            line = fgetl(fid)
            val  = sscanf(line,'%e')
            f.dmerc  = [f.dmerc,  val[1]]
            f.dshear = [f.dshear, val[2]]
            f.dcurr  = [f.dcurr,  val[3]]
            f.dwell  = [f.dwell,  val[4]]
            f.dgeod  = [f.dgeod,  val[5]]
        # end
    # implicitly close the file using the with open (otherwise explicitly)
    except:
        raise
    finally:
        fid.close()
    return f
# end

# ====================================================================== #

def read_vmec_jxbout(filname):
    try:
        fid = fopen(filname, 'r')
#    with open(filname,'r') as fid:
        # fid = fopen(filname,'r')

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
        f.tor_flux = val[0] #, twice?
        f.fnorm    = val[1]
        f.jdotb    = val[2]
        f.bdotv    = val[3]

        line = fgetl(fid)    # percentages
        fgetl(fid)
        fgetl(fid)
        fgetl(fid)

        f.data = _np.zeros( (f.nrad,f.nzeta,f.ntheta,13), dtype=_np.float64)
        for ii in range(nrad): # i=1:nrad:
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

# ====================================================================== #

def read_vmec_netcdf(filname):
    mu0=4*_np.pi*1e-7

    read_netcdf(filname, 'strip', 'flipdim')

    # Now fix named fields so they match the old way of doing things
    f.input_extension = f.inputextension
    f.mgrid_file = char(f.mgridfile)
    f.ierr_vmec = f.ierflag
    f.rmax_surf = f.rmaxsurf
    f.rmin_surf = f.rminsurf
    f.zmax_surf = f.zmaxsurf
    f.ireconstruct = f.lreconlogical
    f.imse = -1
    f.itse = -1
    f.RBtor = f.rbtor
    f.Rmajor = f.Rmajorp
    f.Aminor = f.Aminorp
    f.betatot = f.betatotal
    f.Volume  =f.volumep
    f.VolAvgB = f.volavgB
    f.beta_vol = f.betavol.T
    f.specw = f.specw.T
    f.iasym = f.lasymlogical
    f.freeb = f.lfreeblogical
    f.lfreeb = f.freeb
    f.Itor = f.ctor
    f.Dmerc = f.DMerc
    f.Dwell = f.DWell
    f.Dshear = f.DShear
    f.Dcurr = f.DCurr
    f.Dgeod = f.DGeod

    # Cast some values
    f.ntor = double(f.ntor)
    f.mpol = double(f.mpol)
    f.nfp = double(f.nfp)
    f.ns = float(f.ns)

    # Calculate Currents
    f.currumnc = _np.zeros((f.mnmaxnyq,f.ns), dtype=_np.float64)
    f.currvmnc = _np.zeros((f.mnmaxnyq,f.ns), dtype=_np.float64)

    for ii in range(1,f.ns-1): # i=2:f.ns-1
        f.currumnc[:, ii] = -double(f.xnnyq).T * f.bsubsmns[:, ii] - (f.ns-1) * (f.bsubvmnc[:, ii+1] - f.bsubvmnc[:, ii] )
        f.currvmnc[:, ii] = -double(f.xmnyq).T * f.bsubsmns[:, ii] + (f.ns-1) * (f.bsubumnc[:, ii+1] - f.bsubumnc[:, ii] )
    # end

    f.currumnc[:, 0] = 0.0
    f.currvmnc[:, 0] = 0.0
    for ii in range(f.mnmaxnyq): # i=1:f.mnmaxnyq
        if ( f.xmnyq[ii]==0 ):
            f.currumnc[ii,0] = 2*f.currumnc[ii,1] - f.currumnc[ii,2]
            f.currvmnc[ii,0] = 2*(f.ns-1)*f.bsubumnc[ii,1]
        # end
    # end

    f.currumnc[:, f.ns] = 2*f.currumnc[:, f.ns-1] - f.currumnc[:, f.ns-2]
    f.currvmnc[:, f.ns] = 2*f.currvmnc[:, f.ns-1] - f.currvmnc[:, f.ns-2]
    f.currumnc = f.currumnc/mu0
    f.currvmnc = f.currvmnc/mu0
    if f.iasym:
        f.currumns=_np.zeros( (f.mnmaxnyq,f.ns), dtype=_np.float64)
        f.currvmns=_np.zeros( (f.mnmaxnyq,f.ns), dtype=_np.float64)
        for ii in range(1,f.ns-1): # i=2:f.ns-1
            f.currumns[:, ii] = -double(f.xnnyq).T*f.bsubsmnc[:, ii] - (f.ns-1)*(f.bsubvmns[:, ii+1] - f.bsubvmns[:, ii])
            f.currvmns[:, ii] = -double(f.xmnyq).T*f.bsubsmnc[:, ii] + (f.ns-1)*(f.bsubumns[:, ii+1] - f.bsubumns[:, ii])
        # end
        f.currumns[:, 0] = 0.0
        f.currvmns[:, 0] = 0.0

        for ii in range(f.mnmaxnyq): # i=1:f.mnmaxnyq
            if (f.xmnyq[ii]==0):
                f.currumns[ii,0] = 2*f.currumns[ii,1] - f.currumns[ii,2]
                f.currvmns[ii,0] = 2*(f.ns-1)*f.bsubumns[ii,1]
            # end
        # end
        f.currumns[:, f.ns] = 2*f.currumns[:, f.ns-1] - f.currumns[:, f.ns-2]
        f.currvmns[:, f.ns] = 2*f.currvmns[:, f.ns-1] - f.currvmns[:, f.ns-2]
        f.currumns = f.currumns/mu0
        f.currvmns = f.currvmns/mu0
    # end

    # Remove Renamed Fields
    del f.inputextension # rmfield(self, 'inputextension')
    del f.mgridfile # rmfield(self, 'mgridfile')
    del f.ierflag # rmfield(self, 'ierflag')
    del f.Rmajorp # rmfield(self, 'Rmajorp')
    del f.Aminorp # rmfield(self, 'Aminorp')
    del f.betatotal # rmfield(self, 'betatotal')
    del f.volumep # rmfield(self, 'volumep')
    del f.volavgB # rmfield(self, 'volavgB')
    del f.betavol # rmfield(self, 'betavol')
    del f.lasymlogical # rmfield(self, 'lasymlogical')
    del f.lfreeblogical # rmfield(self, 'lfreeblogical')
    del f.lreconlogical # rmfield(self, 'lreconlogical')
#    f.close()
    return f
# end  def read_netcdf()






