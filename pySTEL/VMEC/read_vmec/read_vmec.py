# -*- coding: utf-8 -*-
"""
    READ_VMEC(filname) This class reads the VMEC wout file.
    This function reads the wout file and returns the data from the file
    in a structure.  Based on the 'readw_only_priv.f' subroutine.
    In addtion to the raw variables, the structure also contains the
    reconstituted Fourier arrays in their full form.
    Currently this function can read VMEC files up to version 8+ (netCDF and
    text). It returns fourier harmonics in the NESCOIL (nu+nv) format.

    Example usage
        data = read_vmec('wout.test')      # Reads VMEC wout file
        mdata = read_vmec('mercier.test')  # Reads VMEC mercier file

    Original version written in MATLAB
        Maintained by: Samuel Lazerson (lazerson@pppl.gov)
        Version:       1.96
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
    1/30/13     Added calculation of chip
    1/10/14     Modified to read 8.51 files
                Uses new methods for calculating J taking into acount the
                odd modes properly.
    3/31/14     Fixed calculation of non-stellarator symmetric J terms.
    3/1/16      Corrected variable names so text files use the new method
                calculation of current densities.

    Original branch point for porting to python:
        03/03/2016     G.M. Weir (IPP-Greifswald) Ported to Python (old version)

    Updated merge in port to python:
        07/12/2019     G.M. Weir updated for general use
"""
# ======================================================================== #
# ======================================================================== #

from __future__ import absolute_import, with_statement, absolute_import, \
                       division, print_function, unicode_literals

import numpy as _np
from utils.vmcutils import finish_import, h2f, h2f_special

try:
#if 1:
    from utils import Struct, ftell, fgetl
except:
    from ..utils import Struct, ftell, fgetl
# end try

try:
    from libstell import read_vmec as read_vmec_libstell
    use_lib_stell = True
except:
    use_lib_stell = False
# end try

__metaclass__ = type

# ======================================================================== #
# ======================================================================== #


class read_vmec(Struct):
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
    """
    # ===================================================================== #
    # ===================================================================== #
    # Defaults
    netcdffile=0
    version=999999.0

    def __init__(self, filname=None, verbose=True):
        self.filname = filname
        # Number of Arguments
        if filname is None:
            if verbose:
                print('read_vmec requires a filname.')
            # end if
            return None
        else:
            if verbose:
                print('Opening %s for reading'%(filname,))
            # end if
        # end if
        success = self._filexists_(verbose=verbose)
        if success:
            self.readfil(verbose=verbose)
        #endif
    #end def __init__

    def _filexists_(self, verbose=True):
        filname = self.filname
        # Check to see the file exists
        # Note that "with open" only works in certain python versions (3+)
        # try, except, finally works in all tested python versions
        success = 0
        try:
            self.fid = open(filname,'r')
            success = 1
        except:
            raise
        finally:
            try: self.fid.close()
            except:  pass
        #end initialization (automatically closes file using "with open," if not already closed)
        return success
    # end def _filexists_

    def readfil(self, verbose=True):
        filname = self.filname

        # Handle File type and library access
        # nfilname = len(filname)
        netcdffile = 0
        if use_lib_stell and filname.find('wout')>-1 or filname.find('.nc')>-1:
            if verbose:
                print('Attempting to open %s using the libstell library'%(filname,))
            # end if
            self.reader = read_vmec_libstell

            self.__dict__.update(self.reader(filname))

        elif filname.find('.nc')>-1 and filname.find('wout')>-1:
            if verbose:     print('netcdf wout file')   # end if
            netcdffile = 1
            try:
                from VMEC.read_vmec.read_vmec_netCDF import read_vmec as read_vmec_nc
            except:
                from .read_vmec_netCDF import read_vmec as read_vmec_nc
            # end try
            self.reader = read_vmec_nc

            self.__dict__.update(self.reader(filname))

        elif filname.find('mercier')>-1:
            if verbose:     print('mercier file')    # end if
            self.reader = self.read_vmec_mercier

            self.__dict__.update(self.reader(filname))
            return
        elif filname.find('jxbout')>-1:
            if verbose:     print('jxbout file')     # end if
            self.reader = self.read_vmec_jxbout

            self.__dict__.update(self.reader(filname))
        else:
            # Assume a VMEC output file in text format
            if verbose:     print('text-based wout file')   # end if
#            with open(filname,'r') as fid:
#            self.reader = self.read_vmec_txt

#            try:
            if 1:
#                from VMEC.read_vmec.read_vmec_txt_matlabVMEC import read_vmec
                from VMEC.read_vmec.read_vmec_txt import read_vmec
#            except:
#                from .read_vmec_txt import read_wout_txt
            # end try
            self.reader = read_vmec
            self.__dict__.update(self.reader(filname))

#            try:
#                self.fid = open(filname,'r')  # Open File
#                # Read First Line and extract version information
#                line = fgetl(self.fid)
#                # line = self.fid.readline()
#                _, self.version = line.split('=')
#                self.version = float(self.version)
#                # Handle unknown versions
#                print(self.version)
#                if (self.version < 0) or (self.version > 8.52):
#                    if verbose:
#                        print('Unknown file type or version.')
#                    # end if
#                    return
#                # end if
#
#                # VMEC files with comma delimited values do exist, in an attempt to
#                # handle them we dynamically create the format specifier for fscanself.
#                start = ftell(self.fid)   # Get the current position to rewind to
#                line = fgetl(self.fid)   # Get the next line
##                start = self.fid.tell()
##                line = self.fid.readline().replace('\n', '')
#                fmt = '%g'
#                if line.find(',') > -1:  # isempty(strfind(line,',')):
#                    fmt += ','
#                    if verbose and self.version < 6.54:
#                        print('Comma delimited file detected!')
#                        print('     Only 6.54+ version supported!')
#                    # end if
#                    return self
#                # end if
#
#                # Go back to just after the version.
#                self.fid.seek(start, 0)
#
#                # Handle versions
#                self.read_vmec_txt(self.fid, fmt)
#
#                # If VMEC threw an error then return what was read
#                if self.ierr_vmec and (self.ierr_vmec != 4):
#                    print('VMEC runtime error detected!')
#                    return self
#                #end
#            except:
#                raise
#            finally:
#                # (automatically closes file using "with open," if not already closed)
#                try:    self.fid.close()
#                except: pass
##                fclose(fid)
#            # end try
#        # end if filename selection
#
        if not netcdffile:
            # This is taken care of in Jonathon's code I believe.
            # TODO:! Double check this!
            if not hasattr(self, 'lrfplogical'):
                self.lrfplogical=0
            # end if

            # === Now comes the geometric manipulations of the VMEC output === #
            #
            # Need to convert various quantities to full mesh
            self.half2fullmesh()

            # Now recompose the Fourier arrays
            self.RecomposeFourierArrays()

            # Create resonance array
            self. CreateResonanceArray()
        else:
            self.finish_import()
        # end if
    # end def readfil

#    def read_vmec_netcdf(self, filname):
#        try:
#            from VMEC.read_vmec_netCDF import read_vmec
#        except:
#            from .read_vmec_netCDF import read_vmec
#        # end try
#        self.reader = read_vmec
#        f = self.reader(filname)
#
#        # Call the initialization function of the super class to store the data
#        if type(f) != dict:   f = f.dict_from_class()   # endif
#        self.__dict__.update(f)
#        return f
#    # end def read_vmec_netcdf
#
#    def read_vmec_txt(self, fid, fmt):
#        try:
#            from VMEC import read_vmec_txt as rvt
##            import read_vmec_txt as rvt
#        except:
#            import VMEC
#            from .VMEC import read_vmec_txt as rvt
#        # end try
#        if (self.version <= 5.10):
#            self.reader = rvt.read_vmec_orig
#        elif (self.version <= 6.05):
#            self.reader = rvt.read_vmec_605
#        elif (self.version <= 6.20):
#            self.reader = rvt.read_vmec_620
#        elif (self.version <= 6.50):
#            self.reader = rvt.read_vmec_650
#        elif (self.version <= 6.95):
#            self.reader = rvt.read_vmec_695
#        elif (self.version <= 8.00):
#            self.reader = rvt.read_vmec_800
#        elif (self.version <= 8.52):
#            self.reader = rvt.read_vmec_847
#        # end if
#        f = self.reader(fid, fmt)
#
#        try:    fid.close()
#        except: pass
#
#        # Call the initialization function of the super class to store the data
#        if type(f) != dict:   f = f.dict_from_class()   # endif
#        self.__dict__.update(f)
#        return
#    # end def read_vmec_txt

    # ===================================================================== #


    def RecomposeFourierArrays(self):
        # Now recompose the Fourier arrays
        msize = _np.max( self.xm )
        self.nfp = _np.float64(self.nfp)
        nsize = _np.max( self.xn / self.nfp ) - _np.min( self.xn / self.nfp )+1

        msize = int(msize)
        nsize = int(nsize)
        self.rbc = _np.zeros( (msize+1, nsize, self.ns), dtype=_np.float64)
        self.zbs = _np.zeros_like( self.rbc )
        self.lbs = _np.zeros_like( self.rbc )

        # Create derivative terms
        # B_R   = B^u*(dR/du) + B^v (dR/dv)
        # B_phi = R*B^v
        # B_Z   = B^u*(dZ/du) + B^v (dZ/dv)
        self.rsc = _np.zeros_like( self.rbc )  # dRmn/ds*cos
        self.rus = _np.zeros_like( self.rbc )  # dRmn/du*sin
        self.rvs = _np.zeros_like( self.rbc )  # dRmn/dv*sin
        self.zss = _np.zeros_like( self.rbc )  # dZmn/ds*sin
        self.zuc = _np.zeros_like( self.rbc )  # dZmn/du*cos
        self.zvc = _np.zeros_like( self.rbc )  # dZmn/dv*cos

        offset = _np.min(self.xn/self.nfp)-1

        # This is a strange line.  #TODO!: is this VMEC Left-hand coord system, or convenience for Sam?
        self.xn *= -1
        for ii in range(self.ns): #ii=1:self.ns
            for jj in range(self.mnmax): #jj=1:self.mnmax
                mm = self.xm[jj]+1
                nn = -offset+self.xn[jj] / self.nfp

                mm = int(mm)-1
                nn = int(nn)-1
                self.rbc[mm, nn, ii] = self.rmnc[jj, ii]
                self.rus[mm, nn, ii] =-self.rmnc[jj, ii]*self.xm[jj]
                self.rvs[mm, nn, ii] =-self.rmnc[jj, ii]*self.xn[jj]
                self.zbs[mm, nn, ii] = self.zmns[jj, ii]
                self.zuc[mm, nn, ii] = self.zmns[jj, ii]*self.xm[jj]
                self.zvc[mm, nn, ii] = self.zmns[jj, ii]*self.xn[jj]
                self.lbs[mm, nn, ii] = self.lmns[jj, ii]
            # end for
        # end for

        # repmat( a, m, n) -> tile( a, (m,n))
        self.rumns = -self.rmnc*_np.tile( _np.atleast_2d(self.xm), (1,self.ns) )  # repmat( self.xm.T, [1 self.ns] )
        self.rvmns = -self.rmnc*_np.tile( _np.atleast_2d(self.xn), (1,self.ns) )  # repmat( self.xn.T, [1 self.ns] )
        self.zumnc =  self.zmns*_np.tile( _np.atleast_2d(self.xm), (1,self.ns) )  # repmat( self.xm.T, [1 self.ns] )
        self.zvmnc =  self.zmns*_np.tile( _np.atleast_2d(self.xn), (1,self.ns) )  # repmat( self.xn.T, [1 self.ns] )

        # Handle Radial Derivatives
        self.rsc = _np.copy(self.rbc)
        self.zss = _np.copy(self.zbs)
        self.rsmnc = _np.copy(self.rmnc)
        self.zsmns = _np.copy(self.zmns)

        # inner points: centered finite differences
        for ii in range(1,self.ns-1): #i=2:self.ns-1
            self.rsmnc[:, ii] = self.rmnc[:, ii+1] - self.rmnc[:, ii-1]
            self.zsmns[:, ii] = self.zmns[:, ii+1] - self.zmns[:, ii-1]
            self.rsc[:, :, ii] = self.rbc[:, :, ii+1] - self.rbc[:, :, ii-1]
            self.zss[:, :, ii] = self.zbs[:, :, ii+1] - self.zbs[:, :, ii-1]
        # end for inner mesh
        self.rsc *= 0.5
        self.zss *= 0.5
        self.rsmnc *= 0.5
        self.zsmns *= 0.5

        # end-points: forward difference in core
        self.rsc[:, :, 0] = self.rbc[:, :, 1] - 2*self.rbc[:, :, 0]
        self.zss[:, :, 0] = self.zbs[:, :, 1] - 2*self.zbs[:, :, 0]
        self.rsmnc[:, 0] = self.rsmnc[:, 1] - 2*self.rsmnc[:, 0]
        self.zsmns[:, 0] = self.zsmns[:, 1] - 2*self.zsmns[:, 0]

        # end-points: backward difference in edge
        self.rsc[:, :, -1] = 2*self.rbc[:, :, -1] - self.rbc[:, :, -2]
        self.zss[:, :, -1] = 2*self.zbs[:, :, -1] - self.zbs[:, :, -2]
        self.rsmnc[:, -1] = 2*self.rsmnc[:, -1] - self.rsmnc[:, -2]
        self.zsmns[:, -1] = 2*self.zsmns[:, -1] - self.zsmns[:, -2]

        # Handle Vector values seperately to deal with nyqyist
        if hasattr(self,'xm_nyq'):
            msize = _np.max( self.xm_nyq )
            nsize = _np.max( self.xn_nyq / self.nfp ) - _np.min( self.xn_nyq / self.nfp )+1

            mnmax = self.mnmax_nyq
            offset = _np.min( self.xn_nyq / self.nfp ) - 1

            self.xn_nyq *= -1
            xn = self.xn_nyq
            xm = self.xm_nyq
            self.mpol = _np.max( self.xm_nyq )              # Do this for VMECplot
            self.ntor = _np.max( self.xn_nyq / self.nfp )   # Do this for VMECplot
        else:
            msize = _np.max( self.xm )
            nsize = _np.max( self.xn / self.nfp ) - _np.min( self.xn / self.nfp )+1
            mnmax = self.mnmax
            offset = _np.min( self.xn / self.nfp )-1
            xn = self.xn
            xm = self.xm
        # end if hasattr()

        # pre-allocate vectors
        msize = int(msize)
        nsize = int(nsize)
        mnmax = int(mnmax)
        self.bc = _np.zeros( (msize+1, nsize, self.ns), dtype=_np.float64)
        self.gc = _np.zeros_like( self.bc )
        self.b_ss = _np.zeros_like( self.bc )
        self.b_uc = _np.zeros_like( self.bc )
        self.b_vc = _np.zeros_like( self.bc )
        self.buc = _np.zeros_like( self.bc )
        self.bvc = _np.zeros_like( self.bc )
        self.currvc = _np.zeros_like( self.bc )
        if (self.version > 8.0):
            self.curruc = _np.zeros( (msize+1,nsize,self.ns), dtype=_np.float64)
        # end if

        for ii in range(self.ns): #ii=1:self.ns
            for jj in range(mnmax): #jj=1:mnmax
                mm =  xm[jj]+1
                nn = -offset+xn[jj] / self.nfp

                mm = int(mm)-1
                nn = int(nn)-1
                self.bc[mm, nn, ii] = self.bmnc[jj, ii]
                self.gc[mm, nn, ii] = self.gmnc[jj, ii]
                self.b_ss[mm, nn, ii] = self.bsubsmns[jj, ii]
                self.b_uc[mm, nn, ii] = self.bsubumnc[jj, ii]
                self.b_vc[mm, nn, ii] = self.bsubvmnc[jj, ii]
                self.buc[mm, nn, ii] = self.bsupumnc[jj, ii]
                self.bvc[mm, nn, ii] = self.bsupvmnc[jj, ii]
                if hasattr(self, 'currvmnc'):
                    self.currvc[mm, nn, ii] = self.currvmnc[jj, ii]
                # end if
#                self.currvc[mm, nn, ii] = self.currvmnc[jj, ii]
                if self.version > 8.0:
                    self.curruc[mm, nn, ii] = self.currumnc[jj, ii]
                #end if
            #end for jj
        #end for ii

        # Note derivative terms not implemented for non-axisymmetric runs
        # Handle non-axisymmetric runs
        if self.iasym==1:
            msize = _np.max( self.xm )
            nsize = _np.max( self.xn / self.nfp ) - _np.min( self.xn / self.nfp )+1
            msize = int(msize)
            nsize = int(nsize)
            self.rbs = _np.zeros( (msize+1,nsize,self.ns), dtype=_np.float64)
            self.zbc = _np.zeros_like( self.rbs )
            self.lbc = _np.zeros_like( self.rbs )

            # Create derivative terms
            # B_R   = B^u*(dR/du) + B^v (dR/dv)
            # B_phi = R*B^v
            # B_Z   = B^u*(dZ/du) + B^v (dZ/dv)
            self.rss = _np.zeros_like( self.rbs )  # dRmn/ds*cos
            self.ruc = _np.zeros_like( self.rbs )  # dRmn/du*sin
            self.rvc = _np.zeros_like( self.rbs )  # dRmn/dv*sin
            self.zsc = _np.zeros_like( self.rbs )  # dZmn/ds*sin
            self.zus = _np.zeros_like( self.rbs )  # dZmn/du*cos
            self.zvs = _np.zeros_like( self.rbs )  # dZmn/dv*cos
            offset   = _np.min( self.xn / self.nfp )-1

            for ii in range(self.ns): #i=1:self.ns:
                for jj in range(int(self.mnmax)): #j=1:self.mnmax:
                    mm =  self.xm[jj]+1
                    nn = -offset+self.xn[jj] / self.nfp

                    mm = int(mm)-1
                    nn = int(nn)-1
                    self.rbs[mm, nn, ii] = self.rmns[jj, ii]
                    self.ruc[mm, nn, ii] =-self.rmns[jj, ii] * self.xm[jj]
                    self.rvc[mm, nn, ii] =-self.rmns[jj, ii] * self.xn[jj]
                    self.zbc[mm, nn, ii] = self.zmnc[jj, ii]
                    self.zus[mm, nn, ii] = self.zmnc[jj, ii] * self.xm[jj]
                    self.zvs[mm, nn, ii] = self.zmnc[jj, ii] * self.xn[jj]
                    self.lbc[mm, nn, ii] = self.lmnc[jj, ii]
                # end for jj
            # end for ii
            self.rumnc = -self.rmns * _np.tile( _np.atleast_2d(self.xm).T, (1, self.ns) )  # repmat(self.xm.T,[1 self.ns])
            self.rvmnc = -self.rmns * _np.tile( _np.atleast_2d(self.xn).T, (1, self.ns) )  # repmat(self.xn.T,[1 self.ns])
            self.zumns =  self.zmnc * _np.tile( _np.atleast_2d(self.xm).T, (1, self.ns) )  # repmat(self.xm.T,[1 self.ns])
            self.zvmns =  self.zmnc * _np.tile( _np.atleast_2d(self.xn).T, (1, self.ns) )  # repmat(self.xn.T,[1 self.ns])

            # Handle Radial Derivatives - centered finite differencing in plasma
            self.rss = _np.copy(self.rbs)
            self.zsc = _np.copy(self.zbc)
            self.rsmns =_np.copy(self.rmns)
            self.zsmnc = _np.copy(self.zmnc)
            for ii in range(1,self.ns-1): #i=2:self.ns-1:
                self.rsmns[:, ii] = self.rmns[:, ii+1] - self.rmns[:, ii-1]
                self.zsmnc[:, ii] = self.zmnc[:, ii+1] - self.zmnc[:, ii-1]
                self.rss[:, :, ii] = self.rbs[:, :, ii+1] - self.rbs[:, :, ii-1]
                self.zsc[:, :, ii] = self.zbc[:, :, ii+1] - self.zbc[:, :, ii-1]
            # end for
            self.rss *= 0.5
            self.zsc *= 0.5
            self.rsmns *= 0.5
            self.zsmnc *= 0.5

            # forward finite difference at core
            self.rss[:, :, 0] = self.rbs[:, :, 1] - 2*self.rbs[:, :, 0]
            self.zsc[:, :, 0] = self.zbc[:, :, 1] - 2*self.zbc[:, :, 0]
            self.rsmns[:, 0] = self.rsmns[:, 1] - 2*self.rsmns[:, 0]
            self.zsmnc[:, 0] = self.zsmnc[:, 1] - 2*self.zsmnc[:, 0]
            # backward finite difference at edge
            self.rss[:, :, -1] = 2*self.rbs[:, :, -1] - self.rbs[:, :, -2]
            self.zsc[:, :, -1] = 2*self.zbc[:, :, -1] - self.zbc[:, :, -2]
            self.rsmns[:, -1] = 2*self.rsmns[:, -1] - self.rsmns[:, -2]
            self.zsmnc[:, -1] = 2*self.zsmnc[:, -1] - self.zsmnc[:, -2]

            # Handle Vector values seperately to deal with nyqyist
            if hasattr(self,'xm_nyq'):
                msize = _np.max( self.xm_nyq )
                nsize = _np.max( self.xn_nyq / self.nfp ) - _np.min( self.xn_nyq / self.nfp ) + 1
                mnmax = self.mnmax_nyq
                offset= _np.min( self.xn_nyq / self.nfp ) - 1
                xn = self.xn_nyq
                xm = self.xm_nyq
            else:
                msize = _np.max( self.xm )
                nsize = _np.max( self.xn / self.nfp ) - _np.min( self.xn / self.nfp )+1
                mnmax = self.mnmax
                offset = _np.min( self.xn /self.nfp )-1
                xn = self.xn
                xm = self.xm
            # end if

            msize = int(msize)
            nsize = int(nsize)
            mnmax = int(mnmax)
            self.bs = _np.zeros((msize+1,nsize,self.ns), dtype=_np.float64)
            self.gs = _np.zeros_like( self.bs )
            self.b_sc = _np.zeros_like( self.bs )
            self.b_us = _np.zeros_like( self.bs )
            self.b_vs = _np.zeros_like( self.bs )
            self.bus = _np.zeros_like( self.bs )
            self.bvs = _np.zeros_like( self.bs )
            self.currvs = _np.zeros_like( self.bs )
            if self.version > 8.0:
                self.currus=_np.zeros( (msize+1, nsize, self.ns), dtype=_np.float64)
            # end if

            for ii in range(self.ns): #i=1:self.ns:
                for jj in range(mnmax): #j=1:mnmax:
                    mm =  xm[jj]+1
                    nn = -offset+xn[jj] / self.nfp

                    mm = int(mm)-1
                    nn = int(nn)-1
                    self.bs[mm, nn, ii] = self.bmns[jj, ii]
                    self.gs[mm, nn, ii] = self.gmns[jj, ii]
                    self.b_sc[mm, nn, ii] = self.bsubsmnc[jj, ii]
                    self.b_us[mm, nn, ii] = self.bsubumns[jj, ii]
                    self.b_vs[mm, nn, ii] = self.bsubvmns[jj, ii]
                    self.bus[mm, nn, ii] = self.bsupumns[jj, ii]
                    self.bvs[mm, nn, ii] =   self.bsupvmns[jj, ii]
                    self.currvs[mm, nn, ii] = self.currvmns[jj, ii]   # check hasattr?
                    if self.version > 8.0:
                        self.currus[mm, nn, ii] = self.currumnc[jj, ii]
                    # end if
                # end for jj
            # end for ii
        # end is axisymmetric

        # Calculate chi and chip
        if self.lrfplogical:
            pass
            # self.chipf = _np.ones((1,self.ns), dtype=_np.float64)
        else:
            self.chipf = self.iotaf*self.phipf
            self.chif = _np.cumsum(self.chipf)
        # end if
        return self
    # end def RecomposeFourierArrays(self)

    # ============================ #

    def CreateResonanceArray(self):
        # Create the Resonance Array
        # self.M = 1:_np.max( self.xm )
        # self.N = _np.min( self.xn( self.xn > 0 ) ) : _np.min( self.xn( self.xn > 0 ) ) : _np.max( self.xn )

        self.N = _np.arange(_np.min(self.xn[self.xn>0]), _np.max( self.xn ), step=_np.min(self.xn[_np.where(self.xn>0)]))
        self.M = _np.linspace(1, _np.max(self.xm), num=len(self.N), endpoint=True)
#        self.M = _np.arange(1, _np.max(self.xm))
#        self.N = _np.arange(_np.min(self.xn[_np.where(self.xn>0)]), _np.max(self.xn), _np.min(self.xn[_np.where(self.xn>0)]))
#        self.M = _np.asarray(range(1, int(_np.max(self.xm))))
#        self.N = _np.asarray(range(int(_np.min( self.xn[_np.where(self.xn>0)])), int(_np.max(self.xn)), int(_np.min(self.xn[_np.where(self.xn>0)]))))

        Minv = 1.0/self.M
        self.iota_res = self.N.T * Minv

        # iota_spline = spline( )
        # iota_spline = pchip( 1:self.ns, _np.abs(self.iotaf) )
        # iota_spline = pchip( _np.linspace(0, self.ns, self.ns), _np.abs(self.iotaf) )
        # self.s_res  = 0.0*self.iota_res
        #for i=1:size(self.s_res,1)
        #    for j=1:size(self.s_res,2)
        #        if (self.iota_res(i,j) > _np.min(_np.abs(self.iotaf))) && (self.iota_res(i,j) < _np.max(_np.abs(self.iotaf)))
        #            f_iota=@(x)ppval(iota_spline,x)-self.iota_res(i,j)
        #            self.s_res(i,j)=fzero(f_iota,50)
        #        #end
        #    #end
        #end

        # ================0 #
        # Fix the flow variables
        if hasattr(self,'protmnc'):
            self.protmnc = self.protmnc / (_np.pi*4.0E-7)
        if hasattr(self,'protrsqmnc'):
            self.protrsqmnc = self.protrsqmnc / (_np.pi*4.0E-7)
        if hasattr(self,'prprmnc'):
            self.prprmnc = self.prprmnc / (_np.pi*4.0E-7)
        if hasattr(self,'protmns'):
            self.protmns = self.protmns / (_np.pi*4.0E-7)
        if hasattr(self,'protrsqmns'):
            self.protrsqmns = self.protrsqmns / (_np.pi*4.0E-7)
        if hasattr(self,'prprmns'):
            self.prprmns = self.prprmns / (_np.pi*4.0E-7)
        # end if

        # ================0 #

        # Calculate the stored energy
        if _np.size(self.vp) == self.ns:
            self.eplasma = 1.5 * _np.pi * _np.pi * _np.sum( self.vp * self.presf ) / self.ns
        else:
            self.eplasma = 1.5 * _np.pi * _np.pi * _np.sum( self.vp * self.pres ) / (self.ns-1)
        # end if

        # Set some defaults
        self.nu = 2*( self.mpol+1 )+6
        self.datatype = 'wout'
        return self
    #end def CreateResonanceArray()

    # ================================================================ #
    # ================================================================ #

    def finish_import(self):
        self.__dict__.update(finish_import(self.dict_from_class()))

    def half2fullmesh(self):
        # Interpolate various quantities to the full mesh
        # The following quantities are on 1D half mesh
        # iotas, mass, pres, beta_vol, buco, bvco, vp, specw, phip, jdotb,
        # bdotgradv
        # The following quantities are on 1D full mesh
        # iotaf, presf, phi, phipf, jcuru, jcurv
        # The following quantities are on 1D full mesh from 1:-1 (MATLAB: 2:ns-1)
        # Dmerc, Dshear, Dwell, Dcurr, Dgeod

        # Now we create full mesh quantities for those that have full mesh names
        if not hasattr(self, 'iotaf'):  # !isfield(self, 'iotaf'):
            self.iotaf = h2f( self.iotas, self.ns)
        if not hasattr(self, 'presf'):  # !isfield(self, 'presf'):
            self.presf = h2f( self.pres, self.ns)
        if not hasattr(self, 'phipf'):  # !isfield(self, 'phipf'):
            self.phipf = h2f(self.phip, self.ns)
        # end if

        # Now the rest of the quantities are remapped to full mesh
        if hasattr(self, 'beta_vol'):
            self.beta_vol = h2f(self.beta_vol, self.ns)
        # end if
        self.buco = h2f(self.buco, self.ns)
        self.bvco = h2f(self.bvco, self.ns)
        self.vp = h2f_special(self.vp, len(self.vp))
        self.overr = h2f(self.overr, self.ns)
        self.specw = h2f(self.specw, self.ns)
        if len(self.jdotb) ==self.ns:
            self.jdotb = h2f(self.jdotb, self.ns)
            self.bdotgradv = h2f(self.bdotgradv, self.ns)
        # end if

        # Check to make sure phi,jcuru and jcurv are the same size as self.ns
        if len(self.phi) != self.ns:
#            temp = _np.zeros( (1, self.ns), dtype=_np.float64)
            temp = _np.zeros( (self.ns,), dtype=_np.float64)
#            temp[0:self.ns] = self.phi
            temp[1:] = _np.copy(self.phi)
            temp[0] = 2*temp[1]-temp[2]
            temp[-1] = 2*temp[-2]-temp[-3]
            self.phi = _np.copy(temp)
        # end if
        if len(self.jcuru) != self.ns:
#            temp = _np.zeros( (1, self.ns), dtype=_np.float64)
            temp = _np.zeros( (self.ns,), dtype=_np.float64)
#            temp[0:self.ns] = self.jcuru
            temp[1:] = _np.copy(self.jcuru)
            temp[0]  = 2*temp[1]-temp[2]
            temp[-1] = 2*temp[-2]-temp[-3]
            self.jcuru = _np.copy(temp)
        # end if
        if len(self.jcurv) != self.ns:
#            temp = _np.zeros( (1,self.ns), dtype=_np.float64)
            temp = _np.zeros( (self.ns,), dtype=_np.float64)
            temp[1:] =_np.copy(self.jcurv)
            temp[0] = 2*temp[1]-temp[2]
            temp[-1] = 2*temp[-2]-temp[-3]
            self.jcurv = _np.copy(temp)
        # end if

        # Put the FLOW arrays on the full mesh
        if hasattr(self, 'pmap'):
#            temp = _np.zeros((1,self.ns), dtype=_np.float64)
            temp = _np.zeros( (self.ns,), dtype=_np.float64)
            temp[1:] = _np.copy(self.pmap[:-1])
            temp[0] = 2*temp[1]-temp[2]
            temp[-1] = 2*temp[-2]-temp[-3]
            self.pmap = _np.copy(temp)
        # end if
        if hasattr(self,'omega'):  # Note this is zero on axis so extrapolate first
#            temp = _np.zeros((1,self.ns), dtype=_np.float64)
            temp = _np.zeros( (self.ns,), dtype=_np.float64)
            temp[1:] = _np.copy(self.omega[:-1])
            temp[1] = 2*temp[2]-temp[3]
            temp[0] = 2*temp[1]-temp[2]
            temp[-1] = 2*temp[-2]-temp[-3]
            self.omega = _np.copy(temp)
        # end if
        if hasattr(self,'tpotb'):  # Note this is zero on axis so extrapolate first
#            temp = _np.zeros((1,self.ns), dtype=_np.float64)
            temp = _np.zeros( (self.ns,), dtype=_np.float64)
            temp[1:] = _np.copy(self.tpotb[:-1])
            temp[1] = 2*temp[2]-temp[3]
            temp[0] = 2*temp[1]-temp[2]
            temp[-1] = 2*temp[-2]-temp[-3]
            self.tpotb = _np.copy(temp)
        # end if

        # Fix the Stability values by making into ns arrays
        if hasattr(self, 'Dmerc') and len(self.Dmerc) != self.ns:
#            temp = _np.zeros( (1, self.ns), dtype=_np.float64)
            temp = _np.zeros( (self.ns,), dtype=_np.float64)
            temp[1:-1] = _np.copy(self.Dmerc)
            self.Dmerc = _np.copy(temp)
            self.Dmerc[0] = 2 * self.Dmerc[1] - self.Dmerc[2]
            self.Dmerc[-1] = 2*self.Dmerc[-2] - self.Dmerc[-3]

            temp[1:-1] = _np.copy(self.Dshear)   # 2:self.ns-1
            self.Dshear = _np.copy(temp)
            self.Dshear[0] = 2*self.Dshear[1] - self.Dshear[2]
            self.Dshear[-1] = 2*self.Dshear[-2] - self.Dshear[-3]

            temp[1:-1] = _np.copy(self.Dwell)
            self.Dwell = _np.copy(temp)
            self.Dwell[0] = 2*self.Dwell[1] - self.Dwell[2]
            self.Dwell[-1] = 2*self.Dwell[-2] - self.Dwell[-3]

            temp[1:-1] = _np.copy(self.Dcurr)
            self.Dcurr = _np.copy(temp)
            self.Dcurr[0] = 2*self.Dcurr[1] - self.Dcurr[2]
            self.Dcurr[-1] = 2*self.Dcurr[-2] - self.Dcurr[-3]

            temp[1:-1] = _np.copy(self.Dgeod)
            self.Dgeod = _np.copy(temp)
            self.Dgeod[0] = 2*self.Dgeod[1] - self.Dgeod[2]
            self.Dgeod[-1] = 2*self.Dgeod[-2] - self.Dgeod[-3]
        # end if

        # ===================================================================== #
        # ================== Now do the matrix values ================ #

#        for name in ['lmns', 'bsupumnc', 'bsupvmnc', 'bsubsmns','bsubumnc', 'bsubvmnc', 'gmnc', 'bmnc']:
#            # First Index (note indexing on vectors is 2:ns when VMEC outputs)
#            self.(name)[:, 0] = 1.5*self.(name)[:, 1] - 0.5*self.(name)[:, 2]
#            for ii in range(1,self.ns-1):  #i=2:self.ns-1
#                # Average
#                self.(name)[:, ii] = 0.5*( self.(name)[:, ii] + self.(name)[:, ii+1] )
#            # end for
#            # Last Index (note indexing on vectors is 2:ns when VMEC outputs)
#            self.(name)[:, -1] = 2.0 * self.(name)[:, -2] - self.(name)[:, -3]
#        # end for

        # First Index (note indexing on vectors is 2:ns when VMEC outputs)
        self.lmns[:, 0] = 1.5*self.lmns[:, 1] - 0.5*self.lmns[:, 2]
        self.bsupumnc[:, 0] = 1.5*self.bsupumnc[:, 1] - 0.5*self.bsupumnc[:, 2]
        self.bsupvmnc[:, 0] = 1.5*self.bsupvmnc[:, 1] - 0.5*self.bsupvmnc[:, 2]
        self.bsubsmns[:, 0] = 1.5*self.bsubsmns[:, 1] - 0.5*self.bsubsmns[:, 2]
        self.bsubumnc[:, 0] = 1.5*self.bsubumnc[:, 1] - 0.5*self.bsubumnc[:, 2]
        self.bsubvmnc[:, 0] = 1.5*self.bsubvmnc[:, 1] - 0.5*self.bsubvmnc[:, 2]
        self.gmnc[:, 0] = 1.5*self.gmnc[:, 1] - 0.5*self.gmnc[:, 2]
        self.bmnc[:, 0] = 1.5*self.bmnc[:, 1] - 0.5*self.bmnc[:, 2]

        # Average
        for ii in range(1,self.ns-1):  #i=2:self.ns-1
            self.lmns[:, ii] = 0.5*( self.lmns[:, ii] + self.lmns[:, ii+1] )
            self.bsupumnc[:, ii] = 0.5*( self.bsupumnc[:, ii] + self.bsupumnc[:, ii+1] )
            self.bsupvmnc[:, ii] = 0.5*( self.bsupvmnc[:, ii] + self.bsupvmnc[:, ii+1] )
            self.bsubsmns[:, ii] = 0.5*( self.bsubsmns[:, ii] + self.bsubsmns[:, ii+1] )
            self.bsubumnc[:, ii] = 0.5*( self.bsubumnc[:, ii] + self.bsubumnc[:, ii+1] )
            self.bsubvmnc[:, ii] = 0.5*( self.bsubvmnc[:, ii] + self.bsubvmnc[:, ii+1] )
            self.gmnc[:, ii] = 0.5*( self.gmnc[:, ii] + self.gmnc[:, ii+1] )
            self.bmnc[:, ii] = 0.5*( self.bmnc[:, ii] + self.bmnc[:, ii+1] )
        # end for

        # Last Index (note indexing on vectors is 2:ns when VMEC outputs)
        self.lmns[:, -1]     = 2.0 * self.lmns[:, -2] - self.lmns[:, -3]
        self.bsupumnc[:, -1] = 2.0 * self.bsupumnc[:, -2] - self.bsupumnc[:, -3]
        self.bsupvmnc[:, -1] = 2.0 * self.bsupvmnc[:, -2] - self.bsupvmnc[:, -3]
        self.bsubsmns[:, -1] = 2.0 * self.bsubsmns[:, -2] - self.bsubsmns[:, -3]
        self.bsubumnc[:, -1] = 2.0 * self.bsubumnc[:, -2] - self.bsubumnc[:, -3]
        self.bsubvmnc[:, -1] = 2.0 * self.bsubvmnc[:, -2] - self.bsubvmnc[:, -3]
        self.gmnc[:, -1] = 2.0 * self.gmnc[:, -2] - self.gmnc[:,-3]
        self.bmnc[:,-1] = 2.0 * self.bmnc[:, -2] - self.bmnc[:, -3]

        # Handle ANI/FLOW Values
        if hasattr(self,'prprmnc'):
            self.prprmnc[:,0] = 1.5*self.prprmnc[:,1] - 0.5*self.prprmnc[:,2]
            for ii in range(1, self.ns-1): # i=2:f.ns-1:
                self.prprmnc[:, ii] = 0.5*(self.prprmnc[:, ii] + self.prprmnc[:,ii+1] )
            # end for
            self.prprmnc[:,-1]= 2.0*self.prprmnc[:,-2] - self.prprmnc[:,-3]
        # end if
        if hasattr(self,'protmnc'):
            self.protmnc[:,0] = 1.5*self.protmnc[:,1] - 0.5*self.protmnc[:,2]
            for ii in range(1, self.ns-1):
                self.protmnc[:, ii] = 0.5*( self.protmnc[:, ii] + self.protmnc[:,ii+1] )
            # end for
            self.protmnc[:,-1]= 2.0*self.protmnc[:,-2] - self.protmnc[:,-3]
        # end if
        if hasattr(self,'protrsqmnc'):
            self.protrsqmnc[:,0] = 1.5*self.protrsqmnc[:,1] - 0.5*self.protrsqmnc[:,2]
            for ii in range(1, self.ns-1): #
                self.protrsqmnc[:, ii] = 0.5 * (self.protrsqmnc[:, ii] + self.protrsqmnc[:,ii+1] )
            # end for
            self.protrsqmnc[:,-1]= 2.0*self.protrsqmnc[:,-2] - self.protrsqmnc[:,-3]
        # end if
        if  not hasattr(self,'iasym'):
            self.iasym=0
        # end if

        if self.iasym >0: # Handle existance of lmnc on half mesh
            self.lmnc[:, 0] =     1.5 *     self.lmnc[:, 1] - 0.5 *     self.lmnc[:, 2] #fixed indices
            self.bsupumns[:, 0] = 1.5 * self.bsupumns[:, 1] - 0.5 * self.bsupumns[:, 2]
            self.bsupvmns[:, 0] = 1.5 * self.bsupvmns[:, 1] - 0.5 * self.bsupvmns[:, 2]
            self.bsubsmnc[:, 0] = 1.5 * self.bsubsmnc[:, 1] - 0.5 * self.bsubsmnc[:, 2]
            self.bsubumns[:, 0] = 1.5 * self.bsubumns[:, 1] - 0.5 * self.bsubumns[:, 2]
            self.bsubvmns[:, 0] = 1.5 * self.bsubvmns[:, 1] - 0.5 * self.bsubvmns[:, 2]
            self.gmns[:, 0] =     1.5 *     self.gmns[:, 1] - 0.5 *     self.gmns[:, 2]
            self.bmns[:, 0] =     1.5 *     self.bmns[:, 1] - 0.5 *     self.bmns[:, 2]

            for ii in range(1,self.ns-1): #i=2:self.ns-1
                self.lmnc[:, ii]=     0.5 * (     self.lmnc[:, ii] +     self.lmnc[:, ii+1] )
                self.bsupumns[:, ii]= 0.5 * ( self.bsupumns[:, ii] + self.bsupumns[:, ii+1] )
                self.bsupvmns[:, ii]= 0.5 * ( self.bsupvmns[:, ii] + self.bsupvmns[:, ii+1] )
                self.bsubsmnc[:, ii]= 0.5 * ( self.bsubsmnc[:, ii] + self.bsubsmnc[:, ii+1] )
                self.bsubumns[:, ii]= 0.5 * ( self.bsubumns[:, ii] + self.bsubumns[:, ii+1] )
                self.bsubvmns[:, ii]= 0.5 * ( self.bsubvmns[:, ii] + self.bsubvmns[:, ii+1] )
                self.gmns[:, ii]    = 0.5 * (     self.gmns[:, ii] +     self.gmns[:, ii+1] )
                self.bmns[:, ii]    = 0.5 * (     self.bmns[:, ii] +     self.bmns[:, ii+1] )
            # end

            self.lmnc[:, -1]  = 2.0 *     self.lmnc[:, -2] -     self.lmnc[:, -3]
            self.bsupumns[:, -1] = 2.0 * self.bsupumns[:, -2] - self.bsupumns[:, -3]
            self.bsupvmns[:, -1] = 2.0 * self.bsupvmns[:, -2] - self.bsupvmns[:, -3]
            self.bsubsmnc[:, -1] = 2.0 * self.bsubsmnc[:, -2] - self.bsubsmnc[:, -3]
            self.bsubumns[:, -1] = 2.0 * self.bsubumns[:, -2] - self.bsubumns[:, -3]
            self.bsubvmns[:, -1] = 2.0 * self.bsubvmns[:, -2] - self.bsubvmns[:, -3]
            self.gmns[:, -1] = 2.0 *     self.gmns[:, -2] -     self.gmns[:, -3]
            self.bmns[:, -1] = 2.0 *     self.bmns[:, -2] -     self.bmns[:, -3]
        # end if
        return self
    # end def half2fullmesh

    # ================================================================ #
    # ================================================================ #

    def __delattr__(self, name):
        delattr(self, name)
        #del self.(name)
    # end def __delattr__()
#end class

# --------------------------------------------------------------------- #

#class vmecdso():
#    def __init__(filename):
#        self.fname = filename
#        self.exists = False
#        self.seekoffset = 0
#        self.filelength = 0
#        if _os.path.exists(filename):
#            self.exists=True
#            self.openFile(filename)
#        # endif
#    # end def
#
#    def openFile(self, filename):
#        try:
#            self.fid = open(filename, 'r')
#
#        self.filelength =
#        return self.fid
#    # end def openFile
#
#    def readFile(self, buffsize, data_offset=0):
#       endpos = buffsize + data_offset
#       if self.filelength >= endpos:
#          ReadBuffSize = buffsize
#       else :
#          ReadBuffSize = self.filelength - data_offset
#      # endif
#       data = _np.zeros((ReadBuffSize,2), order='F')
#       self.fid.seek(data_offset+self.seek_offset)
#       if self.COMM_TYPE == 0:
#          bindata = self.fid.read(ReadBuffSize)
#          intd = struct.unpack(str(ReadBuffSize)+'b',bindata[0:ReadBuffSize])
#       else:
#          bindata = self.fid.read(ReadBuffSize*2)
#          intd = struct.unpack(str(ReadBuffSize)+'h',bindata[0:ReadBuffSize*2])
#       for i in range(len(intd)):
#          data[i,1] = intd[i]*self.VERTICAL_GAIN - self.VERTICAL_OFFSET
#          data[i,0] = (i+data_offset)*self.HORIZ_INTERVAL + self.HORIZ_OFFSET
#       return data
#
#    def closeFile(self, fid):
#        self.fid.close()
#    # end def closeFile

# ======================================================================== #
# ======================================================================== #


if __name__=="__main__":
    import time
    start = time.time()

    import os as _os
    rootdir = _os.path.join('d:/', 'Workshop', 'TRAVIS', 'MagnConfigs', 'W7X')
    if not _os.path.exists(rootdir):
        rootdir = _os.path.join('G:/', 'Workshop', 'TRAVIS_tree', 'MagnConfigs', 'W7X')
    filname = []
    filname.append(_os.path.join(rootdir, 'wout_w7x.1000_1000_1000_1000_+0390_+0000.05.0144.txt'))
#    filname.append(_os.path.join(rootdir, 'wout_w7x.1000_1000_1000_1000_+0000_+0000.01.00.txt'))
#    filname.append(_os.path.join(rootdir, 'wout_w7x.1000_1000_1000_1000_+0390_+0000.05.0144.nc'))
#
#    rootdir = _os.path.join('d:/', 'Workshop', 'TRAVIS', 'MagnConfigs', 'HSX', 'QHS')
#    if not _os.path.exists(rootdir):
#        rootdir = _os.path.join('G:/', 'Workshop', 'TRAVIS_tree', 'MagnConfigs', 'HSX', 'QHS')
#    filname.append(_os.path.join(rootdir, 'wout_QHS_Rstart_1_513_AXIS_v2_VMEC_8p46_2012_03_16_39_52.txt'))
#
#    rootdir = _os.path.join('d:/', 'Workshop', 'TRAVIS', 'MagnConfigs', 'HSX', 'QHS')
#    if not _os.path.exists(rootdir):
#        rootdir = _os.path.join('G:/', 'Workshop', 'TRAVIS_tree', 'MagnConfigs', 'HSX', 'QHS', '8p49')
#    filname.append(_os.path.join(rootdir, 'wout_QHS_Rstart_1_513_AXIS_v2_VMEC_8p46.txt'))
#
#    rootdir = _os.path.join('d:/', 'Workshop', 'TRAVIS', 'MagnConfigs', 'HeliotronJ')
#    if not _os.path.exists(rootdir):
#        rootdir = _os.path.join('G:/', 'Workshop', 'TRAVIS_tree', 'MagnConfigs', 'HeliotronJ')
#    filname.append(_os.path.join(rootdir, 'wout.vm8_hj_HV89TA74TB81AV46IV76_x384m12n12s101bd00id+00'))

    for fil in filname:
#        vmecfile = _os.path.join(rootdir, fil)
        start1 = time.time()
        data = read_vmec(fil)

        tst = data.dict_from_class()
        finish_import(tst)

        end1 = time.time()
        print((fil, end1-start1))
    # end for

    endall = time.time()
    print(('%i files'%(len(filname),), endall-start))
# end if


# ======================================================================== #
# ======================================================================== #




