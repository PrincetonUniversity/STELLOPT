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

try:
    from .Struct import Struct
    from .Spline import Spline
except:
    from libstell.read_vmec import Struct
    from libstell.read_vmec import Spline
# end try

#from pypchip import pchip

# ======================================================================== #
# ======================================================================== #


class read_vmec(object):
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

    def __init__(self, filname = None):
        self.filname = filname
        # Number of Arguments
        if filname is None:
            print('read_vmec requires a filname.')
            return None
        # end if

        if self._filexists_():
            self.readfil()
        #endif
    #end def __init__

    def _filexists_(self):
        filname = self.filname
        # Check to see the file exists
        # Note that "with open" only works in certain python versions (3+)
        # try, except, finally works in all tested python versions
        success = 0
        try:
            self.fid = open(filname,'r+')
#        with open(filname,'r+') as fid:
#            # fid=fopen(filname,'r+')
            if (self.fid < 0):
                print( 'ERROR: Could not find file '+filname)
                # print( fid )
                return success
            # end if
            self.fid.close()
            success = 1
        except:
            pass
        finally:
            self.fid.close()
        #end initialization (automatically closes file using "with open," if not already closed)
        return success
    # end def _filexists_

    def readfil(self):
        filname = self.filname
        # Handle File type
        nfilname = len(filname)
        netcdffile = 0
        if (filname.find('.nc')>-1 and not strfind(filname, 'wout')>-1):
            netcdffile = 1  # This is a netCDF file output by VMEC
        elif filname.find('mercier')>-1:
            print('mercier file')
            f = self.read_vmec_mercier(filname)
            return self
        elif filname.find('jxbout')>-1:
            print('jxbout file')
            f = self.read_vmec_jxbout(filname)
            return self
        else:
            # Assume a VMEC output file in text format
#            with open(filname,'r') as fid:
            try:
                self.fid = fopen(filname,'r')  # Open File

                # Read First Line and extract version information
                line = self.fid.readline()
                _, self.version = line.split('=')
                self.version = float(self.version)
                # Handle unknown versions
                if (self.version < 0) or (self.version > 8.47):
                    print('Unknown file type or version.')
                    return self
                #end

                # VMEC files with comma delimited values do exist, in an attempt to
                # handle them we dynamically create the format specifier for fscanself.
#                start = ftell(fid)   # Get the current position to rewind to
#                line = fgetl(fid)   # Get the next line
                start = self.fid.tell()
                line = self.fid.readline()
                if line.find(',') < 0:  # isempty(strfind(line,',')):
                    fmt = '#g'
                else:
                    fmt = '#g,'
                    print('Comma delimited file detected!')
                    print('     Only 6.54+ version supported!')
                    if self.version < 6.54:
                        return self
                    # end
                # end if
                # Go back to just after the version.
                fseek(self.fid, start, 'bof')

                # Handle versions
                #   pre 5.1, 6.05, 6.20, 6.50, 6.54
                if netcdffile:
                    self.read_vmec_netcdf(filname)
                elif (self.version <= 5.10):
                    self.read_vmec_orig(fid, fmt)
                elif (self.version <= 6.05):
                    self.read_vmec_605(fid, fmt)
                elif (self.version <= 6.20):
                    self.read_vmec_620(fid,fmt)
                elif (self.version <= 6.50):
                    self.read_vmec_650(fid,fmt)
                elif (self.version <= 6.95):
                    self.read_vmec_695(fid,fmt)
                elif (self.version <= 8.00):
                    self.read_vmec_800(fid,fmt)
                elif (self.version <= 8.47):
                    self.read_vmec_847(fid,fmt)
                # end if

                if not netcdffile:
                    self.fid.close()
#                    fclose(fid)  # Close the File
                    self.version = version # add version info to structure
                # end

                # If VMEC through an error then return what was read
                if self.ierr_vmec and (self.ierr_vmec != 4):
                    print('VMEC runtime error detected!')
                    return self
                #end
            except:
                raise
            finally:
                # (automatically closes file using "with open," if not already closed)
                self.fid.close()
#                fclose(fid)
            # end try
        # end if filename selection

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
    # end def readfil

    # ===================================================================== #
    #
    # Public functions within the read_vmec class

    def fread1darray(self, fid, dtyp):
        data1D = _np.fromfile(fid, dtype=dtyp)
        return data1D
    #end def fread1Darray

    # ======================== #

    def fread2darray(self, fid, dtyp):
        data2D = _np.fromfile(fid, dtyp).reshape((-1,2)).T
        return data2D
    #end def fread2Darray

    # ======================== #

    def RecomposeFourierArrays(self):
        # Now recompose the Fourier arrays
        msize = _np.max( self.xm )
        nsize = _np.max( self.xn / self.nfp ) - _np.min( self.xn / self.nfp )+1

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
        self.xn = -self.xn
        for ii in range(self.ns): #ii=1:self.ns
            for jj in range(self.mnmax): #jj=1:self.mnmax
                mm = self.xm[jj]+1
                nn = -offset+self.xn[jj] / self.nfp

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
        self.rumns = -self.rmnc*_np.tile( self.xm.T, (1,self.ns) )  # repmat( self.xm.T, [1 self.ns] )
        self.rvmns = -self.rmnc*_np.tile( self.xn.T, (1,self.ns) )  # repmat( self.xn.T, [1 self.ns] )
        self.zumnc =  self.zmns*_np.tile( self.xm.T, (1,self.ns) )  # repmat( self.xm.T, [1 self.ns] )
        self.zvmnc =  self.zmns*_np.tile( self.xn.T, (1,self.ns) )  # repmat( self.xn.T, [1 self.ns] )

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
        self.rsc[:, :, self.ns] = 2*self.rbc[:, :, self.ns] - self.rbc[:, :, self.ns-1]
        self.zss[:, :, self.ns] = 2*self.zbs[:, :, self.ns] - self.zbs[:, :, self.ns-1]
        self.rsmnc[:, self.ns] = 2*self.rsmnc[:, self.ns] - self.rsmnc[:, self.ns-1]
        self.zsmns[:, self.ns] = 2*self.zsmns[:, self.ns] - self.zsmns[:, self.ns-1]

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
        self.bc = _np.zeros( (msize+1, nsize, self.ns), dtype=_np.float64)
        self.gc = _np.zeros_like( self.bc )
        self.b_ss = _np.zeros_like( self.bc )
        self.b_uc = _np.zeros_like( self.bc )
        self.b_vc = _np.zeros_like( self.bc )
        self.buc = _np.zeros_like( self.bc )
        self.bvc = _np.zeros_like( self.bc )
        self.currvc = _np.zeros_like( self.bc )
        if (self.version > 8.0):
            self.curruc = _np.zeros( (msize,nsize,self.ns), dtype=_np.float64)
        # end if

        for ii in range(self.ns): #ii=1:self.ns
            for jj in range(mnmax): #jj=1:mnmax
                mm =  xm[jj]+1
                nn = -offset+xn[jj] / self.nfp

                self.bc[mm, nn, ii] = self.bmnc[jj, ii]
                self.gc[mm, nn, ii] = self.gmnc[jj, ii]
                self.b_ss[mm, nn, ii] = self.bsubsmns[jj, ii]
                self.b_uc[mm, nn, ii] = self.bsubumnc[jj, ii]
                self.b_vc[mm, nn, ii] = self.bsubvmnc[jj, ii]
                self.buc[mm, nn, ii] = self.bsupumnc[jj, ii]
                self.bvc[mm, nn, ii] = self.bsupvmnc[jj, ii]
                self.currvc[mm, nn, ii] = self.currvmnc[jj, ii]
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
                for jj in range(self.mnmax): #j=1:self.mnmax:
                    mm =  self.xm[jj]+1
                    nn = -offset+self.xn[jj] / self.nfp
                    self.rbs[mm, nn, ii] = self.rmns[jj, ii]
                    self.ruc[mm, nn, ii] =-self.rmns[jj, ii] * self.xm[jj]
                    self.rvc[mm, nn, ii] =-self.rmns[jj, ii] * self.xn[jj]
                    self.zbc[mm, nn, ii] = self.zmnc[jj, ii]
                    self.zus[mm, nn, ii] = self.zmnc[jj, ii] * self.xm[jj]
                    self.zvs[mm, nn, ii] = self.zmnc[jj, ii] * self.xn[jj]
                    self.lbc[mm, nn, ii] = self.lmnc[jj, ii]
                # end for jj
            # end for ii
            self.rumnc = -self.rmns * _np.tile( self.xm.T, (1, self.ns) )  # repmat(self.xm.T,[1 self.ns])
            self.rvmnc = -self.rmns * _np.tile( self.xn.T, (1, self.ns) )  # repmat(self.xn.T,[1 self.ns])
            self.zumns =  self.zmnc * _np.tile( self.xm.T, (1, self.ns) )  # repmat(self.xm.T,[1 self.ns])
            self.zvmns =  self.zmnc * _np.tile( self.xn.T, (1, self.ns) )  # repmat(self.xn.T,[1 self.ns])

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
            self.rss[:, :, 0] = self.rbs[:, :, 1] - 2*self.rbs[:, :, 0] #changed indices
            self.zsc[:, :, 0] = self.zbc[:, :, 1] - 2*self.zbc[:, :, 0]
            self.rsmns[:, 0] = self.rsmns[:, 1] - 2*self.rsmns[:, 0]
            self.zsmnc[:, 0] = self.zsmnc[:, 1] - 2*self.zsmnc[:, 0]
            # backward finite difference at edge
            self.rss[:, :, self.ns] = 2*self.rbs[:, :, self.ns] - self.rbs[:, :, self.ns-1]
            self.zsc[:, :, self.ns] = 2*self.zbc[:, :, self.ns] - self.zbc[:, :, self.ns-1]
            self.rsmns[:, self.ns] = 2*self.rsmns[:, self.ns] - self.rsmns[:, self.ns-1]
            self.zsmnc[:, self.ns] = 2*self.zsmnc[:, self.ns] - self.zsmnc[:, self.ns-1]

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

            self.bs = _np.zeros((msize+1,nsize,self.ns), dtype=_np.float64)
            self.gs = _np.zeros_like( self.bs )
            self.b_sc = _np.zeros_like( self.bs )
            self.b_us = _np.zeros_like( self.bs )
            self.b_vs = _np.zeros_like( self.bs )
            self.bus = _np.zeros_like( self.bs )
            self.bvs = _np.zeros_like( self.bs )
            self.currvs = _np.zeros_like( self.bs )
            if self.version > 8.0:
                self.currus=_np.zeros( (msize,nsize,self.ns), dtype=_np.float64)
            # end if

            for ii in range(self.ns): #i=1:self.ns:
                for jj in range(mnmax): #j=1:mnmax:
                    mm =  xm[jj]+1
                    nn = -offset+xn[jj] / self.nfp
                    self.bs[mm, nn, ii] = self.bmns[jj, ii]
                    self.gs[mm, nn, ii] = self.gmns[jj, ii]
                    self.b_sc[mm, nn, ii] = self.bsubsmnc[jj, ii]
                    self.b_us[mm, nn, ii] = self.bsubumns[jj, ii]
                    self.b_vs[mm, nn, ii] = self.bsubvmns[jj, ii]
                    self.bus[mm, nn, ii] = self.bsupumns[jj, ii]
                    self.bvs[mm, nn, ii] =   self.bsupvmns[jj, ii]
                    self.currvs[mm, nn, ii] = self.currvmns[jj, ii]
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
#        self.M = _np.arange(1, _np.max(self.xm))
#        self.N = _np.arange(_np.min(self.xn[_np.where(self.xn>0)]), _np.max(self.xn), _np.min(self.xn[_np.where(self.xn>0)]))
        self.M = _np.asarray(range(1, _np.max(self.xm)))
        self.N = _np.asarray(range(_np.min( self.xn[_np.where(self.xn>0)]), _np.max(self.xn), _np.min(self.xn[_np.where(self.xn>0)])))

        Minv = 1.0/self.M
        self.iota_res = self.N.T * Minv

        # iota_spline = spline( )
        # iota_spline = pchip( 1:self.ns, _np.abs(self.iotaf) )
        iota_spline = pchip( _np.linspace(0, self.ns, self.ns), _np.abs(self.iotaf) )
        self.s_res  = 0.0*self.iota_res

        #for i=1:size(self.s_res,1)
        #    for j=1:size(self.s_res,2)
        #        if (self.iota_res(i,j) > _np.min(_np.abs(self.iotaf))) && (self.iota_res(i,j) < _np.max(_np.abs(self.iotaf)))
        #            f_iota=@(x)ppval(iota_spline,x)-self.iota_res(i,j)
        #            self.s_res(i,j)=fzero(f_iota,50)
        #        #end
        #    #end
        #end

        # Calculate the stored energy
        self.eplasma = 1.5 * _np.pi * _np.pi * _np.sum( self.vp * self.presf ) / self.ns

        # Set some defaults
        self.nu = 2*( self.mpol+1 )+6
        self.datatype = 'wout'
        return self
    #end def CreateResonanceArray()

    # ---------------------------------------------------------------- #
    # ---------------------------------------------------------------- #

    def read_vmec_orig(self,fid,fmt):
        # For 5.10 and earlier

        # Unit Conversions
        dmu0=(2.0e-7)/_np.pi

        # Read Data
        data = fscanf(fid,'#g',13)
        self.wb    = data[0] #changed indices, twice?
        self.wp    = data[1]
        self.gamma = data[2]
        self.pfac  = data[3]
        self.nfp   = data[4]
        self.ns    = data[5]
        self.mpol  = data[6]
        self.ntor  = data[7]
        self.mnmax = data[8]
        self.itfsq = data[9]
        self.niter = data[10]
        self.iasym = data[11]
        self.ireconstruct[12]

        data = fscanf(fid,'#d',5)
        self.imse    = data[0]  #changed indices, twice?
        self.itse    = data[1]
        self.nbsets  = data[2]
        self.nobd    = data[3]
        self.nextcur = data[4]
        self.nstore_seq = 100

        # Error Check
        if self.ierr_vmec and (self.ierr_vmec != 4):
            print('ierr_vmec >0')
            return None
        #end

        # Read nbfld
        if (self.nbsets > 0):
            self.nbfld = fscanf(fid,'#g',self.nbsets)
        #end

        # Read mgrid filename
        if fmt.find(',')<0:  # isempty(strfind(fmt,',')):
            self.mgrid_file = fscanf(fid, '#s', 1)
            fmt2 = '#g#g'
            fmt3 = '#g#g#g'
            fmt6 = '#g#g#g#g#g#g'
            fmt7 = '#g#g#g#g#g#g#g'
            fmt11 = '#g#g#g#g#g#g#g#g#g#g#g'
            fmt12 = '#g#g#g#g#g#g#g#g#g#g#g#g'
            fmt14 = '#g#g#g#g#g#g#g#g#g#g#g#g#g#g'
            fmt2_11 = '#d#d#g#g#g#g#g#g#g#g#g#g#g'
            fmt2_14 = '#d#d#g#g#g#g#g#g#g#g#g#g#g#g#g#g'
        else:
            self.mgrid_file = fscanf(fid, '#s,', 1)
            fmt2 = '#g,#g'
            fmt3 = '#g,#g,#g'
            fmt6 = '#g,#g,#g,#g,#g,#g'
            fmt7 = '#g,#g,#g,#g,#g,#g,#g'
            fmt11 = '#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g'
            fmt12 = '#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g'
            fmt14 = '#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g'
            fmt2_11 = '#d,#d,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g'
            fmt2_14 = '#d,#d,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g'
        #end

        # Read Arrays
        if self.iasym > 0:
            data1 = fscanf(fid, fmt2_14, [16, self.mnmax])
            data = fscanf(fid, fmt14, [14, self.mnmax*(self.ns-1)])
        else:
            data1 = fscanf(fid,fmt2_11,[13, self.mnmax])
            data = fscanf(fid,fmt11,[11, self.mnmax*(self.ns-1)])
        #end

        # Extract Data from Arrays #changed indices
        self.xm = data1[0, :]
        self.xn = data1[1, :]

        # First reshape data #changed indices
        self.rmnc = [data1[2, :].T, reshape(data[0, :], self.mnmax, self.ns-1)]
        self.zmns = [data1[3, :].T, reshape(data[1, :], self.mnmax, self.ns-1)]
        self.lmns = [data1[4, :].T, reshape(data[2, :], self.mnmax, self.ns-1)]
        self.bmn = [data1[5, :].T, reshape(data[3, :], self.mnmax, self.ns-1)]
        self.gmn = [data1[6, :].T, reshape(data[4, :], self.mnmax, self.ns-1)]
        self.bsubumn = [data1[7, :].T, reshape(data[5, :], self.mnmax, self.ns-1)]
        self.bsubvmn = [data1[8, :].T, reshape(data[6, :], self.mnmax, self.ns-1)]
        self.bsubsmn = [data1[9, :].T, reshape(data[7, :], self.mnmax, self.ns-1)]
        self.bsupumn = [data1[10, :].T, reshape(data[8, :], self.mnmax, self.ns-1)]
        self.bsupvmn = [data1[11, :].T, reshape(data[9, :], self.mnmax, self.ns-1)]
        self.currvmn = [data1[12, :].T, reshape(data[10, :], self.mnmax, self.ns-1)]

        # Read the half-mesh quantities
        data = fscanf(fid, fmt12, [12, self.ns/2])
        self.iotas = data[0, :] #changed indices
        self.mass  = data[1, :]
        self.pres  = data[2, :]
        self.phip  = data[3, :]
        self.buco  = data[4, :]
        self.bvco  = data[5, :]
        self.phi   = data[6, :]
        self.vp    = data[7, :]
        self.overr = data[8, :]
        self.jcuru = data[9, :]
        self.jcurv = data[10, :]
        self.specw = data[11, :]

        data = fscanf(fid, fmt6, 6)
        self.aspect  = data[0] #changed indices
        self.betatot = data[1]
        self.betapol = data[2]
        self.betator = data[3]
        self.betaxis = data[4]
        self.b0      = data[5]

        # Mercier Criterion
        data = fscanf(fid, fmt6, [6, self.ns-2])
        self.Dmerc  = data[0, :] #changed indices
        self.Dshear = data[1, :]
        self.Dwell  = data[2, :]
        self.Dcurr  = data[3, :]
        self.Dgeod  = data[4, :]
        self.equif  = data[5, :]
        if (self.nextcur > 0):
            self.extcur = fscanf(fid, fmt, self.nextcur)
            self.curlabel = fscanf(fid, fmt, self.nextcur)
        # end if

        data = fscanf(fid, fmt2, [2, self.nstore_seq])
        self.sqt  = data(0,:) #changed indices
        self.wdot = data(0,:)

        # Convert from Internal Units to Physical Units
        self.mass  = self.mass/dmu0
        self.pres  = self.pres/dmu0
        self.jcuru = self.jcuru/dmu0
        self.jcurv = self.jcurv/dmu0
        self.jdotb = self.jdotb/dmu0
        self.phi   = -self.phi# Data and MSE Fits

        if (self.ireconstruct > 0):
            if (self.imse >= 2) or (self.itse >0):
                self.twsgt   = fscanf(fid, fmt, 1)
                self.msewgt  = fscanf(fid, fmt, 1)
                self.isnodes = fscanf(fid, fmt, 1)

                data = fscanf(fid, fmt3, [3, self.isnodes])
                self.sknots  = data[0, :] #changed indices
                self.ystark  = data[1, :]
                self.y2stark = data[2, :]
                self.ipnodes = fscanf(fid, fmt, 1)

                data = fscanf(fid, fmt3, [3, self.ipnodes])
                self.pknots = data[0, :] #changed indices
                self.ythom  = data[1, :]
                self.y2thom = data[2, :]

                data = fscanf(fid, fmt7, [7, (2*self.ns)-1])
                self.anglemse = data[0, :] #changed indices
                self.rmid     = data[1, :]
                self.qmid     = data[2, :]
                self.shear    = data[3, :]
                self.presmid  = data[4, :]
                self.alfa     = data[5, :]
                self.curmid   = data[6, :]

                data = fscanf(fid, fmt3, [3, self.imse])
                self.rstark    = data[0, :] #changed indices
                self.datastark = data[1, :]
                self.qmeas     = data[2, :]

                data = fscanf(fid, fmt2, [2, self.itse])
                self.rthom    = data[0, :] #changed indices
                self.datathom = data[1, :]
            #end

            if (self.nobd > 0):
                data = fscanf(fid, fmt3, [3, self.nobd])
                self.dsiext = data[0, :] #changed indices
                self.plflux = data[1, :]
                self.dsiobt = data[2, :]
                self.flmwgt = fscanf(fid, fmt, 1)
            #end

            nbfldn = _np.sum( nbldf(0:self.nbsets) )
            if (nbfldn > 0):
                for nn in range(nbsets): #n=1:nbsets
                    data = fscanf(fid, fmt3, [3 self.nbfld(nn)])
                    self.bcoil[:, nn]  = data[0, :] #changed indices
                    self.plbfld[:, nn] = data[1, :]
                    self.bbc[:, nn]    = data[2, :]
                #end
                self.bcwgt = fscanf(fid, fmt, 1)
            #end
            self.phidiam = fscanf(fid, fmt, 1)
            self.delphid = fscanf(fid, fmt, 1)

            # Read Limiter and Prout Plotting Specs
            self.nsets   = fscanf(fid, fmt, 1)
            self.nparts  = fscanf(fid, fmt, 1)
            self.nlim    = fscanf(fid, fmt, 1)
            self.nsetsn  = fscanf(fid, fmt, self.nsets)
            self.pfcspec = _np.zeros( (self.nparts,_np.max(self.nxsetsn),self.nsets), dtype=_np.float64)

            for kk in range(self.nsets): #k=1:self.nsets
                for jj in range(self.nsetsn[kk]): #j=1:self.nsetsn(k)
                    for ii in range(self.nparts): #i=1:self.nparts
                        self.pfcspec[ii, jj, kk] = fscanf(fid, fmt, 1)
                    #end
                #end
            #end

            self.limitr = fscanf(fid, fmt, self.nlim)
            self.rlim   = _np.zeros( (_np.max(self.limitr),self.nlim), dtype=_np.float64)
            self.zlim   = _np.zeros_like(self.rlim))

            for jj in range(self.nlim): #j=1:self.nlim
                for ii in range(self.limitr[jj]): #i=1:self.limitr(j)
                    data = fscanf(fid,fmt2,2)
                    self.rlim[ii, jj] = data[0] #changed indices
                    self.zlim[ii, jj] = data[1]
                #end
            #end

            self.nrgrid = fscanf(fid, fmt, 1)
            self.nzgrid = fscanf(fid, fmt, 1)

            self.tokid = fscanf(fid, fmt, 1)
            self.rx1   = fscanf(fid, fmt, 1)
            self.rx2   = fscanf(fid, fmt, 1)
            self.zy1   = fscanf(fid, fmt, 1)
            self.zy2   = fscanf(fid, fmt, 1)
            self.conif = fscanf(fid, fmt, 1)
            self.imatch_phiedge = fscanf(fid, fmt, 1)
        #end
        return
    #end def read_vmec_orig

    # ---------------------------------------------------------------- #

    def read_vmec_605(self,fid,fmt):
        # For 6.05

        # Unit Conversions
        dmu0=(2.0e-7)/_np.pi

        # Read Data
        data = fscanf(fid,'#g',6)
        self.wb        = data[0] #fixed indices
        self.wp        = data[1]
        self.gamma     = data[2]
        self.pfac      = data[3]
        self.rmax_surf = data[4]
        self.rmin_surf = data[5]

        data = fscanf(fid,'#d',10)
        self.nfp          = data[0]  #fixed indices
        self.ns           = data[1]
        self.mpol         = data[2]
        self.ntor         = data[3]
        self.mnmax        = data[4]
        self.itfsq        = data[5]
        self.niter        = data[6]
        self.iasym        = data[7]
        self.ireconstruct = data[8]
        self.ierr_vmec    = data[9]

        data = fscanf(fid,'#d',5)
        self.imse       = data[0]  #fixed indices
        self.itse       = data[1]
        self.nbsets     = data[2]
        self.nobd       = data[3]
        self.nextcur    = data[4]
        self.nstore_seq = 100

        # Error Check
        if self.ierr_vmec and (self.ierr_vmec != 4):
            print('ierr_vmec >0')
            return
        #end

        # Read nbfld
        if (self.nbsets > 0):
            self.nbfld = fscanf(fid, '#g', self.nbsets)
        #end

        # Read mgrid filename
        self.mgrid_file = fscanf(fid, '#s', 1)

        # Read Arrays
        if self.iasym > 0:
            data1 = fscanf(fid,'#d#d#g#g#g#g#g#g#g#g#g#g#g#g#g#g',[16 self.mnmax])
            data  = fscanf(fid,'#g#g#g#g#g#g#g#g#g#g#g#g#g#g',[14 self.mnmax*(self.ns-1)])
        else:
            data1 = fscanf(fid,'#d#d#g#g#g#g#g#g#g#g#g#g#g',[13 self.mnmax])
            data  = fscanf(fid,'#g#g#g#g#g#g#g#g#g#g#g',[11 self.mnmax*(self.ns-1)])
        #end

        # Extract Data from Arrays
        self.xm = data1[0, :] #fixed indices
        self.xn = data1[1, :]

        # First reshape data
        self.rmnc = [data1[2, :].T reshape(data[0, :], self.mnmax, self.ns-1)]  #fixed indices
        self.zmns = [data1[3, :].T reshape(data[1, :], self.mnmax, self.ns-1)]
        self.lmns = [data1[4, :].T reshape(data[2, :], self.mnmax, self.ns-1)]
        self.bmnc = [data1[5, :].T reshape(data[3, :], self.mnmax, self.ns-1)]
        self.gmnc = [data1[6, :].T reshape(data[4, :], self.mnmax, self.ns-1)]
        self.bsubumnc = [data1[7, :].T reshape(data[5, :], self.mnmax, self.ns-1)]
        self.bsubvmnc = [data1[8, :].T reshape(data[6, :], self.mnmax, self.ns-1)]
        self.bsubsmns = [data1[9, :].T reshape(data[7, :], self.mnmax, self.ns-1)]
        self.bsupumnc = [data1[10, :].T reshape(data[8, :], self.mnmax, self.ns-1)]
        self.bsupvmnc = [data1[11, :].T reshape(data[9, :], self.mnmax, self.ns-1)]
        self.currvmnc = [data1[12, :].T reshape(data[10, :], self.mnmax, self.ns-1)]

        # Read the half-mesh quantities
        data = fscanf(fid,'#g#g#g#g#g#g#g#g#g#g#g#g',[12 self.ns/2])
        self.iotas = data[0, :]  #fixed indices
        self.mass  = data[1, :]
        self.pres  = data[1, :]
        self.phip  = data[2, :]
        self.buco  = data[3, :]
        self.bvco  = data[4, :]
        self.phi   = data[5, :]
        self.vp    = data[6, :]
        self.overr = data[7, :]
        self.jcuru = data[8, :]
        self.jcurv = data[9, :]
        self.specw = data[10, :]

        data = fscanf(fid,'#g#g#g#g#g#g',6)
        self.aspect  = data[0] #fixed indices
        self.betatot = data[1]
        self.betapol = data[2]
        self.betator = data[3]
        self.betaxis = data[4]
        self.b0      = data[5]

        # Mercier Criterion
        data = fscanf(fid,'#g#g#g#g#g#g',[6 self.ns-2])
        self.Dmerc  = data[0, :] #fixed indices
        self.Dshear = data[1, :]
        self.Dwell  = data[2, :]
        self.Dcurr  = data[3, :]
        self.Dgeod  = data[4, :]
        self.equif  = data[5, :]

        if (self.nextcur > 0):
            self.extcur = fscanf(fid,'#g',self.nextcur)
            self.curlabel = fscanf(fid,'#g',self.nextcur)
        #end

        data = fscanf(fid,'#g#g',[2 self.nstore_seq])
        self.sqt = data(0,:)  #fixed indices??? Mistake here
        self.wdot = data(0,:)

        # Convert from Internal Units to Physical Units
        self.mass =  self.mass / dmu0
        self.pres =  self.pres / dmu0
        self.jcuru =  self.jcuru / dmu0
        self.jcurv =  self.jcurv / dmu0
        self.jdotb =  self.jdotb / dmu0
        self.phi   = -self.phi    # Data and MSE Fits

        if (self.ireconstruct > 0):
            if (self.imse >= 2) or (self.itse >0):
                self.twsgt = fscanf(fid, '#g', 1)
                self.msewgt = fscanf(fid, '#g', 1)
                self.isnodes = fscanf(fid, '#d', 1)

                data = fscanf(fid,'#g#g#g',[3 self.isnodes])
                self.sknots = data[0, :]
                self.ystark = data[1, :]
                self.y2stark = data[2, :]
                self.ipnodes = fscanf(fid, '#d', 1)

                data = fscanf(fid,'#g#g#g',[3 self.ipnodes])
                self.pknots = data[0, :]
                self.ythom = data[1, :]
                self.y2thom = data[2, :]

                data = fscanf(fid,'#g#g#g#g#g#g#g',[7 (2*self.ns)-1])
                self.anglemse = data[0, :]
                self.rmid = data[1, :]
                self.qmid = data[2, :]
                self.shear = data[3, :]
                self.presmid = data[4, :]
                self.alfa = data[5, :]
                self.curmid = data[6, :]

                data = fscanf(fid,'#g#g#g',[3 self.imse])
                self.rstark = data[0, :]
                self.datastark = data[1, :]
                self.qmeas = data[2, :]

                data = fscanf(fid,'#g#g',[2 self.itse])
                self.rthom = data[0, :]
                self.datathom = data[1, :]
            #end

            if (self.nobd > 0):
                data = fscanf(fid,'#g#g#g',[3,self.nobd])
                self.dsiext = data[0, :]
                self.plflux = data[1, :]
                self.dsiobt = data[2, :]
                self.flmwgt = fscanf(fid, '#g', 1)
            #end
            nbfldn=_np.sum(nbldf(1:self.nbsets))

            if (nbfldn > 0):
                for nn in range(nbsets): #n=1:nbsets
                    data = fscanf(fid,'#g#g#g',[3 self.nbfld[nn]])
                    self.bcoil[:, nn] = data[0, :]
                    self.plbfld[:, nn] = data[1, :]
                    self.bbc[:, nn] = data[2, :]
                #end
                self.bcwgt = fscanf(fid, '#g', 1)
            #end
            self.phidiam = fscanf(fid, '#g', 1)
            self.delphid = fscanf(fid, '#g', 1)

            # Read Limiter and Prout Plotting Specs
            self.nsets = fscanf(fid, '#g', 1)
            self.nparts = fscanf(fid, '#g', 1)
            self.nlim = fscanf(fid, '#g', 1)
            self.nsetsn = fscanf(fid, '#g', self.nsets)
            self.pfcspec=_np.zeros( (self.nparts, _np.max(self.nxsetsn), self.nsets), dtype=_np.float64)

            for kk in range(self.nsets): #k=1:self.nsets
                for jj in range(self.nsetsn[kk]): #j=1:self.nsetsn(k)
                    for ii in range(self.nparts): #i=1:self.nparts
                        self.pfcspec[ii, jj, kk] = fscanf(fid, '#g', 1)
                    #end
                #end
            #end

            self.limitr = fscanf(fid,'#g',self.nlim)
            self.rlim = _np.zeros( (_np.max(self.limitr), self.nlim), dtype=_np.float64)
            self.zlim = _np.zeros( (_np.max(self.limitr), self.nlim), dtype=_np.float64)

            for jj in range(self.nlim): #j=1:self.nlim
                for ii in range(self.limitr[jj]): # i=1:self.limitr(j)
                    data = fscanf(fid,'#g#g',2)
                    self.rlim[ii, jj] = data[1]
                    self.zlim[ii, jj] = data[2]
                #end
            #end
            self.nrgrid = fscanf(fid, '#g', 1)
            self.nzgrid = fscanf(fid, '#g', 1)
            self.tokid = fscanf(fid, '#g', 1)
            self.rx1 = fscanf(fid, '#g', 1)
            self.rx2 = fscanf(fid, '#g', 1)
            self.zy1 = fscanf(fid, '#g', 1)
            self.zy2 = fscanf(fid, '#g', 1)
            self.conif = fscanf(fid, '#g', 1)
            self.imatch_phiedge = fscanf(fid, '#g', 1)
        #end
        return
    #end def read_vmec_605(self,fid,fmt)

    # ---------------------------------------------------------------- #

    def read_vmec_620(self,fid,fmt):
        # For 6.20

        # Unit Conversions
        dmu0=(2.0e-7)/_np.pi

        # Read Data
        data = fscanf(fid,'#g',6) #changed indicies
        self.wb = data[0]
        self.wp = data[1]
        self.gamma = data[2]
        self.pfac = data[3]
        self.rmax_surf = data[4]
        self.rmin_surf = data[5]

        data = fscanf(fid,'#d',10) #changed indicies
        self.nfp = data[0]
        self.ns = data[1]
        self.mpol = data[2]
        self.ntor = data[3]
        self.mnmax = data[4]
        self.itfsq = data[5]
        self.niter = data[6]
        self.iasym = data[7]
        self.ireconstruct = data[8]
        self.ierr_vmec = data[9]

        data = fscanf(fid,'#d',5) #changed indicies
        self.imse = data[0]
        self.itse = data[1]
        self.nbsets = data[2]
        self.nobd = data[3]
        self.nextcur = data[4]
        self.nstore_seq=100

        # Error Check
        if (self.ierr_vmec) and (self.ierr_vmec != 4):
            print('ierr_vmec >0')
            return
        #end

        # Read nbfld
        if (self.nbsets > 0):
            self.nbfld = fscanf(fid, '#g', self.nbsets)
        #end

        # Read mgrid filename
        self.mgrid_file = fscanf(fid, '#s', 1)

        # Read Arrays
        if (self.iasym > 0):
            data1 = fscanf(fid, '#d#d#g#g#g#g#g#g#g#g#g#g#g#g#g#g', [16 self.mnmax])
            data = fscanf(fid,'#g#g#g#g#g#g#g#g#g#g#g#g#g#g',[14 self.mnmax*( self.ns-1 )])
        else:
            data1 = fscanf(fid, '#d#d#g#g#g#g#g#g#g#g#g#g#g', [13 self.mnmax])
            data = fscanf(fid,'#g#g#g#g#g#g#g#g#g#g#g',[11 self.mnmax*( self.ns-1 )])
        #end

        # Extract Data from Arrays
        self.xm = data1[0, :] #changed indicies
        self.xn = data1[1, :]

        # First reshape data
        self.rmnc=[data1[2, :].T reshape(data[0, :], self.mnmax, self.ns-1)] #changed indicies
        self.zmns=[data1[3, :].T reshape(data[1, :], self.mnmax, self.ns-1)]
        self.lmns=[data1[4, :].T reshape(data[2, :], self.mnmax, self.ns-1)]
        self.bmn=[data1[5, :].T reshape(data[3, :], self.mnmax, self.ns-1)]
        self.gmn=[data1[6, :].T reshape(data[4, :], self.mnmax, self.ns-1)]
        self.bsubumnc=[data1[7, :].T reshape(data[5, :], self.mnmax, self.ns-1)]
        self.bsubvmnc=[data1[8, :].T reshape(data[6, :], self.mnmax, self.ns-1)]
        self.bsubsmns=[data1[9, :].T reshape(data[7, :], self.mnmax, self.ns-1)]
        self.bsupumnc=[data1[10, :].T reshape(data[8, :], self.mnmax, self.ns-1)]
        self.bsupvmnc=[data1[11, :].T reshape(data[9, :], self.mnmax, self.ns-1)]
        self.currvmnc=[data1[12, :].T reshape(data[10, :], self.mnmax, self.ns-1)]

        # Read the half-mesh quantities
        data = fscanf(fid,'#g#g#g#g#g#g#g#g#g#g#g#g#g',[13 self.ns/2])
        self.iotas = data[0, :]  #changed indicies
        self.mass = data[1, :]
        self.pres = data[2, :]
        self.beta_vol[3, :]
        self.phip = data[4, :]
        self.buco = data[5, :]
        self.bvco = data[6, :]
        self.phi = data[7, :]
        self.vp = data[8, :]
        self.overr = data[9, :]
        self.jcuru = data[10, :]
        self.jcurv = data[11, :]
        self.specw = data[12, :]

        data = fscanf(fid,'#g#g#g#g#g#g',6)
        self.aspect = data[0]  #changed indicies
        self.betatot = data[1]
        self.betapol = data[2]
        self.betator = data[3]
        self.betaxis = data[4]
        self.b0 = data[5]

        self.isigna = fscanf(fid, '#d\n', 1)
        self.input_extension=strtrim(fgetl(fid))

        data = fscanf(fid,'#g',8)
        self.IonLarmor = data[0]  #changed indicies
        self.VolAvgB = data[1]
        self.RBtor0 = data[2]
        self.RBtor = data[3]
        self.Itor = data[4]
        self.Aminor = data[5]
        self.Rmajor = data[6]
        self.Volume = data[7]

        # Mercier Criterion
        data = fscanf(fid,'#g#g#g#g#g#g',[6 self.ns-2])
        self.Dmerc = data[0, :]  #changed indicies
        self.Dshear = data[1, :]
        self.Dwell = data[2, :]
        self.Dcurr = data[3, :]
        self.Dgeod = data[4, :]
        self.equif = data[5, :]

        if (self.nextcur > 0):
            self.extcur = fscanf(fid, '#g', self.nextcur)
            self.curlabel = fscanf(fid, '#g', self.nextcur)
        #end

        data = fscanf(fid,'#g#g',[2 self.nstore_seq])
        self.sqt = data[0, :]  #changed indicies
        self.wdot = data[1, :]

        data = fscanf(fid,'#g#g',[2 self.nstore_seq])
        self.jdotb = data[0, :] # changed indices
        self.bdotgradv = data[1, :]

        # Data and MSE Fits
        if (self.ireconstruct > 0):
            if ((self.imse >= 2) or (self.itse >0)):
                self.twsgt = fscanf(fid, '#g', 1)
                self.msewgt = fscanf(fid, '#g', 1)
                self.isnodes = fscanf(fid, '#d', 1)

                data = fscanf(fid,'#g#g#g',[3 self.isnodes])
                self.sknots = data[0, :] # changed indices
                self.ystark = data[1, :]
                self.y2stark = data[2, :]
                self.ipnodes = fscanf(fid, '#d', 1)

                data = fscanf(fid,'#g#g#g',[3 self.ipnodes])
                self.pknots = data[0, :] # changed indices
                self.ythom = data[1, :]
                self.y2thom = data[2, :]

                data = fscanf(fid,'#g#g#g#g#g#g#g',[7 (2*self.ns)-1])
                self.anglemse = data[0, :] # changed indices
                self.rmid = data[1, :]
                self.qmid = data[2, :]
                self.shear = data[3, :]
                self.presmid = data[4, :]
                self.alfa = data[5, :]
                self.curmid = data[6, :]

                data = fscanf(fid,'#g#g#g',[3 self.imse])
                self.rstark = data[0, :] # changed indices
                self.datastark = data[1, :]
                self.qmeas = data[2, :]

                data = fscanf(fid,'#g#g',[2 self.itse])
                self.rthom = data[0, :] # changed indices
                self.datathom = data[1, :]
            #end

            if (self.nobd > 0):
                data = fscanf(fid,'#g#g#g',[3,self.nobd])
                self.dsiext = data[0, :] # changed indices
                self.plflux = data[1, :]
                self.dsiobt = data[2, :]
                self.flmwgt = fscanf(fid, '#g', 1)
            #end
            nbfldn=_np.sum(nbldf(1:self.nbsets))

            if (nbfldn > 0):
                for nn in range(nbsets): #n=1:nbsets
                    data = fscanf(fid,'#g#g#g',[3 self.nbfld(nn)])
                    self.bcoil[:, nn] = data[0, :] # changed indices
                    self.plbfld[:, nn] = data[1, :]
                    self.bbc[:, nn] = data[2, :]
                #end
                self.bcwgt = fscanf(fid, '#g', 1)
            #end
            self.phidiam = fscanf(fid, '#g', 1)
            self.delphid = fscanf(fid, '#g', 1)

            # Read Limiter and Prout Plotting Specs
            self.nsets = fscanf(fid, '#g', 1)
            self.nparts = fscanf(fid, '#g', 1)
            self.nlim = fscanf(fid, '#g', 1)
            self.nsetsn = fscanf(fid,'#g',self.nsets)
            self.pfcspec=_np.zeros( (self.nparts, _np.max(self.nxsetsn), self.nsets), dtype=_np.float64)
            for kk in range(self.nsets): #k=1:self.nsets
                for jj in range(self.nsetsn[kk]): #j=1:self.nsetsn(k)
                    for ii in range(self.nparts): #i=1:self.nparts
                        self.pfcspec[ii, jj, kk] = fscanf(fid, '#g', 1)
                    #end
                #end
            #end
            self.limitr = fscanf(fid, '#g', self.nlim)
            self.rlim=_np.zeros( (_np.max(self.limitr),self.nlim), dtype=_np.float64)
            self.zlim=_np.zeros( (_np.max(self.limitr),self.nlim), dtype=_np.float64)

            for jj in range(self.nlim): #j=1:self.nlim
                for ii in range(self.limitr[jj]): #i=1:self.limitr(j)
                    data = fscanf(fid,'#g#g',2)
                    self.rlim[ii, jj] = data[0] # changed indices
                    self.zlim[ii, jj] = data[1]
                #end
            #end
            self.nrgrid = fscanf(fid, '#g', 1)
            self.nzgrid = fscanf(fid, '#g', 1)
            self.tokid = fscanf(fid, '#g', 1)
            self.rx1 = fscanf(fid, '#g', 1)
            self.rx2 = fscanf(fid, '#g', 1)
            self.zy1 = fscanf(fid, '#g', 1)
            self.zy2 = fscanf(fid, '#g', 1)
            self.conif = fscanf(fid, '#g', 1)
            self.imatch_phiedge = fscanf(fid, '#g', 1)
        #end
        return
    #end def read_vmec_620()

    # ---------------------------------------------------------------- #

    def read_vmec_650(self, fid, fmt):
        # function f=read_vmec_650(fid, fmt)
        # For 6.50

        # Unit Conversions
        dmu0=(2.0e-7)/_np.pi

        # Read Data
        data = fscanf(fid, '#g', 6)
        self.wb = data[0]
        self.wp = data[1]
        self.gamma = data[2]
        self.pfac = data[3]
        self.rmax_surf = data[4]
        self.rmin_surf = data[5]

        data = fscanf(fid, '#d', 10)
        self.nfp = data[0]
        self.ns = data[1]
        self.mpol = data[2]
        self.ntor = data[3]
        self.mnmax = data[4]
        self.itfsq = data[5]
        self.niter = data[6]
        self.iasym = data[7]
        self.ireconstruct = data[8]
        self.ierr_vmec = data[9]

        data = fscanf(fid, '#d', 6)
        self.imse = data[0]
        self.itse = data[1]
        self.nbsets = data[2]
        self.nobd = data[3]
        self.nextcur = data[4]
        self.nstore_seq = data[5]

        # Error Check
        if (self.ierr_vmec and (self.ierr_vmec != 4)):
            print('ierr_vmec >0')
            return
        #end

        # Read nbfld
        if (self.nbsets > 0):
            self.nbfld = fscanf(fid, '#g', self.nbsets)
        #end

        # Read mgrid filename
        self.mgrid_file = fscanf(fid, '#s', 1)

        # Read Arrays
        if self.iasym > 0:
            data1 = fscanf(fid, '#d#d#g#g#g#g#g#g#g#g#g#g#g#g#g#g', [16 self.mnmax])
            data = fscanf(fid, '#g#g#g#g#g#g#g#g#g#g#g#g#g#g', [14 self.mnmax*(self.ns-1)])
        else:
            data1 = fscanf(fid, '#d#d#g#g#g#g#g#g#g#g#g#g#g', [13 self.mnmax])
            data = fscanf(fid, '#g#g#g#g#g#g#g#g#g#g#g', [11 self.mnmax*(self.ns-1)])
        #end

        # Extract Data from Arrays
        self.xm = data1[0, :]
        self.xn = data1[1, :]

        # First reshape data
        self.rmnc=[data1[2, :].T reshape(data[0, :],self.mnmax,self.ns-1)]
        self.zmns=[data1[3, :].T reshape(data[1, :],self.mnmax,self.ns-1)]
        self.lmns=[data1[4, :].T reshape(data[2, :],self.mnmax,self.ns-1)]
        self.bmnc=[data1[5, :].T reshape(data[3, :],self.mnmax,self.ns-1)]
        self.gmnc=[data1[6, :].T reshape(data[4, :],self.mnmax,self.ns-1)]
        self.bsubumnc=[data1[7, :].T reshape(data[5, :],self.mnmax,self.ns-1)]
        self.bsubvmnc=[data1[8, :].T reshape(data[6, :],self.mnmax,self.ns-1)]
        self.bsubsmns=[data1[9, :].T reshape(data[7, :],self.mnmax,self.ns-1)]
        self.bsupumnc=[data1[10, :].T reshape(data[8, :],self.mnmax,self.ns-1)]
        self.bsupvmnc=[data1[11, :].T reshape(data[9, :],self.mnmax,self.ns-1)]
        self.currvmnc=[data1[12, :].T reshape(data[10, :],self.mnmax,self.ns-1)]

        # Read the half-mesh quantities
        data = fscanf(fid, '#g#g#g#g#g#g#g#g#g#g#g#g#g', [13 self.ns/2])
        self.iotas = data[0, :]
        self.mass = data[1, :]
        self.pres = data[2, :]
        self.beta_vol[3, :]
        self.phip = data[4, :]
        self.buco = data[5, :]
        self.bvco = data[6, :]
        self.phi = data[7, :]
        self.vp = data[8, :]
        self.overr = data[9, :]
        self.jcuru = data[10, :]
        self.jcurv = data[11, :]
        self.specw = data[12, :]

        data = fscanf(fid, '#g#g#g#g#g#g', 6)
        self.aspect = data[0]
        self.betatot = data[1]
        self.betapol = data[2]
        self.betator = data[3]
        self.betaxis = data[4]
        self.b0 = data[5]

        self.isigna = fscanf(fid, '#d\n', 1)
        self.input_extension=strtrim(fgetl(fid))

        data = fscanf(fid, '#g', 8)
        self.IonLarmor = data[0]
        self.VolAvgB = data[1]
        self.RBtor0 = data[2]
        self.RBtor = data[3]
        self.Itor = data[4]
        self.Aminor = data[5]
        self.Rmajor = data[6]
        self.Volume = data[7]

        # Mercier Criterion
        data = fscanf(fid, '#g#g#g#g#g#g', [6 self.ns-2])
        self.Dmerc = data[0, :]
        self.Dshear = data[1, :]
        self.Dwell = data[2, :]
        self.Dcurr = data[3, :]
        self.Dgeod = data[4, :]
        self.equif = data[5, :]

        if (self.nextcur > 0):
            self.extcur = fscanf(fid, '#g', self.nextcur)
            self.curlabel = fscanf(fid, '#g', self.nextcur)
        #end

        data = fscanf(fid, '#g#g', [2 self.nstore_seq])
        self.sqt = data[0, :]
        self.wdot = data[1, :]

        data = fscanf(fid, '#g#g', [2 self.nstore_seq])
        self.jdotb = data[0, :]
        self.bdotgradv = data[1, :]

        # Data and MSE Fits
        if (self.ireconstruct > 0):
            if ((self.imse >= 2) or (self.itse > 0)):
                self.twsgt = fscanf(fid, '#g', 1)
                self.msewgt = fscanf(fid, '#g', 1)
                self.isnodes = fscanf(fid, '#d', 1)

                data = fscanf(fid, '#g#g#g', [3 self.isnodes])
                self.sknots = data[0, :]
                self.ystark = data[1, :]
                self.y2stark = data[2, :]
                self.ipnodes = fscanf(fid, '#d', 1)

                data = fscanf(fid, '#g#g#g', [3 self.ipnodes])
                self.pknots = data[0, :]
                self.ythom = data[1, :]
                self.y2thom = data[2, :]

                data = fscanf(fid, '#g#g#g#g#g#g#g', [7 (2*self.ns)-1])
                self.anglemse = data[0, :]
                self.rmid = data[1, :]
                self.qmid = data[2, :]
                self.shear = data[3, :]
                self.presmid = data[4, :]
                self.alfa = data[5, :]
                self.curmid = data[6, :]

                data = fscanf(fid, '#g#g#g', [3 self.imse])
                self.rstark = data[0, :]
                self.datastark = data[1, :]
                self.qmeas = data[2, :]

                data = fscanf(fid, '#g#g', [2 self.itse])
                self.rthom = data[0, :]
                self.datathom = data[1, :]
            #end

            if (self.nobd > 0):
                data = fscanf(fid, '#g#g#g', [3,self.nobd])
                self.dsiext = data[0, :]
                self.plflux = data[1, :]
                self.dsiobt = data[2, :]
                self.flmwgt = fscanf(fid, '#g', 1)
            #end

            nbfldn = _np.sum(nbldf[:self.nbsets])
            if (nbfldn > 0):
                for nn in range(nbsets):  # for n=1:nbsets
                    data = fscanf(fid, '#g#g#g', [3 self.nbfld[nn]])
                    self.bcoil[:, nn] = data[0, :]
                    self.plbfld[:, nn] = data[1, :]
                    self.bbc[:, nn] = data[2, :]
                #end
                self.bcwgt = fscanf(fid, '#g', 1)
            #end
            self.phidiam = fscanf(fid, '#g', 1)
            self.delphid = fscanf(fid, '#g', 1)

            # Read Limiter and Prout Plotting Specs
            self.nsets = fscanf(fid, '#g', 1)
            self.nparts = fscanf(fid, '#g', 1)
            self.nlim = fscanf(fid, '#g', 1)
            self.nsetsn = fscanf(fid, '#g', self.nsets)
            self.pfcspec=_np.zeros( (self.nparts,_np.max(self.nxsetsn),self.nsets), dtype=_np.float64)
            for kk in range(self.nsets):  # k=1:self.nsets
                for jj in range(self.nsetsn[kk]):  # for j=1:self.nsetsn(k)
                    for ii in range(self.nparts):  # for i=1:self.nparts
                        self.pfcspec[ii, jj, kk] = fscanf(fid, '#g', 1)
                    #end
                #end
            #end

            self.limitr = fscanf(fid, '#g', self.nlim)
            self.rlim=_np.zeros( (_np.max(self.limitr), self.nlim), dtype=_np.float64)
            self.zlim=_np.zeros( (_np.max(self.limitr), self.nlim), dtype=_np.float64)
            for jj in range(self.nlim):  # j=1:self.nlim
                for ii in range(self.limitr[jj]:  # i=1:self.limitr(j)
                    data = fscanf(fid, '#g#g', 2)
                    self.rlim[ii, jj] = data[0]
                    self.zlim[ii, jj] = data[1]
                #end
            #end
            self.nrgrid = fscanf(fid, '#g', 1)
            self.nzgrid = fscanf(fid, '#g', 1)
            self.tokid = fscanf(fid, '#g', 1)
            self.rx1 = fscanf(fid, '#g', 1)
            self.rx2 = fscanf(fid, '#g', 1)
            self.zy1 = fscanf(fid, '#g', 1)
            self.zy2 = fscanf(fid, '#g', 1)
            self.conif = fscanf(fid, '#g', 1)
            self.imatch_phiedge = fscanf(fid, '#g', 1)
        #end
        return
    #end def read_vmec_650()

    # ---------------------------------------------------------------- #

    def read_vmec_695(self, fid, fmt):
        # function f=read_vmec_695(fid,fmt)
        self.lfreeb = 0

        data = fscanf(fid, fmt, 7)
        self.wb = data[0]
        self.wp = data[1]
        self.gamma = data[2]
        self.pfac = data[3]
        self.rmax_surf = data[4]
        self.rmin_surf = data[5]
        self.zmax_surf = data[6]

        data = fscanf(fid, fmt, 10)
        self.nfp = data[0]
        self.ns = data[1]
        self.mpol = data[2]
        self.ntor = data[3]
        self.mnmax = data[4]
        self.itfsq = data[5]
        self.niter = data[6]
        self.iasym = data[7]
        self.ireconstruct = data[8]
        self.ierr_vmec = data[9]

        data = fscanf(fid, fmt, 6)
        self.imse = data[0]
        self.itse = data[1]
        self.nbsets = data[2]
        self.nobd = data[3]
        self.nextcur = data[4]
        self.nstore_seq = data[5]

        # Error Check
        if (self.ierr_vmec and (self.ierr_vmec != 4)):
            print('ierr_vmec:'+ str(self.ierr_vmec))
            return
        #end

        # Read nbfld
        if (self.nbsets > 0):
            self.nbfld = fscanf(fid, fmt, self.nbsets)
        #end

        # Read mgrid filename and setup other format statements
        if fmt.find(',')<0:
            #   self.mgrid_file = fscanf(fid, '#s', 1)
            fmt2 = '#d#d#g#g#g#g#g#g#g#g#g#g#g#g#g#g'
            fmt3 = '#g#g#g#g#g#g#g#g#g#g#g#g#g#g'
            fmt4 = '#d#d#g#g#g#g#g#g#g#g#g#g#g'
            fmt5 = '#g#g#g#g#g#g#g#g#g#g#g'
            fmt6 = '#g#g#g#g#g#g#g#g#g#g#g#g#g'
            fmt7 = '#g#g#g#g#g#g'
            fmt8 = '#g#g'
            fmt9 = '#g#g#g'
            fmt10 = '#g#g#g#g#g#g#g'
        else:
            #    self.mgrid_file = fscanf(fid, '#s,', 1)
            fmt2 = '#d,#d,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g'
            fmt3 = '#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g'
            fmt4 = '#d,#d,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g'
            fmt5 = '#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g'
            fmt6 = '#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g'
            fmt7 = '#g,#g,#g,#g,#g,#g'
            fmt8 = '#g,#g'
            fmt9 = '#g,#g,#g'
            fmt10 = '#g,#g,#g,#g,#g,#g,#g'
        #end

        if fmt.find(',')<0:  # isempty(strfind(fmt,','))
            self.mgrid_file = fscanf(fid,'#s',1)
            fmt2 = '#g#g'
            fmt3 = '#g#g#g'
            fmt6 = '#g#g#g#g#g#g'
            fmt7 = '#g#g#g#g#g#g#g'
            fmt11 = '#g#g#g#g#g#g#g#g#g#g#g'
            fmt12 = '#g#g#g#g#g#g#g#g#g#g#g#g'
            fmt13 = '#g#g#g#g#g#g#g#g#g#g#g#g#g'
            fmt14 = '#g#g#g#g#g#g#g#g#g#g#g#g#g#g'
            fmt2_11 = '#d#d#g#g#g#g#g#g#g#g#g#g#g'
            fmt2_14 = '#d#d#g#g#g#g#g#g#g#g#g#g#g#g#g#g'
        else:
            self.mgrid_file = fscanf(fid, '#s,', 1)
            fmt2 = '#g,#g'
            fmt3 = '#g,#g,#g'
            fmt6 = '#g,#g,#g,#g,#g,#g'
            fmt7 = '#g,#g,#g,#g,#g,#g,#g'
            fmt11 = '#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g'
            fmt12 = '#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g'
            fmt13 = '#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g'
            fmt14 = '#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g'
            fmt2_11 = '#d,#d,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g'
            fmt2_14 = '#d,#d,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g'
        #end

        # Read Arrays
        if self.iasym > 0:
            data1 = fscanf(fid, fmt2_14, [16 self.mnmax])
            data = fscanf(fid, fmt14, [14 self.mnmax*(self.ns-1)])
        else
            data1 = fscanf(fid, fmt2_11, [13 self.mnmax])
            data = fscanf(fid, fmt11, [11 self.mnmax*(self.ns-1)])
        #end

        # Extract Data from Arrays
        self.xm = data1[0, :]
        self.xn = data1[1, :]

        # First reshape data
        self.rmnc=[data1[2, :].T reshape(data[0, :],self.mnmax,self.ns-1)]
        self.zmns=[data1[3, :].T reshape(data[1, :],self.mnmax,self.ns-1)]
        self.lmns=[data1[4, :].T reshape(data[2, :],self.mnmax,self.ns-1)] # On half grid
        self.bmnc=[data1[5, :].T reshape(data[3, :],self.mnmax,self.ns-1)]
        self.gmnc=[data1[6, :].T reshape(data[4, :],self.mnmax,self.ns-1)]
        self.bsubumnc=[data1[7, :].T reshape(data[5, :],self.mnmax,self.ns-1)]
        self.bsubvmnc=[data1[8, :].T reshape(data[6, :],self.mnmax,self.ns-1)]
        self.bsubsmns=[data1[9, :].T reshape(data[7, :],self.mnmax,self.ns-1)]
        self.bsupumnc=[data1[10, :].T reshape(data[8, :],self.mnmax,self.ns-1)]
        self.bsupvmnc=[data1[11, :].T reshape(data[9, :],self.mnmax,self.ns-1)]
        self.currvmnc=[data1[12, :].T reshape(data[10, :],self.mnmax,self.ns-1)]
        if self.iasym > 0:
            self.rmns=[data1[13, :].T reshape(data[11, :],self.mnmax,self.ns-1)]
            self.zmnc=[data1[14, :].T reshape(data[12, :],self.mnmax,self.ns-1)]
            self.lmnc=[data1[15, :].T reshape(data[13, :],self.mnmax,self.ns-1)] # On half grid
        #end

        # Read the half-mesh quantities
        data = fscanf(fid,fmt13,[13 self.ns-1])
        self.iotas = data[0, :]
        self.mass = data[1, :]
        self.pres = data[2, :]
        self.beta_vol = data[3, :]
        self.phip = data[4, :]
        self.buco = data[5, :]
        self.bvco = data[6, :]
        self.phi = data[7, :]
        self.vp = data[8, :]
        self.overr = data[9, :]
        self.jcuru = data[10, :]
        self.jcurv = data[11, :]
        self.specw = data[12, :]

        data = fscanf(fid, fmt, 6)
        self.aspect = data[0]
        self.betatot = data[1]
        self.betapol = data[2]
        self.betator = data[3]
        self.betaxis = data[4]
        self.b0 = data[5]

        self.isigna = fscanf(fid, fmt+'\n', 1)
        self.input_extension = strtrim(fgetl(fid))

        data = fscanf(fid, fmt, 8)
        self.IonLarmor = data[0]
        self.VolAvgB = data[1]
        self.RBtor0 = data[2]
        self.RBtor = data[3]
        self.Itor = data[4]
        self.Aminor = data[5]
        self.Rmajor = data[6]
        self.Volume = data[7]

        # Mercier Criterion
        data = fscanf(fid, fmt6, [6 self.ns-2])
        self.Dmerc = data[0, :]
        self.Dshear = data[1, :]
        self.Dwell = data[2, :]
        self.Dcurr = data[3, :]
        self.Dgeod = data[4, :]
        self.equif = data[5, :]

        self.curlabel=cell(self.nextcur, 1)
        if (self.nextcur > 0):
            self.lfreeb = 1
            self.extcur = fscanf(fid, fmt, self.nextcur)
            fscanf(fid, '\n')
            rem = self.nextcur
            jj = 0
            while rem > 0:
                line = fgetl(fid)
                fscanf(fid,'\n')
                test = line[0]
                index = line.find(test)  # findstr(line,test)
                for ii in range(_np.size(index, axis=1)/2):  # i=1:size(index,2)/2
                    self.curlabel[ii+jj] = strtrim(line[index[2*ii-1]+range(index[2*ii]-1])
                #end
                jj += _np.size(index, axis=1)/2
                rem -= _np.size(index, axis=1)/2
            #end
        #end

        data = fscanf(fid, fmt2, [2 self.nstore_seq])
        self.sqt = data[0, :]
        self.wdot = data[1, :]

        data = fscanf(fid, fmt2, [2 self.nstore_seq])
        self.jdotb = data[0, :]
        self.bdotgradv = data[1, :]

        # No Unit Conversion Necessary
        # Data and MSE Fits
        if (self.ireconstruct > 0):
            if ((self.imse >= 2) or (self.itse > 0)):
                self.twsgt = fscanf(fid, fmt, 1)
                self.msewgt = fscanf(fid, fmt, 1)
                self.isnodes = fscanf(fid, fmt, 1)

                data = fscanf(fid, fmt3, [3 self.isnodes])
                self.sknots = data[0, :]
                self.ystark = data[1, :]
                self.y2stark = data[2, :]
                self.ipnodes = fscanf(fid, fmt, 1)

                data = fscanf(fid, fmt3, [3 self.ipnodes])
                self.pknots = data[0, :]
                self.ythom = data[1, :]
                self.y2thom = data[2, :]

                data = fscanf(fid,fmt7,[7 (2*self.ns)-1])
                self.anglemse = data[0, :]
                self.rmid = data[1, :]
                self.qmid = data[2, :]
                self.shear = data[3, :]
                self.presmid = data[4, :]
                self.alfa = data[5, :]
                self.curmid = data[6, :]

                data = fscanf(fid, fmt3, [3 self.imse])
                self.rstark = data[0, :]
                self.datastark = data[1, :]
                self.qmeas = data[2, :]

                data = fscanf(fid, fmt2, [2 self.itse])
                self.rthom = data[0, :]
                self.datathom = data[1, :]
            #end

            if (self.nobd > 0):
                data = fscanf(fid, fmt3, [3, self.nobd])
                self.dsiext = data[0, :]
                self.plflux = data[1, :]
                self.dsiobt = data[2, :]
                self.flmwgt = fscanf(fid, fmt, 1)
            #end

            nbfldn = _np.sum(nbldf[:self.nbsets])
            if (nbfldn > 0):
                for nn in range(nbsets):  # for n=1:nbsets
                    data = fscanf(fid, fmt3, [3 self.nbfld[nn]])
                    self.bcoil[:, nn] = data[0, :]
                    self.plbfld[:, nn] = data[1, :]
                    self.bbc[:, nn] = data[2, :]
                #end
                self.bcwgt = fscanf(fid, fmt, 1)
            #end
            self.phidiam = fscanf(fid, fmt, 1)
            self.delphid = fscanf(fid, fmt, 1)

            # Read Limiter and Prout Plotting Specs
            self.nsets = fscanf(fid, fmt, 1)
            self.nparts = fscanf(fid, fmt, 1)
            self.nlim = fscanf(fid, fmt, 1)
            self.nsetsn = fscanf(fid, fmt, self.nsets)
            self.pfcspec=_np.zeros( (self.nparts, _np.max(self.nxsetsn), self.nsets), dtype=_np.float64)
            for kk in range(self.nsets):  # k=1:self.nsets
                for jj in range(self.nsetsn[kk]):  # j=1:self.nsetsn(k)
                    for ii in range(self.nparts):  # i=1:self.nparts
                        self.pfcspec[ii, jj, kk] = fscanf(fid, fmt, 1)
                    #end
                #end
            #end

            self.limitr = fscanf(fid, fmt, self.nlim)
            self.rlim = _np.zeros( (_np.max(self.limitr), self.nlim), dtype=_np.float64)
            self.zlim = _np.zeros( (_np.max(self.limitr), self.nlim), dtype=_np.float64)
            for jj in range(self.nlim):  # j=1:self.nlim
                for ii in range(self.limitr[jj]):  # i=1:self.limitr(j)
                    data = fscanf(fid, fmt2, 2)
                    self.rlim[ii, jj] = data[0]
                    self.zlim[ii, jj] = data[1]
                #end
            #end

            self.nrgrid = fscanf(fid, fmt, 1)
            self.nzgrid = fscanf(fid, fmt, 1)
            self.tokid = fscanf(fid, fmt, 1)
            self.rx1 = fscanf(fid, fmt, 1)
            self.rx2 = fscanf(fid, fmt, 1)
            self.zy1 = fscanf(fid, fmt, 1)
            self.zy2 = fscanf(fid, fmt, 1)
            self.conif = fscanf(fid, fmt, 1)
            self.imatch_phiedge = fscanf(fid, fmt, 1)
        #end
        return
    #end def read_vmec_695

    # ---------------------------------------------------------------- #

    def read_vmec_800(self, fid, fmt):
        # function f=read_vmec_800(fid,fmt)
        self.lfreeb = 0
        data = fscanf(fid, fmt, 7)
        self.wb = data[0]
        self.wp = data[1]
        self.gamma = data[2]
        self.pfac = data[3]
        self.rmax_surf = data[4]
        self.rmin_surf = data[5]
        self.zmax_surf = data[6]

        data = fscanf(fid, fmt, 10)
        self.nfp = data[0]
        self.ns = data[1]
        self.mpol = data[2]
        self.ntor = data[3]
        self.mnmax = data[4]
        self.itfsq = data[5]
        self.niter = data[6]
        self.iasym = data[7]
        self.ireconstruct = data[8]
        self.ierr_vmec = data[9]
        self.mnmax_nyq = self.mnmax

        data = fscanf(fid, fmt, 6)
        self.imse = data[0]
        self.itse = data[1]
        self.nbsets = data[2]
        self.nobd = data[3]
        self.nextcur = data[4]
        self.nstore_seq = data[5]

        # Error Check
        if (self.ierr_vmec && (self.ierr_vmec != 4))
            print(strcat('ierr_vmec:',num2str(self.ierr_vmec)))
            return
        #end

        # Read nbfld
        if (self.nbsets > 0)
            self.nbfld = fscanf(fid, fmt, self.nbsets)
        #end

        # Read mgrid filename and setup other format statements
        if fmt.find(',')<0:  # isempty(strfind(fmt,','))
            self.mgrid_file = fscanf(fid,'#s',1)
            fmt2 = '#g#g'
            fmt3 = '#g#g#g'
            fmt6 = '#g#g#g#g#g#g'
            fmt7 = '#g#g#g#g#g#g#g'
            fmt10 = '#g#g#g#g#g#g#g#g#g#g'
            fmt11 = '#g#g#g#g#g#g#g#g#g#g#g'
            fmt12 = '#g#g#g#g#g#g#g#g#g#g#g#g'
            fmt13 = '#g#g#g#g#g#g#g#g#g#g#g#g#g'
            fmt14 = '#g#g#g#g#g#g#g#g#g#g#g#g#g#g'
            fmt2_11 = '#d#d#g#g#g#g#g#g#g#g#g#g#g'
            fmt2_14 = '#d#d#g#g#g#g#g#g#g#g#g#g#g#g#g#g'
        else:
            self.mgrid_file = fscanf(fid, '#s,', 1)
            fmt2 = '#g,#g'
            fmt3 = '#g,#g,#g'
            fmt6 = '#g,#g,#g,#g,#g,#g'
            fmt7 = '#g,#g,#g,#g,#g,#g,#g'
            fmt11 = '#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g'
            fmt12 = '#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g'
            fmt13 = '#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g'
            fmt14 = '#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g'
            fmt2_11 = '#d,#d,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g'
            fmt2_14 = '#d,#d,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g'
        #end

        # Read Arrays
        if self.iasym > 0:
            data1 = fscanf(fid, fmt2_14, [16 self.mnmax])
            data = fscanf(fid, fmt14, [14 self.mnmax*(self.ns-1)])
        else:
            data1 = fscanf(fid, fmt2_11, [13 self.mnmax])
            data = fscanf(fid, fmt11, [11 self.mnmax*(self.ns-1)])
        #end
        # Extract Data from Arrays
        self.xm = data1[0, :]
        self.xn = data1[1, :]

        # First reshape data
        self.rmnc=[data1[2, :].T reshape(data[0, :], self.mnmax,self.ns-1)]
        self.zmns=[data1[3, :].T reshape(data[1, :], self.mnmax,self.ns-1)]
        self.lmns=[data1[4, :].T reshape(data[2, :], self.mnmax,self.ns-1)] # On half grid
        self.bmnc=[data1[5, :].T reshape(data[3, :], self.mnmax,self.ns-1)]
        self.gmnc=[data1[6, :].T reshape(data[4, :], self.mnmax,self.ns-1)]
        self.bsubumnc=[data1[7, :].T reshape(data[5, :], self.mnmax,self.ns-1)]
        self.bsubvmnc=[data1[8, :].T reshape(data[6, :], self.mnmax,self.ns-1)]
        self.bsubsmns=[data1[9, :].T reshape(data[7, :], self.mnmax,self.ns-1)]
        self.bsupumnc=[data1[10, :].T reshape(data[8, :], self.mnmax,self.ns-1)]
        self.bsupvmnc=[data1[11, :].T reshape(data[9, :], self.mnmax,self.ns-1)]
        self.currvmnc=[data1[12, :].T reshape(data[10, :], self.mnmax,self.ns-1)]
        if self.iasym >0:
            self.rmns = [data1[13, :].T reshape(data[11, :], self.mnmax,self.ns-1)]
            self.zmnc = [data1[14, :].T reshape(data[12, :], self.mnmax,self.ns-1)]
            self.lmnc = [data1[15, :].T reshape(data[13, :], self.mnmax,self.ns-1)] # On half grid
        #end

        # Read the full-mesh quantities
        data = fscanf(fid, fmt6, [6 self.ns])
        self.iotaf = data[0, :]
        self.presf = data[1, :]
        self.phipf = data[2, :]
        self.phi = data[3, :]
        self.jcuru = data[4, :]
        self.jcurv = data[5, :]

        # Read the half-mesh quantities
        data = fscanf(fid, fmt10, [10 self.ns-1])
        self.iotas = data[0,:]
        self.mass = data[1,:]
        self.pres = data[2,:]
        self.beta_vol = data[3, :]
        self.phip = data[4, :]
        self.buco = data[5, :]
        self.bvco = data[6, :]
        self.vp = data[7, :]
        self.overr = data[8, :]
        self.specw = data[9, :]

        data = fscanf(fid, fmt, 6)
        self.aspect = data[0]
        self.betatot = data[1]
        self.betapol = data[2]
        self.betator = data[3]
        self.betaxis = data[4]
        self.b0 = data[5]
        self.isigna = fscanf(fid, fmt+'\n', 1)
        self.input_extension = strtrim(fgetl(fid))

        data = fscanf(fid, fmt, 8)
        self.IonLarmor = data[0]
        self.VolAvgB = data[1]
        self.RBtor0 = data[2]
        self.RBtor = data[3]
        self.Itor = data[4]
        self.Aminor = data[5]
        self.Rmajor = data[6]
        self.Volume = data[7]

        # Mercier Criterion
        data = fscanf(fid, fmt6, [6 self.ns-2])
        self.Dmerc = data[0, :]
        self.Dshear = data[1, :]
        self.Dwell = data[2, :]
        self.Dcurr = data[3, :]
        self.Dgeod = data[4, :]
        self.equif = data[5, :]

        self.curlabel = cell(self.nextcur, 1)
        if (self.nextcur > 0):
            self.lfreeb = 1
            self.extcur = fscanf(fid, fmt, self.nextcur)
            fscanf(fid, '\n')
            rem = self.nextcur
            jj = 0
            while rem > 0:
                line = fgetl(fid)
                fscanf(fid, '\n')
                test = line[0]
                index = findstr(line, test)
                for ii in range(_np.size(index, axis=1)/2):  # i=1:size(index, 2)/2
                    self.curlabel[ii+jj] = strtrim(line[index[2*ii-1]+range(index[2*ii])-1])
                #end
                jj += _np.size(index, axis=1)/2
                rem -= _np.size(index, axis=1)/2
            #end
        #end

        data = fscanf(fid, fmt2, [2 self.nstore_seq])
        self.sqt = data[0, :]
        self.wdot = data[1, :]

        data = fscanf(fid, fmt2, [2 self.nstore_seq])
        self.jdotb = data[0, :]
        self.bdotgradv = data[1, :]

        # No Unit Conversion Necessary
        # Data and MSE Fits
        if (self.ireconstruct > 0):
            if ((self.imse >= 2) or (self.itse > 0)):
                self.twsgt = fscanf(fid, fmt, 1)
                self.msewgt = fscanf(fid, fmt, 1)
                self.isnodes = fscanf(fid, fmt, 1)

                data = fscanf(fid, fmt3, [3 self.isnodes])
                self.sknots = data[0, :]
                self.ystark = data[1, :]
                self.y2stark = data[2, :]
                self.ipnodes = fscanf(fid, fmt, 1)

                data = fscanf(fid, fmt3, [3 self.ipnodes])
                self.pknots = data[0, :]
                self.ythom = data[1, :]
                self.y2thom = data[2, :]

                data = fscanf(fid,fmt7,[7 (2*self.ns)-1])
                self.anglemse = data[0, :]
                self.rmid = data[1, :]
                self.qmid = data[2, :]
                self.shear = data[3, :]
                self.presmid = data[4, :]
                self.alfa = data[5, :]
                self.curmid = data[6, :]

                data = fscanf(fid, fmt3, [3 self.imse])
                self.rstark = data[0, :]
                self.datastark = data[1, :]
                self.qmeas = data[2, :]

                data = fscanf(fid, fmt2, [2 self.itse])
                self.rthom = data[0, :]
                self.datathom = data[1, :]
            #end

            if (self.nobd > 0):
                data = fscanf(fid, fmt3, [3,self.nobd])
                self.dsiext = data[0, :]
                self.plflux = data[1, :]
                self.dsiobt = data[2, :]
                self.flmwgt = fscanf(fid, fmt, 1)
            #end

            nbfldn = _np.sum(nbldf[:self.nbsets])
            if (nbfldn > 0):
                for nn in range(nbsets):  # for n=1:nbsets
                    data = fscanf(fid, fmt3, [3 self.nbfld[nn]])
                    self.bcoil[:, nn] = data[0, :]
                    self.plbfld[:, nn] = data[1, :]
                    self.bbc[:, nn] = data[2, :]
                #end
                self.bcwgt = fscanf(fid, fmt, 1)
            #end
            self.phidiam = fscanf(fid, fmt, 1)
            self.delphid = fscanf(fid, fmt, 1)

            # Read Limiter and Prout Plotting Specs
            self.nsets = fscanf(fid, fmt, 1)
            self.nparts = fscanf(fid, fmt, 1)
            self.nlim = fscanf(fid, fmt, 1)
            self.nsetsn = fscanf(fid, fmt, self.nsets)
            self.pfcspec = _np.zeros( (self.nparts, _np.max(self.nxsetsn), self.nsets), dtype=_np.float64)
            for kk in range(self.nsets):  # k=1:self.nsets
                for jj in range(self.nsetns[kk]):  # j=1:self.nsetsn[kk]
                    for ii in range(self.nparts):  # i=1:self.nparts
                        self.pfcspec[ii, jj, kk] = fscanf(fid, fmt, 1)
                    #end
                #end
            #end

            self.limitr = fscanf(fid, fmt, self.nlim)
            self.rlim = _np.zeros( (_np.max(self.limitr), self.nlim), dtype=_np.float64)
            self.zlim = _np.zeros( (_np.max(self.limitr), self.nlim), dtype=_np.float64)
            for jj in range(self.nlim):  # j=1:self.nlim
                for ii in range(self.limitr[jj]):  # i=1:self.limitr(j)
                    data = fscanf(fid, fmt2, 2)
                    self.rlim[ii, jj] = data[0]
                    self.zlim[ii, jj] = data[1]
                #end
            #end

            self.nrgrid = fscanf(fid, fmt, 1)
            self.nzgrid = fscanf(fid, fmt, 1)
            self.tokid = fscanf(fid, fmt, 1)
            self.rx1 = fscanf(fid, fmt, 1)
            self.rx2 = fscanf(fid, fmt, 1)
            self.zy1 = fscanf(fid, fmt, 1)
            self.zy2 = fscanf(fid, fmt, 1)
            self.conif = fscanf(fid, fmt, 1)
            self.imatch_phiedge = fscanf(fid, fmt, 1)
        #end

        return
    #end def read_vmec_800()

    # ---------------------------------------------------------------- #

    def read_vmec_847(self,fid,fmt):
        self.lfreeb = 0

        data = fscanf(fid, fmt, 7)
        self.wb = data[0]
        self.wp = data[1]
        self.gamma = data[2]
        self.pfac = data[3]
        self.rmax_surf = data[4]
        self.rmin_surf = data[5]
        self.zmax_surf = data[6]

        data = fscanf(fid, fmt, 11)
        self.nfp = data[0]
        self.ns = data[1]
        self.mpol = data[2]
        self.ntor = data[3]
        self.mnmax = data[4]
        self.mnmax_nyq = data[5]
        self.itfsq = data[5]
        self.niter = data[6]
        self.iasym = data[8]
        self.ireconstruct = data[9]
        self.ierr_vmec = data[10]

        data = fscanf(fid, fmt, 6)
        self.imse = data[0]
        self.itse = data[1]
        self.nbsets = data[2]
        self.nobd = data[3]
        self.nextcur = data[4]
        self.nstore_seq = data[5]

        # Error Check
        if self.ierr_vmec and (self.ierr_vmec != 4):
            print(strcat('ierr_vmec:',num2str(self.ierr_vmec)))
            return
        #end

        # Read nbfld
        if (self.nbsets > 0):
            self.nbfld = fscanf(fid, fmt, self.nbsets)
        #end

        # Read mgrid filename and setup other format statements
        if fmt.find(',')<0:  # isempty(strfind(fmt,',')):
            self.mgrid_file = fscanf(fid,'#s',1)
            fmt2 = '#g#g'
            fmt3 = '#g#g#g'
            fmt6 = '#g#g#g#g#g#g'
            fmt7 = '#g#g#g#g#g#g#g'
            fmt10 = '#g#g#g#g#g#g#g#g#g#g'
            fmt13 = '#g#g#g#g#g#g#g#g#g#g#g#g#g'
            fmt14 = '#g#g#g#g#g#g#g#g#g#g#g#g#g#g'
            fmt20 = '#g#g#g#g#g#g#g#g#g#g#g#g#g#g#g#g#g#g#g#g'
            fmt2_3_2_7 = '#d#d#g#g#g#d#d#g#g#g#g#g#g#g'
            fmt2_6_2_14 = '#d#d#g#g#g#g#g#g#d#d#g#g#g#g#g#g#g#g#g#g#g#g#g#g'
            fmt2_7 = '#d#d#g#g#g#g#g#g#g'
            fmt2_14 = '#d#d#g#g#g#g#g#g#g#g#g#g#g#g#g#g'
        else:
            self.mgrid_file = fscanf(fid,'#s,',1)
            fmt2 = '#g,#g'
            fmt3 = '#g,#g,#g'
            fmt6 = '#g,#g,#g,#g,#g,#g'
            fmt7 = '#g,#g,#g,#g,#g,#g,#g'
            fmt13 = '#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g,#g'
            fmt2_3 = '#d,#d,#g,#g,#g'
            fmt2_6 = '#d,#d,#g,#g,#g,#g,#g,#g'
        #end

        # Read Arrays
        self.xm = _np.zeros( (1,self.mnmax), dtype=_np.int64)
        self.xn = _np.zeros_like(self.xn)

        self.rmnc = _np.zeros( (self.mnmax,self.ns), dtype=_np.float64)
        self.zmns = _np.zeros_like(self.rmnc)
        self.lmns = _np.zeros_like(self.rmnc)

        self.xm_nyq = _np.zeros( (1,self.mnmax_nyq), dtype_np.int64)
        self.xn_nyq = _np.zeros_like(self.xm_nyq)

        self.bmnc = _np.zeros( (self.mnmax_nyq,self.ns), dtype=_np.float64)
        self.gmnc = _np.zeros_like(self.bmnc)
        self.bsubumnc = _np.zeros_like(self.bmnc)
        self.bsubvmnc = _np.zeros_like(self.bmnc)
        self.bsubsmns = _np.zeros_like(self.bmnc)
        self.bsupumnc = _np.zeros_like(self.bmnc)
        self.bsupvmnc = _np.zeros_like(self.bmnc)

        if self.iasym >0:
            self.rmns = _np.zeros( (self.mnmax,self.ns), dtype=_np.float64)
            self.zmnc = _np.zeros_like(self.rmns)
            self.lmnc = _np.zeros_like(self.rmns)

            self.bmns = _np.zeros( (self.mnmax_nyq,self.ns), dtype=_np.float64)
            self.gmns = _np.zeros_like(self.bmns)
            self.bsubumns = _np.zeros_like(self.bmns)
            self.bsubvmns = _np.zeros_like(self.bmns)
            self.bsubsmnc = _np.zeros_like(self.bmns)
            self.bsupumns = _np.zeros_like(self.bmns)
            self.bsupvmns = _np.zeros_like(self.bmns)
        #end

        for ii in range(self.ns): #i=1:self.ns
            for jj in range(self.mnmax): #j=1:self.mnmax
                if ii==1:
                    self.xm[jj] = fscanf(fid, '#d', 1)
                    self.xn[jj] = fscanf(fid, '#d', 1)
                #end
                self.rmnc[jj, ii] = fscanf(fid, '#g', 1)
                self.zmns[jj, ii] = fscanf(fid, '#g', 1)
                self.lmns[jj, ii] = fscanf(fid, '#g', 1)

                if self.iasym > 0:
                    self.rmns[jj, ii] = fscanf(fid, '#g', 1)
                    self.zmnc[jj, ii] = fscanf(fid, '#g', 1)
                    self.lmnc[jj, ii] = fscanf(fid, '#g', 1)
                #end
            #end

            for jj in range(self.mnmax_nyq): #j=1:self.mnmax_nyq
                if ii==1:
                    self.xm_nyq[jj] = fscanf(fid, '#d', 1)
                    self.xn_nyq[jj] = fscanf(fid, '#d', 1)
                #end

                self.bmnc[jj, ii] = fscanf(fid, '#g', 1)
                self.gmnc[jj, ii] = fscanf(fid, '#g', 1)
                self.bsubumnc[jj, ii] = fscanf(fid, '#g', 1)
                self.bsubvmnc[jj, ii] = fscanf(fid, '#g', 1)
                self.bsubsmns[jj, ii] = fscanf(fid, '#g', 1)
                self.bsupumnc[jj, ii] = fscanf(fid, '#g', 1)
                self.bsupvmnc[jj, ii] = fscanf(fid, '#g', 1)

                if self.iasym > 0:
                    self.bmns[jj, ii] = fscanf(fid, '#g', 1)
                    self.gmns[jj, ii] = fscanf(fid, '#g', 1)
                    self.bsubumns[jj, ii] = fscanf(fid, '#g', 1)
                    self.bsubvmns[jj, ii] = fscanf(fid, '#g', 1)
                    self.bsubsmnc[jj, ii] = fscanf(fid, '#g', 1)
                    self.bsupumns[jj, ii] = fscanf(fid, '#g', 1)
                    self.bsupvmns[jj, ii] = fscanf(fid, '#g', 1)
                #end
            #end
        #end

        self.mnyq = _np.max(self.xm_nyq)
        self.nnyq = _np.max(self.xn_nyq)/self.nfp

        # Calculate the Currents
        self.currumnc = _np.zeros( (self.mnmax_nyq,self.ns), dtype=_np.float64)
        self.currvmnc = _np.zeros_like(self.currumnc)

        for ii in range(1,self.ns-1): #i=2:self.ns-1
            self.currumnc[:, ii] = -self.xn_nyq(:)*self.bsubsmns[:, ii] - (self.bsubvmnc[:, ii+1] - self.bsubvmnc[:, ii-1])/(self.ns-1)
            self.currvmnc[:, ii] = -self.xm_nyq(:)*self.bsubsmns[:, ii] - (self.bsubumnc[:, ii+1] - self.bsubumnc[:, ii-1])/(self.ns-1)
        #end

        self.currvmnc[self.xm_nyq==0, 0] = 2*self.bsubumnc[ self.xm_nyq == 0, 1]/(self.ns-1)
        self.currumnc[self.xm_nyq==0, 0] = 2*self.currumnc[ self.xm_nyq == 0, 1] - self.currumnc( self.xm_nyq == 0, 2]
        self.currvmnc[self.xm_nyq!=0, 0] = 0.0
        self.currumnc[self.xm_nyq!=0, 0] = 0.0
        self.currumnc[:, self.ns] = 2*self.currumnc[:, self.ns-1] - self.currumnc[:, self.ns-2]
        self.currvmnc[:, self.ns] = 2*self.currvmnc[:, self.ns-1] - self.currvmnc[:, self.ns-2]
        self.currumnc = self.currumnc/(4*_np.pi*1e-7)
        self.currvmnc = self.currvmnc/(4*_np.pi*1e-7)

        if self.iasym > 0:
            self.currumns = _np.zeros( (self.mnmax_nyq, self.ns), dtype=_np.float64)
            self.currvmns = zeros_like(self.currumns)
            for ii in range(1, self.ns-1):  # i=2:self.ns-1
                self.currumns[:, ii] = -self.xn_nyq(:)*self.bsubsmnc[:, ii]-(self.bsubvmns[:, ii+1]-self.bsubvmns[:, ii-1])/(self.ns-1)
                self.currvmns[:, ii] = -self.xm_nyq(:)*self.bsubsmnc[:, ii]-(self.bsubumns[:, ii+1]-self.bsubumns[:, ii-1])/(self.ns-1)
            #end

            self.currvmns[self.xm_nyq == 0, 0] = 2*self.bsubumns[self.xm_nyq == 0, 1]/(self.ns-1)
            self.currumns[self.xm_nyq == 0, 0] = 2*self.currumns[self.xm_nyq == 0, 1]-self.currumns[self.xm_nyq == 0, 2]
            self.currvmns[self.xm_nyq != 0, 0] = 0.0
            self.currumns[self.xm_nyq != 0, 0] = 0.0
            self.currumns[:, self.ns] = 2*self.currumns[:, self.ns-1]-self.currumns[:, self.ns-2]
            self.currvmns[:, self.ns] = 2*self.currvmns[:, self.ns-1]-self.currvmns[:, self.ns-2]
            self.currumns = self.currumns/(4*_np.pi*1e-7)
            self.currvmns = self.currvmns/(4*_np.pi*1e-7)
        #end

        # Read the full-mesh quantities
        data = fscanf(fid,fmt6,[6 self.ns])
        self.iotaf = data[0, :]
        self.presf = data[1, :]
        self.phipf = data[2, :]
        self.phi = data[3, :]
        self.jcuru = data[4, :]
        self.jcurv = data[5, :]

        # Read the half-mesh quantities
        data = fscanf(fid,fmt10,[10 self.ns-1])
        self.iotas = data[0, :]
        self.mass = data[1, :]
        self.pres = data[2, :]
        self.beta_vol = data[3, :]
        self.phip = data[4, :]
        self.buco = data[5, :]
        self.bvco = data[6, :]
        self.vp = data[7, :]
        self.overr = data[8, :]
        self.specw = data[9, :]

        data = fscanf(fid, fmt, 6)
        self.aspect = data[0]
        self.betatot = data[1]
        self.betapol = data[2]
        self.betator = data[3]
        self.betaxis = data[4]
        self.b0 = data[5]

        self.isigna = fscanf(fid, fmt+'\n', 1)
        self.input_extension = strtrim(fgetl(fid))

        data = fscanf(fid, fmt, 8)
        self.IonLarmor = data[0]
        self.VolAvgB = data[1]
        self.RBtor0 = data[2]
        self.RBtor = data[3]
        self.Itor = data[4]
        self.Aminor = data[5]
        self.Rmajor = data[6]
        self.Volume = data[7]

        # Mercier Criterion
        data = fscanf(fid,fmt6,[6 self.ns-2])
        self.Dmerc = data[0, :]
        self.Dshear = data[1, :]
        self.Dwell = data[2, :]
        self.Dcurr = data[3, :]
        self.Dgeod = data[4, :]
        self.equif = data[5, :]

        self.curlabel=cell(self.nextcur, 1)
        if (self.nextcur > 0):
            self.lfreeb = 1
            self.extcur = fscanf(fid, fmt, self.nextcur)
            lcurr = strtrim(fscanf(fid, '#s', 1))
            if lcurr == 'T':  # strcmpi(lcurr,'T'):
                fscanf(fid, '\n')
                rem = self.nextcur
                jj = 0
                while rem > 0:
                    line = fgetl(fid)
                    fscanf(fid, '\n')
                    test = line[0]
                    index = line.find(test)  # findstr(line,test)
                    for ii in range(_np.size(index, axis=1)/2):  # i=1:size(index,2)/2
                        self.curlabel[ii+jj] = strtrim(line[index[2*ii-1]+range(index[2*ii])-1)
                    #end
                    jj += _np.size(index, axis=1)/2
                    rem -= _np.size(index, axis=1)/2
                #end
            #end
        #end

        data = fscanf(fid, fmt2, [2 self.nstore_seq])
        self.sqt = data[0, :]
        self.wdot = data[1, :]

        data = fscanf(fid, fmt2, [2 self.nstore_seq])
        self.jdotb = data[0, :]
        self.bdotgradv = data[1, :]

        # No Unit Conversion Necessary
        # Data and MSE Fits
        if (self.ireconstruct > 0):
            if ((self.imse >= 2) or (self.itse >0)):
                self.twsgt = fscanf(fid, fmt, 1)
                self.msewgt = fscanf(fid, fmt, 1)
                self.isnodes = fscanf(fid, fmt, 1)

                data = fscanf(fid, fmt3, [3 self.isnodes])
                self.sknots = data[0, :]
                self.ystark = data[1, :]
                self.y2stark = data[2, :]
                self.ipnodes = fscanf(fid, fmt, 1)

                data = fscanf(fid, fmt3, [3 self.ipnodes])
                self.pknots = data[0, :]
                self.ythom = data[1, :]
                self.y2thom = data[2, :]

                data = fscanf(fid,fmt7,[7 (2*self.ns)-1])
                self.anglemse = data[0, :]
                self.rmid = data[1, :]
                self.qmid = data[2, :]
                self.shear = data[3, :]
                self.presmid = data[4, :]
                self.alfa = data[5, :]
                self.curmid = data[6, :]

                data = fscanf(fid, fmt3, [3 self.imse])
                self.rstark = data[0, :]
                self.datastark = data[1, :]
                self.qmeas = data[2, :]

                data = fscanf(fid, fmt2, [2 self.itse])
                self.rthom = data[0, :]
                self.datathom = data[1, :]
            #end

            if (self.nobd > 0):
                data = fscanf(fid, fmt3, [3,self.nobd])
                self.dsiext = data[0, :]
                self.plflux = data[1, :]
                self.dsiobt = data[2, :]
                self.flmwgt = fscanf(fid, fmt, 1)
            #end

            nbfldn=_np.sum(nbldf[:self.nbsets])
            if (nbfldn > 0):
                for nn in range(nbsets):  # n=1:nbsets
                    data = fscanf(fid, fmt3, [3, self.nbfld[nn]])
                    self.bcoil[:, nn] = data[0, :]
                    self.plbfld[:, nn] = data[1, :]
                    self.bbc[:, nn] = data[2, :]
                #end
                self.bcwgt = fscanf(fid, fmt, 1)
            #end
            self.phidiam = fscanf(fid, fmt, 1)
            self.delphid = fscanf(fid, fmt, 1)

            # Read Limiter and Prout Plotting Specs
            self.nsets = fscanf(fid, fmt, 1)
            self.nparts = fscanf(fid, fmt, 1)
            self.nlim = fscanf(fid, fmt, 1)
            self.nsetsn = fscanf(fid, fmt, self.nsets)
            self.pfcspec = _np.zeros( (self.nparts, _np.max(self.nxsetsn), self.nsets), dtype=_np.float64)
            for kk in range( self.nsets ):  # k=1:self.nsets
                for jj in range( self.nsetsn[kk]):  # j=1:self.nsetsn(k)
                    for ii in range( self.nparts):  # i=1:self.nparts
                        self.pfcspec[ii, jj, kk] = fscanf(fid, fmt, 1)
                    #end
                #end
            #end

            self.limitr = fscanf(fid, fmt,  self.nlim)
            self.rlim = _np.zeros( (_np.max(self.limitr), self.nlim), dtype=_np.float64)
            self.zlim = _np.zeros( (_np.max(self.limitr), self.nlim), dtype=_np.float64)
            for jj in range(self.nlim):  # j=1:self.nlim
                for ii in range(self.limitr[jj]):  # ii=1:self.limitr(jj)
                    data = fscanf(fid, fmt2, 2)
                    self.rlim[ii, jj] = data[0]
                    self.zlim[ii, jj] = data[1]
                #end
            #end

            self.nrgrid = fscanf(fid, fmt, 1)
            self.nzgrid = fscanf(fid, fmt, 1)
            self.tokid = fscanf(fid, fmt, 1)
            self.rx1 = fscanf(fid, fmt, 1)
            self.rx2 = fscanf(fid, fmt, 1)
            self.zy1 = fscanf(fid, fmt, 1)
            self.zy2 = fscanf(fid, fmt, 1)
            self.conif = fscanf(fid, fmt, 1)
            self.imatch_phiedge = fscanf(fid, fmt, 1)
        #end
        self.mgrid_mode=strtrim(fscanf(fid, '#s'))
        return
    #end def read_vmec_847

    # ---------------------------------------------------------------- #
    # ---------------------------------------------------------------- #


    @staticmethod
    def h2f(var,ns): #fixed indices
        # Map quantitiy from half to full grid
        temp = _np.zeros( (1, ns), dtype=_np.float64)
        temp[0]      =  1.5 *   var[0] - 0.5 *    var[1]
        temp[1:ns-1] = 0.5 * ( var[0:ns-2] + var[1:ns-1])
        #for i=2:ns-1
        #    temp(i)= 0.5 * ( var(i) +    var(i+1) )
        #end
        temp[ns] =  1.5 *   var[ns-1] - 0.5 * var[ns-2]
        newval=temp
        return newval
    #end def h2f()

    ########################################################################

    def half2fullmesh(self):
        # Interpolate various quantities to the full mesh
        # The following quantities are on 1D half mesh
        # iotas, mass, pres, beta_vol, buco, bvco, vp, specw, phip, jdotb,
        # bdotgradv
        # The following quantities are on 1D full mesh
        # iotaf, presf, phi, phipf, jcuru, jcurv
        # The following quantities are on 1D full mesh from 2:ns-1
        # Dmerc, Dshear, Dwell, Dcurr, Dgeod

        # Now we create full mesh quantities for those that have full mesh names
        if !isfield(self, 'iotaf'):
            self.iotaf = h2f( self.iotas, self.ns)
        #end

        if !isfield(self, 'presf'):
            self.presf = h2f( self.pres, self.ns)
        #end

        if !isfield(self, 'phipf'):
            self.phipf = h2f(self.phip, self.ns)
        #end

        # Now the rest of the quantities are remapped to full mesh
        self.beta_vol  = h2f(self.beta_vol , self.ns)
        self.buco      = h2f(self.buco     , self.ns)
        self.bvco      = h2f(self.bvco     , self.ns)
        self.vp        = h2f(self.vp       , self.ns)
        self.overr     = h2f(self.overr    , self.ns)
        self.specw     = h2f(self.specw    , self.ns)
        self.jdotb     = h2f(self.jdotb    , self.ns)
        self.bdotgradv = h2f(self.bdotgradv, self.ns)

        # Check to make sure phi,jcuru and jcurv are the same size as self.ns
        if len(self.phi) != self.ns: #fixed indices
            temp = _np.zeros( (1, self.ns), dtype=_np.float64)
            temp[0:self.ns] = self.phi
            temp[0]         = 2*temp[1]-temp[2]
            temp[self.ns]   = 2*temp[self.ns-1]-temp[self.ns-2]
            self.phi = temp
        #end

        if len(self.jcuru) != self.ns:
            temp = _np.zeros( (1, self.ns), dtype=_np.float64) #fixed indices
            temp[0:self.ns] = self.jcuru
            temp[0]       = 2*temp[1]-temp[2]
            temp[self.ns] = 2*temp[self.ns-1]-temp[self.ns-2]
            self.jcuru = temp
        #end

        if len(self.jcurv) != self.ns:
            temp = _np.zeros( (1,self.ns), dtype=_np.float64) #fixed indices
            temp[0:self.ns] = self.jcurv
            temp[0]         = 2*temp[1]-temp[2]
            temp[self.ns]   = 2*temp[self.ns-1]-temp[self.ns-2]
            self.jcurv=temp
        #end

        # Fix the Stability values by making into ns arrays
        if len(self.Dmerc) != self.ns: #fixed indices
            temp = _np.zeros( (1, self.ns), dtype=_np.float64)
            temp[0:self.ns-1] = self.Dmerc
            self.Dmerc=temp
            self.Dmerc[0]       = 2 * self.Dmerc[1] - self.Dmerc[2]
            self.Dmerc[self.ns] = 2 * self.Dmerc[self.ns-1] - self.Dmerc[self.ns-2]

            temp(1:self.ns-1) = self.Dshear
            self.Dshear = temp
            self.Dshear[0]       = 2 * self.Dshear[1] - self.Dshear[2]
            self.Dshear[self.ns] = 2 * self.Dshear[self.ns-1] - self.Dshear[self.ns-2]

            temp(1:self.ns-1) = self.Dwell
            self.Dwell = temp
            self.Dwell[0]       = 2 * self.Dwell[1] - self.Dwell[2]
            self.Dwell[self.ns] = 2 * self.Dwell[self.ns-1] - self.Dwell[self.ns-2]

            temp(1:self.ns-1) = self.Dcurr
            self.Dcurr = temp
            self.Dcurr[0]       = 2 * self.Dcurr[1] - self.Dcurr[2]
            self.Dcurr[self.ns] = 2 * self.Dcurr[self.ns-1] - self.Dcurr[self.ns-2]

            temp(1:self.ns-1) = self.Dgeod
            self.Dgeod = temp
            self.Dgeod[0]       = 2 * self.Dgeod[1] - self.Dgeod[2]
            self.Dgeod[self.ns] = 2 * self.Dgeod[self.ns-1] - self.Dgeod[self.ns-2]
        #end

    # ===================================================================== #

        # Now do the matrix values

        # First Index (note indexing on vectors is 2:ns when VMEC outputs)
        # fixed indices
        self.lmns[:, 0] =     1.5 *     self.lmns[:, 1] - 0.5 *     self.lmns[:, 2]
        self.bsupumnc[:, 0] = 1.5 * self.bsupumnc[:, 1] - 0.5 * self.bsupumnc[:, 2]
        self.bsupvmnc[:, 0] = 1.5 * self.bsupvmnc[:, 1] - 0.5 * self.bsupvmnc[:, 2]
        self.bsubsmns[:, 0] = 1.5 * self.bsubsmns[:, 1] - 0.5 * self.bsubsmns[:, 2]
        self.bsubumnc[:, 0] = 1.5 * self.bsubumnc[:, 1] - 0.5 * self.bsubumnc[:, 2]
        self.bsubvmnc[:, 0] = 1.5 * self.bsubvmnc[:, 1] - 0.5 * self.bsubvmnc[:, 2]
        self.gmnc[:, 0] =     1.5 *     self.gmnc[:, 1] - 0.5 *     self.gmnc[:, 2]
        self.bmnc[:, 0] =     1.5 *     self.bmnc[:, 1] - 0.5 *     self.bmnc[:, 2]

        # Average
        for ii in range(1,self.ns-1) #i=2:self.ns-1
            self.lmns[:, ii] =     0.5 * (     self.lmns[:, ii] +     self.lmns[:, ii+1] )
            self.bsupumnc[:, ii] = 0.5 * ( self.bsupumnc[:, ii] + self.bsupumnc[:, ii+1] )
            self.bsupvmnc[:, ii] = 0.5 * ( self.bsupvmnc[:, ii] + self.bsupvmnc[:, ii+1] )
            self.bsubsmns[:, ii] = 0.5 * ( self.bsubsmns[:, ii] + self.bsubsmns[:, ii+1] )
            self.bsubumnc[:, ii] = 0.5 * ( self.bsubumnc[:, ii] + self.bsubumnc[:, ii+1] )
            self.bsubvmnc[:, ii] = 0.5 * ( self.bsubvmnc[:, ii] + self.bsubvmnc[:, ii+1] )
            self.gmnc[:, ii]     = 0.5 * (     self.gmnc[:, ii] +     self.gmnc[:, ii+1] )
            self.bmnc[:, ii]     = 0.5 * (     self.bmnc[:, ii] +     self.bmnc[:, ii+1] )
        #end

        # Last Index (note indexing on vectors is 2:ns when VMEC outputs)
        self.lmns[:, self.ns]     = 2.0 * self.lmns[:, self.ns-1] - self.lmns[:, self.ns-2]
        self.bsupumnc[:, self.ns] = 2.0 * self.bsupumnc[:, self.ns-1] - self.bsupumnc[:, self.ns-2]
        self.bsupvmnc[:, self.ns] = 2.0 * self.bsupvmnc[:, self.ns-1] - self.bsupvmnc[:, self.ns-2]
        self.bsubsmns[:, self.ns] = 2.0 * self.bsubsmns[:, self.ns-1] - self.bsubsmns[:, self.ns-2]
        self.bsubumnc[:, self.ns] = 2.0 * self.bsubumnc[:, self.ns-1] - self.bsubumnc[:, self.ns-2]
        self.bsubvmnc[:, self.ns] = 2.0 * self.bsubvmnc[:, self.ns-1] - self.bsubvmnc[:, self.ns-2]
        self.gmnc[:, self.ns] = 2.0 * self.gmnc[:, self.ns-1] - self.gmnc[:, self.ns-2]
        self.bmnc[:, self.ns] = 2.0 * self.bmnc[:, self.ns-1] - self.bmnc[:, self.ns-2]

        if self.iasym >0: # Handle existance of lmnc on half mesh
            self.lmnc[:, 0] =     1.5 *     self.lmnc[:, 1] - 0.5 *     self.lmnc[:, 2] #fixed indices
            self.bsupumns[:, 0] = 1.5 * self.bsupumns[:, 1] - 0.5 * self.bsupumns[:, 2]
            self.bsupvmns[:, 0] = 1.5 * self.bsupvmns[:, 1] - 0.5 * self.bsupvmns[:, 2]
            self.bsubsmnc[:, 0] = 1.5 * self.bsubsmnc[:, 1] - 0.5 * self.bsubsmnc[:, 2]
            self.bsubumns[:, 0] = 1.5 * self.bsubumns[:, 1] - 0.5 * self.bsubumns[:, 2]
            self.bsubvmns[:, 0] = 1.5 * self.bsubvmns[:, 1] - 0.5 * self.bsubvmns[:, 2]
            self.gmns[:, 0] =     1.5 *     self.gmns[:, 1] - 0.5 *     self.gmns[:, 2]
            self.bmns[:, 0] =     1.5 *     self.bmns[:, 1] - 0.5 *     self.bmns[:, 2]

            for ii in range(1,self.ns): #i=2:self.ns-1
                self.lmnc[:, ii]=     0.5 * (     self.lmnc[:, ii] +     self.lmnc[:, ii+1] )
                self.bsupumns[:, ii]= 0.5 * ( self.bsupumns[:, ii] + self.bsupumns[:, ii+1] )
                self.bsupvmns[:, ii]= 0.5 * ( self.bsupvmns[:, ii] + self.bsupvmns[:, ii+1] )
                self.bsubsmnc[:, ii]= 0.5 * ( self.bsubsmnc[:, ii] + self.bsubsmnc[:, ii+1] )
                self.bsubumns[:, ii]= 0.5 * ( self.bsubumns[:, ii] + self.bsubumns[:, ii+1] )
                self.bsubvmns[:, ii]= 0.5 * ( self.bsubvmns[:, ii] + self.bsubvmns[:, ii+1] )
                self.gmns[:, ii]    = 0.5 * (     self.gmns[:, ii] +     self.gmns[:, ii+1] )
                self.bmns[:, ii]    = 0.5 * (     self.bmns[:, ii] +     self.bmns[:, ii+1] )
            #end

            self.lmnc[:,    self.ns]  = 2.0 *     self.lmnc[:, self.ns-1] -     self.lmnc[:,self.ns-2]
            self.bsupumns[:, self.ns] = 2.0 * self.bsupumns[:, self.ns-1] - self.bsupumns[:,self.ns-2]
            self.bsupvmns[:, self.ns] = 2.0 * self.bsupvmns[:, self.ns-1] - self.bsupvmns[:,self.ns-2]
            self.bsubsmnc[:, self.ns] = 2.0 * self.bsubsmnc[:, self.ns-1] - self.bsubsmnc[:,self.ns-2]
            self.bsubumns[:, self.ns] = 2.0 * self.bsubumns[:, self.ns-1] - self.bsubumns[:,self.ns-2]
            self.bsubvmns[:, self.ns] = 2.0 * self.bsubvmns[:, self.ns-1] - self.bsubvmns[:,self.ns-2]
            self.gmns[:,    self.ns] = 2.0 *     self.gmns[:, self.ns-1] -     self.gmns[:,self.ns-2]
            self.bmns[:,    self.ns] = 2.0 *     self.bmns[:, self.ns-1] -     self.bmns[:,self.ns-2]
        #end

        return
    #end def half2fullmesh

    ########################################################################

    def read_vmec_mercier(filname):

        with open(filname, 'r') as fid: #fid=fopen(filname,'r')
            fgetl(fid) # First Header
            fgetl(fid) # Second Header ----

            # Read first line
            line = fgetl(fid)
            val  = sscanf(line,'#e')
            self.s     = val[0] #changed indices
            self.phi   = val[1]
            self.iota  = val[2]
            self.shear = val[3]
            self.vp    = val[4]
            self.well  = val[5]
            self.itor  = val[6]
            self.ditor = val[7]
            self.pres  = val[8]
            self.dpres = val[9]
            line = fgetl(fid)

            while line!='':  # !strcmp(line, ''):
                val = sscanf(line,'#e')
                self.s     = _np.vstack( (self.s,     val[0]) #changed indices
                self.phi   = _np.vstack( (self.phi,   val[1])
                self.iota  = _np.vstack( (self.iota,  val[2])
                self.shear = _np.vstack( (self.shear, val[3])
                self.vp    = _np.vstack( (self.vp,    val[4])
                self.well  = _np.vstack( (self.well,  val[5])
                self.itor  = _np.vstack( (self.itor,  val[6])
                self.ditor = _np.vstack( (self.ditor, val[7])
                self.pres  = _np.vstack( (self.pres,  val[8])
                self.dpres = _np.vstack( (self.dpres, val[9])
                line=fgetl(fid)
            #end

            fgetl(fid)
            fgetl(fid)

            # Read first line
            line = fgetl(fid)
            val  = sscanf(line,'#e')
            self.dmerc  = val[0] #changed indices - TWICE?
            self.dshear = val[1]
            self.dcurr  = val[2]
            self.dwell  = val[3]
            self.dgeod  = val[4]

            while !feof(fid):
                line = fgetl(fid)
                val  = sscanf(line,'#e')
                self.dmerc  = [self.dmerc;  val[1]] #changed indices
                self.dshear = [self.dshear; val[2]]
                self.dcurr  = [self.dcurr;  val[3]]
                self.dwell  = [self.dwell;  val[4]]
                self.dgeod  = [self.dgeod;  val[5]]
            #end
        #implicitly close the file using the with open (otherwise explicitly)
        fclose(fid)
        return
    #end

    ########################################################################

    def read_vmec_jxbout(self,filname):
        with open(filname,'r') as fid:
            #fid = fopen(filname,'r')

            fgetl(fid)           # Blank Line
            line = fgetl(fid)    # Radial Poloidal Toroidal points
            val  = sscanf(line,'#*20c #d #*24c #d #*24c #d')

            self.nrad   = val[0] #changed indices
            self.ntheta = val[1]
            self.nzeta  = val[2]

            line = fgetl(fid)     # mpol ntor
            val  = sscanf(line,'#*18c #d #*18c #d')
            self.mpol = val[0] #changed indices
            self.ntor = val[1]
            for ii in range( 13 ): # i=1:13      #Header stuff
                fgetl(fid)
            #end

            line = fgetl(fid)    # Values
            val  = sscanf(line,'#*18c #e #*17c #e #*11c #e #*17c #e')
            self.tor_flux = val[0] #changed indices, twice?
            self.fnorm    = val[1]
            self.jdotb    = val[2]
            self.bdotv    = val[3]

            line = fgetl(fid)    #percentages
            fgetl(fid)
            fgetl(fid)
            fgetl(fid)

            self.data = _np.zeros( (self.nrad,self.nzeta,self.ntheta,13), dtype=_np.float64)
            for ii in range(nrad): #i=1:nrad:
                for jj in range(self.nzeta): #j=1:self.nzeta
                    line=fgetl(fid)    # Angle Information
                    for kk in range(self.ntheta): #k=1:self.ntheta
                        line = fgetl(fid)    # Data
                        self.data[ii, jj, kk, :] = sscanf(line,'#d #12e')
                    #end
                #end
            #end
        #implicitly close the file using the with open (otherwise explicitly)
        fclose(fid)
        return
    #end def read_vmec_jxbout

    #####################################################################

    def read_vmec_netcdf(self,filname):
        mu0=4*_np.pi*1e-7

        read_netcdf(filname, 'strip', 'flipdim')

        # Now fix named fields so they match the old way of doing things
        self.input_extension = self.inputextension
        self.mgrid_file = char(self.mgridfile)
        self.ierr_vmec = self.ierflag
        self.rmax_surf = self.rmaxsurf
        self.rmin_surf = self.rminsurf
        self.zmax_surf = self.zmaxsurf
        self.ireconstruct = self.lreconlogical
        self.imse = -1
        self.itse = -1
        self.RBtor = self.rbtor
        self.Rmajor = self.Rmajorp
        self.Aminor = self.Aminorp
        self.betatot = self.betatotal
        self.Volume  =self.volumep
        self.VolAvgB = self.volavgB
        self.beta_vol = self.betavol.T
        self.specw = self.specw.T
        self.iasym = self.lasymlogical
        self.freeb = self.lfreeblogical
        self.lfreeb = self.freeb
        self.Itor = self.ctor
        self.Dmerc = self.DMerc
        self.Dwell = self.DWell
        self.Dshear = self.DShear
        self.Dcurr = self.DCurr
        self.Dgeod = self.DGeod

        # Cast some values
        self.ntor = double(self.ntor)
        self.mpol = double(self.mpol)
        self.nfp = double(self.nfp)
        self.ns = float(self.ns)

        # Calculate Currents
        self.currumnc = _np.zeros((self.mnmaxnyq,self.ns), dtype=_np.float64)
        self.currvmnc = _np.zeros((self.mnmaxnyq,self.ns), dtype=_np.float64)

        for ii in range(1,self.ns-1): #i=2:self.ns-1
            self.currumnc[:, ii] = -double(self.xnnyq).T * self.bsubsmns[:, ii] - (self.ns-1) * (self.bsubvmnc[:, ii+1] - self.bsubvmnc[:, ii] )
            self.currvmnc[:, ii] = -double(self.xmnyq).T * self.bsubsmns[:, ii] + (self.ns-1) * (self.bsubumnc[:, ii+1] - self.bsubumnc[:, ii] )
        #end

        self.currumnc[:, 0] = 0.0 #fixed indices
        self.currvmnc[:, 0] = 0.0
        for ii in range(self.mnmaxnyq): #i=1:self.mnmaxnyq
            if ( self.xmnyq[ii]==0 ):
                self.currumnc(ii,0) = 2*self.currumnc(ii,1) - self.currumnc(ii,2)  #fixed indices
                self.currvmnc(ii,0) = 2*(self.ns-1)*self.bsubumnc(ii,1)
            #end
        #end

        self.currumnc[:, self.ns] = 2*self.currumnc[:, self.ns-1] - self.currumnc[:, self.ns-2]
        self.currvmnc[:, self.ns] = 2*self.currvmnc[:, self.ns-1] - self.currvmnc[:, self.ns-2]
        self.currumnc = self.currumnc/mu0
        self.currvmnc = self.currvmnc/mu0
        if self.iasym:
            self.currumns=_np.zeros( (self.mnmaxnyq,self.ns), dtype=_np.float64)
            self.currvmns=_np.zeros( (self.mnmaxnyq,self.ns), dtype=_np.float64)
            for ii in range(1,self.ns-1): #i=2:self.ns-1
                self.currumns[:, ii] = -double(self.xnnyq).T*self.bsubsmnc[:, ii] - (self.ns-1)*(self.bsubvmns[:, ii+1] - self.bsubvmns[:, ii])
                self.currvmns[:, ii] = -double(self.xmnyq).T*self.bsubsmnc[:, ii] + (self.ns-1)*(self.bsubumns[:, ii+1] - self.bsubumns[:, ii])
            #end
            self.currumns[:, 0] = 0.0  #fixed indices
            self.currvmns[:, 0] = 0.0

            for ii in range(self.mnmaxnyq): #i=1:self.mnmaxnyq
                if (self.xmnyq[ii]==0):
                    self.currumns(ii,0) = 2*self.currumns(ii,1) - self.currumns(ii,2)  #fixed indices
                    self.currvmns(ii,0) = 2*(self.ns-1)*self.bsubumns(ii,1)
                #end
            #end
            self.currumns[:, self.ns] = 2*self.currumns[:, self.ns-1] - self.currumns[:, self.ns-2]
            self.currvmns[:, self.ns] = 2*self.currvmns[:, self.ns-1] - self.currvmns[:, self.ns-2]
            self.currumns = self.currumns/mu0
            self.currvmns = self.currvmns/mu0
        #end

        # Remove Renamed Fields
        del self.inputextension # rmfield(self, 'inputextension')
        del self.mgridfile # rmfield(self, 'mgridfile')
        del self.ierflag # rmfield(self, 'ierflag')
        del self.Rmajorp # rmfield(self, 'Rmajorp')
        del self.Aminorp # rmfield(self, 'Aminorp')
        del self.betatotal # rmfield(self, 'betatotal')
        del self.volumep # rmfield(self, 'volumep')
        del self.volavgB # rmfield(self, 'volavgB')
        del self.betavol # rmfield(self, 'betavol')
        del self.lasymlogical # rmfield(self, 'lasymlogical')
        del self.lfreeblogical # rmfield(self, 'lfreeblogical')
        del self.lreconlogical # rmfield(self, 'lreconlogical')
        return
    #end def read_netcdf()

    #####################################################################

    def __delattr__(self, name):
        del self.(name)
    # end def __delattr__()
#end class

# --------------------------------------------------------------------- #

class vmecdso():
    def __init__(filename):
        self.fname = filename
        self.exists = False
        self.seekoffset = 0
        self.filelength = 0
        if _os.path.exists(filename):
            self.exists=True
            self.openFile(filename)
        # endif
    # end def

    def openFile(self, filename):
        try:
            self.fid = open(filename, 'r')

        self.filelength =
        return self.fid
    # end def openFile

    def readFile(self, buffsize, data_offset=0):
       endpos = buffsize + data_offset
       if self.filelength >= endpos:
          ReadBuffSize = buffsize
       else :
          ReadBuffSize = self.filelength - data_offset
      # endif
       data = _np.zeros((ReadBuffSize,2), order='F')
       self.fid.seek(data_offset+self.seek_offset)
       if self.COMM_TYPE == 0:
          bindata = self.fid.read(ReadBuffSize)
          intd = struct.unpack(str(ReadBuffSize)+'b',bindata[0:ReadBuffSize])
       else:
          bindata = self.fid.read(ReadBuffSize*2)
          intd = struct.unpack(str(ReadBuffSize)+'h',bindata[0:ReadBuffSize*2])
       for i in range(len(intd)):
          data[i,1] = intd[i]*self.VERTICAL_GAIN - self.VERTICAL_OFFSET
          data[i,0] = (i+data_offset)*self.HORIZ_INTERVAL + self.HORIZ_OFFSET
       return data

    def closeFile(self, fid):
        self.fid.close()
    # end def closeFile
# --------------------------------------------------------------------- #
# --------------------------------------------------------------------- #


#def with_open_finally(fil):
#    try:
#        fid = open(fil,'r+')
#        if fid<1:
#
#            raise ArchiveError
#    except:
#        jsonsignal = None
#        pass
##        raise
#    finally:
#        if resp is not None:
#            resp.close()
#        # endif
    # end try
#    return jsonsignal








