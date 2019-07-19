# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 15:42:44 2019

@author: gawe
"""
# ===================================================================== #
# ===================================================================== #

#This section is to improve compatability between bython 2.7 and python 3x
from __future__ import absolute_import, division, print_function #, unicode_literals
__metaclass__ = type

import os as _os
import sys
import platform
import numpy.ctypeslib as npct
import ctypes as _ct
import utils as _ut

try:
    from libstell import libstell as _ls
except:
    from .libstell import libstell as _ls
# end try

# ===================================================================== #
# ===================================================================== #


class pyvmec_libstell_base(_ut.Struct):
    def __init__(self,fname, **kwargs):
        self.verbose = kwargs.setdefault('verbose', True)
        self.setup_libstell(**kwargs)

        self.filename = fname
        if fname is None:
           print('Needs an equilibrium file')
        else:
           fname = _ct.c_char_p( fname.encode() )
           self.loadlibstell(fname)
        #endif
    #end def init()

    # ================================================================== #

    def setup_libstell(self, **kwargs):
        library_path = kwargs.pop('libstell_path', None)
        use_stellopt = kwargs.pop('use_stellopt_libstell', True)
        libstell_branch = kwargs.pop('branch', 'Release')

        if library_path is None:
            if use_stellopt:
                if self.verbose: print('Attempting to use the STELLOPT library') # end if
                try:
                    # Default to the stellopt environment variable
                    library_path = _os.path.join(_os.environ["STELLOPT_PATH"],
                                                  'LIBSTELL', libstell_branch)
                except KeyError:
                    if self.verbose:
                        print('Please set environment variable STELLOPT_PATH or specify the libstell path')
                    # end if
                    sys.exit()
                # end try
            else:
                # library location in same spot as this module
                library_path = _os.path.dirname(__file__)
                if self.verbose:
                    print('Attempting to set a local library: %s'%(library_path,))
                # end if
            # end if
        elif not _os.path.exists(library_path):
            if self.verbose:
                print('Path to libstell library does not exist, \n'
                    + 'please enter a valid directory to search: \n%s'%(library_path,))
        # end if
        if self.verbose:
            print('Set the libstell library search path to: %s'%(library_path,))
        self.libstell_path = library_path
    # end def


    def loadlibstell( self, fname=None):
        if fname is None:
            fname = _ct.c_char_p(b'w7x-sc1(reduced).bc')
        # end if
        os = platform.system()
        is_64bits = sys.maxsize > 2**32

        if is_64bits:   # 64-bit architecture
            ls_handle = _ct.c_longlong
            if os=='Windows':
                self.libstell = npct.load_library(
                    _os.path.join( self.libstell_path, "libstell64.dll"),".")
            elif os=='Linux':
                libstell = npct.load_library(
                    _os.path.join( self.libstell_path, "libstell64.so"),".")
        else:
            ls_handle = _ct.c_long
            if os=='Windows':
                libstell = npct.load_library(
                    _os.path.join( self.libstell_path, "libstell.dll"),".")
            elif os=='Linux':
                libstell = npct.load_library(
                    _os.path.join( self.libstell_path, "libstell.so"),".")
            # end if
        # end if
        self.libstell = libstell

        # ================================================================== #

        # load the magnetic configuration file
        # @return -- if the function succeeds, the return value is
        # the address of C3dMesh object;  zero otherwise.
        MC = libstell.MCload(fname) # mc is like self in python
        if MC == 0:
            print('libstell: Could not load magnetic configuration')
        #endif
        self.MC = MC

        # use the property methods to load up information from equil.
        self.getamin()
        self.B00*self.amin  # all this does is initialize the value for B00

        return libstell, MC

    # ================================= #

    def unloadmconf(self, iunit=None):
        if iunit is None:
            iunit = self.iunit
        # end def
        _ls.safe_close(iunit)
    # end def

    def loadfile(self, filename):    # , **kwargs):
        self.__dict__.update(_ls.read_vmec(filename))
    # end def loadfile

    def pcurr(self, xx):
        return _ls.pcurr(xx)

    def pmass(self, xx):
        return _ls.pmass(xx)

    def piota(self, xx):
        return _ls.piota(xx)

    def read_indata_namelist(self, iunit, istat=0):
        return _ls.read_indata_namelist(iunit, istat)
    # end def write_indata_namelist

    def write_indata_namelist(self, iunit, istat=0):
        _ls.write_indata_namelist(iunit, istat)
    # end def write_indata_namelist


#    def read_input_file(self, filename):
#        iunit = kwargs.setdefault('iunit', 28) # unit number of VMEC file to open
#        istat = kwargs.setdefault('istat', 0)  # status flag
#        filestat = kwargs.setdefault('filestat', 'old')
#        fileform = kwargs.setdefault('fileform', 'formatted')
#        record_in = kwargs.setdefault('record_in', 1)
#        access_in = kwargs.setdefault('access_in', 'sequential')
#        delim_in = kwargs.setdefault('delim_in', 'none')
#
#        # ============== #
#
#        self.iunit = iunit
#        self.loader_args_in = (iunit, istat, filename, filestat, fileform,
#                               record_in, access_in, delim_in)
#        # ============== #
#
#        strprint = 'attempting to load %s into unit file handle %i'%(filename,iunit,)
#        if not _os.path.exists(filename):
#            strprint += '\n'
#            strprint += '       file does not exist!'
#            strprint += '       I am not going to even to try to load it :('
#        else:
#            if self.verbose:
#                print(strprint)
#            # end if
#            self.istat = _ls.safe_open(*self.loader_args_in)
#
#            strprint = 'File status %i: '%(self.istat,)
#            if self.istat>0:
#                strprint += 'Successfully loaded '
#            else:
#                strprint += 'Failed to load '
#            # end if
#            strprint += '%s into unit file handle %i\n'%(filename,iunit,)
#        # end if
#        if self.verbose:
#            print(strprint)
#        # end if
#    # end def loadmconf

    def unloadlibstell(self, libstell=None, MC=None):
        if libstell is None:  libstell=self.libstell        # endif
        if MC is None:        MC = self.MC                  # endif

        #Free up memory etc.
        mconf.MCfree(MC)
        #unloadlibrary(mconf)
    #end def unloadmconf

#    def __safe_open(self):
#        """
#        Safely open a vmec file.  Note that we have to specify unit handles
#        and status flags for the fortran library
#        """
#
#    # end def safe_close
#
#    # ================================================================== #
#    # ================================================================== #
#
#
#    def __del__(self):
#        self.__safe_close()

    # ================================================================== #
# end class

# ===================================================================== #
# ===================================================================== #
