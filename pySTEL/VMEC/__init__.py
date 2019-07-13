# -*- coding: utf-8 -*-
"""
    This package is designed to read VMEC data in multiple formats

    The best way to do this is to use the compiled library directly from
    STELLOPT. When libstell.so is present on the python path, we use a python
    wrapper written by Jonathon Schilling.
        - The python wrapper library, libstell.py, only works for netcdf files.

    When libstell.so is not on the python path, or we are using text file
    based output, we have to use separate methods.

@author:
    Gavin Weir <gavin.weir@ipp.mpg.de>
    Jonathon Schilling <jons@ipp.mpg.de>
"""
# ======================================================================== #
# ======================================================================== #

from __future__ import absolute_import, with_statement, absolute_import, \
                       division, print_function, unicode_literals

__version__ = "2019.07.12.19"
__all__ = ['read_vmec']

from . import read_vmec
from .read_vmec import read_vmec # analysis:ignore

#from . import CalcFluxCart_VMEC
#from .CalcFluxCart_VMEC import

# ======================================================================== #
# ======================================================================== #







