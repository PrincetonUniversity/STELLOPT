# -*- coding: utf-8 -*-
"""
    A set of utilities required within other modules in this package
@author:
    Gavin Weir <gavin.weir@ipp.mpg.de>
 """
# ======================================================================== #
# ======================================================================== #

from __future__ import absolute_import, with_statement, absolute_import, \
                       division, print_function, unicode_literals

__version__ = "2019.07.12.19"
__all__ = ['Struct', 'Spline', 'utils', 'sscanf']

from . import Struct
from .Struct import Struct  # analysis:ignore

from . import Spline
from .Spline import Spline  # analysis:ignore

from . import utils
from .utils import h2f, cfunct, sfunct # analysis:ignore
from .utils import ftell, fgetl, fscanf, fseek, fgets, findstr # analysis:ignore
from .utils import cylsym_odd, cylsym_even, cart2pol, pol2cart # analysis:ignore
from .utils import lagrange_interpolation, newton, closedpoints, openpoints # analysis:ignore

from . import sscanf
from .sscanf import scanf as sscanf  # analysis:ignore

# ======================================================================== #
# ======================================================================== #







