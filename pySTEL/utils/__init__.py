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
__all__ = ['scanf', 'Struct', 'Spline', 'utils', 'txtutils', 'vmcutils']


from . import Struct
from .Struct import Struct  # analysis:ignore

from . import Spline
from .Spline import Spline  # analysis:ignore

from . import utils
from .utils import cylsym_odd, cylsym_even, cart2pol, pol2cart # analysis:ignore
from .utils import lagrange_interpolation, newton, closedpoints, openpoints # analysis:ignore

from . import txtutils
#from .txtutils import scanf, sscanf, fscanf  # analysis:ignore
#from .txtutils import ftell, fseek, frewind, fgets, fgetl, findstr # analysis:ignore
from txtutils import scanf, sscanf, fscanf, ftell, fseek, frewind, fgets, fgetl, findstr  # analysis:Ignore

from . import vmcutils
from .vmcutils import calc_curr  # analysis:ignore
from .vmcutils import h2f, h2f_special, cfunct, sfunct # analysis:ignore

# ======================================================================== #
# ======================================================================== #







