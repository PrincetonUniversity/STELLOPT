# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 14:34:06 2019

@author: Gavin Weir <gavin.weir@ipp.mpg.de>
"""

# ===================================================================== #
# ===================================================================== #

from __future__ import absolute_import, with_statement, absolute_import, \
                       division, print_function, unicode_literals
import numpy as _np

__metaclass__ = type

# ===================================================================== #
# ===================================================================== #


def cylsym_odd(rr_in, axis=0, exclude_axis=False):
    """
    This function enforces anti-cylindrical symmetry on a vector/array by
    mirroring the negative of the input array about its first index (row for
    multi-dimensional arrays).
    Inputs:
        rr_in - Data array to be mirrored
    Outputs:
        rr - Mirrored data ~ [-fliplr(dat_in), dat_in]
    """
    if axis==0:
        flipper = _np.flipud
    elif axis==1:
        flipper = _np.fliplr
    # end if

    rr = _np.atleast_1d(rr_in).copy()
    if exclude_axis and _np.shape(rr)[axis]>1:
        if axis == 0:
            r1 = rr[1:]
        elif axis==1:
            r1 = rr[:,1:]
        # end if
    else:
        r1 = rr
    # end if

    rr = _np.concatenate((-flipper(r1), rr), axis=axis)
#    rr = _np.concatenate((-_np.flipud(rr), rr), axis=0)
    return rr
# end def cylsym_odd


def cylsym_even(dat_in, axis=0, exclude_axis=False):
    """
    This function enforces cylindrical symmetry on a vector/array by mirroring
    the input array about its first index (row for multi-dimensional arrays).
    Inputs:
        dat_in - Data array to be mirrored
    Outputs:
        dysm - Mirrored data ~ [fliplr(dat_in), dat_in]
    """
    if axis==0:
        flipper = _np.flipud
    elif axis==1:
        flipper = _np.fliplr
    # end if
    dsym = _np.atleast_1d(dat_in).copy()
    if exclude_axis and _np.shape(dsym)[axis]>1:
        if axis == 0:
            r1 = dsym[1:]
        elif axis==1:
            r1 = dsym[:,1:]
        # end if
    else:
        r1 = dsym
    # end if
    dsym = _np.concatenate((flipper(r1), dsym), axis=axis)
    return dsym
# end def cylsym_even

# ===================================================================== #

def cart2pol(x, y, z=None):
    """
    Convert the Cartesian coordinates [x, y] to the Polar radius and
    angle [r, \theta].

    The parameter r is the radial distance, \theta is the polar angle.


    @param vector:  The Cartesian vector [x, y].
    @type vector:   numpy rank-1, 2D array
    @return:        The Polar coordinates [r, \theta].
    @rtype:         numpy rank-1, 2D array
    """
    return _np.arctan2(y, x), _np.sqrt(x*x + y*y), z


def pol2cart(rho, theta, z=None):
    """
    Convert the Polar radius and angle [r, \theta] to the Cartesian
    plane [x, y].

    The parameter r is the radial distance, \theta is the polar angle.


    @param vector:  The Polar coordinates [r, \theta].
    @type vector:   numpy rank-1, 2D array
    @return:        The Cartesian vector [x, y].
    @rtype:         numpy rank-1, 2D array
    """
    return rho*_np.cos(theta), rho*_np.sin(theta), z
# pol2cart

# ========================================================================== #

def _lagrange(xx, il, xm):
    """
    Evaluates the i-th lagrange polynomial at x
    based on grid data xm
    """
    nn = len(xm)-1
    yy = 1.0
    for jj in range(nn+1):
        if il != jj:
            yy *=(xx-xm[jj])/(xm[il]-xm[jj])
        # end if
    # end for
    return yy

def lagrange_interpolation(xx, xm, ym):
    """
    Polynomial interpolation using the lagrange basis polynomials
    xm, ym are the background grid of data points
    xx - is the point to interpolate (xm, ym) => (x, y)

    """
    nn = len(xm)-1
    lagrpoly = _np.array([_lagrange(xx, ii, xm) for ii in range(nn+1)])
    yy = _np.dot(ym, lagrpoly)
    return yy

# ========================================================================== #


def newton(func, funcd, x, TOL=1e-6, verbose=True):   # f(x)=func(x), f'(x)=funcd(x)
    """
    Ubiquitous Newton-Raphson algorithm for solving
        f(x) = 0
    where a root is repeatedly estimated by
        x = x - f(x)/f'(x)
    until |dx|/(1+|x|) < TOL is achieved.  This termination condition is a
    compromise between
        |dx| < TOL,  if x is small
        |dx|/|x| < TOL,  if x is large
    """
    f, fd = func(x), funcd(x)
    count = 0
    while True:
        dx = f / float(fd)
        if abs(dx) < TOL * (1 + abs(x)): return x - dx
        x -= dx
        f, fd = func(x), funcd(x)
        count += 1
        if verbose:  print("newton(%d): x=%s, f(x)=%s" % (count, x, f))  # endif
# end def newton


# ========================================================================== #
# ========================================================================== #


if __name__ == '__main__':
    import doctest
    doctest.testmod()

#    import os as _os
#    if _os.path.exists('utils.testfile'):
#        _os.remove('utils.testfile')
# end if

# ========================================================================== #
# ========================================================================== #

