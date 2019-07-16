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


def closedpoints(func, a, b, TOL=1e-6, verbose=True):         # f(x)=func(x)
    """
    Closed Simpson's rule for
        \int_a^b f(x) dx
    Divide [a,b] iteratively into h, h/2, h/4, h/8, ... step sizes; and,
    for each step size, evaluate f(x) at a+h, a+3h, a+5h, a+7h, ..., b-3h,
    b-h, noting that other points have already been sampled.

    At each iteration step, data are sampled only where necessary so that
    the total data is represented by adding sampled points from all
    previous steps:
        step 1:     h       a---------------b
        step 2:     h/2     a-------^-------b
        step 3:     h/4     a---^-------^---b
        step 4:     h/8     a-^---^---^---^-b
        total:              a-^-^-^-^-^-^-^-b
    So, for step size of h/n, there are n intervals, and the data are
    sampled at the boundaries including the 2 end points.

    If old = Trapezoid formula for an old step size 2h, then Trapezoid
    formula for the new step size h is obtained by
        new = old/2 + h{f(a+h) + f(a+3h) + f(a+5h) + f(a+7h) +...+ f(b-3h)
            + f(b-h)}
    Also, Simpson formula for the new step size h is given by
        simpson = (4 new - old)/3
    """
    h = b - a
    old2 = old = h * (func(a) + func(b)) / 2.0
    count = 0
    while True:
        h /= 2.0
        x, tsum = a + h, 0
        while x < b:
            tsum += func(x)
            x += 2 * h
        new = old / 2.0 + h * tsum
        new2 = (4 * new - old) / 3.0
        if abs(new2 - old2) < TOL * (1 + abs(old2)): return new2
        old = new       # Trapezoid
        old2 = new2     # Simpson
        count += 1
        if verbose: print('closedpoints(%d): Trapezoid=%s, Simpson=%s' % (count, new, new2)) # end if
# end def closedpoints


def openpoints(func, a, b, TOL=1e-6, verbose=True):           # f(x)=func(x)
    """
    Open Simpson's rule (excluding end points) for
        \int_a^b f(x) dx
    Divide [a,b] iteratively into h, h/3, h/9, h/27, ... step sizes; and,
    for each step size, evaluate f(x) at a+h/2, a+2h+h/2, a+3h+h/2,
    a+5h+h/2, a+6h+h/2, ... , b-3h-h/2, b-2h-h/2, b-h/2, noting that other
    points have already been sampled.

    At each iteration step, data are sampled only where necessary so that
    the total data is represented by adding sampled points from all
    previous steps:
        step 1:     h       a-----------------^-----------------b
        step 2:     h/3     a-----^-----------------------^-----b
        step 3:     h/9     a-^-------^---^-------^---^-------^-b
        total:              a-^---^---^---^---^---^---^---^---^-b
    So, for step size of h/n, there are n intervals, and the data are
    sampled at the midpoints.

    If old = Trapezoid formula for an old step size 3h, then Trapezoid
    formula for the new step size h is obtained by
        new = old/3 + h{f(a+h/2) + f(a+2h+h/2) + f(a+3h+h/2) + f(a+5h+h/2)
            + f(a+6h+h/2) +...+ f(b-3h-h/2) + f(b-2h-h/2) + f(b-h/2)}
    Also, Simpson formula for the new step size h is given by
        simpson = (9 new - old)/8
    """
    h = b - a
    old2 = old = h * func((a + b) / 2.0)
    count = 0
    while True:
        h /= 3.0
        x, tsum = a + 0.5 * h, 0
        while x < b:
            tsum += func(x) + func(x + 2 * h)
            x += 3 * h
        new = old / 3.0 + h * tsum
        new2 = (9 * new - old) / 8.0
        if abs(new2 - old2) < TOL * (1 + abs(old2)): return new2
        old = new       # Trapezoid
        old2 = new2     # Simpson
        count += 1
        if verbose: print('openpoints(%d): Trapezoid=%s, Simpson=%s' % (count, new, new2) ) # end if
# end def openpoints



# ========================================================================== #
# ========================================================================== #

if __name__ == '__main__':
    import os as _os
    _os.remove('utils.testfile')
    import doctest
    doctest.testmod()
# end if

# ========================================================================== #
# ========================================================================== #

