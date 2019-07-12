"""
    Created on Wed Jan 11 14:26:46 2017

    pychi_pylab.py
    Michalski
    20090818

    Piecewise cubic Hermite interpolation (monotonic...) in Python

    References:

        Wikipedia:  Monotone cubic interpolation
                    Cubic Hermite spline

    A cubic Hermte spline is a third degree spline with each polynomial of the spline
    in Hermite form.  The Hermite form consists of two control points and two control
    tangents for each polynomial.  Each interpolation is performed on one sub-interval
    at a time (piece-wise).  A monotone cubic interpolation is a variant of cubic
    interpolation that preserves monotonicity of the data to be interpolated (in other
    words, it controls overshoot).  Monotonicity is preserved by linear interpolation
    but not by cubic interpolation.

    Use:

    There are two separate calls, the first call, pchip_init(), computes the slopes that
    the interpolator needs.  If there are a large number of points to compute,
    it is more efficient to compute the slopes once, rather than for every point
    being evaluated.  The second call, pchip_eval(), takes the slopes computed by
    pchip_init() along with X, Y, and a vector of desired 'xnew's and computes a vector
    of 'ynew's.  If only a handful of points is needed, pchip() is a third function
    which combines a call to pchip_init() followed by pchip_eval().
"""

# ========================================================= #
# ========================================================= #

from __future__ import absolute_import, with_statement, absolute_import, \
                       division, print_function, unicode_literals

import pylab as _pylab

# ========================================================= #
# ========================================================= #

def x_is_okay(x, xvec, verbose=False):
    """
    Make sure 'x' and 'xvex' satisfy the conditions for
    running the pchip interpolator
    """
    n = len(x)
    m = len(xvec)
    
    # Make sure 'x' is in sorted order (brute force, but works...)    
    xx = x.copy()
    xx.sort()
    total_matches = (xx == x).sum()
    if total_matches != n:
        if verbose:
            print('*' * 50)
            print('x_is_okay()')
            print('x values were not in sorted order --- aborting')
        # endif
        return False
    # endif

    # Make sure 'x' doesn't have any repeated values
    delta = x[1:] - x[:-1]
    if (delta == 0.0).any():
        if verbose:
            print('*' * 50)
            print('x_is_okay()')
            print('x values were not monotonic--- aborting')
        # endif
        return False
    # endif

    # Check for in-range xvec values (beyond upper edge)
    check = xvec > x[-1]
    if check.any():
        indices = _pylab.compress(check, range(m))

        if verbose:
            print('*' * 50)
            print('x_is_okay()')
            print('Certain xvec values are beyond the upper end of x')
            print('x_max = ', x[-1])
            print('out-of-range xvec = ', xvec[indices])
            print('out-of-range xvec indices = ', indices)
        # endif
        return False
    # end if

    # Second - check for in-range xvec values (beyond lower edge)
    check = xvec< x[0]
    if check.any():
        indices = _pylab.compress(check, range(m))
        if verbose:
            print('*' * 50)
            print('x_is_okay()')
            print('Certain xvec values are beyond the lower end of x')
            print('x_min = ', x[0])
            print('out-of-range xvec = ', xvec[indices])
            print('out-of-range xvec indices = ', indices)
        # end if
        return False
    # endif
    return True

#=========================================================


def pchip_eval(x, y, m, xvec, verbose=False):
    """
     Evaluate the piecewise cubic Hermite interpolant with monoticity preserved

        x = array containing the x-data
        y = array containing the y-data
        m = slopes at each (x,y) point [computed to preserve monotonicity]
        xnew = new 'x' value where the interpolation is desired

        x must be sorted low to high... (no repeats)
        y can have repeated values

     This works with either a scalar or vector of 'xvec'
    """
    n = len(x)
    mm = len(xvec)

    # --------------------------
    # Make sure there aren't problems with the input data
    # --------------------------
    if not x_is_okay(x, xvec, verbose):
        if verbose:
            print('pchip_eval2() - ill formed x vector!!!!!!!!!!!!!')
        # endif

        # Cause a hard crash...
        return #STOP_pchip_eval2
    # endif

    # Find the indices 'k' such that x[k] < xvec < x[k+1]

    # Create 'copies' of 'x' as rows in a mxn 2-dimensional vector
    xx = _pylab.resize(x,(mm,n)).transpose()
    xxx = xx > xvec

    # Compute column by column differences
    z = xxx[:-1,:] - xxx[1:,:]

    # Collapse over rows...
    k = z.argmax(axis=0)

    # Create the Hermite coefficients
    h = x[k+1] - x[k]
    t = (xvec - x[k]) / h[k]

    # Hermite basis functions
    h00 = (2 * t**3) - (3 * t**2) + 1
    h10 =      t**3  - (2 * t**2) + t
    h01 = (-2* t**3) + (3 * t**2)
    h11 =      t**3  -      t**2

    # Compute the interpolated value of 'y'
    ynew = h00*y[k] + h10*h*m[k] + h01*y[k+1] + h11*h*m[k+1]

    return ynew
# end def

#=========================================================


def pchip_init(x, y, verbose=False):
    """
     Evaluate the piecewise cubic Hermite interpolant with monoticity preserved

        x = array containing the x-data
        y = array containing the y-data

        x must be sorted low to high... (no repeats)
        y can have repeated values

        x input conditioning is assumed but not checked
    """

    n = len(x)

    # Compute the slopes of the secant lines between successive points
    delta = (y[1:] - y[:-1]) / (x[1:] - x[:-1])

    # Initialize the tangents at every points as the average of the secants
    m = _pylab.zeros(n, dtype='d')

    # At the endpoints - use one-sided differences
    m[0] = delta[0]
    m[n-1] = delta[-1]

    # In the middle - use the average of the secants
    m[1:-1] = (delta[:-1] + delta[1:]) / 2.0

    # Special case: intervals where y[k] == y[k+1]

    # Setting these slopes to zero guarantees the spline connecting
    # these points will be flat which preserves monotonicity
    indices_to_fix = _pylab.compress((delta == 0.0), range(n))

#    print('zero slope indices to fix = ', indices_to_fix

    for ii in indices_to_fix:
        m[ii]   = 0.0
        m[ii+1] = 0.0
    # end for

    alpha = m[:-1]/delta
    beta  = m[1:]/delta
    dist  = alpha**2 + beta**2
    tau   = 3.0 / _pylab.sqrt(dist)

    # To prevent overshoot or undershoot, restrict the position vector
    # (alpha, beta) to a circle of radius 3.  If (alpha**2 + beta**2)>9,
    # then set m[k] = tau[k]alpha[k]delta[k] and m[k+1] = tau[k]beta[b]delta[k]
    # where tau = 3/sqrt(alpha**2 + beta**2).

    # Find the indices that need adjustment
    over = (dist > 9.0)
    indices_to_fix = _pylab.compress(over, range(n))

#    print('overshoot indices to fix... = ', indices_to_fix

    for ii in indices_to_fix:
        m[ii]   = tau[ii] * alpha[ii] * delta[ii]
        m[ii+1] = tau[ii] * beta[ii]  * delta[ii]

    return m
# end def

# ==================================================================

def pchip(x, y, xnew, verbose=False):
    # Compute the slopes used by the piecewise cubic Hermite interpolator
    m = pchip_init(x, y, verbose)

    # Use these slopes (along with the Hermite basis function) to interpolate
    ynew = pchip_eval(x, y, m, xnew, verbose)
    return ynew
# end def pchip

#==================================================================


def CubicHermiteSpline(x, y, x_new, verbose=False):
    """
     Piecewise Cubic Hermite Interpolation using Catmull-Rom
     method for computing the slopes.

     Note - this only works if delta-x is uniform?
    """

    # Find the two points which 'bracket' 'x_new'
    found_it = False

    for ii in range(len(x)-1):
        if (x[ii] <= x_new) and (x[ii+1] > x_new):
            found_it = True
            break
        # end if
    # end for

    if not found_it:
        if verbose:
            print('')
            print('requested x=<%f> outside X range[%f,%f]' % (x_new, x[0], x[-1]))
        # endif
        return #STOP_CubicHermiteSpline()
    # endif

    # Starting and ending data points
    x0 = x[ii]
    x1 = x[ii+1]

    y0 = y[ii]
    y1 = y[ii+1]

    # Starting and ending tangents (using Catmull-Rom spline method)

    # Handle special cases (hit one of the endpoints...)
    if ii == 0:

        # Hit lower endpoint
        m0 = (y[1] - y[0])
        m1 = (y[2] - y[0]) / 2.0

    elif ii == (len(x) - 2):

        # Hit upper endpoints
        m0 = (y[ii+1] - y[ii-1]) / 2.0
        m1 = (y[ii+1] - y[ii])

    else:

        # Inside the field...
        m0 = (y[ii+1] - y[ii-1])/ 2.0
        m1 = (y[ii+2] - y[ii])  / 2.0
    # end if

    # Normalize to x_new to [0,1] interval
    h = (x1 - x0)
    t = (x_new - x0) / h

    # Compute the four Hermite basis functions
    h00 = ( 2.0 * t**3) - (3.0 * t**2) + 1.0
    h10 = ( 1.0 * t**3) - (2.0 * t**2) + t
    h01 = (-2.0 * t**3) + (3.0 * t**2)
    h11 = ( 1.0 * t**3) - (1.0 * t**2)

    h = 1

    y_new = (h00 * y0) + (h10 * h * m0) + (h01 * y1) + (h11 * h * m1)

    return y_new
# end def

#==============================================================


def main(verbose=True):
    # ======================================================== #
    # Sine wave test
    # ======================================================== #

    # Create a example vector containing a sine wave.
    x = _pylab.arange(30.0)/10.
    y = _pylab.sin(x)

    # Interpolate the data above to the grid defined by 'xvec'
    xvec = _pylab.arange(250.)/100.

    # Initialize the interpolator slopes
    m = pchip_init(x,y,verbose)

    # Call the monotonic piece-wise Hermite cubic interpolator
    yvec2 = pchip_eval(x, y, m, xvec, verbose)

    _pylab.figure(1)
    _pylab.plot(x,y, 'ro')
    _pylab.title('pchip() Sin test code')

    # Plot the interpolated points
    _pylab.plot(xvec, yvec2, 'b')

    # ======================================================== #
    # Step function test...
    # ======================================================== #

    _pylab.figure(2)
    _pylab.title('pchip() step function test')

    # Create a step function (will demonstrate monotonicity)
    x = _pylab.arange(7.0) - 3.0
    y = _pylab.array([-1.0, -1,-1,0,1,1,1])

    # Interpolate using monotonic piecewise Hermite cubic spline
    xvec = _pylab.arange(599.)/100. - 3.0

    # Create the pchip slopes slopes
    m    = pchip_init(x,y,verbose)

    # Interpolate...
    yvec = pchip_eval(x, y, m, xvec, verbose)

    # Call the Scipy cubic spline interpolator
    from scipy.interpolate import interpolate
    function = interpolate.interp1d(x, y, kind='cubic')
    yvec2 = function(xvec)

    # Non-montonic cubic Hermite spline interpolator using
    # Catmul-Rom method for computing slopes...
    yvec3 = []
    for xx in xvec:
        yvec3.append(CubicHermiteSpline(x,y,xx,verbose))
    yvec3 = _pylab.array(yvec3)

    # Plot the results
    _pylab.plot(x,    y,     'ro')
    _pylab.plot(xvec, yvec,  'b')
    _pylab.plot(xvec, yvec2, 'k')
    _pylab.plot(xvec, yvec3, 'g')
    _pylab.xlabel('X')
    _pylab.ylabel('Y')
    _pylab.title('Comparing pypchip() vs. Scipy interp1d() vs. non-monotonic CHS')
    legends = ['Data', 'pypchip()', 'interp1d','CHS']
    _pylab.legend(legends, loc='upper left')

    _pylab.show()
# end def main

# ================================================================== #

if __name__=="__main__":
    main()
# end if


# ================================================================== #
# ================================================================== #





