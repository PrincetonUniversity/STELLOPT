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


#def trapezoidalarea(func, a, b):
#    return 0.5*(b-a)*(func(a)+func(b))
#
#def adaptint(func, a, b, tol=1e-8, recursion=0, max_recursions=1e2):
#    h = b-1.0
#    m = 0.5*(b+1.0)
#    area = 0.0
#    areatot = trapezoidalarea(func, a, b)
#    nextareatot = trapezoidalarea(func, a, m) + trapezoidalarea(func, m, b)
#    err = _np.abs(areatot-nextareatot)
#
#    if err<tol:
#        return areatot
#    elif recursion>max_recursions:
#        print('Maximum recursion depth reached, aborting adaptive integration')
#        return _np.nan
#    else:
#        arealeft = adaptint(func, a, m, 0.5*tol, recursion=recursion+1, max_recursions=max_recursions)
#        arearight = adaptint(func, m, b, 0.5*tol, recursion=recursion+1, max_recursions=max_recursions)
#        area = area + arealeft + arearight
#    return area

def trapezoidal(f, a, b, n):
    """ Vectorized form of trapezoidal function """
    h = float(b-a)/n
    x = _np.linspace(a, b, n+1)
    s = _np.nansum(f(x)) - 0.5*f(a) - 0.5*f(b)   # implicitly replace nan's with zero
    return h*s

def midpoint(f, a, b, n):
    """ Vectorized form of midpoint function """
    h = float(b-a)/n
    x = _np.linspace(a + h/2, b - h/2, n)
    return h*_np.nansum(f(x))   # implicitly replace nan's with zero

def midpoint_double(func, a, b, c, d, nx, ny):
    def gunc(x):
        return midpoint(lambda y: func(x, y), c, d, ny)
    return midpoint(gunc, a, b, nx)

def trapezoid_double(func, a, b, c, d, nx, ny):
    def gunc(x):
        return trapezoidal(lambda y: func(x, y), c, d, ny)
    return trapezoidal(gunc, a, b, nx)


#def _midpoint(f, a, b, n):
#    h = float(b-a)/n
#    result = 0
#    for i in range(n):
#        result += f((a + h/2.0) + i*h)
#    result *= h
#    return result
#
#def _midpoint_double1(func, a, b, c, d, nx, ny):
#    """ doubel integral of func"""
#    hx = (b - a)/float(nx)
#    hy = (d - c)/float(ny)
#    Integ = 0
#    for ii in range(nx):
#        for jj in range(ny):
#            xi = a + hx/2 + ii*hx
#            yj = c + hy/2 + jj*hy
#            Integ += hx*hy*func(xi, yj)
#    return Integ

def adaptive_integration(f, a, b, tol=1e-4, method='midpoint', verbose=False):
    n_limit = 1000000   # Just a choice (used to avoid inf loop)
    n = 2
    if method == 'trapezoidal':
        integral_n  = trapezoidal(f, a, b, n)
        integral_2n = trapezoidal(f, a, b, 2*n)
        diff = abs(integral_2n - integral_n)
        if verbose:
            print(('trapezoidal diff: ', diff))
        # end if
        while (diff > tol) and (n < n_limit):
            integral_n  = trapezoidal(f, a, b, n)
            integral_2n = trapezoidal(f, a, b, 2*n)
            diff = abs(integral_2n - integral_n)
            if verbose:
                print(('trapezoidal diff: ', diff))
            # end if
            n *= 2
        # end while
    elif method == 'midpoint':
        integral_n  = midpoint(f, a, b, n)
        integral_2n = midpoint(f, a, b, 2*n)
        diff = abs(integral_2n - integral_n)
        if verbose:
            print(('midpoint diff: ', diff))
        # end if
        while (diff > tol) and (n < n_limit):
            integral_n  = midpoint(f, a, b, n)
            integral_2n = midpoint(f, a, b, 2*n)
            diff = abs(integral_2n - integral_n)
            if verbose:
                print(('midpoint diff: ', diff))
            # end if
            n *= 2
    else:
        if verbose:
            print('Error - adaptive integration called with unknown par')
        # end if
    # end if
    # Now we check if acceptable n was found or not
    if diff <= tol:   # Success
        if verbose:
            print('The integral computes to: ', integral_2n)
        # end if
        return integral_2n, n
    else:
        return integral_2n, -n   # Return negative n to tell "not found"
    # end if
# end def

# ========================================================================== #
#   test section   / example section

def convergence_rates(f, F, a, b, num_experiments=14):
    expected = F(b) - F(a)
    n = _np.zeros(num_experiments, dtype=int)
    E = _np.zeros(num_experiments)
    r = _np.zeros(num_experiments-1)
    for i in range(num_experiments):
        n[i] = 2**(i+1)
        computed = trapezoidal(f, a, b, n[i])
        E[i] = abs(expected - computed)
        if i > 0:
            r_im1 = _np.log(E[i]/E[i-1])/_np.log(float(n[i])/n[i-1])
            r[i-1] = float('%.2f' % r_im1) # Truncate to two decimals
    return r

def test_midpoint():
    """Check difference between trapezoidal integration and midpoint method"""
    g = lambda y: _np.exp(-y**2.0)
    a = 0
    b = 2
    print('    n        midpoint          trapezoidal          abs.diff')
    for i in range(1, 21):
        n = 2**i
        m = midpoint(g, a, b, n)
        t = trapezoidal(g, a, b, n)
        print('%7d %.16f %.16f %.16f' % (n, m, t, t-m))
    # end for
# end def

def test_midpoint_double():
    """Test that a linear function is integrated exactly."""
    def f(x, y):
        return 2*x + y

    a = 0;  b = 2;  c = 2;  d = 3
    import sympy
    x, y = sympy.symbols('x  y')
    I_expected = sympy.integrate(f(x, y), (x, a, b), (y, c, d))
    # Test three cases: nx < ny, nx = ny, nx > ny
    for nx, ny in (3, 5), (4, 4), (5, 3):
        I_computed1 = midpoint_double(f, a, b, c, d, nx, ny)
        I_computed2 = trapezoid_double(f, a, b, c, d, nx, ny)
        tol = 1E-14
        #print I_expected, I_computed1, I_computed2
        assert abs(I_computed1 - I_expected) < tol
        assert abs(I_computed2 - I_expected) < tol
    # end for
# end def

def test_trapezoidal_linear():
    """Check that linear functions are integrated exactly."""
    f = lambda x: 6*x - 4
    F = lambda x: 3*x**2 - 4*x  # Anti-derivative
    a = 1.2; b = 4.4
    expected = F(b) - F(a)
    tol = 1E-14
    for n in [2, 20, 21]:
        computed = trapezoidal(f, a, b, n)
        error = abs(expected - computed)
        success = error < tol
        msg = 'n=%d, err=%g' % (n, error)
        assert success, msg
    # end for
# end def

def test_trapezoidal_conv_rate():
    """Check empirical convergence rates against the expected -2."""
    v = lambda t: 3*(t**2)*_np.exp(t**3)
    V = lambda t: _np.exp(t**3)
    a = 1.1; b = 1.9
    r = convergence_rates(v, V, a, b, 14)
    print(r)
    tol = 0.01
    msg = str(r[-4:])  # show last 4 estimated rates
    assert (abs(r[-1]) - 2) < tol, msg

# ========================================================================== #

def MonteCarlo_double(f, g, x0, x1, y0, y1, n):
    """
    Monte Carlo integration of f over a domain g>=0, embedded
    in a rectangle [x0,x1]x[y0,y1]. n^2 is the number of
    random points.
    """
    # Draw n**2 random points in the rectangle
    x = _np.random.uniform(x0, x1, n)
    y = _np.random.uniform(y0, y1, n)
    # Compute sum of f values inside the integration domain
    f_mean = 0
    num_inside = 0   # number of x,y points inside domain (g>=0)
    for i in range(len(x)):
        for j in range(len(y)):
            if g(x[i], y[j]) >= 0:
                num_inside += 1
                f_mean += f(x[i], y[j])
    f_mean = f_mean/float(num_inside)
    area = num_inside/float(n**2)*(x1 - x0)*(y1 - y0)
    return area*f_mean

def test_MonteCarlo_double_rectangle_area():
    """Check the area of a rectangle."""
    def g(x, y):
        return (1 if (0 <= x <= 2 and 3 <= y <= 4.5) else -1)

    x0 = 0;  x1 = 3;  y0 = 2;  y1 = 5  # embedded rectangle
    n = 1000
    _np.random.seed(8)      # must fix the seed!
    I_expected = 3.121092  # computed with this seed
    I_computed = MonteCarlo_double(
        lambda x, y: 1, g, x0, x1, y0, y1, n)
    assert abs(I_expected - I_computed) < 1E-14

def test_MonteCarlo_double_circle_r():
    """Check the integral of r over a circle with radius 2."""
    def g(x, y):
        xc, yc = 0, 0  # center
        R = 2          # radius
        return  R**2 - ((x-xc)**2 + (y-yc)**2)

    # Exact: integral of r*r*dr over circle with radius R becomes
    # 2*pi*1/3*R**3
    import sympy
    r = sympy.symbols('r')
    I_exact = sympy.integrate(2*sympy.pi*r*r, (r, 0, 2))
    print(('Exact integral:', I_exact.evalf()))
    x0 = -2;  x1 = 2;  y0 = -2;  y1 = 2
    n = 1000
    _np.random.seed(6)
    I_expected = 16.7970837117376384  # Computed with this seed
    I_computed = MonteCarlo_double(
        lambda x, y: _np.sqrt(x**2 + y**2),
        g, x0, x1, y0, y1, n)
    print('MC approximation %d samples: %.16f' % (n**2, I_computed))
    assert abs(I_expected - I_computed) < 1E-15

# ========================================================================== #
# ========================================================================== #


if __name__ == '__main__':
    import doctest
    doctest.testmod()

#    import os as _os
#    if _os.path.exists('utils.testfile'):
#        _os.remove('utils.testfile')

    test_midpoint()
    test_midpoint_double()
    test_trapezoidal_linear()
    test_trapezoidal_conv_rate()
    test_MonteCarlo_double_rectangle_area()
    test_MonteCarlo_double_circle_r()
# end if

# ========================================================================== #
# ========================================================================== #

