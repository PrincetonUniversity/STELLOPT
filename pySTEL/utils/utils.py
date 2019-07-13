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

    fmt2 = ''.join(['%g'+delim for _ in range(2)]) # '%g%g'
    fmt3 = ''.join(['%g'+delim for _ in range(3)]) # '%g%g%g'
    fmt6 = ''.join(['%g'+delim for _ in range(6)]) # '%g%g%g%g%g%g'
    fmt7 = ''.join(['%g'+delim for _ in range(7)]) # '%g%g%g%g%g%g%g'
    fmt10 = ''.join(['%g'+delim for _ in range(10)]) # '%g%g%g%g%g%g%g%g%g%g'
    fmt11 = ''.join(['%g'+delim for _ in range(11)]) # '%g%g%g%g%g%g%g%g%g%g%g'
    fmt12 = ''.join(['%g'+delim for _ in range(12)]) # '%g%g%g%g%g%g%g%g%g%g%g%g'
    fmt13 = ''.join(['%g'+delim for _ in range(13)]) # '%g%g%g%g%g%g%g%g%g%g%g%g%g'
    fmt14 = ''.join(['%g'+delim for _ in range(14)]) # '%g%g%g%g%g%g%g%g%g%g%g%g%g%g'
    fmt20 = ''.join(['%g'+delim for _ in range(20)]) # '%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g'
    fmt2_3 = (''.join(['%d'+delim for _ in range(2)]) # '%d%d%g%g%g'
           + ''.join(['%g'+delim for _ in range(3)]) )
    fmt2_6 = (''.join(['%d'+delim for _ in range(2)]) # '%d%d%g%g%g%g%g%g'
           + ''.join(['%g'+delim for _ in range(6)]) )
    fmt2_7 = (''.join(['%d'+delim for _ in range(2)]) # '%d%d%g%g%g%g%g%g%g'
            + ''.join(['%g'+delim for _ in range(7)]) )
    fmt2_11 = (''.join(['%d'+delim for _ in range(2)]) # '%d%d%g%g%g%g%g%g%g%g%g%g%g'
           + ''.join(['%g'+delim for _ in range(11)]) )
    fmt2_14 = (''.join(['%d'+delim for _ in range(2)]) # '%d%d%g%g%g%g%g%g%g%g%g%g%g%g%g%g'
             + ''.join(['%g'+delim for _ in range(14)]) )
    fmt2_3_2_7 = (''.join(['%d'+delim for _ in range(2)]) # '%d%d%g%g%g%d%d%g%g%g%g%g%g%g'
                + ''.join(['%g'+delim for _ in range(3)])
                + ''.join(['%d'+delim for _ in range(2)])
                + ''.join(['%g'+delim for _ in range(7)]) )
    fmt2_6_2_14 = (''.join(['%d'+delim for _ in range(2)]) # '%d%d%g%g%g%g%g%g%d%d%g%g%g%g%g%g%g%g%g%g%g%g%g%g'
                + ''.join(['%g'+delim for _ in range(6)])
                + ''.join(['%d'+delim for _ in range(2)])
                + ''.join(['%g'+delim for _ in range(14)]) )


def ftell(fid):
    """
        A python function that mimics ftell from MATLAB
            returns the current file position
    """
    return fid.tell()

def fseek(fid, offset, origin='currentposition'):
    """
        A python function that mimics fseek from MATLAB
            sets the current file position
    inputs:
        fid    - file object identifier
        offset - integer number of bytes to move the current file position
        origin - reference from where to offset the file position

    In text files only seeks relative to the beginning of the file are allowed
    (except for seeking to the very end of the file with seek(0, 2))

    The only valid offset values are those returned from fid.tell or 0
    """
    if type(origin) == type(''):
        if origin.lower().find('currentposition')>-1:   # default
            origin = 1
        elif origin.lower().find('start-of-file')>-1:
            origin = 0
        elif origin.lower().find('end-of-file')>-1:
            origin = 2
        # end if
    # end if
    return fid.seek(offset, origin)

def frewind(fid):
    return fseek(fid, 0, 'start-of-file')

def fgets(fid):
    """
        A python function that mimics fgets from MATLAB
            returns the next line of a file while including the new line character
    """
    return fid.readline()

def fgetl(fid):
    """
        A python function that mimics fgetl from MATLAB
            returns the next line of a file while stripping the new line character
    """
    return fgets(fid).replace('\n','')

def fscanf(fid, dtypes='#g', size=None, offset=None):
    """
    This function tries to mimic fscanf from MATLAB. In MATLAB, fscanf
    reads in formatted data from an open file up to a specified shape.
        fscanf(fileID,formatSpec,sizeA)

    Input:
        fid    - file id generated from fid = open(filename), etc.
        dtypes - data types mimicing the formatSpec keyword in MATLAB's fscanf function
            See the note below
        size   - 1D or 2D shape of returned data, [nrow, ncol] or an integer
            The shape of the output data
            if size = 13, then 13 elements will be read in along a single row
            if size = [3, 2], then 2 elements will be read in for three rows
        offset - The starting point to read from the file
            MATLAB automatically keeps track.
        # delimiter - defaults to white space

    Outputs:
        data - shape is equal to input size, but defaults to (1,length of the row)

    MATLAB cheat-sheet for fscanf

    Formatting
        Integer, signed
            %d - Base 10
            %i - The values in the file determine the base:
                    The default is base 10.
                    If the initial digits are 0x or 0X, then the values are hexadecimal (base 16).
                    If the initial digit is 0, then values are octal (base 8).
            %ld or %li - 64-bit values, base 10, 8, or 16
        Integer, unsigned
            %u - Base 10
            %o - Base 8 (octal)
            %x - Base 16 (hexadecimal)
            %lu, %lo, %lx - 64-bit values, base 10, 8, or 16
        Floating-point number
            %f, %e, %g - Floating-point fields can contain any of the following
                (not case sensitive): Inf, -Inf, NaN, or -NaN.

        String
            %s - string - Read all characters excluding white spaces
            %c - string - Read any single character, including white space.
                    To read multiple characters at a time, specify field width.

    """
#    nptypes = {'i':_np.int, 'd':_np.int, 'ld':_np.int, 'li':_np.int,  # TODO:! actually use correct datatypes
#               'u':_np.uint, 'o':_np.uint, 'x':_np.uint, 'lu':_np.uint, 'lo':_np.uint, 'lx':_np.uint,
#               'f':_np.float, 'e':_np.float, 'g':_np.float,
#               's':_np.str, 'c':_np.str,
#               'bool':_np.bool}
    nptypes = {'i':'int', 'd':'int', 'ld':'int', 'li':'int',  # TODO:! actually use correct datatypes
               'u':'uint', 'o':'uint', 'x':'uint', 'lu':'uint', 'lo':'uint', 'lx':'uint',
               'f':'float', 'e':'float', 'g':'float',
               's':'str', 'c':'str',
               'bool':'bool'}
    # Convert the different datatypes into a list
    delimiter = ' '
    if dtypes.find(',')>-1:   # delimiter hidden in fmt string
        delimiter = ','
        dtypes = dtypes.replace(',', '')
    # end if

    dtypes = list(_np.atleast_1d(dtypes.split('%')))[1:]
    ndt = len(dtypes)
    # ndt = len(dtypes.split('#'))-1

    nrow = 1
    ncol = ndt

    if size is not None:
        size = _np.atleast_1d(size)
        if len(size)>1:
            nrow = size[0]
            ncol = size[1]
        else:
            ncol = size[0]
        # end if
    # end if

#    string_line = False
    if ncol != ndt:
        # Assume that one data type was put in for the entire line regardless of how many columns
#        if len(dtypes)==1 and dtypes[0] =='s':
#            string_line = True
#        else:
        if 1:
#            print(len(dtypes))
            dtypes = dtypes*ncol
#            print(len(dtypes))
            ndt = ncol
        # end if
    # end if

    if offset is not None:
        fid.seek(offset, 0)
    # end if

    dt = _np.dtype(''.join(nptypes[dtypes[jj]] for jj in range(ncol)))
    data = _np.zeros((nrow, ncol), dtype=dt)
#    data = []
    for ii in range(nrow):
        try:
            line = fgetl(fid)
            if len(line)==0:
                continue # end of file has been reached or white space with readline
        except EOFError:
            break
        # end try

        if size == 1 and dtypes[ii]=='s':
            return line
        # end if

        # split up the line into a list (choose where by delimiter)
        line = list(line.split(delimiter))  # now a tuple, now a list

        # Remove empty elements from the list
        line = [line[jj] for jj in range(len(line)) if line[jj] != '']

        data[ii, :] = tuple(line)
#        # Loop over the elements of the list and assign the data type
#        for jj in range(len(line)):
#            line[jj] = nptypes[dtypes[jj]](line[jj])
#        # end for
#        data.append(line)
    # end if
    return data
#    return data, fid.tell()

# ===================================================================== #


# copied from libstell.py by SAL to circumvent libstell import here
def h2f(var,ns):
    temp = _np.zeros((ns,1))
    temp[0] = 1.5 * var[0] - 0.5 * var[1]
    temp[1:ns-1] = 0.5* (var[0:ns-2] + var[1:ns-1])
    temp[ns-1] = 1.5 * var[ns-2] - 0.5 * var[ns-3]
    return temp

def cfunct(theta,zeta,fmnc,xm,xn):
    f=0
    (ns,mn)=fmnc.shape
    lt = len(theta)
    lz = len(zeta)
    mt=_np.matmul(xm,theta.T)
    nz=_np.matmul(xn,zeta.T)
    cosmt=_np.cos(mt)
    sinmt=_np.sin(mt)
    cosnz=_np.cos(nz)
    sinnz=_np.sin(nz)
    f = _np.zeros((ns,lt,lz))

    fmn = _np.ndarray((mn,lt))
    for k in range(ns):
        fmn = _np.broadcast_to(fmnc[k,:],(lt,mn)).T
        fmncosmt=(fmn*cosmt).T
        fmnsinmt=(fmn*sinmt).T
        f[k,:,:]=_np.matmul(fmncosmt, cosnz)-_np.matmul(fmnsinmt, sinnz)
    return f

def sfunct(theta,zeta,fmnc,xm,xn):
    f=0
    (ns,mn)=fmnc.shape
    lt = len(theta)
    lz = len(zeta)
    mt=_np.matmul(xm,theta.T)
    nz=_np.matmul(xn,zeta.T)
    cosmt=_np.cos(mt)
    sinmt=_np.sin(mt)
    cosnz=_np.cos(nz)
    sinnz=_np.sin(nz)
    f = _np.zeros((ns,lt,lz))
    fmn = _np.ndarray((mn,lt))
    for k in range(ns):
        fmn = _np.broadcast_to(fmnc[k,:],(lt,mn)).T
        f[k,:,:]=_np.matmul((fmn*sinmt).T,cosnz)+_np.matmul((fmn*cosmt).T,sinnz)
    return f

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


if __name__ == '__main__':
    import doctest
    doctest.testmod()
# end if

# ========================================================================== #
