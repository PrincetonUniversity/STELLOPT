# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 13:13:58 2019

@author: Gavin Weir <gavin.weir@­ipp.mpg.de>
"""
# ======================================================================== #
# ======================================================================== #

from __future__ import absolute_import, with_statement, absolute_import, \
                       division, print_function, unicode_literals

import numpy as _np

# ======================================================================= #
# ======================================================================= #


class Spline(object):
    def __init__(self, xvar, yvar, xf=None, vary=None, bbox=None, **kwargs):
        self.vary = _np.copy(vary)
        self.yvar = _np.copy(yvar)
        self.xvar = _np.copy(xvar)

#        kwargs.setdefault('n',len(self.xvar))
        kwargs.setdefault('end1', int(0))
        kwargs.setdefault('end2', int(0))
        kwargs.setdefault('slope1', int(0))
        kwargs.setdefault('slope2', int(0))
        self.nmonti = kwargs.pop('nmonti', None)
        self.kwargs = kwargs

        if xf is None:              xf = xvar.copy()    # endif
        self.xf = _np.copy(xf)
        self.bbox = bbox

        self.deriv_out=False
    # end def __init__

    # ========================= #
    # ========================= #

    def derivative(self, deriv_out=True):
        self.deriv_out = deriv_out

    def __call__(self, x, nu=0):
        x = _np.asarray(x)
        # empty input yields empty output
        if x.size == 0:
            return _np.array([])

        yf, dydx = seval(x, self.xvar, self.yvar, self.coeffs[0], self.coeffs[1], self.coeffs[2])
        if nu == 0:
            return yf
        elif nu == 1 or self.deriv_out:
            return dydx
        # end if
    # end def

    def spline(self):
        self.status, b, c, d = _spline(self.xvar[self.msk], self.yvar[self.msk], **self.kwargs)

        if self.status != 0:
            print(self.status)
        coeffs = []
        coeffs.append(b)
        coeffs.append(c)
        coeffs.append(d)
        self.coeffs = coeffs

        self.yf, self.dydx = seval(self.xf, self.xvar, self.yvar,
                  self.coeffs[0], self.coeffs[1], self.coeffs[2])
#        return self.yf, self.dydx
        return self
    # end def

    def spline_and_eval(self):
        if self.nmonti is None:
#            self.status, b, c, d = _spline(self.xvar[self.msk], self.yvar[self.msk], **self.kwargs)
#
#            if self.status != 0:
#                print(self.status)
#            coeffs = []
#            coeffs.append(b)
#            coeffs.append(c)
#            coeffs.append(d)
#            self.coeffs = coeffs
#
#            self.yf, self.dydx = seval(self.xf, self.xvar, self.yvar,
#                      self.coeffs[0], self.coeffs[1], self.coeffs[2])
            self.spline()
            return self.yf, self.dydx
        else:
            self.yf, self.varf, self.dydx, self.vardydx = self.spline_and_eval_bs()
            return self.yf, self.varf, self.dydx, self.vardydx
    # end def

    # ========================= #


    def spline_and_eval_bs(self):
        nmonti = self.nmonti
        self.nmonti = None

        # ============= #

        nxf = len(self.xf)
        ysave = _np.copy(self.yvar)
        vsave = _np.copy(self.vary)

        yvar = _np.atleast_2d(self.yvar).copy()
        vary = _np.atleast_2d(self.vary).copy()
        nsh = _np.shape(yvar)
        nsh = _np.atleast_1d(nsh)
        if nsh[0] == 1:
            yvar = yvar.T
            vary = vary.T
            nsh = _np.flipud(nsh)
        # end if

        # ============= #

#        xvar = self.xvar.copy()
#        xvar = xvar.reshape((nsh[0],), order='C')
#        niterate = len(xvar)
#        niterate *= nmonti
        niterate = nmonti

        dydx = _np.zeros( (niterate, nxf, nsh[1]), dtype=_np.float64)
        yf = _np.zeros( (niterate, nxf, nsh[1]), dtype=_np.float64)
#        cc = -1
        for ii in range(niterate):
            utemp = yvar.copy()
#            vtemp = vary.copy()
            utemp += _np.sqrt(vary)*_np.random.normal(0.0, 1.0, _np.shape(vary))
#            cc += 1
#            if cc >= nsh[0]:
#                cc = 0
#            # end if
#            utemp[cc,:] += _np.sqrt(vary[cc,:])*_np.random.normal(0.0, 1.0, _np.shape(vary[cc,:]))
#            vtemp[cc,:] = (utemp[cc,:]-yvar[cc,:])**2.0

            tmp1 = _np.zeros((nxf, nsh[1]), dtype=utemp.dtype)
            tmp2 = _np.zeros_like(tmp1)
            for jj in range(nsh[1]):
                self.yvar = utemp[:,jj].copy()
                tmp1[:,jj], tmp2[:,jj] = self.spline_and_eval()
            # end for
            yf[ii, :] = tmp1.reshape((nxf, nsh[1]), order='C')
            dydx[ii, :] = tmp2.reshape((nxf, nsh[1]), order='C')
        # end for
        vardydx = _np.var(dydx, axis=0)
        dydx = _np.mean(dydx, axis=0)
        varf = _np.var( yf, axis=0)
        yf = _np.mean(yf, axis=0)
        self.nmonti = nmonti
        self.yvar = ysave
        self.vary = vsave
        return yf, varf, dydx, vardydx

    # ========================= #


    @property
    def bbox(self):
        return self._bbox
    @bbox.setter
    def bbox(self, listvals):
        self._bbox = listvals
        msk = _np.ones((len(self.xvar),), dtype=bool)
        if self._bbox is not None:
            cond = (self.xvar>=self._bbox[0])*(self.xvar<=self._bbox[1])
            if (cond == False).any():
                msk[cond] = False
            # end if
        # end if
        self.msk = msk
#        self.spline()
    @bbox.deleter
    def bbox(self):
        del self._bbox
        del self.msk

    # ========================= #

    def show(self):
        import matplotlib.pyplot as _plt
        _plt.figure()
        _ax1 = _plt.subplot(2,1,1)
        _ax2 = _plt.subplot(2,1,2, sharex=_ax1)
        _ax1.plot(self.xvar, self.yvar, 'ko')
    #    _ax2.plot(self.x, self.dydx, 'ko')

        _ax1.plot(self.xf, self.yf, 'b-')
        _ax1.set_ylabel('function')

        _ax2.plot(self.xf, self.yf, 'b-')
        _ax2.set_xlabel('x')
        _ax2.set_ylabel('deriv')
        _ax1.set_title('Spline functions')
    # ========================= #
# end class

# ======================================================================== #


def spline_bs(xvar, yvar, vary, xf=None, nmonti=300, bbox=None, **kwargs):
    if xf is None:              xf = xvar.copy()    # endif
    if nmonti is None:          nmonti = 300        # endif
    nxf = len(xf)

    if xf is None:
        xf = xvar.copy()
    # endif
    nxf = len(xf)

    # ============= #

    yvar = _np.atleast_2d(yvar)
    vary = _np.atleast_2d(vary)
    nsh = _np.shape(yvar)
    nsh = _np.atleast_1d(nsh)
    if nsh[0] == 1:
        yvar = yvar.T
        vary = vary.T
        nsh = _np.flipud(nsh)
    # end if

    # ============= #

    xvar = xvar.reshape((nsh[0],), order='C')

    niterate = len(xvar)
    niterate *= nmonti

    dydx = _np.zeros( (niterate, nxf, nsh[1]), dtype=_np.float64)
    yf = _np.zeros( (niterate, nxf, nsh[1]), dtype=_np.float64)
    cc = -1
    for ii in range(niterate):
        utemp = yvar.copy()
        vtemp = vary.copy()
        cc += 1
        if cc >= nsh[0]:
            cc = 0
        # end if
        utemp[cc,:] += _np.sqrt(vary[cc,:])*_np.random.normal(0.0, 1.0, _np.shape(vary[cc,:]))
        vtemp[cc,:] = (utemp[cc,:]-yvar[cc,:])**2.0

        tmp1 = _np.zeros((nxf, nsh[1]), dtype=utemp.dtype)
        tmp2 = _np.zeros_like(tmp1)
        for jj in range(nsh[1]):
            SplObj = Spline(xvar, utemp[:,jj], xf, bbox, **kwargs)
            tmp1[:,jj], tmp2[:,jj] = SplObj.spline()

        # end for
            yf[ii, :] = tmp1.reshape((nxf, nsh[1]), order='C')
            dydx[ii, :] = tmp2.reshape((nxf, nsh[1]), order='C')
    # end for

    vardydx = _np.var(dydx, axis=0)
    dydx = _np.mean(dydx, axis=0)
    varf = _np.var( yf, axis=0)
    yf = _np.mean(yf, axis=0)

    return yf, varf, dydx, vardydx

# ======================================================================= #
# ======================================================================= #

#def _spline(x, y, n=int(2), end1=int(0), end2=int(0),
def _spline(x, y, end1=int(0), end2=int(0), slope1=_np.float64(0), slope2=_np.float64(0)):
    """
    *************************************************************
     Evaluate the coefficients for a cubic interpolating spline
    !!Attention
       remove lines bracketed by comments //CStConfig
       if you want to use this function in other applications
    http://www.mech.uq.edu.au/staff/jacobs/nm_lib/cmathsrc/spline.c

    * Evaluate the coefficients b[i], c[i], d[i], i = 0, 1, .. n-1 for
       a cubic interpolating spline

       S(xx) = Y[i] + b[i] * w + c[i] * w**2 + d[i] * w**3
       where w = xx - x[i] and  x[i] <= xx <= x[i+1]

       The n supplied data points are x[i], y[i], i = 0 ... n-1.

    Input :
       n       : The number of data points or knots (n >= 2)
       end1,
       end2    : = 1 to specify the slopes at the end points
                 = 0 to obtain the default conditions
       slope1,
       slope2  : the slopes at the end points x[0] and x[n-1]
                 respectively
       x[]     : the abscissas of the knots in increasing order
       y[]     : the ordinates of the knots

     Output :
       b, c, d : arrays of spline coefficients as defined above
                 (See note 2 for a definition.)
     Return:
                = 0 normal return
                = 1 less than two data points; cannot interpolate
                = 2 x[] are not in ascending order

       This C code written by ...  Peter & Nigel,
       ----------------------      Design Software,
                                   42 Gubberley St,
                                   Kenmore, 4069,
                                   Australia.

       Version ... 1.1, 30 September 1987
       -------     2.0, 6 April 1989    (start with zero subscript)
                                         remove ndim from parameter list
                   2.1, 28 April 1989   (check on x[])
                   2.2, 10 Oct   1989   change number order of matrix
                   2.3, 24 Jun 2019 - Translated to python by GMW

       Notes ...
       -----
       (1) The accompanying function seval() may be used to evaluate the
           spline while deriv will provide the first derivative.
       (2) Using p to denote differentiation
           y[i] = S(X[i])
           b[i] = Sp(X[i])
           c[i] = Spp(X[i])/2
           d[i] = Sppp(X[i])/6  ( Derivative from the right )
       (3) Since the zero elements of the arrays ARE NOW used here,
           all arrays to be passed from the main program should be
           dimensioned at least [n].  These routines will use elements
           [0 .. n-1].
       (4) Adapted from the text
           Forsythe, G.E., Malcolm, M.A. and Moler, C.B. (1977)
           "Computer Methods for Mathematical Computations"
           Prentice Hall
       (5) Note that although there are only n-1 polynomial segments,
           n elements are requird in b, c, d.  The elements b[n-1],
           c[n-1] and d[n-1] are set to continue the last segment
           past x[n-1].
    *----------------------------------------------------------------*/
    """
    # {
    n = len(x)
    nm1 = n-1
    if(n<2):
        status = 1     # // no possible interpolation
        return (status, [], [], [])
    for ii in range(1,n):
        # print((ii, n))
        if (x[ii] <= x[ii-1]):
            status = 2  # //x[] are not in ascending order
            return (status, [], [], [])
        # end if
    # end for

    if n==2:  # /* linear segment only  */
        b = _np.zeros((2,), dtype=_np.float64)
        c = _np.zeros_like(b)
        d = _np.zeros_like(b)
        b[0]=(y[1]-y[0])/(x[1]-x[0])
        c[0]=0
        d[0]=0
        b[1] = _np.copy(b[0])
        c[1] = _np.copy(c[0])
        d[1] = _np.copy(d[0])
        status = 0
        return status, b, c, d
    # end if

    # //CStConfig-beg
    if end1 == 0:
       slope1 = (y[1]-y[0])/(x[1]-x[0])
    if end2 == 0:
       slope2 = (y[n-1]-y[n-2])/(x[n-1]-x[n-2])
    # end if //CStConfig-end

    # //    Set up the symmetric tri-diagonal system
    # //  b = diagonal d = offdiagonal c = right-hand-side
    b = _np.zeros((n,), dtype=_np.float64)
    c = _np.zeros_like(b)
    d = _np.zeros_like(b)

    d[0] = x[1] - x[0]
    c[1] = (y[1] - y[0]) / d[0]
    for ii in range(1, nm1):
        d[ii]   = x[ii+1] - x[ii]
        b[ii]   = 2*(d[ii-1]+d[ii])
        c[ii+1] = (y[ii+1]-y[ii])/d[ii]
        c[ii]   = c[ii+1]-c[ii]
    # end for

    # /* ---- Default End conditions
    #      Third derivatives at x[0] and x[n-1] obtained
    #      from divided differences  */
    b[0]   = -d[0]
    b[nm1] = -d[n-2]
    c[0]   = 0
    c[nm1] = 0
    if (n != 3): #{
        c[0]   = c[2]/(x[3]-x[1]) - c[1]/(x[2]-x[0])
        c[nm1] = c[n-2]/(x[nm1]-x[n-3]) - c[n-3]/(x[n-2]-x[n-4])
        c[0]   = c[0]*d[0]*d[0]/(x[3]-x[0])
        c[nm1] = -c[nm1]*d[n-2]*d[n-2] / (x[nm1]-x[n-4])
    # end if }
    # /* Alternative end conditions -- known slopes */
    if (end1): # {
        b[0] = 2*(x[1]-x[0])
        c[0] = (y[1]-y[0])/(x[1]-x[0]) - slope1
    # end if  }
    if (end2): # {
        b[nm1] = 2*(x[nm1]-x[n-2])
        c[nm1] = slope2 - (y[nm1]-y[n-2])/(x[nm1]-x[n-2])
    # end if }

    # /* Forward elimination */
    for ii in range(1,n): #(i = 1; i < n; ++i) {
        t = (d[ii-1]/b[ii-1])
        b[ii] = b[ii]-t*d[ii-1]
        c[ii] = c[ii]-t*c[ii-1]
    # end for  }

    #  /* Back substitution */
    c[nm1] = c[nm1]/b[nm1]
    for ii in range(n-2, 0-1, -1): # (i=n-2; i>=0; i--)
        c[ii] = (c[ii]-d[ii]*c[ii+1])/b[ii]
    # end for

    # /* c[i] is now the sigma[i] of the text */
    # /* Compute the polynomial coefficients */
    b[nm1] = (y[nm1]-y[n-2])/d[n-2]+d[n-2]*(c[n-2]+2*c[nm1])
    for ii in range(0, nm1): # (i=0; i<nm1; ++i) {
        b[ii] = (y[ii+1]-y[ii])/d[ii] - d[ii]*(c[ii+1]+2*c[ii])
        d[ii] = (c[ii+1]-c[ii])/d[ii]
        c[ii] = 3*c[ii]
    # } end for
    c[nm1] = 3*c[nm1]
    d[nm1] = d[n-2]
    status = 0

#    # //CStConfig-beg
#    ii=n-2
#    b[ii] = slope2
#    d[ii] = c[ii] = 0
#    ii=n-1
#    b[ii] = slope2
#    d[ii] = c[ii] = 0
#    #//CStConfig-end
#    status = 0;
    # }
    return status, b, c, d

# ==========================================================================#
# //***************************************************************************
# //

def seval(xx, x, y, b, c, d):
    xx = _np.atleast_1d(xx)
    nx = len(xx)
    ys = _np.zeros((nx,), dtype=_np.float64)
    dy = _np.zeros_like(ys)
    for ii in range(nx):
        ys[ii], dy[ii] = _seval(xx[ii], x, y, b, c, d)
    # end for
    return ys, dy

def _seval(u, x, y, b, c, d):
    """
    double CStconfig::seval (int n, double u,
                  double x[], double y[],
                  double b[], double c[], double d[],
                  double &yp)

    /*
      Evaluate the cubic spline function
      and the derivative of the cubic spline function

      S(u)  = y[i]+b[i]*w+c[i]*w^2+d[i]*w^3
      S'(u) = b[i]+2*c[i]*w + 3*d[i]*w**2
      where w = u - x[i]
      and   x[i] <= u <= x[i+1]

      If u < x[0]   then i = 0 is used.
      If u > x[n-1] then i = n-1 is used.

      Input :
      n       : The number of data points or knots (n >= 2)
      u       : the abscissa at which the spline is to be evaluated
      x[]     : the abscissas of the knots in increasing order
      y[]     : the ordinates of the knots
      b, c, d : arrays of spline coefficients computed by spline().

      Return:
      seval   : the value of the spline function at u
      yp      : the derivative of the cubic spline function
    */
    """
    # {
    n = len(x)
    i = int(0)        # //Left bracket
    j = int(n)        # //Right bracket
    if (u<=x[0]):
        j=1
    elif (u>=x[j-1]):
        i=j-1
    else:
        while (i+1 < j): # {            //bisection
            k = int((i+j)/2)
            if (u < x[k]):
                j = k
            else:
                i = k
        # } end while
    # end if

    dx = u - x[i]
    s  = y[i] + dx*(b[i]+dx*(c[i]+dx*d[i]))  # // evaluate the spline
    yp = b[i] + dx*(2*c[i]+dx*3*d[i])        #  // evaluate the derivative
    return (s, yp)
    # }

def test_spline_bs(nsegs=None):
    # nsegs will be the number of knot positions in the spline
    if nsegs is None: nsegs = 11  # end if

    xx = _np.linspace(-0.1, 1.3, num=nsegs, endpoint=True)
    yy = 3.3*_np.exp(0.5*(xx-0.6)*(xx-0.6)/0.25/0.25)
    vy = (0.1*_np.random.normal(0.0, 1.0, len(yy))*yy)**2.0
    vy[vy<(0.1*yy.max())**2.0] = (0.1*yy.max())**2.0

    xgrid = _np.linspace(xx.min(), xx.max(), num=100, endpoint=True)

    # first do not specify the endpoint slopes
    kwargs = {}
#    kwargs.setdefault('n',len(xgrid))
    kwargs.setdefault('end1', int(0))
    kwargs.setdefault('end2', int(0))
    Spl1Obj = Spline(xx, yy, xgrid, vy, nmonti=300, **kwargs)
    yspl1, vspl1, dy1dx, vdy1dx = Spl1Obj.spline_and_eval()

    # now specify zero slope at end points
    kwargs['end1'] = int(1)
    kwargs['end2'] = int(1)
    kwargs.setdefault('slope1', int(0))
    kwargs.setdefault('slope2', int(0))
    Spl2Obj = Spline(xx, yy, xgrid, vy, nmonti=300, **kwargs)
    yspl2, vspl2, dy2dx, vdy2dx = Spl2Obj.spline_and_eval()

    # derivative by simple finite differences
    dydx = _np.hstack((_np.atleast_1d(0.5*(yy[1]-2.0*yy[0])/(xx[1]-xx[0])), _np.diff(yy)/_np.diff(xx)))
    vardydx = _np.hstack((_np.atleast_1d((0.25*vy[1]+vy[0])*dydx[0]**2.0), (vy[1:]+vy[:-1])*dydx[1:]**2.0))

    import matplotlib.pyplot as _plt
    _plt.figure()
    _ax1 = _plt.subplot(2,1,1)
    _ax2 = _plt.subplot(2,1,2, sharex=_ax1)
    _ax1.errorbar(xx, yy, yerr=_np.sqrt(vy), fmt='ko')
    _ax2.errorbar(xx, dydx, yerr=_np.sqrt(vardydx), fmt='ko')

    _ax1.errorbar(xgrid, yspl1, yerr=_np.sqrt(vspl1), fmt='b-')
    _ax1.errorbar(xgrid, yspl2, yerr=_np.sqrt(vspl2), fmt='r-')
    _ax1.set_ylabel('function')

    _ax2.errorbar(xgrid, dy1dx, yerr=_np.sqrt(vdy1dx), fmt='b-')
    _ax2.errorbar(xgrid, dy2dx, yerr=_np.sqrt(vdy2dx), fmt='r-')
    _ax2.set_xlabel('x')
    _ax2.set_ylabel('deriv')
    _ax1.set_title('Testing spline functions')
    return xx, yy, yspl1, yspl2

def test_spline(nsegs=None):
    # nsegs will be the number of knot positions in the spline
    if nsegs is None: nsegs = 11  # end if

    xx = _np.linspace(-0.1, 1.3, num=nsegs, endpoint=True)
    yy = 3.3*_np.exp(0.5*(xx-0.6)*(xx-0.6)/0.25/0.25)

    # first do not specify the endpoint slopes
    # status1, b1, c1, d1 = _spline(xx, yy, n=len(xx))
    status1, b1, c1, d1 = _spline(xx, yy)

    # now specify zero slope at end points
#    status2, b2, c2, d2 = _spline(xx, yy, n=len(xx), end1=int(1), end2=int(1),
    status2, b2, c2, d2 = _spline(xx, yy, end1=int(1), end2=int(1),
            slope1=_np.float64(0), slope2=_np.float64(0))

    xgrid = _np.linspace(xx.min(), xx.max(), num=100, endpoint=True)
    yspl1, dy1dx = seval(xgrid, xx, yy, b1, c1, d1)
    yspl2, dy2dx = seval(xgrid, xx, yy, b2, c2, d2)

    xch = _np.linspace(-0.5, 1.5, num=20, endpoint=True)
    yspl3, dy3dx = seval(xch, xx, yy, b1, c1, d1)

    # derivative by simple finite differences
    dydx = _np.hstack((_np.atleast_1d(0.5*(yy[1]-2.0*yy[0])/(xx[1]-xx[0])), _np.diff(yy)/_np.diff(xx)))

    import matplotlib.pyplot as _plt
    _plt.figure()
    _ax1 = _plt.subplot(2,1,1)
    _ax2 = _plt.subplot(2,1,2, sharex=_ax1)
    _ax1.plot(xx, yy, 'ko')
    _ax2.plot(xx, dydx, 'ko')

    _ax1.plot(xgrid, yspl1, 'b-')
    _ax1.plot(xgrid, yspl2, 'r-')
    _ax1.plot(xch, yspl3, 'm-')
    _ax1.set_ylabel('function')

    _ax2.plot(xgrid, dy1dx, 'b-')
    _ax2.plot(xgrid, dy2dx, 'r-')
    _ax2.plot(xch, dy3dx, 'm-')
    _ax2.set_xlabel('x')
    _ax2.set_ylabel('deriv')
    _ax1.set_title('Testing spline functions')
    return xx, yy, yspl1, yspl2

if __name__=="__main__":
    # test it
    test_spline()
    test_spline_bs()
# end if