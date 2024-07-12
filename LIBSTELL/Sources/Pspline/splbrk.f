C--------------------------------------------------------------------
C  SPLBRK -- make a spline with a break (C0 only) at specified locn
C
      SUBROUTINE SPLBRK (IOPT, N, NBRK, X, Y, B, C, D)
C
C  Spline the whole interval, and then respline a subinterval, saving
C  the endpt coeffs from the original spline.
C
C  Result is a spline that is C2 everywhere except C0 only at one
C  interior break point.
C
C  IOPT=0  use SPLAAN for dy/dx=0 at LHS bc
C  IOPT=1  use SPLINE for standard bc
C
C  N -- no. of data pts
C  NBRK -- break point
C
      real X(N),Y(N),B(N),C(N),D(N)
C
      if(iopt.eq.0) then
         call SPLAAN(n, x, y, b, c, d)
      else
         call SPLINE_ez1(n, x, y, b, c, d)
      endif
C
      bsave=b(nbrk)
      csave=c(nbrk)
      dsave=d(nbrk)
C
      call nspline(nbrk, x, y, b, c, d)
C
      b(nbrk)=bsave
      c(nbrk)=csave
      d(nbrk)=dsave
C
      return
      end
