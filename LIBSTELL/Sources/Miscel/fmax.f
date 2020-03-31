      FUNCTION fmax (ax, bx, f, tol)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec) ax, bx, f, tol
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero=0, one=1, two=2,
     1   three = 3, five = 5, p5 = 0.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec):: a,b,c,d,e,eps,xm,p,q,r,tol1,tol2,
     1   u,v,w,fu,fv,fw,fx,x, fmax
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL f
C-----------------------------------------------
c
c  c is the squared inverse of the golden ratio
c
      c = p5*(three - SQRT(five))
c
c  eps is approximately the square root of the relative machine
c  precision.
c
      eps = one
      eps = eps/two
      tol1 = one + eps
 1003 CONTINUE
      IF (tol1 .le. one) GOTO 1002
      eps = eps/two
      tol1 = one + eps
      IF (tol1 .le. one) GOTO 1002
      eps = eps/two
      tol1 = one + eps
      IF (tol1 .le. one) GOTO 1002
      eps = eps/two
      tol1 = one + eps
      IF (tol1 .le. one) GOTO 1002
      eps = eps/two
      tol1 = one + eps
      GOTO 1003
 1002 CONTINUE
      eps = SQRT(eps)
c
c  initialization
c
      a = ax
      b = bx
      v = a + c*(b - a)
      w = v
      x = v
      e = zero
      fx = f(x)
      fv = fx
      fw = fx
c
c  main loop starts here
c
   20 CONTINUE
      xm = p5*(a + b)
      tol1 = eps*ABS(x) + tol/three
      tol2 = two*tol1
c
c  check STOPping criterion
c
      IF (ABS(x - xm) .le. tol2 - p5*(b - a)) GOTO 90
c
c is golden-section necessary
c
      IF (ABS(e) .le. tol1) GOTO 40
c
c  fit parabola
c
      r = (x - w)*(fx - fv)
      q = (x - v)*(fx - fw)
      p = (x - v)*q - (x - w)*r
      q = two*(q - r)
cwie      IF (q.gt.0.0) p = -p
      IF (q < zero) p = -p
      q = ABS(q)
      r = e
      e = d
c
c  is parabola acceptable
c
cwie   30 IF (ABS(p).ge.ABS(0.5*q*r)) GOTO 40
cwie      IF (p.le.q*(a - x)) GOTO 40
cwie      IF (p.ge.q*(b - x)) GOTO 40
      IF (ABS(p) .le. ABS(p5*q*r)) GOTO 40
      IF (p .ge. q*(a - x)) GOTO 40
      IF (p .le. q*(b - x)) GOTO 40
c
c  a parabolic interpolation step
c
      d = p/q
      u = x + d
c
c  f must not be evaluated too CLOSE to ax or bx
c
      IF (u - a < tol2) d = SIGN(tol1,xm - x)
      IF (b - u < tol2) d = SIGN(tol1,xm - x)
      GOTO 50
c
c  a golden-section step
c
   40 CONTINUE
      IF (x .ge. xm) THEN
         e = a - x
      ELSE
         e = b - x
      END IF
      d = c*e
c
c  f must not be evaluated too CLOSE to x
c
   50 CONTINUE
      IF (ABS(d) .ge. tol1) THEN
         u = x + d
      ELSE
         u = x + SIGN(tol1,d)
      END IF
      fu = f(u)
c
c  update  a, b, v, w, and x
c
cwie      IF (fu.gt.fx) GOTO 60
      IF (fu < fx) GOTO 60
      IF (u .ge. x) THEN
         a = x
      ELSE
         b = x
      END IF
      v = w
      fv = fw
      w = x
      fw = fx
      x = u
      fx = fu
      GOTO 20
   60 CONTINUE
      IF (u < x) THEN
         a = u
      ELSE
         b = u
      END IF
cwie      IF (fu.le.fw) GOTO 70
      IF (fu .ge. fw) GOTO 70
      IF (w .eq. x) GOTO 70
cwie      IF (fu.le.fv) GOTO 80
      IF (fu .ge. fv) GOTO 80
      IF (v .eq. x) GOTO 80
      IF (v .eq. w) GOTO 80
      GOTO 20
   70 CONTINUE
      v = w
      fv = fw
      w = u
      fw = fu
      GOTO 20
   80 CONTINUE
      v = u
      fv = fu
      GOTO 20
c
c  END of main loop
c
   90 CONTINUE
      fmax = x

      END FUNCTION fmax
