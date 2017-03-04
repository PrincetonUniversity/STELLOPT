      SUBROUTINE intdy(t, k, yh, nyh, dky, iflag)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER k, nyh, iflag
      REAL(rprec) t
      REAL(rprec), DIMENSION(nyh,*) :: yh
      REAL(rprec), DIMENSION(*) :: dky
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /ls0001/
      COMMON /ls0001/ rowns(209), ccmax, el0, h, hmin, hmxi, hu, rc, tn
     1   , uround, iownd(14), iowns(6), icf, ierpj, iersl, jcur, jstart
     2   , kflag, l, meth, miter, maxord, maxcor, msbp, mxncf, n, nq,
     3   nst, nfe, nje, nqu
      INTEGER   iownd, iowns, icf, ierpj, iersl, jcur, jstart, kflag, l
     1   , meth, miter, maxord, maxcor, msbp, mxncf, n, nq, nst, nfe,
     2   nje, nqu
      REAL(rprec) ::rowns, ccmax, el0, h, hmin, hmxi, hu, rc, tn,
     1   uround
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0
      INTEGER :: ic, j, jb, jb2, jj, jj1, jp1
      REAL(rprec) :: c, r, s, tp
C-----------------------------------------------
clll. optimize
c-----------------------------------------------------------------------
c intdy computes interpolated values of the k-th derivative of the
c dependent variable vector y, and stores it in dky.  this routine
c is called within the package with k = 0 and t = tout, but may
c also be called by the user for any k up to the current order.
c (see detailed instructions in the usage documentation.)
c-----------------------------------------------------------------------
c the computed values in dky are gotten by interpolation using the
c nordsieck history array yh.  this array corresponds uniquely to a
c vector-valued polynomial of degree nqcur or less, and dky is set
c to the k-th derivative of this polynomial at t.
c the formula for dky is..
c              q
c  dky(i)  =  SUM  c(j,k) * (t - tn)**(j-k) * h**(-j) * yh(i,j+1)
c             j=k
c where  c(j,k) = j*(j-1)*...*(j-k+1), q = nqcur, tn = tcur, h = hcur.
c the quantities  nq = nqcur, l = nq+1, n = neq, tn, and h are
c communicated by common.  the above sum is done in reverse order.
c iflag is returned negative if either k or t is out of bounds.
c-----------------------------------------------------------------------
      iflag = 0
      IF (k>=0 .and. k<=nq) THEN
         tp = tn - hu - 100*uround*(tn + hu)
         IF ((t - tp)*(t - tn) > zero) GOTO 90
c
         s = (t - tn)/h
         ic = 1
         IF (k /= 0) THEN
            jj1 = l - k
            DO jj = jj1, nq
               ic = ic*jj
            END DO
         END IF
         c = ic
         dky(:n) = c*yh(:n,l)
         IF (k .ne. nq) THEN
            jb2 = nq - k
            DO jb = 1, jb2
               j = nq - jb
               jp1 = j + 1
               ic = 1
               IF (k /= 0) THEN
                  jj1 = jp1 - k
                  DO jj = jj1, j
                     ic = ic*jj
                  END DO
               END IF
               c = ic
               dky(:n) = c*yh(:n,jp1) + s*dky(:n)
            END DO
            IF (k == 0) RETURN
         END IF
         r = h**(-k)
         dky(:n) = r*dky(:n)
         RETURN
c
      END IF
      CALL xerrwv ('intdy--  k (=i1) illegal      ', 30, 51, 0, 1, k, 0
     1   , 0, zero, zero)
      iflag = -1
      RETURN
   90 CONTINUE
      CALL xerrwv ('intdy--  t (=r1) illegal      ', 30, 52, 0, 0, 0, 0
     1   , 1, t, zero)
      CALL xerrwv (
     1'      t not in interval tcur - hu (= r1) to tcur (=r2)      ',
     2   60, 52, 0, 0, 0, 0, 2, tp, tn)
      iflag = -2

      END SUBROUTINE intdy
