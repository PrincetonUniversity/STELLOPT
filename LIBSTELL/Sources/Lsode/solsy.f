      SUBROUTINE solsy(wm, iwm, x, tem)
      USE stel_kinds
      USE liprec, ONLY: li_gbsl, li_gesl
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, DIMENSION(*) :: iwm
      REAL(rprec), DIMENSION(*) :: wm, x, tem
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
      INTEGER :: i, meband, ml, mu
      REAL(rprec) :: di, hl0, phl0, r
C-----------------------------------------------
clll. optimize
c-----------------------------------------------------------------------
c this routine manages the solution of the linear system arising from
c a chord iteration.  it is called if miter .ne. 0.
c if miter is 1 or 2, it calls sgesl to accomplish this.
c if miter = 3 it updates the coefficient h*el0 in the diagonal
c matrix, and then computes the solution.
c if miter is 4 or 5, it calls sgbsl.
c communication with solsy uses the following variables..
c wm    = real work space containing the inverse diagonal matrix if
c         miter = 3 and the lu decomposition of the matrix otherwise.
c         storage of matrix elements starts at wm(3).
c         wm also contains the following matrix-related data..
c         wm(1) = sqrt(uround) (not used here),
c         wm(2) = hl0, the previous value of h*el0, used if miter = 3.
c iwm   = integer work space containing pivot information, starting at
c         iwm(21), if miter is 1, 2, 4, or 5.  iwm also contains band
c         parameters ml = iwm(1) and mu = iwm(2) if miter is 4 or 5.
c x     = the right-hand side vector on input, and the solution vector
c         on output, of length n.
c tem   = vector of work space of length n, not used in this version.
c iersl = output flag (in common).  iersl = 0 if no trouble occurred.
c         iersl = 1 if a singular matrix arose with miter = 3.
c this routine also uses the common variables el0, h, miter, and n.
c-----------------------------------------------------------------------
      iersl = 0
      SELECT CASE (miter)
      CASE DEFAULT
         CALL li_gesl (wm(3:3+n*n), n, n, iwm(21:21+n), x, 0)
         RETURN
c
      CASE (3)
         phl0 = wm(2)
         hl0 = h*el0
         wm(2) = hl0
         IF (hl0 /= phl0) THEN
            r = hl0/phl0
            DO i = 1, n
               di = 1.0_dp - r*(1.0_dp - 1.0_dp/wm(i+2))
               IF (ABS(di) == 0.0_dp) go to 390
               wm(i+2) = 1.0_dp/di
            END DO
         END IF
         x(:n) = wm(3:n+2)*x(:n)
         RETURN
  390    CONTINUE
         iersl = 1
         RETURN
c
      CASE (4:5)
         ml = iwm(1)
         mu = iwm(2)
         meband = 2*ml + mu + 1
         CALL li_gbsl (wm(3:3+meband*n), meband, n, ml, mu, 
     1                 iwm(21:21+n), x, 0)
         RETURN
      END SELECT

      END SUBROUTINE solsy
