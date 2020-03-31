      SUBROUTINE prepj(neq, y, yh, nyh, ewt, ftem, savf, wm, iwm, 
     1   f, jac)
      USE stel_kinds
      USE liprec, ONLY: li_gbfa, li_gefa
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER nyh
      INTEGER, DIMENSION(*) :: neq, iwm
      REAL(rprec), DIMENSION(*) :: y
      REAL(rprec), DIMENSION(nyh,*) :: yh
      REAL(rprec), DIMENSION(*) :: ewt, ftem, savf, wm
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
      INTEGER :: i, i1, i2, ier, ii, j, j1, jj, lenp, mba, mband, meb1,
     1   meband, ml, ml3, mu, np1
      REAL(rprec) :: con, di, fac, hl0, r, r0, srur, yj, yjj
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      REAL(rprec) , EXTERNAL :: vnorm
      EXTERNAL f, jac 
C-----------------------------------------------
clll. optimize
c-----------------------------------------------------------------------
c prepj is called by stode to compute and process the matrix
c p = i - h*el(1)*j , where j is an approximation to the jacobian.
c here j is computed by the user-supplied routine jac if
c miter = 1 or 4, or by finite differencing if miter = 2, 3, or 5.
c if miter = 3, a diagonal approximation to j is used.
c j is stored in wm and replaced by p.  if miter .ne. 3, p is then
c subjected to lu decomposition in preparation for later solution
c of linear systems with p as coefficient matrix. this is done
c by sgefa if miter = 1 or 2, and by sgbfa if miter = 4 or 5.
c
c in addition to variables described previously, communication
c with prepj uses the following..
c y     = array containing predicted values on entry.
c ftem  = work array of length n (acor in stode).
c savf  = array containing f evaluated at predicted y.
c wm    = real work space for matrices.  on output it contains the
c         inverse diagonal matrix if miter = 3 and the lu decomposition
c         of p if miter is 1, 2 , 4, or 5.
c         storage of matrix elements starts at wm(3).
c         wm also contains the following matrix-related data..
c         wm(1) = sqrt(uround), used in numerical jacobian increments.
c         wm(2) = h*el0, saved for later use if miter = 3.
c iwm   = integer work space containing pivot information, starting at
c         iwm(21), if miter is 1, 2, 4, or 5.  iwm also contains band
c         parameters ml = iwm(1) and mu = iwm(2) if miter is 4 or 5.
c el0   = el(1) (input).
c ierpj = output error flag,  = 0 if no trouble, .gt. 0 if
c         p matrix found to be singular.
c jcur  = output flag = 1 to indicate that the jacobian matrix
c         (or approximation) is now current.
c this routine also uses the common variables el0, h, tn, uround,
c miter, n, nfe, and nje.
c-----------------------------------------------------------------------
      nje = nje + 1
      ierpj = 0
      jcur = 1
      hl0 = h*el0
      GOTO (100,200,300,400,500) miter
c IF miter = 1, CALL jac and multiply by scalar. -----------------------
  100 CONTINUE
      lenp = n*n
      wm(3:lenp+2) = 0.0_dp
      CALL jac (neq, tn, y, 0, 0, wm(3), n)
      con = -hl0
      wm(3:lenp+2) = wm(3:lenp+2)*con
      GOTO 240
c IF miter = 2, make n calls to f to approximate j. --------------------
  200 CONTINUE
      fac = vnorm(n,savf,ewt)
      r0 = 1000.0_dp*ABS(h)*uround*n*fac
      IF (r0 == 0.0_dp) r0 = 1.0_dp
      srur = wm(1)
      j1 = 2
      DO j = 1, n
         yj = y(j)
         r = MAX(srur*ABS(yj),r0/ewt(j))
         y(j) = y(j) + r
         fac = -hl0/r
         CALL f (neq, tn, y, ftem)
         wm(1+j1:n+j1) = (ftem(:n)-savf(:n))*fac
         y(j) = yj
         j1 = j1 + n
      END DO
      nfe = nfe + n
c add identity matrix. -------------------------------------------------
  240 CONTINUE
      j = 3
      np1 = n + 1
      wm(j:(n-1)*np1+j:np1) = wm(j:(n-1)*np1+j:np1) + 1.0_dp
c DO lu decomposition on p. --------------------------------------------
      CALL li_gefa (wm(3:3+n*n), n, n, iwm(21:21+n), ier)
      IF (ier /= 0) ierpj = 1
      RETURN
c IF miter = 3, construct a diagonal approximation to j and p. ---------
  300 CONTINUE
      wm(2) = hl0
      r = el0*0.1_dp
      y(:n) = y(:n) + r*(h*savf(:n)-yh(:n,2))
      CALL f (neq, tn, y, wm(3))
      nfe = nfe + 1
      DO i = 1, n
         r0 = h*savf(i) - yh(i,2)
         di = 0.1_dp*r0 - h*(wm(i+2)-savf(i))
         wm(i+2) = 1.0_dp
         IF (ABS(r0) >= uround/ewt(i)) THEN
            IF (ABS(di) == 0.0_dp) GOTO 330
            wm(i+2) = 0.1_dp*r0/di
         END IF
      END DO
      RETURN
  330 CONTINUE
      ierpj = 1
      RETURN
c IF miter = 4, CALL jac and multiply by scalar. -----------------------
  400 CONTINUE
      ml = iwm(1)
      mu = iwm(2)
      ml3 = ml + 3
      mband = ml + mu + 1
      meband = mband + ml
      lenp = meband*n
      wm(3:lenp+2) = 0.0_dp
      CALL jac (neq, tn, y, ml, mu, wm(ml3), meband)
      con = -hl0
      wm(3:lenp+2) = wm(3:lenp+2)*con
      GOTO 570
c IF miter = 5, make mband calls to f to approximate j. ----------------
  500 CONTINUE
      ml = iwm(1)
      mu = iwm(2)
      mband = ml + mu + 1
      mba = MIN(mband,n)
      meband = mband + ml
      meb1 = meband - 1
      srur = wm(1)
      fac = vnorm(n,savf,ewt)
      r0 = 1000.0_dp*ABS(h)*uround*n*fac
      IF (r0 == 0.0_dp) r0 = 1.0_dp
      DO j = 1, mba
         y(j:n:mband) = y(j:n:mband) + MAX(srur*ABS(y(j:n:mband)),r0/
     1      ewt(j:n:mband))
         CALL f (neq, tn, y, ftem)
         DO jj = j, n, mband
            y(jj) = yh(jj,1)
            yjj = y(jj)
            r = MAX(srur*ABS(yjj),r0/ewt(jj))
            fac = -hl0/r
            i1 = MAX(jj - mu,1)
            i2 = MIN(jj + ml,n)
            ii = jj*meb1 - ml + 2
            wm(ii+i1:i2+ii) = (ftem(i1:i2)-savf(i1:i2))*fac
         END DO
      END DO
      nfe = nfe + mba
c add identity matrix. -------------------------------------------------
  570 CONTINUE
      ii = mband + 2
      wm(ii:(n-1)*meband+ii:meband) = wm(ii:(n-1)*meband+ii:meband) +
     1   1.0_dp
c DO lu decomposition of p. --------------------------------------------
      CALL li_gbfa (wm(3:3+meband*n), meband, n, ml, mu, 
     1              iwm(21:21+n), ier)
      IF (ier /= 0) ierpj = 1

      END SUBROUTINE prepj
