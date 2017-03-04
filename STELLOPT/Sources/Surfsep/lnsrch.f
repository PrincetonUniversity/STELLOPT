      SUBROUTINE lnsrch(n, xold, fold, g, p, x, f, stpmax, check, func)
!  (c) copr. 1986-92 numerical recipes software

      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER n
      REAL(rprec) fold, f, stpmax, func
      LOGICAL check
      REAL(rprec), DIMENSION(n) :: xold, g, p, x
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: alf = 1.e-4_dp
      REAL(rprec), PARAMETER :: tolx = 1.e-7_dp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec) :: a, alam, alam2, alamin, b, disc, f2, fold2,
     1    rhs1, rhs2, slope, sumn, test, tmplam
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      EXTERNAL func
!-----------------------------------------------
!
!     USEs func
!
      check = .FALSE.
!
      sumn = DOT_PRODUCT(p,p)
      sumn = SQRT(sumn)
      IF (sumn > stpmax)  p = p*stpmax/sumn
      slope = DOT_PRODUCT(g,p)
      test = MAXVAL(ABS(p)/MAX(ABS(xold),1._dp))
      alamin = tolx/test
      alam = 1
    1 CONTINUE
      x = xold + alam*p
      f = func(x)
      IF (alam < alamin) THEN
         x = xold
         check = .TRUE.
         RETURN
      ELSE IF (f <= fold + alf*alam*slope) THEN
         RETURN
      ELSE
         IF (alam == 1._dp) THEN
            tmplam = -slope/(2*(f - fold - slope))
         ELSE
            rhs1 = f - fold - alam*slope
            rhs2 = f2 - fold2 - alam2*slope
            a = (rhs1/alam**2 - rhs2/alam2**2)/(alam - alam2)
            b = ((-alam2*rhs1/alam**2) + alam*rhs2/alam2**2)
     1           /(alam - alam2)
            IF (a == 0._dp) THEN
               tmplam = -slope/(2*b)
            ELSE
               disc = b*b - 3*a*slope
               tmplam = ((-b) + SQRT(disc))/(3*a)
            ENDIF
            tmplam = MIN(.5_dp*alam,tmplam)
         ENDIF
      ENDIF
      alam2 = alam
      f2 = f
      fold2 = fold
      alam = MAX(tmplam,.1_dp*alam)
      GOTO 1

      END SUBROUTINE lnsrch
