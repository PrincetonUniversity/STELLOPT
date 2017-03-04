      SUBROUTINE newt(x, n, check)
!  (c) copr. 1986-92 numerical recipes software

!     USEs fdjac,fmin_nr,lnsrch,lubksb,ludcmp

      USE stel_kinds
      USE newtv
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: n
      LOGICAL :: check
      REAL(rprec), DIMENSION(n) :: x
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: maxits = 250
      REAL(rprec), PARAMETER :: tolf = 1.e-4_dp
      REAL(rprec), PARAMETER :: tolMIN = 1.e-6_dp
      REAL(rprec), PARAMETER :: tolx = 1.e-7_dp
      REAL(rprec), PARAMETER :: stpmx = 100
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, its
      INTEGER, DIMENSION(np) :: indx
      REAL(rprec) :: d, den, f, fold, stpmax, test, sumn
      REAL(rprec), DIMENSION(np,np) :: fjac
      REAL(rprec), DIMENSION(np) :: g, p, xold
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      REAL(rprec) , EXTERNAL :: fmin_nr
!-----------------------------------------------
!
      IF (n > np) STOP 'N > NP in NEWT'
      nn = n
      f = fmin_nr(x)
      test = 0
      DO i = 1, n
         IF (ABS(fvec(i)) > test) test = ABS(fvec(i))
      END DO
      IF (test < tolf/100) RETURN
!
      sumn = DOT_PRODUCT(x(:n),x(:n))
!
      stpmax = stpmx*MAX(SQRT(sumn),REAL(n,rprec))
!
      DO its = 1, maxits
         CALL fdjac (n, x, fvec, np, fjac)
         DO i = 1, n
            sumn = SUM(fjac(:n,i)*fvec(:n))
            g(i) = sumn
         END DO
         fold = f
         xold(:n) = x(:n)
         p(:n) = -fvec(:n)
         CALL ludcmp (fjac, n, np, indx, d)
         CALL lubksb (fjac, n, np, indx, p)
         CALL lnsrch (n, xold, fold, g, p, x, f, stpmax, check, fmin_nr)
         test = MAXVAL(ABS(fvec(:n)))
         IF (test < tolf) THEN
            check = .FALSE.
            RETURN
         ENDIF
         IF (check) THEN
            den = MAX(f,.5_dp*n)
            test = MAXVAL(ABS(g(:n))*MAX(ABS(x(:n)),1._dp)/den)
            IF (test < tolmin) THEN
               check = .TRUE.
            ELSE
               check = .FALSE.
            ENDIF
            RETURN
         ENDIF
         test = MAXVAL( ABS(x(:n)-xold(:n))/MAX(ABS(x(:n)),1._dp) )
         IF (test < tolx) RETURN
      END DO
!      WRITE(*,*) 'maxits exceeded in newt'
      check = .TRUE.

      END SUBROUTINE newt
