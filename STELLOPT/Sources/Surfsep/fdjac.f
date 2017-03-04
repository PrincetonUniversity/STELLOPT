      SUBROUTINE fdjac(n, x, fvec, np, df)
!  (c) copr. 1986-92 numerical recipes software

      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER n, np
      REAL(rprec), DIMENSION(n) :: x, fvec
      REAL(rprec), DIMENSION(np,np) :: df
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: nmax = 40
      REAL(rprec), PARAMETER :: eps = 1.e-4_dp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: j
      REAL(rprec) :: h, temp
      REAL(rprec), DIMENSION(nmax) :: f
!-----------------------------------------------
!
      DO j = 1, n
         temp = x(j)
         h = eps*ABS(temp)
         IF (h == 0._dp) h = eps
         x(j) = temp + h
         h = x(j) - temp
         CALL funcv (n, x, f)
         x(j) = temp
         df(:n,j) = (f(:n)-fvec(:n))/h
      END DO

      END SUBROUTINE fdjac
