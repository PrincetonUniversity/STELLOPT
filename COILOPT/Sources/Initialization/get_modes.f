      SUBROUTINE get_modes (mmax, nmax, xm, xn, mnmax)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE stel_constants
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: m, n, mn, mmax, nmax, mnmax
      REAL(rprec) :: xm(*), xn(*)
!-----------------------------------------------

      m = 0
      mn = 0
      DO n = 0, nmax, 3
         mn = mn + 1
         xm(mn) = m
         xn(mn) = n
      END DO

      DO m = 1, mmax
         DO n = -nmax, nmax, 3
            mn = mn + 1
            xm(mn) = m
            xn(mn) = n
         END DO
      END DO

      mnmax = mn

      END SUBROUTINE get_modes
