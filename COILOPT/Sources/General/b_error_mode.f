      SUBROUTINE b_error_mode (m)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE stel_constants
      USE boundary
      USE bnorm_mod
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: ku, kv, n, m, mn
      REAL(rprec) :: sinmn
!-----------------------------------------------

      OPEN (unit=33,file='bmnerr.dat',status='unknown')

      mn = m

      bmn_error(:) = zero
      n = 0
      DO ku = 1, nu
         DO kv = 1, nv
            n = n + 1
            sinmn = SIN (xm_bmn(mn)*(thetab(n) + luv(n))
     1            - xn_bmn(mn)*phib(n))
            bmn_error(n) = bmn(mn)*sinmn
            WRITE (33,100) nfp*phib(n)/twopi, thetab(n)/twopi,
     1         bmn_error(n)
         END DO
         WRITE (33,110)
      END DO

  100 FORMAT (1p,3e15.6)
  110 FORMAT ("")
      CLOSE (33)

      END SUBROUTINE b_error_mode
