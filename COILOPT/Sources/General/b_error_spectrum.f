      SUBROUTINE b_error_spectrum
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE stel_constants
      USE boundary
      USE bnorm_mod
      USE mpi_params                                         !mpi stuff
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: n, mn
      REAL(rprec) :: sinmn, cosmn, dluvdu(nedge)
!-----------------------------------------------
!
!     Want Fourier coefficients of b-error in straight line coordinates
!
!     bmn_err = integral[ sim(mu' - nv') b_error(u',v') du' dv'
!
!     where u' = u + lambda(u),  v' = v. Converting to u,v:
!
!     bmn_err = integral[ sin{m(u+lambda) - nv} b_error(u,v) (1 + d(lambda)/du) du dv
!
!     for m != 0, integrate this by parts to get:
!
!     bmn_err = -1/m integral[ cos{m(u+lambda) - nv} d(b_error)/du du dv ]

!     Compute lamda and d(lambda)/du in u,v space
      DO n = 1, nedge
         luv(n) = zero; dluvdu(n) = 1
         DO mn = 1, mnmax
            sinmn = SIN(xm_b(mn)*thetab(n) - xn_b(mn)*phib(n))
            cosmn = COS(xm_b(mn)*thetab(n) - xn_b(mn)*phib(n))
            luv(n) = luv(n) + lmns_b(mn)*sinmn
            dluvdu(n) = dluvdu(n) + lmns_b(mn)*cosmn*xm_b(mn)
         END DO
      END DO

!     Spectrum of b-error in straight line coordinates
      DO mn = 1, mnmax_bmn
         bmn(mn) = zero
         DO n = 1, nedge
            sinmn = SIN (xm_bmn(mn)*(thetab(n) + luv(n))
     1            -      xn_bmn(mn)*phib(n))
            bmn(mn) = bmn(mn) + b_error(n)*sinmn*dluvdu(n)
         END DO
         bmn(mn) = 2*bmn(mn)/nedge
         IF( myid .eq. master ) THEN                             !mpi stuff
           PRINT 100, NINT(xm_bmn(mn)), NINT(xn_bmn(mn))/3, bmn(mn)
         END IF     !IF( myid .eq. master )                     !mpi stuff
      END DO
  100 FORMAT(2i5,1pe15.5)

      END SUBROUTINE b_error_spectrum
