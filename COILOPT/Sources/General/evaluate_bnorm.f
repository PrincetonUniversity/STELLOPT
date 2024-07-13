      SUBROUTINE evaluate_bnorm
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE stel_constants
      USE boundary
      USE bnorm_mod
      USE tf_coils
      USE safe_open_mod
      USE mpi_params                                         !mpi stuff
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, ku, kv, n, nmatch, iunit=29, ierr
      REAL (rprec) :: pscale, theta, phi, bnrml
!-----------------------------------------------
!     WRITE field normal at plasma edge vs u and v for plot
      IF (myid .EQ. master) CALL safe_open(iunit, ierr, 
     1                      'b_norm_eq.dat', 'unknown', 'formatted')
      pscale = twopi/nfp

      n = 0
      DO ku = 1, nu
          DO kv = 1, nv
              n = n + 1
              phi = phib(n)
              theta = thetab(n)
              bnrml = zero
              DO i=1,mnbn_max
                  bnrml = bnrml + bn_coef(i)*
     1                    SIN(xbn_m(i)*theta + xbn_n(i)*phi)
              END DO
!     bnrml (Tesla) = bnrml*i_pol*mu, where i_pol = Ipol/nfp is the
!     total poloidal current per period obtained from R*Bt = mu*Ipol/twopi
!     (R*BT = R-BTOR(s=1) is output from vmec run producing bnrml)
!     For m3.b15, i_pol = 3.0159*mu (where R*Bt = 1.44 T-m)
              bnormal_match(n) = bnrml*i_pol*mu0
              IF( myid .EQ. master ) WRITE(iunit,1000) phi/pscale,
     1            theta/twopi, bnrml*i_pol*mu0
 1000         FORMAT(1p,3e16.7)
          END DO
          IF( myid .EQ. master ) WRITE(iunit,1010)
 1010     FORMAT("")
      END DO
      CLOSE(iunit)
      nmatch = n

      END SUBROUTINE evaluate_bnorm
