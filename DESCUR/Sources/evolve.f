      SUBROUTINE evolve(g11)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE Vname0
      USE Vname1
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec) g11
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: bmax = 0.15_dp
      REAL(rprec) :: ftest, dtau, otav, b1, fac
C-----------------------------------------------

      ftest = gnorm/g11
      dtau = ABS(one - ftest)
      g11 = gnorm
      otav = dtau/delt
      dtau = delt*otav + 1.e-3_dp
      dtau = MIN(bmax,dtau)
      b1 = one - 0.5_dp*dtau
      fac = one/(one + 0.5_dp*dtau)
      xdot(:n2) = fac*(xdot(:n2)*b1-delt*gvec(:n2))
      xvec(:n2) = xvec(:n2) + xdot(:n2)*delt

      END SUBROUTINE evolve
