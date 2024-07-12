      SUBROUTINE obtain_field_line(zetacn, thetacn, iotac,
     1    nsurf, alpha)
      USE stel_kinds
      USE normalize_data, ONLY: lasym_v                                   ! 110909 RS: logical for asymmetric input (if True)
      USE ballooning_data
      USE general_dimensions
      USE fmesh_quantities
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: nsurf
      REAL(rprec), INTENT(IN) :: zetacn, thetacn, iotac
      REAL(rprec), INTENT(OUT) :: alpha
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: j, lj
      REAL(rprec):: lambda0, ssine, arg, ccosi
!-----------------------------------------------

      lambda0 = 0
      fourier: DO j = 1, mnmax_v                                         ! Fourier invert lambda at initial point (thetacn, zetacn)
         arg = xm_v(j)*thetacn-xn_v(j)*zetacn
         ssine = SIN(arg); ccosi = COS(arg)
         lj = mnmax_v*(nsurf-1)+j
         lambda0 = lambda0 + lmnsf(lj)*ssine
         IF (lasym_v) lambda0 = lambda0 + lmncf(lj)*ccosi                ! 110909  RS: Asymmetric input
      ENDDO fourier

      alpha = thetacn + lambda0 - iotac*zetacn                           ! obtain field line label value

      END SUBROUTINE obtain_field_line
