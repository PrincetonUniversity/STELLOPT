
       SUBROUTINE get_lambda(nsurf, zet, thet, aux1, iotac, fun, 
     1   dfun, dlam)
       USE stel_kinds
       USE normalize_data, ONLY: lasym_v
       USE ballooning_data
       USE general_dimensions
       USE fmesh_quantities
       IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       INTEGER, INTENT(in) :: nsurf
       REAL(rprec), INTENT(in) :: iotac, aux1, thet, zet
       REAL(rprec), INTENT(out) :: fun, dfun, dlam
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER :: j, lj
       REAL(rprec):: lamb, ccosi, ssine, arg
!-----------------------------------------------
       lamb = 0
       dlam = 0
       fourier: DO j = 1, mnmax_v
         arg = xm_v(j)*thet-xn_v(j)*zet
         ssine = sin(arg)
         ccosi = cos(arg)
         lj = mnmax_v*(nsurf-1)+j
         lamb = lamb + lmnsf(lj)*ssine                               ! Fourier invert lambda and d(lambda)/dzeta
         IF (lasym_v) lamb = lamb + lmnsf(lj)*ccosi
         dlam = dlam - xn_v(j)*lmnsf(lj)*ccosi
         IF (lasym_v) dlam = dlam + xn_v(j)*lmncf(lj)*ssine
       END DO fourier
       fun = aux1 + lamb - iotac*zet
       dfun = dlam - iotac

!      END SUBROUTINE fun00
      END
