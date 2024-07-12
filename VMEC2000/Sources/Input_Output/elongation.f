      SUBROUTINE elongation (r1, z1, waist, height)
      USE vmec_main
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), INTENT(out) :: waist(2), height(2)
      REAL(rprec), INTENT(in), DIMENSION(ns,nzeta,ntheta3,0:1) :: r1, z1
      INTEGER :: nv, n1
C-----------------------------------------------
!
!     Compute Waist thickness, Height in phi = 0, pi symmetry planes
!
      n1 = 0
      DO nv = 1, nzeta/2+1
         IF (nv.ne.1 .and. nv.ne.nzeta/2+1) CYCLE
         n1 = n1+1
         waist(n1) = (r1(ns,nv,1,0)       + r1(ns,nv,1,1)) -
     1               (r1(ns,nv,ntheta2,0) + r1(ns,nv,ntheta2,1))
         height(n1) = 2*MAXVAL(ABS(z1(ns,nv,:,0) + z1(ns,nv,:,1)))
      END DO

      END SUBROUTINE elongation
