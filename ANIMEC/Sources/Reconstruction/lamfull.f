      SUBROUTINE lamfull(luef, luof, lue, luo)
      USE vmec_main
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns,*), INTENT(out) :: luef, luof
      REAL(rprec), DIMENSION(ns,nzeta,*), INTENT(in) :: lue, luo
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: p5=0.5_dp, c1p5=1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: lt, lt1
C-----------------------------------------------

C
C     COMPUTES LAMBDA ON FULL RADIAL MESH AT THETA=0, PI
C
      DO lt = 1, 2
         lt1 = 1
         IF (lt .eq. 2) lt1 = ntheta2
         luof(1,lt1) = p5*(c1p5*luo(2,1,1)-p5*luo(3,1,1)+c1p5*luo(2,1,
     1      ntheta2)-p5*luo(3,1,ntheta2))
         luef(1,lt1) = p5*(c1p5*lue(2,1,1)-p5*lue(3,1,1)+c1p5*lue(2,1,
     1      ntheta2)-p5*lue(3,1,ntheta2))
         luof(ns,lt1) = c1p5*luo(ns,1,lt1) - p5*luo(ns1,1,lt1)
         luef(ns,lt1) = c1p5*lue(ns,1,lt1) - p5*lue(ns1,1,lt1)
         luof(2:ns1,lt1) = p5*(luo(3:ns1+1,1,lt1)+luo(2:ns1,1,lt1))
         luef(2:ns1,lt1) = p5*(lue(3:ns1+1,1,lt1)+lue(2:ns1,1,lt1))
      END DO

      END SUBROUTINE lamfull
