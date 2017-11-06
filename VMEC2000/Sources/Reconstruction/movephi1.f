      SUBROUTINE movephi1 (gphifacx)   !wiecode
      USE vmec_main
      USE vsvd
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec) :: gphifacx
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec) :: phifac0
C-----------------------------------------------
!
!     UPDATE PHIEDGE SCALE FACTOR BY THE FORMULA
!     FDOT/F = -2*(1 - Ipexp/Ipcode), WHERE F = PHISCALE
!
      IF (ctor .eq. zero) STOP 'ctor = 0 in movephi1'
      phifac0 = phifac*(currv/ctor)
      gphifacx = rsfac*c1pm2*(phifac0 - phifac)

      END SUBROUTINE movephi1
