      SUBROUTINE bextrema(modb, bmin, bmax, nzeta, ntheta)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: nzeta, ntheta
      REAL(rprec), INTENT(in)  :: modb(nzeta,ntheta)
      REAL(rprec), INTENT(out) :: bmin(ntheta), bmax(ntheta)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ku
C-----------------------------------------------
!
!     Computes MAX, MIN of |B| along v (zeta) between two angle lines (theta = 0,pi)
!
      DO ku = 1,ntheta
         bmin(ku)  = MINVAL(modb(:,ku))
         bmax(ku)  = MAXVAL(modb(:,ku))
      ENDDO

      END SUBROUTINE bextrema
