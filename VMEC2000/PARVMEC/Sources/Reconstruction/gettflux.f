      SUBROUTINE gettflux
      USE vsvd
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec) :: p01=1.e-2_dp, zero = 0.0_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec) :: tfac, phifac0
C-----------------------------------------------

!
!       UPDATE PHIEDGE SCALE FACTOR BY THE FORMULA
!       FDOT/F = -2*(1 - apres/aminor), WHERE F = PHISCALE
!
      tfac = p01*rsfac

      IF (imatch_phiedge .eq. 0) THEN
         IF (aminor .eq. zero) THEN
            phifac0 = phifac
         ELSE
           phifac0 = phifac*(apres/aminor)
         END IF
         gphifac = tfac*(phifac0 - phifac)

      ELSE IF (imatch_phiedge .eq. 2) THEN
         CALL getlim
         gphifac = tfac*phifac*gphifac
      ENDIF

      END SUBROUTINE gettflux
