      SUBROUTINE getangle(rval, zval, angle, rcenter, zcenter)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE Vname0
      USE Vname1
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ntheta, nphi) :: rval, zval, angle
      REAL(rprec), DIMENSION(ntheta) :: rcenter, zcenter
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, j, iterate
      REAL(rprec), DIMENSION(nv) ::
     1  rcos, rsin, zcos, zsin, phiangle
      REAL(rprec) :: xc, yc, dnum, denom, delangle
C-----------------------------------------------
c**********************************************************************
c       Compute angle offset consistent with constraint Z1n = Z1,-n
c       Note: This is done iteratively, since elongation is unknown
c**********************************************************************
      DO i = 1, nphi
         DO j = 1, ntheta
            angle(j,i) = twopi*(j - 1)/REAL(ntheta,rprec)
         END DO
      END DO
      DO iterate = 1, 5
         DO i = 1, nphi
            rcos(i) = zero
            rsin(i) = zero
            zcos(i) = zero
            zsin(i) = zero
            DO j = 1, ntheta
               xc = rval(j,i) - rcenter(i)
               yc = zval(j,i) - zcenter(i)
               rcos(i) = rcos(i) + COS(angle(j,i))*xc
               rsin(i) = rsin(i) + SIN(angle(j,i))*xc
               zcos(i) = zcos(i) + COS(angle(j,i))*yc
               zsin(i) = zsin(i) + SIN(angle(j,i))*yc
            END DO
         END DO
c**********************************************************************
c       Compute new angles starting from offset phiangle(i)
c**********************************************************************
         dnum = zero
         denom = zero
         dnum = SUM(zsin(:nphi))
         denom = SUM(rcos(:nphi))
         IF (denom .ne. zero) THEN
            elongate = dnum/denom
         ELSE
            elongate = 1.E10
         END IF

         delangle = zero
         DO i = 1, nphi
            phiangle(i) = ATAN2(elongate*zcos(i)-rsin(i),
     1                          elongate*zsin(i)+rcos(i))
            delangle = MAX(delangle,ABS(phiangle(i)))
            angle(:,i) = angle(:,i) + phiangle(i)
         END DO
         IF (delangle < 0.02_dp) EXIT
      END DO
      WRITE (*, 1010) elongate, raxis(1), zaxis(1)
      WRITE (3, 1010) elongate, raxis(1), zaxis(1)
      WRITE (*, 1020) ntheta, nphi
      WRITE (3, 1020) ntheta, nphi
      WRITE (*, 1030) mpol1, ntor
      WRITE (3, 1030) mpol1, ntor
 1010 FORMAT(' Average elongation = ',1pe11.4,/,' Raxis = ',1pe12.4,
     1   ' Zaxis = ',1pe12.4)
 1020 FORMAT(' Number of Theta Points Matched = ',i5,/,
     1   ' Number of Phi Planes = ',i5)
 1030 FORMAT(' Max Poloidal Mode Number = ',i5,/,
     1   ' Max Toroidal Mode Number = ',i5)

      END SUBROUTINE getangle
