      SUBROUTINE inguess(rin, zin, angle)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE Vname0
      USE Vname1
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ntheta,nphi) :: rin, zin, angle
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, inside, jskip, j1, j2
C-----------------------------------------------
c**********************************************************************
c       this subroutine obtains initial guesses for the centroid at each
c       toroidal plane. by default, the polar axis is set equal to this
c       centroid.  it is imperative that the polar axis be enclosed by
c       the surface (otherwise the computed theta angles will not span
c       [0,2pi]). for certain complex cross-sections, this simple estimate
c       of the polar axis may fail, and the user must supply values for
c       raxis(i),zaxis(i).  in addition, for non-starlike domains, the
c       points along the curve must be monitonically increasing as a
c       function of the arclength from any fixed point on the curve. this
c       ordering is attempted by the subroutine order, but may fail in
c       general if too few theta points are used.
c
c modified 8/13/97, j. breslau. the subroutine now tries several guesses
c  for the magnetic axis position if the first fails. the mean position
c  of pairs of points is taken until one is found to be inside.
c**********************************************************************
c       COMPUTE CENTROID AND
c       ORDER POINTS ON FLUX SURFACE AT EACH TOROIDAL ANGLE
c**********************************************************************
      PRINT *, 'ORDERING SURFACE POINTS'
      DO i = 1, nphi
         r0n(i) = r0n(i) + SUM(rin(:,i))/REAL(ntheta,rprec)
         z0n(i) = z0n(i) + SUM(zin(:,i))/REAL(ntheta,rprec)
         raxis(i) = r0n(i)
         zaxis(i) = z0n(i)
         CALL order (rin(1,i), zin(1,i), raxis(i), zaxis(i), inside)
         IF (inside .eq. 0) THEN
            jskip = ntheta/2
   40       CONTINUE
            jskip = jskip - 1
            IF (jskip .le. 0) THEN
               PRINT *, 'Could not find internal point'
               WRITE (3, *) 'Could not find internal point'
               STOP
            ENDIF
            j1 = 1
            j2 = j1 + jskip
   50       CONTINUE
            raxis(i) = 0.5_dp*(rin(j1,i)+rin(j2,i))
            zaxis(i) = 0.5_dp*(zin(j1,i)+zin(j2,i))
            CALL order (rin(1,i), zin(1,i), raxis(i), zaxis(i), inside)
            IF (inside .eq. 1) THEN
               r0n(i) = raxis(i)
               z0n(i) = zaxis(i)
            ELSE
               j1 = j1 + 1
               IF (j1 .gt. ntheta) GOTO 40
               j2 = j2 + 1
               IF (j2 .gt. ntheta) j2 = 1
               GOTO 50
            ENDIF
         ENDIF
      END DO
c**********************************************************************
c       COMPUTE OPTIMUM ANGLE OFFSETS FOR M = 1 MODES
c**********************************************************************
      CALL getangle (rin, zin, angle, r0n, z0n)

      END SUBROUTINE inguess
