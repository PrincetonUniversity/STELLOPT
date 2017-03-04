      SUBROUTINE amplitud(rcenter, zcenter, angin, r0c, z0c, rhoc, rhos,
     1   xpts, xin, yin)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE Vname0
      USE Vname1
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec) rcenter, zcenter, r0c, z0c
      REAL(rprec), DIMENSION(ntheta) :: angin, xpts, xin, yin
      REAL(rprec), DIMENSION(0:mrho-1) :: rhoc, rhos
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: mrz
      INTEGER :: m, j
      REAL(rprec) :: xmult, arg, xi, yi, t1, t2,
     1   r1c(mu), r1s(mu), z1c(mu), z1s(mu), tnorm
C-----------------------------------------------
c*****************************************************************
c       This SUBROUTINE assigns initial guesses for angles and
c       Fourier mode amplitudes to the appropriate components of
c       the xvec array
c*****************************************************************
      r0c = rcenter
      z0c = zcenter
      xpts(:ntheta) = angin(:ntheta)

      mrz = mpol-1
      xmult = 2._dp/ntheta

      r1c = zero
      r1s = zero
      z1c = zero
      z1s = zero
      rhoc = zero
      rhos = zero

      DO m = 1, mrz
         DO j = 1, ntheta
            arg = angin(j)
            xi = xmult*(xin(j)-rcenter)
            yi = xmult*(yin(j)-zcenter)
            r1c(m) = r1c(m) + COS(m*arg)*xi
            r1s(m) = r1s(m) + SIN(m*arg)*xi
            z1c(m) = z1c(m) + COS(m*arg)*yi
            z1s(m) = z1s(m) + SIN(m*arg)*yi
         END DO
      END DO

      r10 = SQRT( r1c(1)**2 + r1s(1)**2 + z1c(1)**2 + z1s(1)**2 )
      WRITE(3,'(/,3(a,1pe10.3))')
     1    ' RAXIS = ', rcenter,' ZAXIS = ', zcenter,' R10 = ',r10
      WRITE(*,'(/,3(a,1pe10.3))')
     1    ' RAXIS = ', rcenter,' ZAXIS = ', zcenter,' R10 = ',r10

*
*     INITIALIZE POLAR RHOs for m = 0 thru mrz-1
*
      DO m = 0, mrz-1
         IF( m.le.1 )then
            t1 = t1m(m+1)
            rhoc(m) = 0.5_dp*(r1c(m+1) + z1s(m+1))/t1
            rhos(m) = 0.5_dp*(r1s(m+1) - z1c(m+1))/t1
         ELSE
            t1 = t1m(m+1)
            t2 = t2m(m-1)
            tnorm = 0.5_dp/(t1**2 + t2**2)
            rhoc(m) = tnorm*( (r1c(m+1) + z1s(m+1))*t1
     1                   +  (r1c(m-1) - z1s(m-1))*t2 )
            rhos(m) = tnorm*( (r1s(m+1) - z1c(m+1))*t1
     1                   +  (r1s(m-1) + z1c(m-1))*t2 )
      END IF
      END DO

      END SUBROUTINE amplitud
