      FUNCTION weightfcn(thetaw, phiw, theta0, phi0, wtheta, wphi,
     .   amplw, Nfp)
C EAL 9/00      Form SUM of gaussian weights, in theta, phi plane
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: nfp
      REAL(rprec) :: thetaw, phiw, theta0, phi0,
     1   wtheta, wphi, amplw(3)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec) :: e, pi, twopi, weightfcn
      REAL(rprec) :: c1, c2, aweight
      INTEGER :: i, j
C-----------------------------------------------
      e=0
      Pi=4*ATAN(1._dp)
      twopi = 2*pi
      DO j = 1,3
       c2=phi0 + REAL(j-1,rprec)/twopi/nfp
       DO i = 1,5
        c1=theta0 + (REAL(i-1,rprec)/4)*twopi
        e=e+Amplw(j)*aweight(thetaw,phiw,c1,c2,wtheta,wphi)
       ENDDO
      ENDDO
      weightfcn = e
      RETURN
      END FUNCTION weightfcn
