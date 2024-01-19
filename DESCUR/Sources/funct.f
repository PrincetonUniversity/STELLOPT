      SUBROUTINE funct(r0c, z0c, rhoc, rhos, xpts, gr0c, gz0c, grhoc,
     1                 grhos, gpts, fsq, xin, yin, mrho_in)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE Vname0
      USE Vname1
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER mrho_in
      REAL(rprec) fsq, r0c, z0c, gr0c, gz0c
      REAL(rprec), DIMENSION(ntheta) :: xpts, gpts, xin, yin
      REAL(rprec), DIMENSION(0:mrho_in-1) ::
     1  rhoc, rhos, grhoc, grhos
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: m, mrho1
      REAL(rprec),DIMENSION(ntheta) ::
     1   gcon, gtt, r1, z1, rt1, zt1
      REAL(rprec), DIMENSION(ntheta,0:mrho_in) :: cosa, sina
      REAL(rprec) :: denom, rmc_p, zms_p, rms_p, zmc_p,
     1   t1, t2, tnorm
C-----------------------------------------------
!
!
!     FORCES: dW/dRmn = SUM(i)[ R(i) - RIN(i)]*COS(m*u[i]) ....
!      dW/dZmn = SUM(i)[ Z(i) - ZIN(i)]*SIN(m*u[i]) ....
!      dW/du[i] =    rt(i)*[R(i)-RIN(i)] + zt(i)*[Z(i) - ZIN(i)]
!     THE NORM ON THE ANGLE FORCE (dW/du[i]) FOLLOWS FROM NEWTONS
!     LAW AND IS APPROXIMATELY GTT = RT**2 + ZT**2
!
      denom = zero
      specw = zero
      
      mrho1 = mrho_in-1
      r1(:ntheta) = -xin(:ntheta)
      z1(:ntheta) = -yin(:ntheta)
      rt1(:ntheta) = zero
      zt1(:ntheta) = zero
      cosa(:,0) = one
      cosa(:,1) = COS(xpts)
      sina(:,0) = zero
      sina(:,1) = SIN(xpts)

!
!     COMPUTE CURVE R1 = R - Rin, Z1 = Z - Zin FOR INPUT POINTS
!     NOTE DIMENSIONS: Rm(0:MRHO), Zm(0:MRHO), but RHO(0:MRHO-1)
!
      DO m = 2, mrho_in
         cosa(:,m) = cosa(:,m-1)*cosa(:,1) - sina(:,m-1)*sina(:,1)
         sina(:,m) = sina(:,m-1)*cosa(:,1) + cosa(:,m-1)*sina(:,1)
      END DO

      DO m = 0, mrho_in
         CALL getrz(rmc_p,rms_p,zmc_p,zms_p,r0c,z0c,rhoc,rhos,
     1      m,mrho_in)
!
!        COMPUTE SPECTRAL WIDTH OF CURVE (DIAGNOSTIC)
!
         t2 = rmc_p*rmc_p + zmc_p*zmc_p + rms_p*rms_p + zms_p*zms_p
         denom = denom + t2*xmpq(m,2)
         specw = specw + xmpq(m,1)*t2

         gtt(:ntheta) = rmc_p*cosa(:ntheta,m) + rms_p*sina(:ntheta,m)
         gcon(:ntheta) = zmc_p*cosa(:ntheta,m) + zms_p*sina(:ntheta,m)
         r1(:ntheta) = r1(:ntheta) + gtt(:ntheta)
         z1(:ntheta) = z1(:ntheta) + gcon(:ntheta)
         rt1(:ntheta) = rt1(:ntheta) + dm1(m)*(rms_p*cosa(:ntheta,m)-
     1      rmc_p*sina(:ntheta,m))
         zt1(:ntheta) = zt1(:ntheta) + dm1(m)*(zms_p*cosa(:ntheta,m)-
     1      zmc_p*sina(:ntheta,m))
      END DO

      specw = specw/denom

      gtt(:ntheta) = rt1(:ntheta)**2 + zt1(:ntheta)**2
      gpts(:ntheta) = r1(:ntheta)*rt1(:ntheta)+z1(:ntheta)*zt1(:ntheta)

!
!     COMPUTE MEAN-SQUARE DEVIATION BETWEEN POINTS AND FIT
!
!     t1 = MAXVAL(gtt(:ntheta))
!     gpts(:ntheta) = gpts(:ntheta)/t1
      gpts(:ntheta) = 0.5_dp*gpts(:ntheta)/gtt(:ntheta)
      t1 = MAXVAL(ABS(gpts(:ntheta)))
      t2 = 1.e-3_dp            !Arbitrary threshold on angle motion
      IF (t1.gt.t2) THEN
!     Careful not to make angle points move too rapidly so they DO not cross
        t1 = t2/t1
        gpts(:ntheta) = t1 * gpts(:ntheta)
      END IF
      fsq = 0.5_dp*dnorm*SUM(r1(:ntheta)**2 + z1(:ntheta)**2)
      fsq = SQRT(fsq)/r10

!
!     COMPUTE CURVE FORCES
!     1. Magnetic axis forces
!
      m = 0
      gr0c = dnorm*SUM(cosa(:ntheta,m)*r1(:ntheta))
      gz0c = dnorm*SUM(cosa(:ntheta,m)*z1(:ntheta))

!
!     RHO FORCES
!     (NOTE: RHO for m=0 IS EQUIVALENT TO THE L(v)-TERM => ARCLENGTH, Eq.14 H-B paper)
!
      DO m = 0, mrho1
         IF (m .le. 1) THEN
            t1 = dnorm/t1m(m+1)
            grhoc(m) = t1*SUM(cosa(:ntheta,m+1)*r1(:ntheta) +
     1                     sina(:ntheta,m+1)*z1(:ntheta))
            grhos(m) = t1*SUM(sina(:ntheta,m+1)*r1(:ntheta) -
     1                         cosa(:ntheta,m+1)*z1(:ntheta))
         ELSE
            t1 = t1m(m+1)
            t2 = t2m(m-1)
            tnorm = dnorm/(t1*t1 + t2*t2)
            t1 = t1*tnorm
            t2 = t2*tnorm
            grhoc(m) = SUM((cosa(:ntheta,m+1)*r1(:ntheta) +
     1                      sina(:ntheta,m+1)*z1(:ntheta))*t1 +
     2                     (cosa(:ntheta,m-1)*r1(:ntheta) -
     3                      sina(:ntheta,m-1)*z1(:ntheta))*t2)
            grhos(m) = SUM((sina(:ntheta,m+1)*r1(:ntheta) -
     1                      cosa(:ntheta,m+1)*z1(:ntheta))*t1 +
     2                     (sina(:ntheta,m-1)*r1(:ntheta) +
     3                      cosa(:ntheta,m-1)*z1(:ntheta))*t2)
         ENDIF
      END DO

      grhos(0) = zero        !Enforce constraint on toroidal angle

      gnorm = SUM(grhoc(0:mrho1)*grhoc(0:mrho1) +
     1            grhos(0:mrho1)*grhos(0:mrho1)) +
     2            gr0c**2 + gz0c**2
      gnorm = gnorm/r10**2

      gnorm = gnorm + dnorm*SUM(gpts(:ntheta)*gpts(:ntheta))

      END SUBROUTINE funct
