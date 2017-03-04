      SUBROUTINE scrunch(rin, zin, rbc, zbs, rbs, zbc,
     1   rmnaxis, zmnaxis, ftol, niter, nstep, mexp, mntot)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE Vname0
      USE Vname1
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: niter, nstep, mexp, mntot
      REAL(rprec) :: ftol
      REAL(rprec), DIMENSION(ntheta,nphi) :: rin, zin
      REAL(rprec), DIMENSION(mntot) :: rbc, zbs, rbs, zbc
      REAL(rprec), DIMENSION(*) ::  rmnaxis, zmnaxis
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: pexp = 4, qEXP = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: irst=1, nplane, imodes, iter, modeno, mrho1
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: angle
      REAL(rprec), DIMENSION(2*mpol,nphi) :: result1
      REAL(rprec) :: too_large,
     1   gtrig, fsq, gmin, g11, gout, time_on, time_off
C-----------------------------------------------
c
c       PARAMETER VALUES TO BE SET IN NAME1.INC:
c
c       MU:   MAXIMUM NO. OF POLOIDAL MODES IN FIT
c       NU:   MAXIMUM NO. OF THETA (POLOIDAL) POINTS ALLOWED IN FIT
c       NV:   MAXIMUM NO. OF TOROIDAL PLANES IN FIT
c
!

!
!     CHECK THAT INPUT DIMENSIONS ARE CONSISTENT WITH NAME1 PARAMETERS
!
      deltf = 0.97_dp
      IF (MOD(nphi,2).ne.0 .and. nphi.ne.1 .or. nphi.eq.0) THEN
         PRINT *, ' NPHI > 0 must be EVEN for non-symmetric systems'
         STOP
      ENDIF
      IF (ntheta > nu) THEN
         PRINT *, ' NTHETA input to SCRUNCH is too large'
         STOP
      ENDIF
      IF (nphi > nv) THEN
         PRINT *, ' NPHI input to SCRUNCH is too large'
         STOP
      ENDIF
      IF (mpol > mu) THEN
         PRINT *, ' MPOL input to SCRUNCH is too large'
         STOP
      ENDIF
      IF (nfp .le. 0) THEN
         PRINT *, ' NFP > 0 must be positive and exceed zero'
         STOP
      ENDIF


!
!     INITIALIZE FIXED m,n ARRAYS
!
      OPEN(unit=3, file='outcurve', status='unknown')
      twopi = 8*ATAN(one)
      CALL fixaray (pexp, qexp, mexp)

!
!     COMPUTE INITIAL GUESSES (MUST ORDER ANGLES MONOTONICALLY FOR
!     NON-STARLIKE DOMAINS)
!
      ALLOCATE (angle(ntheta, nphi))
      CALL inguess (rin, zin, angle)

      too_large = 1.E6_dp
      gtrig = 1.E-4_dp
      mrho1 = mrho+1

      ALLOCATE (xvec(n2), gvec(n2), xdot(n2), xstore(n2))
!
!     BEGIN MAIN INTERATION LOOP
!
      CALL second0(time_on)
      DO nplane = 1, nphi
!
!     INITIALIZE M = 0 and M = 1 MODE AMPLITUDES
!
!        STACKING OF XVEC ARRAY (FOR FIXED TOROIDAL PLANE)
!        XVEC(1)              : R0c             (Symmetric components)
!        XVEC(2,mrho1)        : Rhoc            (Symmetric components)
!        XVEC(1+mrho1)        : Z0c             (Asymmetric components)
!        XVEC(2+mrho1,2*mrho1): Rhos            (Asymmetric components)
!        XVEC(2*mrho1+1,n2)   : Theta angle
!
         xstore(:n2) = zero
         xdot(:n2) = zero
         xvec(:n2) = zero

         WRITE (3, 10) nplane
         WRITE (*, 10) nplane
 10      FORMAT(/'                  Fitting toroidal plane # ',i3)

         CALL amplitud (r0n(nplane), z0n(nplane), angle(1,nplane),
     1     xvec, xvec(1+mrho1), xvec(2:mrho1), xvec(2+mrho1:2*mrho1),
     2     xvec(2*mrho1+1:), rin(1,nplane), zin(1,nplane))

      WRITE (3, 20)
      WRITE (*, 20)
 20   FORMAT(/' ITERATIONS    RMS ERROR    FORCE GRADIENT    <M>',
     1        '    MAX m   DELT')

         imodes = mrho
         delt = deltf

         DO iter = 1, niter
            gvec(:n2) = zero
            CALL funct (xvec, xvec(1+mrho1), xvec(2:mrho1),
     1        xvec(2+mrho1:2*mrho1),xvec(2*mrho1+1:), gvec,
     2        gvec(1+mrho1), gvec(2:mrho1), gvec(2+mrho1:2*mrho1),
     3        gvec(2*mrho1+1:), fsq, rin(1,nplane), zin(1,nplane),
     4        imodes)

            SELECT CASE (iter)
            CASE DEFAULT
               gmin = MIN(gmin,gnorm)
               CALL evolve (g11)
!
!     RUDIMENTARY TIME STEP CONTROL
!
               IF (gnorm/gmin .gt. too_large) irst = 2
               IF (irst.eq.2 .or. gmin.eq.gnorm) CALL restart (irst)

               IF (gnorm.lt.gtrig .and. imodes.lt.mrho) THEN
                  imodes = imodes + 1
                  CALL restart (irst)
                  irst = 2
                  delt = delt/.95
               ENDIF
               IF (MOD(iter,nstep).ne.0 .and. gnorm.gt.ftol**2) CYCLE
            CASE (2)
               g11 = gnorm
               CYCLE
            CASE (1)
               gmin = gnorm
               imodes = MIN(4,mrho)
            END SELECT

            gout = SQRT(gnorm)
            modeno = imodes
            IF (iter .eq. 1) modeno = mpol - 1
            IF (qexp .gt. zero) specw = specw**(one/qexp)
            PRINT 110, iter, fsq, gout, specw, modeno, delt
            WRITE (3, 110) iter, fsq, gout, specw, modeno, delt
            IF (gnorm.lt.ftol**2 .and. imodes.eq.mrho
     1         .or. fsq.lt.ftol) EXIT
         END DO

  110    FORMAT(i8,1p,2e16.3,0p,f10.2,i8,1p,e10.2)

!
!     STORE FINAL ANSWER FOR FINAL PHI-TRANSFORM OF RHOC, RHOS
!
         result1(1:2*mpol,nplane) = xvec(1:2*mpol)

      END DO
      CALL second0(time_off)
!
!     OUTPUT LOOP
!
      PRINT 330, mexp, pexp, qexp, time_off-time_on
      WRITE (3, 330) mexp, pexp, qexp, time_off-time_on
  330 FORMAT(/' ANGLE CONSTRAINTS WERE APPLIED ',/,
     1   ' BASED THE POLAR DAMPING EXPONENT (PEXP) = ',i2,/,
     2   ' RM**2 + ZM**2 SPECTRUM COMPUTED WITH P = ',f8.2,
     3   ' AND Q = ',f8.2,/,' TIME: ',1p,e10.2,
     3   ' SEC.'/)
!
!     PERFORM PHI-FOURIER TRANSFORM
!
      CALL fftrans (result1(1,1:nphi), result1(1+mpol,1:nphi),
     1   result1(2:mpol,1:nphi), result1(mpol+2:2*mpol,1:nphi), rbc,
     1   zbs, rbs, zbc, rmnaxis, zmnaxis)

      DEALLOCATE (xvec, gvec, xdot, xstore, angle)

      END SUBROUTINE scrunch
