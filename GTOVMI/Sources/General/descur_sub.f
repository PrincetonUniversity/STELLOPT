      SUBROUTINE descur_sub(g1,g2)
c       THIS IS PROGRAM DESCUR - WHICH USES A STEEPEST DESCENT
c       ALGORITHM TO FIND A LEAST SQUARES APPROXIMATION TO AN
c       ARBITRARY 3-D SPACE CURVE. ANGLE CONSTRAINTS, BASED ON A
c       MINIMUM SPECTRAL WIDTH CRITERION, ARE APPLIED.
c       THE CONSTRAINT IS SATISFIED BY TANGENTIAL VARIATIONS
c       ALONG THE CURVE.
c
c       THE MAIN PROGRAM SETS UP THE INITIAL POINTS TO BE FITS AND
c       THEN CALLS THE SUBROUTINE SCRUNCH, WHICH DOES THE ACTUAL
c       CURVE-FITTING.
c
c***********************************************************************
c       REVISION 1:  January 26, 1989
c                    Compute Fit to individual toroidal planes separately
c                    and THEN Fourier decompose in phi, using optimized
c                    angle representation
c
c       REVISION 2:  July, 1995
c                    Generalized for up-down asymmetric CASE
c
c       REVISION 3:  July, 1997
c                    Converted to Fortran-90
c                    Eliminated Lagrange multiplier constraint
c                    by using an EXPlicit representation
c
c
c       PARAMETER VALUE INTERPRETATIONS:
c
c       MPOL:   ACTUAL NO. OF POLOIDAL MODES USED IN FIT FOR R and Z
c       MRHO:   ACTUAL NO. OF POLOIDAL MODES IN FIT FOR QUASI-POLAR RHO
c       NTHETA: ACTUAL NO. OF THETA (POLOIDAL) POINTS AUSED IN FIT
c       NPHI:   ACTUAL NO. OF TOROIDAL PLANES IN FIT
c
c***********************************************************************
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE Vname0
      USE Vname1
      USE gfile
      USE mapout
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: nphi20 = 1 + nv/2
      INTEGER, PARAMETER :: mntot = mu*(2*nphi20 + 1)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: niter, nstep, itype,
     1   i, k, j, m, n, mexp
      REAL(rprec), DIMENSION(nu*nv) :: rin, zin, LENgth
      REAL(rprec), DIMENSION(mntot) :: rbc, rbs, zbc, zbs
      REAL(rprec), DIMENSION(0:nv) :: rmnaxis, zmnaxis
      REAL(rprec), DIMENSION(0:mu-1,-nphi20:nphi20) ::
     1 rbdy_3d, zbdy_3d
      REAL(rprec), DIMENSION(4) :: rbdy, zbdy, rdshape, zdshape,
     1 rbean, zbean, rbelt, zbelt, rellip, zellip, rsqr, zsqr
      REAL(rprec), DIMENSION(0:3,0:2) :: rheliac, zheliac
      REAL(rprec) :: ftol, pexp, qexp, zeta, texp,
     1   dumr, dumphi, dumz, arg, PSI, XI, DELTA_R,
     2   rmin, rmax, denom
      CHARACTER :: datafile*80
      CHARACTER*1 :: ch_yn
C-----------------------------------------------

      TYPE(GEQDSK) :: g1
      TYPE(RSZS) :: g2
      data rdshape/3.0, 0.991, 0.136, 0./
      data zdshape/0.0, 1.409, -0.118, 0./
      data rbean/3.0, 1.042, 0.502,  - 0.0389/
      data zbean/0.0, 1.339, 0.296,  - 0.0176/
      data rbelt/3.0, 0.453, 0., 0./
      data zbelt/0.0, 0.60, 0., 0.196/
      data rellip/3.0E0, 1.000E0, 0., 0./
      data zellip/0.0, 6.00E0, 0., 0./
      data rheliac/4.115, 0.4753, 0., 0., 0.3225, -0.06208, -0.1358,
     1   -0.04146, 0.0136, 0.0806, -0.0205, 0.0445/
      data zheliac/0., -0.5045, 0., 0., 0.337, 0.0637, 0.1325, -0.04094
     1   , 0.01152, 0.06186, -0.03819,  - 0.02366/
      data rsqr/    3.0, 0.4268, 0., 0.07322/,
     1     zsqr/    0.0, 0.4268, 0., -0.07322/

c**********************************************************************
c       CONTROL DATA THAT CAN BE SPECIFIED BY USER
c**********************************************************************
      data niter/1500/
      data nstep/100/

      data ftol/1.E-5_dbl/
      data pexp/4.0_dbl/
      data qexp/1.0_dbl/
      data mexp/4/

      data datafile/'none'/

      twopi = 8*ATAN(one)

c**********************************************************************
c       CREATE ARRAYS OF CURVE POINTS, RIN(THETA,PHI) , ZIN(THETA,PHI),
c       TO BE FIT TO AN OPTIMIZED FOURIER SERIES BY CALLING "SCRUNCH".
c       NOTE: N IS IN UNITS NORMED TO NFP (NUMBER OF FIELD PERIODS) IN 20 LOOP
c**********************************************************************
c
c     Read in data
c

      HB_Parameter = 1

      nfp = 1

      iTYPE = 0

 1002 CONTINUE
      ntheta = g2%nthet ; nphi = 1 ; nfp = 1
      IF (ntheta .gt. nu) STOP 'ntheta > nu'
      IF (nphi .gt. nv) STOP 'nphi > nv'
      mpol = mu
      mpol1 = mpol - 1
      nphi2 = 1 + nphi/2
      i = 0 ; rin = 0 ; zin = 0
!      rin(1:ntheta)=g2%rs(1:ntheta,1)
!      zin(1:ntheta)=g2%zs(1:ntheta,1)
      rin(1:ntheta)=g2%rs(g2%npsi,1:ntheta)
      zin(1:ntheta)=g2%zs(g2%npsi,1:ntheta)
 4000 CONTINUE
c**********************************************************************
c     PERFORM OPTIMIZING CURVE-FIT
c**********************************************************************
      CALL scrunch (rin, zin, pexp, qexp, rbc, zbs, rbs, zbc, rmnaxis,
     1   zmnaxis, ftol, niter, nstep, mexp)

c**********************************************************************
c     PRINT OPTIMIZED TRANSFORM COEFFICIENTS TO "OUTCURVE" FILE
c     AND STORE RIN, ZIN, RBC, ZBS DATA IN PLOTOUT FILE
c**********************************************************************
      CALL PRINTit(rin,zin,rbc,zbs,rbs,zbc,rmnaxis,zmnaxis,g2)
      CLOSE(10)

      END SUBROUTINE descur_sub
      SUBROUTINE scrunch(rin, zin, pexp, qexp, rbc, zbs, rbs, zbc,
     1   rmnaxis, zmnaxis, ftol, niter, nstep, mexp)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE Vname0
      USE Vname1
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: niter, nstep, mexp
      REAL(rprec) pexp, qexp, ftol
      REAL(rprec), DIMENSION(ntheta,nphi) :: rin, zin
      REAL(rprec), DIMENSION(*) :: rbc, zbs, rbs, zbc,
     1  rmnaxis, zmnaxis
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: irst=1, nplane, imodes, iter, MODeno, mrho1
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: angle
      REAL(rprec), DIMENSION(2*mpol,nphi) :: result
      REAL(rprec) :: too_large,
     1   gtrig, fsq, gmin, g11, gout
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
      deltf = 0.97_dbl
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
CEAL      OPEN(unit=3, file='outcurve', status='unknown')
      twopi = 8*ATAN(one)
      CALL fixaray (pexp, qexp, mexp)

!
!     COMPUTE INITIAL GUESSES (MUST ORDER ANGLES MONOTONICALLY FOR
!     NON-STARLIKE DOMAINS)
!
      ALLOCATE (angle(ntheta, nphi))
      CALL inguess (rin, zin, angle)

      too_large = 1.E6_dbl
      gtrig = 1.E-4_dbl
      mrho1 = mrho+1

      ALLOCATE (xvec(n2), gvec(n2), xdot(n2), xstore(n2))
!
!     BEGIN MAIN INTERATION LOOP
!
      DO nplane = 1, nphi
!
!     INITIALIZE M = 0 and M = 1 MODE AMPLITUDES
!
!        STACKING OF XVEC ARRAY (FOR FIXED TOROIDAL PLANE)
!        XVEC(1)              : R0c              (Symmetric components)
!        XVEC(2,mrho1)        : Rhoc              (Symmetric components)
!        XVEC(1+mrho1)        : Z0c             (Asymmetric components)
!        XVEC(2+mrho1,2*mrho1): Rhos            (Asymmetric components)
!        XVEC(2*mrho1+1,n2)   : Theta angle
!
         xstore(:n2) = zero
         xdot(:n2) = zero
         xvec(:n2) = zero
         r0n = zero ; z0n = zero ; raxis = zero ; zaxis = zero 
CEAL         WRITE (3, 10) nplane
CEAL         WRITE (*, 10) nplane
 10      FORMAT(/'                  Fitting toroidal plane # ',i3)

         CALL amplitud (r0n(nplane), z0n(nplane), angle(1,nplane),
     1     xvec(1), xvec(1+mrho1), xvec(2:mrho1), xvec(2+mrho1:2*mrho1),
     2     xvec(2*mrho1+1:), rin(1,nplane), zin(1,nplane))

CEAL      WRITE (3, 20)
CEAL      WRITE (*, 20)
 20   FORMAT(/' ITERATIONS    RMS ERROR    FORCE GRADIENT    <M>',
     1        '    MAX m   DELT')

         imodes = mrho
         delt = deltf

         DO iter = 1, niter
            gvec(:n2) = zero
            CALL funct (xvec(1), xvec(1+mrho1), xvec(2:mrho1),
     1        xvec(2+mrho1:2*mrho1),xvec(2*mrho1+1:), gvec(1),
     2        gvec(1+mrho1), gvec(2:mrho1), gvec(2+mrho1:2*mrho1),
     3        gvec(2*mrho1+1:), fsq, rin(1,nplane), zin(1,nplane),
     4        imodes)

            SELECT CASE (iter)
            CASE DEFAULT
               gmin = min(gmin,gnorm)
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
               imodes = min(4,mrho)
            END SELECT

            gout = SQRT(gnorm)
            MODeno = imodes
            IF (iter .eq. 1) MODeno = mpol - 1
            IF (qEXP .gt. zero) specw = specw**(one/qexp)
CEAL            PRINT 110, iter, fsq, gout, specw, MODeno, delt
CEAL            WRITE (3, 110) iter, fsq, gout, specw, MODeno, delt
            IF (gnorm.lt.ftol**2 .and. imodes.eq.mrho
     1         .or. fsq.lt.ftol) EXIT
         END DO

  110    FORMAT(i8,1p,2e16.3,0p,f10.2,i8,1p,e10.2)

!
!     STORE FINAL ANSWER FOR FINAL PHI-TRANSFORM OF RHOC, RHOS
!
         result(1:2*mpol,nplane) = xvec(1:2*mpol)

      END DO
!
!     OUTPUT LOOP
!
CEAL      PRINT 330, pexp, qexp, mexp
CEAL      WRITE (3, 330) pexp, qexp, mexp
  330 FORMAT(/' ANGLE CONSTRAINTS WERE APPLIED ',/,
     1   ' BASED ON RM**2 + ZM**2 SPECTRUM WITH P = ',f8.2,' AND Q = ',
     2   f8.2,/,' POLAR DAMPING EXPONENT = ',i2/,' TIME: ',1pe10.2,
     3   ' SEC.'/)
!
!     PERFORM PHI-FOURIER TRANSFORM
!
      CALL fftrans (result(1,1:nphi), result(1+mpol,1:nphi),
     1   result(2:mpol,1:nphi), result(mpol+2:2*mpol,1:nphi), rbc,
     1   zbs, rbs, zbc, rmnaxis, zmnaxis)

      DEALLOCATE (xvec, gvec, xdot, xstore, angle)

      END SUBROUTINE scrunch

      SUBROUTINE fixaray(pexp, qexp, mexp)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE Vname0
      USE Vname1
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: mexp
      REAL(rprec) :: pexp, qexp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: l, ntor, nn0, n, m
      REAL(rprec) :: sgn
C-----------------------------------------------
!
!       This routine stores toroidal and poloidal MODe number arrays
!
!       MPOL = NUMBER OF POLOIDAL MODES USED IN CURVE-FIT
!       NPHI = NUMBER OF EQUALLY-SPACED PHI PLANES PER FIELD
!              PERIOD IN WHICH THE CURVE-FITTING IS PERFORMED
!       IMPORTANT: NPHI MUST BE EVEN FOR NON-SYMMETRIC SYSTEMS
!
!       MPNT = NUMBER OF R AND Z MODES APPEARING IN FIT
!       NTHETA = NUMBER OF THETA POINTS IN EACH PHI PLANE
!       N2   = TOTAL NUMBER OF RESIDUALS IN FSQ (PER PHI PLANE)
!            == 2 (m=0 MODes) + 2*mrho (rhoc, rhos) + ntheta
!       NFP  = NUMBER OF TOROIDAL FIELD PERIODS (IRRELEVANT)
!
      mrho = mpol-1
      nphi2 = nphi/2
      IF (nphi2 .eq. 0) nphi2 = 1
      n2 = 2*(mrho+1) + ntheta

      l = 0
      ntor = max0(1,nphi - 1)
      IF (nphi .eq. 2) ntor = 2

      nn0 = 1 - (ntor + 1)/2
      DO n = 1, ntor
         nn(n) = (nn0 + (n - 1))*nfp
      END DO
      DO m = 1, mpol
         mm(m) = m - 1
         DO n = 1, ntor
            IF (mm(m).ne.0 .or. nn(n).ge.0) THEN
               l = l + 1
               m1(l) = mm(m)
               n1(l) = nn(n)
            ENDIF
         END DO
      END DO
      mpnt = l
      dnorm = 2._dbl/ntheta

      DO m = 0, mpol
         dm1(m) = m
         IF (m .eq. 0) THEN
            xmpq(m,1) = zero
            xmpq(m,2) = zero
         ELSE
            xmpq(m,1) = dm1(m)**(pexp+qexp)
            xmpq(m,2) = dm1(m)**(pexp)
         ENDIF
      END DO

!
!     FOR SECOND DERIVATIVE-TYPE CONSTRAINT, CHOOSE SGN = -one
!
      sgn = one                            !First derivative-TYPE constraint
!     sgn = -one                      !Second derivative-TYPE constraint
      t1m(1) = one
      DO m = 2, mrho
        t1m(m) = (REAL(m-1, rprec)/m)**mexp
      END DO
      DO m = 1,mrho-2
        t2m(m) = sgn*(REAL(m+1, rprec)/m)**mexp
      END DO
      t2m(mrho-1) = zero
      t2m(mrho) = zero

      END SUBROUTINE fixaray
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
c       This SUBROUTINE obtains initial guesses for the centroid at each
c       toroidal plane. By DEFAULT, the polar axis is set equal to this
c       centroid.  It is imperative that the polar axis be enclosed by
c       the surface (otherwise the computed theta angles will not span
c       [0,2pi]). For certain complex cross-sections, this simple estimate
c       of the polar axis may fail, and the USEr must supply values for
c       raxis(i),zaxis(i).  In addition, for non-starlike DOmains, the
c       points along the curve must be monitonically increasing as a
c       FUNCTION of the arclength from ANY fixed poINT on the curve. This
c       ordering is attempted by the SUBROUTINE ORDER, but may fail in
c       general IF too few theta points are USEd.
c
c Modified 8/13/97, J. Breslau. The SUBROUTINE now tries several guesses
c  for the magnetic axis position IF the first fails. The mean position
c  of pairs of points is taken until one is found to be inside.
c**********************************************************************
c       COMPUTE CENTROID AND
c       ORDER POINTS ON FLUX SURFACE AT EACH TOROIDAL ANGLE
c**********************************************************************
c      PRINT *, 'ORDERING SURFACE POINTS'
      r0n = 0 ; z0n = 0
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
               PRINT *, 'Could not find INTernal point'
               WRITE (3, *) 'Could not find INTernal point'
               STOP
            ENDIF
            j1 = 1
            j2 = j1 + jskip
   50       CONTINUE
            raxis(i) = 0.5_dbl*(rin(j1,i)+rin(j2,i))
            zaxis(i) = 0.5_dbl*(zin(j1,i)+zin(j2,i))
            CALL order (rin(1,i), zin(1,i), raxis(i), zaxis(i), inside)
            IF (inside .eq. 1) THEN
               r0n(i) = raxis(i)
               z0n(i) = zaxis(i)
            ELSE
               j1 = j1 + 1
               IF (j1 .gt. ntheta) go to 40
               j2 = j2 + 1
               IF (j2 .gt. ntheta) j2 = 1
               go to 50
            ENDIF
         ENDIF
      END DO
c**********************************************************************
c       COMPUTE OPTIMUM ANGLE OFFSETS FOR M = 1 MODES
c**********************************************************************
      CALL getangle (rin, zin, angle, r0n, z0n)

      END SUBROUTINE inguess
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
c       Fourier MODe amplitudes to the appropriate components of
c       the xvec array
c*****************************************************************
      r0c = rcenter
      z0c = zcenter
      xpts(:ntheta) = angin(:ntheta)

      mrz = mpol-1
      xmult = 2._dbl/ntheta

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
CEAL      WRITE(3,'(/,3(a,1pe10.3))')
CEAL     1    ' ]RAXIS = ', rcenter,' ZAXIS = ', zcenter,' R10 = ',r10
CEAL      WRITE(*,'(/,3(a,1pe10.3))')
CEAL     1    ' RAXIS = ', rcenter,' ZAXIS = ', zcenter,' R10 = ',r10

*
*     INITIALIZE POLAR RHOs for m = 0 thru mrz-1
*
      DO m = 0, mrz-1
         IF( m.le.1 )then
            t1 = t1m(m+1)
            rhoc(m) = 0.5_dbl*(r1c(m+1) + z1s(m+1))/t1
            rhos(m) = 0.5_dbl*(r1s(m+1) - z1c(m+1))/t1
         ELSE
            t1 = t1m(m+1)
            t2 = t2m(m-1)
            tnorm = 0.5_dbl/(t1**2 + t2**2)
            rhoc(m) = tnorm*( (r1c(m+1) + z1s(m+1))*t1
     1                   +  (r1c(m-1) - z1s(m-1))*t2 )
            rhos(m) = tnorm*( (r1s(m+1) - z1c(m+1))*t1
     1                   +  (r1s(m-1) + z1c(m-1))*t2 )
      END IF
      END DO

      END SUBROUTINE amplitud
      SUBROUTINE funct(r0c, z0c, rhoc, rhos, xpts, gr0c, gz0c, grhoc,
     1              grhos, gpts, fsq, xin, yin, mrho_in)
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
      REAL(rprec), DIMENSION(ntheta,0:mrho_in) :: COSa, SINa
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
      COSa(:,0) = one
      COSa(:,1) = COS(xpts)
      SINa(:,0) = zero
      SINa(:,1) = SIN(xpts)

!
!     COMPUTE CURVE R1 = R - Rin, Z1 = Z - Zin FOR INPUT POINTS
!     NOTE DIMENSIONS: Rm(0:MRHO), Zm(0:MRHO), but RHO(0:MRHO-1)
!
      DO m = 2, mrho_in
         COSa(:,m) = COSa(:,m-1)*cosa(:,1) - SINa(:,m-1)*sina(:,1)
         SINa(:,m) = SINa(:,m-1)*cosa(:,1) + COSa(:,m-1)*sina(:,1)
      END DO

      DO m = 0, mrho_in
         CALL getrz(rmc_p,rms_p,zmc_p,zms_p,r0c,z0c,rhoc,rhos,
     1      m,mrho_in)
!
!        COMPUTE SPECTRAL WIDTH OF CURVE
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
      gpts(:ntheta) = 0.5_dbl*gpts(:ntheta)/gtt(:ntheta)
      t1 = MAXVAL(ABS(gpts(:ntheta)))
      t2 = 1.e-3_dbl            !Arbitrary threshold on angle motion
      IF (t1.gt.t2) THEN
!     Careful not to make angle points move too rapidly so they DO not cross
        t1 = t2/t1
        gpts(:ntheta) = t1 * gpts(:ntheta)
      END IF
      fsq = 0.5_dbl*dnorm*SUM(r1(:ntheta)**2 + z1(:ntheta)**2)
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
     1                     SINa(:ntheta,m+1)*z1(:ntheta))
            grhos(m) = t1*SUM(sina(:ntheta,m+1)*r1(:ntheta) -
     1                         COSa(:ntheta,m+1)*z1(:ntheta))
         ELSE
            t1 = t1m(m+1)
            t2 = t2m(m-1)
            tnorm = dnorm/(t1*t1 + t2*t2)
            t1 = t1*tnorm
            t2 = t2*tnorm
            grhoc(m) = SUM((cosa(:ntheta,m+1)*r1(:ntheta) +
     1                      SINa(:ntheta,m+1)*z1(:ntheta))*t1 +
     2                     (cosa(:ntheta,m-1)*r1(:ntheta) -
     3                      SINa(:ntheta,m-1)*z1(:ntheta))*t2)
            grhos(m) = SUM((sina(:ntheta,m+1)*r1(:ntheta) -
     1                      COSa(:ntheta,m+1)*z1(:ntheta))*t1 +
     2                     (sina(:ntheta,m-1)*r1(:ntheta) +
     3                      COSa(:ntheta,m-1)*z1(:ntheta))*t2)
         ENDIF
      END DO

      grhos(0) = zero        !Enforce constraINT on toroidal angle

      gnorm = SUM(grhoc(0:mrho1)*grhoc(0:mrho1) +
     1            grhos(0:mrho1)*grhos(0:mrho1)) +
     2            gr0c**2 + gz0c**2
      gnorm = gnorm/r10**2

      gnorm = gnorm + dnorm*SUM(gpts(:ntheta)*gpts(:ntheta))

      END SUBROUTINE funct
      SUBROUTINE evolve(g11)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE Vname0
      USE Vname1
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec) g11
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: bmax = 0.15_dbl
      REAL(rprec) :: ftest, dtau, otav, b1, fac
C-----------------------------------------------

      ftest = gnorm/g11
      dtau = ABS(one - ftest)
      g11 = gnorm
      otav = dtau/delt
      dtau = delt*otav + 1.e-3_dbl
      dtau = min(bmax,dtau)
      b1 = one - 0.5_dbl*dtau
      fac = one/(one + 0.5_dbl*dtau)
      xdot(:n2) = fac*(xdot(:n2)*b1-delt*gvec(:n2))
      xvec(:n2) = xvec(:n2) + xdot(:n2)*delt

      END SUBROUTINE evolve
      SUBROUTINE fftrans(r0c, z0c, rhoc, rhos, rbc, zbs, rbs, zbc,
     1   rmnaxis, zmnaxis)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE Vname0
      USE Vname1
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(nphi) :: r0c, z0c
      REAL(rprec), DIMENSION(0:mrho-1,nphi) :: rhoc, rhos
      REAL(rprec), DIMENSION(0:mpol-1,-nphi2:nphi2) ::
     1  rbc, zbs, rbs, zbc
      REAL(rprec), DIMENSION(0:nphi2) :: rmnaxis, zmnaxis
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, mn, mREAL, nreal
      REAL(rprec), DIMENSION(nv) :: INTgrate, argi
      REAL(rprec) :: delphi,dn,rmc_p,zms_p,rms_p,zmc_p,
     1  arg,tcosn,tsinn
C-----------------------------------------------
!
!       PERFORM FOURIER TRANSFORM IN phi
!
      delphi = one/nphi
      DO i = 1, nphi
         INTgrate(i) = delphi
         argi(i) = twopi*(i - 1)/REAL(nphi*nfp,rprec)
      END DO

      rbc = 0;   zbs = 0;   rbs = 0;   zbc = 0

      DO mn = 1, mpnt
         mREAL = m1(mn)
         nREAL = n1(mn)/nfp
         dn = REAL(n1(mn))
         DO i = 1, nphi
            CALL getrz(rmc_p,rms_p,zmc_p,zms_p,r0c(i),z0c(i),
     1         rhoc(0,i),rhos(0,i),mREAL,mrho)
            arg = dn*argi(i)
            tcosn = COS(arg)
            tsinn = SIN(arg)
            rbc(mREAL,nreal) = rbc(mreal,nreal) + INTgrate(i)*(tcosn*
     1         rmc_p + tsinn*rms_p)
            zbs(mREAL,nreal) = zbs(mreal,nreal) + INTgrate(i)*(tcosn*
     1         zms_p - tsinn*zmc_p)
            zbc(mREAL,nreal) = zbc(mreal,nreal) + INTgrate(i)*(tcosn*
     1         zmc_p + tsinn*zms_p)
            rbs(mREAL,nreal) = rbs(mreal,nreal) + INTgrate(i)*(tcosn*
     1         rms_p - tsinn*rmc_p)
         END DO
         IF (mreal.eq.0 .and. nreal.eq.0) THEN
            rmnaxis(0) = DOT_PRODUCT(intgrate(:nphi),raxis(:nphi))
            zmnaxis(0) = zero
         ELSE IF (mreal.eq.0 .and. nreal.gt.0) THEN
            rbc(0,nreal) = 2*rbc(0,nreal)
            rbs(0,nreal) = 2*rbs(0,nreal)
            zbc(0,nreal) = 2*zbc(0,nreal)
            zbs(0,nreal) = 2*zbs(0,nreal)
            rmnaxis(nreal) = zero
            zmnaxis(nreal) = zero
            rmnaxis(nreal) = rmnaxis(nreal) + SUM(2*intgrate(:nphi)*
     1         raxis(:nphi)*COS(dn*argi(:nphi)))
            zmnaxis(nreal) = zmnaxis(nreal) - SUM(2*intgrate(:nphi)*
     1         zaxis(:nphi)*SIN(dn*argi(:nphi)))
         ENDIF
      END DO

      END SUBROUTINE fftrans
      SUBROUTINE PRINTit(rin,zin,rbc,zbs,rbs,zbc,rmnaxis,zmnaxis,g2)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE Vname0
      USE Vname1
      USE mapout
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(*) :: rin, zin
      REAL(rprec), DIMENSION(0:mpol - 1,-nphi2:nphi2) ::
     1   rbc, zbs, rbs, zbc
      REAL(rprec), DIMENSION(0:nphi2) :: rmnaxis, zmnaxis
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, m, n, n11
      REAL(rprec) :: tol
      CHARACTER*250 :: form_string
      TYPE(RSZS) :: g2
C-----------------------------------------------
      g2%rbc=0 ; g2%rbs=0 ; g2%zbs=0 ; g2%zbc=0 ;
      g2%rbc(0:mu-1,0)=rbc(0:mu-1,0)
      g2%rbs(0:mu-1,0)=rbs(0:mu-1,0)
      g2%zbc(0:mu-1,0)=zbc(0:mu-1,0)
      g2%zbs(0:mu-1,0)=zbs(0:mu-1,0)

      OPEN(unit=10, file='plotout', status='unknown')
      WRITE (10, 1990) mpol, ntheta, nphi, mpol - 1, nphi2, nfp, mpnt
 1990 FORMAT(7i10)
      DO i = 1, nphi*ntheta
         WRITE (10, 1995) rin(i), zin(i)
      END DO
 1995 FORMAT(1p,2e12.4)

c**********************************************************************
c       This SUBROUTINE PRINTs out data INTo the file "outcurve"
c       The FORMAT of the MODes is compatible with input to VMEC
c**********************************************************************
CEAL      WRITE (3, 10)
   10 FORMAT(/'   MB  NB      RBC         RBS         ',
     1   'ZBC         ZBS        RAXIS       ZAXIS')
      tol = 1.E-6_dbl*ABS(rbc(1,0))
      DO m = 0, mpol - 1
         DO n = -nphi2, nphi2
CEAL            WRITE (10, 2000) rbc(m,n), zbs(m,n), rbs(m,n), zbc(m,n)
            IF (.not.(ABS(rbc(m,n)).lt.tol .and. ABS(zbs(m,n)).lt.tol
     1          .and. ABS(rbs(m,n)).lt.tol .and. ABS(zbc(m,n)).lt.tol))
     2      THEN
               IF (m.eq.0 .and. n.ge.0) THEN
CEAL                  WRITE (3, 30) m, n, rbc(m,n), rbs(m,n), zbc(m,n),
CEAL     1               zbs(m,n), rmnaxis(n), zmnaxis(n)
               ELSE
CEAL                  WRITE (3, 40) m, n, rbc(m,n), rbs(m,n), zbc(m,n),
CEAL     1               zbs(m,n)
               ENDIF
            ENDIF
         END DO
      END DO
   30 FORMAT(i5,i4,1p,6e12.4)
   40 FORMAT(i5,i4,1p,4e12.4)
 2000 FORMAT(1p,4e12.4)

c**********************************************************************
c     WRITE OUT IN FORMAT THAT CAN BE EXTRACTED INTO VMEC
c**********************************************************************
CEAL      WRITE(3, *)
      WRITE(67,fmt='(''  MPOL = '',i2.2)')mpol
      DO m = 0, mpol-1
         DO n = -nphi2, nphi2
            IF (.not.(ABS(rbc(m,n)).lt.tol .and. ABS(zbs(m,n)).lt.tol
     1          .and. ABS(rbs(m,n)).lt.tol .and. ABS(zbc(m,n)).lt.tol))
     2      THEN
               n11 = nphi2/10 + 1
               IF (n .lt. 0) n11 = n11+1
               WRITE(form_string,'(a,4(a,i1,a,i1,a))')"(2x,",
     1"'RBC(',i",n11,",',',i",m/10+1,",') = ',1p,e14.6,3x,",
     2"'RBS(',i",n11,",',',i",m/10+1,",') = ',e14.6,3x,",
     3"'ZBC(',i",n11,",',',i",m/10+1,",') = ',e14.6,3x,",
     4"'ZBS(',i",n11,",',',i",m/10+1,",') = ',e14.6)"
               WRITE(form_string,'(a,4(a,i1,a,i1,a))')"(2x,",
     1"'RBC(',i",n11,",',',i",m/10+1,",') = ',1p,e14.6,3x,",
     2"'RBS(',i",n11,",',',i",m/10+1,",') = ',e14.6,3x,/,4x,",
     3"'ZBC(',i",n11,",',',i",m/10+1,",') = ',e14.6,3x,",
     4"'ZBS(',i",n11,",',',i",m/10+1,",') = ',e14.6)"
               WRITE(67, form_string) n,m,rbc(m,n), n,m,rbs(m,n), 
     1               n,m,zbc(m,n), n,m,zbs(m,n)
c               IF (m.eq.0 .and. n.ge.0) THEN
c                  WRITE (3, 30) m, n, rbc(m,n), rbs(m,n), zbc(m,n),
c     1               zbs(m,n), rmnaxis(n), zmnaxis(n)
            END IF
         END DO
      END DO
 
      END SUBROUTINE PRINTit
      SUBROUTINE getrz(rmc,rms,zmc,zms,r0c,z0c,rhoc,rhos,m,mrho_in)
      USE Vname1
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER m, mrho_in
      REAL(rprec), DIMENSION(0:mrho_in-1) :: rhoc, rhos
      REAL(rprec) :: rmc, rms, zmc, zms, r0c, z0c
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: mrho1
C-----------------------------------------------

      mrho1 = mrho_in - 1
      rhos(0) = zero            !Enforce ConstraINT on toroidal angle

      IF (m .eq. 0) THEN
         rmc = r0c
         zmc = z0c
         rms = zero
         zms = zero
      ELSE IF (m .lt. mrho1) THEN
         rmc = (t1m(m)*rhoc(m-1) + t2m(m)*rhoc(m+1))
         zms = (t1m(m)*rhoc(m-1) - t2m(m)*rhoc(m+1))
         rms = (t1m(m)*rhos(m-1) + t2m(m)*rhos(m+1))
         zmc =-(t1m(m)*rhos(m-1) - t2m(m)*rhos(m+1))
      ELSE                      !Can change highest m constraints here...
         rmc = t1m(m)*rhoc(m-1) * HB_Parameter
         zms = rmc              * HB_Parameter
         rms = t1m(m)*rhos(m-1) * HB_Parameter
         zmc =-rms              * HB_Parameter
      ENDIF

      END SUBROUTINE getrz
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
      REAL(rprec), DIMENSION(nv) :: rcenter, zcenter
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, j, iterate
      REAL(rprec), DIMENSION(nv) ::
     1  rcos, rsin, zcos, zsin, phiangle
      REAL(rprec) :: xc, yc, dnum, denom, delangle
C-----------------------------------------------
c**********************************************************************
c       Compute angle offset consistent with constraINT Z1n = Z1,-n
c       Note: This is DOne iteratively, SINce elongation is unknown
c**********************************************************************
      DO i = 1, nphi
         DO j = 1, ntheta
            angle(j,i) = twopi*(j - 1)/REAL(ntheta,rprec)
         END DO
      END DO
      DO iterate = 1, 5
         DO i = 1, nphi
            rCOS(i) = zero
            rSIN(i) = zero
            zCOS(i) = zero
            zSIN(i) = zero
            DO j = 1, ntheta
               xc = rval(j,i) - rcenter(i)
               yc = zval(j,i) - zcenter(i)
               rCOS(i) = rCOS(i) + COS(angle(j,i))*xc
               rSIN(i) = rSIN(i) + SIN(angle(j,i))*xc
               zCOS(i) = zCOS(i) + COS(angle(j,i))*yc
               zSIN(i) = zSIN(i) + SIN(angle(j,i))*yc
            END DO
         END DO
c**********************************************************************
c       Compute new angles starting from offset phiangle(i)
c**********************************************************************
         dnum = zero
         denom = zero
         dnum = SUM(zSIN(:nphi))
         denom = SUM(rCOS(:nphi))
         IF (denom .ne. zero) THEN
            elongate = dnum/denom
         ELSE
            elongate = 1.E10
         END IF

         delangle = zero
         DO i = 1, nphi
            phiangle(i) = ATAN2(elongate*zCOS(i)-rSIN(i),
     1                          elongate*zSIN(i)+rCOS(i))
            delangle = max(delangle,ABS(phiangle(i)))
            angle(:,i) = angle(:,i) + phiangle(i)
         END DO
         IF (delangle < 0.02_dbl) EXIT
      END DO
CEAL      WRITE (*, 1010) elongate, raxis(1), zaxis(1)
CEAL      WRITE (3, 1010) elongate, raxis(1), zaxis(1)
CEAL      WRITE (*, 1020) ntheta, nphi
CEAL      WRITE (3, 1020) ntheta, nphi
 1010 FORMAT(' Average elongation = ',1pe10.4,/,' Raxis = ',1pe12.4,
     1   ' Zaxis = ',1pe12.4)
 1020 FORMAT(' Number of Theta Points Matched = ',i5,/,
     1   ' Number Phi Planes = ',i5)

      END SUBROUTINE getangle
      SUBROUTINE restart(irst)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE Vname0
      USE Vname1
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER irst
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, SAVE :: nresets = 0
C-----------------------------------------------
!
!       This routine either stores an accepted value of the local solution
!       (irst = 1) or reset the PRESENT solution to a previous value (irst = 2)
!
      SELECT CASE (irst)

      CASE DEFAULT
         xstore(:n2) = xvec(:n2)
         RETURN

      CASE (2)
         xdot(:n2) = zero
         xvec(:n2) = xstore(:n2)
         delt = .95*delt
         irst = 1
         nresets = nresets + 1
         IF (nresets .ge. 100) THEN
            PRINT *, ' Time step reduced 100 times without convergence'
            STOP
         ENDIF
         RETURN
      END SELECT

      END SUBROUTINE restart
      SUBROUTINE order(rval, zval, xaxis, yaxis, inside)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE Vname0
      USE Vname1
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER inside
      REAL(rprec) xaxis, yaxis
      REAL(rprec), DIMENSION(*) :: rval, zval
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, ip1, i1, j, next
      REAL(rprec), DIMENSION(nu) :: tempr, tempz
      REAL(rprec) ::
     1   newdist,olddist,shortest,saver,savez,residue,x,y,dx,dy
C-----------------------------------------------
c**********************************************************************
c       Program ORDER : orders points on a magnetic surface at a
c       fixed toroidal plane and assigns right-handed circulation
c       around flux surface. XAXIS, YAXIS:  Polar-TYPE axis (must lie
c       inside curve to check SIGN of rotation)
c
c Modified 8/13/97, J. Breslau. No longer quits IF the axis is not
c  contained by the points. Instead, RETURNs inside=0 so inguess can
c  try again.
c**********************************************************************
      olddist = 1.E20_dbl
      DO i = 1, ntheta - 1
         ip1 = i + 1
         i1 = i
         shortest = 1.E20_dbl
   15    CONTINUE
         DO j = ip1, ntheta
            IF (i1 > 1) olddist = (rval(i1-1)-rval(j))**2 + (zval(i1-1)-
     1         zval(j))**2
            newdist = (rval(i1)-rval(j))**2 + (zval(i1)-zval(j))**2
            IF (newdist.le.olddist .and. newdist<shortest) THEN
               next = j
               shortest = newdist
            ENDIF
         END DO
c**********************************************************************
c       Swap nearest poINT (next) with current poINT (ip1)
c**********************************************************************
         IF (shortest .ge. 1.E10_dbl) THEN
            SAVEr = rval(i-1)
            rval(i-1) = rval(i)
            rval(i) = SAVEr
            SAVEz = zval(i-1)
            zval(i-1) = zval(i)
            zval(i) = SAVEz
            i1 = i1 - 1
            ip1 = ip1 - 1
            go to 15
         ENDIF
         SAVEr = rval(ip1)
         rval(ip1) = rval(next)
         rval(next) = SAVEr
         SAVEz = zval(ip1)
         zval(ip1) = zval(next)
         zval(next) = SAVEz
      END DO
c**********************************************************************
c       Check that xaxis,yaxis is inside surface and
c       ascertain that the angle rotates counterclockwise
c       using Cauchy theorem in "complex"-plane
c**********************************************************************
      inside = 1
      residue = zero
      DO i = 1, ntheta - 1
         x = 0.5_dbl*(rval(i)+rval(i+1)) - xaxis
         y = 0.5_dbl*(zval(i)+zval(i+1)) - yaxis
         dx = rval(i+1) - rval(i)
         dy = zval(i+1) - zval(i)
         residue = residue + (x*dy - y*dx)/(x**2 + y**2 + 1.E-10_dbl)
      END DO
      x = 0.5_dbl*(rval(1)+rval(ntheta)) - xaxis
      y = 0.5_dbl*(zval(1)+zval(ntheta)) - yaxis
      dx = rval(1) - rval(ntheta)
      dy = zval(1) - zval(ntheta)
      residue = residue + (x*dy - y*dx)/(x**2 + y**2 + 1.E-10_dbl)

      IF (residue < (-0.9_dbl*twopi)) THEN
         DO i = 2, ntheta
            j = ntheta - i + 2
            tempr(i) = rval(j)
            tempz(i) = zval(j)
         END DO
         rval(2:ntheta) = tempr(2:ntheta)
         zval(2:ntheta) = tempz(2:ntheta)
      ELSE IF (ABS(residue) < 0.9_dbl*twopi) THEN
         PRINT *, ' mag. axis not enclosed by bndry; trying again'
         WRITE (3, *) ' mag. axis not enclosed by bndry; trying again'
         inside = 0
c        STOP
      ENDIF

      END SUBROUTINE order
      SUBROUTINE evals_boundary(g2,xp,yp,icount)
      USE Vname0
      USE mapout
      REAL(rprec), DIMENSION(0:mu-1,0:0)  :: rbc,zbs,rbs,zbc 
      REAL*4, DIMENSION(:), ALLOCATABLE  :: rb, zb
      REAL*4, DIMENSION(*) :: xp, yp
      TYPE(RSZS) :: g2
      mpol=mu-1
      ntor=0
      ntheta=g2%nthet
      rbc(0:mu-1,0)=g2%rbc(0:mu-1,0)
      rbs(0:mu-1,0)=g2%rbs(0:mu-1,0)
      zbc(0:mu-1,0)=g2%zbc(0:mu-1,0)
      zbs(0:mu-1,0)=g2%zbs(0:mu-1,0)
      pi=ATAN(1.d0)*4.d0
      pit = 2.d0*pi/(ntheta - 1)
      v=0
      IF( ALLOCATED(rb) ) DEALLOCATE(rb,zb)
      ALLOCATE ( rb(ntheta+1), zb(ntheta+1) )
       icount=0
       DO j=1,ntheta,2
        u=pit*(j-1)
        r=0
        z=0
        n=0
        DO m=0,mpol     
!         DO n=-ntor,ntor
          xxm=m
          arg=xxm*u
          COSa=COS(arg);sina=SIN(arg)
          r=r+rbc(m,n)*cosa+rbs(m,n)*sina
          z=z+zbc(m,n)*cosa+zbs(m,n)*sina
!         ENDdo
        ENDdo
        icount=icount+1
        rb(icount)=r
        zb(icount)=z
       ENDdo
        xp(1:icount)=rb(1:icount)
        yp(1:icount)=zb(1:icount)
       RETURN
      END SUBROUTINE evals_boundary
