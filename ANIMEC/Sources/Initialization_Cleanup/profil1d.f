      SUBROUTINE profil1d(xc, xcdot, lreset)
      USE vmec_main
      USE vmec_params, ONLY: signgs, lamscale, rcc, pdamp
      USE vmec_input, ONLY: lRFP
      USE vsvd, torflux_edge => torflux
      USE vspline
      USE realspace, ONLY: shalf, sqrts
      USE init_geometry, ONLY: lflip
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(neqs2), INTENT(out) :: xc, xcdot
      LOGICAL, INTENT(in) :: lreset
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: c1p5 = 1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i
      REAL(rprec) :: Itor, si, pedge, vpnorm, polflux_edge
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      REAL(rprec), EXTERNAL :: pcurr, pmass, piota, torflux,
     1    torflux_deriv, polflux, polflux_deriv
#ifdef _ANIMEC
     2  , photp, ptrat
#endif
C-----------------------------------------------
!
!                INDEX OF LOCAL VARIABLES
!
!        ai       array of coefficients in phi-series for iota (ncurr=0)
!        ac       array of coefficients in phi-series for the quantity d(Icurv)/ds = toroidal
!                 current density * Vprime, so Icurv(s) = Itor(s) (used for ncurr=1)
!        am       array of coefficients in phi-series for mass (NWT/m**2)
!        iotas    rotational transform , on half radial mesh
!        Icurv    (-)toroidal current inside flux surface (vanishes like s)
!        mass     mass profile on half-grid
!        phiedge  value of real toroidal flux at plasma edge (s=1)
!        phips    toroidal flux (same as phip), one-dimensional array
!        chips    poloidal flux (same as chip), one-dimensional array
!        presf    pressure profile on full-grid, mass/phip**gamma
!        spres_ped value of s beyond which pressure profile is flat (pedestal)
!
!
!     COMPUTE PHIP, IOTA PROFILES ON FULL-GRID
!     COMPUTE MASS PROFILE ON HALF-GRID
!     BY READING INPUT COEFFICIENTS. PRESSURE CONVERTED TO
!     INTERNAL UNITS BY MULTIPLICATION BY mu0 = 4*pi*10**-7
!
      IF (ncurr.EQ.1 .AND. lRFP)STOP 'ncurr=1 inconsistent with lRFP=T!'

      torflux_edge = signgs * phifac * phiedge / twopi
      si = torflux(one)
      IF (si .ne. zero) torflux_edge = torflux_edge/si
      polflux_edge = torflux_edge
      si = polflux(one)
      IF (si .ne. zero) polflux_edge = polflux_edge/si
      r00 = rmn_bdy(0,0,rcc)

      phips(1) = 0
      chips(1) = 0
      icurv(1) = 0

      DO i = 2,ns
         si = hs*(i-c1p5)
!         si = torflux(si)          !SPH060409: Fcn of s, NOT flux
         phips(i) = torflux_edge * torflux_deriv(si)
         chips(i) = torflux_edge * polflux_deriv(si)
         iotas(i) = piota(si)
         icurv(i) = pcurr(si)
      END DO

!
!     Compute lamscale factor for "normalizing" lambda (needed for scaling hessian)
!
      lamscale = SQRT(hs*SUM(phips(2:ns)**2))
      phips(ns+1) = 2*phips(ns)-phips(ns-1)
      IF (lamscale .EQ. 0) STOP 'PHIP == 0: ERROR!'

      IF (lflip) THEN
         iotas = -iotas
         chips = -chips
      END IF

      DO i = 1,ns
         si = hs*(i-1)
!         si = torflux(si)
         iotaf(i) = piota(si)
         phipf(i) = torflux_edge * torflux_deriv(si)
         chipf(i) = torflux_edge * polflux_deriv(si)
      ENDDO
!
!     SCALE CURRENT TO MATCH INPUT EDGE VALUE, CURTOR
!     FACTOR OF SIGNGS NEEDED HERE, SINCE MATCH IS MADE TO LINE
!     INTEGRAL OF BSUBU (IN GETIOTA) ~ SIGNGS * CURTOR
!
      pedge = pcurr(one)
      Itor = 0
      IF (ABS(pedge) .gt. ABS(EPSILON(pedge)*curtor))
     1   Itor = signgs*currv/(twopi*pedge)
      icurv(2:ns) = Itor*icurv(2:ns)

!
!     POSSIBLE PRESSURE PEDESTAL FOR S >= SPRES_PED
!
      spres_ped = ABS(spres_ped)
      IF (.not.lrecon) THEN
        DO i = 2,ns
          si = hs*(i - c1p5)

!         NORMALIZE mass so dV/dPHI (or dV/dPSI) in pressure to mass relation
!         See line 195 of bcovar: pres(2:ns) = mass(2:ns)/vp(2:ns)**gamma

          IF (lRFP) THEN
             vpnorm = polflux_edge * polflux_deriv(si) 
          ELSE
             vpnorm = torflux_edge * torflux_deriv(si)
          END IF
          IF (si .gt. spres_ped) THEN
             pedge = pmass(spres_ped)
          ELSE
!             si = torflux(si)
             pedge = pmass(si)
          END IF
          mass(i) = pedge*(ABS(vpnorm)*r00)**gamma
#ifdef _ANIMEC
!         ANISOTROPIC PRESSURE, Tper/T|| RATIOS
          phot(i) = photp(si)
          tpotb(i)= ptrat(si)
#endif
        END DO

      ELSE
        iotas(:ns) = 0
        iotaf(:ns) = 0
        mass (:ns) = 0
        presf(:ns) = 0
      END IF

      pres(:ns+1) = 0
      xcdot(:neqs2) = 0

#ifdef _ANIMEC
      medge  = pmass (one) * (ABS(phips(ns))*r00)**gamma
      phedg  = photp (one)
#endif
      DO i = 1, ns
         si = hs*ABS(i-1.5_dp)
!         si = torflux(si)                      !SPH060409: shalf = sqrt(s), NOT sqrt(phi(s))!
         shalf(i:nrzt:ns) = SQRT(si)
         si = hs*(i-1)
         sqrts(i:nrzt:ns) = SQRT(si)
         bdamp(i) = 2*pdamp*(1-si)
      END DO

      sqrts(ns:nrzt:ns) = 1     !!Avoid round-off
      shalf(nrzt+1) = 1
      sqrts(nrzt+1) = 1

      DO i = 2,ns
         sm(i) = shalf(i)/sqrts(i)
         sp(i) = shalf(i+1)/sqrts(i)
      ENDDO
      sm(1) = 0
      sp(0) = 0
      sp(1) = sm(2)

      IF (lreset) THEN
        xc(:neqs1) = 0
        IF (lrecon) iresidue = 0
      END IF

      IF (lrecon) THEN
        IF (iresidue .gt. 1) iresidue = 1
!
!       COMPUTE INDEX ARRAY FOR FINDPHI ROUTINE
!
        DO i = 1,ns
          indexr(i)    = ns + 1 - i                    !FINDPHI
          indexr(i+ns) = i + 1                         !FINDPHI
        ENDDO
        indexr(2*ns) = ns
      END IF

      END SUBROUTINE profil1d
