      SUBROUTINE funct3d (lscreen, ier_flag)
      USE vmec_main
      USE vacmod, ONLY: bsqvac, raxis_nestor, zaxis_nestor
      USE vmec_params, ONLY: bad_jacobian_flag, zsc
      USE realspace
      USE vforces
      USE vsvd, ONLY: router, rinner, gphifac, grmse
      USE xstuff
      USE timer_sub
      USE precon2d, ONLY: ictrl_prec2d, lHess_exact
      USE vmec_utils, ONLY: cyl2flx
      USE vparams, ONLY: twopi
      USE totzsp_mod
      USE tomnsp_mod
#ifdef _HBANGLE
      USE angle_constraints, ONLY: getrz, getfrho, xtempa
#endif
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(inout) :: ier_flag
      LOGICAL, INTENT(in) :: lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: l0pi, l, lk, ivacskip
      INTEGER :: nvskip0 = 0
      REAL(dp), DIMENSION(mnmax) ::
     1   rmnc, zmns, lmns, rmns, zmnc, lmnc
      REAL(dp), DIMENSION(:), POINTER :: lu, lv
      REAL(dp) :: presf_ns, delr_mse, delt0
      REAL(dp), EXTERNAL :: pmass
C-----------------------------------------------
!
!     POINTER ALIASES
!
!
!     POINTER ALIASES
!
      lu => czmn;  lv => crmn

      CALL second0 (tfunon)

!     CONVERT ODD M TO 1/SQRT(S) INTERNAL REPRESENTATION
!
#ifdef _HBANGLE
!Overwrites rzl_array, but that is OK since gc = rzl_array in CALL, xc preserved
      xtempa = xc
      CALL getrz(xc)
#endif       
      gc(:neqs2) = xc(:neqs2)*scalxc(:neqs2)

!
!     RIGID BODY SHIFT OF RMNCC(JS.GT.1,0,0) BY DELR_MSE= R00-RAXMSE
!
      IF (lrecon) THEN
         delr_mse = xc(neqs2)
         gc(1:ns) = gc(1:ns) + delr_mse
      ENDIF

!
!     INVERSE FOURIER TRANSFORM TO S,THETA,ZETA SPACE
!     R, Z, AND LAMBDA ARRAYS IN FOURIER SPACE
!     FIRST, DO SYMMETRIC [ F(u,v) = F(-u,-v) ] PIECES
!     ON THE RANGE u = 0,pi  and v = 0,2*pi
!
      CALL second0 (tffton)
      CALL totzsps (gc, r1, ru, rv, z1, zu, zv, lu, lv, rcon, zcon)

!
!     ANTI-SYMMETRIC CONTRIBUTIONS TO INVERSE TRANSFORMS
!
      IF (lasym) THEN

         CALL totzspa (gc, armn, brmn, extra3, azmn, bzmn, extra4, 
     1                 blmn, clmn, extra1, extra2)

!        SUM SYMMETRIC, ANTISYMMETRIC PIECES APPROPRIATELY
!        TO GET R, Z, L, (AND RCON, ZCON) ON FULL RANGE OF u (0 to 2*pi)
         IF (ictrl_prec2d .lt. 2)
     1   CALL symrzl (r1, ru, rv, z1, zu, zv, lu, lv, rcon, zcon, armn,
     2   brmn, extra3, azmn, bzmn, extra4, blmn, clmn, extra1, extra2)

      ENDIF

      CALL second0 (tfftoff)
      IF (ictrl_prec2d .LE. 1) timer(tfft) = 
     1                         timer(tfft) + (tfftoff - tffton)
!
!     IN HESSIAN LOOP, ONLY COMPUTE PERTURBATIONS TO R, Z, AND L FOR A SINGLE (m,n,ntype)
!
      IF (ictrl_prec2d .GE. 2) THEN
#ifdef _HBANGLE
         CALL getrz(xcdot)
#endif       
         xcdot(:neqs2) = xcdot(:neqs2)*scalxc(:neqs2)
         CALL totzsps_hess (xcdot, r1, ru, rv, z1, zu, zv, lu, lv, 
     1                      rcon, zcon)
         IF (lasym) THEN
            CALL totzspa_hess (xcdot, armn, brmn, extra3, azmn, bzmn, 
     1                         extra4, blmn, clmn, extra1, extra2)
            CALL symrzl (r1, ru, rv, z1, zu, zv, lu, lv, rcon, zcon,
     2                   armn, brmn, extra3, azmn, bzmn, extra4, blmn,
     3                   clmn, extra1,extra2)
         END IF
      END IF

      l0pi = ns*(1 + nzeta*(ntheta2 - 1))        !u = pi, v = 0, js = ns
      router = r1(ns,0) + r1(ns,1)
      rinner = r1(l0pi,0) + r1(l0pi,1)
      r00 = r1(1,0)
      z00 = z1(1,0)

!
!     COMPUTE CONSTRAINT RCON, ZCON
!
      DO l = 1,nrzt
#ifndef _HBANGLE
         rcon(l,0) = rcon(l,0) + rcon(l,1)*sqrts(l)
         zcon(l,0) = zcon(l,0) + zcon(l,1)*sqrts(l)
#endif
         ru0(l) = ru(l,0) + ru(l,1)*sqrts(l)
         zu0(l) = zu(l,0) + zu(l,1)*sqrts(l)
      END DO

!
!     COMPUTE RCON0, ZCON0 FOR FIXED BOUNDARY BY SCALING EDGE VALUES
!     SCALE BY POWER OF SQRTS, RATHER THAN USE rcon0 = rcon, etc. THIS
!     PREVENTS A DISCONTINUITY WHEN RESTARTING FIXED BOUNDARY WITH NEW RCON0....
!
!     NOTE: IN ORDER TO MAKE INITIAL CONSTRAINT FORCES SAME FOR FREE/FIXED
!     BOUNDARY, WE SET RCON0,ZCON0 THE SAME INITIALLY, BUT TURN THEM OFF
!     SLOWLY IN FREE-BOUNDARY VACUUM LOOP (BELOW)
!
#ifndef _HBANGLE
      IF (iter2.eq.iter1 .and. ivac.le.0) THEN
!!         rcon0(1:nrzt) = rcon(1:nrzt,0)
!!         zcon0(1:nrzt) = zcon(1:nrzt,0)
         DO l = 1, ns
            rcon0(l:nrzt:ns) = rcon(ns:nrzt:ns,0)*sqrts(l:nrzt:ns)**2
            zcon0(l:nrzt:ns) = zcon(ns:nrzt:ns,0)*sqrts(l:nrzt:ns)**2
         END DO
      ENDIF
#endif
!
!     COMPUTE S AND THETA DERIVATIVE OF R AND Z AND JACOBIAN ON HALF-GRID
!
      CALL jacobian
      IF (irst.eq.2 .and. iequi.eq.0) GOTO 100

!
!     COMPUTE COVARIANT COMPONENTS OF B, MAGNETIC AND KINETIC
!     PRESSURE, AND METRIC ELEMENTS ON HALF-GRID
!
      CALL second0 (tbcovon)
      CALL bcovar (lu, lv, xc(1+2*irzloff+mns*(zsc-1)))    !FIX THIS: last arg used elsewhere!

      CALL second0 (tbcovoff)
      IF (ictrl_prec2d .le. 1) timer(tbcov) = timer(tbcov) 
     1                                      + (tbcovoff - tbcovon)

!     COMPUTE VACUUM MAGNETIC PRESSURE AT PLASMA EDGE
!     NOTE: FOR FREE BOUNDARY RUNS, THE VALUE OF RBTOR=R*BTOR
!     AT THE PLASMA EDGE SHOULD BE ADJUSTED TO APPROXIMATELY
!     EQUAL THE VACUUM VALUE. THIS CAN BE DONE BY CHANGING
!     EITHER PHIEDGE OR THE INITIAL CROSS SECTION ACCORDING
!     TO THE SCALING LAW  R*BTOR .EQ. PHIEDGE/(R1 * Z1).

      IF (lfreeb .and. iter2.gt.1 .and. iequi.eq.0) THEN
         IF (ictrl_prec2d.le.1 .and. (fsqr + fsqz).le.1.e-3_dp) 
     1      ivac = ivac+1   !decreased from e-1 to e-3 - sph12/04
         IF (nvskip0 .eq. 0) nvskip0 = MAX(1, nvacskip)
         IF (ivac .ge. 0) THEN
!           IF INITIALLY ON, MUST TURN OFF rcon0, zcon0 SLOWLY
            IF (ictrl_prec2d .eq. 2) THEN
               rcon0 = 0;  zcon0 = 0
            ELSE IF (ictrl_prec2d .eq. 0) THEN
               rcon0 = 0.9_dp*rcon0;  zcon0 = 0.9_dp*zcon0
            END IF
            CALL second0 (tvacon)
            ivacskip = MOD(iter2-iter1,nvacskip)
            IF (ivac .le. 2) ivacskip = 0

!           EXTEND NVACSKIP AS EQUILIBRIUM CONVERGES
            IF (ivacskip .eq. 0) THEN
               nvacskip = one/MAX(1.e-1_dp, 1.e11_dp*(fsqr+fsqz))
               nvacskip = MAX(nvacskip, nvskip0)
            END IF

!
!           NORMALLY, WHEN COMPUTING THE HESSIAN, IT IS SUFFICIENT TO
!           COMPUTE THE VARIATIONS IN THE "EXACT" SOLUTION, NOT THE ENTIRE
!           FIELD PERIOD SUM. THUS, FOR ictrl_prec2d >= 2, SET ivacskip = 1
!           FOR ictrl_prec2d = 1 (RUN WITH PRECONDITIONER APPLIED), MUST
!           COMPUTE EXACT VACUUM RESPONSE NOW.
!
!           THE EXCEPTION TO THIS IS IF WE ARE TESTING THE HESSIAN (lHess_exact=T), 
!           THEN MUST USE FULL VACUUM CALCULATION TO COMPUTE IT (ivacskip=0)
!
            IF (lHess_exact) THEN
               IF (ictrl_prec2d .ge. 1) ivacskip = 0       !Accurate Hessian
            ELSE IF (ictrl_prec2d .ne. 0) THEN
!               IF (ictrl_prec2d .ge. 2) ivacskip = 1
!               IF (ictrl_prec2d .eq. 1) ivacskip = 0
               IF (ictrl_prec2d .le. 2) ivacskip = 0
               IF (ictrl_prec2d .eq. 3) ivacskip = 1       !Fast vacuum calculation used to compute Hessian
            END IF

!           DO NOT UPDATE THIS WHEN USING PRECONDITIONER: BREAKS TRI-DIAGONAL STRUCTURE
            IF (ictrl_prec2d .eq. 0) THEN
               raxis_nestor(1:nzeta) = r1(1:ns*nzeta:ns,0)
               zaxis_nestor(1:nzeta) = z1(1:ns*nzeta:ns,0)
            END IF

!           NOTE: gc contains correct edge values of r,z,l arrays
!                 convert_sym, convert_asym have been applied to m=1 modes
            CALL convert (rmnc, zmns, lmns, rmns, zmnc, lmnc, gc, ns)
            CALL vacuum (rmnc, rmns, zmns, zmnc, xm, xn, 
     1                   ctor, rbtor, wint, ns, ivacskip, ivac, mnmax,
     2                   ier_flag, lscreen)
            IF (ier_flag .ne. 0) GOTO 100
!
!           RESET FIRST TIME FOR SOFT START
!
            IF (ivac .eq. 1) THEN
               irst = 2;  delt0 = delt
#ifdef _HBANGLE
               gc = xc; xc = xtempa
#endif
               CALL restart_iter(delt0)
#ifdef _HBANGLE
               xc = gc
#endif
               irst = 1
            END IF

!
!           IN CASE PRESSURE IS NOT ZERO AT EXTRAPOLATED EDGE...
!           UNCOMMENT ALL "RPRES" COMMENTS HERE AND IN BCOVAR, FORCES ROUTINES
!           IF NON-VARIATIONAL FORCES ARE DESIRED
!
!           presf_ns = 1.5_dp*pres(ns) - 0.5_dp*pres(ns1)  
!           MUST NOT BREAK TRI-DIAGONAL RADIAL COUPLING: OFFENDS PRECONDITIONER!
            presf_ns = pmass(hs*(ns-1.5_dp))
            IF (presf_ns .ne. zero) 
     1         presf_ns = (pmass(1._dp)/presf_ns) * pres(ns)

            lk = 0
            gcon(:nrzt) = r1(:nrzt,0)+sqrts(:nrzt)*r1(:nrzt,1)
            gcon(1+nrzt) = 0
            DO l = ns, nrzt, ns
               lk = lk + 1
               bsqsav(lk,3) = 1.5_dp*bzmn_o(l) - 0.5_dp*bzmn_o(l-1)
#ifdef _ANIMEC
               gcon(l)     = bsqvac(lk) + pperp_ns(lk)
#else
               gcon(l)     = bsqvac(lk) + presf_ns
#endif 

               rbsq(lk) = gcon(l)*(r1(l,0) + r1(l,1))*ohs
               dbsq(lk) = ABS(gcon(l)-bsqsav(lk,3))  
            END DO
!           
!           COMPUTE m=0,n=0 EDGE "pedestals"
!
!!            alphaR = hs*hs*ard(ns,1)
!!            IF (alphaR .ne. zero) alphaR = 
!!     1         hs*SUM(wint(ns:nrzt:ns)*zu0(ns:nrzt:ns)*rbsq)/alphaR

!!            PRINT *,' alphaR/r1(ns) = ', alphaR/gcon(ns)

            IF (ivac .eq. 1) THEN
               bsqsav(:nznt,1) = bzmn_o(ns:nrzt:ns)
               bsqsav(:nznt,2) = bsqvac(:nznt)
            ENDIF
            CALL second0 (tvacoff)
            IF (ictrl_prec2d .le. 1)
     1         timer(tvac) = timer(tvac) + (tvacoff - tvacon)
         ENDIF
      ENDIF

!
!     COMPUTE CONSTRAINT FORCE
!
#ifdef _HBANGLE
      gcon = 0
      IF (iequi .EQ. 1) GOTO 100
#else
      IF (iequi .NE. 1) THEN
         extra1(:nrzt,0) = (rcon(:nrzt,0) - rcon0(:nrzt))*ru0(:nrzt)
     1                   + (zcon(:nrzt,0) - zcon0(:nrzt))*zu0(:nrzt)
         CALL alias (gcon, extra1(:,0), gc, gc(1+mns), gc(1+2*mns), 
     1               extra1(:,1))
      ELSE 
         IF (lrecon) xc(:ns) = xc(:ns) + delr_mse
         GOTO 100
      END IF
#endif
!
!     COMPUTE MHD FORCES ON INTEGER-MESH
!
      CALL second0 (tforon)
      CALL forces

!
!     FFT TEST
!
!      xc = xc*scalxc
!      CALL tomnsps_t (gc, xc, r1, ru, rv, z1, zu, zv)
!      IF (lasym) CALL tomnspa_t (gc, xc, r1, ru, rv, z1, zu, zv)
!      STOP

!
!     SYMMETRIZE FORCES (in u-v space)
!
      IF (lasym) CALL symforce (armn, brmn, crmn, azmn, bzmn,
     1     czmn, blmn, clmn, rcon, zcon, r1, ru, rv, z1, zu, zv,
     2     extra3, extra4, extra1, extra2)

      CALL second0 (tforoff)
      IF (ictrl_prec2d .le. 1) timer(tfor) = timer(tfor) 
     1                                     + (tforoff - tforon)
!
!     FOURIER-TRANSFORM MHD FORCES TO (M,N)-SPACE
!
      CALL second0 (tffton)
      CALL tomnsps (gc, armn, brmn, crmn, azmn, bzmn, czmn, 
     1              blmn, clmn, rcon, zcon)

      IF (lasym) CALL tomnspa (gc, r1, ru, rv, z1, zu, zv,
     1                         extra3, extra4, extra1, extra2)

      CALL second0 (tfftoff)
      IF (ictrl_prec2d .le. 1) timer(tffi) = timer(tffi) 
     1                                     + (tfftoff - tffton)

!================================================================
!
!     COMPUTE FORCE RESIDUALS (RAW AND PRECONDITIONED)
!
!================================================================
      CALL second0 (treson)

      gc = gc * scalxc    !!IS THIS CORRECT: SPH010214?
#ifdef _HBANGLE
      CALL getfrho(gc)
#endif
      CALL residue (gc, gc(1+irzloff), gc(1+2*irzloff))
!     Force new initial axis guess IF ALLOWED (l_moveaxis=T)
      IF (lmove_axis .and. iter2.eq.1 .and. (fsqr+fsqz+fsql).gt.1.E2_dp)
     1    irst = 4

      CALL second0 (tresoff)
      IF (ictrl_prec2d .le. 1) timer(tres) = timer(tres) 
     1                                     + (tresoff - treson)

      gc(neqs1) = gphifac
      IF (iopt_raxis .gt. 0) gc(neqs2) = grmse

 100  CONTINUE

#ifdef _HBANGLE
      xc = xtempa
#endif
      CALL second0 (tfunoff)
      IF (ictrl_prec2d .le. 1) timer(tfun) = timer(tfun) 
     1                                     + (tfunoff - tfunon)

      END SUBROUTINE funct3d
