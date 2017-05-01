      SUBROUTINE funct3d_par (lscreen, ier_flag)
      USE vmec_main
      USE vacmod, ONLY: bsqvac, raxis_nestor, zaxis_nestor, nuv, nuv3
      USE vmec_params, ONLY: bad_jacobian_flag, zsc, lamscale
      USE vmec_params, ONLY: ntmax, misc_error_flag
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
      USE parallel_include_module
#if defined (SKS)
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
      REAL(rprec), DIMENSION(mnmax) ::
     1   rmnc, zmns, lmns, rmns, zmnc, lmnc
      REAL(rprec), DIMENSION(:,:,:), POINTER :: lu, lv
      REAL(rprec) :: presf_ns, delr_mse, delt0
      REAL(rprec), EXTERNAL :: pmass

      REAL(rprec) :: skston, skstoff, tvaconly_off, tvaconly_on
      INTEGER :: i, j, k, nsmin, nsmax, m
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: bcastbuf
      INTEGER, DIMENSION(4) :: bbuf
C-----------------------------------------------
!
!     POINTER ALIASES
!

      nfunct3d = nfunct3d + 1
      lu => pczmn;  lv => pcrmn

      CALL second0 (tfunon)

!     CONVERT ODD M TO 1/SQRT(S) INTERNAL REPRESENTATION
      ACTIVE1: IF (lactive) THEN 
      DO l=1, 3*ntmax
        DO i=t2lglob, t2rglob
          DO k=0, mpol1
            DO j=0, ntor
              lk = j + (ntor+1)*(k + (mpol1+1)*
     1             ((i-1)+ns*(l-1)))+1
              pgc(lk) = pxc(lk)*pscalxc(lk)
            END DO
          END DO
        END DO
      END DO

!
!     RIGID BODY SHIFT OF RMNCC(JS.GT.1,0,0) BY DELR_MSE= R00-RAXMSE
!

!      IF (lrecon) THEN
!         delr_mse = pxc(neqs2)
!         pgc(1:ns) = pgc(1:ns) + delr_mse
!      ENDIF

!     INVERSE FOURIER TRANSFORM TO S,THETA,ZETA SPACE
!     R, Z, AND LAMBDA ARRAYS IN FOURIER SPACE
!     FIRST, DO SYMMETRIC [ F(u,v) = F(-u,-v) ] PIECES
!     ON THE RANGE u = 0,pi  and v = 0,2*pi
!
      CALL second0 (tffton)
      CALL second0(skston)
      CALL totzsps_par (pgc, pr1, pru, prv, pz1, pzu, pzv, lu, lv, 
     1                  prcon, pzcon)
      CALL second0(skstoff)
      totzsps_time=totzsps_time+(skstoff-skston)

      CALL MPI_BCast(lerror_sam,1,MPI_LOGICAL,0,RUNVMEC_COMM_WORLD,
     1               MPI_ERR)
      IF (lerror_sam) THEN
         ier_flag = bad_jacobian_flag
         irst = 2
         CALL MPI_BCast(ier_flag,1,MPI_INTEGER,0,
     1                  RUNVMEC_COMM_WORLD,MPI_ERR)
         RETURN
      END IF


!     ANTI-SYMMETRIC CONTRIBUTIONS TO INVERSE TRANSFORMS
!
      IF (lasym) THEN
         CALL totzspa_par (pgc, parmn, pbrmn, pextra3, pazmn, pbzmn,
     1                     pextra4, pblmn, pclmn, pextra1, pextra2)
         CALL second0(skstoff)

!        SUM SYMMETRIC, ANTISYMMETRIC PIECES APPROPRIATELY
!        TO GET R, Z, L, (AND RCON, ZCON) ON FULL RANGE OF u (0 to 2*pi)
         IF (ictrl_prec2d .lt. 2)
     1   CALL symrzl_par (pr1, pru, prv, pz1, pzu, pzv, lu, lv, prcon,
     2   pzcon, parmn, pbrmn, pextra3, pazmn, pbzmn, pextra4, pblmn,
     3   pclmn, pextra1, pextra2)
      ENDIF

      CALL second0 (tfftoff)
      IF (ictrl_prec2d .le. 1) timer(tfft) = 
     1                         timer(tfft) + (tfftoff - tffton)

!
!     IN HESSIAN LOOP, ONLY COMPUTE PERTURBATIONS TO R, Z, AND L FOR A SINGLE (m,n,ntype)
!
      ! SKS: The next conditional is skipped for now
      IF (ictrl_prec2d .GE. 2) THEN

         STOP "This conditional has not been parallelized
     1    yet...stopping in funct3d @ 1!"

         xcdot(:neqs2) = xcdot(:neqs2)*pscalxc(:neqs2)
         CALL totzsps_hess (xcdot, r1, ru, rv, z1, zu, zv, lu, lv, 
     1                      rcon, zcon)
         IF (lasym) THEN
            CALL totzspa_hess (xcdot, armn, brmn, extra3, azmn, bzmn, 
     1                         extra4, blmn, clmn, extra1, extra2)
            CALL symrzl_par (r1, ru, rv, z1, zu, zv, lu, lv, rcon, zcon,
     2                   armn, brmn, extra3, azmn, bzmn, extra4, blmn,
     3                   clmn, extra1,extra2)
         END IF
      END IF

!      l0pi = ns*(1 + nzeta*(ntheta2 - 1))        !u = pi, v = 0, js = ns
!      router = r1(ns,0) + r1(ns,1)
!      rinner = r1(l0pi,0) + r1(l0pi,1)
      r00 = pr1(1,1,0)
      z00 = pz1(1,1,0)


!
!     COMPUTE CONSTRAINT RCON, ZCON
!
      nsmin=tlglob; nsmax=t1rglob
      prcon(:,nsmin:nsmax,0) = prcon(:,nsmin:nsmax,0) + 
     1  prcon(:,nsmin:nsmax,1)*psqrts(:,nsmin:nsmax)
      pzcon(:,nsmin:nsmax,0) = pzcon(:,nsmin:nsmax,0) + 
     1  pzcon(:,nsmin:nsmax,1)*psqrts(:,nsmin:nsmax)
      pru0(:,nsmin:nsmax) = pru(:,nsmin:nsmax,0) + 
     1  pru(:,nsmin:nsmax,1)*psqrts(:,nsmin:nsmax)
      pzu0(:,nsmin:nsmax) = pzu(:,nsmin:nsmax,0) + 
     1  pzu(:,nsmin:nsmax,1)*psqrts(:,nsmin:nsmax)

!     COMPUTE RCON0, ZCON0 FOR FIXED BOUNDARY BY SCALING EDGE VALUES
!     SCALE BY POWER OF SQRTS, RATHER THAN USE rcon0 = rcon, etc. THIS
!     PREVENTS A DISCONTINUITY WHEN RESTARTING FIXED BOUNDARY WITH NEW RCON0....
!
!     NOTE: IN ORDER TO MAKE INITIAL CONSTRAINT FORCES SAME FOR FREE/FIXED
!     BOUNDARY, WE SET RCON0,ZCON0 THE SAME INITIALLY, BUT TURN THEM OFF
!     SLOWLY IN FREE-BOUNDARY VACUUM LOOP (BELOW)
!
      nsmin=tlglob; nsmax=t1rglob
      IF (iter2.eq.iter1 .and. ivac.le.0) THEN
        ALLOCATE(bcastbuf(2*nznt))
        bcastbuf(1:nznt)=prcon(:,ns,0)
        bcastbuf(nznt+1:2*nznt)=pzcon(:,ns,0)
        CALL second0(skston)
        CALL MPI_Bcast(bcastbuf,SIZE(bcastbuf),MPI_REAL8,nranks-1,
     1                  NS_COMM,MPI_ERR)
        CALL second0(skstoff)
        broadcast_time = broadcast_time + (skstoff -skston)
        prcon(:,ns,0)=bcastbuf(1:nznt)
        pzcon(:,ns,0)=bcastbuf(nznt+1:2*nznt)
        DEALLOCATE(bcastbuf)
        DO l = nsmin, nsmax
          prcon0(:,l) = prcon(:,ns,0)*psqrts(:,l)**2
          pzcon0(:,l) = pzcon(:,ns,0)*psqrts(:,l)**2
        END DO
      ENDIF

!
!     COMPUTE S AND THETA DERIVATIVE OF R AND Z AND JACOBIAN ON HALF-GRID
!

      CALL second0(skston)
      CALL jacobian_par
      CALL second0(skstoff)
      jacobian_time=jacobian_time+(skstoff-skston)

!     COMPUTE COVARIANT COMPONENTS OF B, MAGNETIC AND KINETIC
!     PRESSURE, AND METRIC ELEMENTS ON HALF-GRID

      CALL second0 (tbcovon)
      CALL bcovar_par (lu, lv, pxc, pxc(1+2*irzloff+mns*(zsc-1)))
      CALL second0 (tbcovoff)
      bcovar_time=bcovar_time+(tbcovoff - tbcovon)

      END IF ACTIVE1

      bbuf(1)=irst; bbuf(2)=iequi; bbuf(3)=ivac; bbuf(4)=iter2
      CALL MPI_BCast(bbuf,4,MPI_INTEGER,0,RUNVMEC_COMM_WORLD,MPI_ERR)
      irst=bbuf(1); iequi=bbuf(2); ivac=bbuf(3); iter2=bbuf(4)
      CALL MPI_BCast(lfreeb,1,MPI_LOGICAL,0,RUNVMEC_COMM_WORLD,MPI_ERR)

      IF (irst.EQ.2 .AND. iequi.EQ.0) GOTO 100

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
         IVAC0: IF (ivac .ge. 0) THEN

!           IF INITIALLY ON, MUST TURN OFF rcon0, zcon0 SLOWLY
            IF (ictrl_prec2d .eq. 2) THEN
               prcon0 = 0;  pzcon0 = 0
            ELSE IF (ictrl_prec2d .eq. 0) THEN
               prcon0 = 0.9_dp*prcon0;  pzcon0 = 0.9_dp*pzcon0
            END IF
            CALL second0 (tvacon)
            ivacskip = MOD(iter2-iter1,nvacskip)
            IF (ivac .le. 2) ivacskip = 0


!           EXTEND NVACSKIP AS EQUILIBRIUM CONVERGES
            IF (ivacskip .EQ. 0) THEN
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
!             IF (ictrl_prec2d .ge. 2) ivacskip = 1
!             IF (ictrl_prec2d .eq. 1) ivacskip = 0
              IF (ictrl_prec2d .le. 2) ivacskip = 0
              IF (ictrl_prec2d .eq. 3) ivacskip = 1       !Fast vacuum calculation used to compute Hessian
           END IF

!           DO NOT UPDATE THIS WHEN USING PRECONDITIONER: BREAKS TRI-DIAGONAL STRUCTURE
           IF (ictrl_prec2d .eq. 0) THEN
              CALL second0(skston)
              raxis_nestor(1:nzeta) = pr1(1:nzeta,1,0)
              zaxis_nestor(1:nzeta) = pz1(1:nzeta,1,0)

              ALLOCATE (bcastbuf(2*nzeta))
              bcastbuf(1:nzeta) = raxis_nestor(1:nzeta)
              bcastbuf(nzeta+1:2*nzeta) = zaxis_nestor(1:nzeta)
              CALL second0(skston)
              CALL MPI_Bcast(bcastbuf,SIZE(bcastbuf),MPI_REAL8,0,
     1                       RUNVMEC_COMM_WORLD,MPI_ERR)
              CALL second0(skstoff)
              broadcast_time = broadcast_time + (skstoff - skston)
              raxis_nestor(1:nzeta) = bcastbuf(1:nzeta)
              zaxis_nestor(1:nzeta) = bcastbuf(nzeta+1:2*nzeta)
              DEALLOCATE (bcastbuf)
              CALL second0(skstoff)
            END IF

!           NOTE: gc contains correct edge values of r,z,l arrays
!                 convert_sym, convert_asym have been applied to m=1 modes

           CALL convert_par(rmnc,zmns,lmns,rmns,zmnc,lmnc,pgc,ns)

           CALL second0(skston)
           ALLOCATE (bcastbuf(2))
           bcastbuf(1)=rbtor
           bcastbuf(2)=ctor
           IF (lactive) THEN
             CALL second0(skston)
             CALL MPI_Bcast(bcastbuf,SIZE(bcastbuf),MPI_REAL8,
     1                        nranks-1,NS_COMM,MPI_ERR)
             CALL second0(skstoff)
             broadcast_time = broadcast_time + (skstoff -skston)
           END IF

           CALL second0(skston)
           IF (vlactive) THEN
             CALL MPI_Bcast(bcastbuf,SIZE(bcastbuf),MPI_REAL8,0,
     1                     VAC_COMM,MPI_ERR)
           END IF
           CALL second0(skstoff)
           broadcast_time = broadcast_time + (skstoff -skston)
           rbtor=bcastbuf(1)
           ctor=bcastbuf(2)
           DEALLOCATE (bcastbuf)
           CALL second0(skstoff)
           broadcast_time = broadcast_time + (skstoff - skston)

           IF (vlactive) THEN
             CALL second0(tvaconly_on)
             CALL vacuum_par (rmnc, rmns, zmns, zmnc, xm, xn, 
     1                       ctor, rbtor, pwint_ns, ns, ivacskip, ivac, 
     2                       mnmax,ier_flag, lscreen)
             CALL second0(tvaconly_off)

           END IF
           IF (vnranks.LT.nranks) THEN
             CALL MPI_Bcast(bsqvac,SIZE(bsqvac),MPI_REAL8,0,
     1                   NS_COMM,MPI_ERR)
           END IF

           IF (ier_flag .ne. 0) RETURN
!
!           RESET FIRST TIME FOR SOFT START
!
           IF (ivac .eq. 1) THEN
             irst = 2;  delt0 = delt
             CALL restart_iter(delt0)
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

           DO l = 1, nznt
             bsqsav(l,3) = 1.5_dp*pbzmn_o(l,ns) -
     1                     0.5_dp*pbzmn_o(l,ns-1)
             pgcon(l,ns) = bsqvac(l) + presf_ns
             rbsq(l) = pgcon(l,ns)*(pr1(l,ns,0) + pr1(l,ns,1))*ohs
             dbsq(l) = ABS(pgcon(l,ns)-bsqsav(l,3))  
           END DO

           IF (ivac .EQ. 1) THEN
             !IF (lActive) THEN
             IF (vlactive) THEN
             bsqsav(:nznt,1) = pbzmn_o(:,ns)
             bsqsav(:nznt,2) = bsqvac(:nznt)
             CALL MPI_Bcast(bsqsav(:,1),nznt,MPI_REAL8,
     1                      nranks-1,NS_COMM,MPI_ERR)
             END IF
           ELSE
             CALL MPI_Bcast(bsqsav(:,1),nznt,MPI_REAL8,
     1                      0,NS_COMM,MPI_ERR)
           END IF

           CALL second0 (tvacoff)
           IF (ictrl_prec2d .le. 1) THEN
             timer(tvac) = timer(tvac) + (tvacoff - tvacon)
             timer_vac(tbcast_vac) = timer_vac(tbcast_vac) +
     1         (tvacoff - tvacon) - (tvaconly_off-tvaconly_on)
           END IF

         ENDIF IVAC0
       ENDIF

!
!     COMPUTE CONSTRAINT FORCE
!
       ACTIVE2: IF (lactive) THEN
         IF (iequi .NE. 1) THEN
           nsmin=tlglob; nsmax=t1rglob
           pextra1(:,nsmin:nsmax,0) = (prcon(:,nsmin:nsmax,0) - 
     1          prcon0(:,nsmin:nsmax))*pru0(:,nsmin:nsmax) +
     2          (pzcon(:,nsmin:nsmax,0) - 
     3          pzcon0(:,nsmin:nsmax))*pzu0(:,nsmin:nsmax)

           CALL second0(skston)
           CALL alias_par (pgcon, pextra1(:,:,0), pgc, pgc(1+mns), 
     1                  pgc(1+2*mns), pextra1(:,:,1))
        CALL second0(skstoff)
        alias_time = alias_time + (skstoff-skston)         
      ELSE 
        IF (lrecon) pxc(:ns) = pxc(:ns) + delr_mse
        GOTO 100
      END IF

!
!     COMPUTE MHD FORCES ON INTEGER-MESH
!
      CALL second0 (tforon)
      skston = tforon
      CALL forces_par
      CALL second0 (skstoff)
      forces_time = forces_time + (skstoff-skston)

!     SYMMETRIZE FORCES (in u-v space)
!

      IF (lasym) CALL symforce_par (parmn, pbrmn, pcrmn, pazmn, pbzmn,
     1     pczmn, pblmn, pclmn, prcon, pzcon, pr1, pru, prv, pz1,
     2     pzu, pzv, pextra3, pextra4, pextra1, pextra2)


      CALL second0 (tforoff)
      IF (ictrl_prec2d .le. 1) timer(tfor) = timer(tfor) 
     1                                     + (tforoff - tforon)
      CALL second0 (tffton)
     
!
!     FOURIER-TRANSFORM MHD FORCES TO (M,N)-SPACE
!

      CALL second0(skston)
         CALL tomnsps_par (pgc, parmn, pbrmn, pcrmn, pazmn, pbzmn,
     1                     pczmn, pblmn, pclmn, prcon, pzcon)
      CALL second0(skstoff)
      tomnsps_time = tomnsps_time + (skstoff-skston)

      IF (lasym) CALL tomnspa_par (pgc, pr1, pru, prv, pz1, pzu, pzv,
     1                             pextra3, pextra4, pextra1, pextra2)

      CALL second0 (tfftoff)
      IF (ictrl_prec2d .le. 1) timer(tffi) = timer(tffi) 
     1                                     + (tfftoff - tffton)

!================================================================
!
!     COMPUTE FORCE RESIDUALS (RAW AND PRECONDITIONED)
!
!================================================================
      CALL second0 (treson)

      DO l=1, 3*ntmax
        DO i=tlglob, t1rglob
          DO k=0, mpol1
            DO j=0, ntor
              lk = j + (ntor+1)*(k +(mpol1+1)*((i-1)+ns*(l-1)))+1
              pgc(lk) = pgc(lk) * pscalxc(lk)
            END DO
          END DO
        END DO
      END DO

      CALL second0 (skston)
      CALL residue_par(pgc, pgc(1+irzloff), pgc(1+2*irzloff))
      CALL second0 (skstoff)
      residue_time = residue_time + (skstoff-skston)

      END IF ACTIVE2

!NEED THIS ON ALL PROCESSORS IN GROUP (NOT JUST ACTIVE ONES) FOR STOPPING CRITERION IN EVOLVE
      IF (gnranks .GT. nranks) THEN
        CALL second0(skston)
        ALLOCATE(bcastbuf(6))
        bcastbuf(1) = fsqr; bcastbuf(2) = fsqr1
        bcastbuf(3) = fsqz; bcastbuf(4) = fsqz1
        bcastbuf(5) = fsql; bcastbuf(6) = fsql1
        CALL second0(skston)
        CALL MPI_Bcast(bcastbuf,SIZE(bcastbuf),MPI_REAL8,0,
     1                 RUNVMEC_COMM_WORLD,MPI_ERR)
        CALL second0(skstoff)
        broadcast_time = broadcast_time + (skstoff -skston)
        fsqr = bcastbuf(1); fsqr1 = bcastbuf(2)
        fsqz = bcastbuf(3); fsqz1 = bcastbuf(4)
        fsql = bcastbuf(5); fsql1 = bcastbuf(6)
        DEALLOCATE(bcastbuf)
        CALL second0(skstoff)
        broadcast_time = broadcast_time + (skstoff -skston)
      ENDIF


!     Force new initial axis guess IF ALLOWED (l_moveaxis=T)
      IF (lmove_axis .and. iter2.eq.1 .and. (fsqr+fsqz+fsql).gt.1.E2_dp)
     1    irst = 4

      CALL second0 (tresoff)
      IF (ictrl_prec2d .le. 1) timer(tres) = timer(tres) 
     1                                     + (tresoff - treson)

      pgc(neqs1) = gphifac
      IF (iopt_raxis .gt. 0) pgc(neqs2) = grmse

 100  CONTINUE

      CALL second0 (tfunoff)
      IF (ictrl_prec2d .le. 1) timer(tfun) = timer(tfun) 
     1                                     + (tfunoff - tfunon)
#endif

      END SUBROUTINE funct3d_par

      SUBROUTINE funct3d (lscreen, ier_flag)
      USE vmec_main
      USE vacmod, ONLY: bsqvac, raxis_nestor, zaxis_nestor
      USE vmec_params, ONLY: bad_jacobian_flag, zsc, lamscale
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
      USE parallel_include_module
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
      REAL(rprec), DIMENSION(mnmax) ::
     1   rmnc, zmns, lmns, rmns, zmnc, lmnc
      REAL(rprec), DIMENSION(:), POINTER :: lu, lv
      REAL(rprec) :: presf_ns, delr_mse, delt0
      REAL(rprec), EXTERNAL :: pmass
#if defined(MPI_OPT)      
      REAL(rprec) :: skston, skstoff
#endif
C-----------------------------------------------
!
!     POINTER ALIASES
!
!      IF (ns.EQ.1024 .AND. iter2.LE.3) THEN
!        WRITE(100+rank,*) "Iteration: ", iter2
!        CALL PrintOutLinearArray (xc,tlglob,trglob,.TRUE.,100)
!        WRITE(200+rank,*) "Iteration: ", iter2
!        CALL PrintOutLinearArray (gc,tlglob,trglob,.TRUE.,200)
!        WRITE(300+rank,*) "Iteration: ", iter2
!        CALL PrintOutLinearArray (scalxc,tlglob,trglob,.TRUE.,300)
!      END IF

      nfunct3d = nfunct3d + 1
      lu => czmn;  lv => crmn

      CALL second0 (tfunon)

!     CONVERT ODD M TO 1/SQRT(S) INTERNAL REPRESENTATION
!

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
      IF (ictrl_prec2d .le. 1) timer(tfft) = 
     1                         timer(tfft) + (tfftoff - tffton)

!
!     IN HESSIAN LOOP, ONLY COMPUTE PERTURBATIONS TO R, Z, AND L FOR A SINGLE (m,n,ntype)
!
      IF (ictrl_prec2d .ge. 2) THEN
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
         rcon(l,0) = rcon(l,0) + rcon(l,1)*sqrts(l)
         zcon(l,0) = zcon(l,0) + zcon(l,1)*sqrts(l)
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
      IF (iter2.eq.iter1 .and. ivac.le.0) THEN
        DO l = 1, ns
          rcon0(l:nrzt:ns) = rcon(ns:nrzt:ns,0)*sqrts(l:nrzt:ns)**2
          zcon0(l:nrzt:ns) = zcon(ns:nrzt:ns,0)*sqrts(l:nrzt:ns)**2
        END DO
      ENDIF
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
      CALL bcovar (lu, lv, xc(1+2*irzloff+mns*(zsc-1)))
      CALL second0 (tbcovoff)
      s_bcovar_time=s_bcovar_time+(tbcovoff - tbcovon)
      
      IF (ictrl_prec2d .le. 1) timer(tbcov) = timer(tbcov) 
     1                                      + (tbcovoff - tbcovon)

!     COMPUTE VACUUM MAGNETIC PRESSURE AT PLASMA EDGE
!     NOTE: FOR FREE BOUNDARY RUNS, THE VALUE OF RBTOR=R*BTOR
!     AT THE PLASMA EDGE SHOULD BE ADJUSTED TO APPROXIMATELY
!     EQUAL THE VACUUM VALUE. THIS CAN BE DONE BY CHANGING
!     EITHER PHIEDGE OR THE INITIAL CROSS SECTION ACCORDING
!     TO THE SCALING LAW  R*BTOR .EQ. PHIEDGE/(R1 * Z1).
#if defined(MPI_OPT)
      CALL second0(skston)
#endif
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

            IF (ier_flag .ne. 0) RETURN
!
!           RESET FIRST TIME FOR SOFT START
!
            IF (ivac .eq. 1) THEN
               irst = 2;  delt0 = delt
               CALL restart_iter(delt0)
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
!            gcon(:nrzt) = r1(:nrzt,0)+sqrts(:nrzt)*r1(:nrzt,1)
!            gcon(1+nrzt) = 0
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
#if defined(MPI_OPT)
      CALL second0(skstoff)
      s_vacuum_time = s_vacuum_time + (skstoff - skston)
#endif

!
!     COMPUTE CONSTRAINT FORCE
!
      IF (iequi .NE. 1) THEN
        extra1(:nrzt,0) = (rcon(:nrzt,0) - rcon0(:nrzt))*ru0(:nrzt)
     1                   + (zcon(:nrzt,0) - zcon0(:nrzt))*zu0(:nrzt)
        CALL alias (gcon, extra1(:,0), gc, gc(1+mns), gc(1+2*mns), 
     1               extra1(:,1))
      ELSE 
        IF (lrecon) xc(:ns) = xc(:ns) + delr_mse
        GOTO 100
      END IF

!
!     COMPUTE MHD FORCES ON INTEGER-MESH
!
      CALL second0 (tforon)
      CALL forces

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

      gc = gc * scalxc

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
      CALL second0 (tfunoff)
      IF (ictrl_prec2d .le. 1) timer(tfun) = timer(tfun) 
     1                                     + (tfunoff - tfunon)

      END SUBROUTINE funct3d
