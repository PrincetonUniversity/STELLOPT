#if defined (SKS)      
      SUBROUTINE bcovar_par (lu, lv, tpxc)
      USE vmec_main, fpsi => bvco, p5 => cp5
      USE vmec_params, ONLY: ns4, signgs, pdamp, lamscale, ntmax
      USE realspace, ONLY: pextra1, pextra2, pextra3, pextra4,
     1                     pguu, pguv, pgvv, pru, pzu,
     2                     pr1, prv, pzv, pshalf, pwint, pz1,
     3                     pru0, pzu0, psqrts
      USE vforces, r12 => parmn_o, ru12 => pazmn_e, gsqrt => pazmn_o,
     1             rs => pbzmn_e, zs => pbrmn_e, zu12 => parmn_e,
     2             bsubu_e => pclmn_e, bsubv_e => pblmn_e, 
     3             bsubu_o => pclmn_o, bsubv_o => pblmn_o,
     4             bsq => pbzmn_o, phipog => pbrmn_o
      USE vsvd, ONLY: phifac, phifsave, imovephi
      USE xstuff, ONLY: pxc
      USE precon2d, ONLY: ictrl_prec2d, lHess_exact,
     1                    ctor_prec2d
      USE fbal
      USE vmec_input, ONLY: nzeta
      USE vmec_dim, ONLY: ntheta3
      USE parallel_include_module
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(nznt,ns,0:1), INTENT(INOUT) :: lu, lv
      REAL(dp), DIMENSION((1+ntor)*(1+mpol1),1:ns,1:2*ntmax), 
     1  INTENT(IN) :: tpxc
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!     GENERALLY, IF TEMPORAL CONVERGENCE IS POOR, TRY TO INCREASE PDAMP (< 1)
!     (STORED IN VMEC_PARAMS)
      REAL(dp), PARAMETER :: c1p5 = (one + p5)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: l, js, ndim
      REAL(dp) :: r2, volume, curpol_temp
#ifndef _HBANGLE
      REAL(dp) :: arnorm, aznorm, tcon_mul
#endif
      REAL(dp) :: bcastton, bcasttoff
      REAL(dp), POINTER, DIMENSION(:,:) :: luu, luv, lvv, tau
      REAL(dp), DIMENSION(:,:), POINTER :: bsupu, bsubuh, 
     1                                      bsupv, bsubvh, r12sq
      LOGICAL :: lctor
      INTEGER :: i, j, k, nsmin, nsmax, istat
      REAL(dp) :: wblocal(ns), wbtotal
      REAL(dp) :: wplocal(ns), wptotal
      REAL(dp) :: vptotal
      REAL(dp) :: fnlocal(ns), fntotal
      REAL(dp) :: fn1local(ns), fn1total
      REAL(dp) :: fnLlocal(ns), fnLtotal
!-----------------------------------------------
      IF (irst.EQ.2 .AND. iequi.EQ.0) RETURN

!
!     POINTER ALIAS ASSIGNMENTS

      tau => pextra1(:,:,1)
      luu => pextra2(:,:,1)  
      luv => pextra3(:,:,1)
      lvv => pextra4(:,:,1)

      bsupu => luu
      bsubuh => bsubu_o
      bsupv => luv
      bsubvh => bsubv_o
      r12sq => bsq

!
!     FOR OPTIMIZATION ON CRAY, MUST USE COMPILER DIRECTIVES TO
!     GET VECTORIZATION OF LOOPS INVOLVING (MORE THAN ONE) POINTER!
!
      nsmin=t1lglob; nsmax=t1rglob
      pguu(:,nsmin:nsmax) = 0
      pguv(:,nsmin:nsmax) = 0
      pgvv(:,nsmin:nsmax) = 0  

!
!     COMPUTE METRIC ELEMENTS GIJ ON HALF MESH
!     FIRST, GIJ = EVEN PART (ON FULL MESH), LIJ = ODD PART (ON FULL MESH)
!     THEN, GIJ(HALF) = < GIJ(even)> + SHALF < GIJ(odd) >

      DO l = nsmin, nsmax
      r12sq(:,l) = psqrts(:,l)*psqrts(:,l)
      pguu(:,l) = pru(:,l,0)*pru(:,l,0) + pzu(:,l,0)*pzu(:,l,0)
     1           + r12sq(:,l)*(pru(:,l,1)*pru(:,l,1) 
     2           + pzu(:,l,1)*pzu(:,l,1))
      luu(:,l) = (pru(:,l,0)*pru(:,l,1) +  pzu(:,l,0)*pzu(:,l,1))*2
      phipog(:,l) = 2*pr1(:,l,0)*pr1(:,l,1) 
      END DO

      IF (lthreed) THEN
         DO l = nsmin, nsmax
         pguv(:,l) = pru(:,l,0) * prv(:,l,0) + pzu(:,l,0)*pzv(:,l,0)
     1             + r12sq(:,l) * (pru(:,l,1)*prv(:,l,1) 
     2             + pzu(:,l,1)*pzv(:,l,1))
         luv(:,l) = pru(:,l,0) * prv(:,l,1) + pru(:,l,1)*prv(:,l,0)
     1            + pzu(:,l,0)*pzv(:,l,1) + pzu(:,l,1)*pzv(:,l,0)
         pgvv(:,l) = prv(:,l,0) * prv(:,l,0) + pzv(:,l,0)*pzv(:,l,0)
     1             + r12sq(:,l) * (prv(:,l,1)*prv(:,l,1)
     2             + pzv(:,l,1)*pzv(:,l,1) )
         lvv(:,l) = (prv(:,l,0) * prv(:,l,1) + pzv(:,l,0)*pzv(:,l,1))*2
         END DO
      END IF

      r12sq(:,nsmin:nsmax) = pr1(:,nsmin:nsmax,0)
     1                * pr1(:,nsmin:nsmax,0) + r12sq(:,nsmin:nsmax)
     2                * pr1(:,nsmin:nsmax,1) * pr1(:,nsmin:nsmax,1)                       

      DO l = t1rglob, MAX(t1lglob,2), -1
        pguu(:,l) = p5*(pguu(:,l) + pguu(:,l-1) + 
     1                  pshalf(:,l)*(luu(:,l) + luu(:,l-1)))
        r12sq(:,l) = p5*(r12sq(:,l)+r12sq(:,l-1)+pshalf(:,l)*  !Comment: r12sq = r12**2
     1                 (phipog(:,l) + phipog(:,l-1)))                      
      END DO

      IF (lthreed) THEN
         DO l = t1rglob, MAX(t1lglob,2), -1
            pguv(:,l) = p5*(pguv(:,l) + pguv(:,l-1) +
     1         pshalf(:,l)*(luv(:,l) + luv(:,l-1)))
            pgvv(:,l) = p5*(pgvv(:,l) + pgvv(:,l-1) +
     1         pshalf(:,l)*(lvv(:,l) + lvv(:,l-1)))
         END DO
      END IF

      pguv(:,1)=0
      pgvv(:,1)=0

      nsmin=tlglob; nsmax=t1rglob 
      DO l=nsmin, nsmax
      tau(:,l) = gsqrt(:,l)
      gsqrt(:,l) = r12(:,l)*tau(:,l) 
      END DO
      gsqrt(:,1) = gsqrt(:,2)

      nsmin=MAX(2,tlglob); nsmax=t1rglob 
      pgvv(:,nsmin:nsmax)= pgvv(:,nsmin:nsmax)
     1                   + r12sq(:,nsmin:nsmax)
      pgvv(:,1)=0

!CATCH THIS AFTER WHERE LINE BELOW phipog = 0
      nsmin=MAX(2,tlglob); nsmax=t1rglob 
      WHERE (gsqrt(:,nsmin:nsmax) .ne. zero) 
     1   phipog(:,nsmin:nsmax) = one/gsqrt(:,nsmin:nsmax)
      phipog(:,1) = 0

      vp(1) = 0;  vp(ns+1) = 0
      DO js = nsmin, nsmax
         vp(js) = signgs*SUM(gsqrt(:,js)*pwint(:,js))
      END DO

!
!     COMPUTE CONTRA-VARIANT COMPONENTS OF B (Bsupu,v) ON RADIAL HALF-MESH
!     TO ACCOMODATE LRFP=T CASES, THE OVERALL PHIP FACTOR (PRIOR TO v8.46)
!     HAS BEEN REMOVED FROM PHIPOG, SO NOW PHIPOG == 1/GSQRT!
!
!     NOTE: LU = LAMU == d(LAM)/du, LV = -LAMV == -d(LAM)/dv COMING INTO THIS ROUTINE
!     WILL ADD CHIP IN CALL TO ADD_FLUXES. THE NET BSUPU, BSUPV ARE (PHIPOG=1/GSQRT AS NOTED ABOVE):
!
!          BSUPU = PHIPOG*(chip + LAMV*LAMSCALE),
!          BSUPV = PHIPOG*(phip + LAMU*LAMSCALE)
!

      nsmin=t1lglob; nsmax=t1rglob

      DO l = nsmin, nsmax
         lu(:,l,:)=lu(:,l,:)*lamscale
         lv(:,l,:)=lv(:,l,:)*lamscale
         lu(:,l,0)=lu(:,l,0)+phipf(l)
      END DO

      nsmin=MAX(2,t1lglob); nsmax=t1rglob
      DO l = nsmin, nsmax
      bsupu(:,l) = p5*phipog(:,l) * (lv(:,l,0) + lv(:,l-1,0) 
     1           + pshalf(:,l)*(lv(:,l,1) + lv(:,l-1,1)))
      bsupv(:,l) = p5*phipog(:,l) * (lu(:,l,0) + lu(:,l-1,0) 
     1           + pshalf(:,l)*(lu(:,l,1) + lu(:,l-1,1)))
      END DO

!v8.49: add ndim points
      IF (rank .EQ. 0) THEN
        bsupu(:,1) =0; ! bsupu(ndim) = 0
        bsupv(:,1) =0; ! bsupv(ndim) = 0
      END IF

!
!     UPDATE IOTA EITHER OF TWO WAYS:
!     1)  FOR ictrl_prec2d = 0, SOLVE THE LINEAR ALGEBRAIC EQUATION <Bsubu> = icurv 
!         FOR iotas (after testing, this is preferred way)
!     2)  FOR ictrl_prec2d > 0, EVOLVE IOTAS IN TIME, USING Force-iota  = <Bsubu> - icurv. 
!         IOTAS IS "STORED" AT LOCATION LAMBDA-SC(0,0) IN XC-ARRAY
!

!     COMPUTE (IF NEEDED) AND ADD CHIP TO BSUPU
#if defined(CHI_FORCE)
      CALL add_fluxes_par(phipog, bsupu, bsupv, ictrl_prec2d.EQ.0)
#else
      CALL add_fluxes_par(phipog, bsupu, bsupv, .TRUE.)
#endif

!
!     COMPUTE LAMBDA FORCE KERNELS (COVARIANT B COMPONENT bsubu,v) ON RADIAL HALF-MESH
!
      nsmin=t1lglob; nsmax=t1rglob
      DO l = nsmin, nsmax
         bsubuh(:,l)=pguu(:,l)*bsupu(:,l) + pguv(:,l)*bsupv(:,l)
         bsubvh(:,l)=pguv(:,l)*bsupu(:,l) + pgvv(:,l)*bsupv(:,l)
      END DO

!v8.49
!
!     COMPUTE MAGNETIC AND KINETIC PRESSURE ON RADIAL HALF-MESH
!
      nsmin=t1lglob; nsmax=t1rglob
      DO l = nsmin, nsmax
         bsq(:,l) = p5*(bsupu(:,l)*bsubuh(:,l) + bsupv(:,l)*bsubvh(:,l))
      END DO

      nsmin=MAX(2,tlglob); nsmax=MIN(ns,t1rglob)
      pres(nsmin:nsmax) = mass(nsmin:nsmax)/vp(nsmin:nsmax)**gamma
      pres(1)=0

      IF (ictrl_prec2d .LE. 1) THEN
      DO l = tlglob, trglob
         wblocal(l) = SUM(pwint(:,l)*gsqrt(:,l) * bsq(:,l))
         wplocal(l) = vp(l)*pres(l)
      END DO

      CALL Gather1XArray(wblocal)
      wbtotal = SUM(wblocal(2:ns))
      CALL Gather1XArray(wplocal)
      wptotal = SUM(wplocal(2:ns))
      wb = hs*ABS(wbtotal)
      wp = hs*wptotal
      END IF

!     ADD KINETIC PRESSURE TO MAGNETIC PRESSURE
      nsmin=tlglob; nsmax=t1rglob
      DO l=nsmin, nsmax
        bsq(:,l) = bsq(:,l) + pres(l)
        lvv(:,l) = phipog(:,l)*pgvv(:,l)
      END DO

!SPH122407-MOVED HERE: COMPUTE LAMBDA FULL MESH FORCES
!     NOTE: bsubu_e is used here ONLY as a temporary array

      nsmin=tlglob; nsmax=MIN(ns-1,trglob)
      bsubv_e(:,nsmin:nsmax) = p5*(lvv(:,nsmin:nsmax)+
     1             lvv(:,nsmin+1:nsmax+1))*lu(:,nsmin:nsmax,0)
      bsubv_e(:,ns) = p5*lvv(:,ns)*lu(:,ns,0)

      nsmin=tlglob; nsmax=t1rglob
      DO l = nsmin, nsmax
      lvv(:,l) = lvv(:,l) * pshalf(:,l)
      bsubu_e(:,l) = pguv(:,l) * bsupu(:,l) !*sigma_an(:nrzt) !sigma_an=1 isotropic
      END DO

      nsmin=tlglob; nsmax=MIN(ns-1,trglob)
      DO l = nsmin, nsmax
      bsubv_e(:,l) = bsubv_e(:,l) + p5*((lvv(:,l)+lvv(:,l+1))*lu(:,l,1)
     1             + bsubu_e(:,l) + bsubu_e(:,l+1))
      END DO
      bsubv_e(:,ns) = bsubv_e(:,ns) 
     1            + p5*(lvv(:,ns)*lu(:,ns,1) + bsubu_e(:,ns))

!
!     COMPUTE AVERAGE FORCE BALANCE AND TOROIDAL/POLOIDAL CURRENTS
!
!WAC: UPDATE buco, bvco AFTER pressure called (Gather buco, bvco in calc_fbal_par
      CALL calc_fbal_par(bsubuh, bsubvh)
      rbtor0= c1p5*fpsi(2)  - p5*fpsi(3)
      rbtor = c1p5*fpsi(ns) - p5*fpsi(ns-1)
!
!     (SPH:08/19/04)
!     MUST AVOID BREAKING TRI-DIAGONAL RADIAL COUPLING AT EDGE WHEN USING PRECONDITIONER
!     CTOR IS PASSED TO VACUUM TO COMPUTE EDGE BSQVAC, SO IT CAN ONLY DEPEND ON NS, NS-1
!     THUS, CTOR ~ buco(ns) WORKS, WITH REMAINDER A FIXED CONSTANT.
!
!     ALSO, IF USING FAST SWEEP IN COMPUTE_BLOCKS, MUST MAKE CTOR CONSTANT
!     TO AVOID BREAKING SYMMETRY OF A+(ns-1) AND B-(ns) HESSIAN ELEMENTS
!
!     TO GET CORRECT HESSIAN, USE THE CTOR=ctor_prec2d +... ASSIGNMENT
!     FOR ictrl_prec2d.ne.0 (replace ictrl_prec2d.gt.1 with ictrl_prec2d.ne.0 in IF test below)
!
!

!     NEXT COMPUTE COVARIANT BSUBV COMPONENT ~ lvv ON FULL RADIAL MESH BY AVERAGING HALF-MESH METRICS
!     NOTE: EDGE VALUES AT JS=NS DOWN BY 1/2
!     THIS IS NEEDED FOR NUMERICAL STABILITY

      IF (lHess_exact) THEN
         lctor = lfreeb .AND. ictrl_prec2d.NE.0      !Yields correct hessian near edge
      ELSE
         lctor = lfreeb .AND. ictrl_prec2d.GT.1      !Yields better accuracy in solution
      END IF

      IF (lctor) THEN       
         IF (ictrl_prec2d .EQ. 2) 
     1       ctor_prec2d = p5*(buco(ns) - buco(ns1))
         ctor = signgs*twopi*(buco(ns)+ctor_prec2d)
      ELSE
         ctor = signgs*twopi*(c1p5*buco(ns) - p5*buco(ns1))
      END IF

!
!     AVERAGE LAMBDA FORCES ONTO FULL RADIAL MESH
!     USE BLENDING FOR bsubv_e FOR NUMERICAL STABILITY NEAR AXIS
!
      nsmin=tlglob; nsmax=t1rglob
      DO l = nsmin, nsmax
         lvv(:,l) = bdamp(l)
      END DO
     
      IF (rank.EQ.0) THEN
        IF (ANY(bsubvh(:,1) .ne. zero)) STOP 'BSUBVH != 0 AT JS=1'
        IF (ANY(bsubuh(:,1) .ne. zero)) STOP 'BSUBUH != 0 AT JS=1'
      END IF

      nsmin=tlglob; nsmax=MIN(trglob,ns-1)
      bsubu_e(:,nsmin:nsmax) = p5*(bsubuh(:,nsmin:nsmax) + 
     1                             bsubuh(:,nsmin+1:nsmax+1))
      IF (trglob .EQ. ns) bsubu_e(:,ns) = p5*bsubuh(:,ns)

      nsmin=tlglob; nsmax=MIN(ns-1,trglob)
      bsubv_e(:,nsmin:nsmax) = bsubv_e(:,nsmin:nsmax)*
     1        lvv(:,nsmin:nsmax) + p5*(1-lvv(:,nsmin:nsmax))*
     2        (bsubvh(:,nsmin:nsmax) + bsubvh(:,nsmin+1:nsmax+1))
      IF (trglob .EQ. ns) bsubv_e(:,ns) = bsubv_e(:,ns)*lvv(:,ns) + 
     1        p5*(1-lvv(:,ns))*bsubvh(:,ns)

!
!     COMPUTE R,Z AND LAMBDA PRE-CONDITIONING MATRIX
!     ELEMENTS AND FORCE NORMS: NOTE THAT lu=>czmn, lv=>crmn externally
!     SO THIS STORES bsupv in czmn_e, bsupu in crmn_e
!
      nsmin=tlglob; nsmax=t1rglob
      IF (iequi .EQ. 1) THEN
         lu(:,nsmin:nsmax,0) = bsupv(:,nsmin:nsmax)
         lv(:,nsmin:nsmax,0) = bsupu(:,nsmin:nsmax)
      END IF

!
!     COMPUTE PRECONDITIONING (1D) AND SCALING PARAMETERS
!     NO NEED TO RECOMPUTE WHEN 2D-PRECONDITIONER ON
!

      IF ((MOD(iter2-iter1,ns4).EQ.0 .AND. iequi.EQ.0) 
     1        .AND. ictrl_prec2d.EQ.0) THEN
         phifsave = phifac

         nsmin=tlglob; nsmax=t1rglob
         phipog(:,nsmin:nsmax) = phipog(:,nsmin:nsmax)
     1                         * pwint (:,nsmin:nsmax)

         CALL lamcal_par(phipog, pguu, pguv, pgvv)

         CALL precondn_par(bsupv,bsq,gsqrt,r12,zs,zu12,pzu(:,:,0),
     1                 pzu(:,:,1),pz1(:,:,1),arm,ard,brm,brd,
     2                 crd,rzu_fac,cos01)

         CALL precondn_par(bsupv,bsq,gsqrt,r12,rs,ru12,pru(:,:,0),
     1                 pru(:,:,1),pr1(:,:,1),azm,azd,bzm,bzd,
     2                 crd,rru_fac,sin01)

         nsmin=MAX(2,tlglob); nsmax=MIN(trglob,ns-1)
         rzu_fac(nsmin:nsmax) = psqrts(1,nsmin:nsmax)*
     1                          rzu_fac(nsmin:nsmax)
         rru_fac(nsmin:nsmax) = psqrts(1,nsmin:nsmax)*
     1                          rru_fac(nsmin:nsmax)
         frcc_fac(nsmin:nsmax) = one/rzu_fac(nsmin:nsmax)  
         rzu_fac(nsmin:nsmax) = rzu_fac(nsmin:nsmax)/2
         fzsc_fac(nsmin:nsmax) =-one/rru_fac(nsmin:nsmax)
         rru_fac(nsmin:nsmax) = rru_fac(nsmin:nsmax)/2

         nsmin=tlglob; nsmax=t1rglob
         pguu(:,nsmin:nsmax) = pguu(:,nsmin:nsmax)*
     1                         r12(:,nsmin:nsmax)**2

         DO l = MAX(2,tlglob), trglob
            fnlocal(l)  = SUM(pguu(:,l)*pwint(:,l))
            fn1local(l) = SUM(tpxc(2:,l,1:ntmax)**2)
     1                  + SUM(tpxc(1:,l,ntmax+1:2*ntmax)**2)
            fnLlocal(l) = SUM((bsubuh(:,l)**2 + bsubvh(:,l)**2)
     1                  * pwint(:,l))*lamscale**2
         END DO

         CALL Gather1XArray(vp);       vptotal = SUM(vp(2:ns))
         CALL Gather1XArray(fnlocal);  fntotal = SUM(fnlocal(2:ns))
         CALL Gather1XArray(fn1local); fn1total= SUM(fn1local(2:ns))
         CALL Gather1XArray(fnLlocal); fnLtotal= SUM(fnLlocal(2:ns))

         volume = hs*vptotal
         r2 = MAX(wb,wp)/volume
         fnorm = one/(fntotal*(r2*r2))
         fnorm1=one/fn1total
         fnormL = one/fnLtotal
!
!        COMPUTE CONSTRAINT FORCE SCALING FACTOR (TCON)
!        OVERRIDE USER INPUT VALUE HERE
!
#ifndef _HBANGLE
         r2 = ns
         tcon0 = MIN(ABS(tcon0), one)                              !!ignore large tcon0 from old-style files
         tcon_mul = tcon0*(1 + r2*(one/60 + r2/(200*120)))

         tcon_mul = tcon_mul/((4*r0scale**2)**2)                   !!Scaling of ard, azd (2*r0scale**2); 
                                                                   !!Scaling of cos**2 in alias (4*r0scale**2)
         tcon = tcon0
         DO js = MAX(2,tlglob), MIN(ns-1,trglob)
           arnorm = SUM(pwint(:,js)*pru0(:,js)**2)
           aznorm = SUM(pwint(:,js)*pzu0(:,js)**2)
           IF (arnorm.eq.zero .or. aznorm.eq.zero)
     1        STOP 'arnorm or aznorm=0 in bcovar'

           tcon(js) = MIN(ABS(ard(js,1)/arnorm),
     1                    ABS(azd(js,1)/aznorm)) * tcon_mul*(32*hs)**2
         END DO
         tcon(ns) = p5*tcon(ns-1)
         IF (lasym) tcon = p5*tcon
#endif
      ENDIF
!
!     COMPUTE COVARIANT BSUBU,V (EVEN, ODD) ON HALF RADIAL MESH
!     FOR FORCE BALANCE AND RETURN (IEQUI=1)
!
!
!     COMPUTE COVARIANT BSUBU,V (EVEN, ODD) ON HALF RADIAL MESH
!     FOR FORCE BALANCE AND RETURN (IEQUI=1)
!

      IF (iequi .EQ. 1) THEN
        nsmin=MAX(tlglob,2); nsmax=MIN(trglob,ns-1)
        DO js = nsmax, nsmin, -1
          bsubvh(:,js) = 2*bsubv_e(:,js) - bsubvh(:,js+1)
        END DO


!     ADJUST <bsubvh> AFTER MESH-BLENDING
        nsmin=MAX(tlglob,2); nsmax=MIN(trglob,ns)
        DO js = nsmin, nsmax
            curpol_temp = fpsi(js) 
     1                  - SUM(bsubvh(:,js)*pwint(:,js))
            bsubvh(:,js) = bsubvh(:,js) + curpol_temp
        END DO

        bsubu_e(:,nsmin:nsmax) = bsubuh(:,nsmin:nsmax)
        bsubv_e(:,nsmin:nsmax) = bsubvh(:,nsmin:nsmax)

        bsubu_o(:,nsmin:nsmax) = pshalf(:,nsmin:nsmax)*
     1                           bsubu_e(:,nsmin:nsmax)
        bsubv_o(:,nsmin:nsmax) = pshalf(:,nsmin:nsmax)*
     1                           bsubv_e(:,nsmin:nsmax)
         RETURN
      END IF

!     MINUS SIGN => HESSIAN DIAGONALS ARE POSITIVE
      nsmin=MAX(tlglob,2); nsmax=trglob
      bsubu_e(:,nsmin:nsmax) = -lamscale*bsubu_e(:,nsmin:nsmax)
      bsubv_e(:,nsmin:nsmax) = -lamscale*bsubv_e(:,nsmin:nsmax)
      bsubu_o(:,nsmin:nsmax)  =
     1          psqrts(:,nsmin:nsmax)*bsubu_e(:,nsmin:nsmax)
      bsubv_o(:,nsmin:nsmax)  =
     1          psqrts(:,nsmin:nsmax)*bsubv_e(:,nsmin:nsmax)

!
!     STORE LU * LV COMBINATIONS USED IN FORCES
!

      nsmin=MAX(tlglob,2); nsmax=t1rglob
      DO l = nsmin, nsmax
      lvv(:,l) = gsqrt(:,l) !*sigma_an(:,l)
      pguu(:,l)  = bsupu(:,l)*bsupu(:,l)*lvv(:,l)
      pguv(:,l)  = bsupu(:,l)*bsupv(:,l)*lvv(:,l)
      pgvv(:,l)  = bsupv(:,l)*bsupv(:,l)*lvv(:,l)
      lv(:,l,0) = bsq(:,l)*tau(:,l)
      lu(:,l,0) = bsq(:,l)*r12(:,l)
      END DO

      END SUBROUTINE bcovar_par
#endif
      
      SUBROUTINE bcovar (lu, lv)
      USE vmec_main, fpsi => bvco, p5 => cp5
      USE vmec_params, ONLY: ns4, signgs, pdamp, lamscale
      USE realspace
      USE vforces, r12 => armn_o, ru12 => azmn_e, gsqrt => azmn_o,
     1             rs => bzmn_e, zs => brmn_e, zu12 => armn_e,
     2             bsubu_e => clmn_e, bsubv_e => blmn_e, 
     3             bsubu_o => clmn_o, bsubv_o => blmn_o,
     4             bsq => bzmn_o, phipog => brmn_o
      USE vsvd, ONLY: phifac, phifsave, imovephi
      USE xstuff, ONLY: xc
      USE precon2d, ONLY: ictrl_prec2d, lHess_exact,
     1                    ctor_prec2d
#ifdef _HBANGLE
      USE angle_constraints, ONLY: precondn_rho, ard2, arm2,            
     1                             azd2, azm2, brd2, brm2, bzd2, bzm2
#endif
      USE fbal
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(nrzt,0:1), INTENT(inout) :: lu, lv
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!     GENERALLY, IF TEMPORAL CONVERGENCE IS POOR, TRY TO INCREASE PDAMP (< 1)
!     (STORED IN VMEC_PARAMS)
      REAL(dp), PARAMETER :: c1p5 = (one + p5)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: l, js, ndim
      REAL(dp) :: r2, volume, curpol_temp
#ifndef _HBANGLE
      REAL(dp) :: arnorm, aznorm, tcon_mul
#endif
      REAL(dp), POINTER, DIMENSION(:) :: luu, luv, lvv, tau
      REAL(dp), DIMENSION(:), POINTER :: bsupu, bsubuh, 
     1                                      bsupv, bsubvh, r12sq
      LOGICAL :: lctor
!-----------------------------------------------
      ndim = 1+nrzt

!
!     POINTER ALIAS ASSIGNMENTS
!   
      tau => extra1(:,1);  luu => extra2(:,1);  
      luv => extra3(:,1);  lvv => extra4(:,1)

      bsupu => luu;  bsubuh => bsubu_o
      bsupv => luv;  bsubvh => bsubv_o
      r12sq => bsq


!
!     FOR OPTIMIZATION ON CRAY, MUST USE COMPILER DIRECTIVES TO
!     GET VECTORIZATION OF LOOPS INVOLVING (MORE THAN ONE) POINTER!
!
      guu(ndim) = 0;  guv = 0;  gvv = 0  

!
!     COMPUTE METRIC ELEMENTS GIJ ON HALF MESH
!     FIRST, GIJ = EVEN PART (ON FULL MESH), LIJ = ODD PART (ON FULL MESH)
!     THEN, GIJ(HALF) = < GIJ(even)> + SHALF < GIJ(odd) >
!

      r12sq(1:nrzt) = sqrts(1:nrzt)*sqrts(1:nrzt)
      guu(1:nrzt)   = ru(1:nrzt,0)*ru(1:nrzt,0) 
     1              + zu(1:nrzt,0)*zu(1:nrzt,0) + r12sq(1:nrzt)*
     2              ( ru(1:nrzt,1)*ru(1:nrzt,1)
     3            +   zu(1:nrzt,1)*zu(1:nrzt,1))

      luu(1:nrzt)   = (ru(1:nrzt,0)*ru(1:nrzt,1) 
     1              +  zu(1:nrzt,0)*zu(1:nrzt,1))*2
      phipog(1:nrzt)= 2*r1(1:nrzt,0)*r1(1:nrzt,1) 

      IF (lthreed) THEN
         guv(1:nrzt)   = ru(1:nrzt,0)*rv(1:nrzt,0)
     1                 + zu(1:nrzt,0)*zv(1:nrzt,0) + r12sq(1:nrzt)*
     2                 ( ru(1:nrzt,1)*rv(1:nrzt,1) 
     3                 + zu(1:nrzt,1)*zv(1:nrzt,1) )
         luv(1:nrzt)   = ru(1:nrzt,0)*rv(1:nrzt,1) 
     1                 + ru(1:nrzt,1)*rv(1:nrzt,0)
     2                 + zu(1:nrzt,0)*zv(1:nrzt,1) 
     3                 + zu(1:nrzt,1)*zv(1:nrzt,0)
         gvv(1:nrzt)   = rv(1:nrzt,0)*rv(1:nrzt,0)
     1                 + zv(1:nrzt,0)*zv(1:nrzt,0) + r12sq(1:nrzt)*
     2                 ( rv(1:nrzt,1)*rv(1:nrzt,1)
     3                 + zv(1:nrzt,1)*zv(1:nrzt,1) )
         lvv(1:nrzt)   =(rv(1:nrzt,0)*rv(1:nrzt,1) 
     1                 + zv(1:nrzt,0)*zv(1:nrzt,1))*2
      END IF

      r12sq(1:nrzt) = r1(1:nrzt,0)*r1(1:nrzt,0) + r12sq(1:nrzt)*
     1                r1(1:nrzt,1)*r1(1:nrzt,1)                       

!DIR$ IVDEP
      DO l = nrzt, 2, -1
         guu(l) = p5*(guu(l) + guu(l-1) + shalf(l)*(luu(l) + luu(l-1)))
         r12sq(l) = p5*(r12sq(l) + r12sq(l-1) + shalf(l)*             !Comment: r12sq = r12**2
     1                (phipog(l) + phipog(l-1)))                      
      END DO

      IF (lthreed) THEN
!DIR$ IVDEP
         DO l = nrzt, 2, -1
            guv(l) = p5*(guv(l) + guv(l-1) +
     1         shalf(l)*(luv(l) + luv(l-1)))
            gvv(l) = p5*(gvv(l) + gvv(l-1) +
     1         shalf(l)*(lvv(l) + lvv(l-1)))
         END DO
      END IF

      tau(1:nrzt) = gsqrt(1:nrzt)
      gsqrt(1:nrzt) = r12(1:nrzt)*tau(1:nrzt)      
      gsqrt(1:nrzt:ns) = gsqrt(2:nrzt:ns)

      gvv(2:nrzt) = gvv(2:nrzt) + r12sq(2:nrzt)

!CATCH THIS AFTER WHERE LINE BELOW   phipog = 0
      WHERE (gsqrt(2:ndim) .ne. zero) 
     1   phipog(2:ndim) = one/gsqrt(2:ndim)
      phipog(1:ndim:ns) = 0

      vp(1) = 0;  vp(ns+1) = 0
      DO js = 2, ns
         vp(js) = signgs*SUM(gsqrt(js:nrzt:ns)*wint(js:nrzt:ns))
      END DO
      IF (iter2 .eq. 1) voli = twopi*twopi*hs*SUM(vp(2:ns))

!
!     COMPUTE CONTRA-VARIANT COMPONENTS OF B (Bsupu,v) ON RADIAL HALF-MESH
!     TO ACCOMODATE LRFP=T CASES, THE OVERALL PHIP FACTOR (PRIOR TO v8.46)
!     HAS BEEN REMOVED FROM PHIPOG, SO NOW PHIPOG == 1/GSQRT!
!
!     NOTE: LU = LAMU == d(LAM)/du, LV = -LAMV == -d(LAM)/dv COMING INTO THIS ROUTINE
!     WILL ADD CHIP IN CALL TO ADD_FLUXES. THE NET BSUPU, BSUPV ARE (PHIPOG=1/GSQRT AS NOTED ABOVE):
!
!          BSUPU = PHIPOG*(chip - LAMV*LAMSCALE),
!          BSUPV = PHIPOG*(phip + LAMU*LAMSCALE)
!
      lu = lu*lamscale
      lv = lv*lamscale

      DO js=1,ns
         lu(js:nrzt:ns,0)=lu(js:nrzt:ns,0)+phipf(js)
      END DO

      bsupu(2:nrzt) = p5*phipog(2:nrzt)*(lv(2:nrzt,0) + lv(1:nrzt-1,0) 
     1              + shalf(2:nrzt)*(lv(2:nrzt,1) + lv(1:nrzt-1,1)))
      bsupv(2:nrzt) = p5*phipog(2:nrzt)*(lu(2:nrzt,0) + lu(1:nrzt-1,0) 
     1              + shalf(2:nrzt)*(lu(2:nrzt,1) + lu(1:nrzt-1,1)))
!v8.49: add ndim points
      bsupu(1)=0;  bsupu(ndim)=0
      bsupv(1)=0;  bsupv(ndim)=0

!
!     UPDATE IOTA EITHER OF TWO WAYS:
!     1)  FOR ictrl_prec2d = 0, SOLVE THE LINEAR ALGEBRAIC EQUATION <Bsubu> = icurv 
!         FOR iotas  
!     2)  FOR ictrl_prec2d > 0, EVOLVE IOTAS IN TIME, USING Force-iota  = <Bsubu> - icurv IN TOMNSP. 
!
!     NEED TO DO IT WAY (#2) TO EASILY COMPUTE THE HESSIAN ELEMENTS DUE TO LAMBDA-VARIATIONS. 
!     IOTAS IS "STORED" AT LOCATION LAMBDA-SC(0,0) IN XC-ARRAY [USE THIS COMPONENT SO IT 
!     WILL WORK EVEN FOR 2D PLASMA], ALTHOUGH ITS VARIATION IS LIKE THAT OF LV-CS(0,0), 
!     WITH N -> 1 IN THE HESSIAN CALCULATION ROUTINES (Compute_Hessian_Flam_lam, etc.)
!

!     COMPUTE (IF NEEDED) AND ADD CHIP TO BSUPU
#if defined(CHI_FORCE)
      CALL add_fluxes(phipog, bsupu, bsupv, ictrl_prec2d.EQ.0)
#else
      CALL add_fluxes(phipog, bsupu, bsupv, .TRUE.)
#endif

!
!     COMPUTE LAMBDA FORCE KERNELS (COVARIANT B COMPONENT bsubu,v) ON RADIAL HALF-MESH
!
      bsubuh(1:nrzt)=guu(1:nrzt)*bsupu(1:nrzt)+guv(1:nrzt)*bsupv(1:nrzt)
      bsubvh(1:nrzt)=guv(1:nrzt)*bsupu(1:nrzt)+gvv(1:nrzt)*bsupv(1:nrzt)
!v8.49
      bsubuh(ndim) = 0; bsubvh(ndim) = 0

!
!     COMPUTE MAGNETIC AND KINETIC PRESSURE ON RADIAL HALF-MESH
!
      bsq(:nrzt) = p5*(bsupu(:nrzt)*bsubuh(:nrzt) 
     1           +     bsupv(:nrzt)*bsubvh(:nrzt))

      wb = hs*ABS(SUM(wint(:nrzt)*gsqrt(:nrzt)*bsq(:nrzt)))

#ifdef _ANIMEC
!SPH: MAKE CALL HERE (bsubX_e are used for scratch arrays)
      CALL an_pressure(bsubu_e, bsubv_e)

!     ADD KINETIC PRESSURE TO MAGNETIC PRESSURE
      bsq(2:nrzt) = bsq(2:nrzt) + pperp(2:nrzt)

!WAC-SPH: MODIFY EFFECTIVE CURRENT K = curl(sigma_an*B)
      phipog(1:nrzt) = phipog(1:nrzt)*sigma_an(1:nrzt)
      bsubuh(1:nrzt) = bsubuh(1:nrzt)*sigma_an(1:nrzt)
      bsubvh(1:nrzt) = bsubvh(1:nrzt)*sigma_an(1:nrzt)

#else
      pres(2:ns) = mass(2:ns)/vp(2:ns)**gamma
      wp = hs*SUM(vp(2:ns)*pres(2:ns))

!     ADD KINETIC PRESSURE TO MAGNETIC PRESSURE
      DO js=2,ns
         bsq(js:nrzt:ns) = bsq(js:nrzt:ns) + pres(js)
      END DO
#endif

!SPH122407-MOVED HERE: COMPUTE LAMBDA FULL MESH FORCES
!     NOTE: bsubu_e is used here ONLY as a temporary array
      lvv = phipog(:ndim)*gvv
      bsubv_e(1:nrzt) = p5*(lvv(1:nrzt)+lvv(2:ndim))*lu(1:nrzt,0)

      lvv = lvv*shalf
      bsubu_e(:nrzt) = guv(:nrzt)*bsupu(:nrzt)*sigma_an(:nrzt)          !sigma_an=1 isotropic
      bsubu_e(ndim) = 0
      bsubv_e(1:nrzt) = bsubv_e(1:nrzt) 
     1            + p5*((lvv(1:nrzt) + lvv(2:ndim))*lu(1:nrzt,1)
     2            +      bsubu_e(1:nrzt) + bsubu_e(2:ndim))

!
!     COMPUTE AVERAGE FORCE BALANCE AND TOROIDAL/POLOIDAL CURRENTS
!
!WAC: UPDATE buco, bvco AFTER pressure called
#ifdef _ANIMEC
      IF (iequi .EQ. 1) papr = pmap*pres/vp
#endif
      CALL calc_fbal(bsubuh, bsubvh)
    
      rbtor0= c1p5*fpsi(2)  - p5*fpsi(3)
      rbtor = c1p5*fpsi(ns) - p5*fpsi(ns-1)
!
!     (SPH:08/19/04)
!     MUST AVOID BREAKING TRI-DIAGONAL RADIAL COUPLING AT EDGE WHEN USING PRECONDITIONER
!     CTOR IS PASSED TO VACUUM TO COMPUTE EDGE BSQVAC, SO IT CAN ONLY DEPEND ON NS, NS-1
!     THUS, CTOR ~ buco(ns) WORKS, WITH REMAINDER A FIXED CONSTANT.
!
!     ALSO, IF USING FAST SWEEP IN COMPUTE_BLOCKS, MUST MAKE CTOR CONSTANT
!     TO AVOID BREAKING SYMMETRY OF A+(ns-1) AND B-(ns) HESSIAN ELEMENTS
!
!     TO GET CORRECT HESSIAN, USE THE CTOR=ctor_prec2d +... ASSIGNMENT
!     FOR ictrl_prec2d.ne.0 (replace ictrl_prec2d.gt.1 with ictrl_prec2d.ne.0 in IF test below)
!
!

!     NEXT COMPUTE COVARIANT BSUBV COMPONENT ~ lvv ON FULL RADIAL MESH BY AVERAGING HALF-MESH METRICS
!     NOTE: EDGE VALUES AT JS=NS DOWN BY 1/2
!     THIS IS NEEDED FOR NUMERICAL STABILITY

      IF (lHess_exact) THEN
         lctor = lfreeb .and. ictrl_prec2d.ne.0      !Yields correct hessian near edge
      ELSE
         lctor = lfreeb .and. ictrl_prec2d.gt.1      !Yields better accuracy in solution
      END IF
      IF (lctor) THEN       
         IF (ictrl_prec2d .eq. 2) 
     1       ctor_prec2d = signgs*twopi*p5*(buco(ns) - buco(ns1))
         ctor = ctor_prec2d + signgs*twopi*buco(ns)
      ELSE
         ctor = signgs*twopi*(c1p5*buco(ns) - p5*buco(ns1))
      END IF

!
!     AVERAGE LAMBDA FORCES ONTO FULL RADIAL MESH
!     USE BLENDING FOR bsubv_e FOR NUMERICAL STABILITY NEAR AXIS
!
      DO l=1,ns
         lvv(l:nrzt:ns) = bdamp(l)
      END DO

      IF (ANY(bsubuh(1:ndim:ns) .ne. zero)) STOP 'BSUBUH != 0 AT JS=1'
      IF (ANY(bsubvh(1:ndim:ns) .ne. zero)) STOP 'BSUBVH != 0 AT JS=1'

      bsubu_e(1:nrzt) = p5*(bsubuh(1:nrzt) + bsubuh(2:ndim))
      bsubv_e(1:nrzt) = bsubv_e(1:nrzt)*lvv(1:nrzt) + p5*(1-lvv(1:nrzt))
     1                *(bsubvh(1:nrzt) + bsubvh(2:ndim))

!
!     COMPUTE R,Z AND LAMBDA PRE-CONDITIONING MATRIX
!     ELEMENTS AND FORCE NORMS: NOTE THAT lu=>czmn, lv=>crmn externally
!     SO THIS STORES bsupv in czmn_e, bsupu in crmn_e
!
      IF (iequi .EQ. 1) THEN
         lu(:nrzt,0) = bsupv(:nrzt)
         lv(:nrzt,0) = bsupu(:nrzt)
      END IF
 
!
!     COMPUTE PRECONDITIONING (1D) AND SCALING PARAMETERS
!     NO NEED TO RECOMPUTE WHEN 2D-PRECONDITIONER ON
!
      IF ((MOD(iter2-iter1,ns4).eq.0 .and. iequi.eq.0) 
     1        .and. ictrl_prec2d.eq.0) THEN
         phifsave = phifac
         phipog(:nrzt) = phipog(:nrzt)*wint(:nrzt)
         CALL lamcal(phipog, guu, guv, gvv)
         CALL precondn(bsupv,bsq,gsqrt,r12,zs,zu12,zu,zu(1,1),
     1                 z1(1,1),arm,ard,brm,brd,
#ifdef _HBANGLE
     2                 arm2, ard2, brm2, brd2,
#endif
     3                 crd,rzu_fac,cos01)
         CALL precondn(bsupv,bsq,gsqrt,r12,rs,ru12,ru,ru(1,1),
     1                 r1(1,1),azm,azd,bzm,bzd,
#ifdef _HBANGLE
     2                 azm2, azd2, bzm2, bzd2,
#endif
     3                 crd,rru_fac,sin01)

         rzu_fac(2:ns-1) = sqrts(2:ns-1)*rzu_fac(2:ns-1)
         rru_fac(2:ns-1) = sqrts(2:ns-1)*rru_fac(2:ns-1)
         frcc_fac(2:ns-1) = one/rzu_fac(2:ns-1);  rzu_fac = rzu_fac/2
         fzsc_fac(2:ns-1) =-one/rru_fac(2:ns-1);  rru_fac = rru_fac/2
#ifdef _HBANGLE
         CALL precondn_rho
#endif


         volume = hs*SUM(vp(2:ns))
         r2 = MAX(wb,wp)/volume
         guu(:nrzt) = guu(:nrzt)*r12(:nrzt)**2                     !R12 from RP in force
         fnorm = one/(SUM(guu(1:nrzt)*wint(1:nrzt))*(r2*r2))       !Norm, unpreconditioned R,Z forces
         fnorm1 = one/SUM(xc(1+ns:2*irzloff)**2)                   !Norm for preconditioned R,Z forces
         fnormL = one/(SUM((bsubuh(1:nrzt)**2 + bsubvh(1:nrzt)**2)
     1            *wint(1:nrzt))*lamscale**2)                                     !Norm for unpreconditioned Lambda force
!         r3 = one/(2*r0scale)**2
!         fnorm2 = one/MAX(SUM(xc(2*irzloff+1:3*irzloff)**2),r3/4) !Norm for preconditioned Lambda force

!
!        COMPUTE CONSTRAINT FORCE SCALING FACTOR (TCON)
!        OVERRIDE USER INPUT VALUE HERE
!
#ifndef _HBANGLE
         r2 = ns
         tcon0 = MIN(ABS(tcon0), one)                              !!ignore large tcon0 from old-style files
         tcon_mul = tcon0*(1 + r2*(one/60 + r2/(200*120)))

         tcon_mul = tcon_mul/((4*r0scale**2)**2)                   !!Scaling of ard, azd (2*r0scale**2); 
                                                                   !!Scaling of cos**2 in alias (4*r0scale**2)

         DO js = 2, ns-1
           arnorm = SUM(wint(js:nrzt:ns)*ru0(js:nrzt:ns)**2)
           aznorm = SUM(wint(js:nrzt:ns)*zu0(js:nrzt:ns)**2)
           IF (arnorm.eq.zero .or. aznorm.eq.zero)
     1        STOP 'arnorm or aznorm=0 in bcovar'

           tcon(js) = MIN(ABS(ard(js,1)/arnorm),
     1                    ABS(azd(js,1)/aznorm)) * tcon_mul*(32*hs)**2
         END DO
         tcon(ns) = p5*tcon(ns-1)
         IF (lasym) tcon = p5*tcon
#endif
      ENDIF

!
!     COMPUTE COVARIANT BSUBU,V (EVEN, ODD) ON HALF RADIAL MESH
!     FOR FORCE BALANCE AND RETURN (IEQUI=1)
!
      IF (iequi .eq. 1) THEN

!         IF (.FALSE.) THEN
         DO js = ns-1,2,-1
            DO l = js, nrzt, ns
               bsubvh(l) = 2*bsubv_e(l) - bsubvh(l+1)
            END DO
         END DO

!     ADJUST <bsubvh> AFTER MESH-BLENDING
         DO js = 2, ns
            curpol_temp = fpsi(js) 
     1                  - SUM(bsubvh(js:nrzt:ns)*wint(js:nrzt:ns))
           DO l = js, nrzt, ns
               bsubvh(l) = bsubvh(l) + curpol_temp
            END DO
         END DO
!         END IF

         bsubu_e(:nrzt) = bsubuh(:nrzt)
         bsubv_e(:nrzt) = bsubvh(:nrzt)

         bsubu_o(:nrzt) = shalf(:nrzt)*bsubu_e(:nrzt)
         bsubv_o(:nrzt) = shalf(:nrzt)*bsubv_e(:nrzt)

         RETURN

      END IF

!     MINUS SIGN => HESSIAN DIAGONALS ARE POSITIVE
      bsubu_e = -lamscale*bsubu_e
      bsubv_e = -lamscale*bsubv_e
      bsubu_o(:nrzt)  = sqrts(:nrzt)*bsubu_e(:nrzt)
      bsubv_o(:nrzt)  = sqrts(:nrzt)*bsubv_e(:nrzt)

!
!     STORE LU * LV COMBINATIONS USED IN FORCES
!
!WAC, SPH122407: sigma_an (=1 for isotropic case)
      lvv(2:nrzt) = gsqrt(2:nrzt)*sigma_an(2:nrzt)
      guu(2:nrzt)  = bsupu(2:nrzt)*bsupu(2:nrzt)*lvv(2:nrzt)
      guv(2:nrzt)  = bsupu(2:nrzt)*bsupv(2:nrzt)*lvv(2:nrzt)
      gvv(2:nrzt)  = bsupv(2:nrzt)*bsupv(2:nrzt)*lvv(2:nrzt)
      lv(2:nrzt,0) = bsq(2:nrzt)*tau(2:nrzt)
      lu(2:nrzt,0) = bsq(2:nrzt)*r12(2:nrzt)

      END SUBROUTINE bcovar

