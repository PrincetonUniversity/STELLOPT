      SUBROUTINE jxbforce(bsupu, bsupv, bsubu, bsubv, bsubs, bsubsu,
     1   bsubsv, gsqrt, bsq, itheta, izeta, brho, sigma_an, ier_flag 
#ifdef _ANIMEC
     2  ,pp1, pp2, ppar, onembc
#endif
     3                   )
      USE safe_open_mod
      USE vmec_main
      USE vmec_params, ONLY: mscale, nscale, signgs, mnyq, nnyq, 
     1                       successful_term_flag
      USE realspace,   ONLY: shalf, wint, guu, guv, gvv, r1, ru, rv,
     1                       zu   , zv  , phip
#ifdef _ANIMEC
     2                      ,pp3 
      USE fbal, ONLY: bimax_ppargrad
#endif
!!#undef NETCDF
#ifdef NETCDF
      USE ezcdf
#endif
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(ns,nznt), INTENT(in) ::
     1  bsupu, bsupv, bsq, gsqrt, sigma_an
#ifdef _ANIMEC
     2 ,ppar, onembc
#endif
      REAL(rprec), DIMENSION(ns,nznt,0:1), TARGET, INTENT(inout) ::
     1  bsubu, bsubv
      REAL(rprec), DIMENSION(ns,nznt), INTENT(inout), TARGET :: bsubs
      REAL(rprec), DIMENSION(ns,nznt), INTENT(out) ::
     1  itheta, brho, izeta
#ifdef _ANIMEC
     1 ,pp1, pp2
#endif
      REAL(rprec), DIMENSION(ns,nznt,0:1) :: bsubsu, bsubsv
      INTEGER, INTENT(in) :: ier_flag
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!RESET lbsubs DEFAULT FLAG TO FALSE TO CAPTURE CURRENT SHEETS!
!      LOGICAL, PARAMETER :: lbsubs = .false.      !!False to use (correct)  bsubs calculation (from metrics)
                                                  !!True  to use (modified) bsubs calculation (from mag. diff. eq.)
!  J Hanson 2014-01-12. Commented out above line. Variable is now declared
!    in module vmec_input, available here through module vmec_main.
!    lbsubs is now a namelist input variable, so user can change.
!      LOGICAL, PARAMETER :: lbsubs = .true.       !!True  to use NEW bsubs calculation (from mag. diff. eq.)
!                                                  !!False to use OLD bsubs calculation (from metrics)
      LOGICAL, PARAMETER :: lprint = .false.      !!Prints out bsubs spectrum to fort.33
      REAL(rprec), PARAMETER :: two=2, p5=0.5_dp, c1p5=1.5_dp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER lk, lz, lt, k, m, js, j, n, injxbout, mparity, nznt1
      INTEGER :: njxbout = jxbout0, nmin, info 
      INTEGER, PARAMETER :: ns_skip = 1, nu_skip = 1, nv_skip = 1
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::
     1    bdotk, bsubuv, bsubvu
      REAL(rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: bsubsmn
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE ::   brhomn,
     1     bsubs3, bsubv3, bsubu3, jxb_gradp, jcrossb, sqrtg3,
     2     bsupv3, bsupu3, jsups3, jsupv3, jsupu3, jdotb_sqrtg
      REAL(rprec), POINTER :: bs1(:), bu1(:,:), bv1(:,:)
      REAL(rprec), DIMENSION(:), ALLOCATABLE     :: kperpu, kperpv, 
     2    sqgb2, sqrtg, kp2, jxb, jxb2, bsupu1, bsupv1, bsubu1, bsubv1,
     3    avforce, aminfor, amaxfor, toroidal_angle, phin, pprim, pprime
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE   :: bsubua, bsubva
      REAL(rprec) ::
     1    bsubsmn1, bsubsmn2, bsubvmn1, bsubvmn2, bsubumn1, bsubumn2,
     1    bsubsmn3, bsubsmn4, bsubvmn3, bsubvmn4, bsubumn3, bsubumn4,
     2    dnorm1, tcos1, tcos2, tsini1, tsini2, tcosi1, tcosi2, 
     3    tcosm1, tcosm2, tcosn1, tcosn2, tsinm1, tsinm2, tsin1, tsin2,
     4    tsinn1, tsinn2, tjnorm, ovp, pnorm, brho00(ns)
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::
     1    bsubu_s, bsubu_a, bsubv_s, bsubv_a
      REAL(rprec), DIMENSION(:), ALLOCATABLE ::
     1    bsubs_s, bsubs_a
      CHARACTER(LEN=100) :: jxbout_file
      CHARACTER(LEN=100) :: legend(13)
      LOGICAL :: lprint_flag
!-----------------------------------------------
#ifdef NETCDF
      CHARACTER(LEN=*), PARAMETER ::
     1  vn_legend = 'legend',
     1  vn_radial_surfaces = 'radial_surfaces',
     1  vn_poloidal_grid_points = 'poloidal_grid_points',
     1  vn_toroidal_grid_points = 'toroidal_grid_points',
     1  vn_mpol = 'mpol',
     1  vn_ntor = 'ntor',
     1  vn_phin = 'phin',
     1  vn_toroidal_angle = 'toroidal_angle',
     1  vn_avforce = 'avforce',
     1  vn_jdotb = 'surf_av_jdotb',
     1  vn_sqg_bdotk = 'sqrt(g)*bdotk',
     1  vn_sqrtg = 'sqrt(g)',
     1  vn_bdotgradv = 'bdotgradv',
     1  vn_amaxfor = 'amaxfor',
     1  vn_aminfor = 'aminfor',
     1  vn_pprime = 'pprime',
     1  vn_jsupu = 'jsupu',
     1  vn_jsupv = 'jsupv',
     1  vn_jsups = 'jsups',
     1  vn_bsupu = 'bsupu',
     1  vn_bsupv = 'bsupv',
     1  vn_jcrossb = 'jcrossb',
     1  vn_jxb_gradp = 'jxb_gradp',
     1  vn_bsubu = 'bsubu',
     1  vn_bsubv = 'bsubv',
     1  vn_bsubs = 'bsubs'
!-----------------------------------------------
#endif
      lprint_flag = (ier_flag.eq.successful_term_flag)
      IF (lprint_flag) THEN
#ifdef NETCDF
      jxbout_file = 'jxbout_'//TRIM(input_extension)//'.nc'

      CALL cdf_open(njxbout,jxbout_file,'w',injxbout)
#else
      jxbout_file = 'jxbout.'//TRIM(input_extension)//'.txt'
      CALL safe_open(njxbout, injxbout, jxbout_file, 'replace',
     1    'formatted')
#endif
#if !defined(_ANIMEC)
      IF (ANY(sigma_an .NE. one)) STOP 'SIGMA_AN != 1'
#endif
      IF (injxbout .ne. 0) THEN
         PRINT *,' Error opening JXBOUT file in jxbforce'
         RETURN
      END IF

!     PROGRAM FOR COMPUTING LOCAL KXB = grad-p FORCE BALANCE
!
!     sigma_an = one + (pperp-ppar)/(2*bsq(1:nrzt))   (=1 for isotropic pressure)
!     K = CURL(sigma_an B) is the "effective" current (=J for isotropic pressure)
!     Compute u (=theta), v (=zeta) derivatives of B sub s
!
      legend(1) = " S = normalized toroidal flux (0 - 1)"
      IF (lasym) THEN
         legend(2) = " U = VMEC poloidal angle (0 - 2*pi, FULL period)"
	ELSE
         legend(2) = " U = VMEC poloidal angle (0 - pi, HALF a period)"
	END IF
      legend(3) = " V = VMEC (geometric) toroidal angle (0 - 2*pi)"
      legend(4) = " SQRT(g') = |SQRT(g-VMEC)| / VOL':" //
     1  " Cylindrical-to-s,u,v Jacobian normed to volume derivative"
      legend(5) = " VOL = Int_s'=0,s Int_u Int_v |SQRT(g_VMEC)| :" //
     1  " plasma volume  enclosed by surface s'=s"
      legend(6) = " VOL' = d(VOL)/ds: differential volume element"
      legend(7) = " Es = SQRT(g') [grad(U) X grad(V)] : covariant" //
     1   " radial unit vector (based on volume radial coordinate)"
      legend(8) = " BSUP{U,V} = sigma_an B DOT GRAD{U,V}:" //
     1   "  contravariant components of B"
      legend(9) = " JSUP{U,V} = SQRT(g') J DOT GRAD{U,V}"
      legend(10)=
     1  " K X B = Es DOT [K X B]: covariant component of K X B force"
      legend(11)= " K * B = K DOT B * SQRT(g')"
      legend(12)= " p' = dp/d(VOL): pressure gradient (based on" //
     1  " volume radial coordinate)"
      legend(13)= " <KSUP{U,V}> = Int_u Int_v [KSUP{U,V}]/dV/ds"

#if !defined(NETCDF)
      WRITE (njxbout,5) (ns1-1)/ns_skip, ntheta3/nu_skip, nzeta/nv_skip,
     1    mpol, ntor, phiedge
 5    FORMAT(/,' Radial surfaces = ',i3, ' Poloidal grid points = ',i3,
     1         ' Toroidal grid points = ',i3,/,
     2         ' Poloidal modes = ',i3,' Toroidal modes = ', i3,
     3         ' Toroidal Flux  = ',1pe12.3)
      WRITE (njxbout, 6) (legend(j), j=1,13)
 6    FORMAT(/,100('-'),/,' LEGEND:',/,100('-'),/,
     1  2(3(a,/),/),5(a,/),/,2(a,/),100('-'),//)
#endif
      ENDIF

      nznt1 = nzeta*ntheta2
      ALLOCATE (avforce(ns),aminfor(ns),amaxfor(ns))
      ALLOCATE (bdotk(ns,nznt), bsubuv(ns,nznt),
     1          bsubvu(ns,nznt), kperpu(nznt), kperpv(nznt), 
     2          sqgb2(nznt), brhomn(0:mnyq,-nnyq:nnyq,0:1),kp2(nznt),
     3          jxb(nznt), jxb2(nznt), bsupu1(nznt), 
     3          bsubua(nznt1,0:1), bsubva(nznt1,0:1), 
     4          bsupv1(nznt), bsubu1(nznt), bsubv1(nznt),
     5          bsubsmn(ns,0:mnyq,-nnyq:nnyq,0:1),
     6          bsubs_s(nznt), bsubs_a(nznt), sqrtg(nznt),
     7          bsubu_s(nznt1,0:1), bsubu_a(nznt1,0:1),
     8          bsubv_s(nznt1,0:1), bsubv_a(nznt1,0:1), stat=j)
      IF (j .ne. 0) STOP 'Allocation error in jxbforce'

!
!     NOTE: bsubuv, bsubvu are used to compute the radial current (should be zero)
!
      bsubsu = 0; bsubsv = 0; bsubuv = 0; bsubvu = 0; bdotk  = 0
      bsubs(1,:) = 0; bsubsmn = 0

      radial: DO js = 1, ns
!
!     Put bsubs on full mesh
!
         IF (js.gt.1 .and. js.lt.ns) THEN     
            bsubs(js,:) = p5*(bsubs(js,:) + bsubs(js+1,:))
         END IF

         bsubu(js,:,1) = bsubu(js,:,1)/shalf(js)
         bsubv(js,:,1) = bsubv(js,:,1)/shalf(js)
         bsubua = 0;   bsubva = 0

!        _s: symmetric in u,v  _a: antisymmetric in u,v on half (ntheta2) interval
         IF (lasym)  THEN
            bs1=>bsubs(js,:); bu1=>bsubu(js,:,:); bv1=>bsubv(js,:,:)
            CALL fsym_fft (bs1, bu1, bv1, bsubs_s, bsubu_s, bsubv_s, 
     1                     bsubs_a, bsubu_a, bsubv_a)
         ELSE
            bsubs_s(:) = bsubs(js,:)
            bsubu_s = bsubu(js,:,:); bsubv_s = bsubv(js,:,:)
         END IF

!
!        FOURIER LOW-PASS FILTER bsubX

         DO m = 0, mpol1
            mparity = MOD(m, 2)
            DO n = 0, ntor
!
!        FOURIER TRANSFORM
!
               dnorm1 = one/r0scale**2
               IF (m .eq. mnyq) dnorm1 = p5*dnorm1
               IF (n.eq.nnyq .and. n.ne.0) dnorm1 = p5*dnorm1
               bsubsmn1 = 0;  bsubsmn2 = 0
               IF (lasym) THEN
                  bsubsmn3 = 0;  bsubsmn4 = 0
               END IF
               bsubumn1 = 0;  bsubumn2 = 0;  bsubvmn1 = 0;  bsubvmn2 = 0
               IF (lasym) THEN
                  bsubumn3 = 0; bsubumn4 = 0; bsubvmn3 = 0; bsubvmn4 = 0
               END IF

               DO k = 1, nzeta
                  lk = k
                  DO j = 1, ntheta2
                     tsini1 = sinmui(j,m)*cosnv(k,n)*dnorm1
                     tsini2 = cosmui(j,m)*sinnv(k,n)*dnorm1
                     tcosi1 = cosmui(j,m)*cosnv(k,n)*dnorm1
                     tcosi2 = sinmui(j,m)*sinnv(k,n)*dnorm1
                     bsubsmn1 = bsubsmn1 + tsini1*bsubs_s(lk)
                     bsubsmn2 = bsubsmn2 + tsini2*bsubs_s(lk)
                     bsubvmn1 = bsubvmn1 + tcosi1*bsubv_s(lk, mparity)
                     bsubvmn2 = bsubvmn2 + tcosi2*bsubv_s(lk, mparity)
                     bsubumn1 = bsubumn1 + tcosi1*bsubu_s(lk, mparity)
                     bsubumn2 = bsubumn2 + tcosi2*bsubu_s(lk, mparity)

                     IF (lasym) THEN
                     bsubsmn3 = bsubsmn3 + tcosi1*bsubs_a(lk)
                     bsubsmn4 = bsubsmn4 + tcosi2*bsubs_a(lk)
                     bsubvmn3 = bsubvmn3 + tsini1*bsubv_a(lk, mparity)
                     bsubvmn4 = bsubvmn4 + tsini2*bsubv_a(lk, mparity)
                     bsubumn3 = bsubumn3 + tsini1*bsubu_a(lk, mparity)
                     bsubumn4 = bsubumn4 + tsini2*bsubu_a(lk, mparity)
                     END IF

                     lk = lk + nzeta

                  END DO
               END DO

!
!              FOURIER INVERSE TRANSFORM
!              Compute on u-v grid (must add symmetric, antisymmetric parts for lasym=T)
! 
               DO k = 1, nzeta
                  lk = k
                  DO j = 1, ntheta2
                     tcos1 = cosmu(j,m)*cosnv(k,n)
                     tcos2 = sinmu(j,m)*sinnv(k,n)
                     bsubua(lk,0) = bsubua(lk,0) + tcos1*bsubumn1 +
     1                  tcos2*bsubumn2
                     bsubva(lk,0) = bsubva(lk,0) + tcos1*bsubvmn1 +
     1                  tcos2*bsubvmn2

                     tcosm1 = cosmum(j,m)*cosnv(k,n)
                     tcosm2 = sinmum(j,m)*sinnv(k,n)
                     bsubsu(js,lk,0) = bsubsu(js,lk,0) +
     1                  tcosm1*bsubsmn1 + tcosm2*bsubsmn2
                     tcosn1 = sinmu(j,m)*sinnvn(k,n)
                     tcosn2 = cosmu(j,m)*cosnvn(k,n)
                     bsubsv(js,lk,0) = bsubsv(js,lk,0) +
     1                  tcosn1*bsubsmn1 + tcosn2*bsubsmn2
                     bsubvu(js,lk) = bsubvu(js,lk) + 
     1                               sinmum(j,m)*cosnv(k,n)*bsubvmn1 +
     2                               cosmum(j,m)*sinnv(k,n)*bsubvmn2
                     bsubuv(js,lk) = bsubuv(js,lk) + 
     1                               cosmu(j,m)*sinnvn(k,n)*bsubumn1 +
     2                               sinmu(j,m)*cosnvn(k,n)*bsubumn2

                     IF (lasym) THEN
                     tsin1 = sinmu(j,m)*cosnv(k,n)
                     tsin2 = cosmu(j,m)*sinnv(k,n)
                     bsubua(lk,1) = bsubua(lk,1) + tsin1*bsubumn3 +
     1                  tsin2*bsubumn4
                     bsubva(lk,1) = bsubva(lk,1) + tsin1*bsubvmn3 +
     1                  tsin2*bsubvmn4

                     tsinm1 = sinmum(j,m)*cosnv(k,n)
                     tsinm2 = cosmum(j,m)*sinnv(k,n)
                     bsubsu(js,lk,1) = bsubsu(js,lk,1) +
     1                   tsinm1*bsubsmn3 + tsinm2*bsubsmn4
                     tsinn1 = cosmu(j,m)*sinnvn(k,n)
                     tsinn2 = sinmu(j,m)*cosnvn(k,n)
                     bsubsv(js,lk,1) = bsubsv(js,lk,1) +
     1                   tsinn1*bsubsmn3 + tsinn2*bsubsmn4
                     bsubvu(js,lk) = bsubvu(js,lk) + 
     1                               cosmum(j,m)*cosnv(k,n)*bsubvmn3 +
     2                               sinmum(j,m)*sinnv(k,n)*bsubvmn4
                     bsubuv(js,lk) = bsubuv(js,lk) + 
     1                               sinmu(j,m)*sinnvn(k,n)*bsubumn3 +
     2                               cosmu(j,m)*cosnvn(k,n)*bsubumn4
                     END IF

                     lk = lk + nzeta

                  END DO
               END DO

!
!              bsubsmn: coefficients of sin(mu)cos(nv), n>=0, cos(mu)sin(nv), n<0 (type=0)
!                                       cos(mu)cos(nv), n>=0, sin(mu)sin(nv), n<0 (type=1, nonzero only for lasym=T)
!
               IF (.not.lprint) CYCLE          !Don't need these except for comparison

               bsubsmn(js,m,n,0) = bsubsmn1
               IF (n .gt. 0) bsubsmn(js,m,-n,0) = bsubsmn2
             
               IF (.not.lasym) CYCLE

               bsubsmn(js,m,n,1) = bsubsmn3
               IF (n .gt. 0) bsubsmn(js,m,-n,0) = bsubsmn4

            END DO
         END DO

         IF (lasym) THEN
!           EXTEND FILTERED bsubu, bsubv TO NTHETA3 MESH
!           NOTE: INDEX 0 - COS(mu-nv) SYMMETRY; 1 - SIN(mu-nv) SYMMETRY
            CALL fext_fft (bsubu(js,:,0), bsubua(:,0), bsubua(:,1))
            CALL fext_fft (bsubv(js,:,0), bsubva(:,0), bsubva(:,1))
         ELSE
            bsubu(js,:,0) = bsubua(:,0)
            bsubv(js,:,0) = bsubva(:,0)
         END IF

      END DO radial

      DEALLOCATE (bsubua, bsubva)

!     EXTEND bsubsu, bsubsv TO NTHETA3 MESH
      IF (lasym) CALL fsym_invfft (bsubsu, bsubsv)

#ifdef _ANIMEC
      CALL bimax_ppargrad(pp1, pp2, gsqrt, ppar, onembc, pres, 
     1                    phot,tpotb)
#endif

!     SKIPS Bsubs Correction - uses Bsubs from metric elements
      IF (.not.lbsubs) GOTO 1500          

!
!     Compute corrected Bsubs coefficients (brhomn) (impacts currents)
!     by solving es dot (KXB - gradp_parallel) = 0 equation for brhomn in REAL SPACE
!     Can be written Bsupu d(bs)/du + Bsupv d(bs)/dv = RHS (jxb below), bs==bsubs
!     brho==sigma B_s, pp1 and pp2 are the Jacobian times the hot particle parallel
!     pressure radial gradient Amplitudes on the full integer mesh
!
      correct_bsubs: DO js = 2, ns-1
         jxb(:) = p5*(gsqrt(js,:) + gsqrt(js+1,:))
         bsupu1(:) = p5*(bsupu(js,:)*gsqrt(js,:)
     1    + bsupu(js+1,:)*gsqrt(js+1,:))
         bsupv1(:) = p5*(bsupv(js,:)*gsqrt(js,:)
     1    + bsupv(js+1,:)*gsqrt(js+1,:))
         brho(js,:) = ohs*
     1   ( bsupu1(:)*(bsubu(js+1,:,0) - bsubu(js,:,0))
     2   + bsupv1(:)*(bsubv(js+1,:,0) - bsubv(js,:,0)))
     3   + (pres(js+1) - pres(js))*ohs*jxb(:)
#ifdef _ANIMEC
!WAC Last two lines of brho contain hot particle parallel pressure gradients
     4   + ohs*((pres(js+1)*phot(js+1) - pres(js)*phot(js)) * pp2(js,:)
     5   +      (tpotb(js+1)           - tpotb(js)      ) * pp1(js,:))
#endif
!
!     SUBTRACT FLUX-SURFACE AVERAGE FORCE BALANCE FROM brho, OTHERWISE
!     LOCAL FORCE BALANCE EQUATION B dot grad(Bs) = brho CAN'T BE SOLVED
!
         brho00(js) = SUM(brho(js,:)*wint(js:nrzt:ns))
         brho(js,:) = brho(js,:) - signgs*jxb(:)*brho00(js)/
     1      (p5*(vp(js) + vp(js+1)))

         jxb(:) = brho(js,:)
         CALL getbsubs (brhomn, jxb, bsupu1, bsupv1, mnyq, nnyq, info)
         IF (info .ne. 0) THEN
            PRINT *, 'Error in GETBRHO: info= ',info, ' js= ',js
         ELSE IF (lprint) THEN
            WRITE (33, *) ' JS = ', js
            IF (lasym) THEN
            WRITE (33, '(a)') 
     1      '   M    N        BSUBS(old)        BSUBS(new)' //
     2      '        BSUBS(old)        BSUBS(new)'
            ELSE
            WRITE (33, *) '  M    N        BSUBS(old)        BSUBS(new)'
            END IF
            DO m = 0, mpol1
               DO n = -ntor, ntor
                  IF (lasym) THEN
                  WRITE(33,1223) m, n, bsubsmn(js,m,n,0), brhomn(m,n,0),
     1                                 bsubsmn(js,m,n,1), brhomn(m,n,1)
                  ELSE
                  WRITE(33,1224) m, n, bsubsmn(js,m,n,0), brhomn(m,n,0)
                  END IF
               END DO
            END DO
         END IF
 1223    FORMAT (i4,1x,i4,4(6x,1p,e12.3))
 1224    FORMAT (i4,1x,i4,2(6x,1p,e12.3))

!
!        Recompute bsubsu,v now using corrected bsubs
!        Store old values (itheta,izeta) for checking force balance later
!
         itheta(js,:) = bsubsu(js,:,0);  izeta (js,:) = bsubsv(js,:,0)

         IF (info .ne. 0) CYCLE
         bsubsu(js,:,:) = 0;   bsubsv(js,:,:) = 0;  bsubs_s = 0
         IF (lasym) bsubs_a = 0;

         DO m = 0, mnyq
            DO n = 0, nnyq
               IF (n .eq. 0) THEN
                  bsubsmn1 = brhomn(m,0,0)
                  bsubsmn2 = 0
               ELSE
                  bsubsmn1 = brhomn(m,n,0)
                  bsubsmn2 = brhomn(m,-n,0)
               END IF

               IF (lasym) THEN
                  IF (n .eq. 0) THEN
                     bsubsmn3 = brhomn(m,0,1)
                     bsubsmn4 = 0
                  ELSE
                     bsubsmn3 = brhomn(m,n,1)
                     bsubsmn4 = brhomn(m,-n,1)
                  END IF
               END IF

               DO k = 1, nzeta
                  lk = k
                  DO j = 1, ntheta2
                     tsin1 = sinmu(j,m)*cosnv(k,n)
                     tsin2 = cosmu(j,m)*sinnv(k,n)
                     bsubs_s(lk) = bsubs_s(lk) + tsin1*bsubsmn1
     1                                         + tsin2*bsubsmn2
                     tcosm1 = cosmum(j,m)*cosnv(k,n)
                     tcosm2 = sinmum(j,m)*sinnv(k,n)
                     bsubsu(js,lk,0) = bsubsu(js,lk,0) +
     1                  tcosm1*bsubsmn1 + tcosm2*bsubsmn2
                     tcosn1 = sinmu(j,m)*sinnvn(k,n)
                     tcosn2 = cosmu(j,m)*cosnvn(k,n)
                     bsubsv(js,lk,0) = bsubsv(js,lk,0) +
     1                  tcosn1*bsubsmn1 + tcosn2*bsubsmn2
                     
                     IF (lasym) THEN
                     tcos1 = cosmu(j,m)*cosnv(k,n)
                     tcos2 = sinmu(j,m)*sinnv(k,n)
                     bsubs_a(lk) = bsubs_a(lk) + tcos1*bsubsmn3
     1                                         + tcos2*bsubsmn4
                     tsinm1 = sinmum(j,m)*cosnv(k,n)
                     tsinm2 = cosmum(j,m)*sinnv(k,n)
                     bsubsu(js,lk,1) = bsubsu(js,lk,1) +
     1                  tsinm1*bsubsmn3 + tsinm2*bsubsmn4
                     tsinn1 = cosmu(j,m)*sinnvn(k,n)
                     tsinn2 = sinmu(j,m)*cosnvn(k,n)
                     bsubsv(js,lk,1) = bsubsv(js,lk,1) +
     1                  tsinn1*bsubsmn3 + tsinn2*bsubsmn4

                     END IF
 
                     lk = lk + nzeta

                  END DO
               END DO
            END DO
         END DO

         IF (lasym) THEN
!           EXTEND TO FULL (ntheta3) u-GRID
            bs1 => bsubs(js,:)
            CALL fext_fft (bs1, bsubs_a, bsubs_s)
         ELSE
            bsubs(js,:) = bsubs_s(:)
         END IF

      END DO correct_bsubs

!     EXTEND bsubsu, bsubsv TO NTHETA3 MESH
      IF (lasym) CALL fsym_invfft (bsubsu, bsubsv)

!
!     CHECK FORCE BALANCE: SQRT(g)*(bsupu*bsubsu + bsupv*bsubsv) = brho
!
      IF (.not.lprint) GOTO 1500

      WRITE (33, '(/,2a,/)') 'ANGLE INDEX       B*grad(Bs)      Frhs',
     1              '          Fold'
      check_fb: DO js = 2, ns-1
         bsupu1(:) = p5*(bsupu(js,:)*gsqrt(js,:)
     1             +     bsupu(js+1,:)*gsqrt(js+1,:))
         bsupv1(:) = p5*(bsupv(js,:)*gsqrt(js,:)
     1             + bsupv(js+1,:)*gsqrt(js+1,:))
         kp2(:) = bsupu1(:)*bsubsu(js,:,0) + bsupv1(:)*bsubsv(js,:,0)
         jxb(:) = bsupu1(:)*itheta(js,:) + bsupv1(:)*izeta(js,:)

         WRITE (33, '(/,a,i4)') 'JS = ',js
         DO lk = 1, nznt
            WRITE(33,1230) lk, brho(js,lk),  kp2(lk),  jxb(lk)
         END DO

      END DO check_fb

 1230 FORMAT (i9,5x, 1p,3e14.4)

 1500 CONTINUE

      DEALLOCATE (bsubs_s, bsubs_a, bsubu_s,
     1            bsubu_a, bsubv_s, bsubv_a, stat=lk)

!
!     Compute end point values for bsubs
!
      bsubs(1,:)  = 2*bsubs(2,:)  - bsubs(3,:)
      bsubs(ns,:) = 2*bsubs(ns,:) - bsubs(ns-1,:)
!
!     Now compute currents on the FULL radial mesh
!     Here:
!
!     Itheta = sqrt(g) * Ksupu
!     Izeta  = sqrt(g) * Ksupv
!     Ksupx  = K dot grad(x)                          x=(u,v)
!     jxb    = (K X B) dot (grad-u X grad-v) sqrt(g)  
!     bdotk  = sigma*sqrt(g)*K dot B
!     kperpx = (B X gradp) dot grad(x) / |B|**2       x=(u,v)
!     sqgb2  = sigma*sqrt(g)*|B|**2
!     sqrtg  = sqrt(g)
!     pprime = d(p||)/dV
!
!     kp2   == |k-perp|**2 = kperpu**2 * guu + 2*kperpu*kperpv*guv + kperpv**2 * gvv
!     This was compared to the alternative expression (agreed very well):
!     |j-perp|**2 = |grad-s|**2 * (dp/ds)**2 / |B|**2
!
!     Note: Multiply currents, pressure by 1/mu0 to get in mks units!
!           TWOPI*TWOPI factor incorporated in vp (thru ovp factor below), so V' = (2pi)**2*vp
!
#ifdef NETCDF
      ALLOCATE(
     1     bsubs3(ns,nzeta,ntheta3), bsubv3(ns,nzeta,ntheta3), 
     2     bsubu3(ns,nzeta,ntheta3), jxb_gradp(ns,nzeta,ntheta3), 
     3     jcrossb(ns,nzeta,ntheta3), bsupv3(ns,nzeta,ntheta3), 
     4     bsupu3(ns,nzeta,ntheta3), jsups3(ns,nzeta,ntheta3), 
     5     jsupv3(ns,nzeta,ntheta3), jsupu3(ns,nzeta,ntheta3),
     6     jdotb_sqrtg(ns,nzeta,ntheta3), sqrtg3(ns,nzeta,ntheta3), 
     7     phin(ns), toroidal_angle(nzeta), stat=j)

      bsubs3=0; bsubv3=0; bsubu3=0; jxb_gradp=0
      jcrossb=0 ; bsupv3=0; bsupu3=0; jsups3=0
      jsupv3=0; jsupu3=0; phin=0; phin(ns)=1
      jdotb_sqrtg=0; sqrtg3=0 
#endif
 
      ALLOCATE (pprime(nznt), pprim(ns),stat=j)
      pprim=0

      avforce=0; aminfor=0; amaxfor=0
      dnorm1 = twopi*twopi

      DO js = 2, ns1
         ovp = two/(vp(js+1) + vp(js))/dnorm1
         tjnorm = ovp*signgs
         sqgb2(:nznt) = sigma_an(js+1,:nznt)*gsqrt(js+1,:nznt)*
     1                  (bsq(js+1,:nznt)- pres(js+1))
     2                + sigma_an(js,:nznt)*gsqrt(js,:nznt)    *
     3                  (bsq(js,:nznt) - pres(js))
!        TAKE THIS OUT: MAY BE POORLY CONVERGED AT THIS POINT....
!         IF (ANY(sqgb2(:nznt)*signgs .le. zero)) 
!     1       STOP ' SQGB2 <= 0 in JXBFORCE'
         pprime(:) = ohs*(pres(js+1)-pres(js))/mu0              !dp/ds here
#ifdef _ANIMEC
!WAC  Last two lines of 'pprime' contain the hot particle parallel pressure
     1 + ohs*((pres(js+1)*phot(js+1) - pres(js)*phot(js))*pp2(js,:nznt) 
     2 +      (tpotb(js+1)           - tpotb(js)      )  *pp1(js,:nznt))
     3  / mu0
#endif
         kperpu(:nznt) = p5*(bsubv(js+1,:nznt,0) + bsubv(js,:nznt,0))*
     1                       pprime(:)/sqgb2
         kperpv(:nznt) =-p5*(bsubu(js+1,:nznt,0) + bsubu(js,:nznt,0))*
     1                       pprime(:)/sqgb2
         kp2(:nznt)=p5*(kperpu**2*(guu(js+1:nrzt:ns) + guu(js:nrzt:ns))
     1          + 2*kperpu*kperpv*(guv(js+1:nrzt:ns) + guv(js:nrzt:ns))
     2          +       kperpv**2*(gvv(js+1:nrzt:ns) + gvv(js:nrzt:ns)))
         itheta(js,:nznt) =  bsubsv(js,:nznt,0) - ohs*
     1                      (bsubv(js+1,:nznt,0) - bsubv(js,:nznt,0))
         izeta(js,:nznt)  = -bsubsu(js,:nznt,0) + ohs*
     1                      (bsubu(js+1,:nznt,0) - bsubu(js,:nznt,0))
         itheta(js,:nznt) = itheta(js,:nznt)/mu0
         izeta(js,:nznt)  = izeta(js,:nznt)/mu0
         sqrtg(:) = p5*(gsqrt(js,:) + gsqrt(js+1,:))
         bsupu1(:nznt) = p5*(bsupu(js+1,:nznt)*gsqrt(js+1,:)
     1                 +     bsupu(js,:nznt)  *gsqrt(js,:)) / sqrtg(:)
         bsupv1(:nznt) = p5*(bsupv(js+1,:nznt)*gsqrt(js+1,:)
     1                 +     bsupv(js,:nznt)  *gsqrt(js,:)) / sqrtg(:)
         bsubu1(:nznt) = p5*(bsubu(js+1,:nznt,0) + bsubu(js,:nznt,0))
         bsubv1(:nznt) = p5*(bsubv(js+1,:nznt,0) + bsubv(js,:nznt,0))
         jxb(:nznt) = ovp*(itheta(js,:nznt) * bsupv1(:nznt)
     1              -      izeta (js,:nznt) * bsupu1(:nznt))
         bdotk(js,:nznt) = itheta(js,:nznt) * bsubu1(:nznt) +
     1                     izeta (js,:nznt) * bsubv1(:nznt)
         pprime(:nznt) = ovp*pprime(:nznt)
         pnorm = one/(ABS(pprime(1)) + EPSILON(pprime(1)))
         amaxfor(js) = MAXVAL(jxb(:nznt)-pprime(:))*pnorm
         aminfor(js) = MINVAL(jxb(:nznt)-pprime(:))*pnorm
         avforce(js) = SUM(wint(2:nrzt:ns)*(jxb(:nznt) - pprime(:)))
         amaxfor(js) = 100*MIN(amaxfor(js),9.999_dp)
         aminfor(js) = 100*MAX(aminfor(js),-9.999_dp)
         pprim(js) = SUM(wint(js:nrzt:ns)*pprime(:))
!        Compute <K dot B>, <B sup v> = signgs*phip
!        jpar2 = <j||**2>, jperp2 = <j-perp**2>,  with <...> = flux surface average

         jdotb(js) = dnorm1*tjnorm*SUM(bdotk(js,:nznt)*wint(2:nrzt:ns)
     1                               / sigma_an(js,:nznt))
         bdotb(js) = dnorm1*tjnorm*SUM(sqgb2(:nznt)*wint(2:nrzt:ns)
     1                               / sigma_an(js,:nznt))
              
         bdotgradv(js) = p5*dnorm1*tjnorm*(phip(js) + phip(js+1))
         jpar2(js) = dnorm1*tjnorm*
     1            SUM(bdotk(js,:nznt)**2*wint(2:nrzt:ns)
     2              /(sigma_an(js,:nznt)*sqgb2(:nznt)))
         jperp2(js)= dnorm1*tjnorm*
     1            SUM(kp2(:nznt)*wint(2:nrzt:ns)*sqrtg(:nznt))

         IF (MOD(js,ns_skip) .eq. 0 .and. lprint_flag) THEN
#ifdef NETCDF
            phin(js) = phi(js)/phi(ns)
            DO lz = 1, nzeta
               toroidal_angle(lz)=REAL(360*(lz-1),rprec)/nzeta
               DO lt = 1, ntheta3
                  lk = lz + nzeta*(lt-1)
C                 lu (js,lz,lt ) =  lt
                  jsupu3 (js,lz,lt) = ovp*itheta(js,lk)
                  jsupv3 (js,lz,lt) = ovp*izeta(js,lk) 
                  jsups3 (js,lz,lt) = ovp*(bsubuv(js,lk)
     1					          -      bsubvu(js,lk))/mu0
                  bsupu3 (js,lz,lt) = bsupu1(lk)
                  bsupv3 (js,lz,lt) = bsupv1(lk)
                  jcrossb (js,lz,lt) = jxb(lk)
                  jxb_gradp (js,lz,lt) = (jxb(lk) - pprime(lk))
                  jdotb_sqrtg (js,lz,lt) = ovp*bdotk(js,lk) 
                  sqrtg3(js,lz,lt) = sqrtg(lk)*ovp
                  bsubu3(js,lz,lt) = bsubu(js,lk,0)
                  bsubv3(js,lz,lt) = bsubv(js,lk,0)
                  bsubs3(js,lz,lt) = bsubs(js,lk)
               END DO
            END DO
#else
            WRITE (njxbout, 200) phi(js), avforce(js), jdotb(js),
     1         bdotgradv(js), pprime(1), one/ovp, 
     2         (twopi**2)*tjnorm*SUM(itheta(js,:)*wint(js:nrzt:ns)),
     3         (twopi**2)*tjnorm*SUM(izeta (js,:)*wint(js:nrzt:ns)),
     4         amaxfor(js), aminfor(js)
            WRITE (njxbout, 90)
            DO lz = 1, nzeta, nv_skip
               WRITE (njxbout, 100) REAL(360*(lz-1),rprec)/nzeta, lz
               DO lt = 1, ntheta3, nu_skip
                  lk = lz + nzeta*(lt - 1)
                  WRITE (njxbout, 110) lt, tjnorm*itheta(js,lk),
     1              tjnorm*izeta(js,lk), ovp*(bsubuv(js,lk) -
     2              bsubvu(js,lk))/mu0, bsupu1(lk), bsupv1(lk),
     3              sqrtg(lk)*ovp, jxb(lk), jxb(lk) - pprime(lk),
     4              ovp*bdotk(js,lk), bsubu(js,lk,0),
     5              bsubv(js,lk,0), bsubs(js,lk)
               END DO
            END DO
#endif
         ENDIF
      END DO

      izeta(1,:nznt) = two*izeta(2,:nznt) - izeta(3,:nznt)           !!Need in wrout
      izeta(ns,:nznt)= two*izeta(ns-1,:nznt) - izeta(ns-2,:nznt)     !!Need in wrout
      jdotb(1) = two*jdotb(2) - jdotb(3)
      jdotb(ns) = two*jdotb(ns-1) - jdotb(ns-2)
      bdotb(1) = two*bdotb(3) - bdotb(2)
      bdotb(ns) = two*bdotb(ns-1) - bdotb(ns-2)
      bdotgradv(1) = two*bdotgradv(2) - bdotgradv(3)
      bdotgradv(ns) = two*bdotgradv(ns-1) - bdotgradv(ns-2)
      jpar2(1)   = 0; jpar2(ns)  = 0; jperp2(1)  = 0; jperp2(ns) = 0
      pprim(1) = 2*pprim(ns-1) - pprim(ns-2)
      pprim(ns) = 2*pprim(ns-1) - pprim(ns-2)

      IF (lprint_flag) THEN
#ifdef NETCDF
      CALL cdf_define(njxbout,vn_legend,legend)
      CALL cdf_define(njxbout,vn_mpol,mpol)
      CALL cdf_define(njxbout,vn_ntor,ntor)
      CALL cdf_define(njxbout,vn_phin,phin)
      CALL cdf_define(njxbout,vn_radial_surfaces,ns)
      CALL cdf_define(njxbout,vn_poloidal_grid_points,ntheta3)
      CALL cdf_define(njxbout,vn_toroidal_grid_points,nzeta)
      CALL cdf_define(njxbout,vn_avforce,avforce)
      CALL cdf_define(njxbout,vn_jdotb,jdotb)
 
      CALL cdf_define(njxbout,vn_sqg_bdotk,jdotb_sqrtg)
      CALL cdf_define(njxbout,vn_sqrtg,sqrtg3)

      CALL cdf_define(njxbout,vn_bdotgradv,bdotgradv)
      CALL cdf_define(njxbout,vn_pprime,pprim)
      CALL cdf_define(njxbout,vn_aminfor,aminfor)
      CALL cdf_define(njxbout,vn_amaxfor,amaxfor)
      CALL cdf_define(njxbout,vn_jsupu,jsupu3)
      CALL cdf_define(njxbout,vn_jsupv,jsupv3)
      CALL cdf_define(njxbout,vn_jsups,jsups3)
      CALL cdf_define(njxbout,vn_bsupu,bsupu3)
      CALL cdf_define(njxbout,vn_bsupv,bsupv3)
      CALL cdf_define(njxbout,vn_jcrossb,jcrossb)
      CALL cdf_define(njxbout,vn_jxb_gradp,jxb_gradp)
      CALL cdf_define(njxbout,vn_bsubu,bsubu3)
      CALL cdf_define(njxbout,vn_bsubv,bsubv3)
      CALL cdf_define(njxbout,vn_bsubs,bsubs3)

      CALL cdf_write(njxbout,vn_legend,legend)
      CALL cdf_write(njxbout,vn_mpol,mpol)
      CALL cdf_write(njxbout,vn_ntor,ntor)
      CALL cdf_write(njxbout,vn_phin,phin)
      CALL cdf_write(njxbout,vn_radial_surfaces,ns)
      CALL cdf_write(njxbout,vn_poloidal_grid_points,ntheta3)
      CALL cdf_write(njxbout,vn_toroidal_grid_points,nzeta)
      CALL cdf_write(njxbout,vn_avforce,avforce)
      CALL cdf_write(njxbout,vn_jdotb,jdotb)

      CALL cdf_write(njxbout,vn_sqg_bdotk,jdotb_sqrtg)
      CALL cdf_write(njxbout,vn_sqrtg,sqrtg3)

      CALL cdf_write(njxbout,vn_bdotgradv,bdotgradv)
      CALL cdf_write(njxbout,vn_pprime,pprim)
      CALL cdf_write(njxbout,vn_aminfor,aminfor)
      CALL cdf_write(njxbout,vn_amaxfor,amaxfor)
      CALL cdf_write(njxbout,vn_jsupu,jsupu3)
      CALL cdf_write(njxbout,vn_jsupv,jsupv3)
      CALL cdf_write(njxbout,vn_jsups,jsups3)
      CALL cdf_write(njxbout,vn_bsupu,bsupu3)
      CALL cdf_write(njxbout,vn_bsupv,bsupv3)
      CALL cdf_write(njxbout,vn_jcrossb,jcrossb)
      CALL cdf_write(njxbout,vn_jxb_gradp,jxb_gradp)
      CALL cdf_write(njxbout,vn_bsubu,bsubu3)
      CALL cdf_write(njxbout,vn_bsubv,bsubv3)
      CALL cdf_write(njxbout,vn_bsubs,bsubs3)
 
      CALL cdf_close(njxbout)

      DEALLOCATE(
     1     bsubs3, bsubv3, bsubu3, jxb_gradp, jcrossb, bsupv3, 
     2     bsupu3, jsups3, jsupv3, jsupu3, jdotb_sqrtg, phin, 
     3     toroidal_angle, sqrtg3, stat=j)

#else
      CLOSE (njxbout)

   90 FORMAT(/"   LU      JSUPU      JSUPV      JSUPS      BSUPU",
     1   "      BSUPV   SQRT(g')     J X B   J X B - p'     J * B",
     2   "      BSUBU      BSUBV      BSUBS   "/)
  100 FORMAT( " TOROIDAL ANGLE (PER PERIOD) = ", f8.3," DEGREES",
     1        " (PLANE #", i3,")")
  110 FORMAT(i5,1p,12e11.3)
  200 FORMAT(/" TOROIDAL FLUX =  ",1p,e12.3,3x,"<J X B - p'> = ",
     1   e12.3,3x,"<J DOT B> = ",e12.3,3x,
     2   "<B DOT GRAD(V)> = ",e12.3,/,
     2   " dp/d(VOL) [p'] = ",e12.3,3x,'d(VOL)/ds    = ',e12.3,3x,
     2   "<JSUPU>   = ",e12.3,3x,"<JSUPV>         = ",e12.3,/,
     3   " MAXIMUM FORCE DEVIATIONS (RELATIVE TO p'): ",sp,0p,f7.2,"%",
     4     3x,f7.2,"%")
#endif

      END IF
      
      DEALLOCATE (kperpu, kperpv, sqgb2, sqrtg, kp2, brhomn, bsubsmn, 
     1    jxb, jxb2, bsupu1, bsupv1, bsubu1, bsubv1, avforce, aminfor, 
     2    amaxfor, pprim, stat=j)
!
!     COMPUTE MERCIER CRITERION
!
      bdotk = mu0*bdotk
      CALL Mercier(gsqrt,bsq,bdotk,iotas,wint,r1,ru,rv,zu,zv,bsubu,
     1             vp,phips,pres,ns,nznt)

      DEALLOCATE (bdotk, bsubuv, bsubvu, pprime, stat=j)

      END SUBROUTINE jxbforce
