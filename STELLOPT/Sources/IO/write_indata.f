      SUBROUTINE write_indata(input_file, istat)
      USE optim, input_file_opt=>input_file
      USE vmec_input
      USE bootsj_input
      USE safe_open_mod
      USE write_array_generic
      USE vacfield_mod, ONLY: write_invac
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(out) :: istat
      CHARACTER(LEN=*) :: input_file
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, DIMENSION(1) :: ins
      INTEGER :: iftol, m, n, i, k, ns0, ns1, iunit,
     1    nextcur_max, max_k
      character(len=8) :: form
      character(len=100) :: temp_write
      character*(*), parameter :: form1 = "(a,i1,a)",
     1    form2 = "(a,i2,a)"
      CHARACTER(len=*), DIMENSION(21), PARAMETER :: comments = (/
     &  'OPTIMIZER (AND MPI) RUN CONTROL PARAMETERS  ',
     &  'LOGICAL VARIABLES: ACCESS TO PHYSICS MODULES',
     &  'SIZE AND FIELD STRENGTH SCALE PARAMETERS    ',
     &  'BOOZER COORDINATE TRANSFORMATION PARAMETERS ',
     &  'SCALAR OPTIMIZATION PARAMETERS              ',
     &  'TRANSPORT OPTIMIZATION PARAMETERS           ',
     &  'BALLOONING MODE OPTIMIZATION PARAMETERS     ',
     &  'BOOTSTRAP CURRENT OPTIMIZATION PARAMETERS   ',
     &  'DKES AND NEO OPTIMIZATION PARAMETERS        ',
     &  'DIAGNOSTICS PARAMETERS  (V3FIT)             ',
     &  'PRESSURE PROFILE PARAMETERS                 ',
     &  'MOTIONAL STARK EFFECT (MSE) DIAGNOSTIC      ',
     &  'MAGNETIC DIAGNOSTICS (DIAGNO)               ',
     &  'LIMITER PARAMETERS                          ',
     &  'VARIABLES ANIMEC/VMEC                       ',
     &  'ELECTRON DENSTITY MEASUREMENTS              ',
     &  'ELECTRON TEMPERATURE MEASUREMENTS           ',
     &  'ION TEMPERATURE MEASUREMENTS                ',
     &  'ELECTRON DENSITY PROFILE                    ',
     &  'ELECTRON TEMPERATURE PROFILE                ',
     &  'ION TEMPERATURE PROFILE                     '
     &   /)
C-----------------------------------------------
      iunit = iout+2
      CALL safe_open(iunit, istat, input_file, 'replace', 'formatted')
      IF (istat .ne. 0) THEN
         istat = -3
         RETURN
      END IF

      ins = MAXLOC(ns_array)
      ns0 = ns_array(ins(1))

      iftol = 1
      DO WHILE(ftol_array(iftol).ne.zero .and. iftol.lt.100)
         iftol = iftol + 1
      END DO

!
!     Compute nextcur_vmec, which could have changed if lcoil_geom = .true.
!
      IF (nextcur_vmec .eq. 0) THEN
         DO nextcur_vmec = SIZE(extcur), 1, -1
            IF (extcur(nextcur_vmec) .ne. zero) EXIT
         END DO
         nextcur_vmec = MAX(0, nextcur_vmec)
      END IF

      DO nextcur_max = SIZE(lextcur), 1, -1
         IF (lextcur(nextcur_max)) EXIT
      END DO
      nextcur_max = MAX(0, nextcur_max)
!
!     Start writing out formatted VMEC data with &INDATA heading
!     so that it can be read in under the indata NAMELIST
!
      WRITE (iunit, '(a7)') '&INDATA'
      WRITE (iunit, '(2x,3a)') "MGRID_FILE = '", TRIM(mgrid_file), "'"
      WRITE (iunit, 1007) 'LFREEB = ', lfreeb
      WRITE (iunit, 1007) 'LFORBAL = ', lforbal
      WRITE (iunit, 1007) 'LASYM = ', lasym
      WRITE (iunit, '(2x,a,1p,e10.2)') 'DELT = ', delt
      WRITE (iunit, '(2x,a,1p,e10.2)') 'TCON0 = ', tcon0
      WRITE (iunit,'(2x,3a)') "PRECON_TYPE = '", TRIM(precon_type),"'" 
      WRITE (iunit,'(2x,a,1p,e14.6)') "PREC2D_THRESHOLD = ", 
     1                                prec2d_threshold
      WRITE (iunit, '(2x,a,i4)') 'NFP = ', nfp
      WRITE (iunit, '(2(2x,a,i4))') 'MPOL = ', mpol, 'NTOR = ', ntor
      IF (ntheta.gt.0) WRITE (iunit, '(2x,a,i4)') 'NTHETA = ',ntheta
      IF (nzeta.gt.0)  WRITE (iunit, '(2x,a,i4)') 'NZETA  = ',nzeta
      IF (mfilter_fbdy.gt.0)
     1    WRITE (iunit, '(2x,a,i4)') 'mfilter_fbdy  = ',mfilter_fbdy
      IF (nfilter_fbdy.gt.0)  
     1    WRITE (iunit, '(2x,a,i4)') 'nfilter_fbdy  = ',nfilter_fbdy

      WRITE (iunit, '(2x,a)') 'NS_ARRAY = '
      WRITE (iunit, 980, err=2000) (ns_array(i),i=1,ins(1))
      WRITE (iunit, '(2x,a,i6)') 'NITER = ', niter
      WRITE (iunit, '(2x,a,i4)') 'NSTEP = ', nstep
      WRITE (iunit, '(2x,a,i4)') 'NVACSKIP = ', nvacskip
      WRITE (iunit, '(2x,a)') 'FTOL_ARRAY ='
      WRITE (iunit, 995, err=2000) (ftol_array(i),i=1,iftol - 1)
      WRITE (iunit, 996) phiedge
      IF (nextcur_vmec .gt. 0)
     1    CALL write_array(iunit, 'EXTCUR', extcur, nextcur_vmec)

      CALL write_rbzb (iunit, istat)
      IF (istat .ne. 0) GOTO 2000

  980 FORMAT(2x,16i5)
  981 FORMAT(2x,a,10(i5))
  990 FORMAT(2x,a,6(1p,e22.14))
  991 FORMAT(1x,a,5(1p,e122.14))
  993 FORMAT(2(2x,a,1p,e22.14))
  994 FORMAT(2x,'GAMMA = ',1p,e14.6)
  995 FORMAT(2x,1p,4e14.6)
  996 FORMAT(2x,'PHIEDGE = ',1p,e22.14)
 1000 FORMAT(5(1p,e22.14))
 1007 FORMAT(2x,a,L1)
 1008 FORMAT(2x,a,L1,2x,a,i2)
!***********************************************************************
! uncomment the next line to create write_indata_only
!      RETURN
!***********************************************************************
!
!     Finish by writing out formatted optimization data with &OPTIMUM
!     so that it can be read in under the optimum NAMELIST
!     Option to add comments (starting with '!') now available
!
      WRITE (iunit, '(a)') '/'
      WRITE (iunit, '(a8)') '&OPTIMUM'
      WRITE (iunit, 100) TRIM(comments(1))
      WRITE (iunit, 990) 'EPSFCN = ', epsfcn
      WRITE (iunit, '(2x,a,i5)') 'NITER_OPT = ', niter_opt
      WRITE (iunit, '(2x,a,i5)') 'NUM_PROCESSORS = ', num_processors
      WRITE (iunit, '(2x,a,i5)') 'NUM_LEVMAR_PARAMS = ',
     1   num_levmar_params
      WRITE (iunit, '(2x,a)') 'NSURF_MASK = '
      ns1 = MIN(ns0, SIZE(nsurf_mask))
      WRITE (iunit, '(10f7.3)', err=2000) (nsurf_mask(i),i=1,ns1)

      WRITE (iunit, '(2x,a,i5)') 'NOPT_ALG = ', nopt_alg
      WRITE (iunit, '(2x,a,i5)') 'NOPT_BOUNDARY = ', nopt_boundary
      WRITE (iunit, 1007) 'LRESET_OPT = ', lreset_opt
      WRITE (iunit, 1007) 'LDIAG_OPT = ', ldiag_opt
      WRITE (iunit, 1007) 'LKEEP_MINS = ', lkeep_mins
c ---------------------------------------------------------------------------
c  SAL 2011 - ANIMEC BASED OPTIMIZATIONS
      WRITE (iunit, 100) TRIM(comments(15))
      WRITE (iunit, 1007) 'LANIMEC = ', lanimec
      IF (lanimec) THEN
         WRITE (iunit, 1007) 'LANI_BCRIT = ', lani_bcrit                   
         WRITE (iunit, 1007) 'LANI_PHOT = ', lani_PHOT                     
         CALL write_array(iunit, 'AH_MASK', ah_mask, 
     1                       size(ah_mask),low_index=0)                
         WRITE (iunit, 1007) 'LANI_TPERP = ', lani_TPERP
         CALL write_array(iunit, 'AT_MASK', at_mask, 
     1                       size(at_mask),low_index=0)
      END IF
      WRITE (iunit, 1007) 'LPHIEDGE = ', lphiedge
      WRITE (iunit, 1007) 'PHIEDGE_DIODE = ',phiedge_diode
      WRITE (iunit, 1007) 'LIOTA_PROF_OPT = ', liota_prof_opt
      WRITE (iunit, 1007) 'LCUR_PROF_OPT = ', lcur_prof_opt
      write (iunit, 1007) 'LCUR_OPT_EDGE0 = ', lcur_opt_edge0
      WRITE (iunit, 1007) 'LEDGE_CURRENT = ', ledge_current
      write (iunit, '(2x,a,11i2)') 'AC_MASK = ',ac_mask                     ! SAL
      WRITE (iunit, 1007) 'LPRES_PROF_OPT = ', lpres_prof_opt
      write (iunit, 1007) 'LPRES_PROF_FIT = ', lpres_prof_fit
      write (iunit, 1007) 'LPRES_OPT_EDGE0 = ', lpres_opt_edge0
      write (iunit, 1007) 'LPRES_OPT_EDGEGR0 = ', lpres_opt_edgegr0
      WRITE (iunit, 1007) 'LRECONP = ', lreconp
      WRITE (iunit, 1007) 'LRECONJ = ', lreconj
      WRITE (iunit, 1008) 'LP1ZERO = ',lp1zero,' KPP = ',kpp
      WRITE (iunit, 1008) 'LJ1ZERO = ',lj1zero,' KJJ = ',kjj
      WRITE (iunit, 100) TRIM(comments(2))
      write (iunit, 1007) 'LBOUNDARY = ', lboundary
      WRITE (iunit, 1007) 'LBMN = ', lbmn
      WRITE (iunit, 1007) 'LJ_STAR = ', lj_star
      WRITE (iunit, 1007) 'LASPECT_MAX = ', laspect_max
      WRITE (iunit, 1007) 'LBETA_MIN = ', lbeta_min
      WRITE (iunit, 1007) 'LJ_INVARIANT = ', lj_invariant
      WRITE (iunit, 1007) 'LBOOTSJ_OPT = ', lbootsj_opt
      WRITE (iunit, 1007) 'LKINK_OPT = ', lkink_opt
      WRITE (iunit, 1007) 'LBALLOON_OPT = ', lballoon_opt
      WRITE (iunit, 1007) 'L_LEGENDRE = ', l_legendre                 !! LEGENDRE
      WRITE (iunit, 1007) 'LDKES_OPT = ', ldkes_opt                   !! RHF
      WRITE (iunit, 1007) 'LNEO_OPT = ', lneo_opt
      WRITE (iunit, 1007) 'LDSUBR_OPT = ', ldsubr_opt
      WRITE (iunit, 1007) 'LORBIT_OPT = ', lorbit_opt
      write (iunit, 1007) 'LDIAGNO_OPT = ', ldiagno_opt
      WRITE (iunit, 1007) 'LNESCOIL_OPT = ', lnescoil_opt
      WRITE (iunit, 1007) 'LCOIL_GEOM = ', lcoil_geom
      write (iunit, 990) 'SIGMA_BERR_AVG = ', sigma_berr_avg
      write (iunit, 990) 'SIGMA_BERR_MAX = ', sigma_berr_max
      WRITE (iunit, 1007) 'LV3POST = ', lv3post
      WRITE (iunit, 1007) 'LVAC_OPT = ', lvac_opt

      i = LBOUND(lfix_ntor,1)
      CALL write_array(iunit, 'LFIX_NTOR', lfix_ntor, SIZE(lfix_ntor),
     1                 low_index=i)
      DO m = 0, mpol1d
         IF (ANY(lfix_rhob(:,m)))
     1      CALL write_array(iunit, 'LFIX_RHOB',
     2           lfix_rhob(-ntord:ntord,m), SIZE(lfix_rhob,1), ndim2=m,
     1           low_index=-ntord)
      END DO

      IF (nextcur_opt .gt. 0) THEN
         CALL write_array(iunit, 'LEXTCUR',
     1                    lextcur, nextcur_vmec)
!         CALL write_array(iunit, 'LEXTCUR',
!     1                    lextcur, nextcur_max)
         CALL write_array(iunit, 'TARGET_EXTCUR',
     1                    target_extcur, nextcur_vmec)
         CALL write_array(iunit, 'SIGMA_EXTCUR',
     1                    sigma_extcur, nextcur_vmec)
         CALL write_array(iunit, 'OH_COEFS',
     1                    oh_coefs, nextcur_vmec)
         WRITE (iunit, 990) 'SIGMA_OH = ', sigma_oh
      ENDIF



      WRITE (iunit, 100) TRIM(comments(3))
      WRITE (iunit, '(3(2x,a,1p,e14.6))')
     1    'R00_OPT = ', r00_opt, 'R00_SCALE = ', r00_scale, 
     2    'B00_SCALE = ', b00_scale

      WRITE (iunit, '(2(2x,a,1p,e14.6))')
     1    'RGRID_MIN = ', rgrid_min, 'RGRID_MAX = ', rgrid_max
      WRITE (iunit, '(2(2x,a,1p,e14.6))')
     1    'ZGRID_MIN = ', zgrid_min, 'ZGRID_MAX = ', zgrid_max
      WRITE (iunit, 100) TRIM(comments(4))
      WRITE (iunit, '(2(2x,a,i4))') 'MBOZ_OPT = ', mboz_opt,
     1                              'NBOZ_OPT = ',nboz_opt

      WRITE (iunit, 100) TRIM(comments(5))
      WRITE (iunit, '(2x,a,1pe14.6)') 'PHIEDGE_MIN = ', phiedge_min
      WRITE (iunit, '(2x,a,1pe14.6)') 'PHIEDGE_MAX = ', phiedge_max
      WRITE (iunit, 990) 'COIL_SEPARATION = ', coil_separation
      WRITE (iunit, 993) 'TARGET_ASPECTRATIO = ', target_aspectratio,
     1                   'SIGMA_ASPECT = ', sigma_aspect
      WRITE (iunit, 993) 'TARGET_BETA = ', target_beta,
     1                   'SIGMA_BETA = ', sigma_beta
      write (iunit, 993) 'TARGET_CURTOR = ', target_curtor,
     1                   'SIGMA_CURTOR = ', sigma_curtor
      write (iunit, 993) 'TARGET_EPLASMA = ', target_eplasma,
     1                   'SIGMA_EPLASMA = ', sigma_eplasma
      write (iunit, 993) 'TARGET_JEDGE = ', target_jedge,
     1                   'SIGMA_JEDGE = ', sigma_jedge
      CALL write_array(iunit, 'TARGET_KINK', target_kink,
     1                 SIZE(target_kink))
      CALL write_array(iunit, 'SIGMA_KINK',sigma_kink,SIZE(sigma_kink))
      WRITE (iunit, 993) 'TARGET_MAXCURRENT = ', target_maxcurrent,
     1                   'SIGMA_MAXCURRENT = ', sigma_maxcurrent
      WRITE (iunit, 993) 'TARGET_RMAX = ', target_rmax,
     1                   'SIGMA_RMAX = ', sigma_rmax
      WRITE (iunit, 993) 'TARGET_RMIN = ', target_rmin,
     1                   'SIGMA_RMIN = ', sigma_rmin
      WRITE (iunit, 993) 'TARGET_ZMAX = ', target_zmax,
     1                   'SIGMA_ZMAX = ', sigma_zmax
      WRITE (iunit, 993) 'TARGET_KAPPA = ', target_kappa,
     1                   'SIGMA_KAPPA = ', sigma_kappa
      WRITE (iunit, 993) 'TARGET_ELLIPTICITY = ', target_ellipticity,
     1                   'SIGMA_ELLIPTICITY = ', sigma_ellipticity
      WRITE (iunit, 993) 'TARGET_FLUXP = ', target_fluxp,
     1                   'SIGMA_FLUXP = ', sigma_fluxp
      WRITE (iunit, 993) 'TARGET_RBTOR = ', target_RBtor,
     1                   'SIGMA_RBTOR = ', sigma_RBtor
      WRITE (iunit, 993) 'TARGET_COIL_COMPLEX = ', target_coil_complex,
     1                   'SIGMA_COIL_COMPLEX = ', sigma_coil_complex
      WRITE (iunit, 993) 'TARGET_COIL_JMAX = ', target_coil_jmax,
     1                   'SIGMA_COIL_JMAX = ', sigma_coil_jmax
      CALL write_array(iunit, 'TARGET_IOTA', target_iota, 
     1                 SIZE(target_iota), low_index=0)
      ns1 = MIN(ns0, SIZE(sigma_iota))
      CALL write_array(iunit, 'SIGMA_IOTA', sigma_iota, ns1)
      CALL write_array(iunit, 'TARGET_IOTA_P', target_iota_p, 
     1                 SIZE(target_iota_p), low_index=0)
      ns1 = MIN(ns0, SIZE(sigma_iota_pmax))
      CALL write_array(iunit,'SIGMA_IOTA_PMAX',sigma_iota_pmax, ns1)
      CALL write_array(iunit,'SIGMA_IOTA_PMIN',sigma_iota_pmin, ns1)
      WRITE (iunit, 993) 'TARGET_IOTA_MIN = ', target_iota_min,
     1                   'SIGMA_IOTA_MIN = ', sigma_iota_min
      WRITE (iunit, 993) 'TARGET_IOTA_MAX = ', target_iota_max,
     1                   'SIGMA_IOTA_MAX = ', sigma_iota_max
      write (iunit, 993) 'TARGET_IOTA_MAX_MIN = ', target_iota_max_min,
     1                   'SIGMA_IOTA_MAX_MIN = ', sigma_iota_max_min
      call write_array(iunit, 'TARGET_WELL', target_well,
     1                 size(target_well), low_index=0)

      WRITE (iunit, 990) 'SIGMA_CURV = ', sigma_curv
      WRITE (iunit, 990) 'SIGMA_BERR_AVE = ', sigma_berr_ave
C      write (iunit, 990) 'SIGMA_PSEUDO = ', sigma_pseudo
C      write (iunit, 990) 'SIGMA_PSEUDO2 = ', sigma_pseudo2
      call write_array(iunit, 'SIGMA_PSEUDO', sigma_pseudo, ns0)
      call write_array(iunit, 'SIGMA_PSEUDO2', sigma_pseudo2, ns0)
      write (iunit, 1007) 'LPSEUDO_SIN = ', lpseudo_sin
      ns1 = MIN(ns0, SIZE(sigma_mercier))
      CALL write_array(iunit, 'SIGMA_MERCIER', sigma_mercier, ns1)
      ns1 = MIN(ns0, SIZE(sigma_vp))
      call write_array(iunit, 'SIGMA_VP', sigma_vp, ns1)
      CALL write_array(iunit, 'N_JAC', n_jac, SIZE(n_jac))
      CALL write_array(iunit, 'M_JAC', m_jac, SIZE(m_jac))
      CALL write_array(iunit, 'SIGMA_JAC', sigma_jac, SIZE(sigma_jac))
      CALL write_array(iunit, 'N_VAC_ISLAND', n_vac_island, 
     1                 SIZE(n_vac_island))
      CALL write_array(iunit, 'M_VAC_ISLAND', m_vac_island, 
     1                 SIZE(m_vac_island))
      CALL write_array(iunit, 'SIGMA_VAC_ISLAND', sigma_vac_island,
     1                 SIZE(sigma_vac_island))

      WRITE (iunit, 100) TRIM(comments(6))
      WRITE (iunit, '(2x,a12,i3,a2,i3,a1)') 'HELICITY = (',
     1  NINT(REAL(helicity)),', ',NINT(AIMAG(helicity)),')'

      ns1 = MIN(ns0, SIZE(sigma_bmin))
      CALL write_array(iunit,'SIGMA_BMIN', sigma_bmin, ns1)
      CALL write_array(iunit,'SIGMA_BMAX', sigma_bmax, ns1)
c      CALL write_array(iunit,'SIGMA_BMN',sigma_bmn, ns1)
      CALL write_array(iunit,'SIGMA_BMN',sigma_bmn, ns1)
      IF( any(sigma_bmn_tgt .lt. bigno)) THEN
         max_k = 0
         DO k = 1, size(sigma_bmn_tgt,2)
            IF( all(sigma_bmn_tgt(1:ns0,k).ge. bigno)) CYCLE
            IF (k .lt. 10) THEN
               form = form1
               WRITE(iunit,'(a,i1,a,i1,a,i1,a,i1)')
     1                     '  N_BMN_TGT(',k,') = ',n_bmn_tgt(k),
     2                     ',  M_BMN_TGT(',k,') = ',m_bmn_tgt(k)
            END IF
            IF (k .ge. 10) THEN
               form = form2
               WRITE(iunit,'(a,i2,a,i2,a,i2,a,i2)')
     1                     '  N_BMN_TGT(',k,') = ',n_bmn_tgt(k),
     2                     ', M_BMN_TGT(',k,') = ',m_bmn_tgt(k)
            END IF
            WRITE (temp_write, form, err=2000) 'SIGMA_BMN_TGT(1:,',
     1                                         k, ')'
!            CALL write_array(iunit, trim(temp_write),
!     1                       sigma_bmn_tgt(1,k), ns0)

            WRITE (temp_write, form, err=2000) 'TARGET_BMN_TGT(1:,',
     1                                         k, ')'
!            CALL write_array(iunit,trim(temp_write),
!     1                       target_bmn_tgt(1,k), ns0)
            max_k = k
         END DO
      END IF
      
      CALL write_array(iunit,'SIGMA_RIPPLE',sigma_ripple, ns1)

      WRITE (iunit, '(2(2x,a,i4))') 'NUMJSTAR = ',NumJstar,
     1   ' NUMJINVARIANT = ', NumJinvariant
      ns1 = MIN(ns0, SIZE(sigma_jstar,1))
      DO k = 1, MIN(NumJstar,SIZE(sigma_jstar,2))
         CALL write_array(iunit,'SIGMA_JSTAR',
     1                    sigma_jstar(1:ns1,k), ns1, k)
      END DO
      ns1 = MIN(ns0, SIZE(sigma_jinvariant,1))
      DO k = 1, MIN(NumJinvariant,SIZE(sigma_jinvariant,2))
         CALL write_array(iunit, 'SIGMA_JINVARIANT',
     1        sigma_jinvariant(1:ns1,k), ns1, k)
      END DO

      WRITE (iunit, '(2(2x,a,i4))') 'NS_JCONF_SRC = ',NS_JConf_Src,
     1   'NS_JCONF_TGT = ',  NS_JConf_Tgt
      WRITE(iunit,981) 'NPITCH_JCONF = ', NPitch_JConf
      WRITE(iunit,990) 'SIGMA_JCONF = ',sigma_jconf

!----------------------------------------------------------------------------------------
!     CODE added by R.SANCHEZ (01/19/99) to WRITE out ballooning-related info.
!     Modified (02/01/99) to include multiple initial position output.
!
      WRITE (iunit,100) TRIM(comments(7))
      IF(lballoon_opt .or. lpres_prof_opt) THEN
        ns1 = MIN(ns0, SIZE(nballoon_mask))
        CALL write_array(iunit,'NBALLOON_MASK',nballoon_mask, ns1)
        CALL write_array(iunit,'TARGET_BALLOON',target_balloon, ns1)
        CALL write_array(iunit,'SIGMA_BALLOON',sigma_balloon, ns1)
        CALL write_array(iunit,'SIGMA_PGRAD',sigma_pgrad, ns1)
        WRITE (iunit, 990) 'SIGMA_PEDGE = ', sigma_pedge
        CALL write_array(iunit,'BAL_ZETA0',bal_zeta0, nini_zeta)
        CALL write_array(iunit,'BAL_THETA0', bal_theta0, nini_theta)
      ENDIF
!-----------------------------------------------------------------------------------------

      WRITE (iunit, 100) TRIM(comments(8))
      IF (lbootsj_opt) THEN
         CALL write_array(iunit,'FBOOT',fboot, 11, low_index=0)
         IF( jboot .ne. 0) WRITE (iunit, '(2x,a,i1)') 'JBOOT = ',jboot
         IF( lseedcur )
     1      CALL write_array(iunit,'ASEEDCUR',aseedcur, 11, low_index=0)

         ns1 = MIN(ns0, SIZE(sigma_bootsj))
         CALL write_array(iunit,'SIGMA_BOOTSJ',sigma_bootsj, ns1)
      END IF

      WRITE (iunit, 100) TRIM(comments(9))                          !RHF
      ns1 = MIN(ns0, SIZE(ndkes_mask))
      CALL write_array(iunit, 'NDKES_MASK', ndkes_mask, ns1)        !RHF
      CALL write_array(iunit, 'DKES_NU', dkes_nu, ns1)              !RHF
      CALL write_array(iunit, 'DKES_EFIELD', dkes_efield, ns1)      !RHF
      CALL write_array(iunit,'SIGMA_DKES', sigma_dkes, ns1)         !RHF

c ---------------------------------------------------------------------------
!  NEO
      IF (lneo_opt) THEN
        ns1 = MIN(ns0, SIZE(nneo_mask))
        CALL write_array(iunit, 'NNEO_MASK', nneo_mask, ns1)
        CALL write_array(iunit, 'SIGMA_NEO', sigma_neo, ns1)
      END IF

      IF (ldsubr_opt) THEN
         ns1 = MIN(ns0, SIZE(sigma_dsubr))
         CALL write_array(iunit, 'SIGMA_DSUBR', sigma_dsubr, ns1)
      END IF

      IF (lorbit_opt)
     1   WRITE (iunit, 990) 'SIGMA_ORBIT = ', sigma_orbit

c ---------------------------------------------------------------------------
c  MCZ & SAL 2011 - Pressure Fitting
      IF( np_prof > 0 .and. any(sigma_p_prof > 0)) then
         WRITE (iunit, 100) TRIM(comments(11))
         WRITE (iunit,981) 'NP_PROF = ', NP_PROF
         WRITE (iunit, 990) 'FACTOR_P_PROF = ', factor_p_prof
         WRITE (iunit, 1007) 'ISOTE = ', isote
         WRITE (iunit, 1007) 'LP_PROF_INCL_EDGE = ', lp_prof_incl_edge
         WRITE (iunit,981) 'PRES_OPT_NMAX = ', pres_opt_nmax
         CALL write_array(iunit,'P_PROF', p_prof, np_prof)
         !CALL write_array(iunit,'NE_PROF', te_prof, np_prof)    
         !CALL write_array(iunit,'TE_PROF', ne_prof, np_prof)    
         CALL write_array(iunit,'SIGMA_P_PROF', sigma_p_prof, np_prof)
         CALL write_array(iunit,'R_P_PROF', r_p_prof, np_prof)    
         CALL write_array(iunit,'Z_P_PROF', z_p_prof, np_prof)
         CALL write_array(iunit,'PHI_P_PROF', phi_p_prof, np_prof)
      ELSE IF (     ANY(sigma_ne_prof > 0)
     1         .or. ANY(sigma_te_prof > 0)
     2         .or. ANY(sigma_ti_prof > 0) ) THEN
               WRITE (iunit, 100) TRIM(comments(11))
               WRITE (iunit, 990) 'FACTOR_P_PROF = ', factor_p_prof
               WRITE (iunit, 990) 'SIGMA_PEDGE = ', sigma_pedge
               CALL write_array(iunit,'SIGMA_PGRAD',sigma_pgrad, ns1)
      END IF
      IF (ANY(sigma_ne_prof > 0)) THEN
         WRITE (iunit, 100) TRIM(comments(16))
         WRITE (iunit,981) 'NNE_PROF = ', nne_prof
         DO m = 1, nne_prof
            WRITE(iunit,'(5(a,i3,a,1p,e13.6,3x))')
     1         '  R_NE_PROF(',m,') =',r_ne_prof(m),
     2         '  PHI_NE_PROF(',m,') =',phi_ne_prof(m),
     3         '  Z_NE_PROF(',m,') =',z_ne_prof(m),
     4         '  NE_PROF(',m,') =',ne_prof(m),
     5         '  SIGMA_NE_PROF(',m,') =',sigma_ne_prof(m)
         END DO
      END IF
      IF (ANY(sigma_te_prof > 0)) THEN
         WRITE (iunit, 100) TRIM(comments(17))
         WRITE (iunit,981) 'NTE_PROF = ', nte_prof
         DO m = 1, nte_prof
            WRITE(iunit,'(5(a,i3,a,1p,e13.6,3x))')
     1         '  R_TE_PROF(',m,') =',r_te_prof(m),
     2         '  PHI_TE_PROF(',m,') =',phi_te_prof(m),
     3         '  Z_TE_PROF(',m,') =',z_te_prof(m),
     4         '  TE_PROF(',m,') =',te_prof(m),
     5         '  SIGMA_TE_PROF(',m,') =',sigma_te_prof(m)
         END DO
      END IF
      IF (ANY(sigma_ti_prof > 0)) THEN
         WRITE (iunit, 100) TRIM(comments(18))
         WRITE (iunit,981) 'NTI_PROF = ', nti_prof
         DO m = 1, nti_prof
            WRITE(iunit,'(5(a,i3,a,1p,e13.6,3x))')
     1         '  R_TI_PROF(',m,') =',r_ti_prof(m),
     2         '  PHI_TI_PROF(',m,') =',phi_ti_prof(m),
     3         '  Z_TI_PROF(',m,') =',z_ti_prof(m),
     4         '  TI_PROF(',m,') =',ti_prof(m),
     5         '  SIGMA_Ti_PROF(',m,') =',sigma_ti_prof(m)
         END DO
      END IF
      
      i = minloc(ne_aux_s(2:),DIM=1)
      IF (i > 4) THEN
         WRITE (iunit, 100) TRIM(comments(19))
         WRITE (iunit, 101) '  NE_AUX_S = ',(ne_aux_s(n), n=1,i)
         WRITE (iunit, 101) '  NE_AUX_F = ',(ne_aux_f(n), n=1,i)
      END IF
      i = minloc(te_aux_s(2:),DIM=1)
      IF (i > 4) THEN
         WRITE (iunit, 100) TRIM(comments(20))
         WRITE (iunit, 101) '  TE_AUX_S = ',(te_aux_s(n), n=1,i)
         WRITE (iunit, 101) '  TE_AUX_F = ',(te_aux_f(n), n=1,i)
      END IF
      i = minloc(ti_aux_s(2:),DIM=1)
      IF (i > 4) THEN
         WRITE (iunit, 100) TRIM(comments(21))
         WRITE (iunit, 101) '  TI_AUX_S = ',(ti_aux_s(n), n=1,i)
         WRITE (iunit, 101) '  TI_AUX_F = ',(ti_aux_f(n), n=1,i)
      END IF
c ---------------------------------------------------------------------------
c  SAL 2011 - MSE Diagnostic
      IF(nmse_cams*nmse_chords > 0) THEN
         WRITE (iunit, 100) TRIM(comments(12))
         WRITE(iunit,981) 'NMSE_CAMS = ',nmse_cams
         WRITE(iunit,981) 'NMSE_CHORDS = ',nmse_chords
         IF (ANY(mse_a1_coef .ne. 0)) THEN
            DO m=1,nmse_cams
               DO n=1,nmse_chords
               WRITE(iunit,'(15(a,i3,a,i2,a,1p,e13.6,3x))')
     1             '  MSE_R(',m,',',n,') = ',mse_r(m,n),
     2             '  MSE_PHI(',m,',',n,') = ',mse_phi(m,n),
     3             '  MSE_Z(',m,',',n,') = ',mse_z(m,n),
     4             '  MSE_A1_COEF(',m,',',n,') = ',mse_a1_coef(m,n),
     4             '  MSE_A2_COEF(',m,',',n,') = ',mse_a2_coef(m,n),
     4             '  MSE_A3_COEF(',m,',',n,') = ',mse_a3_coef(m,n),
     4             '  MSE_A4_COEF(',m,',',n,') = ',mse_a4_coef(m,n),
     4             '  MSE_A5_COEF(',m,',',n,') = ',mse_a5_coef(m,n),
     4             '  MSE_A6_COEF(',m,',',n,') = ',mse_a6_coef(m,n),
     4             '  MSE_A7_COEF(',m,',',n,') = ',mse_a7_coef(m,n),
     4             '  MSE_ER(',m,',',n,')      = ',mse_er(m,n),
     4             '  MSE_EZ(',m,',',n,')      = ',mse_ez(m,n),
     7             '  MSE_VAC(',m,',',n,') = ',mse_vac(m,n),
     8             '  MSE_POL(',m,',',n,') = ',mse_pol(m,n),
     9             '  SIGMA_MSE_POL(',m,',',n,') = ',sigma_mse_pol(m,n)
               END DO
            END DO
         ELSE
            DO m=1,nmse_cams
               DO n=1,nmse_chords
               WRITE(iunit,'(9(a,i3,a,i2,a,1p,e13.6,3x))')
     1             '  MSE_R(',m,',',n,') = ',mse_r(m,n),
     2             '  MSE_PHI(',m,',',n,') = ',mse_phi(m,n),
     3             '  MSE_Z(',m,',',n,') = ',mse_z(m,n),
     4             '  MSE_ALPHA(',m,',',n,') = ',mse_alpha(m,n),
     5             '  MSE_BETA(',m,',',n,') = ',mse_beta(m,n),
     6             '  MSE_THETA(',m,',',n,') = ',mse_theta(m,n),
     7             '  MSE_VAC(',m,',',n,') = ',mse_vac(m,n),
     8             '  MSE_POL(',m,',',n,') = ',mse_pol(m,n),
     9             '  SIGMA_MSE_POL(',m,',',n,') = ',sigma_mse_pol(m,n)
               END DO
            END DO
         END IF
         IF (lmse_er) THEN
            WRITE (iunit, 1007) 'LMSE_ER = ', lmse_er
            i = minloc(er_aux_s(2:),DIM=1)
            IF (i > 4) THEN
               WRITE (iunit, 101) '  ER_AUX_S = ',(er_aux_s(n), n=1,i)
               WRITE (iunit, 101) '  ER_AUX_F = ',(er_aux_f(n), n=1,i)
            END IF
            i = minloc(ez_aux_s(2:),DIM=1)
            IF (i > 4) THEN
               WRITE (iunit, 101) '  EZ_AUX_S = ',(ez_aux_s(n), n=1,i)
               WRITE (iunit, 101) '  EZ_AUX_F = ',(ez_aux_f(n), n=1,i)
            END IF
         END IF
      END IF
c ---------------------------------------------------------------------------
c  MCZ & SAL 2011 -  Diagno simulations of magnetic diagnostics
c
      IF( ldiagno_opt) THEN
         WRITE (iunit, 100) TRIM(comments(13))
         WRITE (iunit,'(2x,3a)') "DIAGNO_CONTROL = '",
     1                            trim(diagno_control),"'"
         WRITE (iunit,'(2x,3a)') "DIAGNO_COIL = '",
     1                            trim(diagno_coil),"'"
         IF( ndiagno_bp > 0) THEN
            CALL write_array(iunit,'SIGMA_DIAGNO_BP', sigma_diagno_bp,
     1                       ndiagno_bp)
            CALL write_array(iunit,'TARGET_DIAGNO_BP',
     1                       target_diagno_bp, ndiagno_bp)
         END IF
         IF( ndiagno_flx > 0) THEN
            CALL write_array(iunit,'SIGMA_DIAGNO_FLX', sigma_diagno_flx,
     1                       ndiagno_flx)
            CALL write_array(iunit,'TARGET_DIAGNO_FLX',
     1                       target_diagno_flx, ndiagno_flx)
         END IF
         IF( ndiagno_seg > 0) THEN
            CALL write_array(iunit,'SIGMA_DIAGNO_SEG', sigma_diagno_seg,
     1                       ndiagno_seg)
            CALL write_array(iunit,'TARGET_DIAGNO_SEG',
     1                       target_diagno_seg, ndiagno_seg)
         END IF
      END IF

c ---------------------------------------------------------------------------
c  MCZ - Soft X-Ray Chords
c     write(iunit,981) 'N_EMIS = ', NP_EMIS
c     if( n_emis > 0 .and. any(sigma_emis > 0)) then
c        call write_array(iunit,'SIGMA_EMIS', sigma_emis, n_emis)
c        write (iunit, 990) 'SIGMA_EMIS_DAMP = ', sigma_emis_damp
c        call write_array(iunit,'AEMIS', aemis, 11)
c        write (iunit,'(2x,3a)') "EMIS_FILE = '",trim(emis_file),"'"
c        call write_array(iunit,'EMIS_CHORD', emis_chord, n_emis)
c     endif

c ---------------------------------------------------------------------------
!  DIAGNOSTICS (V3 FIT)
      WRITE (iunit, 100) TRIM(comments(10))
      WRITE (iunit,200) 'V3RFUN_DIR = ', TRIM(v3rfun_dir)
      WRITE (iunit,200) 'V3POST_IN  = ', TRIM(v3post_in)

      DO k = 1, ndiagno
         WRITE (iunit, 120) k,name_diagno(k),k,data_diagno(k),
     1                      k,sigma_diagno(k)
      END DO
 120  FORMAT ('NAME_DIAGNO(',i3.3,")='",a30,"', DATA_DIAGNO(",i3.3,")=",
     1         1pe12.5,", SIGMA_DIAGNO(",i3.3,")=",1pe12.4)
 200  FORMAT(2x,a,"'",a,"'")

c ---------------------------------------------------------------------------
c  VV Spec
      WRITE (iunit, 100) TRIM(comments(14))
      IF (sigma_vv < bigno .or. sigma_vv_rms < bigno) THEN
         WRITE (iunit,993) 'TARGET_VV = ',target_vv,
     1                     'SIGMA_VV = ',sigma_vv
         WRITE (iunit,993) 'TARGET_VV_RMS = ',  target_vv_rms,
     1                     'SIGMA_VV_RMS = ',sigma_vv_rms
         WRITE (iunit,'(a,i3)' ) '  MPOL_VV = ', mpol_vv
         WRITE (iunit,'(a,i3)' ) '  NTOR_VV = ', ntor_vv
         WRITE (iunit,'(a,i3)' ) '  NU_VV = ', nu_vv
         WRITE (iunit,'(a,i3)' ) '  NV_VV = ', nv_vv
         DO m = 0, mpol_vv
            DO n = -ntor_vv, ntor_vv
               IF (rbc_vv(n,m)==zero .and. zbs_vv(n,m)==zero) CYCLE
               IF (n < -9 .and. m >= 10) THEN
                  WRITE (iunit,'(2(a,i3,a,i2,a,1p,e13.6,3x))',err=2000)
     1                   '  RBC_VV(', n, ',', m, ') = ', rbc_vv(n,m),
     2                   '  ZBS_VV(', n, ',', m, ') = ', zbs_vv(n,m)
               ELSE IF (n < -9 .and. m < 10) THEN
                  WRITE (iunit,'(2(a,i3,a,i1,a,1p,e13.6,3x))',err=2000)
     1                   '  RBC_VV(', n, ',', m, ') = ', rbc_vv(n,m),
     2                   '  ZBS_VV(', n, ',', m, ') = ', zbs_vv(n,m)
               ELSE IF ((n < 0 .or. n > 9).and. m >= 10) THEN
                  WRITE (iunit,'(2(a,i2,a,i2,a,1p,e13.6,3x))',err=2000)
     1                   '  RBC_VV(', n, ',', m, ') = ', rbc_vv(n,m),
     2                   '  ZBS_VV(', n, ',', m, ') = ', zbs_vv(n,m)
               ELSE IF ((n < 0 .or. n > 9).and. m < 10) THEN
                  WRITE (iunit,'(2(a,i2,a,i1,a,1p,e13.6,3x))',err=2000)
     1                   '  RBC_VV(', n, ',', m, ') = ', rbc_vv(n,m),
     2                   '  ZBS_VV(', n, ',', m, ') = ', zbs_vv(n,m)
               ELSE IF ( m < 10 ) THEN
                  WRITE (iunit,'(2(a,i1,a,i1,a,1p,e13.6,3x))',err=2000)
     1                   '  RBC_VV(', n, ',', m, ') = ', rbc_vv(n,m),
     2                   '  ZBS_VV(', n, ',', m, ') = ', zbs_vv(n,m)
               ELSE
                  WRITE (iunit,'(2(a,i1,a,i2,a,1p,e13.6,3x))',err=2000)
     1                   '  RBC_VV(', n, ',', m, ') = ', rbc_vv(n,m),
     2                   '  ZBS_VV(', n, ',', m, ') = ', zbs_vv(n,m)

               ENDIF
            END DO
         END DO

         WRITE (iunit,'(a,L1)' ) '  SHAPEWEIGHT = ', shapeweight
         IF( shapeweight ) THEN
           WRITE (iunit, '(a,1p,3e13.6)') '  PLANES_BW = ',
     1                    (planes_bw(n),n=1,3)
           WRITE (iunit, '(a,1p,3e13.6)') '  AMPLW_BW = ',
     1                    (amplw_bw(n),n=1,3)
           WRITE (iunit, '(a,1p,3e13.6)') '  THETA0_BW = ',
     1                    (theta0_bw(n),n=1,3)
           WRITE (iunit, '(a,1p,e13.6)') '  PHI0_BW = ',phi0_bw
           WRITE (iunit, '(a,1p,e13.6)') '  WTHETA_BW = ', wtheta_bw
           WRITE (iunit, '(a,1p,e13.6)') '  WPHI_BW = ', wphi_bw
         ENDIF
         ! LIMITER MCZ - PPPL
         DO n = 1, size(phi_lim)
            IF( phi_lim(n) < -360 .or. phi_lim(n) > 360 .or.
     1          r_lim(1,n) == 0 .or. r_lim(2,n) == 0) EXIT

            m = 2
            DO i = 3, size(r_lim, 1)
                IF( r_lim(i,n) == 0) EXIT
                m = i
            END DO

            WRITE(temp_write,*) n
            temp_write = adjustl(temp_write)
            WRITE(iunit,'(a,a,a,1pe13.6)') ' PHI_LIM(',
     1               trim(temp_write),') = ',phi_lim(n)
            temp_write = '(1:,'//trim(temp_write)//')'
            CALL write_array(iunit,'R_LIM',r_lim(1:m,n),m,NDIM2=n)
            CALL write_array(iunit,'Z_LIM',z_lim(1:m,n),m,NDIM2=n)
!            CALL write_array(iunit, 'R_LIM'//trim(temp_write),
!     1                       r_lim(1:m,n), m)
!            CALL write_array(iunit, 'Z_LIM'//trim(temp_write),
!     1                       z_lim(1:m,n), m)
         END DO
      ENDIF
      ! Second Limiter Spec
      IF (sigma_bd < bigno .or. sigma_bd_rms < bigno .or.
     1    sigma_bd_max < bigno ) THEN
         WRITE (iunit, '(a,1pe13.6)') '  TARGET_BD = ',target_bd
         WRITE (iunit,'(a,1pe13.6)' ) '  TARGET_BD_RMS = ',target_bd_rms
         WRITE (iunit, '(a,1pe13.6)') '  SIGMA_BD = ',sigma_bd
         WRITE (iunit,'(a,1pe13.6)' ) '  SIGMA_BD_RMS = ',sigma_bd_rms
         WRITE (iunit,'(a,1pe13.6)' ) '  SIGMA_BD_MAX = ',sigma_bd_max

         WRITE (iunit,'(a,i3)' ) '  MPOL_BD = ', mpol_bd
         WRITE (iunit,'(a,i3)' ) '  NTOR_BD = ', ntor_bd
         IF( mpol_bd > 0 .or. ntor_bd > 0) THEN
            WRITE (iunit,'(a,i3)' ) '  NU_BD = ', nu_bd
            WRITE (iunit,'(a,i3)' ) '  NV_BD = ', nv_bd
            DO m = 0, mpol_bd
               DO n = -ntor_bd, ntor_bd
                 IF (rbc_bd(n,m)==zero .and. zbs_bd(n,m)==zero) CYCLE
                 IF (n < -9 .and. m >= 10) THEN
                   WRITE (iunit,'(2(A,i3,A,i2,A,1pe13.6,3x))',err=2000)
     1                   '  RBC_BD(', n, ',', m, ') = ', rbc_bd(n,m),
     2                   '  ZBS_BD(', n, ',', m, ') = ', zbs_bd(n,m)
                 ELSE IF (n < -9 .and. m < 10) THEN
                   WRITE (iunit,'(2(A,i3,A,i1,A,1pe13.6,3x))',err=2000)
     1                   '  RBC_BD(', n, ',', m, ') = ', rbc_bd(n,m),
     2                   '  ZBS_BD(', n, ',', m, ') = ', zbs_bd(n,m)
                 ELSE IF ((n < 0 .or. n > 9).and. m >= 10) THEN
                   WRITE (iunit,'(2(A,i2,A,i2,A,1pe13.6,3x))',err=2000)
     1                   '  RBC_BD(', n, ',', m, ') = ', rbc_bd(n,m),
     2                   '  ZBS_BD(', n, ',', m, ') = ', zbs_bd(n,m)
                 ELSE IF ((n < 0 .or. n > 9).and. m < 10) THEN
                   WRITE (iunit,'(2(A,i2,A,i1,A,1pe13.6,3x))',err=2000)
     1                   '  RBC_BD(', n, ',', m, ') = ', rbc_bd(n,m),
     2                   '  ZBS_BD(', n, ',', m, ') = ', zbs_bd(n,m)
                 ELSE IF ( m < 10 ) THEN
                   WRITE (iunit,'(2(A,i1,A,i1,A,1pe13.6,3x))',err=2000)
     1                   '  RBC_BD(', n, ',', m, ') = ', rbc_bd(n,m),
     2                   '  ZBS_BD(', n, ',', m, ') = ', zbs_bd(n,m)
                 ELSE
                   WRITE (iunit,'(2(A,i1,A,i2,A,1pe13.6,3x))',err=2000)
     1                   '  RBC_BD(', n, ',', m, ') = ', rbc_bd(n,m),
     2                   '  ZBS_BD(', n, ',', m, ') = ', zbs_bd(n,m)
                 ENDIF
               END DO
            END DO
         END IF
      END IF

      WRITE (iunit, '(a)') '/'

  100 FORMAT('!',40('-'),/,'!',5x,a,/,'!',40('-'))
  101 FORMAT(a,(1p,4e22.14))
 1010 FORMAT(8(1p,e16.6))

      WRITE(iunit, nml=bootin, err=2000)
      CALL write_gade_nml(iunit)

      IF (lcoil_geom .or. sigma_berr_avg < bigno .or.
     1    sigma_berr_max < bigno) call write_coilsin (iunit, istat)
      IF (istat .ne. 0) GOTO 2000

      IF (lvac_opt) CALL write_invac (iunit, istat)

      IF (istat .ne. 0) GOTO 2000

      CLOSE (iunit)

      RETURN

!
!     Handle errors here
!
 2000 istat = -5

      END SUBROUTINE write_indata
