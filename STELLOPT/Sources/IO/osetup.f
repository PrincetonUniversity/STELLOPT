      SUBROUTINE osetup(in_file)
      USE optim
      USE vmec_input, ONLY: ns_array, mpol, ntor, gamma, spres_ped,
     1    mgrid_file, lmac, lfreeb, rbc, extcur, raxis_cc, zaxis_cs,
     2    raxis_cs, zaxis_cc, bloat, precon_type, prec2d_threshold
      USE boozer_params, ONLY: xm_bdy, xn_bdy, rmnc_bdy, zmns_bdy,
     1    rmnc_opt, zmns_opt, lmns_opt, rmns_opt, zmnc_opt, lmnc_opt
      USE vparams, ONLY: one, zero, ntord, mpol1d
      USE system_mod
      USE bootsj_input, ONLY: damp_bs
      USE mpi_params                                                     !MPI
      USE safe_open_mod
      USE vacfield_mod, ONLY: coils_file_extension
      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      include 'mpif.h'                                                   !MPI
!DEC$ ENDIF
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      CHARACTER(LEN=*) :: in_file
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
!DEC$ IF DEFINED(MCURIE)
      INTEGER, PARAMETER :: max_process=32
!DEC$ ELSEIF DEFINED(CRAY)
      INTEGER, PARAMETER :: max_process=8
!DEC$ ELSEIF DEFINED(OSF1)
      INTEGER, PARAMETER :: max_process=1
!DEC$ ELSE
      INTEGER, PARAMETER :: max_process=1
!DEC$ ENDIF
      REAL(rprec), PARAMETER :: p360 = 360
      CHARACTER(LEN=*), PARAMETER :: input_ext = 'input.'
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: js, k, k1, k2, ierr, istat, ioutp
      LOGICAL :: lexist, ljacobian, lcoilsdot
      INTEGER, DIMENSION(ini_max):: zetain, thetain                      !COBRA
      CHARACTER(LEN=256) :: stripped_input_file, temp, temp2, temp3
      CHARACTER(LEN=20) :: ch_myid
      CHARACTER(LEN=256) ::v3post_save
C-----------------------------------------------
      CALL getenv('HOME', home_dir)
!DEC$ IF DEFINED (FALCON)
      home_dir = TRIM(home_dir) // '/falcon/'
!DEC$ ELSEIF DEFINED (WIN32)
      home_dir = "D:\bin\"
!DEC$ ELSE
      home_dir = TRIM(home_dir) // '/bin_847/'
!DEC$ ENDIF

!
!     DELETE OUTPUT FILE IF IT ALREADY EXISTS
!

      output_file = 'output.' // TRIM(seq_ext)
      INQUIRE (file=TRIM(output_file), exist=lexist, iostat=istat)
      IF (lexist) THEN
         OPEN(unit=iout, file=output_file)
         CLOSE (iout, status='delete', iostat=istat)
      END IF

      num_levmar_params = MAX(3, max_process/4)           !!Number of lm PARAMETER estimates/jacobian evaluation
      num_processors = max_process                        !!Number of processors to distribute jacobian evaluations over

      CALL initialize_coilsin

!---------------------------------------------------------------------------------
!     Default values for BOOTIN parameters
!
      damp_bs = -1
      rmax_opt = 0; rmin_opt = 0; zmax_opt = 0
      rgrid_max = 0; rgrid_min = 0; zgrid_max = 0; zgrid_min = 0
      nextcur_opt = 0
      nextcur_vmec = 0

!------------------------------------------------------------------------------

      CALL GA_preset           !  initialize the GA parameters
      CALL DE_preset           !  initialize the DE parameters

!------------------------------------------------------------------------------
!
!     Read input (optimum) namelist
!     first, strip any comments from input file (F95-not necessary)
!

!     Produce clean file input_file//'.stripped'
      WRITE (ch_myid, *) myid
      stripped_input_file = TRIM(in_file) // '.stripped' 
     1                      // ADJUSTL(ch_myid)
      CALL strip_comments(in_file,stripped_input_file)                 !MPI
!DEC$ IF DEFINED (MPI_OPT)
!      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr_mpi)                         !MPI
!      IF (ierr_mpi .ne. 0) STOP 'MPI_BARRIER error in STELLOPT Osetup'
!DEC$ ENDIF
!      stripped_input_file = TRIM(in_file) // '.stripped'
!      WRITE (ch_myid, *) myid
!      temp = copy // TRIM(stripped_input_file) // " " //
!     1               TRIM(stripped_input_file) // ADJUSTL(ch_myid)
!      CALL system(TRIM(temp))
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr_mpi)                         !MPI
      IF (ierr_mpi .ne. 0) STOP 'MPI_BARRIER error in STELLOPT Osetup'
!DEC$ ENDIF

!      temp = stripped_input_file
!      stripped_input_file = TRIM(stripped_input_file)//ADJUSTL(ch_myid)
      CALL read_input_opt (stripped_input_file, istat)                   !read in_file and delete it
      IF (istat .ne. 0) THEN
         PRINT *,'Error reading input namelist in STELLOPT',
     1           ' osetup routine: istat = ', istat
         STOP
      END IF

      lcoilp_sep = lfreeb .and. .not.lcoil_geom
      lcoilsdot = .false.

!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr_mpi)
      IF (ierr_mpi .ne. 0) STOP 'MPI_BARRIER error in STELLOPT Osetup'
!DEC$ ENDIF
      temp = remove // TRIM(stripped_input_file)
!      temp = remove // TRIM(temp)
      CALL system (temp)
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr_mpi)
      IF (ierr_mpi .ne. 0) STOP 'MPI_BARRIER error in STELLOPT Osetup'
!DEC$ ENDIF
      
      IF (myid .ne. master) GOTO 100

!
!        CHECK IF PLASMA-COIL SEPARATION MUST BE CALCULATED FOR FREE-BDY, FIXED COIL GEOMETRY
!        (lcoilp_sep = true); IMPLEMENTATION NOT COMPLETE (09/29/03)
!
      IF (lvac_opt) THEN
!
!         Check if coils. file exists.
!
         coils_dot_file = "coils." // TRIM(coils_file_extension)
         temp = ".." // sep // coils_dot_file
         INQUIRE (FILE=temp, EXIST=lvac_opt, iostat=istat)
         IF (istat.ne.0 .or. .not.lvac_opt) THEN
            PRINT *,' No coils-dot file found. Unable to compute' //
     1             ' vacopt chi-sq metrics.'
            STOP
         END IF
         temp = link // ".." // sep // TRIM(coils_dot_file)
     1               // " ." // sep // TRIM(coils_dot_file)
         CALL system(temp)
         lcoilsdot = .true.
         lcoilp_sep = .false.
      END IF

      IF (lcoilp_sep) THEN
         k = INDEX(mgrid_file,sep)          !!If NOT full path, must link to scratch directory
         IF (k .eq. 0) THEN
            temp = link // ".." // sep // TRIM(mgrid_file) 
     1                  // " ." // sep // TRIM(mgrid_file)
            CALL system(temp)
         END IF
         INQUIRE (FILE=mgrid_file, EXIST=lcoilp_sep, iostat=istat)
         IF (istat.ne.0 .or. .not.lcoilp_sep) THEN
            lfreeb = .false.
            PRINT *,' No mgrid file found: running in fixed bdy mode'

         ELSE
!
!           Check if coils. file exists AND has same extension as mgrid file
!
            k = INDEX(mgrid_file,'.',BACK=.TRUE.)
            js = INDEX(mgrid_file,'_')
            IF (js.ne.0 .and. js.lt.k) THEN
               temp = mgrid_file(js+1:k-1)
            ELSE
               temp = mgrid_file(k+1:)
            END IF
            coils_dot_file = 'coils.' // TRIM(temp)
            temp = '..' // sep // coils_dot_file
            INQUIRE (FILE=temp, EXIST=lcoilp_sep, iostat=istat)
            IF (istat.ne.0 .or. .not.lcoilp_sep) THEN
               PRINT *,' No coils-dot file found. Unable to compute',
     1                 ' coils-plasma separation.'
            ELSE
               temp = link // '..' // sep // TRIM(coils_dot_file)
     1                      // ' .' // sep // TRIM(coils_dot_file)
                  CALL system(temp)
            END IF
         END IF
      END IF

!
!        CHECK IF DIAGNOSTICS MATCHING SHOULD BE PERFORMED
!
      ndiagno = 0
      IF (.not.lv3post) GOTO 100

      v3post_save = v3post_in 
      temp = v3post_in
      INQUIRE (FILE=temp, EXIST=lv3post, iostat=istat)
      IF (istat.ne.0 .or. .not.lv3post) THEN
         temp = ".." // sep // v3post_in
         INQUIRE (FILE=temp, EXIST=lv3post, iostat=istat)
         IF (.not.lv3post) THEN
            PRINT *,'Diagnostics input file v3post_in=',
     1            TRIM(v3post_in),' not found.'
            STOP 'Terminating run.'
         END IF
      END IF

!
!     Copy v3post_in to this working directory ('.')
!
      k = INDEX(v3post_in, sep, BACK=.TRUE.)
      v3post_in = v3post_in(k+1:)
      temp = copy // '"' // TRIM(temp) // '"' // " ." 
     1             // sep // TRIM(v3post_in)
!DBG ONLY         PRINT *,'temp=',TRIM(temp)
!DBG ONLY         PRINT *,'v3post_in=',TRIM(v3post_in)
      CALL system(temp)
!
!      Rewrite Diagnostics List with FULL path name
!      (link will not work on MS Windows system, and response
!      function files are too big to copy)
!
      js = 1002
      CALL safe_open (js, istat, v3post_in, 'old', 
     1                'formatted')
      IF (istat .ne. 0) THEN
         lv3post = .false.
         STOP ' Unable to open v3post input file.'
      END IF

      DO WHILE (istat .eq. 0)
         READ(js, '(a)', iostat=istat) temp
         k1 = INDEX(temp,"listin")
         IF (k1 .gt. 0) EXIT
      END DO

      lv3post = (k1 .gt. 0)
      IF (.not.lv3post) THEN
         PRINT *, ' V3POST Diagnostics list file name',
     1            ' not found in v3post input file.'
         STOP
      END IF

! EAL: NEVER use a .LIST file from directory ../ (one above this working dir)
! if it is there, it has virtually ALWAYS been modified to fully
! qualified response function file names and would end up
! path1/path2/file after editing below !!!!!

!     PARSE NAME OF .LIST FILE CONTAINING NAMES OF DIAGNOSTIC RESPONSE FILES
!     LOOK IN v3rfun_dir

      temp2 = temp(k1+1:)
      k = INDEX(temp2,"'")
      temp = temp2(k+1:)
      k = INDEX(temp,"'")
      IF (k .gt. 1) temp = temp(1:k-1)
      temp2 = TRIM(v3rfun_dir) // sep // ADJUSTL(temp)  
      INQUIRE (FILE=temp2,EXIST=lv3post,iostat=istat)
      IF (istat.ne.0 .or. .not.lv3post) 
     1   STOP 'V3POST Diagnostics list file not found.'

!     COPY CONTENTS OF .LIST FILE TO FILE CALLED 'TEMP' IN THIS DIRECTORY
!     STRIPPING ANY PATHS AND PUTTING IN THE V3RFUN_DIR PATH 
!
!     IMPORTANT: THIS ASSUMES ALL THE DIAGNOSTICS ARE IN THE V3RFUN_DIR DIRECTORY
!
      ioutp = js+1
      CALL safe_open (ioutp, k, temp, 'replace', 'formatted')

      ierr = ioutp+1
      CALL safe_open (ierr, k, temp2, 'old', 'formatted')
!     READ THE LIST AND MAKE SURE DIAGNOSTIC FILES EXIST

      DO WHILE (k .eq. 0)
         READ(ierr, '(a)', iostat=k) temp
         IF (k .ne. 0) EXIT
         IF (LEN_TRIM(temp) .le. 0) CYCLE     !Allow blank lines in file
         ndiagno = ndiagno + 1
         k2 = INDEX(temp,' ')
         temp3 = temp(k2+1:)
!     STRIP ANY PATHS FROM EXISTING DIAGNOSTIC NAME
         k1 = SCAN(temp3,sep,BACK=.TRUE.)
         temp2 = temp(1:k2) // TRIM(v3rfun_dir) 
     1                      // sep // temp3(k1+1:)
         WRITE (ioutp, '(a)') TRIM(temp2)
      END DO

      CLOSE (ierr)
      CLOSE (ioutp)
                  
!
!     Copy coil response file (small) to this working directory "."
!
      temp = copy // '"' // TRIM(v3rfun_dir) // '"' // 
     1       sep // 'crfun_*' // ' .' 

      CALL system (temp)

      CLOSE (js)

      v3post_in =  v3post_save  ! this is the one indata will write
 

 100  CONTINUE
!DEC$ IF DEFINED (MPI_OPT)
      IF (lvac_opt .and. .not.lcoilsdot) STOP
      CALL MPI_BCAST(lfreeb, 1, MPI_LOGICAL, master, MPI_COMM_WORLD,
     1                  ierr_mpi)
      CALL MPI_BCAST(lv3post, 1, MPI_LOGICAL, master, MPI_COMM_WORLD,
     1                  ierr_mpi)
      IF (lv3post) THEN
         CALL MPI_BCAST(ndiagno, 1, MPI_INTEGER, master, MPI_COMM_WORLD,
     1                  ierr_mpi)
      END IF
      IF (ierr_mpi .ne. 0) STOP 'MPI_BCAST error #2 in STELLOPT Osetup'
!DEC$ ENDIF
      min_ext = TRIM(seq_ext)
      lniter1 = ((niter_opt.eq.1) .or. lone_step .or. lscale_only)
      IF (.not.lniter1) min_ext = TRIM(min_ext) // '.min'
      min_input_file = input_ext // TRIM(min_ext)
!DEC$ IF DEFINED (NETCDF)
!WILL ADD .nc EXTENSION IN CLEAN-UP ROUTINE      
      min_wout_file = 'wout_' // TRIM(min_ext)
!DEC$ ELSE
      min_wout_file = 'wout.' // TRIM(min_ext)
!DEC$ ENDIF
      min_diagno_file = 'chisq_diagno.' // TRIM(min_ext)

      IF (r00_opt.gt.0 .and. rbc(0,0).ne.zero)
     1    r00_scale = r00_opt/rbc(0,0)
      IF (ABS(r00_scale*b00_scale) .ne. one)
     1   CALL init_scale(r00_scale, b00_scale)

      raxis_old(:,1) = raxis_cc(:);  raxis_old(:,2) = raxis_cs(:)       !Save for writing into restart (.min) file
      zaxis_old(:,1) = zaxis_cs(:);  zaxis_old(:,2) = zaxis_cc(:)
      jcurv_opt = 0;   iota_opt = 0
!
!     COMPATABILITY WITH LPK ADDITIONS
!
      lpres_prof_opt = lpres_prof_opt .or. lpress_opt
      lnescoil_opt = lnescoil_opt .or. lcoil_complex .or. lcoil_opt
      lbootsj_opt = lbootsj_opt .or. lbootstrap
      lballoon_opt = lballoon_opt .or. lbal_opt
      IF (nproc .ge. num_processors) num_processors = nproc
      target_coil_complex = MAX(one, target_coil_complex)
      IF (damp_bs .eq. zero) damp_bs = -1                                !Old-style file: did not mean 0...
      IF (mpol_vv.gt.mpol1d .or. ntor_vv.gt.ntord) STOP
     1    'mpol_vv > mpol1d or ntor_vv > ntord exceed bounds'

      NumJstar = MIN(NumJstard, NumJstar)
      NumJinvariant = MIN(NumJstard, NumJinvariant)
      mpol1_opt = mpol - 1                                               !need in initialize_opt
      mrho1_opt = mpol1_opt - 1                                          
      IF (nopt_boundary .eq. 1) mrho1_opt = mpol1_opt + 1                !Advanced HB representation
      ntor_opt  = ntor
      mnmax_opt = ntor + 1 + (mpol1_opt)*(1 + 2*ntor)

      IF (lvac_opt) lcoil_geom = .false.
      IF (lcoil_geom) lfreeb = .true.

      IF (lfreeb) THEN
         INQUIRE(file=TRIM(mgrid_file), exist=lexist, iostat=istat)
         IF ((istat.ne.0 .or. .not.lexist) .and. .not.lcoil_geom) THEN
            IF (myid.eq.master) PRINT *,' Could not locate mgrid file!'
            lfreeb = .false.
         END IF
         nextcur_opt = count(lextcur)
         IF (lfreeb .and. lnescoil_opt .and. myid.eq.master) PRINT *,
     1   ' NESCOIL Constraints are being used in free-boundary mode!'
      END IF

      IF (lvac_opt) lfreeb = .true.
!
!     Backward compatability
!
      SELECT CASE(TRIM(sym_type))
        CASE ('QA')
           helicity = CMPLX(one,zero)
        CASE ('QH')
           helicity = CMPLX(one,-one)
        CASE ('QP')
           helicity = CMPLX(zero,one)
      END SELECT

!
!     Check consistency of LOGICAL variables with sigmas
!
      IF (ALL(sigma_bmn .ge. bigno) .and. lbmn) THEN
         IF (myid .eq. master)
     1   PRINT *,' LBMN is being set to FALSE; specify SIGMA_BMN array'
         lbmn = .false.
      ENDIF
      IF (ALL(sigma_bootsj .ge. bigno) .and. lbootsj_opt) THEN
         IF (myid .eq. master)
     1   PRINT *,' LBOOTSJ_OPT is being set to FALSE; ',
     2           'specify SIGMA_BOOTSJ array'
         lbootsj_opt = .false.
      ENDIF
      IF (ALL(sigma_balloon .ge. bigno) .and. lballoon_opt) THEN
         IF (myid .eq. master)
     1   PRINT *,' LBALLOON_OPT is being set to FALSE; ',
     2           'specify SIGMA_BALLOON array'
         lballoon_opt = .false.
      ENDIF
      IF (ALL(sigma_dkes .ge. bigno) .and. ldkes_opt) THEN
         IF (myid .eq. master)
     1   PRINT *,' LDKES_OPT is being set to FALSE; ',
     2           'specify SIGMA_DKES array'
         ldkes_opt = .false.
      ENDIF


      IF (REAL(helicity) .lt. zero) helicity = -helicity
      IF (lbmn .and. helicity.eq.CMPLX(zero,zero)) THEN
      IF (myid .eq. master) THEN                                         !MPI
        PRINT *,' Must specify helicity (symmetry TYPE) for lbmn = T'
        PRINT *,' Allowable values for helicity are:'
        PRINT *,'    (1,0)   =>   quasi-axisymmetry'
        PRINT *,'    (0,1)   =>   quasi-poloidal symmetry'
        PRINT *,'    USE with quasi-omnigeneity to reduce 1/R drift'
        PRINT *,'    (l,k)   =>   quasi-helical symmetry'
        PRINT *,'    |B| = F(lu + kv)'
      ENDIF                                                              !MPI
        STOP 'helicity MUST be supplied in OPTIMUM NAMELIST'
      END IF

!      ALLOCATE (rmnc_bdy(mnmax_opt), zmns_bdy(mnmax_opt),
!     1          xm_bdy(mnmax_opt), xn_bdy(mnmax_opt), stat=istat)

      nrad = MAXVAL(ns_array)
!      ALLOCATE (ns_surf(nrad), ns_booz(nrad), ns_ball(nrad),
!     1          ns_neo(nrad), lneed_modB(nrad), stat = js)
      allocate (rmnc_bdy(mnmax_opt), zmns_bdy(mnmax_opt),
     1          xm_bdy(mnmax_opt), xn_bdy(mnmax_opt),
     2          rmnc_opt(mnmax_opt,nrad), zmns_opt(mnmax_opt,nrad),
     3          lmns_opt(mnmax_opt,nrad),
     2          rmns_opt(mnmax_opt,nrad), zmnc_opt(mnmax_opt,nrad),
     3          lmnc_opt(mnmax_opt,nrad),
     4          ns_surf(nrad), ns_booz(nrad), ns_ball(nrad),
     5          ns_neo(nrad), lneed_modB(nrad), stat = istat)
      if (istat.ne.0 ) stop 'Allocation error in OSETUP'
!      IF (istat.ne.0 .or. js.ne.0) STOP 'Allocation error in OSETUP'
!
!     Make old-style, new-style masks conform
!
      CALL update_style (lsurf_mask, nsurf_mask, nrad)
      CALL update_style (ldkes_mask, ndkes_mask, nrad)
      CALL update_style (lballoon_mask, nballoon_mask, nrad)             !VMEC COBRA (RS)
      CALL update_style (lneo_mask, nneo_mask, nrad)                     !NEO

!     IF tailoring the pressure profile: ALL surfaces must be included to compute gradient
      IF (lpres_prof_opt) lballoon_mask = .true.

      lneed_modB = .false.
      lballoon_mask(1) = .false.                !!MUST ignore axis point in COBRA   (RS)
      lballoon_mask(nrad) = .false.             !!MUST ignore last point in COBRA   (RS)

      ns_surf_max = 0                           !!Boozer surfaces based on lsurf_mask
      ns_booz_max = 0                           !!TOTAL NUMBER BOOZER SURFACES
      ns_ball_max = 0                           !!COBRAVMEC   (RS)
      ns_neo_max  = 0
      ljacobian = ANY(ABS(sigma_jac(:)) < bigno)

      DO js = 2, nrad
         IF (lsurf_mask(js)) THEN
           ns_surf_max = ns_surf_max + 1
           ns_surf(ns_surf_max) = js
         END IF
!!    NEO
         IF (lneo_mask(js)) THEN
           ns_neo_max = ns_neo_max + 1
           ns_neo(ns_neo_max) = js
         ENDIF

         IF (lsurf_mask(js) .or. (ldkes_mask(js) .and. ldkes_opt) .or.
     1       ljacobian .or.      (lneo_mask(js)  .and. lneo_opt)) THEN
           ns_booz_max = ns_booz_max + 1
           ns_booz(ns_booz_max) = js
         END IF
         IF(lballoon_mask(js)) THEN                                      !COBRAVMEC   (RS)
           ns_ball_max = ns_ball_max + 1                                 !COBRAVMEC   (RS)
           ns_ball(ns_ball_max) = js                                     !COBRAVMEC   (RS)
         ENDIF                                                           !COBRAVMEC   (RS)
      END DO

!
!     IF USER FORGOT TO SPECIFY LBALLOON_MASK ARRAY, DEFAULT TO NS_SURF...
!
      IF (lballoon_opt .and. ns_ball_max.eq.0) THEN
         ns_ball_max = ns_surf_max
         ns_ball(1:ns_ball_max) = ns_surf(1:ns_surf_max)
         IF (myid .eq. master) THEN                                      !MPI
            PRINT *,' LBALLOON_OPT = TRUE but no nonzero values for',
     1              ' NBALLOON_MASK were found'
            PRINT *,' Using NSURF_MASK for ballooning surfaces'
            PRINT *,' If this is wrong, set LBALLOON_OPT = FALSE,'
            PRINT *,' or set the NBALLOON_MASK array values in',
     1              ' the OPTIMUM NAMELIST'
         ENDIF                                                           !MPI
      END IF
      
!
!   if feeding back on pressure profile, must vary some pressure coefs
!
      if( lpres_prof_opt) then
         if( pres_opt_nmax < 0) then
            print *,'**********************************************'
            print *,'WARNING:'
            print *,' LPRES_PROF_OPT = TRUE but pres_opt_nmax < 0 '
            print *,' Incompatible!!'
            print *,' PRES_OPT_NMAX set to 1'
            print *,'**********************************************'

            pres_opt_nmax = 1
         else if( pres_opt_nmax > 10) then
            pres_opt_nmax = 10

         endif
      endif

!
!   set np_prof on basis of sigma_p_prof
!

      if( np_prof == 0 .and. any(sigma_p_prof<bigno)) then
         np_prof = 2

         do while(any(sigma_p_prof(np_prof:)<bigno) .and.
     1            np_prof < size(sigma_p_prof))
            np_prof = np_prof + 1
         enddo
         np_prof = np_prof - 1
      endif
     
!
!  set ANIMEC default behavior
!
      IF (      lanimec
     1    .and. (.not. lani_bcrit)
     2    .and. (.not. lani_phot)
     3    .and. (.not. lani_tperp)) THEN
         lani_bcrit = .true.
         lani_phot  = .true.
         lani_tperp = .true.
      END IF
      
      IF (       lanimec 
     1    .and.  lani_phot
     2    .and. (.not. ANY(ah_mask .ne. 0))) THEN
         ah_mask = 1
      END IF
      
      IF (       lanimec 
     1    .and.  lani_tperp
     2    .and. (.not. ANY(at_mask .ne. 0))) THEN
         at_mask = 1
      END IF
           
!
!  set AC_MASK default behavior
!
      IF ( lcur_prof_opt .and. (.not. ANY(ac_mask .ne. 0))) THEN
         ac_mask = 1
      END IF
!
!  set ndiagno_seg and ndiagno_flx from respective sigmas
!

      if( ldiagno_opt .and. any(sigma_diagno_seg < bigno)) then
         ndiagno_seg = count(sigma_diagno_seg < bigno)

         do while( any(sigma_diagno_seg(ndiagno_seg:)<bigno) .and.
     1             ndiagno_seg < size(sigma_diagno_seg))
            ndiagno_seg = ndiagno_seg + 1
         enddo
         ndiagno_seg = ndiagno_seg - 1
      else
         ndiagno_seg = 0
      endif

      if( ldiagno_opt .and. any(sigma_diagno_flx < bigno)) then
         ndiagno_flx = count(sigma_diagno_flx < bigno)

         do while( any(sigma_diagno_flx(ndiagno_flx:)<bigno) .and.
     1             ndiagno_flx < size(sigma_diagno_flx))
            ndiagno_flx = ndiagno_flx + 1
         enddo
         ndiagno_flx = ndiagno_flx - 1

         if((.not. lpres_opt_edge0) .and. myid==master) then
            print *,'************************************************'
            print *,'WARNING:'
            print *,'Matching flux loops, but LPRES_OPT_EDGE0 = FALSE'
            print *,'=> will not see all of the plasma pressure'
            print *,'************************************************'
         endif

      else
         ndiagno_flx = 0
      endif

      if( ldiagno_opt .and. any(sigma_diagno_bp < bigno)) then
         ndiagno_bp = count(sigma_diagno_bp < bigno)

         do while( any(sigma_diagno_bp(ndiagno_bp:)<bigno) .and.
     1             ndiagno_bp < size(sigma_diagno_bp))
            ndiagno_bp = ndiagno_bp + 1
         enddo
         ndiagno_bp = ndiagno_bp - 1

      else
         ndiagno_bp = 0
      endif

!     Not necessary (DIAGNO no longer require VMEC DIAGNO file)
!      if( ldiagno_opt .and. ndiagno_seg + ndiagno_flx > 0) then
!         ldiagno = .true.
!      endif

!---------------------------------------------------------------------------------
!   CODE added by R.SANCHEZ (02/01/99):  count number of initial theta values in degrees
!              (NINI_THETA), initial zeta values in degrees (NINI_ZETA) and total number
!              of initial positions (NINI_TOT) in NAMELIST OPTIMUM WHERE
!              ballooning growth rates are to be evaluated.

      IF (lballoon_opt .or. lpres_prof_opt) THEN
         IF(MINVAL(bal_theta0).lt.zero  .or. MAXVAL(bal_theta0).ge.p360
     1   .or. MINVAL(bal_zeta0).lt.zero .or. MAXVAL(bal_zeta0) .ge.p360)
     2   THEN
      IF (myid .eq. master) THEN                                            ! MPI
           PRINT *, 'ALL initial angles (in degrees) MUST be',
     1            ' 0 <= angle < 360 in STELLOPT input file!'
      ENDIF                                                                 ! MPI
           STOP
         END IF
         nini_theta=1
         nini_zeta=1
         thetain(1)=bal_theta0(1)
         zetain(1)=bal_zeta0(1)
         DO k=2, ini_max
            IF (bal_theta0(k) .ge. zero) THEN
               DO js=1, k-1
                  IF(bal_theta0(k).eq.bal_theta0(js)) EXIT
                  IF(js .eq. k-1) THEN
                     nini_theta=nini_theta+1
                     thetain(nini_theta)=bal_theta0(k)
                  ENDIF
               ENDDO
            ENDIF
            IF (bal_zeta0(k) .ge. zero) THEN
               DO js=1, k-1
                  IF (bal_zeta0(k).eq.bal_zeta0(js)) EXIT
                  IF (js .eq. k-1) THEN
                     nini_zeta=nini_zeta+1
                     zetain(nini_zeta)=bal_zeta0(k)
                  ENDIF
              ENDDO
            ENDIF
         ENDDO
         bal_theta0(1:nini_theta)=thetain(1:nini_theta)
         bal_zeta0(1:nini_zeta)=zetain(1:nini_zeta)
         nini_tot=nini_zeta*nini_theta
      ENDIF
!---------------------------------------------------------------------------------
!
!     Maintain Backwards compatability
!
      IF (sigma_rmin .ge. bigno) sigma_rmin = sigma_centering
      IF (sigma_rmax .ge. bigno) sigma_rmax = sigma_centering

      IF (epsfcn .le. zero) THEN
         IF (lreset_opt) THEN
            epsfcn = 1.e-4_dp
         ELSE
            epsfcn = 1.e-3_dp
         ENDIF
      ENDIF

      END SUBROUTINE osetup
