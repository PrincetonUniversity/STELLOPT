!-----------------------------------------------------------------------
!     Subroutine:    stellopt_init
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/26/2012
!     Description:   This subroutine initializes the vector of variables
!                    utilized to optimize the VMEC equilibria.  Note
!                    that for profile coefficients we normalize to the
!                    area under the curve.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_init
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_input_mod
      USE stellopt_targets, ONLY: write_targets
      USE stellopt_vars
      USE equil_utils, ONLY: count_vars
      USE vmec_input
      USE safe_open_mod, ONLY: safe_open
      USE equil_utils, ONLY: profile_norm
      USE vmec_params, ONLY: norm_term_flag, bad_jacobian_flag,&
                             more_iter_flag, jac75_flag, input_error_flag,&
                             phiedge_error_flag, ns_error_flag, &
                             misc_error_flag, successful_term_flag, &
                             restart_flag, readin_flag, timestep_flag, &
                             output_flag, cleanup_flag, reset_jacdt_flag
!                             animec_flag, flow_flag
      USE parallel_vmec_module, ONLY: PARVMEC, gnranks
      USE mpi_params
      USE mpi_inc
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        iunit       File unit number
!----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER ::  i,n,m,ier, iunit,nvar_in,ctrl_dofs,nknots
      INTEGER ::  ictrl(5)
      REAL(rprec) :: norm, delta
      REAL(rprec) :: fvec_temp(1)
      REAL(rprec), DIMENSION(-ntord:ntord,0:mpol1d) :: rbc_temp,zbs_temp
      REAL(rprec), PARAMETER :: norm_fac = 0.5_rprec   ! Used to set bounds for nomalization
!      REAL(rprec), PARAMETER :: pct_domain = 0.05 ! USed to auto-determine domain

      LOGICAL                                      :: ltst
      integer                                      :: color,key, nprocs_total
      CHARACTER(256)                               :: tstr1,tstr2
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      chisq_min = bigno
      ier = 0

      ! Read the OPTIMUM Namelist
      CALL read_stellopt_input(TRIM(id_string),ier,myid)
      !CALL bcast_vars(master,MPI_COMM_STEL,ierr_mpi)
      !IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'stellot_init:bcast_vars',ierr_mpi)

      ! Handle coil geometry
      IF (lcoil_geom) CALL namelist_input_makegrid(id_string)

      ! Handle MPI and shared memory
      CALL stellopt_init_mpi

      ! Send workers to the worker pool subroutine
      IF (myworkid .ne. master) THEN
         ltst  = .false.
         tstr1 = ''
         tstr2 = ''
         ier_paraexe = 0
         CALL stellopt_paraexe(tstr1,tstr2,ltst)
         RETURN
      END IF

      ! Read the Equilibrium input
      CALL tolower(equil_type)
      SELECT CASE (TRIM(equil_type))
         CASE('vmec2000','animec','flow','satire','parvmec','paravmec','vmec2000_oneeq')
              ltst = .false.
              tstr1 = 'parvmec_init'
              id_string = id_string(7:LEN(id_string))
              tstr2 = id_string
              CALL stellopt_paraexe(tstr1,tstr2,ltst)
         CASE('vboot')
              ltst = .false.
              tstr1 = 'parvmec_init'
              id_string = id_string(7:LEN(id_string))
              tstr2 = id_string
              ! Read BOOTSJ NAMELIST
              CALL safe_open(iunit,ier,'input.'//TRIM(id_string),'old','formatted')
              CALL read_namelist (iunit, ier, 'bootin')
              IF (ier < 0 .and. myid == master) THEN
                 WRITE(6,*) '!!!!!!!!!!!!ERRROR!!!!!!!!!!!!!!'
                 WRITE(6,*) '  BOOTIN Namelist not found     '
                 WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                 STOP
              END IF
              CLOSE(iunit)
              CALL stellopt_paraexe(tstr1,tstr2,ltst)
         CASE('test')
              id_string = id_string(7:LEN(id_string))
              ! We don't read indata but use it later so need to do some defaulting
              pcurr_type = 'power_series'
              piota_type = 'power_series'
              pmass_type = 'power_series'
              ph_type = 'power_series'
              pt_type = 'power_series'
      END SELECT
      
      ! Count the number of variables and allocate the array
      iter  = 0
      nvars = 0
      SELECT CASE (TRIM(equil_type))
         CASE('vmec2000','flow','animec','satire','paravmec','parvmec','vboot','vmec2000_oneeq')
              ! Single Values
              IF (lphiedge_opt) nvars = nvars + 1
              IF (lcurtor_opt)  nvars = nvars + 1
              IF (lpscale_opt)  nvars = nvars + 1
              IF (lbcrit_opt)   nvars = nvars + 1
              IF (lmix_ece_opt)   nvars = nvars + 1
              IF (lxics_v0_opt)   nvars = nvars + 1
              ! Arrays
              IF (ANY(lextcur_opt)) nvars = nvars + COUNT(lextcur_opt)
              IF (ANY(laphi_opt)) nvars = nvars + COUNT(laphi_opt)
              ! Profiles
              IF (ANY(lam_opt)) nvars = nvars + count_vars(lam_opt,am,pmass_type)
              IF (ANY(lac_opt)) nvars = nvars + count_vars(lac_opt,ac,pcurr_type)
              IF (ANY(lai_opt)) nvars = nvars + count_vars(lai_opt,ai,piota_type)
              IF (ANY(lah_opt)) nvars = nvars + count_vars(lah_opt,ah,ph_type)
              IF (ANY(lat_opt)) nvars = nvars + count_vars(lat_opt,at,pt_type)
              IF (ANY(lne_opt)) nvars = nvars + count_vars(lne_opt,ne_opt,ne_type)
              IF (ANY(lte_opt)) nvars = nvars + count_vars(lte_opt,te_opt,te_type)
              IF (ANY(lti_opt)) nvars = nvars + count_vars(lti_opt,ti_opt,ti_type)
              IF (ANY(lth_opt)) nvars = nvars + count_vars(lth_opt,th_opt,th_type)
              IF (ANY(lzeff_opt)) nvars = nvars + count_vars(lzeff_opt,zeff_opt,zeff_type)
              IF (ANY(lam_s_opt)) nvars = nvars + COUNT(lam_s_opt)
              IF (ANY(lac_s_opt)) nvars = nvars + COUNT(lac_s_opt)
              IF (ANY(lai_s_opt)) nvars = nvars + COUNT(lai_s_opt)
              IF (ANY(lphi_s_opt)) nvars = nvars + COUNT(lphi_s_opt)
              IF (ANY(lam_f_opt)) nvars = nvars + count_vars(lam_f_opt,am_aux_f,pmass_type)
              IF (ANY(lac_f_opt)) nvars = nvars + count_vars(lac_f_opt,ac_aux_f,pcurr_type)
              IF (ANY(lai_f_opt)) nvars = nvars + count_vars(lai_f_opt,ai_aux_f,piota_type)
              IF (ANY(lah_f_opt)) nvars = nvars + count_vars(lah_f_opt,ah_aux_f,ph_type)
              IF (ANY(lat_f_opt)) nvars = nvars + count_vars(lat_f_opt,at_aux_f,pt_type)
              IF (ANY(lne_f_opt)) nvars = nvars + count_vars(lne_f_opt,ne_aux_f,ne_type)
              IF (ANY(lte_f_opt)) nvars = nvars + count_vars(lte_f_opt,te_aux_f,te_type)
              IF (ANY(lti_f_opt)) nvars = nvars + count_vars(lti_f_opt,ti_aux_f,ti_type)
              IF (ANY(lth_f_opt)) nvars = nvars + count_vars(lth_f_opt,th_aux_f,th_type)
              IF (ANY(lzeff_f_opt)) nvars = nvars + count_vars(lzeff_f_opt,zeff_aux_f,zeff_type)
              IF (ANY(lemis_xics_f_opt)) nvars = nvars + count_vars(lemis_xics_f_opt,emis_xics_f,emis_xics_type)
              IF (ANY(lphi_f_opt)) nvars = nvars + count_vars(lphi_f_opt,phi_aux_f,phi_type)
              IF (ANY(lbeamj_f_opt)) nvars = nvars + count_vars(lbeamj_f_opt,beamj_aux_f,beamj_type)
              IF (ANY(lbootj_f_opt)) nvars = nvars + count_vars(lbootj_f_opt,bootj_aux_f,bootj_type)
              ! Boundary Parameters
              DO n = 0, ntord
                 IF (laxis_opt(n)) THEN
                    nvars = nvars + 1
                    IF (n/=0) nvars = nvars + 1
                    !IF (lasym) THEN
                    !   nvars = nvars + 1
                    !   IF (n/=0) nvars = nvars + 1
                    !END IF
                 END IF
              END DO
              IF (lbound_opt(0,0)) THEN
                 nvars = nvars + 1
                 IF (lasym) nvars = nvars + 1
              END IF
              DO n = -ntord, ntord
                 DO m = 0, mpol1d
                    IF (m==0 .and. n<=0) CYCLE
                    IF (lbound_opt(n,m)) THEN
                       nvars = nvars + 2
                       IF (lasym) nvars = nvars + 2
                    END IF
                 END DO
              END DO
              DO n = -ntord, ntord
                 DO m = 0, mpol1d
                    IF (lmode_opt(n,m)) THEN
                       nvars = nvars + 1
                       IF (lasym) nvars = nvars + 1
                    END IF
                 END DO
              END DO
              DO n = -ntord, ntord
                 DO m = 0, mpol1d
                    IF (lrho_opt(n,m) .and. (m /= 0 .or. n >= 0)) THEN
                       nvars = nvars + 1
                    !   IF (lasym) nvars = nvars + 1
                    END IF
                 END DO
              END DO
              DO n = -ntord, ntord
                 DO m = -mpol1d, mpol1d
                    IF (ldeltamn_opt(n,m) .and. .not.(n == 0 .and. m==0)) THEN
                       nvars = nvars + 1
                    !   IF (lasym) nvars = nvars + 1
                    END IF
                 END DO
              END DO
              ctrl_dofs = 3                 !x,y,z at each point
              IF (lwindsurf) ctrl_dofs = 2  ! u,v  at each point
              DO n = LBOUND(lcoil_spline,DIM=1), UBOUND(lcoil_spline,DIM=1)
                 ! Actual no. of knots for coil spline n, less two to enforce periodicity of f,f'
                 nknots = COUNT(coil_splinesx(n,:) >= 0.0) - 2

                 ! First ctrl of modular loses one dof (u or z).
                 IF ((coil_type(n).eq.'M').AND.lcoil_spline(n,1)) nvars = nvars - 1

                 DO m = 1,nknots-4
                    IF (lcoil_spline(n,m)) nvars = nvars + ctrl_dofs
                 END DO
              END DO
              ier = 0

              ! REGCOIL options
              IF (lregcoil_winding_surface_separation_opt) nvars = nvars + 1
              IF (lregcoil_current_density_opt) nvars = nvars + 1
              DO m = -mpol_rcws, mpol_rcws
                 DO n = -ntor_rcws, ntor_rcws
                    IF (lregcoil_rcws_rbound_c_opt(m,n)) THEN
                       nvars = nvars + 1
                    END IF
                    IF (lregcoil_rcws_rbound_s_opt(m,n)) THEN
                       nvars = nvars + 1
                    END IF
                    IF (lregcoil_rcws_zbound_c_opt(m,n)) THEN
                       nvars = nvars + 1
                    END IF
                    IF (lregcoil_rcws_zbound_s_opt(m,n)) THEN
                       nvars = nvars + 1
                    END IF
                 END DO
              END DO
         CASE('spec')
         CASE('test')
            IF (lxval_opt)  nvars = nvars + 1
            IF (lyval_opt)  nvars = nvars + 1
            IF (ANY(lRosenbrock_X_opt)) nvars = nvars + COUNT(lRosenbrock_X_opt)
      END SELECT

      ! Allocate Arrays
      ALLOCATE(vars(1:nvars),var_dex(1:nvars),diag(1:nvars),&
               vars_min(1:nvars),vars_max(1:nvars),arr_dex(1:nvars,2))
      vars = 0.0; var_dex(:) = 0; diag = -1.0; arr_dex(:,:) = 0
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER( MPI_COMM_STEL, ierr_mpi )                   ! MPI
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BARRIER_ERR,'stellot_init',ierr_mpi)
!DEC$ ENDIF
      ! Read the Equilibrium Namelist and initalize the var arrays
      ! Initialize nvar_in to 0
      nvar_in=0
      SELECT CASE (TRIM(equil_type))
         CASE('vmec2000','animec','flow','satire','paravmec','parvmec','vboot','vmec2000_oneeq')
              ! Set some defaults
              phiedge_old = phiedge
              IF (ncurr /= 0 .and. ANY(lai_opt)) lai_opt(:) = .false.
              IF (ncurr /= 0 .and. ANY(lai_f_opt)) lai_f_opt(:) = .false.
              IF (ncurr /= 1 .and. lcurtor_opt) lcurtor_opt = .false.
              IF (ncurr /= 1 .and. ANY(lac_opt)) lac_opt(:) = .false.
              IF (ncurr /= 1 .and. ANY(lac_f_opt)) lac_f_opt(:) = .false.
              IF (ncurr /= 1 .and. ANY(lbeamj_f_opt)) lbeamj_f_opt(:) = .false.
              IF (ncurr /= 1 .and. ANY(lbootj_f_opt)) lbootj_f_opt(:) = .false.
              IF (.not. lfreeb) lno_restart = .true.
!DEC$ IF DEFINED (MPI_OPT)
              CALL MPI_BARRIER( MPI_COMM_STEL, ierr_mpi )                   ! MPI
              IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BARRIER_ERR,'stellot_init',ierr_mpi)
!DEC$ ENDIF
              ! Now count
              IF (lregcoil_winding_surface_separation_opt) THEN
                 IF (lauto_domain) THEN
                    regcoil_winding_surface_separation_min = &
                        regcoil_winding_surface_separation - &
                        ABS(pct_domain*regcoil_winding_surface_separation)
                    regcoil_winding_surface_separation_max = &
                        regcoil_winding_surface_separation + &
                        ABS(pct_domain*regcoil_winding_surface_separation)
                 END IF
                 nvar_in = nvar_in + 1
                 vars(nvar_in) = regcoil_winding_surface_separation
                 vars_min(nvar_in) = regcoil_winding_surface_separation_min
                 vars_max(nvar_in) = regcoil_winding_surface_separation_max
                 var_dex(nvar_in) = iregcoil_winding_surface_separation
                 diag(nvar_in)    = dregcoil_winding_surface_separation_opt
                 arr_dex(nvar_in,1) = 1
              END IF
              IF (lregcoil_current_density_opt) THEN
                 IF (lauto_domain) THEN
                    regcoil_current_density_min = &
                        regcoil_current_density - &
                        ABS(pct_domain*regcoil_current_density)
                    regcoil_current_density_max = &
                        regcoil_current_density + &
                        ABS(pct_domain*regcoil_current_density)
                 END IF
                 nvar_in = nvar_in + 1
                 vars(nvar_in) = regcoil_current_density
                 vars_min(nvar_in) = regcoil_current_density_min
                 vars_max(nvar_in) = regcoil_current_density_max
                 var_dex(nvar_in) = iregcoil_current_density
                 diag(nvar_in)    = dregcoil_current_density_opt
                 arr_dex(nvar_in,1) = 1
              END IF
              IF (ANY(lregcoil_rcws_rbound_c_opt) ) THEN
                 DO m = -mpol_rcws,mpol_rcws
                    DO n = -ntor_rcws,ntor_rcws
                       ! IF (m==0 .and. n<=0) CYCLE
                       IF (lregcoil_rcws_rbound_c_opt(m,n)) THEN
                          IF (lauto_domain) THEN
                             regcoil_rcws_rbound_c_min(m,n) = regcoil_rcws_rbound_c(m,n) - ABS(pct_domain*regcoil_rcws_rbound_c(m,n))
                             regcoil_rcws_rbound_c_max(m,n) = regcoil_rcws_rbound_c(m,n) + ABS(pct_domain*regcoil_rcws_rbound_c(m,n))
                          END IF
                          nvar_in = nvar_in + 1
                          vars(nvar_in) = regcoil_rcws_rbound_c(m,n)
                          vars_min(nvar_in) = regcoil_rcws_rbound_c_min(m,n)
                          vars_max(nvar_in) = regcoil_rcws_rbound_c_max(m,n)
                          var_dex(nvar_in) = iregcoil_rcws_rbound_c
                          diag(nvar_in)    = dregcoil_rcws_rbound_c_opt(m,n)
                          arr_dex(nvar_in,1) = m
                          arr_dex(nvar_in,2) = n
                       END IF
                    END DO
                 END DO
              END IF
              IF (ANY(lregcoil_rcws_rbound_s_opt) ) THEN
                 DO m = -mpol_rcws,mpol_rcws
                    DO n = -ntor_rcws,ntor_rcws
                       ! IF (m==0 .and. n<=0) CYCLE
                       IF (lregcoil_rcws_rbound_s_opt(m,n)) THEN
                          IF (lauto_domain) THEN
                             regcoil_rcws_rbound_s_min(m,n) = regcoil_rcws_rbound_s(m,n) - ABS(pct_domain*regcoil_rcws_rbound_s(m,n))
                             regcoil_rcws_rbound_s_max(m,n) = regcoil_rcws_rbound_s(m,n) + ABS(pct_domain*regcoil_rcws_rbound_s(m,n))
                          END IF
                          nvar_in = nvar_in + 1
                          vars(nvar_in) = regcoil_rcws_rbound_s(m,n)
                          vars_min(nvar_in) = regcoil_rcws_rbound_s_min(m,n)
                          vars_max(nvar_in) = regcoil_rcws_rbound_s_max(m,n)
                          var_dex(nvar_in) = iregcoil_rcws_rbound_s
                          diag(nvar_in)    = dregcoil_rcws_rbound_s_opt(m,n)
                          arr_dex(nvar_in,1) = m
                          arr_dex(nvar_in,2) = n
                       END IF
                    END DO
                 END DO
              END IF
              IF (ANY(lregcoil_rcws_zbound_c_opt) ) THEN
                 DO m = -mpol_rcws,mpol_rcws
                    DO n = -ntor_rcws,ntor_rcws
                       ! IF (m==0 .and. n<=0) CYCLE
                       IF (lregcoil_rcws_zbound_c_opt(m,n)) THEN
                          IF (lauto_domain) THEN
                             regcoil_rcws_zbound_c_min(m,n) = regcoil_rcws_zbound_c(m,n) - ABS(pct_domain*regcoil_rcws_zbound_c(m,n))
                             regcoil_rcws_zbound_c_max(m,n) = regcoil_rcws_zbound_c(m,n) + ABS(pct_domain*regcoil_rcws_zbound_c(m,n))
                          END IF
                          nvar_in = nvar_in + 1
                          vars(nvar_in) = regcoil_rcws_zbound_c(m,n)
                          vars_min(nvar_in) = regcoil_rcws_zbound_c_min(m,n)
                          vars_max(nvar_in) = regcoil_rcws_zbound_c_max(m,n)
                          var_dex(nvar_in) = iregcoil_rcws_zbound_c
                          diag(nvar_in)    = dregcoil_rcws_zbound_c_opt(m,n)
                          arr_dex(nvar_in,1) = m
                          arr_dex(nvar_in,2) = n
                       END IF
                    END DO
                 END DO
              END IF
              IF (ANY(lregcoil_rcws_zbound_s_opt) ) THEN
                 DO m = -mpol_rcws,mpol_rcws
                    DO n = -ntor_rcws,ntor_rcws
                       ! IF (m==0 .and. n<=0) CYCLE
                       IF (lregcoil_rcws_zbound_s_opt(m,n)) THEN
                          IF (lauto_domain) THEN
                             regcoil_rcws_zbound_s_min(m,n) = regcoil_rcws_zbound_s(m,n) - ABS(pct_domain*regcoil_rcws_zbound_s(m,n))
                             regcoil_rcws_zbound_s_max(m,n) = regcoil_rcws_zbound_s(m,n) + ABS(pct_domain*regcoil_rcws_zbound_s(m,n))
                          END IF
                          nvar_in = nvar_in + 1
                          vars(nvar_in) = regcoil_rcws_zbound_s(m,n)
                          vars_min(nvar_in) = regcoil_rcws_zbound_s_min(m,n)
                          vars_max(nvar_in) = regcoil_rcws_zbound_s_max(m,n)
                          var_dex(nvar_in) = iregcoil_rcws_zbound_s
                          diag(nvar_in)    = dregcoil_rcws_zbound_s_opt(m,n)
                          arr_dex(nvar_in,1) = m
                          arr_dex(nvar_in,2) = n
                       END IF
                    END DO
                 END DO
              END IF
              IF (lphiedge_opt) THEN
                 IF (lauto_domain) THEN
                    phiedge_min = phiedge - ABS(pct_domain*phiedge)
                    phiedge_max = phiedge + ABS(pct_domain*phiedge)
                 END IF
                 nvar_in = nvar_in + 1
                 vars(nvar_in) = phiedge
                 vars_min(nvar_in) = phiedge_min
                 vars_max(nvar_in) = phiedge_max
                 var_dex(nvar_in) = iphiedge
                 diag(nvar_in)    = dphiedge_opt
                 arr_dex(nvar_in,1) = 1
              END IF
              IF (lcurtor_opt) THEN
                 IF (lauto_domain) THEN
                    curtor_min = curtor - ABS(pct_domain*curtor)
                    curtor_max = curtor + ABS(pct_domain*curtor)
                 END IF
                 nvar_in = nvar_in + 1
                 vars(nvar_in) = curtor
                 vars_min(nvar_in) = curtor_min
                 vars_max(nvar_in) = curtor_max
                 var_dex(nvar_in) = icurtor
                 diag(nvar_in)    = dcurtor_opt
                 arr_dex(nvar_in,1) = 1
              END IF
              IF (lpscale_opt) THEN
                 IF (lauto_domain) THEN
                    pscale_min = pres_scale - ABS(pct_domain*pres_scale)
                    pscale_max = pres_scale + ABS(pct_domain*pres_scale)
                 END IF
                 nvar_in = nvar_in + 1
                 vars(nvar_in) = pres_scale
                 vars_min(nvar_in) = MAX(pscale_min,0.0_rprec)
                 vars_max(nvar_in) = MAX(pscale_max,0.0_rprec)
                 var_dex(nvar_in) = ipscale
                 diag(nvar_in)    = dpscale_opt
                 arr_dex(nvar_in,1) = 1
              END IF
              IF (lbcrit_opt) THEN
                 IF (lauto_domain) THEN
                    bcrit_min = bcrit - ABS(pct_domain*bcrit)
                    bcrit_max = bcrit + ABS(pct_domain*bcrit)
                 END IF
                 nvar_in = nvar_in + 1
                 vars(nvar_in) = bcrit
                 vars_min(nvar_in) = bcrit_min
                 vars_max(nvar_in) = bcrit_max
                 var_dex(nvar_in) = ibcrit
                 diag(nvar_in)    = dbcrit_opt
                 arr_dex(nvar_in,1) = 1
              END IF
              IF (lmix_ece_opt) THEN
                 IF (lauto_domain) THEN
                    mix_ece_min = MAX(mix_ece - ABS(pct_domain*mix_ece),0.0)
                    mix_ece_max = MIN(mix_ece + ABS(pct_domain*mix_ece),1.0)
                 END IF
                 nvar_in = nvar_in + 1
                 vars(nvar_in) = mix_ece
                 vars_min(nvar_in) = mix_ece_min
                 vars_max(nvar_in) = mix_ece_max
                 var_dex(nvar_in) = imixece
                 diag(nvar_in)    = dmix_ece_opt
                 arr_dex(nvar_in,1) = 1
              END IF
              IF (lxics_v0_opt) THEN
                 IF (lauto_domain) THEN
                    xics_v0_min = MAX(xics_v0 - ABS(pct_domain*xics_v0),5.0E3)
                    xics_v0_max = MIN(xics_v0 + ABS(pct_domain*xics_v0),5.0E3)
                 END IF
                 nvar_in = nvar_in + 1
                 vars(nvar_in) = xics_v0
                 vars_min(nvar_in) = xics_v0_min
                 vars_max(nvar_in) = xics_v0_max
                 var_dex(nvar_in) = ixics_v0
                 diag(nvar_in)    = dxics_v0_opt
                 arr_dex(nvar_in,1) = 1
              END IF
              IF (ANY(lextcur_opt)) THEN
                 DO i = LBOUND(lextcur_opt,DIM=1), UBOUND(lextcur_opt,DIM=1)
                    IF (lextcur_opt(i)) THEN
                       IF (lauto_domain) THEN
                          extcur_min(i) = extcur(i) - ABS(pct_domain*extcur(i))
                          extcur_max(i) = extcur(i) + ABS(pct_domain*extcur(i))
                       END IF
                       nvar_in = nvar_in + 1
                       vars(nvar_in) = extcur(i)
                       vars_min(nvar_in) = extcur_min(i)
                       vars_max(nvar_in) = extcur_max(i)
                       var_dex(nvar_in) = iextcur
                       diag(nvar_in)    = dextcur_opt(i)
                       arr_dex(nvar_in,1) = i
                    END IF
                 END DO
              END IF
              IF (ANY(laphi_opt)) THEN
                 norm = profile_norm(aphi,'power_series')
                 IF (norm /= 0) THEN
                    nvar_in = nvar_in + 1
                    vars(nvar_in) = norm
                    vars_min(nvar_in) = norm - abs(norm_fac*norm)
                    vars_max(nvar_in) = norm + abs(norm_fac*norm)
                    var_dex(nvar_in) = iaphi
                    diag(nvar_in)    = 1.0_rprec
                    arr_dex(nvar_in,2) = norm_dex
                    aphi = aphi / norm
                    aphi_min = aphi_min/norm
                    aphi_max = aphi_max/norm
                 END IF
                 DO i = LBOUND(laphi_opt,DIM=1), UBOUND(laphi_opt,DIM=1)
                    IF (laphi_opt(i)) THEN
                       IF (lauto_domain) THEN
                          aphi_min(i) = aphi(i) - ABS(pct_domain*aphi(i))
                          aphi_max(i) = aphi(i) + ABS(pct_domain*aphi(i))
                       END IF
                       nvar_in = nvar_in + 1
                       vars(nvar_in) = aphi(i)
                       vars_min(nvar_in) = aphi_min(i)
                       vars_max(nvar_in) = aphi_max(i)
                       var_dex(nvar_in) = iaphi
                       diag(nvar_in)    = daphi_opt(i)
                       arr_dex(nvar_in,1) = i
                    END IF
                 END DO
              END IF
              IF (ANY(lam_opt)) THEN
                 norm = profile_norm(am,pmass_type)
                 IF (norm /= 0) THEN
                    nvar_in = nvar_in + 1
                    vars(nvar_in) = norm
                    vars_min(nvar_in) = norm - abs(norm_fac*norm)
                    vars_max(nvar_in) = norm + abs(norm_fac*norm)
                    var_dex(nvar_in) = iam
                    diag(nvar_in)    = 1.0_rprec
                    arr_dex(nvar_in,2) = norm_dex
                    am = am / norm
                    am_min = am_min/norm
                    am_max = am_max/norm
                 END IF
                 DO i = LBOUND(lam_opt,DIM=1), UBOUND(lam_opt,DIM=1)
                    IF (lam_opt(i)) THEN
                       IF (lauto_domain) THEN
                          am_min(i) = am(i) - ABS(pct_domain*am(i))
                          am_max(i) = am(i) + ABS(pct_domain*am(i))
                       END IF
                       nvar_in = nvar_in + 1
                       vars(nvar_in) = am(i)
                       vars_min(nvar_in) = am_min(i)
                       vars_max(nvar_in) = am_max(i)
                       var_dex(nvar_in) = iam
                       diag(nvar_in)    = dam_opt(i)
                       arr_dex(nvar_in,1) = i
                    END IF
                 END DO
              END IF
              ! Quick note on AC, sometimes we don't want to normalize
              IF (ANY(lac_opt)) THEN
                 norm = profile_norm(ac,pcurr_type)
                 IF (norm /= 0) THEN
                    nvar_in = nvar_in + 1
                    vars(nvar_in) = norm
                    vars_min(nvar_in) = norm - abs(norm_fac*norm)
                    vars_max(nvar_in) = norm + abs(norm_fac*norm)
                    var_dex(nvar_in) = iac
                    diag(nvar_in)    = 1.0_rprec
                    arr_dex(nvar_in,2) = norm_dex
                    ac = ac / norm
                    ac_min = ac_min/norm
                    ac_max = ac_max/norm
                 END IF
                 DO i = LBOUND(lac_opt,DIM=1), UBOUND(lac_opt,DIM=1)
                    IF (lac_opt(i)) THEN
                       IF (lauto_domain) THEN
                          ac_min(i) = ac(i) - ABS(pct_domain*ac(i))
                          ac_max(i) = ac(i) + ABS(pct_domain*ac(i))
                       END IF
                       nvar_in = nvar_in + 1
                       vars(nvar_in) = ac(i)
                       vars_min(nvar_in) = ac_min(i)
                       vars_max(nvar_in) = ac_max(i)
                       var_dex(nvar_in) = iac
                       diag(nvar_in)    = dac_opt(i)
                       arr_dex(nvar_in,1) = i
                    END IF
                 END DO
              END IF
              IF (ANY(lai_opt)) THEN
                 norm = profile_norm(ai,piota_type)
                 IF (norm /= 0) THEN
                    nvar_in = nvar_in + 1
                    vars(nvar_in) = norm
                    vars_min(nvar_in) = norm - abs(norm_fac*norm)
                    vars_max(nvar_in) = norm + abs(norm_fac*norm)
                    var_dex(nvar_in) = iai
                    diag(nvar_in)    = 1.0_rprec
                    arr_dex(nvar_in,2) = norm_dex
                    ai = ai / norm
                    ai_min = ai_min/norm
                    ai_max = ai_max/norm
                 END IF
                 DO i = LBOUND(lai_opt,DIM=1), UBOUND(lai_opt,DIM=1)
                    IF (lai_opt(i)) THEN
                       IF (lauto_domain) THEN
                          ai_min(i) = ai(i) - ABS(pct_domain*ai(i))
                          ai_max(i) = ai(i) + ABS(pct_domain*ai(i))
                       END IF
                       nvar_in = nvar_in + 1
                       vars(nvar_in) = ai(i)
                       vars_min(nvar_in) = ai_min(i)
                       vars_max(nvar_in) = ai_max(i)
                       var_dex(nvar_in) = iai
                       diag(nvar_in)    = dai_opt(i)
                       arr_dex(nvar_in,1) = i
                    END IF
                 END DO
              END IF
              IF (ANY(lah_opt)) THEN
                 norm = profile_norm(ah,ph_type)
                 IF (norm /= 0) THEN
                    nvar_in = nvar_in + 1
                    vars(nvar_in) = norm
                    vars_min(nvar_in) = norm - abs(norm_fac*norm)
                    vars_max(nvar_in) = norm + abs(norm_fac*norm)
                    var_dex(nvar_in) = iah
                    diag(nvar_in)    = 1.0_rprec
                    arr_dex(nvar_in,2) = norm_dex
                    ah = ah / norm
                    ah_min = ah_min/norm
                    ah_max = ah_max/norm
                 END IF
                 DO i = LBOUND(lah_opt,DIM=1), UBOUND(lah_opt,DIM=1)
                    IF (lah_opt(i)) THEN
                       IF (lauto_domain) THEN
                          ah_min(i) = ah(i) - ABS(pct_domain*ah(i))
                          ah_max(i) = ah(i) + ABS(pct_domain*ah(i))
                       END IF
                       nvar_in = nvar_in + 1
                       vars(nvar_in) = ah(i)
                       vars_min(nvar_in) = ah_min(i)
                       vars_max(nvar_in) = ah_max(i)
                       var_dex(nvar_in) = iah
                       diag(nvar_in)    = dah_opt(i)
                       arr_dex(nvar_in,1) = i
                    END IF
                 END DO
              END IF
              IF (ANY(lat_opt)) THEN
                 norm = profile_norm(at,pt_type)
                 IF (norm /= 0) THEN
                    nvar_in = nvar_in + 1
                    vars(nvar_in) = norm
                    vars_min(nvar_in) = norm - abs(norm_fac*norm)
                    vars_max(nvar_in) = norm + abs(norm_fac*norm)
                    var_dex(nvar_in) = iat
                    diag(nvar_in)    = 1.0_rprec
                    arr_dex(nvar_in,2) = norm_dex
                    at = at / norm
                    at_min = at_min/norm
                    at_max = at_max/norm
                 END IF
                 DO i = LBOUND(lat_opt,DIM=1), UBOUND(lat_opt,DIM=1)
                    IF (lat_opt(i)) THEN
                       IF (lauto_domain) THEN
                          at_min(i) = at(i) - ABS(pct_domain*at(i))
                          at_max(i) = at(i) + ABS(pct_domain*at(i))
                       END IF
                       nvar_in = nvar_in + 1
                       vars(nvar_in) = at(i)
                       vars_min(nvar_in) = at_min(i)
                       vars_max(nvar_in) = at_max(i)
                       var_dex(nvar_in) = iat
                       diag(nvar_in)    = dat_opt(i)
                       arr_dex(nvar_in,1) = i
                    END IF
                 END DO
              END IF
              IF (ANY(lne_opt)) THEN
                 norm = profile_norm(ne_opt,ne_type)
                 IF (norm /= 0) THEN
                    nvar_in = nvar_in + 1
                    vars(nvar_in) = norm
                    vars_min(nvar_in) = norm - abs(norm_fac*norm)
                    vars_max(nvar_in) = norm + abs(norm_fac*norm)
                    var_dex(nvar_in) = ine
                    diag(nvar_in)    = 1.0_rprec
                    arr_dex(nvar_in,2) = norm_dex
                    ne_opt = ne_opt / norm
                    ne_min = ne_min/norm
                    ne_max = ne_max/norm
                 END IF
                 DO i = LBOUND(lne_opt,DIM=1), UBOUND(lne_opt,DIM=1)
                    IF (lne_opt(i)) THEN
                       IF (lauto_domain) THEN
                          ne_min(i) = ne_opt(i) - ABS(pct_domain*ne_opt(i))
                          ne_max(i) = ne_opt(i) + ABS(pct_domain*ne_opt(i))
                       END IF
                       nvar_in = nvar_in + 1
                       vars(nvar_in) = ne_opt(i)
                       vars_min(nvar_in) = ne_min(i)
                       vars_max(nvar_in) = ne_max(i)
                       var_dex(nvar_in) = ine
                       diag(nvar_in)    = dne_opt(i)
                       arr_dex(nvar_in,1) = i
                    END IF
                 END DO
              END IF
              
              IF (ANY(lzeff_opt)) THEN
                 norm = profile_norm(zeff_opt,zeff_type)
                 IF (norm /= 0) THEN
                    nvar_in = nvar_in + 1
                    vars(nvar_in) = norm
                    vars_min(nvar_in) = norm - abs(norm_fac*norm)
                    vars_max(nvar_in) = norm + abs(norm_fac*norm)
                    var_dex(nvar_in) = izeff
                    diag(nvar_in)    = 1.0_rprec
                    arr_dex(nvar_in,2) = norm_dex
                    zeff_opt = zeff_opt / norm
                    zeff_min = zeff_min/norm
                    zeff_max = zeff_max/norm
                 END IF
                 DO i = LBOUND(lzeff_opt,DIM=1), UBOUND(lzeff_opt,DIM=1)
                    IF (lzeff_opt(i)) THEN
                       IF (lauto_domain) THEN
                          zeff_min(i) = zeff_opt(i) - ABS(pct_domain*zeff_opt(i))
                          zeff_max(i) = zeff_opt(i) + ABS(pct_domain*zeff_opt(i))
                       END IF
                       nvar_in = nvar_in + 1
                       vars(nvar_in) = zeff_opt(i)
                       vars_min(nvar_in) = zeff_min(i)
                       vars_max(nvar_in) = zeff_max(i)
                       var_dex(nvar_in) = izeff
                       diag(nvar_in)    = dzeff_opt(i)
                       arr_dex(nvar_in,1) = i
                    END IF
                 END DO
              END IF
              
              IF (ANY(lte_opt)) THEN
                 norm = profile_norm(te_opt,te_type)
                 IF (norm /= 0.0) THEN
                    nvar_in = nvar_in + 1
                    vars(nvar_in) = norm
                    vars_min(nvar_in) = norm - abs(norm_fac*norm)
                    vars_max(nvar_in) = norm + abs(norm_fac*norm)
                    var_dex(nvar_in) = ite
                    diag(nvar_in)    = 1.0_rprec
                    arr_dex(nvar_in,2) = norm_dex
                    te_opt = te_opt / norm
                    te_min = te_min/norm
                    te_max = te_max/norm
                 END IF
                 DO i = LBOUND(lte_opt,DIM=1), UBOUND(lte_opt,DIM=1)
                    IF (lte_opt(i)) THEN
                       IF (lauto_domain) THEN
                          te_min(i) = te_opt(i) - ABS(pct_domain*te_opt(i))
                          te_max(i) = te_opt(i) + ABS(pct_domain*te_opt(i))
                       END IF
                       nvar_in = nvar_in + 1
                       vars(nvar_in) = te_opt(i)
                       vars_min(nvar_in) = te_min(i)
                       vars_max(nvar_in) = te_max(i)
                       var_dex(nvar_in) = ite
                       diag(nvar_in)    = dte_opt(i)
                       arr_dex(nvar_in,1) = i
                    END IF
                 END DO
              END IF
              IF (ANY(lti_opt)) THEN
                 norm = profile_norm(ti_opt,ti_type)
                 IF (norm /= 0) THEN
                    nvar_in = nvar_in + 1
                    vars(nvar_in) = norm
                    vars_min(nvar_in) = norm - abs(norm_fac*norm)
                    vars_max(nvar_in) = norm + abs(norm_fac*norm)
                    var_dex(nvar_in) = iti
                    diag(nvar_in)    = 1.0_rprec
                    arr_dex(nvar_in,2) = norm_dex
                    ti_opt = ti_opt / norm
                    ti_min = ti_min/norm
                    ti_max = ti_max/norm
                 END IF
                 DO i = LBOUND(lti_opt,DIM=1), UBOUND(lti_opt,DIM=1)
                    IF (lti_opt(i)) THEN
                       IF (lauto_domain) THEN
                          ti_min(i) = ti_opt(i) - ABS(pct_domain*ti_opt(i))
                          ti_max(i) = ti_opt(i) + ABS(pct_domain*ti_opt(i))
                       END IF
                       nvar_in = nvar_in + 1
                       vars(nvar_in) = ti_opt(i)
                       vars_min(nvar_in) = ti_min(i)
                       vars_max(nvar_in) = ti_max(i)
                       var_dex(nvar_in) = iti
                       diag(nvar_in)    = dti_opt(i)
                       arr_dex(nvar_in,1) = i
                    END IF
                 END DO
              END IF
              IF (ANY(lth_opt)) THEN
                 norm = profile_norm(th_opt,th_type)
                 IF (norm /=0) THEN
                    nvar_in = nvar_in + 1
                    vars(nvar_in) = norm
                    vars_min(nvar_in) = norm - abs(norm_fac*norm)
                    vars_max(nvar_in) = norm + abs(norm_fac*norm)
                    var_dex(nvar_in) = ith
                    diag(nvar_in)    = 1.0_rprec
                    arr_dex(nvar_in,2) = norm_dex
                    th_opt = th_opt / norm
                    th_min = th_min/norm
                    th_max = th_max/norm
                 END IF
                 DO i = LBOUND(lth_opt,DIM=1), UBOUND(lth_opt,DIM=1)
                    IF (lth_opt(i)) THEN
                       IF (lauto_domain) THEN
                          th_min(i) = th_opt(i) - ABS(pct_domain*th_opt(i))
                          th_max(i) = th_opt(i) + ABS(pct_domain*th_opt(i))
                       END IF
                       nvar_in = nvar_in + 1
                       vars(nvar_in) = th_opt(i)
                       vars_min(nvar_in) = th_min(i)
                       vars_max(nvar_in) = th_max(i)
                       var_dex(nvar_in) = ith
                       diag(nvar_in)    = dth_opt(i)
                       arr_dex(nvar_in,1) = i
                    END IF
                 END DO
              END IF
              
              IF (ANY(lam_s_opt)) THEN
                 DO i = LBOUND(lam_s_opt,DIM=1), UBOUND(lam_s_opt,DIM=1)
                    IF (lam_s_opt(i)) THEN
                       nvar_in = nvar_in + 1
                       vars(nvar_in) = am_aux_s(i)
                       vars_min(nvar_in) = 0.0_rprec
                       vars_max(nvar_in) = 1.0_rprec
                       var_dex(nvar_in) = iam_aux_s
                       diag(nvar_in)    = dam_s_opt(i)
                       arr_dex(nvar_in,1) = i
                    END IF
                 END DO
              END IF
              IF (ANY(lam_f_opt)) THEN
                 norm = profile_norm(am_aux_f,pmass_type)
                 IF (norm /=0) THEN
                    nvar_in = nvar_in + 1
                    vars(nvar_in) = norm
                    vars_min(nvar_in) = norm - abs(norm_fac*norm)
                    vars_max(nvar_in) = norm + abs(norm_fac*norm)
                    var_dex(nvar_in) = iam_aux_f
                    diag(nvar_in)    = 1.0_rprec
                    arr_dex(nvar_in,2) = norm_dex
                    am_aux_f = am_aux_f / norm
                    am_f_min = am_f_min/norm
                    am_f_max = am_f_max/norm
                 END IF
                 DO i = LBOUND(lam_f_opt,DIM=1), UBOUND(lam_f_opt,DIM=1)
                    IF (lam_f_opt(i)) THEN
                       IF (lauto_domain) THEN
                          am_f_min(i) = am_aux_f(i) - ABS(pct_domain*am_aux_f(i))
                          am_f_max(i) = am_aux_f(i) + ABS(pct_domain*am_aux_f(i))
                       END IF
                       nvar_in = nvar_in + 1
                       vars(nvar_in) = am_aux_f(i)
                       vars_min(nvar_in) = am_f_min(i)
                       vars_max(nvar_in) = am_f_max(i)
                       var_dex(nvar_in) = iam_aux_f
                       diag(nvar_in)    = dam_f_opt(i)
                       arr_dex(nvar_in,1) = i
                    END IF
                 END DO
              END IF
              IF (ANY(lac_s_opt)) THEN
                 DO i = LBOUND(lac_s_opt,DIM=1), UBOUND(lac_s_opt,DIM=1)
                    IF (lac_s_opt(i)) THEN
                       nvar_in = nvar_in + 1
                       vars(nvar_in) = ac_aux_s(i)
                       vars_min(nvar_in) = 0.0_rprec
                       vars_max(nvar_in) = 1.0_rprec
                       var_dex(nvar_in) = iac_aux_s
                       diag(nvar_in)    = dac_s_opt(i)
                       arr_dex(nvar_in,1) = i
                    END IF
                 END DO
              END IF
              IF (ANY(lac_f_opt)) THEN
                 norm = profile_norm(ac_aux_f,pcurr_type)
                 IF (norm /=0) THEN
                    nvar_in = nvar_in + 1
                    vars(nvar_in) = norm
                    vars_min(nvar_in) = norm - abs(norm_fac*norm)
                    vars_max(nvar_in) = norm + abs(norm_fac*norm)
                    var_dex(nvar_in) = iac_aux_f
                    diag(nvar_in)    = 1.0_rprec
                    arr_dex(nvar_in,2) = norm_dex
                    ac_aux_f = ac_aux_f / norm
                    ac_f_min = ac_f_min/norm
                    ac_f_max = ac_f_max/norm
                 END IF
                 DO i = LBOUND(lac_f_opt,DIM=1), UBOUND(lac_f_opt,DIM=1)
                    IF (lac_f_opt(i)) THEN
                       IF (lauto_domain) THEN
                          ac_f_min(i) = ac_aux_f(i) - ABS(pct_domain*ac_aux_f(i))
                          ac_f_max(i) = ac_aux_f(i) + ABS(pct_domain*ac_aux_f(i))
                       END IF
                       nvar_in = nvar_in + 1
                       vars(nvar_in) = ac_aux_f(i)
                       vars_min(nvar_in) = ac_f_min(i)
                       vars_max(nvar_in) = ac_f_max(i)
                       var_dex(nvar_in) = iac_aux_f
                       diag(nvar_in)    = dac_f_opt(i)
                       arr_dex(nvar_in,1) = i
                    END IF
                 END DO
              END IF
              IF (ANY(lbeamj_f_opt)) THEN
                 norm = profile_norm(beamj_aux_f,beamj_type)
                 IF (norm /=0) THEN
                    nvar_in = nvar_in + 1
                    vars(nvar_in) = norm
                    vars_min(nvar_in) = norm - abs(norm_fac*norm)
                    vars_max(nvar_in) = norm + abs(norm_fac*norm)
                    var_dex(nvar_in) = ibeamj_aux_f
                    diag(nvar_in)    = 1.0_rprec
                    arr_dex(nvar_in,2) = norm_dex
                    beamj_aux_f = beamj_aux_f / norm
                    beamj_f_min = beamj_f_min/norm
                    beamj_f_max = beamj_f_max/norm
                 END IF
                 DO i = LBOUND(lbeamj_f_opt,DIM=1), UBOUND(lbeamj_f_opt,DIM=1)
                    IF (lbeamj_f_opt(i)) THEN
                       IF (lauto_domain) THEN
                          beamj_f_min(i) = beamj_aux_f(i) - ABS(pct_domain*beamj_aux_f(i))
                          beamj_f_max(i) = beamj_aux_f(i) + ABS(pct_domain*beamj_aux_f(i))
                       END IF
                       nvar_in = nvar_in + 1
                       vars(nvar_in) = beamj_aux_f(i)
                       vars_min(nvar_in) = beamj_f_min(i)
                       vars_max(nvar_in) = beamj_f_max(i)
                       var_dex(nvar_in) = ibeamj_aux_f
                       diag(nvar_in)    = dbeamj_f_opt(i)
                       arr_dex(nvar_in,1) = i
                    END IF
                 END DO
              END IF
              IF (ANY(lbootj_f_opt)) THEN
                 norm = profile_norm(bootj_aux_f,bootj_type)
                 IF (norm /=0) THEN
                    nvar_in = nvar_in + 1
                    vars(nvar_in) = norm
                    vars_min(nvar_in) = norm - abs(norm_fac*norm)
                    vars_max(nvar_in) = norm + abs(norm_fac*norm)
                    var_dex(nvar_in) = ibootj_aux_f
                    diag(nvar_in)    = 1.0_rprec
                    arr_dex(nvar_in,2) = norm_dex
                    bootj_aux_f = bootj_aux_f / norm
                    bootj_f_min = bootj_f_min/norm
                    bootj_f_max = bootj_f_max/norm
                 END IF
                 DO i = LBOUND(lbootj_f_opt,DIM=1), UBOUND(lbootj_f_opt,DIM=1)
                    IF (lbootj_f_opt(i)) THEN
                       IF (lauto_domain) THEN
                          bootj_f_min(i) = bootj_aux_f(i) - ABS(pct_domain*bootj_aux_f(i))
                          bootj_f_max(i) = bootj_aux_f(i) + ABS(pct_domain*bootj_aux_f(i))
                       END IF
                       nvar_in = nvar_in + 1
                       vars(nvar_in) = bootj_aux_f(i)
                       vars_min(nvar_in) = bootj_f_min(i)
                       vars_max(nvar_in) = bootj_f_max(i)
                       var_dex(nvar_in) = ibootj_aux_f
                       diag(nvar_in)    = dbootj_f_opt(i)
                       arr_dex(nvar_in,1) = i
                    END IF
                 END DO
              END IF
              IF (ANY(lemis_xics_f_opt)) THEN
                 norm = profile_norm(emis_xics_f,emis_xics_type)
                 IF (norm /=0) THEN
                    nvar_in = nvar_in + 1
                    vars(nvar_in) = norm
                    vars_min(nvar_in) = norm - abs(norm_fac*norm)
                    vars_max(nvar_in) = norm + abs(norm_fac*norm)
                    var_dex(nvar_in) = iemis_xics_f
                    diag(nvar_in)    = 1.0_rprec
                    arr_dex(nvar_in,2) = norm_dex
                    emis_xics_f = emis_xics_f / norm
                    emis_xics_f_min = emis_xics_f_min/norm
                    emis_xics_f_max = emis_xics_f_max/norm
                 END IF
                 DO i = LBOUND(lemis_xics_f_opt,DIM=1), UBOUND(lemis_xics_f_opt,DIM=1)
                    IF (lemis_xics_f_opt(i)) THEN
                       IF (lauto_domain) THEN
                          emis_xics_f_min(i) = emis_xics_f(i) - ABS(pct_domain*emis_xics_f(i))
                          emis_xics_f_max(i) = emis_xics_f(i) + ABS(pct_domain*emis_xics_f(i))
                       END IF
                       nvar_in = nvar_in + 1
                       vars(nvar_in) = emis_xics_f(i)
                       vars_min(nvar_in) = emis_xics_f_min(i)
                       vars_max(nvar_in) = emis_xics_f_max(i)
                       var_dex(nvar_in) = iemis_xics_f
                       diag(nvar_in)    = demis_xics_f_opt(i)
                       arr_dex(nvar_in,1) = i
                    END IF
                 END DO
              END IF
              IF (ANY(lai_s_opt)) THEN
                 DO i = LBOUND(lai_s_opt,DIM=1), UBOUND(lai_s_opt,DIM=1)
                    IF (lai_s_opt(i)) THEN
                       nvar_in = nvar_in + 1
                       vars(nvar_in) = ai_aux_s(i)
                       vars_min(nvar_in) = 0.0_rprec
                       vars_max(nvar_in) = 1.0_rprec
                       var_dex(nvar_in) = iai_aux_s
                       diag(nvar_in)    = dai_s_opt(i)
                       arr_dex(nvar_in,1) = i
                    END IF
                 END DO
              END IF
              IF (ANY(lai_f_opt)) THEN
                 norm = profile_norm(ai_aux_f,piota_type)
                 IF (norm /=0) THEN
                    nvar_in = nvar_in + 1
                    vars(nvar_in) = norm
                    vars_min(nvar_in) = norm - abs(norm_fac*norm)
                    vars_max(nvar_in) = norm + abs(norm_fac*norm)
                    var_dex(nvar_in) = iai_aux_f
                    diag(nvar_in)    = 1.0_rprec
                    arr_dex(nvar_in,2) = norm_dex
                    ai_aux_f = ai_aux_f / norm
                    ai_f_min = ai_f_min/norm
                    ai_f_max = ai_f_max/norm
                 END IF
                 DO i = LBOUND(lai_f_opt,DIM=1), UBOUND(lai_f_opt,DIM=1)
                    IF (lai_f_opt(i)) THEN
                       IF (lauto_domain) THEN
                          ai_f_min(i) = ai_aux_f(i) - ABS(pct_domain*ai_aux_f(i))
                          ai_f_max(i) = ai_aux_f(i) + ABS(pct_domain*ai_aux_f(i))
                       END IF
                       nvar_in = nvar_in + 1
                       vars(nvar_in) = ai_aux_f(i)
                       vars_min(nvar_in) = ai_f_min(i)
                       vars_max(nvar_in) = ai_f_max(i)
                       var_dex(nvar_in) = iai_aux_f
                       diag(nvar_in)    = dai_f_opt(i)
                       arr_dex(nvar_in,1) = i
                    END IF
                 END DO
              END IF
              IF (ANY(lphi_s_opt)) THEN
                 DO i = LBOUND(lphi_s_opt,DIM=1), UBOUND(lphi_s_opt,DIM=1)
                    IF (lphi_s_opt(i)) THEN
                       nvar_in = nvar_in + 1
                       vars(nvar_in) = phi_aux_s(i)
                       vars_min(nvar_in) = 0.0_rprec
                       vars_max(nvar_in) = 1.0_rprec
                       var_dex(nvar_in) = iphi_aux_s
                       diag(nvar_in)    = dphi_s_opt(i)
                       arr_dex(nvar_in,1) = i
                    END IF
                 END DO
              END IF
              IF (ANY(lphi_f_opt)) THEN
                 norm = profile_norm(phi_aux_f,phi_type)
                 IF (norm /=0) THEN
                    nvar_in = nvar_in + 1
                    vars(nvar_in) = norm
                    vars_min(nvar_in) = norm - abs(norm_fac*norm)
                    vars_max(nvar_in) = norm + abs(norm_fac*norm)
                    var_dex(nvar_in) = iphi_aux_f
                    diag(nvar_in)    = 1.0_rprec
                    arr_dex(nvar_in,2) = norm_dex
                    phi_aux_f = phi_aux_f / norm
                    phi_f_min = phi_f_min/norm
                    phi_f_max = phi_f_max/norm
                 END IF
                 DO i = LBOUND(lphi_f_opt,DIM=1), UBOUND(lphi_f_opt,DIM=1)
                    IF (lphi_f_opt(i)) THEN
                       IF (lauto_domain) THEN
                          phi_f_min(i) = phi_aux_f(i) - ABS(pct_domain*phi_aux_f(i))
                          phi_f_max(i) = phi_aux_f(i) + ABS(pct_domain*phi_aux_f(i))
                       END IF
                       nvar_in = nvar_in + 1
                       vars(nvar_in) = phi_aux_f(i)
                       vars_min(nvar_in) = phi_f_min(i)
                       vars_max(nvar_in) = phi_f_max(i)
                       var_dex(nvar_in) = iphi_aux_f
                       diag(nvar_in)    = dphi_f_opt(i)
                       arr_dex(nvar_in,1) = i
                    END IF
                 END DO
              END IF
              IF (ANY(lne_f_opt)) THEN
                 norm = profile_norm(ne_aux_f,ne_type)
                 IF (norm /=0) THEN
                    nvar_in = nvar_in + 1
                    vars(nvar_in) = norm
                    vars_min(nvar_in) = norm - abs(norm_fac*norm)
                    vars_max(nvar_in) = norm + abs(norm_fac*norm)
                    var_dex(nvar_in) = ine_aux_f
                    diag(nvar_in)    = 1.0_rprec
                    arr_dex(nvar_in,2) = norm_dex
                    ne_aux_f = ne_aux_f / norm
                    ne_f_min = ne_f_min/norm
                    ne_f_max = ne_f_max/norm
                 END IF
                 DO i = LBOUND(lne_f_opt,DIM=1), UBOUND(lne_f_opt,DIM=1)
                    IF (lne_f_opt(i)) THEN
                       IF (lauto_domain) THEN
                          ne_f_min(i) = ne_aux_f(i) - ABS(pct_domain*ne_aux_f(i))
                          ne_f_max(i) = ne_aux_f(i) + ABS(pct_domain*ne_aux_f(i))
                       END IF
                       nvar_in = nvar_in + 1
                       vars(nvar_in) = ne_aux_f(i)
                       vars_min(nvar_in) = MAX(ne_f_min(i),0.0_rprec)
                       vars_max(nvar_in) = ne_f_max(i)
                       var_dex(nvar_in) = ine_aux_f
                       diag(nvar_in)    = dne_f_opt(i)
                       arr_dex(nvar_in,1) = i
                    END IF
                 END DO
              END IF
              IF (ANY(lzeff_f_opt)) THEN
                 norm = profile_norm(zeff_aux_f,zeff_type)
                 IF (norm /=0) THEN
                    nvar_in = nvar_in + 1
                    vars(nvar_in) = norm
                    vars_min(nvar_in) = norm - abs(norm_fac*norm)
                    vars_max(nvar_in) = norm + abs(norm_fac*norm)
                    var_dex(nvar_in) = izeff_aux_f
                    diag(nvar_in)    = 1.0_rprec
                    arr_dex(nvar_in,2) = norm_dex
                    zeff_aux_f = zeff_aux_f / norm
                    zeff_f_min = zeff_f_min/norm
                    zeff_f_max = zeff_f_max/norm
                 END IF
                 DO i = LBOUND(lzeff_f_opt,DIM=1), UBOUND(lzeff_f_opt,DIM=1)
                    IF (lzeff_f_opt(i)) THEN
                       IF (lauto_domain) THEN
                          zeff_f_min(i) = zeff_aux_f(i) - ABS(pct_domain*zeff_aux_f(i))
                          zeff_f_max(i) = zeff_aux_f(i) + ABS(pct_domain*zeff_aux_f(i))
                       END IF
                       nvar_in = nvar_in + 1
                       vars(nvar_in) = zeff_aux_f(i)
                       vars_min(nvar_in) = MAX(zeff_f_min(i),0.0_rprec)
                       vars_max(nvar_in) = zeff_f_max(i)
                       var_dex(nvar_in) = izeff_aux_f
                       diag(nvar_in)    = dzeff_f_opt(i)
                       arr_dex(nvar_in,1) = i
                    END IF
                 END DO
              END IF
              IF (ANY(lte_f_opt)) THEN
                 norm = profile_norm(te_aux_f,te_type)
                 IF (norm /=0) THEN
                    nvar_in = nvar_in + 1
                    vars(nvar_in) = norm
                    vars_min(nvar_in) = norm - abs(norm_fac*norm)
                    vars_max(nvar_in) = norm + abs(norm_fac*norm)
                    var_dex(nvar_in) = ite_aux_f
                    diag(nvar_in)    = 1.0_rprec
                    arr_dex(nvar_in,2) = norm_dex
                    te_aux_f = te_aux_f / norm
                    te_f_min = te_f_min/norm
                    te_f_max = te_f_max/norm
                 END IF
                 DO i = LBOUND(lte_f_opt,DIM=1), UBOUND(lte_f_opt,DIM=1)
                    IF (lte_f_opt(i)) THEN
                       IF (lauto_domain) THEN
                          te_f_min(i) = te_aux_f(i) - ABS(pct_domain*te_aux_f(i))
                          te_f_max(i) = te_aux_f(i) + ABS(pct_domain*te_aux_f(i))
                       END IF
                       nvar_in = nvar_in + 1
                       vars(nvar_in) = te_aux_f(i)
                       vars_min(nvar_in) = MAX(te_f_min(i),0.0_rprec)
                       vars_max(nvar_in) = te_f_max(i)
                       var_dex(nvar_in) = ite_aux_f
                       diag(nvar_in)    = dte_f_opt(i)
                       arr_dex(nvar_in,1) = i
                    END IF
                 END DO
              END IF
              IF (ANY(lti_f_opt)) THEN
                 norm = profile_norm(ti_aux_f,ti_type)
                 IF (norm /=0) THEN
                    nvar_in = nvar_in + 1
                    vars(nvar_in) = norm
                    vars_min(nvar_in) = norm - abs(norm_fac*norm)
                    vars_max(nvar_in) = norm + abs(norm_fac*norm)
                    var_dex(nvar_in) = iti_aux_f
                    diag(nvar_in)    = 1.0_rprec
                    arr_dex(nvar_in,2) = norm_dex
                    ti_aux_f = ti_aux_f / norm
                    ti_f_min = ti_f_min/norm
                    ti_f_max = ti_f_max/norm
                 END IF
                 DO i = LBOUND(lti_f_opt,DIM=1), UBOUND(lti_f_opt,DIM=1)
                    IF (lti_f_opt(i)) THEN
                       IF (lauto_domain) THEN
                          ti_f_min(i) = ti_aux_f(i) - ABS(pct_domain*ti_aux_f(i))
                          ti_f_max(i) = ti_aux_f(i) + ABS(pct_domain*ti_aux_f(i))
                       END IF
                       nvar_in = nvar_in + 1
                       vars(nvar_in) = ti_aux_f(i)
                       vars_min(nvar_in) = MAX(ti_f_min(i),0.0_rprec)
                       vars_max(nvar_in) = ti_f_max(i)
                       var_dex(nvar_in) = iti_aux_f
                       diag(nvar_in)    = dti_f_opt(i)
                       arr_dex(nvar_in,1) = i
                    END IF
                 END DO
              END IF
              IF (ANY(lth_f_opt)) THEN
                 norm = profile_norm(th_aux_f,th_type)
                 IF (norm /=0) THEN
                    nvar_in = nvar_in + 1
                    vars(nvar_in) = norm
                    vars_min(nvar_in) = norm - abs(norm_fac*norm)
                    vars_max(nvar_in) = norm + abs(norm_fac*norm)
                    var_dex(nvar_in) = ith_aux_f
                    diag(nvar_in)    = 1.0_rprec
                    arr_dex(nvar_in,2) = norm_dex
                    th_aux_f = th_aux_f / norm
                    th_f_min = th_f_min/norm
                    th_f_max = th_f_max/norm
                 END IF
                 DO i = LBOUND(lth_f_opt,DIM=1), UBOUND(lth_f_opt,DIM=1)
                    IF (lth_f_opt(i)) THEN
                       IF (lauto_domain) THEN
                          th_f_min(i) = th_aux_f(i) - ABS(pct_domain*th_aux_f(i))
                          th_f_max(i) = th_aux_f(i) + ABS(pct_domain*th_aux_f(i))
                       END IF
                       nvar_in = nvar_in + 1
                       vars(nvar_in) = th_aux_f(i)
                       vars_min(nvar_in) = MAX(th_f_min(i),0.0_rprec)
                       vars_max(nvar_in) = th_f_max(i)
                       var_dex(nvar_in) = ith_aux_f
                       diag(nvar_in)    = dth_f_opt(i)
                       arr_dex(nvar_in,1) = i
                    END IF
                 END DO
              END IF
              IF (ANY(lah_f_opt)) THEN
                 norm = profile_norm(ah_aux_f,ph_type)
                 IF (norm /=0) THEN
                    nvar_in = nvar_in + 1
                    vars(nvar_in) = norm
                    vars_min(nvar_in) = norm - abs(norm_fac*norm)
                    vars_max(nvar_in) = norm + abs(norm_fac*norm)
                    var_dex(nvar_in) = iah_aux_f
                    diag(nvar_in)    = 1.0_rprec
                    arr_dex(nvar_in,2) = norm_dex
                    ah_aux_f = ah_aux_f / norm
                    ah_f_min = ah_f_min/norm
                    ah_f_max = ah_f_max/norm
                 END IF
                 DO i = LBOUND(lah_f_opt,DIM=1), UBOUND(lah_f_opt,DIM=1)
                    IF (lah_f_opt(i)) THEN
                       IF (lauto_domain) THEN
                          ah_f_min(i) = ah_aux_f(i) - ABS(pct_domain*ah_aux_f(i))
                          ah_f_max(i) = ah_aux_f(i) + ABS(pct_domain*ah_aux_f(i))
                       END IF
                       nvar_in = nvar_in + 1
                       vars(nvar_in) = ah_aux_f(i)
                       vars_min(nvar_in) = ah_f_min(i)
                       vars_max(nvar_in) = ah_f_max(i)
                       var_dex(nvar_in) = iah_aux_f
                       diag(nvar_in)    = dah_f_opt(i)
                       arr_dex(nvar_in,1) = i
                    END IF
                 END DO
              END IF
              IF (ANY(lat_f_opt)) THEN
                 norm = profile_norm(at_aux_f,pt_type)
                 IF (norm /=0) THEN
                    nvar_in = nvar_in + 1
                    vars(nvar_in) = norm
                    vars_min(nvar_in) = norm - abs(norm_fac*norm)
                    vars_max(nvar_in) = norm + abs(norm_fac*norm)
                    var_dex(nvar_in) = iat_aux_f
                    diag(nvar_in)    = 1.0_rprec
                    arr_dex(nvar_in,2) = norm_dex
                    at_aux_f = at_aux_f / norm
                    at_f_min = at_f_min/norm
                    at_f_max = at_f_max/norm
                 END IF
                 DO i = LBOUND(lat_f_opt,DIM=1), UBOUND(lat_f_opt,DIM=1)
                    IF (lat_f_opt(i)) THEN
                       IF (lauto_domain) THEN
                          at_f_min(i) = at_aux_f(i) - ABS(pct_domain*at_aux_f(i))
                          at_f_max(i) = at_aux_f(i) + ABS(pct_domain*at_aux_f(i))
                       END IF
                       nvar_in = nvar_in + 1
                       vars(nvar_in) = ah_aux_f(i)
                       vars_min(nvar_in) = at_f_min(i)
                       vars_max(nvar_in) = at_f_max(i)
                       var_dex(nvar_in) = iat_aux_f
                       diag(nvar_in)    = dat_f_opt(i)
                       arr_dex(nvar_in,1) = i
                    END IF
                 END DO
              END IF
              IF (ANY(laxis_opt)) THEN
                 n = UBOUND(laxis_opt,1)
                 DO n = LBOUND(laxis_opt,1), UBOUND(laxis_opt,1)
                    IF (laxis_opt(n)) THEN
                       nvar_in = nvar_in + 1
                       vars(nvar_in) = raxis_cc(n)
                       IF (lauto_domain) THEN
                          norm = MIN(raxis_cc(n),raxis_cs(n))
                          raxis_min(n) = norm - ABS(pct_domain*norm)
                          norm = MAX(raxis_cc(n),raxis_cs(n))
                          raxis_max(n) = norm + ABS(pct_domain*norm)
                       END IF
                       vars_min(nvar_in) = raxis_min(n)
                       vars_max(nvar_in) = raxis_max(n)
                       var_dex(nvar_in)  = iraxis_cc
                       diag(nvar_in)     = daxis_opt(n)
                       arr_dex(nvar_in,1) = n
                       IF (n /= 0) THEN
                          nvar_in = nvar_in + 1
                          vars(nvar_in) = zaxis_cs(n)
                          IF (lauto_domain) THEN
                             norm = MIN(zaxis_cs(n),zaxis_cc(n))
                             zaxis_min(n) = norm - ABS(pct_domain*norm)
                             norm = MAX(zaxis_cs(n),zaxis_cc(n))
                             zaxis_max(n) = norm + ABS(pct_domain*norm)
                          END IF
                          vars_min(nvar_in) = zaxis_min(n)
                          vars_max(nvar_in) = zaxis_max(n)
                          var_dex(nvar_in)  = izaxis_cs
                          diag(nvar_in)     = daxis_opt(n)
                          arr_dex(nvar_in,1) = n
                       END IF
                       IF (lasym) THEN
                          IF (n /= 0) THEN
                             nvar_in = nvar_in + 1
                             vars(nvar_in) = raxis_cs(n)
                             vars_min(nvar_in) = raxis_min(n)
                             vars_max(nvar_in) = raxis_max(n)
                             var_dex(nvar_in)  = iraxis_cs
                             diag(nvar_in)     = daxis_opt(n)
                             arr_dex(nvar_in,1) = n
                          END IF
                          nvar_in = nvar_in + 1
                          vars(nvar_in) = zaxis_cc(n)
                          vars_min(nvar_in) = zaxis_min(n)
                          vars_max(nvar_in) = zaxis_max(n)
                          var_dex(nvar_in)  = izaxis_cc
                          diag(nvar_in)     = daxis_opt(n)
                          arr_dex(nvar_in,1) = n
                       END IF
                    END IF
                 END DO
              END IF
              IF (ANY(lmode_opt)) THEN
                 DO n = LBOUND(lmode_opt,1), UBOUND(lmode_opt,1)
                    DO m = LBOUND(lmode_opt,2), UBOUND(lmode_opt,2)
                       IF (lmode_opt(n,m)) THEN
                          nvar_in = nvar_in + 1
                          vars(nvar_in) = 0.5*(rbc(n,m)+zbs(n,m))
                          IF (lauto_domain) THEN
                             bound_min(n,m) = vars(nvar_in) - ABS(pct_domain*vars(nvar_in))
                             bound_max(n,m) = vars(nvar_in) + ABS(pct_domain*vars(nvar_in))
                          END IF
                          vars_min(nvar_in) = bound_min(n,m)
                          vars_max(nvar_in) = bound_max(n,m)
                          var_dex(nvar_in)  = imodemn
                          diag(nvar_in)     = dbound_opt(n,m)
                          arr_dex(nvar_in,1) = n
                          arr_dex(nvar_in,2) = m
                       END IF
                    END DO
                 END DO
              END IF
              IF (ANY(lrho_opt)) THEN
                 m = UBOUND(lrho_opt,2)
                 n = UBOUND(lrho_opt,1)
                 rhobc = 0.0_rprec
                 CALL convert_boundary(rbc,zbs,rhobc,mpol1d,ntord,rho_exp)
                 DO n = LBOUND(lrho_opt,1), UBOUND(lrho_opt,1)
                    DO m = LBOUND(lrho_opt,2), UBOUND(lrho_opt,2)
                       IF (lrho_opt(n,m) .and. (m /= 0 .or. n >= 0)) THEN
                          nvar_in = nvar_in + 1
                          vars(nvar_in) = rhobc(n,m)
                          IF (lauto_domain) THEN
                             bound_min(n,m) = rhobc(n,m) - ABS(pct_domain*rhobc(n,m))
                             bound_max(n,m) = rhobc(n,m) + ABS(pct_domain*rhobc(n,m))
                          END IF
                          vars_min(nvar_in) = bound_min(n,m)
                          vars_max(nvar_in) = bound_max(n,m)
                          var_dex(nvar_in)  = irhobc
                          diag(nvar_in)     = drho_opt(n,m)
                          arr_dex(nvar_in,1) = n
                          arr_dex(nvar_in,2) = m
                       ELSE
                          lrho_opt(n,m) = .FALSE.
                       END IF
                    END DO
                 END DO
              END IF
              IF (ANY(ldeltamn_opt)) THEN
                 deltamn = 0.0_rprec
                 rbc_temp = rbc
                 zbs_temp = zbs
                 CALL convert_boundary_PG(rbc_temp,zbs_temp,deltamn,mpol1d,ntord)
                 DO n = LBOUND(ldeltamn_opt,1), UBOUND(ldeltamn_opt,1)
                    DO m = LBOUND(ldeltamn_opt,2), UBOUND(ldeltamn_opt,2)
                       IF (ldeltamn_opt(n,m) .and. .not.(n == 0 .and. m==0)) THEN
                          nvar_in = nvar_in + 1
                          vars(nvar_in) = deltamn(n,m)
                          IF (lauto_domain) THEN
                             delta_min(n,m) = deltamn(n,m) - ABS(pct_domain*deltamn(n,m))
                             delta_max(n,m) = deltamn(n,m) + ABS(pct_domain*deltamn(n,m))
                          END IF
                          vars_min(nvar_in) = delta_min(n,m)
                          vars_max(nvar_in) = delta_max(n,m)
                          var_dex(nvar_in)  = ideltamn
                          diag(nvar_in)     = ddeltamn_opt(n,m)
                          arr_dex(nvar_in,1) = n
                          arr_dex(nvar_in,2) = m
                       END IF
                    END DO
                 END DO
              END IF
              IF (ANY(lbound_opt)) THEN
                 IF (lbound_opt(0,0)) THEN
                    IF (lauto_domain) THEN
                       rbc_min(0,0) = rbc(0,0) - ABS(pct_domain*rbc(0,0))
                       rbc_max(0,0) = rbc(0,0) + ABS(pct_domain*rbc(0,0))
                    END IF
                    nvar_in = nvar_in + 1
                    vars(nvar_in) = rbc(0,0)
                    vars_min(nvar_in) = rbc_min(0,0)
                    vars_max(nvar_in) = rbc_max(0,0)
                    var_dex(nvar_in) = ibound_rbc
                    diag(nvar_in)    = dbound_opt(0,0)
                    arr_dex(nvar_in,1) = 0
                    arr_dex(nvar_in,2) = 0
                    IF (lasym) THEN
                       IF (lauto_domain) THEN
                          zbc_min(0,0) = zbc(0,0) - ABS(pct_domain*zbc(0,0))
                          zbc_max(0,0) = zbc(0,0) + ABS(pct_domain*zbc(0,0))
                       END IF
                       nvar_in = nvar_in + 1
                       vars(nvar_in) = zbc(0,0)
                       vars_min(nvar_in) = zbc_min(0,0)
                       vars_max(nvar_in) = zbc_max(0,0)
                       var_dex(nvar_in) = ibound_zbc
                       diag(nvar_in)    = dbound_opt(0,0)
                       arr_dex(nvar_in,1) = 0
                       arr_dex(nvar_in,2) = 0
                    END IF
                 END IF
                 DO n = LBOUND(lbound_opt,1), UBOUND(lbound_opt,1)
                    DO m = 0, UBOUND(lbound_opt,2)
                       IF (m==0 .and. n<=0) CYCLE
                       IF (lbound_opt(n,m)) THEN
                          IF (lauto_domain) THEN
                             rbc_min(n,m) = rbc(n,m) - ABS(pct_domain*rbc(n,m))
                             rbc_max(n,m) = rbc(n,m) + ABS(pct_domain*rbc(n,m))
                          END IF
                          nvar_in = nvar_in + 1
                          vars(nvar_in) = rbc(n,m)
                          vars_min(nvar_in) = rbc_min(n,m)
                          vars_max(nvar_in) = rbc_max(n,m)
                          var_dex(nvar_in) = ibound_rbc
                          diag(nvar_in)    = dbound_opt(n,m)
                          arr_dex(nvar_in,1) = n
                          arr_dex(nvar_in,2) = m
                          IF (lauto_domain) THEN
                             zbs_min(n,m) = zbs(n,m) - ABS(pct_domain*zbs(n,m))
                             zbs_max(n,m) = zbs(n,m) + ABS(pct_domain*zbs(n,m))
                          END IF
                          nvar_in = nvar_in + 1
                          vars(nvar_in) = zbs(n,m)
                          vars_min(nvar_in) = zbs_min(n,m)
                          vars_max(nvar_in) = zbs_max(n,m)
                          var_dex(nvar_in) = ibound_zbs
                          diag(nvar_in)    = dbound_opt(n,m)
                          arr_dex(nvar_in,1) = n
                          arr_dex(nvar_in,2) = m
                          IF (lasym) THEN
                             IF (lauto_domain) THEN
                                rbs_min(n,m) = rbs(n,m) - ABS(pct_domain*rbs(n,m))
                                rbs_max(n,m) = rbs(n,m) + ABS(pct_domain*rbs(n,m))
                             END IF
                             nvar_in = nvar_in + 1
                             vars(nvar_in) = rbs(n,m)
                             vars_min(nvar_in) = rbs_min(n,m)
                             vars_max(nvar_in) = rbs_max(n,m)
                             var_dex(nvar_in) = ibound_rbs
                             diag(nvar_in)    = dbound_opt(n,m)
                             arr_dex(nvar_in,1) = n
                             arr_dex(nvar_in,2) = m
                             IF (lauto_domain) THEN
                                zbc_min(n,m) = zbc(n,m) - ABS(pct_domain*zbc(n,m))
                                zbc_max(n,m) = zbc(n,m) + ABS(pct_domain*zbc(n,m))
                             END IF
                             nvar_in = nvar_in + 1
                             vars(nvar_in) = zbc(n,m)
                             vars_min(nvar_in) = zbc_min(n,m)
                             vars_max(nvar_in) = zbc_max(n,m)
                             var_dex(nvar_in) = ibound_zbc
                             diag(nvar_in)    = dbound_opt(n,m)
                             arr_dex(nvar_in,1) = n
                             arr_dex(nvar_in,2) = m
                          END IF
                       END IF
                    END DO
                 END DO
              END IF
              IF (ANY(lcoil_spline)) THEN
                 DO n = LBOUND(lcoil_spline,1), UBOUND(lcoil_spline,1)
                    ! Actual no. of knots for coil spline n, less two to enforce periodicity of f,f'
                    nknots = COUNT(coil_splinesx(n,:) >= 0.0) - 2

                    DO m = 1,nknots-4
                       IF (lcoil_spline(n,m)) THEN
                          IF (lauto_domain) THEN
                             coil_splinefx_min(n,m) = coil_splinefx(n,m) - ABS(pct_domain*coil_splinefx(n,m))
                             coil_splinefx_max(n,m) = coil_splinefx(n,m) + ABS(pct_domain*coil_splinefx(n,m))
                             coil_splinefy_min(n,m) = coil_splinefy(n,m) - ABS(pct_domain*coil_splinefy(n,m))
                             coil_splinefy_max(n,m) = coil_splinefy(n,m) + ABS(pct_domain*coil_splinefy(n,m))
                             coil_splinefz_min(n,m) = coil_splinefz(n,m) - ABS(pct_domain*coil_splinefz(n,m))
                             coil_splinefz_max(n,m) = coil_splinefz(n,m) + ABS(pct_domain*coil_splinefz(n,m))
                          END IF

                          IF ((m > 1).OR.(coil_type(n).NE.'M').OR.(.NOT.lwindsurf)) THEN !u0 fixed for mod on ws
                             nvar_in = nvar_in + 1
                             vars(nvar_in) = coil_splinefx(n,m)
                             vars_min(nvar_in) = coil_splinefx_min(n,m)
                             vars_max(nvar_in) = coil_splinefx_max(n,m)
                             var_dex(nvar_in) = icoil_splinefx
                             diag(nvar_in)    = dcoil_spline(n,m)
                             arr_dex(nvar_in,1) = n
                             arr_dex(nvar_in,2) = m
                          END IF

                          nvar_in = nvar_in + 1
                          vars(nvar_in) = coil_splinefy(n,m)
                          vars_min(nvar_in) = coil_splinefy_min(n,m)
                          vars_max(nvar_in) = coil_splinefy_max(n,m)
                          var_dex(nvar_in) = icoil_splinefy
                          diag(nvar_in)    = dcoil_spline(n,m)
                          arr_dex(nvar_in,1) = n
                          arr_dex(nvar_in,2) = m

                          ! z gets ignored if winding surface is present;
                          !  z0 is held fixed for modular coils.
                          IF ((.NOT.lwindsurf).AND.((m > 1).OR.(coil_type(n).NE.'M'))) THEN
                             nvar_in = nvar_in + 1
                             vars(nvar_in) = coil_splinefz(n,m)
                             vars_min(nvar_in) = coil_splinefz_min(n,m)
                             vars_max(nvar_in) = coil_splinefz_max(n,m)
                             var_dex(nvar_in) = icoil_splinefz
                             diag(nvar_in)    = dcoil_spline(n,m)
                             arr_dex(nvar_in,1) = n
                             arr_dex(nvar_in,2) = m
                          END IF
                       END IF
                    END DO !m
                 END DO !n
              END IF
              ier = -327
              CALL stellopt_prof_to_vmec('init',ier)
         CASE('test')
              IF (lxval_opt) THEN
                 IF (lauto_domain) THEN
                    xval_min = phiedge - ABS(pct_domain*xval)
                    xval_max = phiedge + ABS(pct_domain*xval)
                 END IF
                 nvar_in = nvar_in + 1
                 vars(nvar_in) = xval
                 vars_min(nvar_in) = xval_min
                 vars_max(nvar_in) = xval_max
                 var_dex(nvar_in) = ixval
                 diag(nvar_in)    = dxval_opt
                 arr_dex(nvar_in,1) = 1
              END IF
              IF (lyval_opt) THEN
                 IF (lauto_domain) THEN
                    yval_min = phiedge - ABS(pct_domain*yval)
                    yval_max = phiedge + ABS(pct_domain*yval)
                 END IF
                 nvar_in = nvar_in + 1
                 vars(nvar_in) = yval
                 vars_min(nvar_in) = yval_min
                 vars_max(nvar_in) = yval_max
                 var_dex(nvar_in) = iyval
                 diag(nvar_in)    = dyval_opt
                 arr_dex(nvar_in,1) = 1
              END IF
                 DO m = 1,rosenbrock_dim
                    IF (lRosenbrock_X_opt(m)) THEN
                       IF (lauto_domain) THEN
                         Rosenbrock_X_min(m) = Rosenbrock_X(m) - ABS(pct_domain*Rosenbrock_X(m))
                         Rosenbrock_X_max(m) = Rosenbrock_X(m) + ABS(pct_domain*Rosenbrock_X(m))
                       END IF
                       nvar_in = nvar_in + 1
                       vars(nvar_in) = Rosenbrock_X(m)
                       vars_min(nvar_in) = Rosenbrock_X_min(m)
                       vars_max(nvar_in) = Rosenbrock_X_max(m)
                       var_dex(nvar_in) = iRosenbrock_X
                       diag(nvar_in)    = dRosenbrock_X_opt(m)
                       arr_dex(nvar_in,1) = m
                    END IF
                 END DO
      END SELECT
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER( MPI_COMM_STEL, ierr_mpi )                   ! MPI
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BARRIER_ERR,'stellot_init',ierr_mpi)
!DEC$ ENDIF

      ! Now initalize the targets
      mtargets = 1
      CALL stellopt_load_targets(mtargets,fvec_temp,ier,-1)          ! Count
      ALLOCATE(vals(mtargets),targets(mtargets),sigmas(mtargets),target_dex(mtargets))
      mtargets = 1
      CALL stellopt_load_targets(mtargets,fvec_temp,ier,-2)          ! Load index
      IF (ier /=0) CALL handle_err(NAMELIST_READ_ERR,'IN STELLOPT_LOAD_TARGETS',ier)
      IF (lverb) THEN
         WRITE(6,*) '-----  Optimization  -----'
         WRITE(6,*) '   =======VARS======='
         DO i = 1, nvars
            CALL write_vars(6,var_dex(i),arr_dex(i,1),arr_dex(i,2))
         END DO
         ! Check H-B represenation
         IF (ANY(lrho_opt)) THEN
            rbc_temp = rbc
            zbs_temp = zbs
            CALL unique_boundary(rbc_temp,zbs_temp,rhobc,mpol1d,ntord,mpol-1,ntor,mpol-1,rho_exp)
            delta = 0
            DO m = 0, mpol
               DO n = -ntor, ntor
                  IF (rbc(n,m) /= 0) delta = delta + (rbc(n,m)-rbc_temp(n,m))**2/rbc(n,m)**2
                  IF (zbs(n,m) /= 0) delta = delta + (zbs(n,m)-zbs_temp(n,m))**2/zbs(n,m)**2
               END DO
            END DO
            delta = sqrt(delta)
            WRITE(6,'(A,F7.2,A)')'   == Accuracy of conversion = ',100*(1-delta),'%  =='
            iunit = 12; ier = 0
            CALL safe_open(iunit,ier,'rhomn.txt','unknown','formatted')
            DO m = LBOUND(lrho_opt,2), UBOUND(lrho_opt,2)
               DO n = LBOUND(lrho_opt,1), UBOUND(lrho_opt,1)
                  IF (rhobc(n,m) /= 0.0_rprec) WRITE(iunit,*) n,m,rhobc(n,m)
               END DO
            END DO
            CLOSE(iunit)
         END IF
         ! Check PG represenation
         IF (ANY(ldeltamn_opt)) THEN
            rbc_temp = rbc
            zbs_temp = zbs
            CALL unique_boundary_PG(rbc_temp,zbs_temp,deltamn,ntord,mpol1d,mpol,ntor)
            delta = 0
            DO m = 0, mpol
               DO n = -ntor, ntor
                  IF (rbc(n,m) /= 0) delta = delta + (rbc(n,m)-rbc_temp(n,m))**2/rbc(n,m)**2
                  IF (zbs(n,m) /= 0) delta = delta + (zbs(n,m)-zbs_temp(n,m))**2/zbs(n,m)**2
               END DO
            END DO
            delta = sqrt(delta)
            WRITE(6,'(A,F7.2,A)')'   == Accuracy of conversion = ',100*(1-delta),'%  =='
            iunit = 12; ier = 0
            CALL safe_open(iunit,ier,'deltamn.txt','unknown','formatted')
            DO m = LBOUND(ldeltamn_opt,2), UBOUND(ldeltamn_opt,2)
               DO n = LBOUND(ldeltamn_opt,1), UBOUND(ldeltamn_opt,1)
                  IF (deltamn(n,m) /= 0.0_rprec) WRITE(iunit,*) n,m,deltamn(n,m)
                  !WRITE(iunit,*) n,m,deltamn(n,m)
               END DO
            END DO
            CLOSE(iunit)
         END IF
         WRITE(6,*) '   ======TARGETS====='
         m = 0
         DO i = 1, mtargets
            IF (target_dex(i) == m) CYCLE
            CALL write_targets(6,target_dex(i))
            m=target_dex(i)
         END DO
         WRITE(6,*) '   =================='
         WRITE(6,*) '   Number of Parameters: ',nvars
         WRITE(6,*) '   Number of Targets:    ',mtargets
         IF (lno_restart) WRITE(6,*) '   !!!! EQUILIBRIUM RESTARTING NOT UTILIZED !!!!'
      END IF

      IF (myid == master) THEN
         ! Write the var_labels file
         iunit = 12; ier = 0
         CALL safe_open(iunit,ier,'var_labels','unknown','formatted')
         WRITE(iunit,'(I4.4)') nvars
         DO i = 1, nvars
            CALL write_vars(iunit,var_dex(i),arr_dex(i,1),arr_dex(i,2))
         END DO
         WRITE(iunit,'(I8.8)') mtargets
         DO i = 1, mtargets
            CALL write_targets(iunit,target_dex(i))
         END DO
         CLOSE(iunit)
      END IF
      
      
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER( MPI_COMM_STEL, ierr_mpi )                   ! MPI
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BARRIER_ERR,'stellot_init',ierr_mpi)
!DEC$ ENDIF
      
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_init
