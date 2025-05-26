!-----------------------------------------------------------------------
!     Subroutine:    thrift_init
!     Authors:       L. van Ham
!     Date:          11/XX/2022
!     Description:   This subroutine initialzies the code for performing
!                    a run.
!-----------------------------------------------------------------------
      SUBROUTINE thrift_init
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE thrift_runtime
      USE thrift_input_mod
      USE thrift_vars
      USE thrift_profiles_mod
      USE diagno_input_mod, ONLY:   read_diagno_input
      USE penta_interface_mod, ONLY:   init_penta_input, &
                                       read_penta_run_params_namelist
      USE thrift_plasma_solver_mod, ONLY: initialize_plasma_solver, Nt_total_plasma_solver, &
      dt_plasma_solver, time_plasma_grid, N_plasma_steps_per_THRIFT_step
      USE thrift_equil, ONLY : eq_Aminor, eq_phiedge, vp_spl, bsq_spl, bcs1
      USE safe_open_mod
      USE mpi_params
      USE mpi_inc
      USE mpi_sharmem
      USE EZspline
      USE EZspline_obj
#if defined(LHDF5)
      USE ez_hdf5
#endif
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        ltst        logical for supressing screen output
!        tstr1/2     String for calling paraexe
!-----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL        :: ltst
      INTEGER        :: ier, i, iunit, ntimesteps_restart, ns_restart, k
      CHARACTER(256) :: tstr1,tstr2
      REAL(rprec)    :: dt, tend_restart
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: temp2d
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: temp1d
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: temp3d
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      ! Read the Input Namelist for THRIFT
      IF (lverb) WRITE(6,'(A)') '----- THRIFT Input Parameters -----'
      CALL init_thrift_input

      ! Read the THRIFT input
      IF (lvmec) THEN
         CALL read_thrift_input('input.' // TRIM(id_string),ier)
      END IF

      ! Read diagno file
      IF (ldiagno) THEN
         CALL read_diagno_input('input.' // TRIM(id_string),ier)
      END IF

      ! Output to screen
      IF (lverb) THEN 
         WRITE(6,'(A)') '   FILE:             input.' // TRIM(id_string)
         WRITE(6,'(A)') '   BOOTSTRAP MODEL:  ' // TRIM(bootstrap_type)
         WRITE(6,'(A)') '   ETAPAR MODEL:  ' // TRIM(etapar_type)
         IF (leccd) WRITE(6,'(A)') '   ECCD MODEL:       ' // TRIM(eccd_type)
         IF (lnbcd) WRITE(6,'(A)') '   NBCD MODEL:       ' // TRIM(nbcd_type)
         WRITE(6,'(A)') ''
         WRITE(6,'(A11,I5)') '   NRHO:   ', nrho
         WRITE(6,'(A11,I5)') '   NS:     ', nsj
         WRITE(6,'(A11,I5)') '   NT:     ', ntimesteps
         WRITE(6,'(A11,F8.4)') '   TSTART: ', tstart
         WRITE(6,'(A11,F8.4)') '   TEND:   ', tend
         IF (lvmec_reset) WRITE(6,'(A)') '   VMEC RESET FEATURE ON.'
      END IF

      ! Grid allocations
      CALL mpialloc(THRIFT_RHO,       nrho,    myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_rho)
      CALL mpialloc(THRIFT_RHOFULL, nrho+2,    myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_rhofull)
      CALL mpialloc(THRIFT_S,          nsj,    myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_s)
      CALL mpialloc(THRIFT_SNOB,     nsj-2,    myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_snob)
      CALL mpialloc(THRIFT_T,   ntimesteps,    myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_t)
      CALL mpialloc(THRIFT_PHIEDGE, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_phiedge)

      ! Current densities
      CALL mpialloc(THRIFT_J,         nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_j)
      CALL mpialloc(THRIFT_JBOOT,     nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_jboot)
      CALL mpialloc(THRIFT_JPLASMA,   nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_jplasma)
      CALL mpialloc(THRIFT_JECCD,     nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_jeccd)
      CALL mpialloc(THRIFT_JNBCD,     nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_jnbcd)
      CALL mpialloc(THRIFT_JOHMIC,    nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_johmic)
      CALL mpialloc(THRIFT_JSOURCE,   nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_jsource)     
      ! Total currents
      CALL mpialloc(THRIFT_I,       nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_i)
      CALL mpialloc(THRIFT_IBOOT,   nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_iboot)
      CALL mpialloc(THRIFT_IPLASMA, nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_iplasma)
      CALL mpialloc(THRIFT_IECCD,   nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_ieccd)
      CALL mpialloc(THRIFT_INBCD,   nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_inbcd)
      CALL mpialloc(THRIFT_IOHMIC,  nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_iohmic)
      CALL mpialloc(THRIFT_ISOURCE, nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_isource)
      CALL mpialloc(THRIFT_UGRID,   nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_ugrid)  
      
      ! Profile variables
      CALL mpialloc(THRIFT_ETAPARA,nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_etapara)  
      CALL mpialloc(THRIFT_PPRIME ,nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_pprime)  
      CALL mpialloc(THRIFT_P      ,nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_p)  
    
      ! Magnetic variables
      CALL mpialloc(THRIFT_S11,    nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_s11)     
      CALL mpialloc(THRIFT_S12,    nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_s12)   
      CALL mpialloc(THRIFT_BAV,    nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_bav)     
      CALL mpialloc(THRIFT_BSQAV,  nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_bsqav)    
      CALL mpialloc(THRIFT_IOTA,   nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_iota)    
      CALL mpialloc(THRIFT_AMINOR, nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_aminor)     
      CALL mpialloc(THRIFT_RMAJOR, nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_rmajor)  
      CALL mpialloc(THRIFT_VP,     nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_vp)      
      CALL mpialloc(THRIFT_BVAV,   nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_bvav)   
      CALL mpialloc(THRIFT_BETATOT,   ntimesteps,   myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_betatot)

      ! Electric field
      CALL mpialloc(THRIFT_EPARB,  nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_eparb)
      CALL mpialloc(THRIFT_ER,     nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_er) 

      ! ABCD
      CALL mpialloc(THRIFT_COEFF_A, nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_coeff_a)
      CALL mpialloc(THRIFT_COEFF_B, nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_coeff_b)
      CALL mpialloc(THRIFT_COEFF_C, nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_coeff_c)
      CALL mpialloc(THRIFT_COEFF_D, nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_coeff_d)
      CALL mpialloc(THRIFT_COEFF_BP,nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_coeff_bp)
      CALL mpialloc(THRIFT_COEFF_CP,nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_coeff_cp)
      CALL mpialloc(THRIFT_COEFF_DP,nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_coeff_dp)
      
      ! Alphas
      CALL mpialloc(THRIFT_ALPHA1,   nsj-2, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_alpha1)
      CALL mpialloc(THRIFT_ALPHA2,   nsj-2, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_alpha2)
      CALL mpialloc(THRIFT_ALPHA3,   nsj-2, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_alpha3)
      CALL mpialloc(THRIFT_ALPHA4,   nsj-2, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_alpha4)
      
      ! System of equations
      CALL mpialloc(THRIFT_MATLD,    nsj-1, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_matld)
      CALL mpialloc(THRIFT_MATMD,      nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_matmd)
      CALL mpialloc(THRIFT_MATUD,    nsj-1, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_matud)
      CALL mpialloc(THRIFT_MATRHS,     nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_matrhs)   
      
      ! Restart arrays
      CALL mpialloc(UGRID_RESTART,   nsj, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_ugrid_restart)
      CALL mpialloc(J_RESTART,       nsj, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_j_restart)

      ! Read the Bootstrap input
      CALL tolower(bootstrap_type)
      SELECT CASE (TRIM(bootstrap_type))
         CASE('bootsj')
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
         CASE('dkespenta')
            ier = 0
            CALL init_penta_input
            CALL read_penta_run_params_namelist('input.'//TRIM(id_string),ier)
            IF (ier < 0 .and. myid == master) THEN
               WRITE(6,*) '!!!!!!!!!!!!ERRROR!!!!!!!!!!!!!!'
               WRITE(6,*) '  RUN_PARAMS Namelist not found     '
               WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               STOP
            END IF
         CASE('sfincs')
      END SELECT

      ! Check that tend > tstart
      IF(tend < tstart .and. lverb) THEN 
         WRITE(6,*) '!!!!!!!!!!!!ERRROR!!!!!!!!!!!!!!'
         WRITE(6,*) '          tend < tstart         '
         WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         STOP
      ENDIF
      
      ! Define grids
      IF( ntimesteps==1 ) THEN 
         dt = 0.0_rprec
      ELSE IF( ntimesteps > 1) THEN 
         dt = (tend-tstart)/(ntimesteps-1)
      ELSE
         IF(lverb) THEN
            WRITE(6,*) '!!!!!!!!!!!!ERRROR!!!!!!!!!!!!!!'
            WRITE(6,*) '          ntimesteps < 1        '
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            STOP
         END IF
      END IF

      ! Read restart file
      IF (lrestart_from_file) THEN
         UGRID_RESTART = 0.0
         IF (lverb) THEN 
            WRITE(6,'(A)') '----- Reading Restart File -----'
            WRITE(6,'(A)')  '   FILE: '//TRIM(restart_filename)
         END IF
         ! Read file 
         IF (myid_sharmem == master) THEN
            CALL open_hdf5(TRIM(restart_filename),fid,ier,LCREATE=.false.)
            IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,TRIM(restart_filename),ier)

            CALL read_scalar_hdf5(fid,'ntimesteps',ier,INTVAR=ntimesteps_restart)
            IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'ntimesteps',ier)

            CALL read_scalar_hdf5(fid,'nssize',ier,INTVAR=ns_restart)
            IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'nssize',ier)

            !Check that ns_restart is equal to current ns
            IF(ns_restart /= nsj) THEN
               WRITE(6,*) '!!!!!!!!!!!!ERRROR!!!!!!!!!!!!!!'
               WRITE(6,*) '  ns_restart different from nsj '
               WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               STOP
            ENDIF

            ALLOCATE(temp2d(ns_restart,ntimesteps_restart),temp1d(ntimesteps_restart))

            CALL read_var_hdf5(fid,'THRIFT_T',ntimesteps_restart,ier,DBLVAR=temp1d)
            IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'THRIFT_T',ier)
            tend_restart = temp1d(ntimesteps_restart)

            IF(lverb) WRITE(6,'(A17,F8.4)') '   TEND_RESTART: ', tend_restart

            ! Check tstart > tend_restart
            IF(tstart < tend_restart) THEN 
               WRITE(6,*) '!!!!!!!!!!!!ERRROR!!!!!!!!!!!!!!'
               WRITE(6,*) '          tstart < tend_resart  '
               WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               STOP
            ENDIF

            dt_first_iter = tstart - tend_restart
            IF(ntimesteps == 1) dt = dt_first_iter

            CALL read_var_hdf5(fid,'THRIFT_UGRID',ns_restart,ntimesteps_restart,ier,DBLVAR=temp2d)
            IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'THRIFT_UGRID',ier)
            UGRID_RESTART = temp2d(:,ntimesteps_restart)

            CALL read_var_hdf5(fid,'THRIFT_J',ns_restart,ntimesteps_restart,ier,DBLVAR=temp2d)
            IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'THRIFT_J',ier)
            J_RESTART = temp2d(:,ntimesteps_restart)

            CALL read_var_hdf5(fid,'eq_Aminor',ier,DBLVAR=eq_Aminor)
            IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'eq_Aminor',ier)

            CALL read_var_hdf5(fid,'THRIFT_PHIEDGE',ntimesteps_restart,ier,DBLVAR=temp1d)
            IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'THRIFT_PHIEDGE',ier)
            eq_phiedge = temp1d(ntimesteps_restart)

            CALL read_var_hdf5(fid,'THRIFT_VP',ns_restart,ntimesteps_restart,ier,DBLVAR=temp2d)
            IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'THRIFT_VP',ier)

            ! dV/dPhi Spline (Volume derivative)
            bcs1=(/ 0, 0/)
            IF (EZspline_allocated(vp_spl)) CALL EZspline_free(vp_spl,ier)
            CALL EZspline_init(vp_spl,ns_restart,bcs1,ier)
            IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'thrift_init: vp_spl',ier)
            vp_spl%isHermite = 0
            FORALL (k=1:ns_restart) vp_spl%x1(k) = sqrt(DBLE(k-1)/DBLE(ns_restart-1))
            CALL EZspline_setup(vp_spl,temp2d(:,ntimesteps_restart)/eq_phiedge,ier,EXACT_DIM=.true.)
            IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'thrift_init: vp_spl',ier)

            CALL read_var_hdf5(fid,'THRIFT_BSQAV',ns_restart,ntimesteps_restart,ier,DBLVAR=temp2d)
            IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'THRIFT_BSQAV',ier)

            ! Bsq Spline
            bcs1=(/ 0, 0/)
            IF (EZspline_allocated(bsq_spl)) CALL EZspline_free(bsq_spl,ier)
            CALL EZspline_init(bsq_spl,ns_restart,bcs1,ier)
            IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'thrift_init: bsq_spl',ier)
            bsq_spl%isHermite = 0
            FORALL (k=1:ns_restart) bsq_spl%x1(k) = sqrt(DBLE(k-1)/DBLE(ns_restart-1))
            CALL EZspline_setup(bsq_spl,temp2d(:,ntimesteps_restart),ier,EXACT_DIM=.true.)
            IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'thrift_init: bsq_spl',ier)

            DEALLOCATE(temp2d,temp1d)
            
            !Close the HDF5 file
            CALL close_hdf5(fid,ier)
            IF (ier /= 0) CALL handle_err(HDF5_CLOSE_ERR,TRIM(restart_filename),ier)
         END IF
      END IF

      CALL MPI_BCAST(dt,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_BCAST(dt_first_iter,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)

      IF(solve_plasma_equations) THEN
         ! Check dt_plasma_solver and ajust it
         IF( dt_plasma_solver .GT. dt .AND. ntimesteps > 1) THEN
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!ERROR!!!!!!!!!!!!!!'
            WRITE(6,*) '   dt_plasma_solver < dt_THRIFT        '
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            STOP
         END IF

         IF( lrestart_from_file ) THEN
            ! check dt_first_iter is equal to dt
            IF( ABS(dt_first_iter-dt) > 1.0E-12_rprec) THEN
               WRITE(6,*) '!!!!!!!!!!!!!!!!!!! ERROR!!!!!!!!!!!!!!!!!!!!!!!'
               WRITE(6,*) '              NOT POSSIBLE TO HAVE:             '
               WRITE(6,*) '   dt_first_iter != dt with plasma_solver ON    '
               WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               STOP
            END IF
            ! Need to add the plasma steps before tstart of THRIFT
            Nt_total_plasma_solver = 1 + ntimesteps*NINT(dt/dt_plasma_solver)
            dt_plasma_solver = dt / NINT(dt/dt_plasma_solver)
            IF(lverb) PRINT *, 'dt_plasma_solver adjusted to ', dt_plasma_solver
            N_plasma_steps_per_THRIFT_step = NINT(dt/dt_plasma_solver)
            ALLOCATE(time_plasma_grid(Nt_total_plasma_solver))
            FORALL(i = 1:Nt_total_plasma_solver) time_plasma_grid(i) = tend_restart + (i-1)*dt_plasma_solver
         ELSE
            Nt_total_plasma_solver = 1 + (ntimesteps-1)*NINT(dt/dt_plasma_solver)
            dt_plasma_solver = dt / NINT(dt/dt_plasma_solver)
            IF(lverb) PRINT *, 'dt_plasma_solver adjusted to ', dt_plasma_solver
            N_plasma_steps_per_THRIFT_step = NINT(dt/dt_plasma_solver)
            ALLOCATE(time_plasma_grid(Nt_total_plasma_solver))
            FORALL(i = 1:Nt_total_plasma_solver) time_plasma_grid(i) = tstart + (i-1)*dt_plasma_solver
         END IF
      END IF

      ! Now setup the profiles (plasma profiles if not solving plasma eqs; external source profiles if solving plasma eqs.)
      IF(solve_plasma_equations) THEN
         CALL initialize_plasma_solver((TRIM(prof_string))) 
      ELSE
         CALL read_thrift_profh5(TRIM(prof_string))
      ENDIF

      ! Allocate particle and heat fluxes (do it here because nion_prof only now available)
      CALL mpialloc(THRIFT_GNEO,   nion_prof+1, nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_gneo) 
      CALL mpialloc(THRIFT_QNEO,   nion_prof+1, nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_qneo)
      ! Allocate densities, temperatures and pressures
      CALL mpialloc(THRIFT_DENS,   nion_prof+1, nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_dens)
      CALL mpialloc(THRIFT_TEMP,   nion_prof+1, nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_temp)
      CALL mpialloc(THRIFT_PRESS,  nion_prof+1, nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_press)    
      ! Restart vars
      IF(lrestart_from_file) THEN
         CALL mpialloc(DENS_RESTART,   nion_prof+1, ns_restart, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_dens_restart)
         CALL mpialloc(TEMP_RESTART,   nion_prof+1, ns_restart, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_temp_restart)
      END IF

      IF(lrestart_from_file .AND. solve_plasma_equations .AND. myid_sharmem == master) THEN
         CALL open_hdf5(TRIM(restart_filename),fid,ier,LCREATE=.false.)
         IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,TRIM(restart_filename),ier)

         ALLOCATE(temp3d(nion_prof+1,ns_restart,ntimesteps_restart))
         ! Read density of all species at last time step
         CALL read_var_hdf5(fid,'THRIFT_DENS',nion_prof+1,ns_restart,ntimesteps_restart,ier,DBLVAR=temp3d)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'THRIFT_DENS',ier)
         DENS_RESTART = temp3d(:,:,ntimesteps_restart)
         ! Read temperature of all species at last time step
         CALL read_var_hdf5(fid,'THRIFT_TEMP',nion_prof+1,ns_restart,ntimesteps_restart,ier,DBLVAR=temp3d)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'THRIFT_TEMP',ier)
         TEMP_RESTART = temp3d(:,:,ntimesteps_restart)
         DEALLOCATE(temp3d)

         !Close the HDF5 file
         CALL close_hdf5(fid,ier)
         IF (ier /= 0) CALL handle_err(HDF5_CLOSE_ERR,TRIM(restart_filename),ier)
      END IF

      IF (myid_sharmem == master) THEN
        FORALL(i = 1:nrho) THRIFT_RHO(i) = DBLE(i-0.5)/DBLE(nrho) ! (half) rho grid
        FORALL(i = 1:nsj)  THRIFT_S(i)   = DBLE(i-1)/DBLE(nsj-1)  ! (full)  s  grid
        FORALL(i = 1:ntimesteps) THRIFT_T(i) = tstart + (i-1)*dt  !       time grid
      END IF

      ! Extra variables (used in debugging process)
      THRIFT_RHOFULL(1) = 0.0
      THRIFT_RHOFULL(2:nrho+1) = THRIFT_RHO
      THRIFT_RHOFULL(nrho+2) = 1.0
      THRIFT_SNOB = THRIFT_S(2:nsj-1)

      ! Split off workers
      CALL thrift_init_mpisubgroup
      IF (myworkid .ne. master) THEN
         ltst  = .false.
         tstr1 = ''
         tstr2 = ''
         ier_paraexe = 0
         CALL thrift_paraexe(tstr1,tstr2,ltst)
         RETURN
      END IF
      ! - From this point on only the main thread of each run executes

      ! Initialize the equilbrium code
      IF (lvmec) THEN
         ltst = .false.
         tstr1 = 'parvmec_init'
         tstr2 = id_string
         CALL thrift_paraexe(tstr1,tstr2,ltst)
      END IF

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_init

