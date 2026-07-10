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
#if defined(LHDF5)
      USE ez_hdf5
#endif
      USE EZspline
      USE EZspline_obj
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        ltst        logical for supressing screen output
!        tstr1/2     String for calling paraexe
!-----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL        :: ltst
      INTEGER        :: ier, i, iunit, ntimesteps_ecrh, nbeams
      CHARACTER(256) :: tstr1,tstr2
      REAL(rprec)    :: dt, tend_restart
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: total_current_edge_bc,time_grid_edge_bc
      INTEGER        :: total_current_edge_bc_ntimesteps
      TYPE(EZspline1_r8) :: temp_spl
      INTEGER :: bcs0(2)
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

      ! ECCD power (saved when using TRAVIS)
      CALL mpialloc(THRIFT_DPECRHDV,  nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_dpecrhdv)
      CALL mpialloc(THRIFT_PECRH,  nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_pecrh)
      
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
      IF (lrestart_from_file) CALL thrift_restart

      CALL MPI_BCAST(dt,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_BCAST(dt_first_iter,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
      IF (lrestart_from_file .AND. ntimesteps == 1) dt = dt_first_iter
      IF (lrestart_from_file) tend_restart = tstart - dt_first_iter

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

      ! Check that the number of ion species in the restart file matches the profiles file
      IF (lrestart_from_file .AND. nion_prof_restart /= nion_prof) THEN
         WRITE(6,*) '!!!!!!!!!!!!ERRROR!!!!!!!!!!!!!!'
         WRITE(6,*) '  nion_prof mismatch: restart file has ', nion_prof_restart, &
                    ' ion species but profiles file has ', nion_prof
         WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         STOP
      END IF

      ! Allocate particle and heat fluxes (do it here because nion_prof only now available)
      CALL mpialloc(THRIFT_GNEO,   nion_prof+1, nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_gneo) 
      CALL mpialloc(THRIFT_QNEO,   nion_prof+1, nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_qneo)
      ! Allocate densities, temperatures and pressures
      CALL mpialloc(THRIFT_DENS,   nion_prof+1, nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_dens)
      CALL mpialloc(THRIFT_TEMP,   nion_prof+1, nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_temp)
      CALL mpialloc(THRIFT_PRESS,  nion_prof+1, nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_press)    
      CALL mpialloc(THRIFT_FAST_ALPHAS_DENS,    nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_fast_alphas_dens)

      IF (myid_sharmem == master) THEN
        FORALL(i = 1:nrho) THRIFT_RHO(i) = DBLE(i-0.5)/DBLE(nrho) ! (half) rho grid
        FORALL(i = 1:nsj)  THRIFT_S(i)   = DBLE(i-1)/DBLE(nsj-1)  ! (full)  s  grid
        FORALL(i = 1:ntimesteps) THRIFT_T(i) = tstart + (i-1)*dt  !       time grid
      END IF

      ! Read PECRH_AUX_T and PECRH_AUX_F in case they exist in profiles file
      IF(leccd) THEN
         IF (myid_sharmem == master) THEN
            !
            nbeams = 0
            DO i = 1,nsys
               IF (ANY(antennaposition_ecrh(i,:) .ne. 0)) nbeams = nbeams + 1
            END DO
            !
            CALL open_hdf5(TRIM(prof_string),fid,ier,LCREATE=.false.)
            IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,TRIM(prof_string),ier)
            !
            IF( TRIM(power_type) .EQ. 'read_from_file' ) THEN
               CALL read_scalar_hdf5(fid,'ecrh_ntimesteps',ier,INTVAR=ntimesteps_ecrh)
               IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'ecrh_ntimesteps',ier)
               CALL read_scalar_hdf5(fid,'ecrh_ngyrotrons',ier,INTVAR=ngyrotrons)
               IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'ecrh_ngyrotrons',ier)
               !
               ! Let's check that nygrotrons is the same as number of antennas in input file
               IF(nbeams .ne. ngyrotrons) THEN
                  WRITE(6,*) '!!!!!!!!!!!!ERRROR!!!!!!!!!!!!!!'
                  WRITE(6,*) '  Number of gyrotrons in profiles file different from number of beams in input file '
                  WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                  STOP
               END IF
               !
               ALLOCATE(PECRH_AUX_T(ngyrotrons,ntimesteps_ecrh),PECRH_AUX_F(ngyrotrons,ntimesteps_ecrh))
               !
               CALL read_var_hdf5(fid,'PECRH_AUX_T',ngyrotrons,ntimesteps_ecrh,ier,DBLVAR=PECRH_AUX_T)
               IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'PECRH_AUX_T',ier)
               CALL read_var_hdf5(fid,'PECRH_AUX_F',ngyrotrons,ntimesteps_ecrh,ier,DBLVAR=PECRH_AUX_F)
               IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'PECRH_AUX_F',ier)
               !
               CALL close_hdf5(fid,ier)   
            ELSE IF( TRIM(power_type) .EQ. 'read_from_namelist') THEN
               ! When power read from namelist, assumes the same power at all times
               ntimesteps_ecrh = ntimesteps
               ngyrotrons = nbeams
               ALLOCATE(PECRH_AUX_T(ngyrotrons,ntimesteps_ecrh),PECRH_AUX_F(ngyrotrons,ntimesteps_ecrh))
               DO i=1,ngyrotrons
                  PECRH_AUX_T(i,:) = THRIFT_T
                  PECRH_AUX_F(i,:) = power_ecrh(i)
               END DO
            ELSE
               WRITE(6,*) '  power_type MUST BE read_from_file OR read_from_namelist '
               FLUSH(6)
               STOP
            END IF
         END IF
         ! CALL barrier??
         CALL MPI_BCAST(ntimesteps_ecrh,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(ngyrotrons,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
         IF( .NOT. ALLOCATED(PECRH_AUX_T)) ALLOCATE(PECRH_AUX_T(ngyrotrons,ntimesteps_ecrh))
         IF( .NOT. ALLOCATED(PECRH_AUX_F)) ALLOCATE(PECRH_AUX_F(ngyrotrons,ntimesteps_ecrh))
         CALL MPI_BCAST(PECRH_AUX_T,ngyrotrons*ntimesteps_ecrh,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(PECRH_AUX_F,ngyrotrons*ntimesteps_ecrh,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
      END IF

      ! If dirichlet boundary condition, need to read experimental current as a function of time from file
      IF(trim(edge_bc_type) == 'dirichlet') THEN
         WRITE(6,*) 'EDGE BC: dirichlet (read from file) '
         !
         CALL open_hdf5(TRIM(prof_string),fid,ier,LCREATE=.false.)
         IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,TRIM(prof_string),ier)
         !
         CALL read_scalar_hdf5(fid,'total_current_edge_bc_ntimesteps',ier,INTVAR=total_current_edge_bc_ntimesteps)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'total_current_edge_bc_ntimesteps',ier)
         !
         ALLOCATE(total_current_edge_bc(total_current_edge_bc_ntimesteps),time_grid_edge_bc(total_current_edge_bc_ntimesteps))
         !
         CALL read_var_hdf5(fid,'total_current_edge_bc',total_current_edge_bc_ntimesteps,ier,DBLVAR=total_current_edge_bc)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'total_current_edge_bc',ier)
         !
         CALL read_var_hdf5(fid,'time_grid_edge_bc',total_current_edge_bc_ntimesteps,ier,DBLVAR=time_grid_edge_bc)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'time_grid_edge_bc',ier)
         !
         CALL close_hdf5(fid,ier)
         
         ! Now create spline and evaluate THRIFT_DIRICHLET_EDGE_BC at THRIFT_T
         CALL mpialloc(THRIFT_DIRICHLET_EDGE_BC, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_dirichlet_edge_bc)
         bcs0=(/ 0, 0/)
         CALL EZspline_init(temp_spl,total_current_edge_bc_ntimesteps,bcs0,ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'init: edge BC spline',ier)
         temp_spl%x1          = time_grid_edge_bc
         temp_spl%isHermite   = 1
         CALL EZspline_setup(temp_spl,total_current_edge_bc,ier,EXACT_DIM=.true.)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'setup: edge BC spline',ier)
         CALL EZspline_interp(temp_spl,ntimesteps,THRIFT_T,THRIFT_DIRICHLET_EDGE_BC,ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'interp: edge BC spline',ier)
         CALL EZspline_free(temp_spl,ier)

         DEALLOCATE(total_current_edge_bc,time_grid_edge_bc)

      ELSEIF(trim(edge_bc_type) == 'robin') THEN
         WRITE(6,*) 'EDGE BC: robin'

      ELSE
         WRITE(6,*) '!!!!!!!!!!!! ERROR !!!!!!!!!!!!!!'
         WRITE(6,*) '  edge_bc_type must be either dirichlet or robin'
         WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         STOP
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

      ! Run VMEC+booz_xform+dkes using restart profiles so DKES_D** are
      ! available at the first thrift_penta call when add_NEO=.true.
      IF (lrestart_from_file .AND. solve_plasma_equations) CALL thrift_restart_equil

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_init

