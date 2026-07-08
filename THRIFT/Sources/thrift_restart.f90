!-----------------------------------------------------------------------
!     Subroutine:    thrift_restart
!     Authors:       A. Coelho
!     Date:          07/2026
!     Description:   Reads the restart file and initialises all restart
!                    arrays.  Must be called after UGRID_RESTART and
!                    J_RESTART have been mpialloc'd in thrift_init.
!-----------------------------------------------------------------------
      SUBROUTINE thrift_restart
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE thrift_runtime
      USE thrift_vars
      USE thrift_equil, ONLY: eq_Aminor, eq_phiedge
      USE thrift_globals, ONLY: nsj, tstart, solve_plasma_equations
      USE mpi_params
      USE mpi_inc
      USE mpi_sharmem
#if defined(LHDF5)
      USE ez_hdf5
#endif
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER     :: ier, ns_restart, ntimesteps_restart
      REAL(rprec) :: tend_restart
      REAL(rprec), ALLOCATABLE :: temp1d(:), temp2d(:,:), temp3d(:,:,:)
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      UGRID_RESTART = 0.0_rprec

      IF (lverb) THEN
         WRITE(6,'(A)') '----- Reading Restart File -----'
         WRITE(6,'(A)') '   FILE: '//TRIM(restart_filename)
      END IF

      ! Read dimension scalars on master to know sizes before mpialloc
      IF (myid_sharmem == master) THEN
         CALL open_hdf5(TRIM(restart_filename),fid,ier,LCREATE=.false.)
         IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,TRIM(restart_filename),ier)

         CALL read_scalar_hdf5(fid,'ntimesteps',ier,INTVAR=ntimesteps_restart)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'ntimesteps',ier)

         CALL read_scalar_hdf5(fid,'nssize',ier,INTVAR=ns_restart)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'nssize',ier)

         IF (ns_restart /= nsj) THEN
            WRITE(6,*) '!!!!!!!!!!!!ERRROR!!!!!!!!!!!!!!'
            WRITE(6,*) '  ns_restart different from nsj '
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            STOP
         END IF

         CALL read_scalar_hdf5(fid,'nion_prof',ier,INTVAR=nion_prof_restart)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'nion_prof',ier)
      END IF

      ! Broadcast nion_prof_restart so all ranks can participate in mpialloc
      CALL MPI_BCAST(nion_prof_restart, 1, MPI_INTEGER, master, MPI_COMM_MYWORLD, ierr_mpi)

      ! Allocate restart density/temperature/Er arrays (all ranks)
      CALL mpialloc(DENS_RESTART,           nion_prof_restart+1, nsj, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_dens_restart)
      CALL mpialloc(TEMP_RESTART,           nion_prof_restart+1, nsj, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_temp_restart)
      CALL mpialloc(DENS_FAST_ALPHAS_RESTART,                nsj, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_dens_fast_alphas_restart)
      CALL mpialloc(ER_RESTART,                              nsj, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_er_restart)

      ! Read all restart arrays on master (file is still open)
      IF (myid_sharmem == master) THEN
         ALLOCATE(temp1d(ntimesteps_restart), temp2d(nsj, ntimesteps_restart))

         CALL read_var_hdf5(fid,'THRIFT_T',ntimesteps_restart,ier,DBLVAR=temp1d)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'THRIFT_T',ier)
         tend_restart = temp1d(ntimesteps_restart)
         IF (lverb) WRITE(6,'(A17,F8.4)') '   TEND_RESTART: ', tend_restart

         IF (tstart < tend_restart) THEN
            WRITE(6,*) '!!!!!!!!!!!!ERRROR!!!!!!!!!!!!!!'
            WRITE(6,*) '          tstart < tend_resart  '
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            STOP
         END IF

         dt_first_iter = tstart - tend_restart

         CALL read_var_hdf5(fid,'THRIFT_UGRID',nsj,ntimesteps_restart,ier,DBLVAR=temp2d)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'THRIFT_UGRID',ier)
         UGRID_RESTART = temp2d(:,ntimesteps_restart)

         CALL read_var_hdf5(fid,'THRIFT_J',nsj,ntimesteps_restart,ier,DBLVAR=temp2d)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'THRIFT_J',ier)
         J_RESTART = temp2d(:,ntimesteps_restart)

         CALL read_var_hdf5(fid,'eq_Aminor',ier,DBLVAR=eq_Aminor)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'eq_Aminor',ier)

         CALL read_var_hdf5(fid,'THRIFT_PHIEDGE',ntimesteps_restart,ier,DBLVAR=temp1d)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'THRIFT_PHIEDGE',ier)
         eq_phiedge = temp1d(ntimesteps_restart)

         DEALLOCATE(temp1d, temp2d)

         IF (solve_plasma_equations) THEN
            ALLOCATE(temp3d(nion_prof_restart+1,nsj,ntimesteps_restart), &
                     temp2d(nsj,ntimesteps_restart))

            CALL read_var_hdf5(fid,'THRIFT_DENS',nion_prof_restart+1,nsj,ntimesteps_restart,ier,DBLVAR=temp3d)
            IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'THRIFT_DENS',ier)
            DENS_RESTART = temp3d(:,:,ntimesteps_restart)

            CALL read_var_hdf5(fid,'THRIFT_FAST_ALPHAS_DENS',nsj,ntimesteps_restart,ier,DBLVAR=temp2d)
            IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'THRIFT_FAST_ALPHAS_DENS',ier)
            DENS_FAST_ALPHAS_RESTART = temp2d(:,ntimesteps_restart)

            CALL read_var_hdf5(fid,'THRIFT_TEMP',nion_prof_restart+1,nsj,ntimesteps_restart,ier,DBLVAR=temp3d)
            IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'THRIFT_TEMP',ier)
            TEMP_RESTART = temp3d(:,:,ntimesteps_restart)

            CALL read_var_hdf5(fid,'THRIFT_ER',nsj,ntimesteps_restart,ier,DBLVAR=temp2d)
            IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'THRIFT_ER',ier)
            ER_RESTART = temp2d(:,ntimesteps_restart)

            DEALLOCATE(temp3d, temp2d)
         END IF

         CALL close_hdf5(fid,ier)
         IF (ier /= 0) CALL handle_err(HDF5_CLOSE_ERR,TRIM(restart_filename),ier)
      END IF

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_restart

!-----------------------------------------------------------------------
!     Subroutine:    thrift_restart_equil
!     Authors:       A. Coelho
!     Date:          07/2026
!     Description:   Reconstructs the VMEC equilibrium from restart
!                    arrays and (if solve_plasma_equations .AND. add_NEO)
!                    runs booz_xform+dkes so that DKES_D** are available
!                    for thrift_penta at the first plasma iteration.
!                    Must be called after thrift_init_mpisubgroup
!                    (requires thrift_paraexe).
!-----------------------------------------------------------------------
      SUBROUTINE thrift_restart_equil
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE thrift_runtime
      USE thrift_vars
      USE thrift_equil, ONLY: eq_phiedge, ns_eq
      USE thrift_globals, ONLY: add_NEO, solve_plasma_equations
      USE booz_params, ONLY: lsurf_boz
      USE vmec_input, ONLY: am_aux_s, am_aux_f, ac_aux_s, ac_aux_f, &
                            pmass_type, pcurr_type, pres_scale, ncurr, curtor
      USE EZspline
      USE EZspline_obj
      USE stel_tools
      USE mpi_params
      USE mpi_inc
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: i, ier
      INTEGER :: bcs0(2)
      REAL(rprec), ALLOCATABLE :: p_restart(:), I_restart(:), s_temp(:)
      TYPE(EZspline1_r8) :: p_spl, i_spl
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      ! Build total pressure and enclosed current from restart arrays
      ALLOCATE(p_restart(nsj), I_restart(nsj))
      DO i = 1, nsj
         p_restart(i) = SUM(DENS_RESTART(:,i) * TEMP_RESTART(:,i)) * e_charge
      END DO
      I_restart = eq_phiedge / mu0 * UGRID_RESTART

      IF (lvmec) THEN
         bcs0 = (/ 0, 0 /)
         ALLOCATE(s_temp(n_eq))
         FORALL(i = 1:n_eq) s_temp(i) = DBLE(i-1)/DBLE(n_eq-1)

         ! Pressure profile -> AM_AUX_S/F
         CALL EZspline_init(p_spl, nsj, bcs0, ier)
         p_spl%x1        = THRIFT_S
         p_spl%isHermite = 1
         CALL EZspline_setup(p_spl, p_restart, ier, EXACT_DIM=.true.)
         PMASS_TYPE = 'akima_spline'
         PRES_SCALE = one
         DO i = 1, n_eq
            CALL EZspline_interp(p_spl, s_temp(i), AM_AUX_F(i), ier)
            AM_AUX_S(i) = s_temp(i)
         END DO
         CALL EZspline_free(p_spl, ier)

         ! Current profile dI/ds -> AC_AUX_S/F
         CALL EZspline_init(i_spl, nsj, bcs0, ier)
         i_spl%x1        = THRIFT_S
         i_spl%isHermite = 1
         CALL EZspline_setup(i_spl, I_restart, ier, EXACT_DIM=.true.)
         PCURR_TYPE = 'akima_spline_ip'
         NCURR = 1
         CALL EZspline_derivative1_array_r8(i_spl, 1, n_eq, s_temp, AC_AUX_F(1:n_eq), ier)
         CURTOR = I_restart(nsj)
         AC_AUX_S(1:n_eq) = s_temp
         CALL EZspline_free(i_spl, ier)

         DEALLOCATE(s_temp)
      END IF

      DEALLOCATE(p_restart, I_restart)

      IF (lverb) WRITE(6,'(A)') '----- Running equilibrium for restart -----'

      ! Run VMEC and load equilibrium splines (sets iota, vp, bsq, bu, bv, phip, ...)
      lscreen_subcodes = .TRUE.
      proc_string = TRIM(TRIM(id_string) // '.000_000')
      CALL thrift_run_equil
      lscreen_subcodes = .FALSE.

      ! Run booz_xform and dkes so DKES_D** are ready for thrift_penta
      IF (solve_plasma_equations .AND. add_NEO) THEN
         IF (ALLOCATED(lsurf_boz)) DEALLOCATE(lsurf_boz)
         ALLOCATE(lsurf_boz(ns_eq))
         lsurf_boz = .FALSE.
         DO i = 1, DKES_NS_MAX
            IF (DKES_K(i) < 1) CYCLE
            lsurf_boz(DKES_K(i)) = .TRUE.
         END DO
         IF (lverb) WRITE(6,'(A)') '----- Running booz_xform for restart -----'
         CALL thrift_paraexe('booz_xform', proc_string, lscreen_subcodes)
         IF (lverb) WRITE(6,'(A)') '----- Running dkes for restart -----'
         CALL thrift_paraexe('dkes', proc_string, lscreen_subcodes)
      END IF

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_restart_equil
