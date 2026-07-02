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
      USE thrift_equil, ONLY: eq_Aminor, eq_phiedge, vp_spl, bsq_spl, bcs1
      USE thrift_globals, ONLY: nsj, tstart, solve_plasma_equations
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
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER     :: ier, ns_restart, ntimesteps_restart, k
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

      ! Allocate restart density/temperature arrays (all ranks)
      CALL mpialloc(DENS_RESTART,           nion_prof_restart+1, nsj, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_dens_restart)
      CALL mpialloc(TEMP_RESTART,           nion_prof_restart+1, nsj, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_temp_restart)
      CALL mpialloc(DENS_FAST_ALPHAS_RESTART,                nsj, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_dens_fast_alphas_restart)

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

         CALL read_var_hdf5(fid,'THRIFT_VP',nsj,ntimesteps_restart,ier,DBLVAR=temp2d)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'THRIFT_VP',ier)

         bcs1 = (/ 0, 0 /)
         IF (EZspline_allocated(vp_spl)) CALL EZspline_free(vp_spl,ier)
         CALL EZspline_init(vp_spl,nsj,bcs1,ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'thrift_restart: vp_spl',ier)
         vp_spl%isHermite = 0
         FORALL (k=1:nsj) vp_spl%x1(k) = sqrt(DBLE(k-1)/DBLE(nsj-1))
         CALL EZspline_setup(vp_spl,temp2d(:,ntimesteps_restart)/eq_phiedge,ier,EXACT_DIM=.true.)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'thrift_restart: vp_spl',ier)

         CALL read_var_hdf5(fid,'THRIFT_BSQAV',nsj,ntimesteps_restart,ier,DBLVAR=temp2d)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'THRIFT_BSQAV',ier)

         bcs1 = (/ 0, 0 /)
         IF (EZspline_allocated(bsq_spl)) CALL EZspline_free(bsq_spl,ier)
         CALL EZspline_init(bsq_spl,nsj,bcs1,ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'thrift_restart: bsq_spl',ier)
         bsq_spl%isHermite = 0
         FORALL (k=1:nsj) bsq_spl%x1(k) = sqrt(DBLE(k-1)/DBLE(nsj-1))
         CALL EZspline_setup(bsq_spl,temp2d(:,ntimesteps_restart),ier,EXACT_DIM=.true.)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'thrift_restart: bsq_spl',ier)

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
