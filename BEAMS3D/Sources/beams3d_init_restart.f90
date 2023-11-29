! -----------------------------------------------------------------------
! Module:        beams3d_init_restart
! Authors:       S. Lazerson (lazerson@pppl.gov)
! Date:          03/27/2012
! Description:   This subroutine loads data from a previous run.
! For now we define local variables to control how
! this is done.  In the future we will let the user
! do this interactively.
! -----------------------------------------------------------------------
SUBROUTINE beams3d_init_restart
   ! -----------------------------------------------------------------------
   ! Libraries
   ! -----------------------------------------------------------------------
   USE stel_kinds, ONLY: rprec
   USE beams3d_runtime
   USE beams3d_grid
   USE beams3d_lines
   USE beams3d_physics_mod, ONLY: beams3d_MODB
#if defined(LHDF5)
   USE ez_hdf5
#endif
   USE mpi_params
   USE mpi_inc
   ! -----------------------------------------------------------------------
   ! Local Variables
   ! ier            Error Flag
   ! npoinc_extract Which save state to extract from file.
   ! -----------------------------------------------------------------------
   IMPLICIT NONE
   logical :: lplasma_old, ldepo_old, lfusion_old
   integer :: i, k, ier, npoinc_extract, npoinc_save, state_flag
   LOGICAL, DIMENSION(:), ALLOCATABLE :: lgc2fo_old
   INTEGER, DIMENSION(:), ALLOCATABLE :: beam2, start_dex
   real(rprec) :: vpartmax, B_help, version_old, s_fullorbit
   REAL(rprec), DIMENSION(3) :: q
   REAL(rprec), DIMENSION(:), ALLOCATABLE :: mass2, charge2, Zatom2, &
      weight2
   ! -----------------------------------------------------------------------
   ! Begin Subroutine
   ! -----------------------------------------------------------------------
   IF (lverb) THEN
      WRITE(6,'(A)')  '----- Reading Restart File -----'
      WRITE(6,'(A)')  '   FILE: '// TRIM(restart_string)
   END IF

   IF (myworkid == master) THEN
      ! Save quantities
      lfusion_old = .FALSE.
      npoinc_save = npoinc
      ! Read the data
      CALL open_hdf5(TRIM(restart_string), fid, ier, LCREATE = .false.)
      IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR, TRIM(restart_string), ier)
      CALL read_scalar_hdf5(fid,'nparticles', ier, INTVAR = nparticles)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'nparticles', ier)
      CALL read_scalar_hdf5(fid,'npoinc', ier, INTVAR = npoinc)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'npoinc', ier)
      CALL read_scalar_hdf5(fid,'VERSION', ier, DBLVAR = version_old)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'version_old', ier)
      IF (version_old >= 2.8) THEN
         CALL read_scalar_hdf5(fid,'lfusion', ier, BOOVAR = lfusion_old)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'lfusion', ier)
      END IF
      IF (lverb) THEN
         WRITE(6,'(A,I8)') '   NPARTICLES_OLD: ', nparticles
         WRITE(6,'(A,I8)') '   NPOINC_OLD: ', npoinc
         IF (lfusion_old) WRITE(6,'(A,I8)') '   FUSION RUN DETECTED'
      END IF
      ! CALL read_scalar_hdf5(fid,'partvmax',ier,DBLVAR=vpartmax)
      ! partvmax = MAX(partvmax,vpartmax)
      ! IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'partvmax',ier)
      IF (ALLOCATED(mass)) DEALLOCATE(mass)
      IF (ALLOCATED(charge)) DEALLOCATE(charge)
      IF (ALLOCATED(Zatom)) DEALLOCATE(charge)
      IF (ALLOCATED(beam)) DEALLOCATE(beam)
      IF (ALLOCATED(weight)) DEALLOCATE(weight)
      IF (ALLOCATED(end_state)) DEALLOCATE(end_state)
      IF (ALLOCATED(R_lines)) DEALLOCATE(R_lines)
      IF (ALLOCATED(PHI_lines)) DEALLOCATE(PHI_lines)
      IF (ALLOCATED(Z_lines)) DEALLOCATE(Z_lines)
      IF (ALLOCATED(vr_lines)) DEALLOCATE(vr_lines)
      IF (ALLOCATED(vphi_lines)) DEALLOCATE(vphi_lines)
      IF (ALLOCATED(vz_lines)) DEALLOCATE(vz_lines)
      IF (ALLOCATED(moment_lines)) DEALLOCATE(moment_lines)
      IF (ALLOCATED(vll_lines)) DEALLOCATE(vll_lines)
      IF (ALLOCATED(neut_lines)) DEALLOCATE(neut_lines)
      ALLOCATE(mass2(nparticles), charge2(nparticles), Zatom2(nparticles),&
         beam2(nparticles), weight2(nparticles), end_state(nparticles), lgc2fo_old(nparticles))
      lgc2fo_old(:) = .TRUE.
      ALLOCATE(R_lines(0:npoinc, nparticles), Z_lines(0:npoinc, nparticles), PHI_lines(0:npoinc, nparticles),&
         vll_lines(0:npoinc, nparticles), neut_lines(0:npoinc, nparticles), moment_lines(0:npoinc, nparticles),&
         S_lines(0:npoinc, nparticles), B_lines(0:npoinc, nparticles), &
         vr_lines(0:npoinc, nparticles), vphi_lines(0:npoinc, nparticles), vz_lines(0:npoinc, nparticles))
      CALL read_var_hdf5(fid,'mass', nparticles, ier, DBLVAR = mass2)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'mass2', ier)
      CALL read_var_hdf5(fid,'charge', nparticles, ier, DBLVAR = charge2)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'charge2', ier)
      CALL read_var_hdf5(fid,'Zatom', nparticles, ier, DBLVAR = Zatom2)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'Zatom2', ier)
      CALL read_var_hdf5(fid,'Weight', nparticles, ier, DBLVAR = weight2)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'weight2', ier)
      CALL read_var_hdf5(fid,'Beam', nparticles, ier, INTVAR = beam2)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'beam2', ier)
      CALL read_var_hdf5(fid,'end_state', nparticles, ier, INTVAR = end_state)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'end_state', ier)
      CALL read_var_hdf5(fid,'R_lines', npoinc + 1, nparticles, ier, DBLVAR = R_lines)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'R_lines', ier)
      CALL read_var_hdf5(fid,'Z_lines', npoinc + 1, nparticles, ier, DBLVAR = Z_lines)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'Z_lines', ier)
      CALL read_var_hdf5(fid,'PHI_lines', npoinc + 1, nparticles, ier, DBLVAR = PHI_lines)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'PHI_lines', ier)
      CALL read_var_hdf5(fid,'vll_lines', npoinc + 1, nparticles, ier, DBLVAR = vll_lines)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'vll_lines', ier)
      CALL read_var_hdf5(fid,'neut_lines', npoinc + 1, nparticles, ier, BOOVAR = neut_lines)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'neut_lines', ier)
      CALL read_var_hdf5(fid,'moment_lines', npoinc + 1, nparticles, ier, DBLVAR = moment_lines)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'moment_lines', ier)
      CALL read_var_hdf5(fid,'S_lines', npoinc + 1, nparticles, ier, DBLVAR = S_lines)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'S_lines', ier)
      CALL read_var_hdf5(fid,'B_lines', npoinc + 1, nparticles, ier, DBLVAR = B_lines)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'B_lines', ier)
      IF (version_old >= 4.0) THEN
         CALL read_var_hdf5(fid,'vr_lines', npoinc + 1, nparticles, ier, DBLVAR = vr_lines)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'vr_lines', ier)
         CALL read_var_hdf5(fid,'vphi_lines', npoinc + 1, nparticles, ier, DBLVAR = vphi_lines)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'vphi_lines', ier)
         CALL read_var_hdf5(fid,'vz_lines', npoinc + 1, nparticles, ier, DBLVAR = vz_lines)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'vz_lines', ier)
      END IF
      CALL close_hdf5(fid, ier)
      IF (ier /= 0) CALL handle_err(HDF5_CLOSE_ERR,'beams3d_'// TRIM(restart_string) //'.h5', ier)


      ! Helper for where to start loading particles
      ALLOCATE(start_dex(nparticles))
      start_dex = 0 ! Default

      ! Decide what to do
      ! IF lfusion_old then start from 0
      ! IF ldepo run then start from 2
      ! ELSE Start from wall_hit
      ldepo_old = .false.
      state_flag = 0
      IF (ANY(end_state == 3)) ldepo_old = .true.
      IF (lfusion_old) THEN
         WRITE(6,'(A)') '   Detected old fusion run! '
         end_state = 0
         state_flag = 0
         IF (lplasma_only) THEN
            WHERE(S_lines(0,:) >= 1) end_state = - 1
         END IF
      ELSEIF (ldepo_old) THEN
         WRITE(6,'(A)') '   Detected old deposition run! '
         ! Only orbiting particles
         state_flag = 0
         ! Use ionization point unless outside FO radius
         start_dex = 2
         s_fullorbit = SIGN(rho_fullorbit * rho_fullorbit, rho_fullorbit)
         WHERE(S_lines(1,:) >= s_fullorbit)
         start_dex = 1
         lgc2fo_old = .FALSE.
         END WHERE
         ! IF plasma run only consider particles born inside equilibrium
         IF (lplasma_only) THEN
            WRITE(6,'(A)') '   Detected old plasma only run! '
            WHERE((S_lines(1,:) >= 1) .and. (start_dex == 1)) end_state = - 1
            WHERE((S_lines(2,:) >= 1) .and. (start_dex == 2)) end_state = - 1
         END IF
      ELSE
         WRITE(6,'(A)') '   Detected old restart run! Using wall strikes for restart. '
         state_flag = 2
         DO i = 1, nparticles
            start_dex(i) = COUNT(R_lines(:, i) > 0) - 1 ! Note indexed from 0
         END DO
      END IF
      k = COUNT(end_state == state_flag)

      IF (lverb) THEN
         WRITE(6,'(A,I8)') '   NPARTICLES_RESTART: ', k
      END IF

      ! Allocate the particles
      ALLOCATE(  R_start(k), phi_start(k), Z_start(k), &
         vr_start(k), vphi_start(k), vz_start(k), &
         mass(k), charge(k), &
         mu_start(k), Zatom(k), t_end(k), vll_start(k), &
         beam(k), weight(k), lgc2fo_start(k) )

      ! Now fill the arrays downselecting for non-shinethrough particles
      k = 1
      DO i = 1, nparticles
         IF (end_state(i) /= state_flag) CYCLE
         npoinc_extract = start_dex(i)
         R_start(k)   = R_lines(npoinc_extract, i)
         Z_start(k)   = Z_lines(npoinc_extract, i)
         phi_start(k) = PHI_lines(npoinc_extract, i)
         vll_start(k) = vll_lines(npoinc_extract, i)
         vr_start(k) = vr_lines(npoinc_extract, i)
         vphi_start(k) = vphi_lines(npoinc_extract, i)
         vz_start(k)  = vz_lines(npoinc_extract, i)
         mu_start(k)  = moment_lines(npoinc_extract, i) * B_lines(npoinc_extract, i)
         mass(k)      = mass2(i)
         charge(k)   = charge2(i)
         Zatom(k)    = Zatom2(i)
         beam(k)     = beam2(i)
         weight(k)   = weight2(i)
         lgc2fo_start(k) = lgc2fo_old(i)
         t_end(k)    = MAXVAL(t_end_in)
         q = (/ R_start(k), phi_start(k), Z_start(k) /)
         CALL beams3d_MODB(q, B_help)
         mu_start(k) = mu_start(k) / B_help
         k = k + 1
      END DO
      DEALLOCATE(R_lines, Z_lines, PHI_lines, vll_lines, moment_lines, &
         neut_lines, end_state, S_lines, B_lines, vr_lines, vphi_lines, vz_lines)
      DEALLOCATE(mass2, charge2, Zatom2, beam2, weight2, start_dex, lgc2fo_old)

      ! Restore quantities
      nparticles = k - 1
      npoinc = npoinc_save
      nbeams = MAXVAL(beam)
      IF (lverb) THEN
         WRITE(6,'(A,I6)') '   # of Beams: ', nbeams
      END IF
   END IF

#if defined(MPI_OPT)
   CALL MPI_BARRIER(MPI_COMM_BEAMS, ierr_mpi)
   CALL MPI_BCAST(nparticles, 1, MPI_INTEGER, master, MPI_COMM_BEAMS, ierr_mpi)
   CALL MPI_BCAST(nbeams, 1, MPI_INTEGER, master, MPI_COMM_BEAMS, ierr_mpi)
   CALL MPI_BCAST(partvmax, 1, MPI_REAL8, master, MPI_COMM_BEAMS, ierr_mpi)
   IF (myworkid /= master) THEN
      ALLOCATE(  R_start(nparticles), phi_start(nparticles), Z_start(nparticles), &
         vr_start(nparticles), vphi_start(nparticles), vz_start(nparticles), &
         mass(nparticles), charge(nparticles), &
         mu_start(nparticles), Zatom(nparticles), t_end(nparticles), vll_start(nparticles), &
         beam(nparticles), weight(nparticles), lgc2fo_start(nparticles) )
   END IF
   CALL MPI_BCAST(lgc2fo_start, nparticles, MPI_LOGICAL, master, MPI_COMM_BEAMS, ierr_mpi)
   CALL MPI_BCAST(mu_start, nparticles, MPI_REAL8, master, MPI_COMM_BEAMS, ierr_mpi)
   CALL MPI_BCAST(t_end, nparticles, MPI_REAL8, master, MPI_COMM_BEAMS, ierr_mpi)
   CALL MPI_BCAST(mass, nparticles, MPI_REAL8, master, MPI_COMM_BEAMS, ierr_mpi)
   CALL MPI_BCAST(charge, nparticles, MPI_REAL8, master, MPI_COMM_BEAMS, ierr_mpi)
   CALL MPI_BCAST(Zatom, nparticles, MPI_REAL8, master, MPI_COMM_BEAMS, ierr_mpi)
   CALL MPI_BCAST(weight, nparticles, MPI_REAL8, master, MPI_COMM_BEAMS, ierr_mpi)
   CALL MPI_BCAST(R_start, nparticles, MPI_REAL8, master, MPI_COMM_BEAMS, ierr_mpi)
   CALL MPI_BCAST(phi_start, nparticles, MPI_REAL8, master, MPI_COMM_BEAMS, ierr_mpi)
   CALL MPI_BCAST(Z_start, nparticles, MPI_REAL8, master, MPI_COMM_BEAMS, ierr_mpi)
   CALL MPI_BCAST(vr_start, nparticles, MPI_REAL8, master, MPI_COMM_BEAMS, ierr_mpi)
   CALL MPI_BCAST(vphi_start, nparticles, MPI_REAL8, master, MPI_COMM_BEAMS, ierr_mpi)
   CALL MPI_BCAST(vz_start, nparticles, MPI_REAL8, master, MPI_COMM_BEAMS, ierr_mpi)
   CALL MPI_BCAST(vll_start, nparticles, MPI_REAL8, master, MPI_COMM_BEAMS, ierr_mpi)
   CALL MPI_BCAST(beam, nparticles, MPI_INTEGER, master, MPI_COMM_BEAMS, ierr_mpi)
#endif

   RETURN

   ! -----------------------------------------------------------------------
   ! End Subroutine
   ! -----------------------------------------------------------------------
END SUBROUTINE beams3d_init_restart
