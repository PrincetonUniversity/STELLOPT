!-----------------------------------------------------------------------
!     Module:        beams3d_init_ascot5_endstate
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          01/23/2026
!     Description:   This subroutine loads particle data from an ASCOT5
!                    endstate run.
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_init_ascot5_endstate
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE beams3d_globals, ONLY: a5_run_name, a5_marker_name
      USE beams3d_runtime
      USE beams3d_grid
      USE beams3d_lines
#if defined(LHDF5)
      USE ez_hdf5
#endif
      USE mpi_sharmem
      USE mpi_params
      USE mpi_inc
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          npoinc_extract Which save state to extract from file.
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: i, k, ier
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: temp2
      CHARACTER(10)  :: marker_id
      INTEGER :: MPI_COMM_LOCAL
      DOUBLE PRECISION, PARAMETER :: dalton    = 1.66053906892E-27 ! AMU [kg]
      DOUBLE PRECISION, PARAMETER :: e_charge  = 1.60217662E-19    ![C]
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------


      ! DEALLOCATE the start variables
#if defined(MPI_OPT)
      CALL mpidealloc(R_start, win_R_start)
      CALL mpidealloc(PHI_start, win_PHI_start)
      CALL mpidealloc(Z_start, win_Z_start)
      CALL mpidealloc(vr_start, win_vr_start)
      CALL mpidealloc(vphi_start, win_vphi_start)
      CALL mpidealloc(vz_start, win_vz_start)
      CALL mpidealloc(mass, win_mass)
      CALL mpidealloc(charge, win_charge)
      CALL mpidealloc(mu_start, win_mu_start)
      CALL mpidealloc(Zatom, win_Zatom)
      CALL mpidealloc(t_end, win_t_end)
      CALL mpidealloc(vll_start, win_vll_start)
      CALL mpidealloc(beam, win_beam)
      CALL mpidealloc(weight, win_weight)
      CALL mpidealloc(lgc2fo_start, win_lgc2fo_start)
      CALL mpidealloc(end_state, win_end_state)
#endif

      IF (lverb) THEN
         WRITE(6,'(A)')  '----- Reading Endstates from ASCOT5 File -----'
         WRITE(6,'(A)')  '               FILE: '//TRIM(ascot_endstate_file)
      END IF

      ! Load the array size information from the old file
      ! Now the strcuture of the files is
      IF (myworkid == master) THEN
         ! Read the data
         CALL open_hdf5(TRIM(ascot_endstate_file),fid,ier,LCREATE=.false.)
         IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,TRIM(restart_string),ier)
         ! First get the marker ID for the run
         !CALL read_group_attscalar_hdf5(fid,'/results/'//TRIM(a5_run_name),'qid_marker',ier,STRVAR=marker_id)
         !IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'qid_marker',ier)
         !WRITE(6,*) 'id: '//TRIM(marker_id)
         ! Now get the numer of markers
         !WRITE(marker_name,'(A,i10.10)') 'prt_',marker_id
         !marker_name = 'prt_'//TRIM(marker_id)
         CALL read_scalar_hdf5(fid,'/marker/'//TRIM(a5_marker_name)//'/n',ier,INTVAR=nparticles)
         nparticles = nparticles*2
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'nparticles',ier)
         IF (lverb) THEN
            WRITE(6,'(A,A)')  '      RUN_NAME: ', TRIM(a5_run_name)
            WRITE(6,'(A,A)')  '   MARKER_NAME: ', TRIM(a5_marker_name)
            WRITE(6,'(A,I8)') '    NPARTICLES: ', nparticles
         END IF
      END IF      
#if defined(MPI_OPT)
      CALL MPI_BCAST(nparticles, 1, MPI_INTEGER, master, MPI_COMM_BEAMS,ierr_mpi)

      ! Allocate the particles
      CALL mpialloc(end_state,    nparticles, myid_sharmem, 0, MPI_COMM_SHARMEM, win_end_state)
      CALL mpialloc(R_start,      nparticles, myid_sharmem, 0, MPI_COMM_SHARMEM, win_R_start)
      CALL mpialloc(PHI_start,    nparticles, myid_sharmem, 0, MPI_COMM_SHARMEM, win_PHI_start)
      CALL mpialloc(Z_start,      nparticles, myid_sharmem, 0, MPI_COMM_SHARMEM, win_Z_start)
      CALL mpialloc(vr_start,     nparticles, myid_sharmem, 0, MPI_COMM_SHARMEM, win_vr_start)
      CALL mpialloc(vphi_start,   nparticles, myid_sharmem, 0, MPI_COMM_SHARMEM, win_vphi_start)
      CALL mpialloc(vz_start,     nparticles, myid_sharmem, 0, MPI_COMM_SHARMEM, win_vz_start)
      CALL mpialloc(mass,         nparticles, myid_sharmem, 0, MPI_COMM_SHARMEM, win_mass)
      CALL mpialloc(charge,       nparticles, myid_sharmem, 0, MPI_COMM_SHARMEM, win_charge)
      CALL mpialloc(mu_start,     nparticles, myid_sharmem, 0, MPI_COMM_SHARMEM, win_mu_start)
      CALL mpialloc(Zatom,        nparticles, myid_sharmem, 0, MPI_COMM_SHARMEM, win_Zatom)
      CALL mpialloc(t_end,        nparticles, myid_sharmem, 0, MPI_COMM_SHARMEM, win_t_end)
      CALL mpialloc(vll_start,    nparticles, myid_sharmem, 0, MPI_COMM_SHARMEM, win_vll_start)
      CALL mpialloc(beam,         nparticles, myid_sharmem, 0, MPI_COMM_SHARMEM, win_beam)
      CALL mpialloc(weight,       nparticles, myid_sharmem, 0, MPI_COMM_SHARMEM, win_weight)
      CALL mpialloc(lgc2fo_start, nparticles, myid_sharmem, 0, MPI_COMM_SHARMEM, win_lgc2fo_start)
#endif
      IF (myworkid == master) THEN
         ALLOCATE(temp2(nparticles))
         beam = 1
         end_state = 0
         lgc2fo_start = .false.
         t_end        = MAXVAL(t_end_in)
         CALL read_var_hdf5(fid,'/results/'//TRIM(a5_run_name)//'/endstate/mass',nparticles,ier,DBLVAR=mass)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'mass',ier)
         mass = mass * dalton
         CALL read_var_hdf5(fid,'/results/'//TRIM(a5_run_name)//'/endstate/charge',nparticles,ier,DBLVAR=charge)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'charge',ier)
         charge = charge * e_charge
         CALL read_var_hdf5(fid,'/results/'//TRIM(a5_run_name)//'/endstate/znum',nparticles,ier,DBLVAR=Zatom)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'Zatom',ier)
         CALL read_var_hdf5(fid,'/results/'//TRIM(a5_run_name)//'/endstate/weight',nparticles,ier,DBLVAR=weight)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'weight',ier)
         CALL read_var_hdf5(fid,'/results/'//TRIM(a5_run_name)//'/endstate/rprt',nparticles,ier,DBLVAR=temp2)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'rprt',ier)
         R_start = temp2
         CALL read_var_hdf5(fid,'/results/'//TRIM(a5_run_name)//'/endstate/zprt',nparticles,ier,DBLVAR=temp2)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'zprt',ier)
         Z_start = temp2
         CALL read_var_hdf5(fid,'/results/'//TRIM(a5_run_name)//'/endstate/phiprt',nparticles,ier,DBLVAR=temp2)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'phiprt',ier)
         PHI_start = temp2 * pi2 / 360.0
         CALL read_var_hdf5(fid,'/results/'//TRIM(a5_run_name)//'/endstate/ppar',nparticles,ier,DBLVAR=temp2)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'ppar',ier)
         vll_start = temp2/mass
         CALL read_var_hdf5(fid,'/results/'//TRIM(a5_run_name)//'/endstate/mu',nparticles,ier,DBLVAR=temp2)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'mu',ier)
         mu_start = temp2 * e_charge
         CALL read_var_hdf5(fid,'/results/'//TRIM(a5_run_name)//'/endstate/prprt',nparticles,ier,DBLVAR=temp2)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'prprt',ier)
         vr_start = temp2/mass
         CALL read_var_hdf5(fid,'/results/'//TRIM(a5_run_name)//'/endstate/pphiprt',nparticles,ier,DBLVAR=temp2)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'pphiprt',ier)
         vphi_start = temp2/mass
         CALL read_var_hdf5(fid,'/results/'//TRIM(a5_run_name)//'/endstate/pzprt',nparticles,ier,DBLVAR=temp2)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'pzprt',ier)
         vz_start = temp2/mass
         CALL close_hdf5(fid,ier)
         IF (ier /= 0) CALL handle_err(HDF5_CLOSE_ERR,'beams3d_'//TRIM(restart_string)//'.h5',ier)
         DEALLOCATE(temp2)
      END IF

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_BEAMS,ierr_mpi)
      i = MPI_UNDEFINED
      IF (myid_sharmem == master) i = 0
      CALL MPI_COMM_SPLIT( MPI_COMM_BEAMS,i,myworkid,MPI_COMM_LOCAL,ierr_mpi)
      IF (myid_sharmem == master) THEN
         CALL MPI_BCAST(     R_start, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(   PHI_start, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(     Z_start, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(   vll_start, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(    mu_start, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(    vr_start, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(  vphi_start, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(    vz_start, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(       t_end, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(        mass, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(      charge, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(       Zatom, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(      weight, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(        beam, nparticles, MPI_INTEGER, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(lgc2fo_start, nparticles, MPI_LOGICAL, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_COMM_FREE(MPI_COMM_LOCAL,ierr_mpi)
      END IF
      nbeams = 1
#endif

      RETURN

!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE beams3d_init_ascot5_endstate
