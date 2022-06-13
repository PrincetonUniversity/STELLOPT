!-----------------------------------------------------------------------
!     Module:        BEAMS3D_IMAS_MODULE
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          09/22/2021
!     Description:   Provides IMAS interfaces for running BEAMS3D from 
!                    a variety of sources.
!-----------------------------------------------------------------------
MODULE BEAMS3D_IMAS_MODULE
#if defined(IMAS)
  !---------------------------------------------------------------------
  !     Libraries
  !---------------------------------------------------------------------
  USE ids_schemas

  !---------------------------------------------------------------------
  !     MODULE VARIABLES
  !        IMAS_EQUI : Equilibrium input
  !        IMAS_MACH : Machine input
  !        IMAS_PROF : Profile input
  !---------------------------------------------------------------------
  TYPE(ids_equilibrium) :: IMAS_EQUI
  TYPE(ids_equilibrium) :: IMAS_MACH
  TYPE(ids_equilibrium) :: IMAS_PROF
   
!-----------------------------------------------------------------------
!     SUBROUTINE:        VMEC_IMAS
!     PURPOSE:           RUNS VMEC ASSUMING XML DEFINES RUN LIKE THE
!                        NAMELIST
!-----------------------------------------------------------------------
SUBROUTINE BEAMS3D_IMAS_NBI(EQ_IN, MACHINE_IN, PROFILES_IN, INDATA, &
                        status_code, status_message)
  !---------------------------------------------------------------------
  !     Libraries
  !---------------------------------------------------------------------
  USE ids_schemas
  USE beams3d_runtime
  USE wall_mod, ONLY: wall_free
  USE mpi_params
  USE mpi_inc
#if defined(LHDF5)
  USE hdf5
#endif

  !---------------------------------------------------------------------
  !     INPUT/OUTPUT VARIABLES
  !        EQ_IN : Equilibrium input
  !        MACHINE_IN : Machine input
  !        PROFILES_IN : Profile input
  !        STATUS_CODE : STATUS Flag
  !        STATUS_MESSAGE : Text Message
  !---------------------------------------------------------------------
  IMPLICIT NONE
  TYPE(ids_equilibrium), INTENT(IN) :: EQ_IN
  TYPE(ids_equilibrium), INTENT(IN) :: MACHINE_IN
  TYPE(ids_equilibrium), INTENT(IN) :: PROFILES_IN
  TYPE(ids_parameters_input), INTENT(IN) :: INDATA
  INTEGER, INTENT(OUT) :: status_code
  CHARACTER(LEN=:), POINTER, INTENT(OUT) :: status_message

  !---------------------------------------------------------------------
  !     SUBROUTINE VARIABLES
  !---------------------------------------------------------------------
  LOGICAL :: lmpi_flag
  INTEGER :: impi_flag

  !---------------------------------------------------------------------
  !     BEGIN EXECUTION
  !---------------------------------------------------------------------

  !----  MPI initialisation ----
#if defined(MPI_OPT)
  CALL MPI_initialized(lmpi_flag, impi_flag)
  if (.not. lmpi_flag)   call MPI_INIT(impi_flag)
#endif

  !----  Handle the MPI Communicators ----
#if defined(MPI_OPT)
  CALL MPI_COMM_DUP( MPI_COMM_WORLD, MPI_COMM_BEAMS, ierr_mpi)
  CALL MPI_COMM_RANK(MPI_COMM_BEAMS, myworkid, ierr_mpi) ! MPI
  CALL MPI_COMM_SIZE(MPI_COMM_BEAMS, nprocs_beams, ierr_mpi) ! MPI
  CALL MPI_COMM_SPLIT_TYPE(MPI_COMM_BEAMS, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, MPI_COMM_SHARMEM, ierr_mpi)
  CALL MPI_COMM_RANK(MPI_COMM_SHARMEM, myid_sharmem, ierr_mpi)
  CALL MPI_COMM_SIZE(MPI_COMM_SHARMEM, nshar, ierr_mpi) ! MPI
  CALL MPI_GET_VERSION(vmajor,vminor,ierr_mpi)
  CALL MPI_GET_LIBRARY_VERSION(mpi_lib_name,liblen,ierr_mpi)
  ! Now we set some info
  CALL MPI_INFO_CREATE(mpi_info_beams3d, ierr_mpi)
  CALL MPI_INFO_SET(mpi_info_beams3d, "IBM_largeblock_io", "true",    ierr_mpi)
  CALL MPI_INFO_SET(mpi_info_beams3d, "stripping_unit",    "1048576", ierr_mpi)
  CALL MPI_INFO_SET(mpi_info_beams3d, "romio_ds_read",     "disable", ierr_mpi)
  CALL MPI_INFO_SET(mpi_info_beams3d, "romio_ds_write",    "disable", ierr_mpi)
#endif

  !----  HDF5  ----
#if defined(LHDF5)
  CALL H5GET_LIBVERSION_F(h5major, h5minor, h5rel, ier)
  h5par = 0
#endif
#if defined(HDF5_PAR)
  h5par = 1
#endif

  !----  Constants  ----
  pi = 4.0 * ATAN(1.0)
  pi2 = 8.0 * ATAN(1.0)
  invpi2 = 1./pi2
  mu0 = (16.0E-7) * ATAN(1.0)
  to3 = REAL(2)/REAL(3)

  !----  Defaults  ----
  lverb        = .FALSE.
  lvmec        = .FALSE.
  lpies        = .FALSE.
  lspec        = .FALSE.
  leqdsk       = .FALSE.
  lcoil        = .FALSE.
  lmgrid       = .FALSE.
  lascot       = .FALSE.
  lascotfl     = .FALSE.
  lascot4      = .FALSE.
  lbbnbi       = .FALSE.
  lraw         = .FALSE.
  lvessel      = .FALSE.
  lvac         = .FALSE.
  lrestart_grid     = .FALSE.
  lrestart_particles     = .FALSE.
  lbeam_simple = .FALSE.
  lhitonly           = .FALSE. 
  lplasma_only = .FALSE.
  lbeam        = .FALSE.
  lread_input  = .FALSE.
  lcollision   = .FALSE.
  lw7x   = .FALSE.
  lrandomize = .FALSE.
  lsuzuki = .FALSE.
  id_string    = TRIM(file_str)
  coil_string  = ''
  mgrid_string = ''
  vessel_string = ''
  restart_string = ''
  eqdsk_string = ''

  !----  Initialize the code  ----
  CALL beams3d_init

  !----  Follow particles  ----
  CALL beams3d_follow

  !----  Write ouput  ----
  CALL beams3d_write('TRAJECTORY_PARTIAL')
  IF (lascot) THEN
     IF (lascotfl) THEN
        CALL beams3d_write_ascoth5('FIELDLINES')
     ELSE
        CALL beams3d_write_ascoth5('MARKER')
     END IF
  END IF
  IF (lascot4) CALL beams3d_write_ascoth4('MARKER')

  !----  Diagnostics ----
  CALL beams3d_diagnostics

  !----  Clean up memory ----
  CALL beams3d_free(MPI_COMM_SHARMEM)
  CALL wall_free(ier,MPI_COMM_SHARMEM)
              
  !----  Free local communicators ----
  CALL MPI_COMM_FREE(MPI_COMM_SHARMEM,ierr_mpi)
  CALL MPI_COMM_FREE(MPI_COMM_BEAMS,ierr_mpi)
  CALL MPI_INFO_FREE(mpi_info_beams3d, ierr_mpi)

  !---------------------------------------------------------------------
  !     END EXECUTION - Don't touch below here
  !---------------------------------------------------------------------

  !----  MPI Finalisation ----
#if defined(MPI_OPT)
  call MPI_finalized(lmpi_flag, impi_flag)
  IF (.NOT. lmpi_flag)   CALL MPI_Finalize(impi_flag)
#endif
  
END SUBROUTINE BEAMS3D_IMAS

#endif
END MODULE BEAMS3D_IMAS_MODULE