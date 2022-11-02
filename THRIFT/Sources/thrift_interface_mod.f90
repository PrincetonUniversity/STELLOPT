!-----------------------------------------------------------------------
!     Module:        THRIFT_INTERFACE_MOD
!     Authors:       L. van Ham, S. Lazerson
!     Date:          11/XX/2022
!     Description:   Module provides routines for handling various
!                    initialization tasks.
!-----------------------------------------------------------------------
MODULE THRIFT_INTERFACE_MOD
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE thrift_runtime
      USE mpi_params
      USE mpi_inc
#if defined(LHDF5)
      USE hdf5
#endif

!-----------------------------------------------------------------------
!     Module Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: vmajor, vminor, liblen
      INTEGER :: h5major, h5minor, h5rel, h5par
      INTEGER :: mpi_info_thrift
      CHARACTER(LEN=MPI_MAX_LIBRARY_VERSION_STRING) :: mpi_lib_name
      
!-----------------------------------------------------------------------
!     Subroutines
!         thrift_init_mpi:         MPI initialization
!         thrift_init_HDF5:        HDF5 initialization
!         thrift_init_constants:   Constants initialization
!         thrift_init_commandline: Handle the command line
!-----------------------------------------------------------------------
CONTAINS

      SUBROUTINE thrift_init_mpi
      IMPLICIT NONE
      myworkid = master
#if defined(MPI_OPT)
      CALL MPI_INIT(ierr_mpi) ! MPI
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_INIT_ERR, 'thrift_main', ierr_mpi)
      CALL MPI_COMM_DUP( MPI_COMM_WORLD, MPI_COMM_THRIFT, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_RANK_ERR, 'thrift_main', ierr_mpi)
      CALL MPI_COMM_RANK(MPI_COMM_THRIFT, myworkid, ierr_mpi) ! MPI
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_RANK_ERR, 'thrift_main', ierr_mpi)
      CALL MPI_COMM_SIZE(MPI_COMM_THRIFT, nprocs_thrift, ierr_mpi) ! MPI
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_SIZE_ERR, 'thrift_main', ierr_mpi)
      CALL thrift_init_mpi_split(MPI_COMM_THRIFT)
      CALL MPI_GET_VERSION(vmajor,vminor,ierr_mpi)
      CALL MPI_GET_LIBRARY_VERSION(mpi_lib_name,liblen,ierr_mpi)
      ! Now we set some info
      CALL MPI_INFO_CREATE(mpi_info_thrift, ierr_mpi)
      CALL MPI_INFO_SET(mpi_info_thrift, "IBM_largeblock_io", "true",    ierr_mpi)
      CALL MPI_INFO_SET(mpi_info_thrift, "stripping_unit",    "1048576", ierr_mpi)
      CALL MPI_INFO_SET(mpi_info_thrift, "romio_ds_read",     "disable", ierr_mpi)
      CALL MPI_INFO_SET(mpi_info_thrift, "romio_ds_write",    "disable", ierr_mpi)
#endif
      RETURN
      END SUBROUTINE thrift_init_mpi

      SUBROUTINE thrift_init_mpi_split(comm)
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: comm
#if defined(MPI_OPT)
      CALL MPI_COMM_SPLIT_TYPE(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, MPI_COMM_SHARMEM, ierr_mpi)
      CALL MPI_COMM_RANK(MPI_COMM_SHARMEM, myid_sharmem, ierr_mpi)
#endif
      END SUBROUTINE thrift_init_mpi_split

      SUBROUTINE thrift_cleanup
      IMPLICIT NONE
      INTEGER :: ier
      ! Clean up
      ier = 0
      !CALL thrift_free(MPI_COMM_SHARMEM)
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_THRIFT, ierr_mpi)
      IF (ierr_mpi /= 0) CALL handle_err(MPI_BARRIER_ERR, 'thrift_main', ierr_mpi)
      ierr_mpi=0
      CALL MPI_INFO_FREE(mpi_info_thrift, ierr_mpi)
      CALL MPI_FINALIZE(ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_FINE_ERR, 'thrift_main', ierr_mpi)
#endif
      IF (lverb) WRITE(6, '(A)') '----- THRIFT DONE -----'
      RETURN
      END SUBROUTINE

      SUBROUTINE thrift_init_hdf5
      IMPLICIT NONE
      INTEGER :: ier
#if defined(LHDF5)
      CALL H5GET_LIBVERSION_F(h5major, h5minor, h5rel, ier)
      h5par = 0
#endif
#if defined(HDF5_PAR)
      h5par = 1
#endif
      RETURN
      END SUBROUTINE thrift_init_hdf5

      SUBROUTINE thrift_init_constants
      IMPLICIT NONE
      pi = 4.0 * ATAN(1.0)
      pi2 = 8.0 * ATAN(1.0)
      invpi2 = 1./pi2
      mu0 = (16.0E-7) * ATAN(1.0)
      to3 = REAL(2)/REAL(3)
      lverb = .true.
      lread_input = .true.
      RETURN
      END SUBROUTINE thrift_init_constants

      SUBROUTINE thrift_init_commandline
      IMPLICIT NONE
      INTEGER :: numargs, i, ier
      INTEGER, parameter :: arg_len = 256
      CHARACTER*(arg_len) :: arg1
      CHARACTER*(arg_len), allocatable, dimension(:) :: args
      lverb = .false.
      IF (myworkid == master) THEN
        numargs = 0
        i = 0
        arg1 = ''
        limas = .false.
        lverb = .true.
        lvmec = .false.
        id_string = ''

        ! First Handle the input arguments
        CALL GETCARG(1, arg1, numargs)
        ALLOCATE(args(numargs))
        ! Cycle through Arguments
        i = 1
        DO WHILE (i <= numargs)
            call GETCARG(i, args(i), numargs)
            select case (args(i))
            case ("-noverb") ! No Verbose Output
                lverb = .false.
            case ("-vmec")
                i = i + 1
                lvmec = .true.
                CALL GETCARG(i, id_string, numargs)
            case ("-prof")
                i = i + 1
                lvmec = .true.
                CALL GETCARG(i, prof_string, numargs)
            case ("-help", "-h") ! Output Help message
                write(6, *) ' Beam MC Code'
                write(6, *) ' Usage: xthrift <options>'
                write(6, *) '    <options>'
                write(6, *) '     -vmec ext:     VMEC input/wout extension'
                write(6, *) '     -prof file:    Profile file'
                write(6, *) '     -noverb:       Supress all screen output'
                write(6, *) '     -help:         Output help message'
            end select
            i = i + 1
        END DO
        DEALLOCATE(args)
      END IF
    ! Broadcast variables
#if defined(MPI_OPT)
      CALL MPI_BCAST(id_string, 256, MPI_CHARACTER, master, MPI_COMM_THRIFT, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'thrift_main', ierr_mpi)
      CALL MPI_BCAST(prof_string, 256, MPI_CHARACTER, master, MPI_COMM_THRIFT, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'thrift_main', ierr_mpi)
      CALL MPI_BCAST(lvmec, 1, MPI_LOGICAL, master, MPI_COMM_THRIFT, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'thrift_main', ierr_mpi)
      CALL MPI_BCAST(limas, 1, MPI_LOGICAL, master, MPI_COMM_THRIFT, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'thrift_main', ierr_mpi)
#endif
      RETURN
      END SUBROUTINE thrift_init_commandline

      SUBROUTINE thrift_output_header
      IMPLICIT NONE
      INTEGER :: nshar
#if defined(MPI_OPT)
      CALL MPI_COMM_SIZE(MPI_COMM_SHARMEM, nshar, ierr_mpi) ! MPI
#endif
      IF (lverb) THEN
        WRITE(6, '(a,f5.2)') 'THRIFT Version ', THRIFT_VERSION
#if defined(LHDF5)
        IF (h5par > 0) THEN
           WRITE(6,'(A)')      '-----  HDF5 (Parallel) Parameters  -----'
        ELSE
           WRITE(6,'(A)')      '-----  HDF5 Parameters  -----'
        ENDIF
        WRITE(6,'(A,I2,2(A,I2.2))')  '   HDF5_version:  ', h5major,'.',h5minor,' release: ',h5rel
#endif
        WRITE(6,'(A)')      '-----  MPI Parameters  -----'
        WRITE(6,'(A,I2,A,I2.2)')  '   MPI_version:  ', vmajor,'.',vminor
        WRITE(6,'(A,A)')  '   ', TRIM(mpi_lib_name(1:liblen))
        WRITE(6,'(A,I8)')  '   Nproc_total:  ', nprocs_thrift
        WRITE(6,'(A,3X,I5)')  '   Nproc_shared: ', nshar
      END IF
      RETURN
      END SUBROUTINE thrift_output_header

END MODULE THRIFT_INTERFACE_MOD
