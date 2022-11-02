!-----------------------------------------------------------------------
!     Module:        thrift_runtime
!     Authors:       L. van Ham
!     Date:          11/XX/2022
!     Description:   This module contains various runtime parameters
!                    for the THRIFT code.  It also contains the
!                    handle_error subroutine.
!     v0.00 11/XX/22 - Generally used to track major version information
!-----------------------------------------------------------------------
MODULE thrift_runtime
    !-----------------------------------------------------------------------
    !     Libraries
    !-----------------------------------------------------------------------
    USE stel_kinds, ONLY: rprec
    USE mpi_params
    USE EZspline
    !-----------------------------------------------------------------------
    !     Module Variables
    !          lverb         Logical to control screen output
    !----------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER, PARAMETER :: MPI_CHECK = 0
    INTEGER, PARAMETER :: FILE_OPEN_ERR = 1
    INTEGER, PARAMETER :: ALLOC_ERR = 11
    INTEGER, PARAMETER :: NAMELIST_READ_ERR = 12
    INTEGER, PARAMETER :: BAD_INPUT_ERR = 13
    INTEGER, PARAMETER :: BAD_BEAMDEX_ERR = 14
    INTEGER, PARAMETER :: VMEC_INPUT_ERR = 2
    INTEGER, PARAMETER :: VMEC_WOUT_ERR = 21
    INTEGER, PARAMETER :: MGRID_ERR = 22
    INTEGER, PARAMETER :: INIT_BEAMS_ERR = 23
    INTEGER, PARAMETER :: ADAS_ERR = 24
    INTEGER, PARAMETER :: BOOZER_ERR = 31
    INTEGER, PARAMETER :: EZSPLINE_ERR = 4
    INTEGER, PARAMETER :: NETCDF_OPEN_ERR = 5
    INTEGER, PARAMETER :: PIES_NETCDF_READ_ERR = 51
    INTEGER, PARAMETER :: HDF5_ERR = 6
    INTEGER, PARAMETER :: HDF5_OPEN_ERR = 61
    INTEGER, PARAMETER :: HDF5_WRITE_ERR = 62
    INTEGER, PARAMETER :: HDF5_READ_ERR = 63
    INTEGER, PARAMETER :: HDF5_CLOSE_ERR = 69
    INTEGER, PARAMETER :: NAG_ERR = 70
    INTEGER, PARAMETER :: D02CJF_ERR = 71
    INTEGER, PARAMETER :: D05ADF_ERR = 72
    INTEGER, PARAMETER :: C05AJF_ERR = 73
    INTEGER, PARAMETER :: LSODE_ERR = 791
    INTEGER, PARAMETER :: RKH68_ERR = 792
    INTEGER, PARAMETER :: MPI_ERR = 80
    INTEGER, PARAMETER :: MPI_INIT_ERR = 801
    INTEGER, PARAMETER :: MPI_RANK_ERR = 802
    INTEGER, PARAMETER :: MPI_SIZE_ERR = 803
    INTEGER, PARAMETER :: MPI_BARRIER_ERR = 81
    INTEGER, PARAMETER :: MPI_SEND_ERR = 821
    INTEGER, PARAMETER :: MPI_RECV_ERR = 822
    INTEGER, PARAMETER :: MPI_REDU_ERR = 823
    INTEGER, PARAMETER :: MPI_BCAST_ERR = 83
    INTEGER, PARAMETER :: MPI_FINE_ERR = 89

    INTEGER, PARAMETER :: MAXPARTICLES = 2**18
    INTEGER, PARAMETER :: MAXBEAMS = 32
    INTEGER, PARAMETER :: MAXPROFLEN = 512
    INTEGER, PARAMETER :: NION = 4

    DOUBLE PRECISION, PARAMETER :: one           = 1.0D0 ! 1.0

    LOGICAL :: lverb, lvmec, lread_input, limas
    INTEGER :: nprocs_thrift, nparallel_runs, mboz, nboz, ier_paraexe, &
               mytimestep
    REAL(rprec) :: pi, pi2, invpi2, mu0, to3
    CHARACTER(256) :: id_string, prof_string, bootstrap_type, proc_string

    REAL(rprec), PARAMETER :: THRIFT_VERSION = 0.00 
    !-----------------------------------------------------------------------
    !     Subroutines
    !          handle_err  Controls Program Termination
    !-----------------------------------------------------------------------
CONTAINS

    SUBROUTINE handle_err(error_num, string_val, ierr)
        USE mpi_params
        USE mpi_inc
        IMPLICIT NONE
        INTEGER, INTENT(in) :: error_num
        INTEGER, INTENT(in) :: ierr
        INTEGER, ALLOCATABLE :: error_array(:)
        CHARACTER(*), INTENT(in) :: string_val
        IF (error_num .ne. MPI_CHECK) WRITE(6, *) '!!!!! ERROR !!!!!'

        IF (error_num .eq. FILE_OPEN_ERR) THEN
            WRITE(6, *) '  THRIFT COULD NOT OPEN A FILE.'
            WRITE(6, *) '  FILENAME: ', TRIM(string_val)
            WRITE(6, *) '  IERR:     ', ierr
        ELSEIF (error_num .eq. VMEC_INPUT_ERR) THEN
            WRITE(6, *) '  THRIFT COULD NOT READ THE VMEC INDATA NAMELIST'
            WRITE(6, *) '  FILENAME: ', TRIM(string_val)
            WRITE(6, *) '  IERR:     ', ierr
        ELSEIF (error_num .eq. VMEC_WOUT_ERR) THEN
            WRITE(6, *) '  THRIFT COULD NOT READ THE VMEC WOUT FILE'
            WRITE(6, *) '  FILENAME: ', TRIM(string_val)
            WRITE(6, *) '  IERR:     ', ierr
        ELSEIF (error_num .eq. ALLOC_ERR) THEN
            WRITE(6, *) '  THRIFT ENCOUNTERED AN ALLOCATION ERROR'
            WRITE(6, *) '  VARIABLES: ', TRIM(string_val)
            WRITE(6, *) '  IERR:      ', ierr
        ELSEIF (error_num .eq. EZSPLINE_ERR) THEN
            WRITE(6, *) '  THRIFT ENCOUNTERED AN EZSPLINE ERROR'
            WRITE(6, *) '  ROUTINE/VAR: ', TRIM(string_val)
            WRITE(6, *) '  IERR:        ', ierr
            CALL EZspline_error(ierr)
        ELSEIF (error_num .eq. MGRID_ERR) THEN
            WRITE(6, *) '  THRIFT ENCOUNTERED AN ERROR READING THE MGRID FILE'
            WRITE(6, *) '  VARIABLES: ', TRIM(string_val)
            WRITE(6, *) '  IERR:      ', ierr
        ELSEIF (error_num .eq. INIT_BEAMS_ERR) THEN
            WRITE(6, *) '  THRIFT ENCOUNTERED AN ERROR INITIALIZING THE BEAMS'
            WRITE(6, *) '  VARIABLES: ', TRIM(string_val)
            WRITE(6, *) '  IERR:      ', ierr
        ELSEIF (error_num .eq. ADAS_ERR) THEN
            WRITE(6, *) '  THRIFT ENCOUNTERED AN ERROR CALLING ADAS'
            WRITE(6, *) '  VARIABLES: ', TRIM(string_val)
            WRITE(6, *) '  IERR:      ', ierr
        ELSEIF (error_num .eq. BAD_BEAMDEX_ERR) THEN
            WRITE(6, *) '  THRIFT ENCOUNTERED AN NAMELIST ERROR'
            WRITE(6, *) '    -BEAMLET USED BUT NO DEX_BEAM SET!'
        ELSEIF (error_num .eq. NETCDF_OPEN_ERR) THEN
            WRITE(6, *) '  THRIFT ENCOUNTERED AN ERROR OPENING A NETCDF FILE'
            WRITE(6, *) '  VARIABLES: ', TRIM(string_val)
            WRITE(6, *) '  IERR:      ', ierr
        ELSEIF (error_num .eq. HDF5_OPEN_ERR) THEN
            WRITE(6, *) '  THRIFT ENCOUNTERED AN ERROR OPENING AN HDF5 FILE'
            WRITE(6, *) '  FILENAME:  ', TRIM(string_val)
            WRITE(6, *) '  IERR:      ', ierr
        ELSEIF (error_num .eq. HDF5_WRITE_ERR) THEN
            WRITE(6, *) '  THRIFT ENCOUNTERED AN ERROR WRITING TO AN HDF5 FILE'
            WRITE(6, *) '  VARNAME:   ', TRIM(string_val)
            WRITE(6, *) '  IERR:      ', ierr
        ELSEIF (error_num .eq. HDF5_READ_ERR) THEN
            WRITE(6, *) '  THRIFT ENCOUNTERED AN ERROR READING FROM AN HDF5 FILE'
            WRITE(6, *) '  VARNAME:   ', TRIM(string_val)
            WRITE(6, *) '  IERR:      ', ierr
        ELSEIF (error_num .eq. HDF5_CLOSE_ERR) THEN
            WRITE(6, *) '  THRIFT ENCOUNTERED AN ERROR CLOSING AN HDF5 FILE'
            WRITE(6, *) '  FILENAME:  ', TRIM(string_val)
            WRITE(6, *) '  IERR:      ', ierr
        ELSEIF (error_num .eq. NAMELIST_READ_ERR) THEN
            WRITE(6, *) '  THRIFT ENCOUNTERED AN ERROR READING A NAMELIST'
            WRITE(6, *) '  ', TRIM(string_val)
            WRITE(6, *) '  IERR:      ', ierr
        ELSEIF (error_num .eq. D02CJF_ERR) THEN
            WRITE(6, *) '  THRIFT ENCOUNTERED A NAG ERROR (D02CJF)'
            WRITE(6, *) '     CALLING FUNCTION ', TRIM(string_val)
            WRITE(6, *) '     IERR:      ', ierr
            IF (ierr .eq. 1) THEN
                WRITE(6, *) '     CHECK INPUT PARAMETERS!'
            ELSEIF (ierr .eq. 2) THEN
                WRITE(6, *) '     VALUE OF TOL PREVENTS INTEGRATION!'
            ELSEIF (ierr .eq. 3) THEN
                WRITE(6, *) '     TOL TOO SMALL FOR FIRST STEP!'
            ELSEIF (ierr .eq. 4) THEN
                WRITE(6, *) '     XSOL HAS NOT BEEN RESET OR IS BEHIND X!'
            ELSEIF (ierr .eq. 5) THEN
                WRITE(6, *) '     XSOL HAS NOT BEEN RESET OR IS BEHIND X!'
            ELSEIF (ierr .eq. 6) THEN
                WRITE(6, *) '     ROOT COULD NOT BE FOUND!'
            ELSEIF (ierr .eq. 7) THEN
                WRITE(6, *) '     A SERIOUS NAG ERROR HAS OCCURED!'
            ELSE
                WRITE(6, *) '     UNKNOWN NAG ERROR!'
            END IF
        ELSEIF (error_num .eq. D05ADF_ERR) THEN
            WRITE(6, *) '  THRIFT ENCOUNTERED A NAG ERROR (D05ADF)'
            WRITE(6, *) '     CALLING FUNCTION ', TRIM(string_val)
            WRITE(6, *) '     IERR:      ', ierr
            IF (ierr .eq. 1) THEN
                WRITE(6, *) '     EPS <= 0 or A=B or F(A)xF(B)>0'
            ELSEIF (ierr .eq. 2) THEN
                WRITE(6, *) '     EPS IS TOO SMALL!'
            ELSEIF (ierr .eq. 3) THEN
                WRITE(6, *) '     POLE OF F FOUND!'
            ELSEIF (ierr .eq. 4) THEN
                WRITE(6, *) '     A SERIOUS NAG ERROR HAS OCCURED!'
            ELSE
                WRITE(6, *) '     UNKNOWN NAG ERROR!'
            END IF
        ELSEIF (error_num .eq. C05AJF_ERR) THEN
            WRITE(6, *) '  THRIFT ENCOUNTERED A NAG ERROR (C05AJF)'
            WRITE(6, *) '     CALLING FUNCTION ', TRIM(string_val)
            WRITE(6, *) '     IERR:      ', ierr
            IF (ierr .eq. 1) THEN
                WRITE(6, *) '     EPS <= 0 or NFMAX <= 0'
            ELSEIF (ierr .eq. 2) THEN
                WRITE(6, *) '     BAD SCALE FACTOR!'
            ELSEIF (ierr .eq. 3) THEN
                WRITE(6, *) '     NO ZERO OR ACCURACY TOO HIGH!'
            ELSEIF (ierr .eq. 4) THEN
                WRITE(6, *) '     MORE THAN NFMAX CALLS HAVE BEEN MADE!'
            ELSEIF (ierr .eq. 5) THEN
                WRITE(6, *) '     A SERIOUS NAG ERROR HAS OCCURED!'
            ELSE
                WRITE(6, *) '     UNKNOWN NAG ERROR!'
            END IF
        ELSEIF (error_num .eq. LSODE_ERR) THEN
            WRITE(6, *) '  THRIFT ENCOUNTERED AN LSODE ERROR'
            WRITE(6, *) '     CALLING FUNCTION ', TRIM(string_val)
            WRITE(6, *) '     IERR:      ', ierr
            IF (ierr .eq. - 1) THEN
                WRITE(6, *) '     EXCESS WORK DONE (check mf)!'
            ELSEIF (ierr .eq. - 2) THEN
                WRITE(6, *) '     TOLLERANCE TOO SMALL!'
            ELSEIF (ierr .eq. - 3) THEN
                WRITE(6, *) '     ILLEGAL INPUT DETECTED!'
            ELSEIF (ierr .eq. - 4) THEN
                WRITE(6, *) '     REPEATED ERROR TEST FAILURES!'
            ELSEIF (ierr .eq. - 5) THEN
                WRITE(6, *) '     REPEATED CONVERGENCE TEST FAILURES!'
            ELSEIF (ierr .eq. - 6) THEN
                WRITE(6, *) '     ERROR WEIGHT BECAME ZERO!'
            ELSE
                WRITE(6, *) '     UNKNOWN LSODE ERROR!'
            END IF
        ELSEIF (error_num .eq. RKH68_ERR) THEN
            WRITE(6, *) '  THRIFT ENCOUNTERED A RKH68 ERROR'
            WRITE(6, *) '  ROUTINE:   ', TRIM(string_val)
            WRITE(6, *) '  IERR:      ', ierr
        ELSEIF (error_num .eq. NAG_ERR) THEN
            WRITE(6, *) '  YOUR MACHINE DOES NOT SUPPORT NAG'
            WRITE(6, *) '  ROUTINE:   ', TRIM(string_val)
            WRITE(6, *) '  IERR:      ', ierr
        ELSEIF (error_num .eq. HDF5_ERR) THEN
            WRITE(6, *) '  THRIFT ENCOUNTERED AN HDF ERROR'
            WRITE(6, *) '  ROUTINE:   ', TRIM(string_val)
            WRITE(6, *) '  IERR:      ', ierr
        ELSEIF (error_num .eq. BOOZER_ERR) THEN
            WRITE(6, *) '  ERROR LOADING BOOZER FILE'
            WRITE(6, *) '  ROUTINE:   ', TRIM(string_val)
            WRITE(6, *) '  IERR:      ', ierr
        ELSEIF (error_num .eq. MPI_REDU_ERR) THEN
            WRITE(6, *) '  MPI REDUCE ERROR'
            WRITE(6, *) '  ROUTINE:   ', TRIM(string_val)
            WRITE(6, *) '  IERR:      ', ierr
        ELSEIF (error_num .eq. MPI_FINE_ERR) THEN
            WRITE(6, *) '  MPI FINALIZE ERROR'
            WRITE(6, *) '  ROUTINE:   ', TRIM(string_val)
            WRITE(6, *) '  IERR:      ', ierr
            CALL FLUSH(6)
            PRINT *,myworkid,' calling STOP'
            CALL FLUSH(6)
            STOP
        ELSEIF (error_num .eq. MPI_CHECK) THEN
        ELSE
            WRITE(6, *) '  THRIFT ENCOUNTERED AN UNKNOWN ERROR'
            WRITE(6, *) '  STRING: ', TRIM(string_val)
            WRITE(6, *) '  ierr:   ', ierr
        END IF
        CALL FLUSH(6)
#if defined(MPI_OPT)
!        CALL MPI_BARRIER(MPI_COMM_THRIFT,ierr_mpi)
!        ALLOCATE(error_array(1:nprocs_beams))
!        ierr_mpi = 0
!        CALL MPI_ALLGATHER(error_num,1,MPI_INTEGER,error_array,1,MPI_INTEGER,MPI_COMM_THRIFT,ierr_mpi)
!        ierr_mpi = 0
!        CALL MPI_BARRIER(MPI_COMM_THRIFT,ierr_mpi)
!        IF (ANY(error_array .ne. 0)) CALL MPI_FINALIZE(ierr_mpi)
!        DEALLOCATE(error_array)
        RETURN
#else
        IF (error_num .eq. MPI_CHECK) RETURN
#endif
        PRINT *,myworkid,' calling STOP'
        CALL FLUSH(6)
        STOP
    END SUBROUTINE handle_err

END MODULE thrift_runtime
