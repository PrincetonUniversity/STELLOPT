!-----------------------------------------------------------------------
!     Module:        torlines_runtime
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/2/2011
!     Description:   This module contains various runtime parameters
!                    for the TORLINES code.  It also contains the
!                    handle_error subroutine.
!     v1.0  6/2/15   - Parallelization and boundary extrapolation
!                    - emc3 routines in testing
!                    - Vessel targeting buggy and not yet working
!     v1.1  6/18/15  - Adjusted array references in
!                      torlines_init_external for GCC compliance
!                    - Fixed various GCC related bugs
!                    - MGRID bugs fixed
!     v1.12  6/18/15 - open changed to safe_open
!                    - Addressed extcur bug
!     v1.20  7/15/15 - Corrected vessel limiting
!                    - Fixed glitch causing VC to not be used when
!                      compiling with GCC.
!                    - Fixed glitch in using LSODE for backwards
!                      fieldline tracing.
!     v1.21  7/16/15 - Fixed glitch in vacuum metric elements
!-----------------------------------------------------------------------
      MODULE torlines_runtime
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds
      USE EZspline
!-----------------------------------------------------------------------
!     Module Variables
!          lverb        Logical to control screen output
!          ldone        Logical to control code completion
!          lrestart     Logical to control restarting from netCDF
!          lwout        Logical to control starting from a VMEC wout file.
!          lasym        Logical to conrol non-up/down symmetry
!          lfreeb       Logical to free boundary run
!          iter         Iteration number
!          nextcur      Number of external current systems
!          extsurfs     Number of extrapolated surfaces beyond VMEC surface
!          lastsurf     Last closed flux surface
!          vsurf        Location of VMEC surface in PIES coordinates
!          isurf        Dummy index used for passing current surface to subroutines
!          id_string    Filename of input file
!          mgrid_file   Filename of background field file
!          extcur       External currents for MGRID calculation
!----------------------------------------------------------------------
      IMPLICIT NONE
      
      
      INTEGER, PARAMETER ::  FILE_OPEN_ERR     = 1
      INTEGER, PARAMETER ::  ALLOC_ERR         = 11
      INTEGER, PARAMETER ::  NAMELIST_READ_ERR = 12
      INTEGER, PARAMETER ::  VMEC_INPUT_ERR    = 2
      INTEGER, PARAMETER ::  VMEC_WOUT_ERR     = 21
      INTEGER, PARAMETER ::  MGRID_ERR         = 22
      INTEGER, PARAMETER ::  EZSPLINE_ERR      = 4
      INTEGER, PARAMETER ::  NETCDF_OPEN_ERR      = 5
      INTEGER, PARAMETER ::  PIES_NETCDF_READ_ERR = 51
      INTEGER, PARAMETER ::  HDF5_ERR       = 6
      INTEGER, PARAMETER ::  HDF5_OPEN_ERR  = 61
      INTEGER, PARAMETER ::  HDF5_WRITE_ERR = 62
      INTEGER, PARAMETER ::  HDF5_CLOSE_ERR = 69
      INTEGER, PARAMETER ::  D02CJF_ERR = 71
      INTEGER, PARAMETER ::  D05ADF_ERR = 72
      INTEGER, PARAMETER ::  C05AJF_ERR = 73
      INTEGER, PARAMETER ::  E04CCF_ERR = 75
      INTEGER, PARAMETER ::  LSODE_ERR  = 791
      INTEGER, PARAMETER ::  NAG_ERR  = 792
      INTEGER, PARAMETER ::  MPI_ERR            = 80
      INTEGER, PARAMETER ::  MPI_INIT_ERR       = 801
      INTEGER, PARAMETER ::  MPI_RANK_ERR       = 802
      INTEGER, PARAMETER ::  MPI_SIZE_ERR       = 803
      INTEGER, PARAMETER ::  MPI_DUP_ERR        = 804
      INTEGER, PARAMETER ::  MPI_BARRIER_ERR    = 81
      INTEGER, PARAMETER ::  MPI_SEND_ERR       = 821
      INTEGER, PARAMETER ::  MPI_RECV_ERR       = 822
      INTEGER, PARAMETER ::  MPI_BCAST_ERR      = 83
      INTEGER, PARAMETER ::  MPI_FINE_ERR       = 89
      
      INTEGER, PARAMETER ::  MAXLINES   = 65536
      
      
      LOGICAL         :: lverb, ldone, lrestart, lasym, lfreeb, &
                         lcoil, lmgrid, lvmec, lpies, lspec, lvac,&
                         lvessel, lraw, lemc3, lauto, lvc_field
      INTEGER         :: nextcur, npoinc, nruntype, num_hcp, vsurf,&
                         nexternal, nprocs_torlines
      REAL(rprec)     :: mu, dphi, follow_tol, pi, pi2, mu0, &
                         delta_hc, eps1, eps2, eps3
      REAL(rprec), DIMENSION(MAXLINES)     :: r_start, phi_start, &
                                              z_start, phi_end, &
                                              r_hc, z_hc, phi_hc
      REAL(rprec), ALLOCATABLE :: extcur(:)
      CHARACTER(256)  :: id_string, mgrid_string, coil_string, &
                         vessel_string, int_type, restart_string
      
      REAL(rprec), PARAMETER :: TORLINES_VERSION = 2.0
!-----------------------------------------------------------------------
!     Subroutines
!          handle_error  Controls Program Termination
!-----------------------------------------------------------------------
      CONTAINS
      
      SUBROUTINE handle_err(error_num,string_val,ierr)
      IMPLICIT NONE
      INTEGER,INTENT(in)      :: error_num
      INTEGER,INTENT(in)      :: ierr
      CHARACTER(*),INTENT(in) :: string_val
      WRITE(6,*) '!!!!! ERROR !!!!!'
      
      IF (error_num == FILE_OPEN_ERR) THEN
            WRITE(6,*) '  TORLINES COULD NOT OPEN A FILE.'
            WRITE(6,*) '  FILENAME: ',TRIM(string_val)
            WRITE(6,*) '  IERR:     ',ierr
      ELSEIF (error_num == VMEC_INPUT_ERR) THEN
            WRITE(6,*) '  TORLINES COULD NOT READ THE VMEC INDATA NAMELIST'
            WRITE(6,*) '  FILENAME: ',TRIM(string_val)
            WRITE(6,*) '  IERR:     ',ierr
      ELSEIF (error_num == VMEC_WOUT_ERR) THEN
            WRITE(6,*) '  TORLINES COULD NOT READ THE VMEC WOUT FILE'
            WRITE(6,*) '  FILENAME: ',TRIM(string_val)
            WRITE(6,*) '  IERR:     ',ierr
      ELSEIF (error_num == ALLOC_ERR) THEN
            WRITE(6,*) '  TORLINES ENCOUNTERED AN ALLOCATION ERROR'
            WRITE(6,*) '  VARIABLES: ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num == EZSPLINE_ERR) THEN
            WRITE(6,*) '  TORLINES ENCOUNTERED AN EZSPLINE ERROR'
            WRITE(6,*) '  ROUTINE/VAR: ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
            CALL EZspline_error(ierr)
      ELSEIF (error_num == MGRID_ERR) THEN
            WRITE(6,*) '  TORLINES ENCOUNTERED AN ERROR READING THE MGRID FILE'
            WRITE(6,*) '  VARIABLES: ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num == NETCDF_OPEN_ERR) THEN
            WRITE(6,*) '  TORLINES ENCOUNTERED AN ERROR OPENING A NETCDF FILE'
            WRITE(6,*) '  VARIABLES: ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num .eq. HDF5_OPEN_ERR) THEN
            WRITE(6,*) '  FIELDLINES ENCOUNTERED AN ERROR OPENING AN HDF5 FILE'
            WRITE(6,*) '  FILENAME:  ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num .eq. HDF5_WRITE_ERR) THEN
            WRITE(6,*) '  FIELDLINES ENCOUNTERED AN ERROR WRITING TO AN HDF5 FILE'
            WRITE(6,*) '  VARNAME:   ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num .eq. HDF5_CLOSE_ERR) THEN
            WRITE(6,*) '  FIELDLINES ENCOUNTERED AN ERROR CLOSING AN HDF5 FILE'
            WRITE(6,*) '  FILENAME:  ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num == NAMELIST_READ_ERR) THEN
            WRITE(6,*) '  TORLINES ENCOUNTERED AN ERROR READING A NAMELIST'
            WRITE(6,*) '  ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num == D02CJF_ERR) THEN
            WRITE(6,*) '  TORLINES ENCOUNTERED A NAG ERROR (D02CJF)'
            WRITE(6,*) '     CALLING FUNCTION ',TRIM(string_val)
            WRITE(6,*) '     IERR:      ',ierr
            IF (ierr == 1) THEN
               WRITE(6,*) '     CHECK INPUT PARAMETERS!'
            ELSEIF (ierr == 2) THEN
               WRITE(6,*) '     VALUE OF TOL PREVENTS INTEGRATION!'
            ELSEIF (ierr == 3) THEN
               WRITE(6,*) '     TOL TOO SMALL FOR FIRST STEP!'
            ELSEIF (ierr == 4) THEN
               WRITE(6,*) '     XSOL HAS NOT BEEN RESET OR IS BEHIND X!'
            ELSEIF (ierr == 5) THEN
               WRITE(6,*) '     XSOL HAS NOT BEEN RESET OR IS BEHIND X!'
            ELSEIF (ierr == 6) THEN
               WRITE(6,*) '     ROOT COULD NOT BE FOUND!'
            ELSEIF (ierr == 7) THEN
               WRITE(6,*) '     A SERIOUS NAG ERROR HAS OCCURED!'
            ELSE
               WRITE(6,*) '     UNKNOWN NAG ERROR!'
            END IF
      ELSEIF (error_num == D05ADF_ERR) THEN
            WRITE(6,*) '  TORLINES ENCOUNTERED A NAG ERROR (D05ADF)'
            WRITE(6,*) '     CALLING FUNCTION ',TRIM(string_val)
            WRITE(6,*) '     IERR:      ',ierr
            IF (ierr == 1) THEN
               WRITE(6,*) '     EPS <= 0 or A=B or F(A)xF(B)>0'
            ELSEIF (ierr == 2) THEN
               WRITE(6,*) '     EPS IS TOO SMALL!'
            ELSEIF (ierr == 3) THEN
               WRITE(6,*) '     POLE OF F FOUND!'
            ELSEIF (ierr == 4) THEN
               WRITE(6,*) '     A SERIOUS NAG ERROR HAS OCCURED!'
            ELSE
               WRITE(6,*) '     UNKNOWN NAG ERROR!'
            END IF
      ELSEIF (error_num == C05AJF_ERR) THEN
            WRITE(6,*) '  TORLINES ENCOUNTERED A NAG ERROR (C05AJF)'
            WRITE(6,*) '     CALLING FUNCTION ',TRIM(string_val)
            WRITE(6,*) '     IERR:      ',ierr
            IF (ierr == 1) THEN
               WRITE(6,*) '     EPS <= 0 or NFMAX <= 0'
            ELSEIF (ierr == 2) THEN
               WRITE(6,*) '     BAD SCALE FACTOR!'
            ELSEIF (ierr == 3) THEN
               WRITE(6,*) '     NO ZERO OR ACCURACY TOO HIGH!'
            ELSEIF (ierr == 4) THEN
               WRITE(6,*) '     MORE THAN NFMAX CALLS HAVE BEEN MADE!'
            ELSEIF (ierr == 5) THEN
               WRITE(6,*) '     A SERIOUS NAG ERROR HAS OCCURED!'
            ELSE
               WRITE(6,*) '     UNKNOWN NAG ERROR!'
            END IF
      ELSEIF (error_num == E04CCF_ERR) THEN
            WRITE(6,*) '  TORLINES ENCOUNTERED A NAG ERROR (E04CCF)'
            WRITE(6,*) '     CALLING FUNCTION ',TRIM(string_val)
            WRITE(6,*) '     IERR:      ',ierr
            IF (ierr == 1) THEN
               WRITE(6,*) '     ON ENTRY, N < 1'
               WRITE(6,*) '     OR        TOL < MACHINE PRECISION'
               WRITE(6,*) '     OR        IW /= N+1'
               WRITE(6,*) '     OR        MAXCAL < 1'
            ELSEIF (ierr == 2) THEN
               WRITE(6,*) '     MAXCAL FUNCTION EVALUATIONS REACHED'
               RETURN
            END IF
      ELSEIF (error_num == HDF5_ERR) THEN
            WRITE(6,*) '  TORLINES ENCOUNTERED AN HDF ERROR'
            WRITE(6,*) '  ROUTINE:   ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num == NAG_ERR) THEN
            WRITE(6,*) '  NAG NOT INSTALLED ON THIS MACHINE'
            WRITE(6,*) '  SWITCH INT_TYPE TO LSODE AND RERUN'
            WRITE(6,*) '  ROUTINE:   ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSE
           WRITE(6,*) '  TORLINES ENCOUNTERED A GENERIC ERROR'
           WRITE(6,*) '  STRING: ',TRIM(string_val)
           WRITE(6,*) '  ierr:   ',ierr
      END IF
      CALL FLUSH(6)
      STOP
      END SUBROUTINE handle_err
      
      END MODULE torlines_runtime
