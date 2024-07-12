!-----------------------------------------------------------------------
!     Module:        spec_runtime
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          04/06/2012
!     Description:   This module contains various runtime parameters
!                    for the VMEC2SPEC code.  It also contains the
!                    handle_error subroutine.
!-----------------------------------------------------------------------
      MODULE spec_runtime
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
!          lmake_coils  Logical to control creation of a coil_data file.
!          iter         Iteration number
!          nextcur      Number of external current systems
!          extsurfs     Number of extrapolated surfaces beyond VMEC surface
!          lastsurf     Last closed flux surface
!          vsurf        Location of VMEC surface in SPEC coordinates
!          m_new        SPEC input m
!          n_new        SPEC input n
!          free_override Overrides Fixed or Free Boundary selction in VMEC
!          id_string    Filename of input file
!          mgrid_file   Filename of background field file
!          extcur       External currents for MGRID calculation
!          coils_file   Coils File for conversion
!----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL         :: lverb, lwout, lasym, lfreeb, lflipped
      INTEGER         :: nextcur
      INTEGER         :: m_new, n_new
      CHARACTER(256)  :: id_string
      REAL(rprec), ALLOCATABLE :: extcur(:)
      
      INTEGER, PARAMETER ::  FILE_OPEN_ERR     = 1
      INTEGER, PARAMETER ::  ALLOC_ERR         = 11
      INTEGER, PARAMETER ::  NAMELIST_READ_ERR = 12
      INTEGER, PARAMETER ::  BAD_INPUT_ERR     = 13
      INTEGER, PARAMETER ::  VMEC_INPUT_ERR    = 2
      INTEGER, PARAMETER ::  VMEC_WOUT_ERR     = 21
      INTEGER, PARAMETER ::  MGRID_ERR         = 22
      INTEGER, PARAMETER ::  EZSPLINE_ERR      = 4
      INTEGER, PARAMETER ::  NETCDF_OPEN_ERR      = 5
      INTEGER, PARAMETER ::  SPEC_NETCDF_READ_ERR = 51
      INTEGER, PARAMETER ::  HDF5_ERR       = 6
      INTEGER, PARAMETER ::  D02CJF_ERR = 71
      INTEGER, PARAMETER ::  D05ADF_ERR = 72
      INTEGER, PARAMETER ::  C05AJF_ERR = 73
      
      REAL(rprec), PARAMETER :: VMEC2SPEC_VERSION = 1.5
      REAL(rprec), PARAMETER :: mu0 = 1.25663706144E-06
!-----------------------------------------------------------------------
!     Subroutines
!          handle_error  Controls Program Termination
!-----------------------------------------------------------------------
      CONTAINS
      
      SUBROUTINE handle_error(error_num,string_val,ierr)
      IMPLICIT NONE
      INTEGER,INTENT(in)      :: error_num
      INTEGER,INTENT(in)      :: ierr
      CHARACTER(*),INTENT(in) :: string_val
      WRITE(6,*) '!!!!! ERROR !!!!!'
      
      IF (error_num == FILE_OPEN_ERR) THEN
            WRITE(6,*) '  VMEC2SPEC COULD NOT OPEN A FILE.'
            WRITE(6,*) '  FILENAME: ',TRIM(string_val)
            WRITE(6,*) '  IERR:     ',ierr
      ELSEIF (error_num == VMEC_INPUT_ERR) THEN
            WRITE(6,*) '  VMEC2SPEC COULD READ THE VMEC INDATA NAMELIST'
            WRITE(6,*) '  FILENAME: ',TRIM(string_val)
            WRITE(6,*) '  IERR:     ',ierr
      ELSEIF (error_num == VMEC_WOUT_ERR) THEN
            WRITE(6,*) '  VMEC2SPEC COULD READ THE VMEC WOUT FILE'
            WRITE(6,*) '  FILENAME: ',TRIM(string_val)
            WRITE(6,*) '  IERR:     ',ierr
      ELSEIF (error_num == ALLOC_ERR) THEN
            WRITE(6,*) '  VMEC2SPEC ENCOUNTERED AN ALLOCATION ERROR'
            WRITE(6,*) '  VARIABLES: ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num == EZSPLINE_ERR) THEN
            WRITE(6,*) '  VMEC2SPEC ENCOUNTERED AN EZSPLINE ERROR'
            WRITE(6,*) '  ROUTINE/VAR: ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
            CALL EZspline_error(ierr)
      ELSEIF (error_num == MGRID_ERR) THEN
            WRITE(6,*) '  VMEC2SPEC ENCOUNTERED AN ERROR READING THE MGRID FILE'
            WRITE(6,*) '  VARIABLES: ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num == NETCDF_OPEN_ERR) THEN
            WRITE(6,*) '  VMEC2SPEC ENCOUNTERED AN ERROR OPENING A NETCDF FILE'
            WRITE(6,*) '  VARIABLES: ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num == NAMELIST_READ_ERR) THEN
            WRITE(6,*) '  VMEC2SPEC ENCOUNTERED AN ERROR READING A NAMELIST'
            WRITE(6,*) '  ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num == D02CJF_ERR) THEN
            WRITE(6,*) '  VMEC2SPEC ENCOUNTERED A NAG ERROR (D02CJF)'
            WRITE(6,*) '     CALLING FUNCTION',TRIM(string_val)
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
            WRITE(6,*) '  VMEC2SPEC ENCOUNTERED A NAG ERROR (D05ADF)'
            WRITE(6,*) '     CALLING FUNCTION',TRIM(string_val)
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
            WRITE(6,*) '  VMEC2SPEC ENCOUNTERED A NAG ERROR (C05AJF)'
            WRITE(6,*) '     CALLING FUNCTION',TRIM(string_val)
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
      ELSEIF (error_num == HDF5_ERR) THEN
            WRITE(6,*) '  VMEC2SPEC ENCOUNTERED AN HDF ERROR'
            WRITE(6,*) '  ROUTINE:   ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSE
           WRITE(6,*) '  VMEC2SPEC ENCOUNTERED A GENERIC ERROR'
           WRITE(6,*) '  STRING: ',TRIM(string_val)
           WRITE(6,*) '  ierr:   ',ierr
      END IF
      CALL FLUSH(6)
      STOP
      END SUBROUTINE handle_error
      
      END MODULE spec_runtime
