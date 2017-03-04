!-----------------------------------------------------------------------
!     Module:        stelltran_runtime
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          08/21/2015
!     Description:   This module contains various runtime parameters
!                    for the STELLTRAN code.  It also contains the
!                    handle_error subroutine.
!     Version Notes:
!         1.00     - Initial Release 
!-----------------------------------------------------------------------
      MODULE stelltran_runtime
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE EZspline
!-----------------------------------------------------------------------
!     Module Variables
!          lverb         Logical to control screen output
!          pi            PI
!          pi2           2*PI
!          mu0           4*PI*10^-7
!          id_string     Equilibrium Filename
!----------------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER, PARAMETER ::  FILE_OPEN_ERR     = 1
      INTEGER, PARAMETER ::  FILE_EXIST_ERR    = 10
      INTEGER, PARAMETER ::  ALLOC_ERR         = 11
      INTEGER, PARAMETER ::  NAMELIST_READ_ERR = 12
      INTEGER, PARAMETER ::  BAD_INPUT_ERR     = 13
      INTEGER, PARAMETER ::  VMEC_INPUT_ERR    = 2
      INTEGER, PARAMETER ::  VMEC_WOUT_ERR     = 21
      INTEGER, PARAMETER ::  MGRID_ERR         = 22
      INTEGER, PARAMETER ::  VMEC_RUN_ERR      = 23
      INTEGER, PARAMETER ::  EZSPLINE_ERR           = 4
      INTEGER, PARAMETER ::  NETCDF_OPEN_ERR      = 5
      INTEGER, PARAMETER ::  PIES_NETCDF_READ_ERR = 51
      INTEGER, PARAMETER ::  HDF5_ERR       = 6
      INTEGER, PARAMETER ::  HDF5_OPEN_ERR  = 61
      INTEGER, PARAMETER ::  HDF5_WRITE_ERR = 62
      INTEGER, PARAMETER ::  HDF5_CLOSE_ERR = 69
      INTEGER, PARAMETER ::  NAG_ERR    = 70
      INTEGER, PARAMETER ::  D02CJF_ERR = 71
      INTEGER, PARAMETER ::  D05ADF_ERR = 72
      INTEGER, PARAMETER ::  C05AJF_ERR = 73
      INTEGER, PARAMETER ::  LSODE_ERR  = 791
      INTEGER, PARAMETER ::  RKH68_ERR  = 792
      INTEGER, PARAMETER ::  MPI_ERR            = 800
      INTEGER, PARAMETER ::  MPI_INIT_ERR       = 801
      INTEGER, PARAMETER ::  MPI_RANK_ERR       = 802
      INTEGER, PARAMETER ::  MPI_SIZE_ERR       = 803
      INTEGER, PARAMETER ::  MPI_BARRIER_ERR    = 810
      INTEGER, PARAMETER ::  MPI_SEND_ERR       = 821
      INTEGER, PARAMETER ::  MPI_RECV_ERR       = 822
      INTEGER, PARAMETER ::  MPI_BCAST_ERR      = 830
      INTEGER, PARAMETER ::  MPI_FINE_ERR       = 890

      REAL(rprec), PARAMETER :: STELLTRAN_VERSION = 0.01
      
      LOGICAL                  :: lverb, lrestart, lfirst_pass,&
                                  lveryfirst_pass
      REAL(rprec)              :: pi, pi2, mu0, fit_tol
      CHARACTER(256)           :: id_string, proc_string, &
                                  proc_string_old, db_file, &
                                  equil_type, run_type, vessel_ecrh,&
                                  mirror_ecrh, antennatype_ecrh,&
                                  targettype_ecrh
      
!-----------------------------------------------------------------------
!     Subroutines
!          handle_err  Controls Program Termination
!-----------------------------------------------------------------------
      CONTAINS
      
      SUBROUTINE handle_err(error_num,string_val,ierr)
      IMPLICIT NONE
      INTEGER,INTENT(in)      :: error_num
      INTEGER,INTENT(in)      :: ierr
      CHARACTER(*),INTENT(in) :: string_val
      INTEGER                 :: ierr2
   
      WRITE(6,*) '!!!!! ERROR !!!!!'
      
      IF (error_num .eq. FILE_OPEN_ERR) THEN
            WRITE(6,*) '  STELLOPT COULD NOT OPEN A FILE.'
            WRITE(6,*) '  FILENAME: ',TRIM(string_val)
            WRITE(6,*) '  IERR:     ',ierr
      ELSEIF (error_num .eq. FILE_EXIST_ERR) THEN
            WRITE(6,*) '  STELLOPT COULD NOT FIND A FILE'
            WRITE(6,*) '  FILENAME: ',TRIM(string_val)
      ELSEIF (error_num .eq. VMEC_INPUT_ERR) THEN
            WRITE(6,*) '  STELLOPT COULD READ THE VMEC INDATA NAMELIST'
            WRITE(6,*) '  FILENAME: ',TRIM(string_val)
            WRITE(6,*) '  IERR:     ',ierr
      ELSEIF (error_num .eq. VMEC_WOUT_ERR) THEN
            WRITE(6,*) '  STELLOPT COULD READ THE VMEC WOUT FILE'
            WRITE(6,*) '  FILENAME: ',TRIM(string_val)
            WRITE(6,*) '  IERR:     ',ierr
      ELSEIF (error_num .eq. VMEC_RUN_ERR) THEN
            WRITE(6,*) '  STELLOPT COULD NOT RUN VMEC'
            WRITE(6,*) '  MESSAGE:  ',TRIM(string_val)
            WRITE(6,*) '  IERR:     ',ierr
      ELSEIF (error_num .eq. ALLOC_ERR) THEN
            WRITE(6,*) '  STELLOPT ENCOUNTERED AN ALLOCATION ERROR'
            WRITE(6,*) '  VARIABLES: ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num .eq. EZSPLINE_ERR) THEN
            WRITE(6,*) '  STELLOPT ENCOUNTERED AN EZSPLINE ERROR'
            WRITE(6,*) '  ROUTINE/VAR: ',TRIM(string_val)
            WRITE(6,*) '  IERR:        ',ierr
            CALL EZspline_error(ierr)
      ELSEIF (error_num .eq. MGRID_ERR) THEN
            WRITE(6,*) '  STELLOPT ENCOUNTERED AN ERROR READING THE MGRID FILE'
            WRITE(6,*) '  VARIABLES: ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num .eq. NETCDF_OPEN_ERR) THEN
            WRITE(6,*) '  STELLOPT ENCOUNTERED AN ERROR OPENING A NETCDF FILE'
            WRITE(6,*) '  VARIABLES: ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num .eq. HDF5_OPEN_ERR) THEN
            WRITE(6,*) '  STELLOPT ENCOUNTERED AN ERROR OPENING AN HDF5 FILE'
            WRITE(6,*) '  FILENAME:  ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num .eq. HDF5_WRITE_ERR) THEN
            WRITE(6,*) '  STELLOPT ENCOUNTERED AN ERROR WRITING TO AN HDF5 FILE'
            WRITE(6,*) '  VARNAME:   ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num .eq. HDF5_CLOSE_ERR) THEN
            WRITE(6,*) '  STELLOPT ENCOUNTERED AN ERROR CLOSING AN HDF5 FILE'
            WRITE(6,*) '  FILENAME:  ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num .eq. NAMELIST_READ_ERR) THEN
            WRITE(6,*) '  STELLOPT ENCOUNTERED AN ERROR READING A NAMELIST'
            WRITE(6,*) '  ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num .eq. D02CJF_ERR) THEN
            WRITE(6,*) '  STELLOPT ENCOUNTERED A NAG ERROR (D02CJF)'
            WRITE(6,*) '     CALLING FUNCTION ',TRIM(string_val)
            WRITE(6,*) '     IERR:      ',ierr
            IF (ierr .eq. 1) THEN
               WRITE(6,*) '     CHECK INPUT PARAMETERS!'
            ELSEIF (ierr .eq. 2) THEN
               WRITE(6,*) '     VALUE OF TOL PREVENTS INTEGRATION!'
            ELSEIF (ierr .eq. 3) THEN
               WRITE(6,*) '     TOL TOO SMALL FOR FIRST STEP!'
            ELSEIF (ierr .eq. 4) THEN
               WRITE(6,*) '     XSOL HAS NOT BEEN RESET OR IS BEHIND X!'
            ELSEIF (ierr .eq. 5) THEN
               WRITE(6,*) '     XSOL HAS NOT BEEN RESET OR IS BEHIND X!'
            ELSEIF (ierr .eq. 6) THEN
               WRITE(6,*) '     ROOT COULD NOT BE FOUND!'
            ELSEIF (ierr .eq. 7) THEN
               WRITE(6,*) '     A SERIOUS NAG ERROR HAS OCCURED!'
            ELSE
               WRITE(6,*) '     UNKNOWN NAG ERROR!'
            END IF
      ELSEIF (error_num .eq. D05ADF_ERR) THEN
            WRITE(6,*) '  STELLOPT ENCOUNTERED A NAG ERROR (D05ADF)'
            WRITE(6,*) '     CALLING FUNCTION ',TRIM(string_val)
            WRITE(6,*) '     IERR:      ',ierr
            IF (ierr .eq. 1) THEN
               WRITE(6,*) '     EPS <= 0 or A=B or F(A)xF(B)>0'
            ELSEIF (ierr .eq. 2) THEN
               WRITE(6,*) '     EPS IS TOO SMALL!'
            ELSEIF (ierr .eq. 3) THEN
               WRITE(6,*) '     POLE OF F FOUND!'
            ELSEIF (ierr .eq. 4) THEN
               WRITE(6,*) '     A SERIOUS NAG ERROR HAS OCCURED!'
            ELSE
               WRITE(6,*) '     UNKNOWN NAG ERROR!'
            END IF
      ELSEIF (error_num .eq. C05AJF_ERR) THEN
            WRITE(6,*) '  STELLOPT ENCOUNTERED A NAG ERROR (C05AJF)'
            WRITE(6,*) '     CALLING FUNCTION ',TRIM(string_val)
            WRITE(6,*) '     IERR:      ',ierr
            IF (ierr .eq. 1) THEN
               WRITE(6,*) '     EPS <= 0 or NFMAX <= 0'
            ELSEIF (ierr .eq. 2) THEN
               WRITE(6,*) '     BAD SCALE FACTOR!'
            ELSEIF (ierr .eq. 3) THEN
               WRITE(6,*) '     NO ZERO OR ACCURACY TOO HIGH!'
            ELSEIF (ierr .eq. 4) THEN
               WRITE(6,*) '     MORE THAN NFMAX CALLS HAVE BEEN MADE!'
            ELSEIF (ierr .eq. 5) THEN
               WRITE(6,*) '     A SERIOUS NAG ERROR HAS OCCURED!'
            ELSE
               WRITE(6,*) '     UNKNOWN NAG ERROR!'
            END IF
      ELSEIF (error_num .eq. LSODE_ERR) THEN
            WRITE(6,*) '  STELLOPT ENCOUNTERED AN LSODE ERROR'
            WRITE(6,*) '     CALLING FUNCTION ',TRIM(string_val)
            WRITE(6,*) '     IERR:      ',ierr
            IF (ierr .eq. -1) THEN
               WRITE(6,*) '     EXCESS WORK DONE (check mf)!'
            ELSEIF (ierr .eq. -2) THEN
               WRITE(6,*) '     TOLLERANCE TOO SMALL!'
            ELSEIF (ierr .eq. -3) THEN
               WRITE(6,*) '     ILLEGAL INPUT DETECTED!'
            ELSEIF (ierr .eq. -4) THEN
               WRITE(6,*) '     REPEATED ERROR TEST FAILURES!'
            ELSEIF (ierr .eq. -5) THEN
               WRITE(6,*) '     REPEATED CONVERGENCE TEST FAILURES!'
            ELSEIF (ierr .eq. -6) THEN
               WRITE(6,*) '     ERROR WEIGHT BECAME ZERO!'
            ELSE
               WRITE(6,*) '     UNKNOWN LSODE ERROR!'
            END IF
      ELSEIF (error_num .eq. RKH68_ERR) THEN
            WRITE(6,*) '  STELLOPT ENCOUNTERED A RKH68 ERROR'
            WRITE(6,*) '  ROUTINE:   ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num .eq. NAG_ERR) THEN
            WRITE(6,*) '  YOUR MACHINE DOES NOT SUPPORT NAG'
            WRITE(6,*) '  ROUTINE:   ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num .eq. HDF5_ERR) THEN
            WRITE(6,*) '  STELLOPT ENCOUNTERED AN HDF ERROR'
            WRITE(6,*) '  ROUTINE:   ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num <= MPI_ERR) THEN
            WRITE(6,*) '  STELLOPT ENCOUNTERED AN MPI ERROR'
            WRITE(6,*) '  ROUTINE:   ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSE
           WRITE(6,*) '  STELLOPT ENCOUNTERED AN UNKNOWN ERROR'
           WRITE(6,*) '  STRING: ',TRIM(string_val)
           WRITE(6,*) '  ierr:   ',ierr
      END IF
      CALL FLUSH(6)
!DEC$ IF DEFINED (MPI_OPT)
      ierr2 = 0
      !CALL stellopt_paraexe('exit','exit',.false.)
      CALL MPI_FINALIZE(ierr2)   
!DEC$ ENDIF
      STOP
      END SUBROUTINE handle_err
      
      END MODULE stelltran_runtime
