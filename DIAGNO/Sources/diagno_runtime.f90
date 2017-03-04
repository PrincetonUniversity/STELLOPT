!-----------------------------------------------------------------------
!     Module:        diagno_runtime
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/02/2012
!     Description:   This module contains various runtime parameters
!                    for the DIAGNO code.  It also contains the
!                    handle_error subroutine.
!
!     Notes:
!       v3.11        Properly handles curtor printout.
!       v3.12        Added BPROBE_TURNS to namelist
!       v3.13        Adjusted ordering of BPROBE output to give accurate
!                    values, adjusted form of rogowski output
!       v3.20        Fixed error in vacuum integrated fields
!       v3.25        Vaccum currents now selectable (luse_extcur)
!                    Output format now ES22.12E3.
!       v3.30        Added -bench benchmarking capability
!                    bprobe_turns increased to 2048 elements for ITER
!       v3.31        Fixed bug in EXTCUR array not reading all elements
!                    in diagno_init_coil.f90
!       v3.32        Fixed bug which caused lfreeb in VMEC to default
!                    to true.
!       v3.40        Modified Rogowski deffinition to handle W7-X style
!                    sensors with multiple unconnected traces in a
!                    given sensor.
!       v3.41        Added LSKIP_FLUX and LSKIP_ROGO to skip
!                    sensor calculations to speed up STELLOPT runs
!                    where not all sensors are used.
!                    Added NFP factor back in to VC
!                    Added EQ_SGNS to handle PHIEDGE sign issues
!-----------------------------------------------------------------------
      MODULE diagno_runtime
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE EZspline
!-----------------------------------------------------------------------
!     Module Variables
!          lverb          Logical to control screen output
!          lvmec          Logical to indicate VMEC Equilibria
!          lpies          Logical to indicate PIES Equilibria
!          lspec          Logical to indicate SPEC Equilibria
!          lcoil          Logical to indicate coil file
!          lvac           Logical to indicate vacuum fields only
!          lmut           Logical to control mutual induction calc
!          luse_mut       Logical to control use of mutal induction (not fully implemented)
!          lvc_field      Logical to control use of VIRTUAL Casing
!          nextcur        Number of external currents
!          int_step       Number of sub steps to take for diagnostic integration
!          nu             Number of polodial gridpoints
!          nv             Number of toroidal gridpoints
!          nfp_diagno     Number of field periods
!          flux_turns     Scalar integer multiple for flux loops
!          vc_adapt_tol   Absolute tollerance for adaptive integration
!          vc_adapt_rel   Relative tollerance for adaptive integration
!          units          Scale factor for units
!          phiedge        Total Enclosed toroidal flux
!          id_string      Equilibrium Filename
!          coil_string    Coil Filename
!          flux_diag_file Fluxloop specification file
!          bprobes_file   B-Field diagnostics specification file
!          mirnov_file    Mirnov Array specification file
!          seg_rog_file   Segmented Rogowski Coil specification file
!          int_type       Integration type ('midpoint','simpson','bode')
!          flux_mut_file  Mutual inductance file
!          bfield_points_file  B-Field at a point diagnostic
!          extcur         External currents for MGRID calculation
!----------------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER, PARAMETER ::  FILE_OPEN_ERR     = 1
      INTEGER, PARAMETER ::  ALLOC_ERR         = 11
      INTEGER, PARAMETER ::  NAMELIST_READ_ERR = 12
      INTEGER, PARAMETER ::  BAD_INPUT_ERR     = 13
      INTEGER, PARAMETER ::  VMEC_INPUT_ERR    = 2
      INTEGER, PARAMETER ::  VMEC_WOUT_ERR     = 21
      INTEGER, PARAMETER ::  MGRID_ERR         = 22
      INTEGER, PARAMETER ::  EZSPLINE_ERR      = 4
      INTEGER, PARAMETER ::  NETCDF_OPEN_ERR      = 5
      INTEGER, PARAMETER ::  PIES_NETCDF_READ_ERR = 51
      INTEGER, PARAMETER ::  HDF5_ERR       = 6
      INTEGER, PARAMETER ::  D02CJF_ERR = 71
      INTEGER, PARAMETER ::  D05ADF_ERR = 72
      INTEGER, PARAMETER ::  C05AJF_ERR = 73
      INTEGER, PARAMETER ::  MAXLINES   = 256
      
      LOGICAL         :: lverb, lvmec, lpies, lspec, lcoil, lvac, &
                         lrphiz, lmut, luse_mut, lvc_field
      LOGICAL         :: luse_extcur(512),lskip_flux(512),lskip_rogo(512)
      INTEGER         :: nextcur, int_step, nu, nv, nfp_diagno, eq_sgns
      REAL(rprec)     :: flux_turns(512), bprobe_turns(2048)
      REAL(rprec)     :: vc_adapt_tol, units, phiedge, vc_adapt_rel
      REAL(rprec), ALLOCATABLE :: extcur(:)
      CHARACTER(256)  :: id_string, flux_diag_file, bprobes_file, &
                         mirnov_file, seg_rog_file, bfield_points_file,&
                         int_type, coil_string, flux_mut_file
                         
      REAL(rprec), PARAMETER :: DIAGNO_VERSION = 3.41
      REAL(rprec), PARAMETER ::      pi = 3.14159265358979312D+00
      REAL(rprec), PARAMETER ::     pi2 = 6.28318530717958623
      REAL(rprec), PARAMETER ::  onerad = 1.74532925199432955E-02
      REAL(rprec), PARAMETER ::     mu0 = 1.25663706143591729E-06
      REAL(rprec), PARAMETER :: sqrtmu0 = 1.12099824327958572E-03
!-----------------------------------------------------------------------
!     Subroutines
!          handle_error  Controls Program Termination
!-----------------------------------------------------------------------
      CONTAINS
      
      SUBROUTINE handle_error(error_num,string_val,ierr)
      IMPLICIT NONE
      INTEGER,INTENT(in)      :: error_num
      INTEGER,INTENT(inout)      :: ierr
      CHARACTER(*),INTENT(in) :: string_val
      WRITE(6,*) '!!!!! ERROR !!!!!'
      
      IF (error_num .eq. FILE_OPEN_ERR) THEN
            WRITE(6,*) '  DIAGNO COULD NOT OPEN A FILE.'
            WRITE(6,*) '  FILENAME: ',TRIM(string_val)
            WRITE(6,*) '  IERR:     ',ierr
      ELSEIF (error_num .eq. VMEC_INPUT_ERR) THEN
            WRITE(6,*) '  DIAGNO COULD READ THE VMEC INDATA NAMELIST'
            WRITE(6,*) '  FILENAME: ',TRIM(string_val)
            WRITE(6,*) '  IERR:     ',ierr
      ELSEIF (error_num .eq. VMEC_WOUT_ERR) THEN
            WRITE(6,*) '  DIAGNO COULD READ THE VMEC WOUT FILE'
            WRITE(6,*) '  FILENAME: ',TRIM(string_val)
            WRITE(6,*) '  IERR:     ',ierr
      ELSEIF (error_num .eq. ALLOC_ERR) THEN
            WRITE(6,*) '  DIAGNO ENCOUNTERED AN ALLOCATION ERROR'
            WRITE(6,*) '  VARIABLES: ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num .eq. EZSPLINE_ERR) THEN
            WRITE(6,*) '  DIAGNO ENCOUNTERED AN EZSPLINE ERROR'
            WRITE(6,*) '  ROUTINE/VAR: ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
            CALL EZspline_error(ierr)
      ELSEIF (error_num .eq. MGRID_ERR) THEN
            WRITE(6,*) '  DIAGNO ENCOUNTERED AN ERROR READING THE MGRID FILE'
            WRITE(6,*) '  VARIABLES: ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num .eq. NETCDF_OPEN_ERR) THEN
            WRITE(6,*) '  DIAGNO ENCOUNTERED AN ERROR OPENING A NETCDF FILE'
            WRITE(6,*) '  VARIABLES: ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num .eq. NAMELIST_READ_ERR) THEN
            WRITE(6,*) '  DIAGNO ENCOUNTERED AN ERROR READING A NAMELIST'
            WRITE(6,*) '  ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num .eq. D02CJF_ERR) THEN
            WRITE(6,*) '  DIAGNO ENCOUNTERED A NAG ERROR (D02CJF)'
            WRITE(6,*) '     CALLING FUNCTION',TRIM(string_val)
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
            WRITE(6,*) '  DIAGNO ENCOUNTERED A NAG ERROR (D05ADF)'
            WRITE(6,*) '     CALLING FUNCTION',TRIM(string_val)
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
            WRITE(6,*) '  DIAGNO ENCOUNTERED A NAG ERROR (C05AJF)'
            WRITE(6,*) '     CALLING FUNCTION',TRIM(string_val)
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
      ELSEIF (error_num .eq. HDF5_ERR) THEN
            WRITE(6,*) '  DIAGNO ENCOUNTERED AN HDF ERROR'
            WRITE(6,*) '  ROUTINE:   ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSE
           WRITE(6,*) '  DIAGNO ENCOUNTERED A GENERIC ERROR'
           WRITE(6,*) '  STRING: ',TRIM(string_val)
           WRITE(6,*) '  ierr:   ',ierr
      END IF
      CALL FLUSH(6)
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_FINALIZE(ierr)   
!DEC$ ENDIF
      STOP
      END SUBROUTINE handle_error
      
      END MODULE diagno_runtime
