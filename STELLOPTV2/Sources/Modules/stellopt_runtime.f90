!-----------------------------------------------------------------------
!     Module:        stellopt_runtime
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/24/2012
!     Description:   This module contains various runtime parameters
!                    for the STELLOPT code.  It also contains the
!                    handle_error subroutine.
!     Version Notes:
!         2.11     - First stable version.
!         2.12     - Added support for calculation of the bootstrap
!                    current density separately from the total current.
!                  - Now supports all profile functions in VMEC.
!                  - Added functionality to control which current groups
!                    are used to calculate vac_mse signal.
!         2.20     - Robust version of LMDIF implemented
!                  - Profile functions normalization error in fcn fixed.
!                  - Added support for line integrated electron density.
!                  - Added support for Faraday rotation.
!         2.21     - Parameterized Beam and Bootstrap current profiles.
!                  - Added RBtor target
!                  - Added R0 target
!                  - Temporary input file written for robustness.
!                  - Glitches in robust LMDIF corrected
!         2.22     - NE_AUX_F and NE_OPT now normalized to (1E18 [m^-3])
!         2.23     - Addressed bug in get_equil_B(cyl)
!                  - MSE now compares total mse signals
!         2.24     - Faraday Rotation Diagnostic Implemented
!                  - Line integrals modified to handle C-shape plasma
!                  - Line integral bugs fixed (bad points now ignored)
!                  - Headers and version now printed in stellopt. file
!         2.25     - Added Refit Feature
!                  - Ne in stellopt. file output in [m^-3] units.
!                  - Added Z0 target
!                  - Addressed bug in CHISQ_SEPARATRIX
!         2.26     - Added Target Beta Poloidal
!                  - Added Target Beta Toroidal
!                  - Added Target Curvature (kertosis)
!                  - Added Target Resonant. Jacobian
!                  - Fixed Refit so it does restarts
!         2.27     - GIST incorporated into stellopt_txport
!                  - refit_param added to input namelist
!                  - Added error message regarding NE normalization
!                  - Addressed ERROR where fixed boundary calcs didn't
!                    report vmec error flag and would crash.
!                  - Adjusted unique_boundary and convert_boundary so
!                    exponent is passed via the rho_exp value below.
!         2.28     - rho_exp moved to input namelist (default = 2)
!                  - Fixed bug in chisq_iota for targeting s values.
!                  - Added auto_domain option
!                  - Added line integrated quantities to refit
!                  - Tested J_bootstrap refitting
!                  - Added npopulation to input namelist
!                  - Fixed glitch in line integrated ne refitting
!                  - Fixed bug in chisq_separatrix
!                  - Fixed Error with zeff_norm in stellopt_fcn
!         2.29     - Adjusted how Ne is normalized in chisq_ne
!                  - Made prof_to_vmec more robust
!                  - get_equil_XX modified to val intent(INOUT) from OUT
!                  - Modified how profiles are written to input files.
!                  - DIAG now set properly if mode = 2 in LMDIF optimization
!         2.30     - DKES Optimization now implemented
!         2.33     - AT/AH spline optimization now implemented
!                  - FLOW and ANIMEC now implemented
!                  - V_phi targeting now implemented
!                  - Support for Differential Evolution now implemented
!                  - Levenberg Routine Modified to get indexes on
!                    output files correct.
!                  - Stepopt routine now correctly outputs to xvec.dat
!                  - LAXIS_OPT option implemented for axis optimization
!                  - RBC/ZBS updating bug addressed
!         2.34     - Temporary backport to v8.47
!                  - chisq_ne normalization adjusted so if no targets are
!                    found (s<0.05) then s_min target is utilized.
!                  - Adjusted TXPORT screen output to include the word
!                    turbulence.
!                  - Fixed glitch with DKES where NU_DKES was defaulted
!                    to 0.0, now defaults to 0.01
!         2.35     - Now iterfaced to VMEC v8.51
!                  - Glitch identified in calculation of DKES at ns=2
!                    sigma defaulted to bigno to skip.
!         2.36     - Added nprobes to stellopt_targets so the number
!                    will be globally defined....makes modification
!                    easier...note diagno_runtime still only use 2048
!                  - Addressed bug where lfreeb = false and no boundary
!                    harmonics varied would cause code to crash.
!                  - Now compiles on NERSC (HOPPER) using CRAY
!                  - Addressed bug where txport_proxy not written to
!                    input files correctly.
!         2.40     - Added support for GENE linear (serial) computation
!                    of turbulent transport.
!                  - Added support for parallel codes.
!                  - Added support for GENE lienar (parallel)
!                    computation of turbluent transport.
!                  - Added TEM proxy prox_tem_proll
!         2.41     - Addressed bug in STELLOPT_TXPORT where GENE
!                    was only being run for one flux tube.
!                  - Various syntax fixes in LIBSTELL for support
!                    of GNU Fortran compilers on OSX
!                  - Now possible to compile STELLOPTV2 on OSX
!         2.42     - Addressed bug in TEM Proxy
!                  - Compiles with PGI compilers at PPPL and on Hopper
!         2.43     - Added ability to do hyperspace plane maps
!                  - RBC/ZBS_MIN/MAX variables now defined
!                  - stellopt_init.f90 reworked
!                  - Fixed glitch where files were moved instead of
!                    being copied which caused DIAGNO targets to crash.
!         2.44     - Modified txport calcualtion from
!                      kp1 = L2 - dpdx/2
!                          to
!                      kp1 = L2 - dpdx/2/Bhat
!                  - Fixed bugs in how magnetic diagnostics were handled
!                  - Addressed bugs in init regarding normalizations
!                    and min/max values when triangulate is used.
!         2.45     - Fixed minus sign glitch in tem_overlap
!                  - Fixed -m to m glitch in deltamn array size.
!         2.46     - Fixed glitch in DE2_interface.
!                  - Added 'one_iter' optimizaiton type
!                  - TEM proxy 'tem_bounce_tau' added
!                  - Energetic Particle Optimization introduced
!         2.47     - Added 'eval_xvec' optimization type
!                  - Switch from equil_utils to stel_tools
!                  - Corrected issues with DIAGNO and added ability
!                    to skip diagnostics.
!         2.48     - Added Kink Stability Target with TERPSICHORE
!                  - Added ECE Reflectrometry Target with TRAVIS
!                  - Removed old proxies from Turblent Transport
!                  - Modified default behavior of Zeff
!         2.49     - Upgraded interface to TERPSICHORE
!                  - Debugged interface to BEAMS3D for PPPL/Edison
!                  - Adjusted deffinition of alpha in TXPORT
!                  - Introduced Quasi-linear GENE transport runs
!                  - Added Vacuum Iota Targeting
!                  - PARVMEC Interfaced
!                  - Added NOPTIMIZERS option
!         2.50     - PARALLEL VMEC Working (no restart yet)
!                  - PARALLEL Boozer Working
!                  - PARALLEL BOOTSJ Working
!                  - Saving of Boozer files enabled
!                  - Helicity target weighting corrected.
!                  - Non-stellopt parameters properly handled by BOOTSJ
!                  - VBOOT Added as Equilibrium Option (VMEC+BOOTSJ)
!         2.51     - Elongation Targeting corrected (Kappa)
!-----------------------------------------------------------------------
      MODULE stellopt_runtime
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
      INTEGER, PARAMETER ::  CWS_READ_ERR      = 14
      INTEGER, PARAMETER ::  BAD_CWS_ERR       = 15
      INTEGER, PARAMETER ::  KNOT_MISMATCH_ERR = 16
      INTEGER, PARAMETER ::  KNOT_DEF_ERR      = 17
      INTEGER, PARAMETER ::  KNOT_ORDER_ERR    = 18
      INTEGER, PARAMETER ::  KNOT_CONST_ERR    = 19
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
      
      LOGICAL                  :: lverb, lkeep_mins, lneed_output, lrestart,&
                                  lrefit, lno_restart, lauto_domain, lparallel,&
                                  ltriangulate, lcoil_geom
      INTEGER                  :: nvars, mtargets, iter, mode, iunit_out,&
                                  cr_strategy, rho_exp, npopulation, noptimizers,&
                                  ier_paraexe
      INTEGER, ALLOCATABLE     :: var_dex(:),target_dex(:)
      INTEGER, ALLOCATABLE     :: arr_dex(:,:)
      REAL(rprec)              :: pi, pi2, mu0, ftol, xtol, gtol, epsfcn,&
                                  factor, chisq_min, refit_param, pct_domain
      REAL(rprec), ALLOCATABLE :: vars(:),targets(:),sigmas(:),vals(:),&
                                  diag(:),vars_min(:),vars_max(:)
      CHARACTER(256)           :: id_string, opt_type, proc_string, &
                                  proc_string_old, screen_str, xvec_file
      LOGICAL :: lcentered_differences
      
      REAL(rprec), PARAMETER :: STELLOPT_VERSION = 2.70
      
      REAL(rprec), PARAMETER :: bigno = 1.0E+10
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
            WRITE(6,*) '  SUBROUTINE: ',TRIM(string_val)
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
      ELSEIF (error_num .eq. CWS_READ_ERR) THEN
            WRITE(6,*) '  STELLOPT ENCOUNTERED AN ERROR READING A WINDING SURFACE'
            WRITE(6,*) '  FILENAME: ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num .eq. BAD_CWS_ERR) THEN
            WRITE(6,*) '  STELLOPT ENCOUNTERED A WINDING SURFACE ERROR'
            WRITE(6,*) '  ROUTINE:   ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num .eq. KNOT_MISMATCH_ERR) THEN
            WRITE(6,*) '  STELLOPT ENCOUNTERED A COIL COORDINATE SPLINE KNOT COUNT MISMATCH'
            WRITE(6,*) '  ROUTINE:   ',TRIM(string_val)
            WRITE(6,*) '  CTRL PT COUNT:      ',ierr
      ELSEIF (error_num .eq. KNOT_DEF_ERR) THEN
            WRITE(6,*) '  STELLOPT ENCOUNTERED A COIL SPLINE WITH LESS THAN FOUR KNOTS'
            WRITE(6,*) '  ROUTINE:   ',TRIM(string_val)
            WRITE(6,*) '  KNOT COUNT:      ',ierr
      ELSEIF (error_num .eq. KNOT_ORDER_ERR) THEN
            WRITE(6,*) '  STELLOPT ENCOUNTERED A COIL SPLINE WITH DESCENDING KNOTS'
            WRITE(6,*) '  ROUTINE:   ',TRIM(string_val)
            WRITE(6,*) '  KNOT:      ',ierr
      ELSEIF (error_num .eq. KNOT_CONST_ERR) THEN
            WRITE(6,*) '  FIRST AND LAST FOUR COIL SPLINE KNOTS MUST BE IDENTICAL'
            WRITE(6,*) '  ROUTINE:   ',TRIM(string_val)
            WRITE(6,*) '  COIL:      ',ierr
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
      
      END MODULE stellopt_runtime
