!-----------------------------------------------------------------------
!     Module:        fieldlines_runtime
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          02/21/2012
!     Description:   This module contains various runtime parameters
!                    for the FIELDLINES code.  It also contains the
!                    handle_error subroutine.
!     v1.1 12/20/12  - Added -full modifier for axis finding and grid
!                      refinement.
!     v1.2 01/03/12  - Added Homocline finder
!     v1.21 03/10/14 - Fixed bug in EXTCUR array not reading all elements
!                      in fieldlines_init_coil.f90
!     v1.21 03/10/14 - Cleaned up VC and added grid masking arrays
!     v1.22 05/02/14 - Added support for wall_mod over vessel_mod
!                    - Modified all tracing routines to use
!                      out_fieldlines_nag routine so we don't duplicate
!                      efforts.
!                    - Added diffusion operator to out_fieldlines_nag,
!                      removed old ones from fblin routines.
!     v1.25 07/31/14 - Fixed binary output 
!     v1.26 08/25/14 - Addressed LSODE initial step error.
!     v1.27 08/25/14 - Addressed memory leak when running in batch mode.
!     v1.28 12/18/14 - Added -raw flag for using scaled currents.
!     v1.29 01/23/15 - Added -emc3 flag for EMC3 compatible output.
!     v1.30 03/10/15 - Fixed RKH68 and LSODE routines to work properly.
!                    - Added support for LSODE and RKH68 in unstable
!                      manifold solver. (untested as of this release)
!     v1.31 03/18/15 - Reverse flag now properly works.
!     v1.32 04/18/16 - Fixed issue with ATOL in follow_single.f90
!                    - Added lwall_trans option
!                    - MU Should properly work
!                    - Fixed glitch in reverse option
!     v1.40 08/30/16 - Switch away from que type data transfer to 
!                      gatherv statements.
!                    - ModB along fieldline now saved.
!                    - MU benchmarked by B. Isralei
!                    - Wall strike points saved.
!     v1.41 03/02/17 - Various edits for long runs
!                    - Field construction improved for large machines
!     v1.75 10/30/20 - Store iota_axis in file.
!                    - Added restarting from FIELDLINES run for FLD
!     v1.80 02/27/21 - Added interface to HINT magnetic fields
!-----------------------------------------------------------------------
      MODULE fieldlines_runtime
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE EZspline
!-----------------------------------------------------------------------
!     Module Variables
!          lverb         Logical to control screen output
!          lvmec         Logical to indicate VMEC Equilibria
!          lpies         Logical to indicate PIES Equilibria
!          lspec         Logical to indicate SPEC Equilibria
!          lcoil         Logical to indicate coil file
!          lmgrid        Logical to indicate MAKEGRID File
!          lmu           Logical to control fieldline diffusion
!          lvessel       Logical to indicate limiting surface
!          lvac          Logical to indicate vacuum fields only
!          lrestart      Logical to indicate a restart of FIELDLINES
!          laxis_i       Logical to control putting a coil at VMEC input axis
!          lraw          Logical to control use of raw currents in coils
!          nextcur       Number of external current systems
!          npoinc        Number of poincare sections
!          mu            Fieldline diffusion coefficient (mu=sqrt(D*tau*2)
!          dphi          Toroidal angle stepsize
!          follow_tol    Fieldlines following tollerance
!          pi            PI
!          pi2           2*PI
!          mu0           4*PI*10^-7
!          id_string     Equilibrium Filename
!          mgrid_string  MAKEGRID Filename
!          coil_string   Coil Filename
!          vessel_string Limiting Surface Filename
!          extcur        External currents for MGRID calculation
!          r(/z/phi)_hc  Starting points for homoclinic tangle calc
!          num_hcp       Number of points for homocline plot (half)
!          delta_hc      Length of initial homocline line (half length)
!----------------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER, PARAMETER ::  FILE_OPEN_ERR     = 1
      INTEGER, PARAMETER ::  ALLOC_ERR         = 11
      INTEGER, PARAMETER ::  NAMELIST_READ_ERR = 12
      INTEGER, PARAMETER ::  BAD_INPUT_ERR     = 13
      INTEGER, PARAMETER ::  VMEC_INPUT_ERR    = 2
      INTEGER, PARAMETER ::  VMEC_WOUT_ERR     = 21
      INTEGER, PARAMETER ::  MGRID_ERR         = 22
      INTEGER, PARAMETER ::  WALL_ERR          = 23
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
      
      INTEGER, PARAMETER ::  runtype_old       = 0
      INTEGER, PARAMETER ::  runtype_advanced  = 1
      INTEGER, PARAMETER ::  runtype_full      = 327
      INTEGER, PARAMETER ::  runtype_norun     = 328
      
      INTEGER, PARAMETER ::  MAXLINES   = 2**18
      INTEGER, PARAMETER ::  NLOCAL = 128  ! Number of local processors
      
      LOGICAL         :: lverb, lvmec, lpies, lspec, lcoil, lmgrid, &
                         lmu, lvessel, lvac, lrestart, laxis_i, &
                         ladvanced, lauto, lplasma_only, lbfield_only,&
                         lreverse, lhitonly, lafield_only, lraw, lemc3, &
                         lerror_field, lwall_trans, ledge_start, lnescoil,&
                         lmodb, lfield_start, lhint
      INTEGER         :: nextcur, npoinc, nruntype, num_hcp, &
                         nprocs_fieldlines, line_select
      REAL(rprec)     :: mu, dphi, follow_tol, pi, pi2, mu0, delta_hc, iota0
      REAL(rprec), DIMENSION(MAXLINES)     :: r_start, phi_start, &
                                              z_start, phi_end, &
                                              r_hc, z_hc, phi_hc
      REAL(rprec), ALLOCATABLE :: extcur(:)
      CHARACTER(256)  :: id_string, mgrid_string, coil_string, &
                         vessel_string, int_type, restart_string, &
                         nescoil_string
      
      REAL(rprec), PARAMETER :: FIELDLINES_VERSION = 1.75
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
      WRITE(6,*) '!!!!! ERROR !!!!!'
      
      IF (error_num .eq. FILE_OPEN_ERR) THEN
            WRITE(6,*) '  FIELDLINES COULD NOT OPEN A FILE.'
            WRITE(6,*) '  FILENAME: ',TRIM(string_val)
            WRITE(6,*) '  IERR:     ',ierr
      ELSEIF (error_num .eq. VMEC_INPUT_ERR) THEN
            WRITE(6,*) '  FIELDLINES COULD READ THE VMEC INDATA NAMELIST'
            WRITE(6,*) '  FILENAME: ',TRIM(string_val)
            WRITE(6,*) '  IERR:     ',ierr
      ELSEIF (error_num .eq. VMEC_WOUT_ERR) THEN
            WRITE(6,*) '  FIELDLINES COULD READ THE VMEC WOUT FILE'
            WRITE(6,*) '  FILENAME: ',TRIM(string_val)
            WRITE(6,*) '  IERR:     ',ierr
      ELSEIF (error_num .eq. ALLOC_ERR) THEN
            WRITE(6,*) '  FIELDLINES ENCOUNTERED AN ALLOCATION ERROR'
            WRITE(6,*) '  VARIABLES: ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num .eq. EZSPLINE_ERR) THEN
            WRITE(6,*) '  FIELDLINES ENCOUNTERED AN EZSPLINE ERROR'
            WRITE(6,*) '  ROUTINE/VAR: ',TRIM(string_val)
            WRITE(6,*) '  IERR:        ',ierr
            CALL EZspline_error(ierr)
      ELSEIF (error_num .eq. MGRID_ERR) THEN
            WRITE(6,*) '  FIELDLINES ENCOUNTERED AN ERROR READING THE MGRID FILE'
            WRITE(6,*) '  VARIABLES: ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num .eq. NETCDF_OPEN_ERR) THEN
            WRITE(6,*) '  FIELDLINES ENCOUNTERED AN ERROR OPENING A NETCDF FILE'
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
      ELSEIF (error_num .eq. NAMELIST_READ_ERR) THEN
            WRITE(6,*) '  FIELDLINES ENCOUNTERED AN ERROR READING A NAMELIST'
            WRITE(6,*) '  ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num .eq. WALL_ERR) THEN
            WRITE(6,*) '  WALL_MOD ENCOUTERED AN ERROR LOADING THE WALL'
            WRITE(6,*) '  ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num .eq. D02CJF_ERR) THEN
            WRITE(6,*) '  FIELDLINES ENCOUNTERED A NAG ERROR (D02CJF)'
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
            WRITE(6,*) '  FIELDLINES ENCOUNTERED A NAG ERROR (D05ADF)'
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
            WRITE(6,*) '  FIELDLINES ENCOUNTERED A NAG ERROR (C05AJF)'
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
            WRITE(6,*) '  FIELDLINES ENCOUNTERED AN LSODE ERROR'
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
            WRITE(6,*) '  FIELDLINES ENCOUNTERED A RKH68 ERROR'
            WRITE(6,*) '  ROUTINE:   ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num .eq. NAG_ERR) THEN
            WRITE(6,*) '  YOUR MACHINE DOES NOT SUPPORT NAG'
            WRITE(6,*) '  ROUTINE:   ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSEIF (error_num .eq. HDF5_ERR) THEN
            WRITE(6,*) '  FIELDLINES ENCOUNTERED AN HDF ERROR'
            WRITE(6,*) '  ROUTINE:   ',TRIM(string_val)
            WRITE(6,*) '  IERR:      ',ierr
      ELSE
           WRITE(6,*) '  FIELDLINES ENCOUNTERED AN UNKNOWN ERROR'
           WRITE(6,*) '  STRING: ',TRIM(string_val)
           WRITE(6,*) '  ierr:   ',ierr
      END IF
      CALL FLUSH(6)
#if defined(MPI_OPT)
      CALL MPI_FINALIZE(ierr)   
#endif
      STOP
      END SUBROUTINE handle_err

#if defined(MPI_OPT)
    SUBROUTINE FIELDLINES_TRANSMIT_2DDBL(n1,n2,m1,m2,data_in,n1_gbl,n2_gbl,id,root,COMM_local,ier)
    USE stel_kinds, ONLY: rprec
    USE mpi
    IMPLICIT NONE
    INTEGER, INTENT(in)           :: n1,n2,m1,m2,id,root,COMM_local
    INTEGER, INTENT(in)           :: n1_gbl, n2_gbl
    REAL(rprec), INTENT(inout)    :: data_in(n1:n2,m1:m2)
    INTEGER, INTENT(inout)        :: ier
    INTEGER                       :: nt_gbl
    DOUBLE PRECISION, ALLOCATABLE :: buffer_temp(:,:)

    IF (ier <0) RETURN
    nt_gbl=(n2_gbl-n1_gbl+1)*(m2-m1+1)
    ALLOCATE(buffer_temp(n1_gbl:n2_gbl,m1:m2))
    buffer_temp = 0
    buffer_temp(n1:n2,m1:m2) = data_in(n1:n2,m1:m2)
    IF (id == root) THEN
       CALL MPI_REDUCE(MPI_IN_PLACE,buffer_temp,nt_gbl,MPI_DOUBLE_PRECISION,MPI_SUM,&
                       root,COMM_local,ier)
       data_in = buffer_temp
    ELSE
         CALL MPI_REDUCE(buffer_temp,buffer_temp,nt_gbl,MPI_DOUBLE_PRECISION,MPI_SUM,&
                      root,COMM_local,ier)
    END IF
    DEALLOCATE(buffer_temp)
    ier = 0
    CALL MPI_BARRIER(COMM_local, ier)
    RETURN
    END SUBROUTINE FIELDLINES_TRANSMIT_2DDBL
#endif
      
      END MODULE fieldlines_runtime
