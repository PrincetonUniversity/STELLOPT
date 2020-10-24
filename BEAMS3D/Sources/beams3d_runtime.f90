!-----------------------------------------------------------------------
!     Module:        beams3d_runtime
!     Authors:       S. Lazerson (lazerson@pppl.gov), M. McMillan (matthew.mcmillan@my.wheaton.edu)
!     Date:          06/20/2012
!     Description:   This module contains various runtime parameters
!                    for the BEAMS3D code.  It also contains the
!                    handle_error subroutine.
!     v1.05 03/10/15 - Fixed RKH68 and LSODE routines to work properly.
!     v1.10 09/28/15 - Optimized calling splines in physics and fpart.
!                    - Added shine through and current drive output
!                    - Corrected ioninzation glitch
!                    - Neutral beams may start outside vacuum field box.
!                    - Diagnostics added for optimization.
!     v1.11 03/29/16 - Fixed bug with -lflux incrementing i.
!     v1.20 06/03/16 - Beam power incorporated along with ndot profile calculation (real units).
!     v1.20 08/10/16 - More efficient array passing methods implemented.
!                    - Max stepsize set to 0.05 m
!                    - Diagnostics fully working
!                    - Added Wall output including number of strikes per cell.
!     v1.50 08/29/16 - Switch away from que type data transfer to gatherv statements.
!                    - Collisional slowing down operators check for accuracy.
!                    - S and U coordinate now defaulted on.  -flux depreicated
!     v1.51 09/13/16 - Substep size implemnted for neutral deposition
!     v1.52 11/22/16 - Added ability to model W7-X injector geometry
!     v2.00 05/06/19 - Shared Memory model implemented
!     v2.01 08/21/19 - Added ASCOT Interface
!     v2.50 01/31/20 - Added Heating, Deposition, and Dist FUNCTION
!                    - Tested ASCOT4 and ASCOT5 interfaces
!                    - Now output wall and shinethrough heat flux
!     v2.70 06/09/20 - Converged on format of dV/drho factor
!                    - ASCOT5 interface updated
!                    - DIST5D implemented
!     v2.80 09/24/20 - DIST5D Normalized to phase space volume
!                    - J and dense now calculated in diagnostics
!-----------------------------------------------------------------------
MODULE beams3d_runtime
    !-----------------------------------------------------------------------
    !     Libraries
    !-----------------------------------------------------------------------
    USE stel_kinds, ONLY: rprec
    USE mpi_params
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
    !          lrestart      Logical to indicate a restart of BEAMS3D
    !          laxis_i       Logical to control putting a coil at VMEC input axis
    !          nextcur       Number of external current systems
    !          npoinc        Number of poincare sections
    !          mu            Trajectory diffusion coefficient
    !          dt            Time angle stepsize
    !          follow_tol    Trajectory following tollerance
    !          pi            PI
    !          pi2           2*PI
    !          mu0           4*PI*10^-7
    !          id_string     Equilibrium Filename
    !          mgrid_string  MAKEGRID Filename
    !          coil_string   Coil Filename
    !          vessel_string Limiting Surface Filename
    !          extcur        External currents for MGRID calculation
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

    DOUBLE PRECISION, PARAMETER :: one           = 1.0D0 ! 1.0

    LOGICAL :: lverb, lvmec, lpies, lspec, lcoil, lmgrid, &
               lvessel, lvac, lrestart_grid, lrestart_particles, lneut, &
               lbeam, lhitonly, lread_input, lplasma_only, lraw,&
               ldepo, lbeam_simple, ldebug, lcollision, lw7x, lsuzuki, &
               lascot, lascot4, lbbnbi, lvessel_beam, lascotfl, lrandomize, &
               lfusion, lfusion_alpha
    INTEGER :: nextcur, npoinc, nbeams, nparticles_start, nprocs_beams
    INTEGER, DIMENSION(MAXBEAMS) :: Dex_beams
    INTEGER, ALLOCATABLE :: beam(:)
    REAL(rprec) :: dt, follow_tol, pi, pi2, invpi2, mu0, to3, dt_save, &
                   ne_scale, te_scale, ti_scale, zeff_scale, fusion_scale, &
                   lendt_m
    REAL(rprec), DIMENSION(MAXBEAMS) :: Adist_beams, Asize_beams, Div_beams, E_beams, mass_beams, &
                                        charge_beams, Zatom_beams, P_beams
    REAL(rprec), DIMENSION(MAXBEAMS, 2) :: r_beams, z_beams, phi_beams
    REAL(rprec), DIMENSION(MAXPROFLEN) :: TE_AUX_S, TE_AUX_F, NE_AUX_S, NE_AUX_F, TI_AUX_S, TI_AUX_F,&
                                            POT_AUX_S, POT_AUX_F, ZEFF_AUX_S, ZEFF_AUX_F
    REAL(rprec), DIMENSION(MAXPARTICLES) :: r_start_in, phi_start_in, z_start_in, vll_start_in, &
                                            & mu_start_in, charge_in, Zatom_in, mass_in, t_end_in
    REAL(rprec), ALLOCATABLE :: R_start(:), phi_start(:), Z_start(:), vll_start(:), v_neut(:,:), mu_start(:), &
                                & mass(:), charge(:), Zatom(:), t_end(:), weight(:)
    REAL(rprec), ALLOCATABLE :: extcur(:)
    CHARACTER(LEN=10) ::  qid_str_saved ! For ASCOT5
    CHARACTER(256) :: id_string, mgrid_string, coil_string, &
    vessel_string, int_type, restart_string, bbnbi_string

    REAL(rprec), PARAMETER :: BEAMS3D_VERSION = 2.80
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
            WRITE(6, *) '  BEAMS3D COULD NOT OPEN A FILE.'
            WRITE(6, *) '  FILENAME: ', TRIM(string_val)
            WRITE(6, *) '  IERR:     ', ierr
        ELSEIF (error_num .eq. VMEC_INPUT_ERR) THEN
            WRITE(6, *) '  BEAMS3D COULD NOT READ THE VMEC INDATA NAMELIST'
            WRITE(6, *) '  FILENAME: ', TRIM(string_val)
            WRITE(6, *) '  IERR:     ', ierr
        ELSEIF (error_num .eq. VMEC_WOUT_ERR) THEN
            WRITE(6, *) '  BEAMS3D COULD NOT READ THE VMEC WOUT FILE'
            WRITE(6, *) '  FILENAME: ', TRIM(string_val)
            WRITE(6, *) '  IERR:     ', ierr
        ELSEIF (error_num .eq. ALLOC_ERR) THEN
            WRITE(6, *) '  BEAMS3D ENCOUNTERED AN ALLOCATION ERROR'
            WRITE(6, *) '  VARIABLES: ', TRIM(string_val)
            WRITE(6, *) '  IERR:      ', ierr
        ELSEIF (error_num .eq. EZSPLINE_ERR) THEN
            WRITE(6, *) '  BEAMS3D ENCOUNTERED AN EZSPLINE ERROR'
            WRITE(6, *) '  ROUTINE/VAR: ', TRIM(string_val)
            WRITE(6, *) '  IERR:        ', ierr
            CALL EZspline_error(ierr)
        ELSEIF (error_num .eq. MGRID_ERR) THEN
            WRITE(6, *) '  BEAMS3D ENCOUNTERED AN ERROR READING THE MGRID FILE'
            WRITE(6, *) '  VARIABLES: ', TRIM(string_val)
            WRITE(6, *) '  IERR:      ', ierr
        ELSEIF (error_num .eq. INIT_BEAMS_ERR) THEN
            WRITE(6, *) '  BEAMS3D ENCOUNTERED AN ERROR INITIALIZING THE BEAMS'
            WRITE(6, *) '  VARIABLES: ', TRIM(string_val)
            WRITE(6, *) '  IERR:      ', ierr
        ELSEIF (error_num .eq. ADAS_ERR) THEN
            WRITE(6, *) '  BEAMS3D ENCOUNTERED AN ERROR CALLING ADAS'
            WRITE(6, *) '  VARIABLES: ', TRIM(string_val)
            WRITE(6, *) '  IERR:      ', ierr
        ELSEIF (error_num .eq. BAD_BEAMDEX_ERR) THEN
            WRITE(6, *) '  BEAMS3D ENCOUNTERED AN NAMELIST ERROR'
            WRITE(6, *) '    -BEAMLET USED BUT NO DEX_BEAM SET!'
        ELSEIF (error_num .eq. NETCDF_OPEN_ERR) THEN
            WRITE(6, *) '  BEAMS3D ENCOUNTERED AN ERROR OPENING A NETCDF FILE'
            WRITE(6, *) '  VARIABLES: ', TRIM(string_val)
            WRITE(6, *) '  IERR:      ', ierr
        ELSEIF (error_num .eq. HDF5_OPEN_ERR) THEN
            WRITE(6, *) '  BEAMS3D ENCOUNTERED AN ERROR OPENING AN HDF5 FILE'
            WRITE(6, *) '  FILENAME:  ', TRIM(string_val)
            WRITE(6, *) '  IERR:      ', ierr
        ELSEIF (error_num .eq. HDF5_WRITE_ERR) THEN
            WRITE(6, *) '  BEAMS3D ENCOUNTERED AN ERROR WRITING TO AN HDF5 FILE'
            WRITE(6, *) '  VARNAME:   ', TRIM(string_val)
            WRITE(6, *) '  IERR:      ', ierr
        ELSEIF (error_num .eq. HDF5_READ_ERR) THEN
            WRITE(6, *) '  BEAMS3D ENCOUNTERED AN ERROR READING FROM AN HDF5 FILE'
            WRITE(6, *) '  VARNAME:   ', TRIM(string_val)
            WRITE(6, *) '  IERR:      ', ierr
        ELSEIF (error_num .eq. HDF5_CLOSE_ERR) THEN
            WRITE(6, *) '  BEAMS3D ENCOUNTERED AN ERROR CLOSING AN HDF5 FILE'
            WRITE(6, *) '  FILENAME:  ', TRIM(string_val)
            WRITE(6, *) '  IERR:      ', ierr
        ELSEIF (error_num .eq. NAMELIST_READ_ERR) THEN
            WRITE(6, *) '  BEAMS3D ENCOUNTERED AN ERROR READING A NAMELIST'
            WRITE(6, *) '  ', TRIM(string_val)
            WRITE(6, *) '  IERR:      ', ierr
        ELSEIF (error_num .eq. D02CJF_ERR) THEN
            WRITE(6, *) '  BEAMS3D ENCOUNTERED A NAG ERROR (D02CJF)'
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
            WRITE(6, *) '  BEAMS3D ENCOUNTERED A NAG ERROR (D05ADF)'
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
            WRITE(6, *) '  BEAMS3D ENCOUNTERED A NAG ERROR (C05AJF)'
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
            WRITE(6, *) '  BEAMS3D ENCOUNTERED AN LSODE ERROR'
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
            WRITE(6, *) '  BEAMS3D ENCOUNTERED A RKH68 ERROR'
            WRITE(6, *) '  ROUTINE:   ', TRIM(string_val)
            WRITE(6, *) '  IERR:      ', ierr
        ELSEIF (error_num .eq. NAG_ERR) THEN
            WRITE(6, *) '  YOUR MACHINE DOES NOT SUPPORT NAG'
            WRITE(6, *) '  ROUTINE:   ', TRIM(string_val)
            WRITE(6, *) '  IERR:      ', ierr
        ELSEIF (error_num .eq. HDF5_ERR) THEN
            WRITE(6, *) '  BEAMS3D ENCOUNTERED AN HDF ERROR'
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
            WRITE(6, *) '  BEAMS3D ENCOUNTERED AN UNKNOWN ERROR'
            WRITE(6, *) '  STRING: ', TRIM(string_val)
            WRITE(6, *) '  ierr:   ', ierr
        END IF
        CALL FLUSH(6)
#if defined(MPI_OPT)
        CALL MPI_BARRIER(MPI_COMM_BEAMS,ierr_mpi)
        ALLOCATE(error_array(1:nprocs_beams))
        ierr_mpi = 0
        CALL MPI_ALLGATHER(error_num,1,MPI_INTEGER,error_array,1,MPI_INTEGER,MPI_COMM_BEAMS,ierr_mpi)
        ierr_mpi = 0
        CALL MPI_BARRIER(MPI_COMM_BEAMS,ierr_mpi)
        IF (ANY(error_array .ne. 0)) CALL MPI_FINALIZE(ierr_mpi)
        DEALLOCATE(error_array)
        RETURN
#else
        IF (error_num .eq. MPI_CHECK) RETURN
#endif
        PRINT *,myworkid,' calling STOP'
        CALL FLUSH(6)
        STOP
    END SUBROUTINE handle_err

#if defined(MPI_OPT)
    SUBROUTINE BEAMS3D_TRANSMIT_2DDBL(n1,n2,m1,m2,data_in,nproc,mnum,moffsets,id,root,COMM_local,ier)
    USE stel_kinds, ONLY: rprec
    USE mpi_inc
    IMPLICIT NONE
    INTEGER, INTENT(in)           :: n1,n2,m1,m2,nproc,id,root,COMM_local
    INTEGER, INTENT(in)           :: mnum(nproc), moffsets(nproc)
    REAL(rprec), INTENT(inout)    :: data_in(n1:n2,m1:m2)
    INTEGER, INTENT(inout)        :: ier
    INTEGER, PARAMETER            :: ndims=2
    INTEGER, PARAMETER            :: sstart(2) = (/0,0/) ! Starting offsets
    INTEGER                       :: dbl_size, localsize, ARRAY_SEND_TYPE, RESIZED_ARRAY_SEND_TYPE
    INTEGER                       :: asize(ndims), ssize(ndims), mrec(nproc)
    INTEGER(KIND=MPI_ADDRESS_KIND):: low_bound,extent
    DOUBLE PRECISION, ALLOCATABLE :: buffer_temp(:,:)

    IF (ier <0) RETURN
    mrec = 1
    ! Size of the data array
    ssize(1) = n2-n1+1
    ssize(2) = m2-m1+1
    localsize = mnum(id+1) ! Number of elements in block of data
    ALLOCATE(buffer_temp(ssize(1),ssize(2)))  ! Buffer is the local array (global for master)
    buffer_temp(1:ssize(1),1:ssize(2)) = data_in(n1:n2,m1:m2)
    asize    = ssize  ! Size of Global array (BCAST from master later)
    ! Set asize for everyone
    CALL MPI_BCAST(asize, 2, MPI_INTEGER, root, COMM_local, ier); IF (ier <0) RETURN
    ! Define the subarray
    CALL MPI_TYPE_CREATE_SUBARRAY(ndims,asize,ssize,sstart,MPI_ORDER_FORTRAN,&
                                  MPI_DOUBLE_PRECISION,ARRAY_SEND_TYPE,ier); IF (ier <0) RETURN
    CALL MPI_TYPE_COMMIT(ARRAY_SEND_TYPE,ier); IF (ier <0) RETURN
    ! Define Resized Type array
    CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, dbl_size,ier); IF (ier<0) RETURN
    low_bound = 0
    extent = dbl_size
    CALL MPI_TYPE_CREATE_RESIZED(ARRAY_SEND_TYPE,low_bound,extent,RESIZED_ARRAY_SEND_TYPE,ier); IF (ier<0) RETURN
    CALL MPI_TYPE_COMMIT(RESIZED_ARRAY_SEND_TYPE,ier); IF (ier <0) RETURN
    ! Now gather the variable to master
    IF (id == root) THEN
       localsize = PRODUCT(ssize)
       mrec(1)   = localsize
       CALL MPI_GATHERV(MPI_IN_PLACE,1,MPI_DOUBLE_PRECISION,&
                        data_in,mrec, moffsets,RESIZED_ARRAY_SEND_TYPE,&
                        root,COMM_local,ier)
       !data_in(n1:n2,m1:m2) = buffer_temp(n1:n2,m1:m2)  ! Make sure we pass back the full array.
    ELSE
       CALL MPI_GATHERV(buffer_temp,localsize,MPI_DOUBLE_PRECISION,&
                        buffer_temp,mrec, moffsets,RESIZED_ARRAY_SEND_TYPE,&
                        root,COMM_local,ier)
    END IF
    ! Free the types
    CALL MPI_TYPE_FREE(ARRAY_SEND_TYPE,ier); IF (ier <0) RETURN
    CALL MPI_TYPE_FREE(RESIZED_ARRAY_SEND_TYPE,ier); IF (ier <0) RETURN
    DEALLOCATE(buffer_temp)
    ier = 0
    CALL MPI_BARRIER(COMM_local, ier)
    RETURN
    END SUBROUTINE BEAMS3D_TRANSMIT_2DDBL

    SUBROUTINE BEAMS3D_TRANSMIT_2DLOG(n1,n2,m1,m2,data_in,nproc,mnum,moffsets,id,root,COMM_local,ier)
    USE stel_kinds, ONLY: rprec
    USE mpi_inc
    IMPLICIT NONE
    INTEGER, INTENT(in)           :: n1,n2,m1,m2,nproc,id,root,COMM_local
    INTEGER, INTENT(in)           :: mnum(nproc), moffsets(nproc)
    LOGICAL, INTENT(inout)    :: data_in(n1:n2,m1:m2)
    INTEGER, INTENT(inout)        :: ier
    INTEGER, PARAMETER            :: ndims=2
    INTEGER, PARAMETER            :: sstart(2) = (/0,0/) ! Starting offsets
    INTEGER                       :: log_size, localsize, ARRAY_SEND_TYPE, RESIZED_ARRAY_SEND_TYPE
    INTEGER                       :: asize(ndims), ssize(ndims), mrec(nproc)
    INTEGER(KIND=MPI_ADDRESS_KIND):: low_bound,extent
    LOGICAL, ALLOCATABLE :: buffer_temp(:,:)

    IF (ier <0) RETURN
    mrec = 1
    ! Size of the data array
    ssize(1) = n2-n1+1
    ssize(2) = m2-m1+1
    localsize = mnum(id+1) ! Number of elements in block of data
    ALLOCATE(buffer_temp(ssize(1),ssize(2)))  ! Buffer is the local array (global for master)
    buffer_temp(1:ssize(1),1:ssize(2)) = data_in(n1:n2,m1:m2)
    asize    = ssize  ! Size of Global array (BCAST from master later)
    ! Set asize for everyone
    CALL MPI_BCAST(asize, 2, MPI_INTEGER, root, COMM_local, ier); IF (ier <0) RETURN
    ! Define the subarray
    CALL MPI_TYPE_CREATE_SUBARRAY(ndims,asize,ssize,sstart,MPI_ORDER_FORTRAN,&
                                  MPI_LOGICAL,ARRAY_SEND_TYPE,ier); IF (ier <0) RETURN
    CALL MPI_TYPE_COMMIT(ARRAY_SEND_TYPE,ier); IF (ier <0) RETURN
    ! Define Resized Type array
    CALL MPI_TYPE_SIZE(MPI_LOGICAL, log_size,ier); IF (ier<0) RETURN
    low_bound = 0
    extent = log_size
    CALL MPI_TYPE_CREATE_RESIZED(ARRAY_SEND_TYPE,low_bound,extent,RESIZED_ARRAY_SEND_TYPE,ier); IF (ier<0) RETURN
    CALL MPI_TYPE_COMMIT(RESIZED_ARRAY_SEND_TYPE,ier); IF (ier <0) RETURN
    ! Now gather the variable to master
    IF (id == root) THEN
       localsize = PRODUCT(ssize)
       mrec(1)   = localsize
       CALL MPI_GATHERV(MPI_IN_PLACE,1,MPI_LOGICAL,&
                        data_in,mrec, moffsets,RESIZED_ARRAY_SEND_TYPE,&
                        root,COMM_local,ier)
       !data_in(n1:n2,m1:m2) = buffer_temp(n1:n2,m1:m2)  ! Make sure we pass back the full array.
    ELSE
       CALL MPI_GATHERV(buffer_temp,localsize,MPI_LOGICAL,&
                        buffer_temp,mrec, moffsets,RESIZED_ARRAY_SEND_TYPE,&
                        root,COMM_local,ier)
    END IF
    ! Free the types
    CALL MPI_TYPE_FREE(ARRAY_SEND_TYPE,ier); IF (ier <0) RETURN
    CALL MPI_TYPE_FREE(RESIZED_ARRAY_SEND_TYPE,ier); IF (ier <0) RETURN
    DEALLOCATE(buffer_temp)
    ier = 0
    CALL MPI_BARRIER(COMM_local, ier)
    RETURN
    END SUBROUTINE BEAMS3D_TRANSMIT_2DLOG
#endif

END MODULE beams3d_runtime
