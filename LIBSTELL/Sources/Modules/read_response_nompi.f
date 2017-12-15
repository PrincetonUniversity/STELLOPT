!*******************************************************************************
!  File read_response.f
!      Module for use in program v3post
!      Encapsulates the reading of the netCDF files
!         for the coil responses and the plasma responses.
!      Response data is generated n code v3rfun.
!      This module is a companion to module write_response
!      NB.  Geometric data in the cdf files IS NOT read back, only response data is read.
!-------------------------------------------------------------------------------
!   DEPENDENCIES
!-------------------------------------------------------------------------------
!
!    This module uses the following modules:
!       stel_kinds
!       stel_constants
!       bsc_cdf
!       response_arrays ! No Longer JDH 6.28.03
!       ezcdf
!       bsc
!
!-------------------------------------------------------------------------------
!   CHANGE HISTORY
!-------------------------------------------------------------------------------
!
!  Initial Coding - Ed Lazarus 12.12.2002
! 01.11.2003 SPH - added istat error handler, in case file can not be opened
! 04.11.2003 SPH - added mmode arg to check "Raw" or "Scaled" mode used for extcur array
! 06.27.2003 JDH - Add kp_store to TYPE prfun. Add lots of comments and white space
! 06.28.2003 JDH - Eliminated get_response. Only used once, changed call in v3post.f
! 06.29.2003 JDH - Eliminated cdf_vctpot_read. In cdf_prfun_read, the vector potential 
!   arrays a_r, a_f, a_z that are pointer components of the derived structure prfun,
!   now have space allocated within cdf_prfun_read. 
! 07.08.2003 SPH - Added logical flag ldim_only to cdf_prfun_read so only dimensions are read
!   in, NOT a-pointer components (needed for initializations).
! 08.15.2003 JDH - Revised logic in cdf_prfun_read so that dimensional information
!   (including s_name) would ALWAYS be read in. Other small changes.
!-------------------------------------------------------------------------------
!   USAGE
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!   COMMENTS
!-------------------------------------------------------------------------------
!
!
!*******************************************************************************
!  MODULE read_response
!    
! SECTION I. VARIABLE DECLARATIONS
! SECTION II. INTERFACE BLOCKS
! SECTION III.  COIL RESPONSE FUNCTION READ
! SECTION IV. PLASMA RESPONSE FUNCTION READ
!*******************************************************************************

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

      MODULE read_response_nompi

! work_arrays avoids premature size definition of arrays
      USE stel_kinds
      USE stel_constants
#if defined(NETCDF)
      USE bsc_cdf, ONLY: vn_c_type, vn_s_name, vn_l_name,
     1  vn_current, vn_raux, vn_xnod, vn_ehnod, vn_rcirc,
     1  vn_xcent, vn_enhat
!      USE response_arrays
      USE ezcdf
#endif
      
!  Define derived type clresfun - coil response function
      TYPE clresfun
         INTEGER(iprec) :: n_field_cg
         INTEGER(iprec) :: n_diagn_c
         REAL(rprec), DIMENSION(:,:), POINTER :: rdiag_coilg => null()
      END TYPE clresfun

!  Define derived type prfun
!  Contains vector potential data for a single magnetic diagnostic
!  TYPE prfun is 247+8*size(a) bytes (+ size of kp_store, JDH)
      TYPE prfun
         CHARACTER (len=63) :: cdffil
         CHARACTER (len=120) :: name_diagnostic_dot
         CHARACTER (len=30) :: idrfun
         INTEGER(iprec) :: idc
         INTEGER(iprec) :: ir
         INTEGER(iprec) :: jz
         INTEGER(iprec) :: kp    
         INTEGER(iprec) :: kp_store    
         INTEGER(iprec) :: n_field_periods
         REAL(rprec) :: rmin
         REAL(rprec) :: rmax       
         REAL(rprec) :: zmin
         REAL(rprec) :: zmax
         REAL(rprec), DIMENSION(:,:,:), POINTER :: a_r => null()
         REAL(rprec), DIMENSION(:,:,:), POINTER :: a_f => null()
         REAL(rprec), DIMENSION(:,:,:), POINTER :: a_z => null()
         CHARACTER (len=30) :: s_name
         CHARACTER (len=8) :: c_type
         LOGICAL :: lstell_sym
      END TYPE prfun

! vn_ are the variable names loaded into NETCDF
      CHARACTER (LEN=*), PARAMETER :: 
     &  vn_rdiag_coilg = 'rdiag_coilg',                                        &
     &  vn_n_field_cg = 'n_field_cg',                                          &
     &  vn_n_diagn_c = 'n_diagn_c',                                            &
     &  vn_inductance_coilg ='external_inductance'
      CHARACTER (LEN=*), PARAMETER ::                                          &
     &  vn_idrfun = 'idrfun',                                                  &
     &  vn_name_diagnostic_dot = 'name_diagnostic_dot',                        &
     &  vn_lstell_sym = 'lstell_sym',                                          &
     &  vn_ir = 'ir',                                                          &
     &  vn_jz = 'jz',                                                          &
     &  vn_kp = 'kp',                                                          & 
     &  vn_kp_store = 'kp_store',                                              & 
     &  vn_idc = 'idc',                                                        &
     &  vn_cdffil = 'cdffil',                                                  &
     &  vn_a_r = 'a_r',                                                        &
     &  vn_a_f = 'a_f',                                                        &
     &  vn_a_z = 'a_z',                                                        & 
     &  vn_plasma_response = 'plasma_response',                                & 
     &  vn_n_field_periods = 'n_field_periods',                                &
     &  vn_rmin = 'rmin',                                                      &
     &  vn_rmax = 'rmax',                                                      &
     &  vn_zmin = 'zmin',                                                      &
     &  vn_zmax = 'zmax'
       

!  Other variables
!      INTEGER :: nprfun = 0

#if defined(NETCDF)
!*******************************************************************************
! SECTION II. INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------

      CONTAINS

!*******************************************************************************
! SECTION III. COIL RESPONSE FUNCTION READ
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE cdf_crfun_read_nompi(cdffil, crf, istat, mmode)

! cdf_crfun_read reads back the (sensor,coil) vector potential table
! (ncrfun, cdffil) are unit number and file name. Called by get_response.
! 01.11.2003 SPH 
! added istat error handler, in case file can not be opened
! 04.11.2003 SPH
! added mmode arg to check "Raw" or "Scaled" mode used for extcur array
! ASSUMES that V3POST output for vn_rdiag_coilg is in "Raw" mode
!  JDH 06.28.2003. Added coil response function argument (derived type)
!  SPH 05.12.2005  Added MPI enabling code
! this version for xv3post => NO MPI! EAL
      IMPLICIT NONE

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER (LEN=*), INTENT(in) :: cdffil
      TYPE(clresfun), INTENT(out)   :: crf
      INTEGER, INTENT(out)        :: istat
      CHARACTER(LEN=1), OPTIONAL, INTENT(in) :: mmode
!  cdffil      name of the netcdf file that contains the coil response function
!  crf         Derived Type clresfun variable - contains the coil response function
!  istat       status variable
!  mmode       variable to specify mode for "Scaled' or "Raw"
!              in "S" mode, B-fields are "unit" current responses, so EXTCUR
!              has units of A. The responses are obtained from SUM(inductance * EXTCUR)
!              in "R" mode, B-fields are the true fields corresponding to
!              the currents in the coils-dot file, so EXTCUR ~ unity dimensionless
!              multiplier. Thus, the responses are SUM(rdiag_coil * EXTCUR)

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(LEN=1) :: mode = 'N'
      INTEGER          :: ncrfun, nwprocs
      LOGICAL          :: bReadIO
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
      IF (PRESENT(mmode)) mode = mmode
      nwprocs = 0
      istat = 0

      bReadIO = .true.

      IF (bReadIO) THEN
         CALL cdf_open(ncrfun, cdffil, 'r', istat)
         IF (istat .ne. 0) STOP 'ERROR OPENING CDFFIL!'
         CALL cdf_read(ncrfun, vn_n_field_cg, crf%n_field_cg)
         CALL cdf_read(ncrfun, vn_n_diagn_c, crf%n_diagn_c)
      END IF
      
      IF (ASSOCIATED(crf%rdiag_coilg)) DEALLOCATE(crf%rdiag_coilg)
      ALLOCATE(crf%rdiag_coilg(crf%n_diagn_c,crf%n_field_cg),                  &
     &   stat = istat)
      IF (istat .ne. 0) THEN
          WRITE(6,*) 'In cdf_crfun_read, istat = ', istat
          RETURN
       END IF

      IF (bReadIO) THEN
         IF (mode .eq. 'R' .or. mode .eq. 'N') THEN
            CALL cdf_read(ncrfun, vn_rdiag_coilg, crf%rdiag_coilg)
         ELSE
            CALL cdf_read(ncrfun, vn_inductance_coilg, crf%rdiag_coilg)
         END IF
         CALL cdf_close(ncrfun, istat)
      END IF


      END SUBROUTINE cdf_crfun_read_nompi

!*******************************************************************************
! SECTION IV. PLASMA RESPONSE FUNCTION READ
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------

      SUBROUTINE cdf_prfun_read_nompi(cdffil, pl_str, istat, ldim_only)
      
!  Subroutine cdf_pr_fun reads a NETCDF file, and puts all the plasma
!  response function information into the defined type plrfun.

! cdf_prfun_read reads NETCDF file of vector potential for single sensor
! (nprfun, cdffil, idc, pl_str) unit#, filename, sensor#, structure for
! returning data. Called by get_response.
! SPH added istat error handler, in case file can not be opened
! SPH (05/12/05) added MPI_ logic
      IMPLICIT NONE
!      USE bsc

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      CHARACTER (LEN=*), INTENT(in) :: cdffil
      TYPE (prfun), INTENT(inout)   :: pl_str
      INTEGER, INTENT(out)          :: istat
      LOGICAL, INTENT(in), OPTIONAL :: ldim_only

!  cdffil      name of the plasma response netcdf file
!  pl_str      type prfun variable where the plasma response is stored
!  istat       integer status variable
!  ldim_only   logical, true if only need to read the dimension information. 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!      INTEGER         :: dims(3) 
!      CHARACTER(LEN=5)     :: xtype
      LOGICAL         :: ldim_only_local, lfile
      INTEGER         :: nprfun, nwprocs
      LOGICAL         :: bReadIO
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

      ldim_only_local = .false.
      IF (PRESENT(ldim_only)) ldim_only_local = ldim_only
      nwprocs = 0
      istat = 0

      bReadIO = .true.

      IF (bReadIO) THEN
         INQUIRE(file=cdffil, exist=lfile)
         IF (.not. lfile) PRINT*,"file=",cdffil
         IF (.not. lfile) STOP 'cdf file not found in CDF_PRFUN_READ!'
!         DO WHILE (istat .ne. 0)
         CALL cdf_open(nprfun, cdffil, 'r', istat)
!         END DO

!  Read the scalar (dim) information
         CALL cdf_read(nprfun,vn_s_name,pl_str%s_name)
         CALL cdf_read(nprfun,vn_kp,pl_str%kp)
      
         CALL cdf_read(nprfun,vn_kp_store,pl_str%kp_store,istat) ! note istat
         CALL cdf_read(nprfun,vn_lstell_sym,pl_str%lstell_sym)
! Coding to take care of files generated before the stellarator symmetry
! was correctly implemented.
         IF (istat .ne. 0) THEN
            pl_str%kp_store = pl_str%kp                          
            pl_str%lstell_sym = .false.
         END IF

         CALL cdf_read(nprfun,vn_jz,pl_str%jz)
         CALL cdf_read(nprfun,vn_ir,pl_str%ir)
         CALL cdf_read(nprfun,vn_rmin,pl_str%rmin)
         CALL cdf_read(nprfun,vn_rmax,pl_str%rmax)
         CALL cdf_read(nprfun,vn_zmin,pl_str%zmin)
         CALL cdf_read(nprfun,vn_zmax,pl_str%zmax)
         CALL cdf_read(nprfun,vn_cdffil,pl_str%cdffil)
         CALL cdf_read(nprfun,vn_idc,pl_str%idc)
         CALL cdf_read(nprfun,vn_name_diagnostic_dot,                          &
     &                 pl_str%name_diagnostic_dot)   
         CALL cdf_read(nprfun,vn_n_field_periods,pl_str%n_field_periods)
         CALL cdf_read(nprfun,vn_idrfun,pl_str%idrfun)
      END IF


      IF (ldim_only_local) THEN
         IF (bReadIO) CALL cdf_close(nprfun, istat)
         RETURN
      END IF
         
      IF (ASSOCIATED(pl_str%a_r)) DEALLOCATE(pl_str%a_r)
      IF (ASSOCIATED(pl_str%a_f)) DEALLOCATE(pl_str%a_f)
      IF (ASSOCIATED(pl_str%a_z)) DEALLOCATE(pl_str%a_z)
      ALLOCATE(pl_str%a_r(pl_str%ir,pl_str%jz,pl_str%kp_store))
      ALLOCATE(pl_str%a_f(pl_str%ir,pl_str%jz,pl_str%kp_store))
      ALLOCATE(pl_str%a_z(pl_str%ir,pl_str%jz,pl_str%kp_store))

      IF (bReadIO) THEN
         CALL cdf_read(nprfun,vn_a_r,pl_str%a_r)
         CALL cdf_read(nprfun,vn_a_f,pl_str%a_f)
         CALL cdf_read(nprfun,vn_a_z,pl_str%a_z)
         CALL cdf_close(nprfun, istat)
      END IF


      END SUBROUTINE cdf_prfun_read_nompi
#endif
      END MODULE read_response_nompi
