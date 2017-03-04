!     SPH010908 REPLACE INTEGER(iprec) with INTEGER
!*******************************************************************************
!  File signal_cdf.f
!  Contains the module signal_cdf
!    Module for defining variables and writing NETCDF files, and reading
!    netCDF files with the derived types signal_desc, signal_data, and 
!    signal_mrf (from the signal_T module).
!    
!    Information about the  EZcdf module is at:
!       http://w3.pppl.gov/NTCC/EZcdf/
!
!-------------------------------------------------------------------------------
!   DEPENDENCIES
!-------------------------------------------------------------------------------
!
!    This module uses the following modules:
!       stel_kinds
!       signal_T
!       diagnostic_T
!       diagnostic_cdf
!       ezcdf
!
!-------------------------------------------------------------------------------
!   CHANGE HISTORY
!-------------------------------------------------------------------------------
!
!  See Section X at the end of the module.
!
!-------------------------------------------------------------------------------
!   USAGE
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!   COMMENTS
!
! 1) All cdf_define calls must be completed before the first cdf_write call.
!-------------------------------------------------------------------------------
!
!*******************************************************************************

!*******************************************************************************
!  MODULE signal_cdf
!    
! SECTION I. VARIABLE DECLARATIONS
! SECTION II. INTERFACE BLOCKS
! SECTION III. DEFINE SUBROUTINES
! SECTION IV. WRITE SUBROUTINES
! SECTION V. READ SUBROUTINES
! SECTION VI AUXILIARY FUNCTIONS AND SUBROUTINES

! SECTION X.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************

      MODULE signal_cdf

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

!-------------------------------------------------------------------------------
!  Type declarations - lengths of reals, integers, and complexes.
!  Frequently used mathematical constants, lots of extra precision.
!-------------------------------------------------------------------------------
      USE stel_kinds, only : rprec, iprec, cprec
      USE stel_constants, only : pi, twopi, one, zero

!-------------------------------------------------------------------------------
!  Modules to USE
!-------------------------------------------------------------------------------

      USE bsc_T
      USE signal_T
      USE diagnostic_T
      USE diagnostic_cdf
      USE ezcdf
      USE v3_utilities

!-------------------------------------------------------------------------------
!  Implicit None comes after USE statements, before other declarations
!-------------------------------------------------------------------------------
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Lengths of Character Variables
!-------------------------------------------------------------------------------
      INTEGER, PARAMETER, PRIVATE :: type_len=10      
      INTEGER, PARAMETER, PRIVATE :: sn_len=30      
      INTEGER, PARAMETER, PRIVATE :: ln_len=80      
      INTEGER, PARAMETER, PRIVATE :: units_len=30      

!-------------------------------------------------------------------------------
!  Variable Names for netCDF. Make them Private.
!-------------------------------------------------------------------------------

      CHARACTER (LEN=*), PRIVATE, PARAMETER ::                                 &
     &  vn_s_type = 'signal_desc_s_type',                                      &
     &  vn_s_name = 'signal_desc_s_name',                                      &
     &  vn_l_name = 'signal_desc_l_name',                                      &
     &  vn_units = 'signal_desc_units',                                        &         
     &  vn_diag_s_name = 'signal_desc_diag_s_name'
     
      CHARACTER (LEN=64), PRIVATE ::                                           &
     &  vn_s_type_use,                                                         &
     &  vn_s_name_use,                                                         &
     &  vn_l_name_use,                                                         &
     &  vn_units_use,                                                          &         
     &  vn_diag_s_name_use

      CHARACTER (LEN=*), PRIVATE, PARAMETER ::                                 &
     &  vn_sd_type = 'signal_data_sd_type',                                    &
     &  vn_desc_s_name = 'signal_data_desc_s_name',                            &
     &  vn_data = 'signal_data_data',                                          &
     &  vn_sigma = 'signal_data_sigma'

      CHARACTER (LEN=64), PRIVATE ::                                           &
     &  vn_sd_type_use,                                                        &
     &  vn_desc_s_name_use,                                                    &
     &  vn_data_use,                                                           &
     &  vn_sigma_use

!*******************************************************************************
! SECTION II. INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Generic Define
!-------------------------------------------------------------------------------
      INTERFACE signal_cdf_define
         MODULE PROCEDURE signal_cdf_define_desc,                              &
     &                    signal_cdf_define_data
         END INTERFACE

!-------------------------------------------------------------------------------
!  Generic Write
!-------------------------------------------------------------------------------
      INTERFACE signal_cdf_write
         MODULE PROCEDURE signal_cdf_write_desc,                               &
     &                    signal_cdf_write_data
         END INTERFACE

!-------------------------------------------------------------------------------
!  Generic Read
!-------------------------------------------------------------------------------
      INTERFACE signal_cdf_read
         MODULE PROCEDURE signal_cdf_read_desc,                                &
     &                    signal_cdf_read_data
         END INTERFACE
!-------------------------------------------------------------------------------

      CONTAINS
!*******************************************************************************
! SECTION III. DEFINE SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE signal_cdf_define_desc(this,iou,prefix)
!  Subroutine to do the appropriate netCDF definition calls for a signal_desc
!  
!  For the pointer to the signal_diag, the s_name of the 
!  signal_diag is written in to the cdf file. (vn_diag_s_name)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE (signal_desc), INTENT (in)          :: this
      INTEGER, INTENT(in)                      :: iou
      CHARACTER (len=*), INTENT(in), OPTIONAL  :: prefix

!  this        signal_desc - this is what gets defined and written.
!  iou         i/o unit number of the netCDF file. 
!  prefix      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'signal_cdf_define_desc: '
      CHARACTER(len=32) :: prefix_use

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

! Define the prefix to actually use
      IF (PRESENT(prefix)) THEN
         prefix_use = TRIM(ADJUSTL(prefix))
      ELSE
         prefix_use = ' '
      ENDIF

! Define all vn_--_use variable names
      CALL signal_cdf_defvn_desc(prefix_use)
         
! Define Components
      CALL cdf_define(iou, TRIM(vn_s_type_use), this % s_type)
      CALL cdf_define(iou, TRIM(vn_s_name_use), this % s_name)
      CALL cdf_define(iou, TRIM(vn_l_name_use), this % l_name)
      CALL cdf_define(iou, TRIM(vn_units_use), this % units)
      IF (ASSOCIATED(this % diag)) THEN
         CALL cdf_define(iou, TRIM(vn_diag_s_name_use),                        &
     &          this % diag % s_name)
      ELSE
         CALL cdf_define(iou, TRIM(vn_diag_s_name_use),'null()')
      END IF 
            
      RETURN
      
      END SUBROUTINE signal_cdf_define_desc
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE signal_cdf_define_data(this,iou,prefix)
!  Subroutine to do the appropriate netCDF definition calls for a signal_data
!  
!  For the pointer to the signal_desc, the s_name of the signal_desc is
!  written in to the cdf file. (vn_desc_s_name)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE (signal_data), INTENT (in)          :: this
      INTEGER, INTENT(in)                      :: iou
      CHARACTER (len=*), INTENT(in), OPTIONAL  :: prefix

!  this        signal_data - this is what gets defined and written.
!  iou         i/o unit number of the netCDF file. 
!  prefix      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'signal_cdf_define_data: '
      CHARACTER(len=32) :: prefix_use

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

! Define the prefix to actually use
      IF (PRESENT(prefix)) THEN
         prefix_use = TRIM(ADJUSTL(prefix))
      ELSE
         prefix_use = ' '
      ENDIF

! Define all vn_--_use variable names
      CALL signal_cdf_defvn_data(prefix_use)
         
! Define Components
      CALL cdf_define(iou, TRIM(vn_sd_type_use), this % sd_type)
      IF (ASSOCIATED(this % desc)) THEN
         CALL cdf_define(iou, TRIM(vn_desc_s_name_use),                   &
     &          this % desc % s_name)
      ELSE
         CALL cdf_define(iou, TRIM(vn_desc_s_name_use),'null()')
      END IF 
      CALL cdf_define(iou, TRIM(vn_data_use), this % data)
      CALL cdf_define(iou, TRIM(vn_sigma_use), this % sigma)
      
      RETURN
      
      END SUBROUTINE signal_cdf_define_data

!*******************************************************************************
! SECTION III. WRITE SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE signal_cdf_write_desc(this,iou,prefix)
!  Subroutine to do the appropriate netCDF write calls for a signal_desc
!  
!  For the pointer to the signal_diag, the s_name of the 
!  signal_diag is written in to the cdf file. (vn_diag_s_name)
!
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE (signal_desc), INTENT (in)            :: this
      INTEGER, INTENT (in)                       :: iou
      CHARACTER (len=*), INTENT (in), OPTIONAL   :: prefix

!  this        signal_desc - this is what gets written.
!  iou         i/o unit number of the netCDF file. 
!  prefix      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'signal_cdf_write_desc: '
      CHARACTER(len=32) :: prefix_use

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
! Define the prefix to actually use
      IF (PRESENT(prefix)) THEN
         prefix_use = TRIM(ADJUSTL(prefix))
      ELSE
         prefix_use = ' '
      ENDIF

! Define all vn_--_use variable names
      CALL signal_cdf_defvn_desc(prefix_use)

! Write Components
      CALL cdf_write(iou, TRIM(vn_s_type_use), this % s_type)
      CALL cdf_write(iou, TRIM(vn_s_name_use), this % s_name)
      CALL cdf_write(iou, TRIM(vn_l_name_use), this % l_name)
      CALL cdf_write(iou, TRIM(vn_units_use), this % units)
      IF (ASSOCIATED(this % diag)) THEN
         CALL cdf_write(iou, TRIM(vn_diag_s_name_use),                         &
     &          this % diag % s_name)
      ELSE
         CALL cdf_write(iou, TRIM(vn_diag_s_name_use),'null()')
      END IF 
      
      RETURN
      
      END SUBROUTINE signal_cdf_write_desc
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE signal_cdf_write_data(this,iou,prefix)
!  Subroutine to do the appropriate netCDF write calls for a signal_desc
!  
!  For the pointer to the signal_desc, the s_name of the signal_desc is
!  written in to the cdf file. (vn_desc_s_name)
!
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE (signal_data), INTENT (in)            :: this
      INTEGER, INTENT (in)                       :: iou
      CHARACTER (len=*), INTENT (in), OPTIONAL   :: prefix

!  this        signal_data - this is what gets written.
!  iou         i/o unit number of the netCDF file. 
!  prefix      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'signal_cdf_write_data: '
      CHARACTER(len=32) :: prefix_use

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
! Define the prefix to actually use
      IF (PRESENT(prefix)) THEN
         prefix_use = TRIM(ADJUSTL(prefix))
      ELSE
         prefix_use = ' '
      ENDIF

! Define all vn_--_use variable names
      CALL signal_cdf_defvn_data(prefix_use)

! Write Components
      CALL cdf_write(iou, TRIM(vn_sd_type_use), this % sd_type)
      IF (ASSOCIATED(this % desc)) THEN
         CALL cdf_write(iou, TRIM(vn_desc_s_name_use),                         &
     &          this % desc % s_name)
      ELSE
         CALL cdf_write(iou, TRIM(vn_desc_s_name_use),'null()')
      END IF 
      CALL cdf_write(iou, TRIM(vn_data_use), this % data)
      CALL cdf_write(iou, TRIM(vn_sigma_use), this % sigma)
      
      RETURN
      
      END SUBROUTINE signal_cdf_write_data

!*******************************************************************************
! SECTION IV. READ SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE signal_cdf_read_desc(this,iou,prefix,psname)
!  Subroutine to do the appropriate netCDF read calls for a signal_desc
!  
!  For the pointer to the diagnostic_desc, the s_name of the diagnostic_desc is
!  read from to the cdf file (vn_diag_s_name). It is communicated to the
!  outside of this subroutine via the optional argument psname.

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE (signal_desc), INTENT (inout)         :: this
      INTEGER, INTENT (in)                       :: iou
      CHARACTER (len=*), INTENT (in), OPTIONAL   :: prefix
      CHARACTER (len=*), INTENT (out), OPTIONAL  :: psname

!  this        signal_desc - this is what gets defined and written.
!  iou         i/o unit number of the netCDF file. 
!  prefix      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names
!  psname      character - the s_name of the signal_diag this % diag
!              This subroutine can't point to that signal_diag. That has to
!              be done outside  of this subroutine.
!              If this % diag was not associated, then psname = 'null()'
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'signal_cdf_read_desc: '
      CHARACTER(len=32) :: prefix_use
      INTEGER, DIMENSION(3) :: dimlens
      CHARACTER(len=64) :: psname_use
      
      CHARACTER (len=type_len)         :: s_type
      CHARACTER (len=sn_len)           :: s_name                                 
      CHARACTER (len=ln_len)           :: l_name
      CHARACTER (len=units_len)        :: units                                 
      TYPE (diagnostic_desc), TARGET   :: diag 

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

! Define the prefix to actually use
      IF (PRESENT(prefix)) THEN
         prefix_use = TRIM(ADJUSTL(prefix))
      ELSE
         prefix_use = ' '
      ENDIF

! Define all vn_--_use variable names
      CALL signal_cdf_defvn_desc(prefix_use)
         
! Read Components common to all d_types
! Note: Read in to variables local to this subroutine.
      CALL cdf_read(iou, TRIM(vn_s_type_use), s_type)
      CALL cdf_read(iou, TRIM(vn_s_name_use), s_name)
      CALL cdf_read(iou, TRIM(vn_l_name_use), l_name)
      CALL cdf_read(iou, TRIM(vn_units_use), units)
      CALL cdf_read(iou, TRIM(vn_diag_s_name_use), psname_use)

!  Create the signal_desc, this

!  JDH 12-11-2004. Have problem with null() argument here
!      diag => null()

!  JDH 12-13-2004. Changed diag from a pointer to just a local structure (with
!  the TARGET attribute)
!  It has NOT been constructed. The correct pointer assignment will have to 
!  be made later. This is NOT clean code. It should be fixed up.
!  [Optional argument, and a subroutine to do the pointer assignment]

      CALL signal_construct(this,s_type,s_name,l_name, units, diag)
      
      IF (PRESENT(psname)) THEN
         psname = psname_use
      ENDIF
      
      RETURN
      
      END SUBROUTINE signal_cdf_read_desc
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE signal_cdf_read_data(this,iou,prefix,psname)
!  Subroutine to do the appropriate netCDF read calls for a signal_data
!  
!  For the pointer to the signal_desc, the s_name of the signal_desc is
!  read from to the cdf file (vn_desc_s_name). It is communicated to the
!  outside of this subroutine via the optional argument psname.

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE (signal_data), INTENT (inout)         :: this
      INTEGER, INTENT (in)                       :: iou
      CHARACTER (len=*), INTENT (in), OPTIONAL   :: prefix
      CHARACTER (len=*), INTENT (out), OPTIONAL  :: psname

!  this        signal_data - this is what gets defined and written.
!  iou         i/o unit number of the netCDF file. 
!  prefix      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names
!  psname      character - the s_name of the signal_desc this % desc
!              This subroutine can't point to that signal_desc. That has to
!              be done outside  of this subroutine.
!              If this % desc was not associated, then psname = 'null()'
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(len=*), PARAMETER  :: sub_name =                               &
     &  'signal_cdf_read_data: '
      CHARACTER(len=32)            :: prefix_use
      INTEGER, DIMENSION(3) :: dimlens
      INTEGER               :: n_data, n_sigma
      INTEGER               :: ier1
      CHARACTER(len=64)            :: psname_use
      
      CHARACTER(len=type_len)            :: sd_type
      TYPE (signal_desc), POINTER        :: desc => null()
      REAL(rprec), DIMENSION(:), POINTER :: data => null()
      REAL(rprec), DIMENSION(:), POINTER :: sigma => null()

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

! Define the prefix to actually use
      IF (PRESENT(prefix)) THEN
         prefix_use = TRIM(ADJUSTL(prefix))
      ELSE
         prefix_use = ' '
      ENDIF

! Define all vn_--_use variable names
      CALL signal_cdf_defvn_data(prefix_use)
         
! Read Components
! Note: Read in to variables local to this subroutine.
! Arrays require inquiry regarding size, and allocation before actual reading.

      CALL cdf_read(iou, TRIM(vn_sd_type_use), sd_type)
      
      CALL cdf_read(iou, TRIM(vn_desc_s_name_use),psname_use)

      CALL cdf_inquire(iou, TRIM(vn_data_use),dimlens)
      n_data = dimlens(1)
      CALL assert_eq(1,dimlens(2),dimlens(3),                                  &
     &   sub_name // 'Unexpected data dimensions')
      ALLOCATE(data(n_data),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc data')
      CALL cdf_read(iou, TRIM(vn_data_use), data)

      CALL cdf_inquire(iou, TRIM(vn_sigma_use),dimlens)
      n_sigma = dimlens(1)
      CALL assert_eq(1,dimlens(2),dimlens(3),                                  &
     &   sub_name // 'Unexpected sigma dimensions')
      ALLOCATE(sigma(n_data),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc sigma')
      CALL cdf_read(iou, TRIM(vn_sigma_use), sigma)

! Create the signal_data, this
      desc => null()
      CALL signal_data_construct(this,desc,sd_type,data,sigma)
      
      IF (PRESENT(psname)) THEN
         psname = psname_use
      ENDIF

! Deallocate the local pointers
      DEALLOCATE(data,sigma,STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'dealloc')
      
      RETURN
      
      END SUBROUTINE signal_cdf_read_data


!*******************************************************************************
! SECTION IV. AUXILIARY FUNCTIONS AND SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE signal_cdf_defvn_desc(prefix_use)
!  Subroutine to do define the character variable names for a signal_desc,
!  using the prefix. All the vn_ variables are module variables, and so do not
!   need to be declared here

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER (len=*), INTENT (in)   :: prefix_use

!  prefix_use      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'signal_cdf_defvn_desc: '

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

! Define all variable names
      vn_s_type_use = signal_cdf_mknam(prefix_use,vn_s_type)
      vn_s_name_use = signal_cdf_mknam(prefix_use,vn_s_name)
      vn_l_name_use = signal_cdf_mknam(prefix_use,vn_l_name)
      vn_units_use = signal_cdf_mknam(prefix_use,vn_units)
      vn_diag_s_name_use = signal_cdf_mknam(prefix_use,                        &
     &   vn_diag_s_name)
      
      RETURN
      
      END SUBROUTINE signal_cdf_defvn_desc
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE signal_cdf_defvn_data(prefix_use)
!  Subroutine to do define the character variable names for a signal_desc,
!  using the prefix. All the vn_ variables are module variables, and so do not
!   need to be declared here

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER (len=*), INTENT (in)   :: prefix_use

!  prefix_use      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'signal_cdf_defvn_data: '

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

! Define all variable names
      vn_sd_type_use = signal_cdf_mknam(prefix_use,vn_sd_type)
      vn_desc_s_name_use = signal_cdf_mknam(prefix_use,vn_desc_s_name)
      vn_data_use = signal_cdf_mknam(prefix_use,vn_data)
      vn_sigma_use = signal_cdf_mknam(prefix_use,vn_sigma)
      
      RETURN
      
      END SUBROUTINE signal_cdf_defvn_data
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      FUNCTION signal_cdf_mknam(c1,c2)
! A simple function to help in the generation of names

!-----------------------------------------------
!   F u n c t i o n   N a m e
!-----------------------------------------------
      CHARACTER(LEN=64) signal_cdf_mknam

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER(LEN=*), INTENT (in) :: c1,c2

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
      IF (LEN_TRIM(c1) .eq. 0) THEN
         signal_cdf_mknam = TRIM(c2)
      ELSE
         signal_cdf_mknam = ADJUSTL(TRIM(c1) // '_' // TRIM(c2))
      ENDIF
      
      RETURN
       
      END FUNCTION signal_cdf_mknam

!-----------------------------------------------
!-----------------------------------------------
!*******************************************************************************
! SECTION X.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!  
!!  09-02-2004 JDH - Initial Coding, based on bsc_cdf
!
!  JDH 12-11-2004
!     Modified _desc and _mrf subroutines, to correspond with changes in
!     signal_T. Still have problems with null() arguments.
!
!  JDH 07-01-2005
!     Did signal destroy of local mrf in signal_cdf_read_desc, to avoid
!     memory leakage.
!
!  JDH 2007-06-12
!     Eliminated _mrf. 
!     Not sure of details with diagnostic_desc component, pointer, target, etc.
!
!  JDH 2008-01-19 - Added => null() to pointer declarations

      END MODULE signal_cdf
