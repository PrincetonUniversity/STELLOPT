!     SPH010908 REPLACE INTEGER(iprec) with INTEGER
!*******************************************************************************
!  File diagnostic_cdf.f
!  Contains the module diagnostic_cdf
!    Module for defining variables and writing netCDF files, and reading
!    netCDF files with the derived types diagnostic_desc
!    (from the diagnostic_T module).
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
!       stel_constants
!       mddc_T
!       mddc_cdf
!       diagnostic_T
!       ezcdf
!       v3_utilities
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
!  MODULE diagnostic_cdf
!    
! SECTION I. VARIABLE DECLARATIONS
! SECTION II. INTERFACE BLOCKS
! SECTION III. DEFINE SUBROUTINES
! SECTION IV. WRITE SUBROUTINES
! SECTION V. READ SUBROUTINES
! SECTION VI. AUXILIARY FUNCTIONS AND SUBROUTINES

! SECTION X.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************

      MODULE diagnostic_cdf

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

      USE diagnostic_T
      USE mddc_T
      USE mddc_cdf
      USE ezcdf
      USE v3_utilities

!-------------------------------------------------------------------------------
!  Implicit None comes after USE statements, before other declarations
!-------------------------------------------------------------------------------
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Make type declarations and constants Private, so there are no conflicts.
!-------------------------------------------------------------------------------
      PRIVATE rprec, iprec, cprec, pi, twopi, one, zero

!-------------------------------------------------------------------------------
!  Variable Names for netCDF. Make them Private.
!-------------------------------------------------------------------------------

      CHARACTER (LEN=*), PRIVATE, PARAMETER ::                                 &
     &  vn_d_type = 'diagnostic_desc_d_type',                                  &
     &  vn_s_name = 'diagnostic_desc_s_name',                                  &
     &  vn_l_name = 'diagnostic_desc_l_name',                                  &
     &  vn_units = 'diagnostic_desc_units',                                    &         
     &  vn_sigma_default = 'diagnostic_desc_sigma_default'            

      CHARACTER (LEN=64), PRIVATE ::                                           &
     &  vn_d_type_use,                                                         &
     &  vn_s_name_use,                                                         &
     &  vn_l_name_use,                                                         &
     &  vn_units_use,                                                          &         
     &  vn_sigma_default_use                

      CHARACTER (LEN=*), PRIVATE, PARAMETER ::                                 &
     &  vn_desc_s_name = 'diagnostic_data_desc_s_name'

      CHARACTER (LEN=64), PRIVATE ::                                           &
     &  vn_desc_s_name_use

!-------------------------------------------------------------------------------
!  Lengths of Character Variables
!-------------------------------------------------------------------------------
      INTEGER, PARAMETER, PRIVATE :: type_len=10      
      INTEGER, PARAMETER, PRIVATE :: sn_len=30      
      INTEGER, PARAMETER, PRIVATE :: ln_len=80      
      INTEGER, PARAMETER, PRIVATE :: units_len=30      

!*******************************************************************************
! SECTION II. INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Generic Define
!-------------------------------------------------------------------------------
      INTERFACE diagnostic_cdf_define
         MODULE PROCEDURE diagnostic_cdf_define_desc
         END INTERFACE

!-------------------------------------------------------------------------------
!  Generic Write
!-------------------------------------------------------------------------------
      INTERFACE diagnostic_cdf_write
         MODULE PROCEDURE diagnostic_cdf_write_desc
         END INTERFACE

!-------------------------------------------------------------------------------
!  Generic Read
!-------------------------------------------------------------------------------
      INTERFACE diagnostic_cdf_read
         MODULE PROCEDURE diagnostic_cdf_read_desc
         END INTERFACE
!-------------------------------------------------------------------------------

      CONTAINS
!*******************************************************************************
! SECTION III. DEFINE SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE diagnostic_cdf_define_desc(this,iou,prefix)
!  Subroutine to do the appropriate netCDF definition calls for a diagnostic_desc
!  
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE (diagnostic_desc), INTENT (in)      :: this
      INTEGER, INTENT(in)                      :: iou
      CHARACTER (len=*), INTENT(in), OPTIONAL  :: prefix

!  this        diagnostic_desc - this is what gets defined and written.
!  iou         i/o unit number of the netCDF file. 
!  prefix      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'diagnostic_cdf_define_desc: '
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
      CALL diagnostic_cdf_defvn_desc(prefix_use)
         
! Define Components common to all d_types
      CALL cdf_define(iou, TRIM(vn_d_type_use), this % d_type)
      CALL cdf_define(iou, TRIM(vn_s_name_use), this % s_name)
      CALL cdf_define(iou, TRIM(vn_l_name_use), this % l_name)
      CALL cdf_define(iou, TRIM(vn_units_use), this % units)
      CALL cdf_define(iou, TRIM(vn_sigma_default_use),                         &
     &   this % sigma_default)

! Particular coding, depending on d_type

      SELECT CASE (this % d_type)
      
      CASE ('mddc') ! Magnetic Diagnostic
         CALL mddc_cdf_define(this % mddc,iou,prefix)

      CASE DEFAULT
         CALL err_fatal(sub_name // 'unrecognized d_type: ',                   &
     &      char = this % d_type)

      END SELECT
      
      RETURN
      
      END SUBROUTINE diagnostic_cdf_define_desc

!*******************************************************************************
! SECTION III. WRITE SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE diagnostic_cdf_write_desc(this,iou,prefix)
!  Subroutine to do the appropriate netCDF write calls for a diagnostic_desc
!  
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE (diagnostic_desc), INTENT (in)        :: this
      INTEGER, INTENT (in)                       :: iou
      CHARACTER (len=*), INTENT (in), OPTIONAL   :: prefix

!  this        diagnostic_desc - this is what gets written.
!  iou         i/o unit number of the netCDF file. 
!  prefix      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'diagnostic_cdf_write_desc: '
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
      CALL diagnostic_cdf_defvn_desc(prefix_use)

! Write Components
      CALL cdf_write(iou, TRIM(vn_d_type_use), this % d_type)
      CALL cdf_write(iou, TRIM(vn_s_name_use), this % s_name)
      CALL cdf_write(iou, TRIM(vn_l_name_use), this % l_name)
      CALL cdf_write(iou, TRIM(vn_units_use), this % units)
      CALL cdf_write(iou, TRIM(vn_sigma_default_use),                          &
     &   this % sigma_default)

! Particular coding, depending on d_type

      SELECT CASE (this % d_type)
      
      CASE ('mddc') ! Magnetic Diagnostic
         CALL mddc_cdf_write(this % mddc,iou,prefix)

      CASE DEFAULT
         CALL err_fatal(sub_name // 'unrecognized d_type: ',                   &
     &      char = this % d_type)

      END SELECT
      
      RETURN
      
      END SUBROUTINE diagnostic_cdf_write_desc

!*******************************************************************************
! SECTION IV. READ SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE diagnostic_cdf_read_desc(this,iou,prefix)
!  Subroutine to do the appropriate netCDF read calls for a diagnostic_desc
!  
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE (diagnostic_desc), INTENT (inout)        :: this
      INTEGER, INTENT (in)                       :: iou
      CHARACTER (len=*), INTENT (in), OPTIONAL   :: prefix

!  this        diagnostic_desc - this is what gets defined and written.
!  iou         i/o unit number of the netCDF file. 
!  prefix      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'diagnostic_cdf_read_desc: '
      CHARACTER(len=32) :: prefix_use
      INTEGER, DIMENSION(3) :: dimlens
      
      CHARACTER (len=type_len)   :: d_type
      CHARACTER (len=sn_len)     :: s_name                                 
      CHARACTER (len=ln_len)     :: l_name
      CHARACTER (len=units_len)  :: units                                 
      REAL(rprec)                :: sigma_default
      TYPE (mddc_desc)           :: mddc

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
      CALL diagnostic_cdf_defvn_desc(prefix_use)
         
! Read Components common to all d_types
! Note: Read in to variables local to this subroutine.
      CALL cdf_read(iou, TRIM(vn_d_type_use), d_type)
      CALL cdf_read(iou, TRIM(vn_s_name_use), s_name)
      CALL cdf_read(iou, TRIM(vn_l_name_use), l_name)
      CALL cdf_read(iou, TRIM(vn_units_use), units)
      CALL cdf_read(iou, TRIM(vn_sigma_default_use),                           &
     &   sigma_default)

! Particular coding, depending on d_type
      SELECT CASE (d_type)
      
      CASE ('mddc') ! Magnetic Diagnostic
         CALL mddc_cdf_read(mddc,iou,prefix)

!...............Create the diagnostic_desc, this...............!
         CALL diagnostic_construct(this,d_type,s_name,l_name,                     &
     &                             units,sigma_default,mddc)

      CASE DEFAULT
         CALL err_fatal(sub_name // 'unrecognized d_type: ',                   &
     &      char = d_type)

      END SELECT



!  Destroy the local mddc, to avoid memory leakage
      CALL mddc_destroy(mddc)
      
      RETURN
      
      END SUBROUTINE diagnostic_cdf_read_desc

!*******************************************************************************
! SECTION IV. AUXILIARY FUNCTIONS AND SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE diagnostic_cdf_defvn_desc(prefix_use)
!  Subroutine to do define the character variable names for a diagnostic_desc,
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
     &  'diagnostic_cdf_defvn_desc: '

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

! Define all variable names
      vn_d_type_use = diagnostic_cdf_mknam(prefix_use,vn_d_type)
      vn_s_name_use = diagnostic_cdf_mknam(prefix_use,vn_s_name)
      vn_l_name_use = diagnostic_cdf_mknam(prefix_use,vn_l_name)
      vn_units_use = diagnostic_cdf_mknam(prefix_use,vn_units)
      vn_sigma_default_use = diagnostic_cdf_mknam(prefix_use,                  &
     &   vn_sigma_default)
      
      RETURN
      
      END SUBROUTINE diagnostic_cdf_defvn_desc

!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      FUNCTION diagnostic_cdf_mknam(c1,c2)
! A simple function to help in the generation of names

!-----------------------------------------------
!   F u n c t i o n   N a m e
!-----------------------------------------------
      CHARACTER(LEN=64) diagnostic_cdf_mknam

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER(LEN=*), INTENT (in) :: c1,c2

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
      IF (LEN_TRIM(c1) .eq. 0) THEN
         diagnostic_cdf_mknam = TRIM(c2)
      ELSE
         diagnostic_cdf_mknam = ADJUSTL(TRIM(c1) // '_' // TRIM(c2))
      ENDIF
      
      RETURN
       
      END FUNCTION diagnostic_cdf_mknam

!-----------------------------------------------
!-----------------------------------------------
!*******************************************************************************
! SECTION X.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!  
!  09-02-2004 JDH - Initial Coding, based on bsc_cdf
!
!  JDH 12-11-2004
!     Modified _desc subroutines, to correspond with changes in diagnostic_T.
!
!  JDH 07-01-2005
!     Added bsc_destroy of local bsc_coil mdcoil in subroutine 
!     diagnostic_cdf_read_desc, to avoid memory leaks.
!
!  JDH 2007-06-11 - Altered coding to reflect definition of mddc_desc derived type
!
!  JDH 2008-01-19 - Added => null() to pointer declarations
!
!  JDH 2009-06-15
!    Eliminated diagnostic_data derived type - not needed.

      END MODULE diagnostic_cdf
