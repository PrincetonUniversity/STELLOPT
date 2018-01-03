!     SPH010908 REPLACE INTEGER(iprec) with INTEGER
!*******************************************************************************
!  File mddc_cdf.f
!  Contains the module mddc_cdf
!    Module for defining variables and writing netCDF files, and reading
!    netCDF files with the derived types mddc_desc and mddc_mrf
!    (from the mddc_T module).
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
!       bsc
!       bsc_cdf
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
!  MODULE mddc_cdf
!    
! SECTION I. VARIABLE DECLARATIONS
! SECTION II. INTERFACE BLOCKS
! SECTION III. DEFINE SUBROUTINES
! SECTION IV. WRITE SUBROUTINES
! SECTION V. READ SUBROUTINES
! SECTION VI. AUXILIARY FUNCTIONS AND SUBROUTINES

! SECTION X.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************

      MODULE mddc_cdf

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
      USE bsc_cdf, only : bsc_cdf_define_coil, bsc_cdf_write_coil,
     &   bsc_cdf_read_coil
      USE mddc_T
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
     &  vn_s_name = 'mddc_desc_s_name',                                        &            
     &  vn_l_name = 'mddc_desc_l_name',                                        &            
     &  vn_units = 'mddc_desc_units',                                          &            
     &  vn_sigma_default = 'mddc_desc_sigma_default',                          &            
     &  vn_l_mdcoil_def = 'mddc_desc_l_mdcoil_def',                            &            
     &  vn_mddc_type = 'mddc_desc_mddc_type',                                  &            
     &  vn_flux_factor = 'mddc_desc_flux_factor'

      CHARACTER (LEN=64), PRIVATE ::                                           &
     &  vn_s_name_use,                                                         &
     &  vn_l_name_use,                                                         &
     &  vn_units_use,                                                          &         
     &  vn_sigma_default_use,                                                  &            
     &  vn_l_mdcoil_def_use,                                                   &            
     &  vn_mddc_type_use,                                                      &            
     &  vn_flux_factor_use                   

      CHARACTER (LEN=*), PRIVATE, PARAMETER ::                                 &
     &  vn_desc_s_name = 'mddc_data_desc_s_name'

      CHARACTER (LEN=64), PRIVATE ::                                           &
     &  vn_desc_s_name_use

      CHARACTER (LEN=*), PRIVATE, PARAMETER ::                                 &
     &  vn_code_name = 'mddc_mrf_code_name',                                   &
     &  vn_code_version = 'mddc_mrf_code_version',                             &
     &  vn_date_run = 'mddc_mrf_date_run',                                     &         
     &  vn_field_coils_id = 'mddc_mrf_field_coils_id',                         &         
     &  vn_n_field_cg = 'mddc_mrf_n_field_cg',                                 &         
     &  vn_rdiag_coilg_1 = 'mddc_mrf_rdiag_coilg_1',                           &         
     &  vn_extcur_mg = 'mddc_mrf_extcur_mg',                                   &         
     &  vn_ir = 'mddc_mrf_ir',                                                 &         
     &  vn_jz = 'mddc_mrf_jz',                                                 &         
     &  vn_kp = 'mddc_mrf_kp',                                                 &         
     &  vn_kp_store = 'mddc_mrf_kp_store',                                     &         
     &  vn_rmin = 'mddc_mrf_rmin',                                             &         
     &  vn_rmax = 'mddc_mrf_rmax',                                             &         
     &  vn_zmin = 'mddc_mrf_zmin',                                             &         
     &  vn_zmax = 'mddc_mrf_zmax',                                             &         
     &  vn_n_field_periods = 'mddc_mrf_n_field_periods',                       &         
     &  vn_lstell_sym = 'mddc_mrf_lstell_sym',                                 &         
     &  vn_a_r = 'mddc_mrf_a_r',                                               &         
     &  vn_a_f = 'mddc_mrf_a_f',                                               &         
     &  vn_a_z = 'mddc_mrf_a_z',                                               &
     &  vn_use_con_shell = 'mddc_mrf_use_con_shell',                           &
     &  vn_a_s_r = 'mddc_mrf_a_s_r',                                           &
     &  vn_a_s_f = 'mddc_mrf_a_s_f',                                           &
     &  vn_a_s_z = 'mddc_mrf_a_s_z',                                           &
     &  vn_kp_shell = 'mddc_mrf_kp_shell',                                     &
     &  vn_kp_shell_store = 'mddc_mrf_kp_shell_store'

      CHARACTER (LEN=64), PRIVATE ::                                           &
     &  vn_code_name_use,                                                      &
     &  vn_code_version_use,                                                   &
     &  vn_date_run_use,                                                       &
     &  vn_field_coils_id_use,                                                 &
     &  vn_n_field_cg_use,                                                     &
     &  vn_rdiag_coilg_1_use,                                                  &
     &  vn_extcur_mg_use,                                                      &
     &  vn_ir_use,                                                             &
     &  vn_jz_use,                                                             &
     &  vn_kp_use,                                                             &
     &  vn_kp_store_use,                                                       &
     &  vn_rmin_use,                                                           &
     &  vn_rmax_use,                                                           &
     &  vn_zmin_use,                                                           &
     &  vn_zmax_use,                                                           &
     &  vn_n_field_periods_use,                                                &
     &  vn_lstell_sym_use,                                                     &
     &  vn_a_r_use,                                                            &
     &  vn_a_f_use,                                                            &
     &  vn_a_z_use,                                                            &
     &  vn_use_con_shell_use,                                                  &
     &  vn_a_s_r_use,                                                          &
     &  vn_a_s_f_use,                                                          &
     &  vn_a_S_z_use,                                                          &
     &  vn_kp_shell_use,                                                       &
     &  vn_kp_shell_store_use

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
      INTERFACE mddc_cdf_define
         MODULE PROCEDURE mddc_cdf_define_desc,                                &
     &                    mddc_cdf_define_mrf
         END INTERFACE

!-------------------------------------------------------------------------------
!  Generic Write
!-------------------------------------------------------------------------------
      INTERFACE mddc_cdf_write
         MODULE PROCEDURE mddc_cdf_write_desc,                                 &
     &                    mddc_cdf_write_mrf
         END INTERFACE

!-------------------------------------------------------------------------------
!  Generic Read
!-------------------------------------------------------------------------------
      INTERFACE mddc_cdf_read
         MODULE PROCEDURE mddc_cdf_read_desc,                                  &
     &                    mddc_cdf_read_mrf
         END INTERFACE
!-------------------------------------------------------------------------------

      CONTAINS
!*******************************************************************************
! SECTION III. DEFINE SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE mddc_cdf_define_desc(this,iou,prefix)
!  Subroutine to do the appropriate netCDF definition calls for a mddc_desc
!  
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE (mddc_desc), INTENT (in)            :: this
      INTEGER, INTENT(in)                      :: iou
      CHARACTER (len=*), INTENT(in), OPTIONAL  :: prefix

!  this        mddc_desc - this is what gets defined and written.
!  iou         i/o unit number of the netCDF file. 
!  prefix      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'mddc_cdf_define_desc: '
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
      CALL mddc_cdf_defvn_desc(prefix_use)
         
! Scalar Components
      CALL cdf_define(iou, TRIM(vn_s_name_use), this % s_name)
      CALL cdf_define(iou, TRIM(vn_l_name_use), this % l_name)
      CALL cdf_define(iou, TRIM(vn_units_use), this % units)
      CALL cdf_define(iou, TRIM(vn_sigma_default_use),                         &
     &   this % sigma_default)
      CALL cdf_define(iou, TRIM(vn_l_mdcoil_def_use),                          &
     &      this % l_mdcoil_def)
      CALL cdf_define(iou, TRIM(vn_mddc_type_use),                             &
     &      this % mddc_type)
      CALL cdf_define(iou, TRIM(vn_flux_factor_use),                           &
     &      this % flux_factor)

!  bsc_coil and mrf
      IF (this % l_mdcoil_def) THEN
            CALL bsc_cdf_define_coil(this % mdcoil,iou,prefix)
      END IF 
     
      CALL mddc_cdf_define(this % mrf,iou,prefix)   
      
      RETURN
      
      END SUBROUTINE mddc_cdf_define_desc
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE mddc_cdf_define_mrf(this,iou,prefix)
!  Subroutine to do the appropriate netCDF definition calls for a mddc_mrf

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE (mddc_mrf), INTENT (in)           :: this
      INTEGER, INTENT(in)                      :: iou
      CHARACTER (len=*), INTENT(in), OPTIONAL  :: prefix

!  this        mddc_mrf - this is what gets defined and written.
!  iou         i/o unit number of the netCDF file. 
!  prefix      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'mddc_cdf_define_mrf: '
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
      CALL mddc_cdf_defvn_mrf(prefix_use)
         
! Define Components
      CALL cdf_define(iou, TRIM(vn_code_name_use), this % code_name)
      CALL cdf_define(iou, TRIM(vn_code_version_use),                          &
     &   this % code_version)
      CALL cdf_define(iou, TRIM(vn_date_run_use), this % date_run)
      CALL cdf_define(iou, TRIM(vn_field_coils_id_use),                        &
     &   this % field_coils_id)
      CALL cdf_define(iou, TRIM(vn_n_field_cg_use),                            &
     &   this % n_field_cg)
      CALL cdf_define(iou, TRIM(vn_rdiag_coilg_1_use),                         &
     &   this % rdiag_coilg_1)
      CALL cdf_define(iou, TRIM(vn_extcur_mg_use), this % extcur_mg)
      CALL cdf_define(iou, TRIM(vn_ir_use), this % ir)
      CALL cdf_define(iou, TRIM(vn_jz_use), this % jz)
      CALL cdf_define(iou, TRIM(vn_kp_use), this % kp)
      CALL cdf_define(iou, TRIM(vn_kp_store_use), this % kp_store)
      CALL cdf_define(iou, TRIM(vn_rmin_use), this % rmin)
      CALL cdf_define(iou, TRIM(vn_rmax_use), this % rmax)
      CALL cdf_define(iou, TRIM(vn_zmin_use), this % zmin)
      CALL cdf_define(iou, TRIM(vn_zmax_use), this % zmax)
      CALL cdf_define(iou, TRIM(vn_n_field_periods_use),                       &
     &   this % n_field_periods)
      CALL cdf_define(iou, TRIM(vn_lstell_sym_use),                            &
     &   this % lstell_sym)
      CALL cdf_define(iou, TRIM(vn_a_r_use), this % a_r)
      CALL cdf_define(iou, TRIM(vn_a_f_use), this % a_f)
      CALL cdf_define(iou, TRIM(vn_a_z_use), this % a_z)

      CALL cdf_define(iou, TRIM(vn_use_con_shell_use),                         &
     &                this%use_con_shell)
      IF (this%use_con_shell) THEN
         CALL cdf_define(iou, TRIM(vn_kp_shell_use), this % kp_shell)
         CALL cdf_define(iou, TRIM(vn_a_s_r_use), this % a_s_r)
         CALL cdf_define(iou, TRIM(vn_a_s_f_use), this % a_s_f)
         CALL cdf_define(iou, TRIM(vn_a_s_z_use), this % a_s_z)
      END IF

      RETURN

      END SUBROUTINE mddc_cdf_define_mrf

!*******************************************************************************
! SECTION III. WRITE SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE mddc_cdf_write_desc(this,iou,prefix)
!  Subroutine to do the appropriate netCDF write calls for a mddc_desc
!  
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE (mddc_desc), INTENT (in)        :: this
      INTEGER, INTENT (in)                       :: iou
      CHARACTER (len=*), INTENT (in), OPTIONAL   :: prefix

!  this        mddc_desc - this is what gets written.
!  iou         i/o unit number of the netCDF file. 
!  prefix      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'mddc_cdf_write_desc: '
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
      CALL mddc_cdf_defvn_desc(prefix_use)

! Scalar Components
      CALL cdf_write(iou, TRIM(vn_s_name_use), this % s_name)
      CALL cdf_write(iou, TRIM(vn_l_name_use), this % l_name)
      CALL cdf_write(iou, TRIM(vn_units_use), this % units)
      CALL cdf_write(iou, TRIM(vn_sigma_default_use),                         &
     &   this % sigma_default)
      CALL cdf_write(iou, TRIM(vn_l_mdcoil_def_use),                          &
     &      this % l_mdcoil_def)
      CALL cdf_write(iou, TRIM(vn_mddc_type_use),                             &
     &      this % mddc_type)
      CALL cdf_write(iou, TRIM(vn_flux_factor_use),                           &
     &      this % flux_factor)

!  bsc_coil and mrf
      IF (this % l_mdcoil_def) THEN
            CALL bsc_cdf_write_coil(this % mdcoil,iou,prefix)
      END IF 
     
      CALL mddc_cdf_write(this % mrf,iou,prefix)   
      
      RETURN
      
      END SUBROUTINE mddc_cdf_write_desc
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE mddc_cdf_write_mrf(this,iou,prefix)
!  Subroutine to do the appropriate netCDF write calls for a mddc_mrf

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE (mddc_mrf), INTENT (in)           :: this
      INTEGER, INTENT(in)                      :: iou
      CHARACTER (len=*), INTENT(in), OPTIONAL  :: prefix

!  this        mddc_mrf - this is what gets defined and written.
!  iou         i/o unit number of the netCDF file. 
!  prefix      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'mddc_cdf_write_mrf: '
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
      CALL mddc_cdf_defvn_mrf(prefix_use)
         
! Write Components
      CALL cdf_write(iou, TRIM(vn_code_name_use), this % code_name)
      CALL cdf_write(iou, TRIM(vn_code_version_use),                           &
     &   this % code_version)
      CALL cdf_write(iou, TRIM(vn_date_run_use), this % date_run)
      CALL cdf_write(iou, TRIM(vn_field_coils_id_use),                         &
     &   this % field_coils_id)
      CALL cdf_write(iou, TRIM(vn_n_field_cg_use),                             &
     &   this % n_field_cg)
      CALL cdf_write(iou, TRIM(vn_rdiag_coilg_1_use),                          &
     &   this % rdiag_coilg_1)
      CALL cdf_write(iou, TRIM(vn_extcur_mg_use), this % extcur_mg)
      CALL cdf_write(iou, TRIM(vn_ir_use), this % ir)
      CALL cdf_write(iou, TRIM(vn_jz_use), this % jz)
      CALL cdf_write(iou, TRIM(vn_kp_use), this % kp)
      CALL cdf_write(iou, TRIM(vn_kp_store_use), this % kp_store)
      CALL cdf_write(iou, TRIM(vn_rmin_use), this % rmin)
      CALL cdf_write(iou, TRIM(vn_rmax_use), this % rmax)
      CALL cdf_write(iou, TRIM(vn_zmin_use), this % zmin)
      CALL cdf_write(iou, TRIM(vn_zmax_use), this % zmax)
      CALL cdf_write(iou, TRIM(vn_n_field_periods_use),                        &
     &   this % n_field_periods)
      CALL cdf_write(iou, TRIM(vn_lstell_sym_use),                             &
     &   this % lstell_sym)
      CALL cdf_write(iou, TRIM(vn_a_r_use), this % a_r)
      CALL cdf_write(iou, TRIM(vn_a_f_use), this % a_f)
      CALL cdf_write(iou, TRIM(vn_a_z_use), this % a_z)

      CALL cdf_write(iou, TRIM(vn_use_con_shell_use),                          &
     &               this%use_con_shell)
      IF (this%use_con_shell) THEN
         CALL cdf_write(iou, TRIM(vn_kp_shell_use), this % kp_shell)
         CALL cdf_write(iou, TRIM(vn_kp_shell_store_use),                      &
     &                  this % kp_shell_store)

         CALL cdf_write(iou, TRIM(vn_a_s_r_use), this % a_s_r)
         CALL cdf_write(iou, TRIM(vn_a_s_f_use), this % a_s_f)
         CALL cdf_write(iou, TRIM(vn_a_s_z_use), this % a_s_z)
      END IF

      RETURN

      END SUBROUTINE mddc_cdf_write_mrf

!*******************************************************************************
! SECTION IV. READ SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE mddc_cdf_read_desc(this,iou,prefix)
!  Subroutine to do the appropriate netCDF read calls for a mddc_desc
!  
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE (mddc_desc), INTENT (inout)        :: this
      INTEGER, INTENT (in)                       :: iou
      CHARACTER (len=*), INTENT (in), OPTIONAL   :: prefix

!  this        mddc_desc - this is what gets defined and written.
!  iou         i/o unit number of the netCDF file. 
!  prefix      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'mddc_cdf_read_desc: '
      CHARACTER(len=32) :: prefix_use
      INTEGER, DIMENSION(3) :: dimlens
      
      CHARACTER (len=sn_len)    :: s_name                                 
      CHARACTER (len=ln_len)    :: l_name
      CHARACTER (len=units_len) :: units                                 
      CHARACTER (len=30)        :: mddc_type                                 
      LOGICAL                   :: l_mdcoil_def
      REAL(rprec)               :: sigma_default
      REAL(rprec)               :: flux_factor
      TYPE (bsc_coil)           :: mdcoil
      TYPE (mddc_mrf)           :: mrf      

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
      CALL mddc_cdf_defvn_desc(prefix_use)
         
! Read Scalar pieces
! Note: Read in to variables local to this subroutine.
      CALL cdf_read(iou, TRIM(vn_s_name_use), s_name)
      CALL cdf_read(iou, TRIM(vn_l_name_use), l_name)
      CALL cdf_read(iou, TRIM(vn_units_use), units)

      CALL cdf_read(iou, TRIM(vn_sigma_default_use), sigma_default)
      CALL cdf_read(iou, TRIM(vn_l_mdcoil_def_use),l_mdcoil_def)
      CALL cdf_read(iou, TRIM(vn_mddc_type_use),mddc_type)
      CALL cdf_read(iou, TRIM(vn_flux_factor_use),flux_factor)

!  Read Derived Types
      IF (l_mdcoil_def) THEN
         CALL bsc_cdf_read_coil(mdcoil,iou,prefix)
      ENDIF
      
      CALL mddc_cdf_read(mrf,iou,prefix)   
      

! Create the mddc_desc, this
      CALL mddc_desc_construct(this,s_name,l_name,units,                       &
     &   sigma_default,mddc_type,mdcoil,mrf,flux_factor)

!  Destroy the local bsc_coil mdcoil, and mrf to avoid memory leakage
      CALL bsc_destroy(mdcoil)
      CALL mddc_destroy(mrf)
      
      RETURN
      
      END SUBROUTINE mddc_cdf_read_desc
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE mddc_cdf_read_mrf(this,iou,prefix)
!  Subroutine to do the appropriate netCDF read calls for a mddc_mrf
!  
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE (mddc_mrf), INTENT (inout)          :: this
      INTEGER, INTENT (in)                       :: iou
      CHARACTER (len=*), INTENT (in), OPTIONAL   :: prefix

!  this        mddc_mrf - this is what gets defined and written.
!  iou         i/o unit number of the netCDF file. 
!  prefix      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(len=32)            :: prefix_use
      INTEGER, DIMENSION(3) :: dimlens
      INTEGER, DIMENSION(2) :: dimlens2
      
      CHARACTER(len=80)                          :: code_name
      CHARACTER(len=80)                          :: code_version
      CHARACTER(len=80)                          :: date_run
      CHARACTER(len=80)                          :: field_coils_id
      INTEGER                             :: n_field_cg    
      REAL(rprec), DIMENSION(:), ALLOCATABLE     :: rdiag_coilg_1
      REAL(rprec), DIMENSION(:), ALLOCATABLE     :: extcur_mg
      INTEGER                             :: ir    
      INTEGER                             :: jz    
      INTEGER                             :: kp    
      INTEGER                             :: kp_store    
      REAL(rprec)                                :: rmin
      REAL(rprec)                                :: rmax       
      REAL(rprec)                                :: zmin
      REAL(rprec)                                :: zmax
      INTEGER                             :: n_field_periods
      LOGICAL                                    :: lstell_sym
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: a_r
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: a_f
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: a_z

      LOGICAL                                    :: use_con_shell
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE   :: a_s_r
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE   :: a_s_f
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE   :: a_s_z
      INTEGER                                    :: kp_shell
      INTEGER                                    :: kp_shell_store

! Declare local variables
      INTEGER           :: ir1, ir2, ir3, if1, if2, if3, iz1,           &
     &    iz2, iz3
      INTEGER           :: ier1, ier2, ier3
      INTEGER                             :: n_field_cg_decl    

      CHARACTER(len=*), PARAMETER  :: sub_name =                               &
     &  'mddc_cdf_read_data: '

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
      CALL mddc_cdf_defvn_mrf(prefix_use)
         
! Read Components
! Note: Read in to variables local to this subroutine.
! Arrays require inquiry regarding size, and allocation before actual reading.

      CALL cdf_read(iou, TRIM(vn_code_name_use), code_name)
      CALL cdf_read(iou, TRIM(vn_code_version_use),code_version)
      CALL cdf_read(iou, TRIM(vn_date_run_use), date_run)
      CALL cdf_read(iou, TRIM(vn_field_coils_id_use),field_coils_id)
      CALL cdf_read(iou, TRIM(vn_n_field_cg_use),n_field_cg)

      CALL cdf_inquire(iou, TRIM(vn_rdiag_coilg_1_use),dimlens)
      n_field_cg_decl = dimlens(1)
      CALL assert_eq(0,dimlens(2),dimlens(3),                                  &
     &   sub_name // 'Unexpected rdiag_coilg_1 dimensions')
      CALL assert_eq(n_field_cg,n_field_cg_decl,                               &
     &   sub_name // 'Disagreement rdiag_coilg_1 dimensions')
      ALLOCATE(rdiag_coilg_1(n_field_cg),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc rdiag_coilg_1')
      CALL cdf_read(iou, TRIM(vn_rdiag_coilg_1_use), rdiag_coilg_1)

      CALL cdf_inquire(iou, TRIM(vn_extcur_mg_use),dimlens)
      n_field_cg_decl = dimlens(1)
      CALL assert_eq(0,dimlens(2),dimlens(3),                                  &
     &   sub_name // 'Unexpected extcur_mg dimensions')
      CALL assert_eq(n_field_cg,n_field_cg_decl,                               &
     &   sub_name // 'Disagreement extcur_mg dimensions')
      ALLOCATE(extcur_mg(n_field_cg),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc extcur_mg')
      CALL cdf_read(iou, TRIM(vn_extcur_mg_use), extcur_mg)

      CALL cdf_read(iou, TRIM(vn_ir_use), ir)
      CALL cdf_read(iou, TRIM(vn_jz_use), jz)
      CALL cdf_read(iou, TRIM(vn_kp_use), kp)
      CALL cdf_read(iou, TRIM(vn_kp_store_use), kp_store)
      CALL cdf_read(iou, TRIM(vn_rmin_use), rmin)
      CALL cdf_read(iou, TRIM(vn_rmax_use), rmax)
      CALL cdf_read(iou, TRIM(vn_zmin_use), zmin)
      CALL cdf_read(iou, TRIM(vn_zmax_use), zmax)
      CALL cdf_read(iou, TRIM(vn_n_field_periods_use), n_field_periods)
      CALL cdf_read(iou, TRIM(vn_lstell_sym_use), lstell_sym)

      CALL cdf_inquire(iou, TRIM(vn_a_r_use),dimlens)
      ir1 = dimlens(1)
      ir2 = dimlens(2)
      ir3 = dimlens(3)
      CALL assert_eq(ir,ir1,                                                   &
     &   sub_name // 'Disagreement ir a_r dimensions')
      CALL assert_eq(jz,ir2,                                                   &
     &   sub_name // 'Disagreement jz a_r dimensions')
      CALL assert_eq(kp_store,ir3,                                             &
     &   sub_name // 'Disagreement kp_store a_r dimensions')
      ALLOCATE(a_r(ir1,ir2,ir3),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc a_r')
      CALL cdf_read(iou, TRIM(vn_a_r_use), a_r)

      CALL cdf_inquire(iou, TRIM(vn_a_f_use),dimlens)
      if1 = dimlens(1)
      if2 = dimlens(2)
      if3 = dimlens(3)
      CALL assert_eq(ir,if1,                                                   &
     &   sub_name // 'Disagreement ir a_f dimensions')
      CALL assert_eq(jz,if2,                                                   &
     &   sub_name // 'Disagreement jz a_f dimensions')
      CALL assert_eq(kp_store,if3,                                             &
     &   sub_name // 'Disagreement kp_store a_f dimensions')
      ALLOCATE(a_f(if1,if2,if3),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc a_f')
      CALL cdf_read(iou, TRIM(vn_a_f_use), a_f)

      CALL cdf_inquire(iou, TRIM(vn_a_z_use),dimlens)
      iz1 = dimlens(1)
      iz2 = dimlens(2)
      iz3 = dimlens(3)
      CALL assert_eq(ir,iz1,                                                   &
     &   sub_name // 'Disagreement ir a_z dimensions')
      CALL assert_eq(jz,iz2,                                                   &
     &   sub_name // 'Disagreement jz a_z dimensions')
      CALL assert_eq(kp_store,iz3,                                             &
     &   sub_name // 'Disagreement kp_store a_z dimensions')
      ALLOCATE(a_z(iz1,iz2,iz3),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc a_z')
      CALL cdf_read(iou, TRIM(vn_a_z_use), a_z)

! Read
      SELECT CASE (code_version)

         CASE ('MRC 2014-09-28')
            CALL cdf_read(iou, TRIM(vn_use_con_shell_use),                     &
     &                    use_con_shell)

            IF (use_con_shell) THEN
               CALL cdf_read(iou, TRIM(vn_kp_shell_use), kp_shell)

               CALL cdf_inquire(iou, TRIM(vn_a_r_use), dimlens2)
               ir1 = dimlens2(1)
               ir2 = dimlens2(2)
               ALLOCATE(a_s_r(ir1,ir2))
               CALL cdf_read(iou, TRIM(vn_a_s_r_use), a_s_r)

               CALL cdf_inquire(iou, TRIM(vn_a_f_use), dimlens2)
               if1 = dimlens2(1)
               if2 = dimlens2(2)
               ALLOCATE(a_s_r(if1,if2))
               CALL cdf_read(iou, TRIM(vn_a_s_f_use), a_s_f)

               CALL cdf_inquire(iou, TRIM(vn_a_z_use), dimlens2)
               iz1 = dimlens2(1)
               iz2 = dimlens2(2)
               ALLOCATE(a_s_r(iz1,iz2))
               CALL cdf_read(iou, TRIM(vn_a_s_z_use), a_s_z)
            END IF

         CASE DEFAULT
            use_con_shell = .false.

      END SELECT

! Create the mddc_mrf, this
      CALL mddc_mrf_construct(this,code_name,code_version,                     &
     &  date_run,field_coils_id,rdiag_coilg_1,extcur_mg,kp,                    &
     &  rmin,rmax,zmin,zmax,n_field_periods,lstell_sym,a_r,a_f,a_z,            &
     &  use_con_shell, a_s_r, a_s_f, a_s_z, kp_shell)

! Deallocate the local pointers
      DEALLOCATE(rdiag_coilg_1,extcur_mg,STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'dealloc 1')
       DEALLOCATE(a_r,a_f,a_z,STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'dealloc 2')

      IF (ALLOCATED(a_s_r)) DEALLOCATE(a_s_r)
      IF (ALLOCATED(a_s_f)) DEALLOCATE(a_s_f)
      IF (ALLOCATED(a_s_z)) DEALLOCATE(a_s_z)

      RETURN
      
      END SUBROUTINE mddc_cdf_read_mrf

!*******************************************************************************
! SECTION IV. AUXILIARY FUNCTIONS AND SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE mddc_cdf_defvn_desc(prefix_use)
!  Subroutine to do define the character variable names for a mddc_desc,
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
     &  'mddc_cdf_defvn_desc: '

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

! Define all variable names
      vn_s_name_use = mddc_cdf_mknam(prefix_use,vn_s_name)
      vn_l_name_use = mddc_cdf_mknam(prefix_use,vn_l_name)
      vn_units_use = mddc_cdf_mknam(prefix_use,vn_units)
      vn_sigma_default_use = mddc_cdf_mknam(prefix_use,                        &
     &   vn_sigma_default)
      vn_l_mdcoil_def_use = mddc_cdf_mknam(prefix_use,                         &
     &   vn_l_mdcoil_def)
      vn_mddc_type_use = mddc_cdf_mknam(prefix_use,                            &
     &   vn_mddc_type)
      vn_flux_factor_use = mddc_cdf_mknam(prefix_use,                          &
     &   vn_flux_factor)
      
      RETURN
      
      END SUBROUTINE mddc_cdf_defvn_desc
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE mddc_cdf_defvn_mrf(prefix_use)
!  Subroutine to do define the character variable names for a mddc_mrf,
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
     &  'mddc_cdf_defvn_mrf: '

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

! Define all variable names
      vn_code_name_use = mddc_cdf_mknam(prefix_use,vn_code_name)
      vn_code_version_use = mddc_cdf_mknam(prefix_use,                         &
     &   vn_code_version)
      vn_date_run_use = mddc_cdf_mknam(prefix_use,vn_date_run)
      vn_field_coils_id_use = mddc_cdf_mknam(prefix_use,                       &
     &   vn_field_coils_id)
      vn_n_field_cg_use = mddc_cdf_mknam(prefix_use,                           &
     &   vn_n_field_cg)
      vn_rdiag_coilg_1_use = mddc_cdf_mknam(prefix_use,                        &
     &   vn_rdiag_coilg_1)
      vn_extcur_mg_use = mddc_cdf_mknam(prefix_use,vn_extcur_mg)
      vn_ir_use = mddc_cdf_mknam(prefix_use,vn_ir)
      vn_jz_use = mddc_cdf_mknam(prefix_use,vn_jz)
      vn_kp_use = mddc_cdf_mknam(prefix_use,vn_kp)
      vn_kp_store_use = mddc_cdf_mknam(prefix_use,vn_kp_store)
      vn_rmin_use = mddc_cdf_mknam(prefix_use,vn_rmin)
      vn_rmax_use = mddc_cdf_mknam(prefix_use,vn_rmax)
      vn_zmin_use = mddc_cdf_mknam(prefix_use,vn_zmin)
      vn_zmax_use = mddc_cdf_mknam(prefix_use,vn_zmax)
      vn_n_field_periods_use = mddc_cdf_mknam(prefix_use,                      &
     &   vn_n_field_periods)
      vn_lstell_sym_use = mddc_cdf_mknam(prefix_use,vn_lstell_sym)
      vn_a_r_use = mddc_cdf_mknam(prefix_use,vn_a_r)
      vn_a_f_use = mddc_cdf_mknam(prefix_use,vn_a_f)
      vn_a_z_use = mddc_cdf_mknam(prefix_use,vn_a_z)
      vn_use_con_shell_use = mddc_cdf_mknam(prefix_use,vn_use_con_shell)
      vn_a_s_r_use = mddc_cdf_mknam(prefix_use,vn_a_s_r)
      vn_a_s_f_use = mddc_cdf_mknam(prefix_use,vn_a_s_f)
      vn_a_s_z_use = mddc_cdf_mknam(prefix_use,vn_a_s_z)
      vn_kp_shell_use = mddc_cdf_mknam(prefix_use,vn_kp_shell)
      vn_kp_shell_store_use = mddc_cdf_mknam(prefix_use,                       &
     &                                       vn_kp_shell_store)

      RETURN

      END SUBROUTINE mddc_cdf_defvn_mrf
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      FUNCTION mddc_cdf_mknam(c1,c2)
! A simple function to help in the generation of names

!-----------------------------------------------
!   F u n c t i o n   N a m e
!-----------------------------------------------
      CHARACTER(LEN=64) mddc_cdf_mknam

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER(LEN=*), INTENT (in) :: c1,c2

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
      IF (LEN_TRIM(c1) .eq. 0) THEN
         mddc_cdf_mknam = TRIM(c2)
      ELSE
         mddc_cdf_mknam = ADJUSTL(TRIM(c1) // '_' // TRIM(c2))
      ENDIF
      
      RETURN
       
      END FUNCTION mddc_cdf_mknam

!-----------------------------------------------
!-----------------------------------------------
!*******************************************************************************
! SECTION X.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!
! JDH 2007-06-11. First version of mddc_cdf.
!    Based on diagnostic_cdf
!
!  JDH 2008-01-19 - Added => null() to pointer declarations
!
!------------------ (Below comments for diagnostic_cdf) ------------------------
!  09-02-2004 JDH - Initial Coding, based on bsc_cdf
!
!  JDH 12-11-2004
!     Modified _desc subroutines, to correspond with changes in diagnostic_T.
!
!  JDH 07-01-2005
!     Added bsc_destroy of local bsc_coil mdcoil in subroutine 
!     diagnostic_cdf_read_desc, to avoid memory leaks.
!       
!  JDH 2009-06-15
!    Eliminated mddc_data derived type - not needed.


      END MODULE mddc_cdf
