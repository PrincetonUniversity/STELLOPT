!*******************************************************************************
!  File bsc_cdf.f
!  Contains the module bsc_cdf.f
!    Module for defining variables and writing netCDF files with the
!    derived types bsc_coil and bsc_coilcoll, from the bsc (Biot-Savart Coil) module
!
!-------------------------------------------------------------------------------
!   DEPENDENCIES
!-------------------------------------------------------------------------------
!
!    This module uses the following modules:
!       stel_kinds
!       bsc
!       ezcdf
!
!-------------------------------------------------------------------------------
!   CHANGE HISTORY
!-------------------------------------------------------------------------------
!
!  12.13.2002 - Initial Coding - Ed Lazarus, 
!  12.16.2002 - JDH Initial Comments, limit to bsc_cdf subroutines
!  12.17.2002 - JDH - return to using stel_kinds. Eliminate some unused variables.
!  12.18.2002 - JDH - Eliminated identifier. Added prefix. Made _coilcoll routines
!     call the _coil routines.
!  09.11.2003 - JDH - Added 'fil_rogo'wski information
!  09-27-2004 - JDH - Added bsc_cdf_read_coil subroutine. Modified coding to be more
!     consistent with structure of diagnostic_cdf and signal_cdf.
!
!-------------------------------------------------------------------------------
!   USAGE
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!   COMMENTS
!-------------------------------------------------------------------------------
!
!*******************************************************************************

!*******************************************************************************
!  MODULE bsc_cdf
!    
! SECTION I.   VARIABLE DECLARATIONS
! SECTION II.  INTERFACE BLOCKS
! SECTION III. DEFINITION SUBROUTINES
! SECTION IV.  WRITING SUBROUTINES
! SECTION V.   READING SUBROUTINES
! SECTION VI.  AUXILLIARY FUNCTIONS
!*******************************************************************************

      MODULE bsc_cdf

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

!-------------------------------------------------------------------------------
!  Type declarations - lengths of reals, integers, and complexes.
!-------------------------------------------------------------------------------

      USE stel_kinds
      USE stel_constants
      USE bsc_T 
      USE ezcdf
      USE v3_utilities

!-------------------------------------------------------------------------------
!  Implicit None comes after USE statements, before other declarations
!-------------------------------------------------------------------------------
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variable Names for netCDF
!-------------------------------------------------------------------------------

      CHARACTER (LEN=*),  PARAMETER :: 
     &  vn_c_type = 'c_type',                                                  &
     &  vn_s_name = 's_name',                                                  &
     &  vn_l_name = 'l_name',                                                  &
     &  vn_current = 'current',                                                &         
     &  vn_raux = 'raux',                                                      &            
     &  vn_xnod = 'xnod',                                                      &            
     &  vn_ehnod = 'ehnod',                                                    &    
     &  vn_rcirc = 'rcirc',                                                    &        
     &  vn_xcent = 'xcent',                                                    &        
     &  vn_enhat = 'enhat',                                                    &
     &  vn_ave_n_area = 'ave_n_area'

      CHARACTER (LEN=64), PRIVATE :: 
     &  vn_c_type_use,                                                         &
     &  vn_s_name_use,                                                         &
     &  vn_l_name_use,                                                         &
     &  vn_current_use,                                                        &         
     &  vn_raux_use,                                                           &            
     &  vn_xnod_use,                                                           &            
     &  vn_ehnod_use,                                                          &    
     &  vn_rcirc_use,                                                          &        
     &  vn_xcent_use,                                                          &        
     &  vn_enhat_use,                                                          &
     &  vn_ave_n_area_use

!*******************************************************************************
! SECTION II. INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------

      CONTAINS
!*******************************************************************************
! SECTION III. DEFINITION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_cdf_define_coil(this,lunit,prefix)
!  Subroutine to do the appropriate netCDF definition calls for a bsc_coil

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE (bsc_coil), INTENT (in) :: this
      INTEGER                      :: lunit
      CHARACTER (len=*)            :: prefix

!  this        bsc_coil - this is the coils that gets defined.
!  lunit       i/o unit number
!  prefix      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(len=32) :: prefix_use
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

! Define the prefix to actually use
      prefix_use = TRIM(ADJUSTL(prefix))

! Define all vn_--_use variable names
      CALL bsc_cdf_defvn_coil(prefix_use)
         
! Define Components common to all c_types
      CALL cdf_define(lunit, TRIM(vn_c_type_use), this%c_type)
      CALL cdf_define(lunit, TRIM(vn_s_name_use), this%s_name)
      CALL cdf_define(lunit, TRIM(vn_l_name_use), this%l_name)
      CALL cdf_define(lunit, TRIM(vn_current_use), this%current)
      CALL cdf_define(lunit, TRIM(vn_raux_use), this%raux)

! Particular coding, depending on c_type

      SELECT CASE (this%c_type)
      
      CASE ('fil_loop','floop') ! Filamentary Loop Variables
         IF (ASSOCIATED(this%xnod)) THEN
            CALL cdf_define(lunit, TRIM(vn_xnod_use), this%xnod)
         END IF ! this%xnod ASSOCIATED

      CASE ('fil_circ', 'fcirc') ! Filamentary Circle Variables
         CALL cdf_define(lunit, TRIM(vn_rcirc_use), this%rcirc)
         CALL cdf_define(lunit, TRIM(vn_xcent_use), this%xcent(1:3))
         CALL cdf_define(lunit, TRIM(vn_enhat_use), this%enhat(1:3))

      CASE ('fil_rogo') ! Rogowskis
         IF (ASSOCIATED(this%xnod)) THEN
            CALL cdf_define(lunit, TRIM(vn_xnod_use), this%xnod)
         END IF ! this%xnod ASSOCIATED
         CALL cdf_define(lunit, TRIM(vn_ave_n_area_use),                       &
     &         this%ave_n_area)
      
      END SELECT
      
      END SUBROUTINE bsc_cdf_define_coil

!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_cdf_define_coilcoll(this,lunit)
!  Subroutine to  do the appropriate netCDF definition calls for a bsc_coilcoll
!  To avoid duplicate names in the netCDF files, the variable names will
!  have a prefix added on.

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE (bsc_coilcoll), INTENT (in) :: this
      INTEGER                          :: lunit

!  this        bsc_coilcoll - this is the coilcoll that gets defined.
!  lunit       i/o unit number
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, n, ncoild
      INTEGER dimlens(2)
      CHARACTER(len=40) nowname
      CHARACTER(len=40) :: prefix

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
!  Not sure of the reason for the next IF test. JDH
!  
      IF (this%s_name .eq. ' ') THEN
         WRITE(*,*) 'this%s_name = one blank. bsc_cdf_define_coilcoll'
         WRITE(*,*) ' is returning'
         RETURN
      END IF
      
      ncoild = this%ncoil

!  Next loop could be augmented to make sure that the prefixes are unique.

      DO i = 1,ncoild  ! Loop over coils in the coilcoll
         prefix = this%coils(i)%s_name
         CALL bsc_cdf_define_coil(this%coils(i),lunit,prefix)
      END DO

      END SUBROUTINE bsc_cdf_define_coilcoll

!*******************************************************************************
! SECTION IV. WRITING SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_cdf_write_coil(this,lunit,prefix)
!  Subroutine to  do the appropriate netCDF definition calls for a bsc_coil

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE (bsc_coil), INTENT (in) :: this
      INTEGER                      :: lunit
      CHARACTER (len=*)            :: prefix

!  this        bsc_coil - this is the coil that gets written.
!  lunit       i/o unit number
!  prefix      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(len=32) :: prefix_use
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

! Define the prefix to actually use
      prefix_use = TRIM(ADJUSTL(prefix))

! Define all vn_--_use variable names
      CALL bsc_cdf_defvn_coil(prefix_use)
         
! Write Components common to all c_types
      CALL cdf_write(lunit, TRIM(vn_c_type_use), this%c_type)
      CALL cdf_write(lunit, TRIM(vn_s_name_use), this%s_name)
      CALL cdf_write(lunit, TRIM(vn_l_name_use), this%l_name)
      CALL cdf_write(lunit, TRIM(vn_current_use), this%current)
      CALL cdf_write(lunit, TRIM(vn_raux_use), this%raux)

! Particular coding, depending on c_type

      SELECT CASE (this%c_type)
      
      CASE ('fil_loop','floop') ! Filamentary Loop Variables
         IF (ASSOCIATED(this%xnod)) THEN
            CALL cdf_write(lunit, TRIM(vn_xnod_use), this%xnod)
         END IF ! this%xnod ASSOCIATED

      CASE ('fil_circ', 'fcirc') ! Filamentary Circle Variables
         CALL cdf_write(lunit, TRIM(vn_rcirc_use), this%rcirc)
         CALL cdf_write(lunit, TRIM(vn_xcent_use), this%xcent(1:3))
         CALL cdf_write(lunit, TRIM(vn_enhat_use), this%enhat(1:3))

      CASE ('fil_rogo') ! Rogowskis
         IF (ASSOCIATED(this%xnod)) THEN
            CALL cdf_write(lunit, TRIM(vn_xnod_use), this%xnod)
         END IF ! this%xnod ASSOCIATED
         CALL cdf_write(lunit, TRIM(vn_ave_n_area_use),                        &
     &       this%ave_n_area)      
      
      END SELECT
      
      END SUBROUTINE bsc_cdf_write_coil

!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_cdf_write_coilcoll(this,lunit)
!  Subroutine to  do the appropriate netCDF definition calls for a bsc_coilcoll
!  To avoid duplicate names in the netCDF files, the variable names will
!  have a prefix added on.

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE (bsc_coilcoll), INTENT (in) :: this
      INTEGER                          :: lunit

!  this        bsc_coilcoll - this is the coilcoll that gets written.
!  lunit       i/o unit number
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, n, ncoild
      INTEGER dimlens(2)
      CHARACTER(LEN=40) nowname
      CHARACTER (len=40) :: prefix

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
!  Not sure of the reason for the next IF test. JDH
!  
      IF (this%s_name .eq. ' ') THEN
         WRITE(*,*) 'this%s_name = one blank. bsc_cdf_write_coilcoll'
         WRITE(*,*) ' is returning'
         RETURN
      END IF
      
      ncoild = this%ncoil

!  Next loop could be augmented to make sure that the prefixes are unique.

      DO i = 1,ncoild  ! Loop over coils in the coilcoll
         prefix = this%coils(i)%s_name
         CALL bsc_cdf_write_coil(this%coils(i),lunit,prefix)
      END DO

      END SUBROUTINE bsc_cdf_write_coilcoll

!*******************************************************************************
! SECTION V. READING SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_cdf_read_coil(this,iou,prefix)
!  Subroutine to do the appropriate netCDF read calls for a bsc_coil
!  

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE (bsc_coil), INTENT (inout)            :: this
      INTEGER, INTENT (in)                       :: iou
      CHARACTER (len=*), INTENT (in), OPTIONAL   :: prefix

!  this        bsc_coil - this is what gets defined and written.
!  iou         i/o unit number of the netCDF file. 
!  prefix      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(len=*), PARAMETER  :: sub_name =                               &
     &  'bsc_cdf_read_coil: '
      CHARACTER(len=32)            :: prefix_use
      INTEGER, DIMENSION(3) :: dimlens
      INTEGER               :: ier1, n2
      
      CHARACTER (len=8) :: c_type
      CHARACTER (len=30) :: s_name                                 
      CHARACTER (len=80) :: l_name
      REAL(rprec) :: eps_sq
      REAL(rprec) :: current
      REAL(rprec) :: raux
      REAL(rprec) :: rcirc
      REAL(rprec) :: ave_n_area
      REAL(rprec), DIMENSION(3) :: xcent, enhat
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: xnod 

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
      CALL bsc_cdf_defvn_coil(prefix_use)
         
! Read Components
! Note: Read in to variables local to this subroutine.
! Arrays require inquiry regarding size, and allocation before actual reading.

      CALL cdf_read(iou, TRIM(vn_c_type_use), c_type)
      CALL cdf_read(iou, TRIM(vn_s_name_use),s_name)
      CALL cdf_read(iou, TRIM(vn_l_name_use),l_name)
      CALL cdf_read(iou, TRIM(vn_current_use),current)      
      CALL cdf_read(iou, TRIM(vn_raux_use),raux)

      SELECT CASE (TRIM(c_type))
      
      CASE ('fil_loop','floop') ! Filamentary Loop Variables
         CALL cdf_inquire(iou, TRIM(vn_xnod_use),dimlens)
         ALLOCATE(xnod(dimlens(1),dimlens(2)),STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'alloc xnod')
         CALL assert_eq(3,dimlens(1),sub_name // 'bad xnod dim')
         CALL cdf_read(iou, TRIM(vn_xnod_use), xnod)
! Take into account logic in bsc_construct for loop closure
         IF (dimlens(2) .ge. 4) THEN
            n2 = dimlens(2) - 1
         ELSE
            n2 = dimlens(2)
         ENDIF
         CALL bsc_construct(this,c_type,s_name,l_name,current,                 &
     &   xnod(1:3,1:n2),raux=raux)

      CASE ('fil_circ', 'fcirc') ! Filamentary Circle Variables
         CALL cdf_read(iou, TRIM(vn_rcirc_use),rcirc)
         CALL cdf_read(iou, TRIM(vn_xcent_use),xcent)
         CALL cdf_read(iou, TRIM(vn_enhat_use),enhat)
         CALL bsc_construct(this,c_type,s_name,l_name,current,                 &
     &   rcirc = rcirc,xcent = xcent,enhat = enhat,raux = raux)

      CASE ('fil_rogo') ! Rogowskis
         CALL cdf_inquire(iou, TRIM(vn_xnod_use),dimlens)
         ALLOCATE(xnod(dimlens(1),dimlens(2)),STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'alloc xnod')
         CALL cdf_read(iou, TRIM(vn_xnod_use), xnod)
         CALL cdf_read(iou, TRIM(vn_ave_n_area_use),ave_n_area)
         CALL bsc_construct(this,c_type,s_name,l_name,current,                 &
     &      xnod,raux = raux,anturns = one,xsarea = ave_n_area)
      
      END SELECT

! Deallocate the local allocatable space
      IF (ALLOCATED(xnod)) THEN
         DEALLOCATE(xnod,STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'dealloc xnod')
      END IF
      
      RETURN
      
      END SUBROUTINE bsc_cdf_read_coil

!*******************************************************************************
! SECTION VI. AUXILLIARY FUNCTIONS
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_cdf_defvn_coil(prefix_use)
!  Subroutine to do define the character variable names for a bsc_coil,
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
     &  'bsc_cdf_defvn_coil: '

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

! Define all variable names
      vn_c_type_use = bsc_cdf_mknam(prefix_use,vn_c_type)
      vn_s_name_use = bsc_cdf_mknam(prefix_use,vn_s_name)
      vn_l_name_use = bsc_cdf_mknam(prefix_use,vn_l_name)
      vn_current_use = bsc_cdf_mknam(prefix_use,vn_current)
      vn_raux_use = bsc_cdf_mknam(prefix_use,vn_raux)
      vn_xnod_use = bsc_cdf_mknam(prefix_use,vn_xnod)
      vn_rcirc_use = bsc_cdf_mknam(prefix_use,vn_rcirc)
      vn_xcent_use = bsc_cdf_mknam(prefix_use,vn_xcent)
      vn_enhat_use = bsc_cdf_mknam(prefix_use,vn_enhat)
      vn_ave_n_area_use = bsc_cdf_mknam(prefix_use,vn_ave_n_area)
      
      RETURN
      
      END SUBROUTINE bsc_cdf_defvn_coil
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      
      FUNCTION bsc_cdf_mknam(c1,c2)
!  A simple function to help in the generation of names

!-----------------------------------------------
!   F u n c t i o n   N a m e
!-----------------------------------------------
      CHARACTER(LEN=40) bsc_cdf_mknam

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER(LEN=*), INTENT (in) :: c1,c2

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
      IF (LEN_TRIM(c1) .eq. 0) THEN
         bsc_cdf_mknam = TRIM(c2)
      ELSE
         bsc_cdf_mknam = ADJUSTL(TRIM(c1) // '_' // TRIM(c2))
      ENDIF

      RETURN
       
      END FUNCTION bsc_cdf_mknam

!-----------------------------------------------
!-----------------------------------------------

      END MODULE bsc_cdf
