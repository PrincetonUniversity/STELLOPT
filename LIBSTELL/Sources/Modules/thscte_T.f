
!*******************************************************************************
!  File thscte_T.f
!  Contains module thscte_T
!  Defines derived-types: thscte_desc
!  A type of Diagnostic - Thomspon Scattering Te diagnostic
!
!*******************************************************************************
!  MODULE sxr_T
!    (SXR Type Definition, for the V3FIT code)
! SECTION I.    VARIABLE DECLARATIONS
! SECTION II.   DERIVED-TYPE DECLARATIONS
! SECTION III.  INTERFACE BLOCKS
! SECTION IV.   CONSTRUCTION SUBROUTINES
! SECTION V.    DESTRUCTION SUBROUTINES
! SECTION VI.   ASSIGNMENT SUBROUTINES
! SECTION VII.  OUTPUT SUBROUTINES
! SECTION VIII. PRIVATE ROUTINES USED IN SXR_T
!
! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
      MODULE thscte_T

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

!-------------------------------------------------------------------------------
!  Type declarations - lengths of reals, integers, and complexes.
!  Frequently used mathematical constants, lots of extra precision.
!-------------------------------------------------------------------------------

      USE stel_kinds , only : rprec
      USE stel_constants, only : pi, zero
      USE safe_open_mod    !from LIBSTELL/MODULES
     
!-------------------------------------------------------------------------------
!  Use Statements for other structures, V3 Utilities
!-------------------------------------------------------------------------------
!     USE v3_utilities

!-------------------------------------------------------------------------------
!  Implicit None comes after USE statements, before other declarations
!-------------------------------------------------------------------------------
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Make type declarations and constants Private, so there are no conflicts.
!-------------------------------------------------------------------------------
      PRIVATE rprec, pi, zero

!-------------------------------------------------------------------------------
!  Lengths of Character Variables
!-------------------------------------------------------------------------------

!------------------------------------------------------------------------------ 
      INTEGER,PARAMETER  :: chord_name_len=30           
!------------------------------------------------------------------------------      

!*******************************************************************************
! SECTION II. DERIVED-TYPE DECLARATIONS
!     thscte Description:
!       thscte_desc  
!     Type of diagnostic specified by  % d_type = 'thscte'.
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!  Declare type thscte_desc
!         chord_name            - character, chord name
!         xcart(3)              - Cartesian position vector - where Te is measured
!         thsc_type             - 't' - Temperature
!                               - 'd' - Density
!                               - 'p' - Pressure
!-------------------------------------------------------------------------------
      TYPE thscte_desc
        CHARACTER(LEN=chord_name_len) :: chord_name
        REAL(rprec), DIMENSION(3) :: xcart
        CHARACTER(LEN=1) :: thsc_type
      END TYPE thscte_desc

!*******************************************************************************
! SECTION III. INTERFACE BLOCKS
!*******************************************************************************

      CONTAINS
!*******************************************************************************
! SECTION IV. CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a thscte_desc
!
!  For d_type = 'thscte' (Thomson Scattering Te)
!-------------------------------------------------------------------------------
      SUBROUTINE thscte_desc_construct(this, chord_name, xcart,                &
     &                                 chord_type)

         IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Argument Declarations
!-------------------------------------------------------------------------------
         TYPE (thscte_desc), INTENT(inout)        :: this
         CHARACTER(LEN=chord_name_len),INTENT(in) :: chord_name
         REAL(rprec), DIMENSION(3), INTENT(in)    :: xcart
         CHARACTER(LEN=1),INTENT(in)              :: chord_type
 
!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------

!  Assignments    
         this % chord_name = chord_name
         this % xcart = xcart
         this % thsc_type = chord_type
      
      END SUBROUTINE thscte_desc_construct

!*******************************************************************************
! SECTION V. DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Destroy an thscte_desc
!
!  ARGUMENT
!  this     - an thscte_desc
!-------------------------------------------------------------------------------
      SUBROUTINE thscte_desc_destroy(this)
      
         TYPE (thscte_desc),INTENT(inout) :: this
         
         this % chord_name = ''
         this % xcart = zero

      
      END SUBROUTINE thscte_desc_destroy


!*******************************************************************************
! SECTION VI.   ASSIGNMENT SUBROUTINES
!
!   These are not needed because the intrinsic assignments work
!*******************************************************************************

!*******************************************************************************
! SECTION VII.  OUTPUT SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Write out the contents  of a thscte_desc
!  if iou and filaname are present - write to file
!  if iou and filename are not present - write to stdout (screen)
!
! THIS NEEDS MODIFYING TO BE ABLE TO APPEND RECORDS AND NOT OVERWRITE FILES
!-------------------------------------------------------------------------------

      SUBROUTINE thscte_desc_write(this,iounit,filename)
      IMPLICIT NONE
!-------------------------------------------------------------------------------
! Arguments
! this     - sxr_chord
! iou      - output io unit number
! filename - output file name
!-------------------------------------------------------------------------------     

      TYPE (thscte_desc),INTENT(in)    :: this
      INTEGER, OPTIONAL,INTENT(in)     :: iounit
      CHARACTER*300,OPTIONAL,INTENT(in)  :: filename
!-------------------------------------------------------------------------------
! Local Variables
! iou  - iounit to use
! istat - status of file opening
!-------------------------------------------------------------------------------     
     
      INTEGER :: iou = 6
      INTEGER :: istat = 0      !status of safe_open call
      
      IF (PRESENT(iounit).AND.PRESENT(filename)) THEN
         iou=iounit
         CALL safe_open(iou,istat,filename,'replace','formatted')
         WRITE(iou,*) 'chord name  - ', this % chord_name 
         WRITE(iou,*) 'position -', this % xcart
      ELSE
         WRITE(*,*)'chord name  - ',this % chord_name 
         WRITE(*,*)'position -',this % xcart
      END IF
      
      END SUBROUTINE thscte_desc_write
!*******************************************************************************
! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!
!  TO DO 2011-10-23
! 1)    Change name of chord_name component to something more appropriate
!
!  2011-10-23 JDH
!    First version of module, based on sxrch_T
!
!
      END MODULE thscte_T
