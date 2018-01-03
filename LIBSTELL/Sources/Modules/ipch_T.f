
!*******************************************************************************
!  File ipch_T.f
!  Contains module ipch_T
!  Defines derived-types: ipch_desc
!  A type of Diagnostic - Interferometry-Polarimetry chordal diagnostic
!
!*******************************************************************************
!  MODULE ip_T
!    (ip Type Definition, for the V3FIT code)
! SECTION I.    VARIABLE DECLARATIONS
! SECTION II.   DERIVED-TYPE DECLARATIONS
! SECTION III.  INTERFACE BLOCKS
! SECTION IV.   CONSTRUCTION SUBROUTINES
! SECTION V.    DESTRUCTION SUBROUTINES
! SECTION VI.   ASSIGNMENT SUBROUTINES
! SECTION VII.  OUTPUT SUBROUTINES
! SECTION VIII. PRIVATE ROUTINES USED IN ip_T
!
! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
      MODULE ipch_T

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
      USE integration_path
     
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
!     ipCH Description:
!       ipch_desc  
!     Type of diagnostic specified by  % d_type = 'ipch'.
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!  Declare type ipch_desc
!         chord_name            - character, chord name
!         xcart_i(3)            - Cartesian position vector, start of chord (meters)
!         xcart_f(3)            - Cartesian position vector, end of chord (meters)
!         ip_type               - 'i' - interferometry
!                               - 'p' - polarimetry
!         wavelength            _ wavelength of the beam
!-------------------------------------------------------------------------------
      TYPE ipch_desc
        CHARACTER(LEN=chord_name_len) :: chord_name
        CHARACTER(LEN=1)              :: ip_type
        REAL(rprec)                   :: wavelength
        LOGICAL                       :: inDegrees
        TYPE(vertex), POINTER         :: chordPath => null()
      END TYPE ipch_desc

!*******************************************************************************
! SECTION III. INTERFACE BLOCKS
!*******************************************************************************

      CONTAINS
!*******************************************************************************
! SECTION IV. CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a ipch_desc
!
!  For d_type = 'ipch' (interferometry-polarimetry chord)
!-------------------------------------------------------------------------------
      SUBROUTINE ipch_desc_construct(this, chord_name,                         &
     &                               xcart_i, xcart_f,                         &
     &                               ip_type, wavelength,                      &
     &                               inDegrees)

         IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Argument Declarations
!-------------------------------------------------------------------------------
         TYPE (ipch_desc), INTENT(inout)          :: this
         CHARACTER(LEN=chord_name_len),INTENT(in) :: chord_name
         REAL(rprec), DIMENSION(3), INTENT(in)    :: xcart_i
         REAL(rprec), DIMENSION(3), INTENT(in)    :: xcart_f
         CHARACTER(LEN=1), INTENT(in)             :: ip_type
         REAL(rprec), INTENT(in)                  :: wavelength
         LOGICAL, INTENT(in)                      :: inDegrees
 
!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------

!  Assignments    
         this % chord_name = chord_name
         this % ip_type = ip_type
         this % wavelength = wavelength
         this % inDegrees = inDegrees

!  Must NULL out the vertices array or else it will point to the last
!  integration_path created in memory

         CALL path_append_vertex(this%chordPath, xcart_i)
         CALL path_append_vertex(this%chordPath, xcart_f)

      END SUBROUTINE ipch_desc_construct

!*******************************************************************************
! SECTION V. DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Destroy an ipch_desc
!
!  ARGUMENT
!  this     - an ipch_desc
!-------------------------------------------------------------------------------
      SUBROUTINE ipch_desc_destroy(this)
      
         TYPE (ipch_desc),INTENT(inout) :: this
         
         this % chord_name = ''
         this % ip_type = ' '
         this % wavelength = zero

         CALL path_destruct(this%chordPath)

      END SUBROUTINE ipch_desc_destroy


!*******************************************************************************
! SECTION VI.   ASSIGNMENT SUBROUTINES
!
!   These are not needed because the intrinsic assignments work
!*******************************************************************************

!*******************************************************************************
! SECTION VII.  OUTPUT SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Write out the contents  of a ipch_desc
!  if iou and filaname are present - write to file
!  if iou and filename are not present - write to stdout (screen)
!
! THIS NEEDS MODIFYING TO BE ABLE TO APPEND RECORDS AND NOT OVERWRITE FILES
!-------------------------------------------------------------------------------

      SUBROUTINE ipch_desc_write(this,iounit,filename)
      IMPLICIT NONE
!-------------------------------------------------------------------------------
! Arguments
! this     - ip_chord
! iou      - output io unit number
! filename - output file name
!-------------------------------------------------------------------------------     

      TYPE (ipch_desc),INTENT(in)    :: this
      INTEGER, OPTIONAL,INTENT(in)   :: iounit
      CHARACTER*300,OPTIONAL,INTENT(in)  :: filename
!-------------------------------------------------------------------------------
! Local Variables
! iou  - iounit to use
! istat - status of file opening
!-------------------------------------------------------------------------------     
      REAL(rprec), DIMENSION(3)    :: xcart_i
      REAL(rprec), DIMENSION(3)    :: xcart_f

      INTEGER :: iou = 6
      INTEGER :: istat = 0      !status of safe_open call

      xcart_i = this%chordPath%position
      xcart_f = this%chordPath%next%position

      IF (PRESENT(iounit).AND.PRESENT(filename)) THEN
         iou=iounit
         CALL safe_open(iou,istat,filename,'replace','formatted')
         WRITE(iou,*) 'chord name  - ', this % chord_name 
         WRITE(iou,*) 'start position -', xcart_i
         WRITE(iou,*) 'end position   -', xcart_f
         WRITE(iou,*) 'ip_type   -', this % ip_type
         WRITE(iou,*) 'wavelength -', this % wavelength
      ELSE
         WRITE(*,*) 'chord name  - ',this % chord_name 
         WRITE(*,*) 'start position -', xcart_i
         WRITE(*,*) 'end position   -', xcart_f
         WRITE(*,*) 'ip_type   -', this % ip_type
         WRITE(*,*) 'wavelength -', this % wavelength
      END IF
      
      END SUBROUTINE ipch_desc_write
!*******************************************************************************
! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!
!  GJH 2009-08-18. First version of sxr_T. Copied and edited from mddc_T
!
!  GJH 2010-01-22  Added measurement units to sxr_chords
!
!  JDH 2011-08-01
!    Refactor sxrc -> sxrch
!
!  JDH 2011-08-29 - Added these changes from GJH:
!  GJH 2010-09-08
!    Changed Ro,Zo,Phio to Ri,Zi,Phii 
!    removed camera_type
!    removed position units
!    added chord_num 
!             in sxrch_desc 
!
!  JDH 2011-09-06
!    Modified R2x in sxrch_desc_construct
!
!  2011-09-08 JDH
!   Significant modification and code elimination. Just ID, start and end
!   positions, and calibration constant. Creation from camera description
!   must now be done elsewhere.

!  2011-10-17 JDH
!    Further simplification. Now just start and end position (cartesian) a
!    and chord name (longer - 30 characters)

!  2012-03-15 JDH
!    First version for ipch. Copied and edited from sxrch, added ip_type
     
      END MODULE ipch_T
