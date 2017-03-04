!*******************************************************************************
!  File limiter_iso_T.f
!  Contains module limiter_iso_T
!  Defines derived-types: limiter_iso

!*******************************************************************************
!  MODULE limiter_iso_T
!    (limiter_iso Type Definition, for the V3FIT code)
! SECTION I.    VARIABLE DECLARATIONS
! SECTION II.   DERIVED-TYPE DECLARATIONS
! SECTION III.  INTERFACE BLOCKS
! SECTION IV.   CONSTRUCTION SUBROUTINES
! SECTION V.    DESTRUCTION SUBROUTINES
! SECTION VI.   ASSIGNMENT SUBROUTINES
! SECTION VII.  OUTPUT SUBROUTINES
! SECTION VIII. FUNCTION EVALUATION SUBROUTINES
! SECTION XII.  AUXILIARY SUBROUTINES
! SECTION XIII. DEBUGGING SUBROUTINES
! SECTION XV.   DUPLICATE CODING FOR TESTING
! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
      MODULE limiter_iso_T

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

!-------------------------------------------------------------------------------
!  Type declarations - lengths of reals, integers, and complexes.
!  Frequently used mathematical constants, lots of extra precision.
!-------------------------------------------------------------------------------
      USE stel_kinds, only : rprec, cprec
      USE stel_constants, only : pi, twopi, one, zero
     
!-------------------------------------------------------------------------------
!  Use Statements for V3 Utilities
!-------------------------------------------------------------------------------
      USE v3_utilities

!-------------------------------------------------------------------------------
!  Implicit None comes after USE statements, before other declarations
!-------------------------------------------------------------------------------
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Make type declarations and constants Private, so there are no conflicts.
!-------------------------------------------------------------------------------
      PRIVATE rprec, cprec, pi, twopi, one, zero


!*******************************************************************************
! SECTION II. DERIVED-TYPE DECLARATIONS
! 1)   Limiter iso:
!         limiter_iso  
!     Contains the data to define a scalar function of position, so that the
!     iso-contour f=0 of the function corresponds to a geometric limit to the
!     plasma
!        Also contains data to specify the minimum number poloidal points on the
!     s=1 surface to use, and the toroidal planes on which to s=1 surface will
!     be evaluated
!
!    Second attempt:
!     e = SUM_over_i(0,4)_j(0,4) [
!       arz(i,j) (r - rc)^i (z - zc)^j]
!     f = e / |grad(e)|
!
!  Note 
!   1) For now, axisymmetric
!   2) Easy to put in circles, ellipses, and planes
!   3) Ffnction f is approximately distance.
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Declare type limiter_iso
!       arz(i,j) (r - rc)^i (z - zc)^j]
!     f = e / |grad(e)|
!   arz         Real array (0:4,0:4), coefficients for function
!   rc          Real, offset for function
!   zc          Real, offset for function!   zc          Real array (4), coefficients for function
!   numin       Integer, minimum number of poloidal angles for s=1 surface
!   vgrid       Real array (POINTER), values of toroidal angle at which
!                  to compute s=1 surface.
!      ^ USE AS ALLOCATABLE ARRAY
!-------------------------------------------------------------------------------
      TYPE limiter_iso
         REAL(rprec), DIMENSION(0:4,0:4)       :: arz
         REAL(rprec)                           :: rc, zc
         INTEGER                               :: numin
         REAL(rprec), DIMENSION(:), POINTER    :: vgrid => null()
       END TYPE limiter_iso

!*******************************************************************************
! SECTION III. INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Assignment for structures
!-------------------------------------------------------------------------------
      INTERFACE ASSIGNMENT (=)
         MODULE PROCEDURE limiter_iso_assign, limiter_iso_assign_a
      END INTERFACE

!-------------------------------------------------------------------------------
!  End of interface blocks
!-------------------------------------------------------------------------------
      CONTAINS
!*******************************************************************************
! SECTION IV. CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a limiter_iso
!-------------------------------------------------------------------------------
      SUBROUTINE limiter_iso_construct(this,arz,rc,zc,numin,vgrid)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (limiter_iso), INTENT(inout)                     :: this
      REAL(rprec), DIMENSION(0:4,0:4), INTENT(in), OPTIONAL :: arz
      REAL(rprec), INTENT(in), OPTIONAL                     :: rc
      REAL(rprec), INTENT(in), OPTIONAL                     :: zc
      INTEGER, INTENT(in), OPTIONAL                         :: numin
      REAL(rprec), DIMENSION(:), INTENT(in), OPTIONAL       :: vgrid

!  Declare local variables
      INTEGER :: ier1
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'limiter_iso_construct: '

!  Start of executable code
      IF (PRESENT(arz)) THEN
         this % arz = arz
      ELSE
!  Default case is a circle of radius 1.
         this % arz = zero
         this % arz(2,0) = one
         this % arz(0,2) = one
         this % arz(0,0) = - one
      ENDIF

      IF (PRESENT(rc)) THEN
         this % rc = rc
      ELSE
         this % rc = one
      ENDIF

      IF (PRESENT(zc)) THEN
         this % zc = zc
      ELSE
         this % zc = zero
      ENDIF

      IF (PRESENT(numin)) THEN
         this % numin = numin
      ELSE
         this % numin = 1
      ENDIF
      
      IF (PRESENT(vgrid)) THEN
         ALLOCATE(this % vgrid(1:SIZE(vgrid)),STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'alloc')
         this % vgrid = vgrid
      ELSE
         ALLOCATE(this % vgrid(1),STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'alloc')
         this % vgrid(1) = zero
      ENDIF
      
      END SUBROUTINE limiter_iso_construct
      
!*******************************************************************************
! SECTION V. DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Destroy a limiter_iso
!-------------------------------------------------------------------------------
      SUBROUTINE limiter_iso_destroy(this)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (limiter_iso), INTENT(inout) :: this

!  Declare local variables
      INTEGER :: ier1
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'limiter_iso_destroy: '

!  Start of executable code

!  Get rid of all components
      this % arz = zero
      this % rc = zero
      this % zc = zero
      this % numin = 0
      IF (ASSOCIATED(this % vgrid)) THEN
         DEALLOCATE(this % vgrid,STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'deallocate error')
      ENDIF

      END SUBROUTINE limiter_iso_destroy

!*******************************************************************************
! SECTION VI. ASSIGNMENT SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Assignment for limiter_iso
!-------------------------------------------------------------------------------
      SUBROUTINE limiter_iso_assign(left,right)

      IMPLICIT NONE
      
!  Declare Arguments 
      TYPE (limiter_iso), INTENT (inout) :: left
      TYPE (limiter_iso), INTENT (in) :: right
      
!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'limiter_iso_assign: '
      CHARACTER (len=*), PARAMETER :: err_mess1 =                              &
     & 'left-right vgrid are the same?. FIX IT'
      INTEGER :: ier1
         
!  Start of executable code

!  Check to see if the 'use as allocatable array' pointers are pointing
!  to the same location
      CALL assert(.not.ASSOCIATED(left % vgrid,right % vgrid),                 &
     &   sub_name // err_mess1)

!  Destroy left
      CALL limiter_iso_destroy(left)
      
      left % arz = right % arz
      left % rc = right % rc
      left % zc = right % zc
      left % numin = right % numin

!  Allocate space for the 'use as allocatable array' pointer variables
      ALLOCATE(left % vgrid(1:SIZE(right % vgrid)),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'Allocation error')

      left % vgrid = right % vgrid
         
      END SUBROUTINE limiter_iso_assign

!-------------------------------------------------------------------------------
!  Assignment for array of limiter_iso
!-------------------------------------------------------------------------------
      SUBROUTINE limiter_iso_assign_a(left,right)

      IMPLICIT NONE
      
!  Declare Arguments 
      TYPE (limiter_iso), DIMENSION(:), INTENT (inout) :: left
      TYPE (limiter_iso), DIMENSION(:), INTENT (in) :: right
      
!  Declare local variables
      INTEGER :: n_left, n_right, i
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'limiter_iso_assign_a: '
      CHARACTER (len=*), PARAMETER :: err_mess1 =                              &
     & 'left-right array lengths are not the same.'
         
!  Start of executable code

      n_left = SIZE(left)
      n_right = SIZE(right)
      CALL assert_eq(n_left,n_right,sub_name // err_mess1)
      DO i = 1,n_left
         left(i) = right(i)
      END DO
         
      END SUBROUTINE limiter_iso_assign_a

!*******************************************************************************
! SECTION VII.  OUTPUT SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Write out the contents of a limiter_iso
!-------------------------------------------------------------------------------

      SUBROUTINE limiter_iso_write(this,identifier,unit,verbose)
      IMPLICIT NONE

!  Declare Arguments
      TYPE (limiter_iso), INTENT (in) :: this
      CHARACTER (len=*), INTENT(in), OPTIONAL :: identifier
      INTEGER, INTENT(in), OPTIONAL :: unit
      INTEGER, INTENT(in), OPTIONAL :: verbose
!  identifier   character variable, also written out
!  unit         I/O unit number to write to
!  verbose      integer, to specify verbosity level of write

!  Declare local variables and constants
      INTEGER :: iv_default = 1
      INTEGER :: iv
      INTEGER :: iou_default = 6
      INTEGER :: iou
      CHARACTER (len=60) :: id

!  Declare Format array
      CHARACTER(len=*), PARAMETER, DIMENSION(7) :: fmt1 = (/                   &
     & '(" start limiter_iso write, called with id = ",a)         ',           &
     & '(" arz = ",5(3x,es12.5))                                  ',           &
     & '(" rc = ",3x,es12.5)                                      ',           &
     & '(" zc = ",3x,es12.5)                                      ',           &
     & '(" numin = ",i6)                                          ',           &
     & '(" vgrid = ",5(3x,es12.5))                                ',           &
     & '(" end limiter_iso write, called with id = ",a)           '            &
     &  /) 

!  start of executable code
!  Check for arguments present
      IF (PRESENT(identifier)) THEN
         id = identifier
      ELSE
         id = ' '
      END IF

      IF (PRESENT(unit)) THEN
         iou = unit
      ELSE
         iou = iou_default
      END IF

      IF (PRESENT(verbose)) THEN
         iv = verbose
      ELSE
         iv = iv_default
      END IF
      
!  Select Case of Verbosity Level
!  Will need Select Case on s_type
      SELECT CASE(iv)
      CASE( :0)  ! VERY Terse
         WRITE(iou,*) this % arz
         WRITE(iou,*) this % rc
         WRITE(iou,*) this % zc
         WRITE(iou,*) this % numin
         WRITE(iou,*) this % vgrid

      
      CASE(1:)    ! Default, more verbose
         WRITE(iou,fmt1(1)) id
         WRITE(iou,fmt1(2)) this % arz
         WRITE(iou,fmt1(3)) this % rc
         WRITE(iou,fmt1(4)) this % zc
         WRITE(iou,fmt1(5)) this % numin
         WRITE(iou,fmt1(6)) this % vgrid
         WRITE(iou,fmt1(7)) id
      
      END SELECT

      END SUBROUTINE limiter_iso_write

!*******************************************************************************
! SECTION VIII. FUNCTION EVALUATION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Evaluate the limiter function
!-------------------------------------------------------------------------------
      SUBROUTINE limiter_iso_f_cyl(this,rpz_arg,fval)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (limiter_iso), INTENT(in)          :: this
      REAL(rprec), DIMENSION(3), INTENT(in)   :: rpz_arg
      REAL(rprec), INTENT(out)                :: fval

!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'limiter_iso_f_cyl: '
      REAL(rprec)                 :: e, er, ez, grad_e
      INTEGER                     :: i,j

!  Start of executable code

      e = zero
      DO i = 0,4
      DO j = 0,4
         e = e + 
     &      this % arz(i,j) * ((rpz_arg(1) - this % rc) ** i) *                &
     &         ((rpz_arg(3) - this % zc) ** j)
      END DO
      END DO

      er = zero
      DO i = 1,4
      DO j = 0,4
         er = er + 
     &      this % arz(i,j) * i * ((rpz_arg(1) - this % rc) ** (i-1)) *        &
     &         ((rpz_arg(3) - this % zc) ** j)
      END DO
      END DO
      
      ez = zero
      DO i = 0,4
      DO j = 1,4
         ez = ez + 
     &      this % arz(i,j) * j * ((rpz_arg(1) - this % rc) ** i) *            &
     &         ((rpz_arg(3) - this % zc) ** (j - 1))
      END DO
      END DO
      
      grad_e = MAX(1.E-12,SQRT(er * er + ez * ez))
      
      fval = e / grad_e
      
      END SUBROUTINE limiter_iso_f_cyl
      
!*******************************************************************************
! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!  
!  JDH 2009-01-14
!     Modifying signal_T.f to get limiter_function_T
!
!  JDH 2009-01-31
!     Changed name to limiter_iso_T. Added numin and vgrid components.
!
!  JDH 2009-04-10
!     Changed details of function parameterization. (second attempt)
!
      END MODULE limiter_iso_T
