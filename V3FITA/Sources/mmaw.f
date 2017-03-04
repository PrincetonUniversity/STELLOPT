!     SPH: INTEGER(iprec) -> INTEGER
!*******************************************************************************
! SECTION I. MODULE mmaw
!*******************************************************************************
      MODULE mmaw
!-------------------------------------------------------------------------------
!  PURPOSE
!-------------------------------------------------------------------------------
!
!    This module is to help write output files from Fortran that can EASILY
!      be read in to Mathematica, with a GET["filename"] or a copy and paste.
!
!  The big issue is that Mathematica uses *^ to indicate an exponent, and a
!  simple Fortran write uses E to indicate an exponent.
!  This also takes care of the {} necessary to get arrays read in correctly.
!
!-------------------------------------------------------------------------------
!   CHANGE HISTORY
!-------------------------------------------------------------------------------
!
!  See Section V, at end of file
!
!-------------------------------------------------------------------------------
!   USAGE
!-------------------------------------------------------------------------------
!
!  1) Pick the iou number (Default is 6)
!     mmaw_iou = 57
!
!  2) Write an integer (works for scalar, or one-dimensional array)
!     CALL mmaw_wint('myint',myint(1:13))
!
!  3) Write a real (works for scalar, one, two, or three-dimensional arrays)
!     CALL mmaw_wreal('my_x',my_x(1:2,1:12))
!
!  4) Write a logical (scalar only)
!     CALL mmaw_wlog('my_log',my_log)
!
!  5) Change the number of digits in the floating point output (Default is 12)
!     mmaw_rdg = 7
!
!-------------------------------------------------------------------------------


      IMPLICIT NONE
!*******************************************************************************
!    SubSection A     VARIABLE DECLARATIONS
!*******************************************************************************

!-------------------------------------------------------------------------------
!  Type declarations - lengths of reals, integers, and complexes.
!-------------------------------------------------------------------------------
      INTEGER, PARAMETER :: rprec = SELECTED_REAL_KIND(12,100)
      INTEGER, PARAMETER :: iprec = SELECTED_INT_KIND(8)
      
!-------------------------------------------------------------------------------
!  Make type declarations and constants Private, so there are no conflicts.
!-------------------------------------------------------------------------------
      PRIVATE rprec, iprec
!-------------------------------------------------------------------------------
!  I/O Unit Number
!  Number of digits for reals
!-------------------------------------------------------------------------------
      INTEGER :: mmaw_iou=6
      INTEGER :: mmaw_rdg=12

!*******************************************************************************
!    SubSection B     INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Generic Write
!-------------------------------------------------------------------------------

      INTERFACE mmaw_wreal
         MODULE PROCEDURE mmaw_wreal0i, mmaw_wreal1i, mmaw_wreal1in,           &
     &      mmaw_wreal2i, mmaw_wreal2in, mmaw_wreal3i, mmaw_wreal3in
      END INTERFACE

      INTERFACE mmaw_wint
         MODULE PROCEDURE mmaw_wint0i, mmaw_wint1i, mmaw_wint1in
      END INTERFACE

!      END INTERFACE
!      INTERFACE mmaw_write
!         MODULE PROCEDURE mmaw_write_int, mmaw_write_real,                     &
!     &                    mmaw_write_real_a1, mmaw_write_real_a2,              &
!     &                    mmaw_write_real_a3
!      END INTERFACE

      CONTAINS
!*******************************************************************************
!    SubSection C     WRITING SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Write a comment
!-------------------------------------------------------------------------------
      SUBROUTINE mmaw_write_comment(cstring)

      IMPLICIT NONE

!  Declare Arguments 
      CHARACTER(len=*), INTENT(in)   :: cstring

!  Declare local variables

!  Start of executable code
      WRITE(mmaw_iou,'("(* ",a," *)")') cstring
      
      END SUBROUTINE mmaw_write_comment

!-------------------------------------------------------------------------------
!  Write an integer
!-------------------------------------------------------------------------------
      SUBROUTINE mmaw_wint0i(cname, iarg)

      IMPLICIT NONE

!  Declare Arguments 
      CHARACTER(len=*), INTENT(in)   :: cname
      INTEGER, INTENT(in)     :: iarg

!  Declare local variables
      CHARACTER(len=24)   :: cvalue

!  Start of executable code
      WRITE(mmaw_iou,'(a," = ",i0)') cname, iarg
      
      END SUBROUTINE mmaw_wint0i
!-------------------------------------------------------------------------------
!  Write an integer array, one index, no name
!-------------------------------------------------------------------------------
      SUBROUTINE mmaw_wint1i(iarg)

      IMPLICIT NONE

!  Declare Arguments 
      INTEGER, INTENT(in), DIMENSION(:)     :: iarg

!  Declare local variables
      CHARACTER(len=26)   :: cvalue
      INTEGER :: i, imin, imax

!  Start of executable code
      imin = lbound(iarg,1)
      imax = ubound(iarg,1)
      DO i = imin,imax
         WRITE(cvalue,'(i0)') iarg(i)
         IF (i .eq. imin) THEN
            cvalue = '{' // cvalue
         ENDIF
         IF (i .eq. imax) THEN
            cvalue = TRIM(cvalue) // '}'
         ELSE
            cvalue = TRIM(cvalue) // ','
         ENDIF
         WRITE(mmaw_iou,'(a)') cvalue
      END DO
      
      END SUBROUTINE mmaw_wint1i
!-------------------------------------------------------------------------------
!  Write an integer array, one index, with name
!-------------------------------------------------------------------------------
      SUBROUTINE mmaw_wint1in(cname,iarg)

      IMPLICIT NONE

!  Declare Arguments 
      CHARACTER(len=*), INTENT(in)                :: cname
      INTEGER, INTENT(in), DIMENSION(:)    :: iarg

!  Declare local variables

!  Start of executable code
      WRITE(mmaw_iou,'(a,a)') cname, '='
      CALL mmaw_wint1i(iarg)

      END SUBROUTINE mmaw_wint1in
!-------------------------------------------------------------------------------
!  Write a logical
!-------------------------------------------------------------------------------
      SUBROUTINE mmaw_wlog(cname, larg)

      IMPLICIT NONE

!  Declare Arguments 
      CHARACTER(len=*), INTENT(in)   :: cname
      LOGICAL, INTENT(in)            :: larg

!  Declare local variables

!  Start of executable code
      IF (larg) THEN
         WRITE(mmaw_iou,'(a," = ",a)') cname, 'True'
      ELSE
         WRITE(mmaw_iou,'(a," = ",a)') cname, 'False'
      ENDIF
      
      END SUBROUTINE mmaw_wlog
!-------------------------------------------------------------------------------
!  Write a real, no index
!-------------------------------------------------------------------------------
      SUBROUTINE mmaw_wreal0i(cname, xarg)

      IMPLICIT NONE

!  Declare Arguments 
      CHARACTER(len=*), INTENT(in)   :: cname
      REAL(rprec), INTENT(in)        :: xarg

!  Declare local variables
      CHARACTER(len=24)   :: cvalue

!  Start of executable code
      CALL mmaw_r2s(xarg,mmaw_rdg,cvalue)
      WRITE(mmaw_iou,'(a," = ",a)') cname, cvalue
      
      END SUBROUTINE mmaw_wreal0i
!-------------------------------------------------------------------------------
!  Write a real, one index, no name
!-------------------------------------------------------------------------------
      SUBROUTINE mmaw_wreal1i(xarg)

      IMPLICIT NONE

!  Declare Arguments 
      REAL(rprec), INTENT(in), DIMENSION(:)  :: xarg

!  Declare local variables
      CHARACTER(len=26)   :: cvalue
      INTEGER :: i, imin, imax

!  Start of executable code
      imin = lbound(xarg,1)
      imax = ubound(xarg,1)
      DO i = imin,imax
         CALL mmaw_r2s(xarg(i),mmaw_rdg,cvalue)
         IF (i .eq. imin) THEN
            cvalue = '{' // cvalue
         ENDIF
         IF (i .eq. imax) THEN
            cvalue = TRIM(cvalue) // '}'
         ELSE
            cvalue = TRIM(cvalue) // ','
         ENDIF
         WRITE(mmaw_iou,'(a)') cvalue
      END DO

      END SUBROUTINE mmaw_wreal1i
!-------------------------------------------------------------------------------
!  Write a real, one index, with name
!-------------------------------------------------------------------------------
      SUBROUTINE mmaw_wreal1in(cname,xarg)

      IMPLICIT NONE

!  Declare Arguments 
      CHARACTER(len=*), INTENT(in)             :: cname
      REAL(rprec), INTENT(in), DIMENSION(:)  :: xarg

!  Declare local variables

!  Start of executable code
      WRITE(mmaw_iou,'(a,a)') cname, '='
      CALL mmaw_wreal1i(xarg)

      END SUBROUTINE mmaw_wreal1in
!-------------------------------------------------------------------------------
!  Write a real, two indices, no name
!-------------------------------------------------------------------------------
      SUBROUTINE mmaw_wreal2i(xarg)

      IMPLICIT NONE

!  Declare Arguments 
      REAL(rprec), INTENT(in), DIMENSION(:,:)  :: xarg

!  Declare local variables
      INTEGER :: i, imin, imax

!  Start of executable code
      imin = lbound(xarg,1)
      imax = ubound(xarg,1)
      DO i = imin,imax
         IF (i .eq. imin) THEN
            WRITE(mmaw_iou,'(a)') '{'
         ELSE
            WRITE(mmaw_iou,'(a)') ','
         ENDIF
         CALL mmaw_wreal1i(xarg(i,:))
      END DO
      WRITE(mmaw_iou,'(a)') '}'
      
      END SUBROUTINE mmaw_wreal2i
!-------------------------------------------------------------------------------
!  Write a real, two indices, with name
!-------------------------------------------------------------------------------
      SUBROUTINE mmaw_wreal2in(cname,xarg)

      IMPLICIT NONE

!  Declare Arguments 
      CHARACTER(len=*), INTENT(in)             :: cname
      REAL(rprec), INTENT(in), DIMENSION(:,:)  :: xarg

!  Declare local variables

!  Start of executable code
      WRITE(mmaw_iou,'(a,a)') cname, '='
      CALL mmaw_wreal2i(xarg)

      END SUBROUTINE mmaw_wreal2in
!-------------------------------------------------------------------------------
!  Write a real, three indices, no name
!-------------------------------------------------------------------------------
      SUBROUTINE mmaw_wreal3i(xarg)

      IMPLICIT NONE

!  Declare Arguments 
      REAL(rprec), INTENT(in), DIMENSION(:,:,:)  :: xarg

!  Declare local variables
      INTEGER :: i, imin, imax

!  Start of executable code
      imin = lbound(xarg,1)
      imax = ubound(xarg,1)
      DO i = imin,imax
         IF (i .eq. imin) THEN
            WRITE(mmaw_iou,'(a)') '{'
         ELSE
            WRITE(mmaw_iou,'(a)') ','
         ENDIF
         CALL mmaw_wreal2i(xarg(i,:,:))
      END DO
      WRITE(mmaw_iou,'(a)') '}'
      
      END SUBROUTINE mmaw_wreal3i
!-------------------------------------------------------------------------------
!  Write a real, three indices, with name
!-------------------------------------------------------------------------------
      SUBROUTINE mmaw_wreal3in(cname,xarg)

      IMPLICIT NONE

!  Declare Arguments 
      CHARACTER(len=*), INTENT(in)               :: cname
      REAL(rprec), INTENT(in), DIMENSION(:,:,:)  :: xarg

!  Declare local variables

!  Start of executable code
      WRITE(mmaw_iou,'(a,a)') cname, '='
      CALL mmaw_wreal3i(xarg)

      END SUBROUTINE mmaw_wreal3in
      
!*******************************************************************************
!    SubSection D     AUXILLIARY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Convert a real to a string
!     This replaces the D+nn or E+nn FORTRAN format with Mathematica's
!     "`*^" floating exponent indicator
!-------------------------------------------------------------------------------
      SUBROUTINE mmaw_r2s(xarg,idigits,cstring)

      IMPLICIT NONE

!  Declare Arguments 
      REAL(rprec), INTENT(in)      :: xarg
      INTEGER, INTENT(in)   :: idigits   ! number of digits past decimal.
      CHARACTER(len=*), INTENT(out):: cstring

!  Declare local variables
      INTEGER   :: idigits_use     ! Digits after decimal to print out
      INTEGER   :: iexp            ! exponent
      INTEGER   :: iint            ! length of integer format
      REAL(rprec)      :: xmant           ! mantissa
      REAL(rprec)      :: my_tiny=1.d-300 ! Tiny number, to check for zero
      CHARACTER(len=100) :: temp_string	  ! Characters of floating point
      INTEGER   :: id, ilen

!  Start of executable code
      idigits_use = min(14,max(0,idigits))      ! Make sure between 0 and 14
      IF (ABS(xarg) .lt. my_tiny) THEN
!       Small enough that consider zero
         temp_string = ' 0.00000000000000`*^+0'
      ELSE
!      Find exponent and mantissa
         iexp = FLOOR(LOG10(ABS(xarg)))
         xmant = xarg * 10.D0 ** (-iexp)
         WRITE(temp_string,'(f17.14,"`*^",sp,i0)') xmant, iexp
      ENDIF
!  Write to temporary string. 
!  Note format uses 'sp', - print optional '+'
!  Note format uses i0 - prints integer to shortest allowable space
         

!  Trim down the printed mantissa
!  Note that the value is NOT rounded.
      id = idigits_use + 3
      ilen = len_trim(temp_string)
      temp_string = temp_string(1:id) //  temp_string(18:ilen) 
      
      cstring = TRIM(temp_string)
      
      END SUBROUTINE mmaw_r2s

!*******************************************************************************
!    SubSection E     COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!  
!  JDH 05-19-2006
!     Start on module
!
!  JDH 05-22-2006
!     Quick and Easy - one real per line.
!
!  JDH 05-23-2006
!     Add logical
!
!  JDH 10-06-2006
!     Pulled out of the wout2m code, for stand-alone use. Some editing.
           
      END MODULE mmaw

!*******************************************************************************
