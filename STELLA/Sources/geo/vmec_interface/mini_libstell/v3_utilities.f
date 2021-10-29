!  JDH 08-12-2004. First version. Some subroutines modeled after Numerical Recipes
!  nrutil subroutines

!*******************************************************************************
!  File v3_utilitlies.f
!  Contains module v3_utilitlies
!  Utility Functions for V3FIT

!*******************************************************************************
!  MODULE v3_utilitlies
!    (Utilities, for the V3FIT code)
! SECTION I.    Variable Declarations
! SECTION II.   Interface Blocks
! SECTION III.  Error Trapping
! SECTION IV.   Input-Output Utilities
!*******************************************************************************

      MODULE v3_utilities

      IMPLICIT NONE
!*******************************************************************************
! SECTION I. Variable Declarations
!*******************************************************************************

!-------------------------------------------------------------------------------
!  Type declarations - lengths of reals, integers, and complexes.
!  Frequently used mathematical constants, lots of extra precision.
!  Make type declarations and constants Private, so there are no conflicts.
!-------------------------------------------------------------------------------
      INTEGER, PARAMETER :: rprec = SELECTED_REAL_KIND(12,100)
      INTEGER, PARAMETER :: iprec = SELECTED_INT_KIND(8)
      INTEGER, PARAMETER :: cprec = KIND((1.0_rprec,1.0_rprec))
      REAL(rprec), PARAMETER :: pi=3.14159265358979323846264338328_rprec
      REAL(rprec), PARAMETER :: twopi=6.28318530717958647692528677_rprec
      REAL(rprec), PARAMETER :: one = 1.0_rprec
      REAL(rprec), PARAMETER :: zero = 0.0_rprec

!      USE stell_kinds, only : rprec, iprec, cprec
!      USE stell_constants, only: pi, twopi, one, zero
      PRIVATE rprec, iprec, cprec, pi, twopi, one, zero

!-------------------------------------------------------------------------------
!  JDH 08-13-04. Perhaps add variable for error IO unit numnber, and print there also
!-------------------------------------------------------------------------------

!*******************************************************************************
! SECTION II. INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!  assert, with varying numbers of arguments
!-------------------------------------------------------------------------------
      INTERFACE assert
         MODULE PROCEDURE assert1,assert2,assert3,assert4,assert_v
      END INTERFACE

!-------------------------------------------------------------------------------
!  assert_eq, with varying numbers of arguments
!-------------------------------------------------------------------------------
      INTERFACE assert_eq
         MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
      END INTERFACE

      CONTAINS

!*******************************************************************************
! SECTION III. Error Trapping
!*******************************************************************************

!-------------------------------------------------------------------------------
!  Assert
!    This subroutine is modeled after 'assert' in the Numerical Recipes in F90.
!    The last argument is a string, the first argument(s) are logicals
!    If any of the logicals are false, an error message is printed, and the
!    subroutine either stops or continues.
!    The optional argument err_class determines if the error is fatal, or just
!    a warning. The default action is fatal (execution stops). If the argument
!    err_class is present, and its first character is 'w' or 'W', then execution
!    continues.
!    Different versions have 1, 2, 3, 4, or a vector of logicals.
!  Assert, 1 logical
!-------------------------------------------------------------------------------
      SUBROUTINE assert1(n1,string,err_class)
!  Argument Declaration
      CHARACTER(LEN=*), INTENT(IN) :: string
      LOGICAL, INTENT(IN) :: n1
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: err_class

!  Local variable declarations
      CHARACTER(LEN=*), PARAMETER :: fmt='(1x,a,/,1x,a)'
      CHARACTER(LEN=1) :: first_char
      LOGICAL :: lfatal

!  Start of executable code
      IF (.not. n1) THEN
         WRITE (*,fmt) 'error: an assertion failed with this tag:',            &
     &      string
         WRITE(*,*) ' n1 = ',n1
!  Is error Fatal or Warning?
         lfatal = .TRUE.
         IF (PRESENT(err_class)) THEN
            first_char = err_class(1:1)
            IF ((first_char .eq. 'w') .or. (first_char .eq. 'W'))
     &         lfatal = .FALSE.
         ENDIF
         IF (lfatal) THEN
            STOP 'program terminated by assert1'
         ELSE
            WRITE(*,*) ' end of warning error message from assert1'
         END IF
      END IF
      END SUBROUTINE assert1

!-------------------------------------------------------------------------------
!  Assert, 2 logicals
!-------------------------------------------------------------------------------
      SUBROUTINE assert2(n1,n2,string,err_class)
!  Argument Declaration
      CHARACTER(LEN=*), INTENT(IN) :: string
      LOGICAL, INTENT(IN) :: n1, n2
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: err_class

!  Local variable declarations
      CHARACTER(LEN=*), PARAMETER :: fmt='(1x,a,/,1x,a)'
      CHARACTER(LEN=1) :: first_char
      LOGICAL :: lfatal

!  Start of executable code
      IF (.not. (n1 .and. n2)) THEN
         WRITE (*,fmt) 'error: an assertion failed with this tag:',            &
     &      string
         WRITE(*,*) ' n1, n2 = ',n1, n2
!  Is error Fatal or Warning?
         lfatal = .TRUE.
         IF (PRESENT(err_class)) THEN
            first_char = err_class(1:1)
            IF ((first_char .eq. 'w') .or. (first_char .eq. 'W'))
     &         lfatal = .FALSE.
         ENDIF
         IF (lfatal) THEN
            STOP 'program terminated by assert2'
         ELSE
            WRITE(*,*) ' end of warning error message from assert2'
         END IF
      END IF
      END SUBROUTINE assert2

!-------------------------------------------------------------------------------
!  Assert, 3 logicals
!-------------------------------------------------------------------------------
      SUBROUTINE assert3(n1,n2,n3,string,err_class)
!  Argument Declaration
      CHARACTER(LEN=*), INTENT(IN) :: string
      LOGICAL, INTENT(IN) :: n1,n2,n3
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: err_class

!  Local variable declarations
      CHARACTER(LEN=*), PARAMETER :: fmt='(1x,a,/,1x,a)'
      CHARACTER(LEN=1) :: first_char
      LOGICAL :: lfatal

!  Start of executable code
      IF (.not. (n1 .and. n2 .and. n3)) THEN
         WRITE (*,fmt) 'error: an assertion failed with this tag:',              &
     &      string
         WRITE(*,*) ' n1, n2, n3 = ',n1, n2, n3
!  Is error Fatal or Warning?
         lfatal = .TRUE.
         IF (PRESENT(err_class)) THEN
            first_char = err_class(1:1)
            IF ((first_char .eq. 'w') .or. (first_char .eq. 'W'))
     &         lfatal = .FALSE.
         ENDIF
         IF (lfatal) THEN
            STOP 'program terminated by assert3'
         ELSE
            WRITE(*,*) ' end of warning error message from assert3'
         END IF
      END IF
      END SUBROUTINE assert3

!-------------------------------------------------------------------------------
!  Assert, 4 logicals
!-------------------------------------------------------------------------------
      SUBROUTINE assert4(n1,n2,n3,n4,string,err_class)
!  Argument Declaration
      CHARACTER(LEN=*), INTENT(IN) :: string
      LOGICAL, INTENT(IN) :: n1,n2,n3,n4
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: err_class

!  Local variable declarations
      CHARACTER(LEN=*), PARAMETER :: fmt='(1x,a,/,1x,a)'
      CHARACTER(LEN=1) :: first_char
      LOGICAL :: lfatal

!  Start of executable code
      IF (.not. (n1 .and. n2 .and. n3 .and. n4)) THEN
         WRITE (*,fmt) 'error: an assertion failed with this tag:',              &
     &      string
         WRITE(*,*) ' n1, n2, n3, n4 = ',n1, n2, n3, n4
!  Is error Fatal or Warning?
         lfatal = .TRUE.
         IF (PRESENT(err_class)) THEN
            first_char = err_class(1:1)
            IF ((first_char .eq. 'w') .or. (first_char .eq. 'W'))
     &         lfatal = .FALSE.
         ENDIF
         IF (lfatal) THEN
            STOP 'program terminated by assert4'
         ELSE
            WRITE(*,*) ' end of warning error message from assert4'
         END IF
      END IF
      END SUBROUTINE assert4

!-------------------------------------------------------------------------------
!  Assert, vector of logicals
!-------------------------------------------------------------------------------
      SUBROUTINE assert_v(n,string,err_class)
!  Argument Declaration
      CHARACTER(LEN=*), INTENT(IN) :: string
      LOGICAL, DIMENSION(:), INTENT(IN) :: n
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: err_class
      
!  Local variable declarations
      CHARACTER(LEN=*), PARAMETER :: fmt='(1x,a,/,1x,a)'
      INTEGER(iprec) :: ntot, nfalse, ifirst, i
      CHARACTER(LEN=1) :: first_char
      LOGICAL :: lfatal

!  Start of executable code
      IF (.not. all(n)) THEN
         WRITE (*,fmt) 'error: an assertion failed with this tag:',              &
     &      string
         ntot = SIZE(n)
         WRITE (*,*) ntot, ' logicals in array. indices of .F. are:'
         nfalse = 0
         DO i = 1,ntot
            IF(.not. n(i)) THEN
               WRITE(*,*) i
               nfalse = nfalse + 1
            END IF
         END DO
         WRITE (*,*) nfalse, ' logicals are false'
!  Is error Fatal or Warning?
         lfatal = .TRUE.
         IF (PRESENT(err_class)) THEN
            first_char = err_class(1:1)
            IF ((first_char .eq. 'w') .or. (first_char .eq. 'W'))
     &         lfatal = .FALSE.
         ENDIF
         IF (lfatal) THEN
            STOP 'program terminated by assert_v'
         ELSE
            WRITE(*,*) ' end of warning error message from assert_v'
         END IF
      END IF
      END SUBROUTINE assert_v

!-------------------------------------------------------------------------------
!  Assert_eq
!    This subroutine is modeled after 'assert_eq' in the Numerical Recipes in F90.
!    The last argument is a string, the first argument are integers
!    If any of the integers is different from the first integer, then
!    an error message is printed, and the subroutine either STOPs or continues.
!    The optional argument err_class determines if the error is fatal, or just
!    a warning. The default action is fatal (execution stops). If the argument
!    err_class is present, and its first character is 'w' or 'W', then execution
!    continues.
!    Different versions have 2, 3, 4, or a vector of integer arguments.
!    Note: NR had this as a FUNCTION. I have changed it to a subroutine.
!  Assert_eq, 2 integer arguments
!-------------------------------------------------------------------------------
      SUBROUTINE assert_eq2(n1,n2,string,err_class)
!  Argument Declaration
      CHARACTER(LEN=*), INTENT(IN) :: string
      INTEGER, INTENT(IN) :: n1,n2
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: err_class

!  Local variable declarations
      CHARACTER(LEN=*), PARAMETER :: fmt='(1x,a,/,1x,a)'
      CHARACTER(LEN=1) :: first_char
      LOGICAL :: lfatal

!  Start of executable code
      IF (.not.(n1 == n2)) THEN
         WRITE (*,fmt) 'error: an assert_eq failed with this tag:',              &
     &      string
         WRITE (*,*) ' n1, n2 = ',n1, n2
!  Is error Fatal or Warning?
         lfatal = .TRUE.
         IF (PRESENT(err_class)) THEN
            first_char = err_class(1:1)
            IF ((first_char .eq. 'w') .or. (first_char .eq. 'W'))
     &         lfatal = .FALSE.
         ENDIF
         IF (lfatal) THEN
            STOP 'program terminated by assert_eq2'
         ELSE
            WRITE(*,*) ' end of warning error message from assert_eq2'
         END IF
      END IF
      END SUBROUTINE assert_eq2

!-------------------------------------------------------------------------------
!  Assert_eq, 3 integer arguments
!-------------------------------------------------------------------------------
      SUBROUTINE assert_eq3(n1,n2,n3,string,err_class)
!  Argument Declaration
      CHARACTER(LEN=*), INTENT(IN) :: string
      INTEGER, INTENT(IN) :: n1,n2,n3
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: err_class

!  Local variable declarations
      CHARACTER(LEN=*), PARAMETER :: fmt='(1x,a,/,1x,a)'
      CHARACTER(LEN=1) :: first_char
      LOGICAL :: lfatal

!  Start of executable code
      IF (.not.(n1 == n2 .and. n2 == n3)) THEN
         WRITE (*,fmt) 'error: an assert_eq failed with this tag:',              &
     &      string
         WRITE (*,*) ' n1, n2, n3 = ',n1, n2, n3
!  Is error Fatal or Warning?
         lfatal = .TRUE.
         IF (PRESENT(err_class)) THEN
            first_char = err_class(1:1)
            IF ((first_char .eq. 'w') .or. (first_char .eq. 'W'))
     &         lfatal = .FALSE.
         ENDIF
         IF (lfatal) THEN
            STOP 'program terminated by assert_eq3'
         ELSE
            WRITE(*,*) ' end of warning error message from assert_eq3'
         END IF
      END IF
      END SUBROUTINE assert_eq3

!-------------------------------------------------------------------------------
!  Assert_eq, 4 integer arguments
!-------------------------------------------------------------------------------
      SUBROUTINE assert_eq4(n1,n2,n3,n4,string,err_class)
!  Argument Declaration
      CHARACTER(LEN=*), INTENT(IN) :: string
      INTEGER, INTENT(IN) :: n1,n2,n3,n4
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: err_class

!  Local variable declarations
      CHARACTER(LEN=*), PARAMETER :: fmt='(1x,a,/,1x,a)'
      CHARACTER(LEN=1) :: first_char
      LOGICAL :: lfatal

!  Start of executable code
      IF (.not.(n1 == n2 .and. n2 == n3 .and. n3 == n4)) THEN
         WRITE (*,fmt) 'error: an assert_eq failed with this tag:',              &
     &      string
         WRITE (*,*) ' n1, n2, n3, n4 = ',n1, n2, n3, n4
!  Is error Fatal or Warning?
         lfatal = .TRUE.
         IF (PRESENT(err_class)) THEN
            first_char = err_class(1:1)
            IF ((first_char .eq. 'w') .or. (first_char .eq. 'W'))
     &         lfatal = .FALSE.
         ENDIF
         IF (lfatal) THEN
            STOP 'program terminated by assert_eq4'
         ELSE
            WRITE(*,*) ' end of warning error message from assert_eq4'
         END IF
      END IF
      END SUBROUTINE assert_eq4

!-------------------------------------------------------------------------------
!  Assert_eq, a vector of integers
!-------------------------------------------------------------------------------
      SUBROUTINE assert_eqn(nn,string,err_class)
!  Argument Declaration
      CHARACTER(LEN=*), INTENT(IN) :: string
      INTEGER, DIMENSION(:), INTENT(IN) :: nn
      INTEGER :: ntot, nne, i
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: err_class

!  Local variable declarations
      CHARACTER(LEN=*), PARAMETER :: fmt='(1x,a,/,1x,a)'
      CHARACTER(LEN=1) :: first_char
      LOGICAL :: lfatal

!  Start of executable code
      IF (.not.(all(nn(2:) == nn(1)))) THEN
         WRITE (*,fmt) 'error: an assert_eq failed with this tag:',              &
     &      string
         ntot = SIZE(nn)
         WRITE (*,*) ntot, ' integers in the array.'
         WRITE (*,*) nn(1), ' is the first value in the array'
         WRITE(*,*) ' index    value  (of those .ne. to first)'
         nne = 0
         DO i = 2,ntot
            IF (nn(i) .ne. nn(1)) THEN
               WRITE (*,*) i,nn(i)
               nne = nne + 1
            END IF
         END DO
         WRITE (*,*) ' There are ',nne, ' integers .ne. to first'
!  Is error Fatal or Warning?
         lfatal = .TRUE.
         IF (PRESENT(err_class)) THEN
            first_char = err_class(1:1)
            IF ((first_char .eq. 'w') .or. (first_char .eq. 'W'))
     &         lfatal = .FALSE.
         ENDIF
         IF (lfatal) THEN
            STOP 'program terminated by assert_eqn'
         ELSE
            WRITE(*,*) ' end of warning error message from assert_eqn'
         END IF
      END IF
      END SUBROUTINE assert_eqn

!-------------------------------------------------------------------------------
!  SPH: Changed int to INTEGER (from INTEGER(iprec))
!  Subroutine err_fatal
!    This subroutine is for fatal errors.
!    The first argument is a string, and there are optional arguments for
!    character, real, integer and logical variables.
!    (The user should use KEYWORDs for the optional arguments)
!    When called, err_fatal prints the first argument, and all of the
!    optional arguments, and then STOPs the program.
!-------------------------------------------------------------------------------
      SUBROUTINE err_fatal(string, char, real1, int, log)
!  Argument Declaration
      CHARACTER(LEN=*), INTENT(IN)            :: string
      CHARACTER(LEN=*), OPTIONAL, INTENT(IN)  :: char
      REAL(rprec), OPTIONAL, INTENT(IN)       :: real1
      INTEGER, OPTIONAL, INTENT(IN)           :: int
      LOGICAL, OPTIONAL, INTENT(IN)           :: log

!  Local variable declarations
      CHARACTER(LEN=*), PARAMETER :: fmt='(1x,a,/,1x,a)'

!  Start of executable code
      WRITE (*,fmt) 'FATAL ERROR: ',string

      IF(PRESENT(char)) THEN
         WRITE (*,*) ' char argument = ', char
      END IF
      
      IF(PRESENT(real1)) THEN
         WRITE (*,*) ' real argument = ', real1
      END IF
      
      IF(PRESENT(int)) THEN
         WRITE (*,*) ' integer argument = ', int
      END IF
      
      IF(PRESENT(log)) THEN
         WRITE (*,*) ' logical argument = ', log
      END IF

      STOP 'program terminated by err_fatal'
      END SUBROUTINE err_fatal

!-------------------------------------------------------------------------------
! Warning
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Subroutine err_warn
!    This subroutine is for warnings. Execution is not stopped.
!    The first argument is a string, and there are optional arguments for
!    character, real, integer and logical variables.
!    (The user should use KEYWORDs for the optional arguments)
!    When called, err_warn prints the first argument, and all of the
!    optional arguments, and then returns.
!-------------------------------------------------------------------------------
      SUBROUTINE err_warn(string, char, real, int, log)
!  Argument Declaration
      CHARACTER(LEN=*), INTENT(IN)            :: string
      CHARACTER(LEN=*), OPTIONAL, INTENT(IN)  :: char
      REAL(rprec), OPTIONAL, INTENT(IN)       :: real
      INTEGER, OPTIONAL, INTENT(IN)           :: int
      LOGICAL, OPTIONAL, INTENT(IN)           :: log

!  Local variable declarations
      CHARACTER(LEN=*), PARAMETER :: fmt='(1x,a,/,1x,a)'

!  Start of executable code
      WRITE (*,fmt) 'WARNING ERROR: ',string

      IF(PRESENT(char)) THEN
         WRITE (*,*) ' char argument = ', char
      END IF
      
      IF(PRESENT(real)) THEN
         WRITE (*,*) ' real argument = ', real
      END IF
      
      IF(PRESENT(int)) THEN
         WRITE (*,*) ' integer argument = ', int
      END IF
      
      IF(PRESENT(log)) THEN
         WRITE (*,*) ' logical argument = ', log
      END IF

      RETURN
      END SUBROUTINE err_warn


!*******************************************************************************
! SECTION IV.   Input-Output Utilities
!*******************************************************************************
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      SUBROUTINE svdproducts(b,svprod,numzeros)

! This subroutine takes as its first argument an arbitrary real matrix B.
! It's purpose is to find the variable in column-space that is
! "most redundant"

! It returns a vector of singular value products, computed by eliminating
! individual columns of B, and finding the product of the Singular Values (SV) 
! for the resulting column-reduced matrix.

! For example, suppose a matrix B has 13 columns, and svprod(7) is the maximum
! of all the SV products. Thus eliminating column 7 yields the reduced matrix
! with the largest SV product. Variable 7 in the column space is then
! the "most redundant" variable.

!  Note that when B is a normalized Jacobian matrix, a surface of
!  constant-chi-squared in the column-space is an 
!  ellipsoid, and the volume of the ellipsoid is proportional to the 
!  reciprocal of the product of the singular values.

! There is a problem if there are TWO (or more) variables that essentially
! give ZERO singular values. (For example, if B has 13 columns, and 
! columns 7 and 9 are both zero.) Then the SV-product array will all be zero, and
! won't be able to distinguish the "most redundant" variable. 
! A solution to this problem is to keep track of the number of zero SVs.
! These are stored in the numzeros array.
! Then the column number with the MINIMUM number of zero singular values
! will be (one of) the "most-redundant" variables.

!  The length of the SV product vector must be equal to the number of columns of B.
!  The length of the numzeros vector must be equal to the number of columns of B
!  JDH 2008-02-21, 24
!  JDH 2008-05-19 - CHanged for SV ratios to SV product value

      USE stel_kinds
      
      IMPLICIT NONE

!  Declare Arguments
      REAL(rprec), DIMENSION(:,:), INTENT(inout) :: b
      REAL(rprec), DIMENSION(:), INTENT(inout) :: svprod
      INTEGER, DIMENSION(:), INTENT(inout) :: numzeros

!  Declare local variables
      INTEGER, DIMENSION(2) :: dimlens
      INTEGER :: nrowb, ncolb, nrow, ncol, minrowcol, l_work_svd, ier1 
      INTEGER :: jelim, info_svd
      REAL(rprec)                              :: svp
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: c
      REAL(rprec), PARAMETER                   :: tiny = 1.e-14_rprec
      REAL(rprec), PARAMETER                 :: verytiny = 1.e-28_rprec
      REAL(rprec), DIMENSION(:), ALLOCATABLE   :: work_svd, svec
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'svdproducts: '

!  Start of executable code

!  Array sizes
      dimlens = SHAPE(b)
      nrowb = dimlens(1)
      ncolb = dimlens(2)
      CALL assert(nrowb .ge. 2,sub_name // 'nrowb too small')
      CALL assert(ncolb .ge. 2,sub_name // 'ncolb too small')
      CALL assert_eq(ncolb,SIZE(svprod),sub_name // 'svprod wrong size')
      CALL assert_eq(ncolb,SIZE(numzeros),sub_name //                          &
     &   'numzeros wrong size')

!  Allocate space for the number-columns-reduced-by-one arrays
      nrow = nrowb
      ncol = ncolb - 1
      minrowcol = min(nrow,ncol)
      l_work_svd = 5 * MAX(nrow,ncol)
      ALLOCATE(work_svd(l_work_svd), stat = ier1)
      CALL assert_eq(0,ier1,sub_name // 'Allocate work_svd')
      ALLOCATE(c(nrow,ncol), stat = ier1)
      CALL assert_eq(0,ier1,sub_name // 'Allocate c')
      ALLOCATE(svec(minrowcol), stat = ier1)
      CALL assert_eq(0,ier1,sub_name // 'Allocate svec')

!  Loop over the column to be eliminated
      DO jelim = 1,ncolb
         c(1:nrow,1:jelim-1) = b(1:nrow,1:jelim-1)
         c(1:nrow,jelim:ncolb-1) = b(1:nrow,jelim+1:ncolb)

! Next line commented out by MJL 2017-08-15 since not needed for mini_libstell.
!         CALL dgesvd('None','None',nrow,ncol,c,nrow,                           &
!     &      svec,c,nrow,c,ncol,work_svd,l_work_svd,info_svd)
         CALL assert_eq(0,info_svd,sub_name // 'dgesvd problem')
         svp = PRODUCT(svec)
         numzeros(jelim) = COUNT(svec < tiny)
         svprod(jelim) = svp
      END DO
      RETURN
      END SUBROUTINE svdproducts
      
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      SUBROUTINE most_redundant(a,ncol_array,svprod_array,j_col_elim)

!  This subroutine orders the variables in the column-space of a matrix A from
!  "most redundant" to "least redundant". Redundancy is measured by an increase
!  in the product of singular values when the column is eliminated.

!  Arguments (Input)
!    A               Input: An array with nrow rows and ncol columns
!    ncol_array      Integer array (ncol): number of columns of A remaining
!    svprod_array    Real array (ncol): singular value product
!    j_col_elim      Integer array (ncol): index of the column of A 
!                    (in the ORIGINAL ORDERING) that was eliminated most recently
!    
!  The subroutine makes repeated calls to the related subroutine svdproducts.
!    
!  JDH 2008-02-21
!  JDH 2008-05-19 - Changed from SV ratios to SV products.

      USE stel_kinds
      
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Declare Arguments
!-------------------------------------------------------------------------------
      REAL(rprec), DIMENSION(:,:), INTENT(inout) :: a
      REAL(rprec), DIMENSION(:), INTENT(inout) :: svprod_array
      INTEGER, DIMENSION(:), INTENT(inout) :: ncol_array, j_col_elim 

!-------------------------------------------------------------------------------
!  Declare local variables
!-------------------------------------------------------------------------------
      REAL(rprec), PARAMETER :: tiny = 1.e-14_rprec
      REAL(rprec), PARAMETER :: verytiny = 1.e-28_rprec
!   Array dimensions, error status, indices, etc
      INTEGER, DIMENSION(2) :: dimlens
      INTEGER :: nrowa, ncola, nrow, ncol, minrowcol, ier1, ncolb
      INTEGER :: i, ncolelim, minnz, k, kk
!   For calling the SVD subroutine
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: work_svd, svec
      INTEGER :: l_work_svd, info_svd
!   Arrays for the call to svdproducts
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: atemp, b
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: svprod
      INTEGER, DIMENSION(:), ALLOCATABLE :: numzeros
!   Keeping track of the j-k indices, locations
      INTEGER, DIMENSION(1) :: k_minnz_a, k_maxprod_a
      INTEGER, DIMENSION(:), ALLOCATABLE :: j_now
!  
      REAL(rprec) :: maxprod
      
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'most_redundant: '

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------

!  Array sizes
      dimlens = SHAPE(a)
      nrowa = dimlens(1)
      ncola = dimlens(2)
      CALL assert(nrowa .ge. 2,sub_name // 'nrowa too small')
      CALL assert(ncola .ge. 2,sub_name // 'ncola too small')
      CALL assert_eq(ncola,SIZE(ncol_array),sub_name //                        &
     &   'ncol_array wrong size')
      CALL assert_eq(ncola,SIZE(svprod_array),sub_name //                      &
     &   'svprod_array wrong size')
      CALL assert_eq(ncola,SIZE(j_col_elim),sub_name //                        &
     &   'j_col_elim wrong size')

!  Allocations for local arrays
      nrow = nrowa
      ncol = ncola 
      minrowcol = min(nrow,ncol)
      l_work_svd = 5 * MAX(nrow,ncol)
      ALLOCATE(work_svd(l_work_svd), stat = ier1)
      CALL assert_eq(0,ier1,sub_name // 'Allocate work_svd')
      ALLOCATE(atemp(nrow,ncol),b(nrow,ncol), stat = ier1)
      CALL assert_eq(0,ier1,sub_name // 'Allocate atemp, b')
      ALLOCATE(svec(minrowcol), stat = ier1)
      CALL assert_eq(0,ier1,sub_name // 'Allocate svec')
      ALLOCATE(svprod(ncola), stat = ier1)
      CALL assert_eq(0,ier1,sub_name // 'Allocate svprod')
      ALLOCATE(numzeros(ncola), stat = ier1)
      CALL assert_eq(0,ier1,sub_name // 'Allocate numzeros')
      ALLOCATE(j_now(ncola), stat = ier1)
      CALL assert_eq(0,ier1,sub_name // 'Allocate j_now')

!  All Columns - do SVD and find ratio
      b = a
! Next line commented out by MJL 2017-08-15 since not needed for mini_libstell.
!      CALL dgesvd('None','None',nrow,ncol,b,nrow,                              &
!     &   svec,b,nrow,b,ncol,work_svd,l_work_svd,info_svd)
      CALL assert_eq(0,info_svd,sub_name // 'dgesvd problem')
      
      ncol_array(1) = ncola
      svprod_array(1) = PRODUCT(svec)
      j_col_elim(1) = 0

!  Initialize array of j-indices
!  Used to convert from k-index (current column number, B matrix) 
!  to j-index (original column number, A matrix)
      DO i = 1,ncola
         j_now(i) = i
      END DO
      
!  Loop over the number of columns to be eliminated
      atemp = a
      DO ncolelim = 1,ncola-2
         ncolb = ncola + 1 - ncolelim
         b(1:nrow,1:ncolb) = atemp(1:nrow,1:ncolb)
         CALL svdproducts(b(1:nrow,1:ncolb),svprod(1:ncolb),                   &
     &      numzeros(1:ncolb))
         minnz = MINVAL(numzeros(1:ncolb))
         k_minnz_a = MINLOC(numzeros(1:ncolb))
         maxprod = MAXVAL(svprod(1:ncolb))
         k_maxprod_a = MAXLOC(svprod(1:ncolb))
         IF (minnz .gt. 0) THEN
            k = k_minnz_a(1)
         ELSE
            k = k_maxprod_a(1)
         ENDIF
!         WRITE(*,*) 'products', svprod(1:ncolb)
!         WRITE(*,*) 'numzeros', numzeros(1:ncolb)
         ncol_array(ncolelim + 1) = ncola - ncolelim
         svprod_array(ncolelim + 1) = maxprod
         j_col_elim(ncolelim + 1) = j_now(k)
         DO kk = k,ncola - ncolelim + 1
            j_now(kk) = j_now(kk + 1)
            atemp(1:nrow,kk) = atemp(1:nrow,kk + 1)
         END DO
      END DO

!  Define last elements. Make o choice between last two variables.
!  Put the last two variables into j_col_elim and
!  into the (real) svprod_array.
      ncol_array(ncola) = 1
      svprod_array(ncola) = j_now(2)
      j_col_elim(ncola) = j_now(1)
      
      RETURN
      END SUBROUTINE most_redundant


      END MODULE v3_utilities
