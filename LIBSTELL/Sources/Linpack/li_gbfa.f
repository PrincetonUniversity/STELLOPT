      SUBROUTINE sgbfa1 (ABD, LDA, N, ML, MU, IPVT, INFO)
      USE LIPREC, ONLY: WP => SP
      IMPLICIT NONE
!-----------------------------------------------
!   D U M M Y   A R G U M E N T S
!-----------------------------------------------
      INTEGER LDA, N, ML, MU, INFO
      INTEGER, DIMENSION(N) :: IPVT
      REAL(WP), DIMENSION(*) :: ABD

      CALL sgbfa (ABD, LDA, N, ML, MU, IPVT, INFO)

      END SUBROUTINE sgbfa1 


      SUBROUTINE dgbfa1 (ABD, LDA, N, ML, MU, IPVT, INFO)
      USE LIPREC, ONLY: WP => DP
      IMPLICIT NONE
!-----------------------------------------------
!   D U M M Y   A R G U M E N T S
!-----------------------------------------------
      INTEGER LDA, N, ML, MU, INFO
      INTEGER, DIMENSION(N) :: IPVT
      REAL(WP), DIMENSION(*) :: ABD

      CALL dgbfa (ABD, LDA, N, ML, MU, IPVT, INFO)

      END SUBROUTINE dgbfa1 

#ifndef CRAY
      SUBROUTINE sgbfa (ABD, LDA, N, ML, MU, IPVT, INFO)
      USE LIPREC, ONLY: WP => SP
      IMPLICIT NONE
!-----------------------------------------------
!   D U M M Y   A R G U M E N T S
!-----------------------------------------------
      INTEGER LDA, N, ML, MU, INFO
      INTEGER, DIMENSION(N) :: IPVT
      REAL(WP), DIMENSION(LDA,N) :: ABD
!-----------------------------------------------
!   L O C A L   V A R I A B L E S
!-----------------------------------------------
      REAL(WP), PARAMETER :: ZERO = 0, ONE = 1
      INTEGER :: I0, J, JU, JZ, J0, J1, K, KP1, L, LM, M, MM, NM1
      REAL(WP) :: T
      INTEGER :: ISAMAX1(1)
!-----------------------------------------------
      EXTERNAL SSCAL, SAXPY
!
!     sgbfa factors a REAL band matrix by elimination.
!
!     sgbfa is usually called by sgbco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!
!     on entry
!
!        abd     real(lda, n)
!                contains the matrix in band storage.  the columns
!                of the matrix are stored in the columns of  abd  and
!                the diagonals of the matrix are stored in rows
!                ml+1 through 2*ml+mu+1 of  abd .
!                see the comments below for details.
!
!        lda     integer
!                the leading dimension of the array  abd .
!                lda must be .ge. 2*ml + mu + 1 .
!
!        n       integer
!                the order of the original matrix.
!
!        ml      integer
!                number of diagonals below the main diagonal.
!                0 .le. ml .lt. n .
!
!        mu      integer
!                number of diagonals above the main diagonal.
!                0 .le. mu .lt. n .
!                more efficient if  ml .le. mu .
!     on return
!
!        abd     an upper triangular matrix in band storage and
!                the multipliers which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        info    integer
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that sgbsl will divide by zero if
!                     called.  use  rcond  in sgbco for a reliable
!                     indication of singularity.
!
!     band storage
!
!           if  a  is a band matrix, the following program segment
!           will set up the input.
!
!                   ml = (band width below the diagonal)
!                   mu = (band width above the diagonal)
!                   m = ml + mu + 1
!                   do 20 j = 1, n
!                      i1 = max(1, j-mu)
!                      i2 = min(n, j+ml)
!                      do 10 i = i1, i2
!                         k = i - j + m
!                         abd(k,j) = a(i,j)
!                10    continue
!                20 continue
!
!           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
!           in addition, the first  ml  rows in  abd  are used for
!           elements generated during the triangularization.
!           the total number of rows needed in  abd  is  2*ml+mu+1 .
!           the  ml+mu by ml+mu  upper left triangle and the
!           ml by ml  lower right triangle are not referenced.
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas saxpy ,sscal
!     fortran max,min
!
!     internal variables
!
!
      m = ml + mu + 1
      info = 0
!
!     zero initial fill-in columns
!
      j0 = mu + 2
      j1 = MIN(n,m) - 1
      DO jz = j0, j1
         i0 = m + 1 - jz
         abd(i0:ml,jz) = zero
      END DO
      jz = j1
      ju = 0
!
!     gaussian elimination with partial pivoting
!
      nm1 = n - 1
      IF (nm1 .ge. 1) THEN
         DO k = 1, nm1
            kp1 = k + 1
!
!        zero next fill-in column
!
            jz = jz + 1
            IF (jz .le. n) THEN
               abd(:ml,jz) = zero
            END IF
!
!        find l = pivot INDEX
!
            lm = MIN(ml,n - k)
            isamax1 = MAXLOC(ABS(abd(m:m+lm,k)))
            l = isamax1(1) + m - 1
!           l = isamax(lm + 1,abd(m,k),1) + m - 1
            ipvt(k) = l + k - m
!
!        zero pivot implies this column already triangularized
!
            IF (abd(l,k) .ne. zero) THEN
!
!           interchange if necessary
!
               IF (l .ne. m) THEN
                  t = abd(l,k)
                  abd(l,k) = abd(m,k)
                  abd(m,k) = t
               END IF
!
!           compute multipliers
!
               t = -one/abd(m,k)
               CALL sscal (lm, t, abd(m+1,k), 1)
!
!           row elimination with column indexing
!
               ju = MIN(MAX(ju,mu + ipvt(k)),n)
               mm = m
               DO j = kp1, ju
                  l = l - 1
                  mm = mm - 1
                  t = abd(l,j)
                  IF (l .ne. mm) THEN
                     abd(l,j) = abd(mm,j)
                     abd(mm,j) = t
                  END IF
                  CALL saxpy (lm, t, abd(m+1,k), 1, abd(mm+1,j), 1)
               END DO
            ELSE
               info = k
            END IF
         END DO
      END IF
      ipvt(n) = n
      IF (abd(m,n) .eq. zero) info = n

      END SUBROUTINE sgbfa

      SUBROUTINE dgbfa (ABD, LDA, N, ML, MU, IPVT, INFO)
      USE LIPREC, ONLY: WP => DP
      IMPLICIT NONE
!-----------------------------------------------
!   D U M M Y   A R G U M E N T S
!-----------------------------------------------
      INTEGER LDA, N, ML, MU, INFO
      INTEGER, DIMENSION(N) :: IPVT
      REAL(WP), DIMENSION(LDA,N) :: ABD
!-----------------------------------------------
!   L O C A L   V A R I A B L E S
!-----------------------------------------------
      REAL(WP), PARAMETER :: ZERO = 0, ONE = 1
      INTEGER :: I0, J, JU, JZ, J0, J1, K, KP1, L, LM, M, MM, NM1
      REAL(WP) :: T
      INTEGER :: ISAMAX1(1)
!-----------------------------------------------
      EXTERNAL DSCAL, DAXPY
!
!     sgbfa factors a REAL band matrix by elimination.
!
!     sgbfa is usually called by sgbco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!
!     on entry
!
!        abd     real(lda, n)
!                contains the matrix in band storage.  the columns
!                of the matrix are stored in the columns of  abd  and
!                the diagonals of the matrix are stored in rows
!                ml+1 through 2*ml+mu+1 of  abd .
!                see the comments below for details.
!
!        lda     integer
!                the leading dimension of the array  abd .
!                lda must be .ge. 2*ml + mu + 1 .
!
!        n       integer
!                the order of the original matrix.
!
!        ml      integer
!                number of diagonals below the main diagonal.
!                0 .le. ml .lt. n .
!
!        mu      integer
!                number of diagonals above the main diagonal.
!                0 .le. mu .lt. n .
!                more efficient if  ml .le. mu .
!     on return
!
!        abd     an upper triangular matrix in band storage and
!                the multipliers which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        info    integer
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that sgbsl will divide by zero if
!                     called.  use  rcond  in sgbco for a reliable
!                     indication of singularity.
!
!     band storage
!
!           if  a  is a band matrix, the following program segment
!           will set up the input.
!
!                   ml = (band width below the diagonal)
!                   mu = (band width above the diagonal)
!                   m = ml + mu + 1
!                   do 20 j = 1, n
!                      i1 = max(1, j-mu)
!                      i2 = min(n, j+ml)
!                      do 10 i = i1, i2
!                         k = i - j + m
!                         abd(k,j) = a(i,j)
!                10    continue
!                20 continue
!
!           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
!           in addition, the first  ml  rows in  abd  are used for
!           elements generated during the triangularization.
!           the total number of rows needed in  abd  is  2*ml+mu+1 .
!           the  ml+mu by ml+mu  upper left triangle and the
!           ml by ml  lower right triangle are not referenced.
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas saxpy ,sscal
!     fortran max,min
!
!     internal variables
!
      m = ml + mu + 1
      info = 0
!
!     zero initial fill-in columns
!
      j0 = mu + 2
      j1 = MIN(n,m) - 1
      DO jz = j0, j1
         i0 = m + 1 - jz
         abd(i0:ml,jz) = zero
      END DO
      jz = j1
      ju = 0
!
!     gaussian elimination with partial pivoting
!
      nm1 = n - 1
      IF (nm1 .ge. 1) THEN
         DO k = 1, nm1
            kp1 = k + 1
!
!        zero next fill-in column
!
            jz = jz + 1
            IF (jz .le. n) THEN
               abd(:ml,jz) = zero
            END IF
!
!        find l = pivot INDEX
!
            lm = MIN(ml,n - k)
            isamax1 = MAXLOC(ABS(abd(m:m+lm,k)))
            l = isamax1(1) + m - 1
!           l = isamax(lm + 1,abd(m,k),1) + m - 1
            ipvt(k) = l + k - m
!
!        zero pivot implies this column already triangularized
!
            IF (abd(l,k) .ne. zero) THEN
!
!           interchange if necessary
!
               IF (l .ne. m) THEN
                  t = abd(l,k)
                  abd(l,k) = abd(m,k)
                  abd(m,k) = t
               END IF
!
!           compute multipliers
!
               t = -one/abd(m,k)
               CALL dscal (lm, t, abd(m+1,k), 1)
!
!           row elimination with column indexing
!
               ju = MIN(MAX(ju,mu + ipvt(k)),n)
               mm = m
               DO j = kp1, ju
                  l = l - 1
                  mm = mm - 1
                  t = abd(l,j)
                  IF (l .ne. mm) THEN
                     abd(l,j) = abd(mm,j)
                     abd(mm,j) = t
                  END IF
                  CALL daxpy (lm, t, abd(m+1,k), 1, abd(mm+1,j), 1)
               END DO
            ELSE
               info = k
            END IF
         END DO
      END IF
      ipvt(n) = n
      IF (abd(m,n) .eq. zero) info = n

      END SUBROUTINE dgbfa
#endif
