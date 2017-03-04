!DEC$ IF .NOT.DEFINED (CRAY)
      SUBROUTINE sgbfa (ABD, LDA, N, ML, MU, IPVT, INFO)
      USE LIPREC, ONLY: WP => SP
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER LDA, N, ML, MU, INFO
      INTEGER, DIMENSION(N) :: IPVT
      REAL(WP), DIMENSION(LDA,N) :: ABD
C-----------------------------------------------
C   L O C A L   V A R I A B L E S
C-----------------------------------------------
      REAL(WP), PARAMETER :: ZERO = 0, ONE = 1
      INTEGER :: I0, J, JU, JZ, J0, J1, K, KP1, L, LM, M, MM, NM1
      REAL(WP) :: T
      INTEGER :: ISAMAX1(1)
C-----------------------------------------------
      EXTERNAL SSCAL, SAXPY
c
c     sgbfa factors a REAL band matrix by elimination.
c
c     sgbfa is usually called by sgbco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c
c     on entry
c
c        abd     real(lda, n)
c                contains the matrix in band storage.  the columns
c                of the matrix are stored in the columns of  abd  and
c                the diagonals of the matrix are stored in rows
c                ml+1 through 2*ml+mu+1 of  abd .
c                see the comments below for details.
c
c        lda     integer
c                the leading dimension of the array  abd .
c                lda must be .ge. 2*ml + mu + 1 .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c                0 .le. ml .lt. n .
c
c        mu      integer
c                number of diagonals above the main diagonal.
c                0 .le. mu .lt. n .
c                more efficient if  ml .le. mu .
c     on return
c
c        abd     an upper triangular matrix in band storage and
c                the multipliers which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that sgbsl will divide by zero if
c                     called.  use  rcond  in sgbco for a reliable
c                     indication of singularity.
c
c     band storage
c
c           if  a  is a band matrix, the following program segment
c           will set up the input.
c
c                   ml = (band width below the diagonal)
c                   mu = (band width above the diagonal)
c                   m = ml + mu + 1
c                   do 20 j = 1, n
c                      i1 = max(1, j-mu)
c                      i2 = min(n, j+ml)
c                      do 10 i = i1, i2
c                         k = i - j + m
c                         abd(k,j) = a(i,j)
c                10    continue
c                20 continue
c
c           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
c           in addition, the first  ml  rows in  abd  are used for
c           elements generated during the triangularization.
c           the total number of rows needed in  abd  is  2*ml+mu+1 .
c           the  ml+mu by ml+mu  upper left triangle and the
c           ml by ml  lower right triangle are not referenced.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas saxpy ,sscal
c     fortran max,min
c
c     internal variables
c
c
c
      m = ml + mu + 1
      info = 0
c
c     zero initial fill-in columns
c
      j0 = mu + 2
      j1 = MIN(n,m) - 1
      DO jz = j0, j1
         i0 = m + 1 - jz
         abd(i0:ml,jz) = zero
      END DO
      jz = j1
      ju = 0
c
c     gaussian elimination with partial pivoting
c
      nm1 = n - 1
      IF (nm1 .ge. 1) THEN
         DO k = 1, nm1
            kp1 = k + 1
c
c        zero next fill-in column
c
            jz = jz + 1
            IF (jz .le. n) THEN
               abd(:ml,jz) = zero
            END IF
c
c        find l = pivot INDEX
c
            lm = MIN(ml,n - k)
            isamax1 = MAXLOC(ABS(abd(m:m+lm,k)))
            l = isamax1(1) + m - 1
!           l = isamax(lm + 1,abd(m,k),1) + m - 1
            ipvt(k) = l + k - m
c
c        zero pivot implies this column already triangularized
c
            IF (abd(l,k) .ne. zero) THEN
c
c           interchange if necessary
c
               IF (l .ne. m) THEN
                  t = abd(l,k)
                  abd(l,k) = abd(m,k)
                  abd(m,k) = t
               END IF
c
c           compute multipliers
c
               t = -one/abd(m,k)
               CALL sscal (lm, t, abd(m+1,k), 1)
c
c           row elimination with column indexing
c
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
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER LDA, N, ML, MU, INFO
      INTEGER, DIMENSION(N) :: IPVT
      REAL(WP), DIMENSION(LDA,N) :: ABD
C-----------------------------------------------
C   L O C A L   V A R I A B L E S
C-----------------------------------------------
      REAL(WP), PARAMETER :: ZERO = 0, ONE = 1
      INTEGER :: I0, J, JU, JZ, J0, J1, K, KP1, L, LM, M, MM, NM1
      REAL(WP) :: T
      INTEGER :: ISAMAX1(1)
C-----------------------------------------------
      EXTERNAL DSCAL, DAXPY
c
c     sgbfa factors a REAL band matrix by elimination.
c
c     sgbfa is usually called by sgbco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c
c     on entry
c
c        abd     real(lda, n)
c                contains the matrix in band storage.  the columns
c                of the matrix are stored in the columns of  abd  and
c                the diagonals of the matrix are stored in rows
c                ml+1 through 2*ml+mu+1 of  abd .
c                see the comments below for details.
c
c        lda     integer
c                the leading dimension of the array  abd .
c                lda must be .ge. 2*ml + mu + 1 .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c                0 .le. ml .lt. n .
c
c        mu      integer
c                number of diagonals above the main diagonal.
c                0 .le. mu .lt. n .
c                more efficient if  ml .le. mu .
c     on return
c
c        abd     an upper triangular matrix in band storage and
c                the multipliers which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that sgbsl will divide by zero if
c                     called.  use  rcond  in sgbco for a reliable
c                     indication of singularity.
c
c     band storage
c
c           if  a  is a band matrix, the following program segment
c           will set up the input.
c
c                   ml = (band width below the diagonal)
c                   mu = (band width above the diagonal)
c                   m = ml + mu + 1
c                   do 20 j = 1, n
c                      i1 = max(1, j-mu)
c                      i2 = min(n, j+ml)
c                      do 10 i = i1, i2
c                         k = i - j + m
c                         abd(k,j) = a(i,j)
c                10    continue
c                20 continue
c
c           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
c           in addition, the first  ml  rows in  abd  are used for
c           elements generated during the triangularization.
c           the total number of rows needed in  abd  is  2*ml+mu+1 .
c           the  ml+mu by ml+mu  upper left triangle and the
c           ml by ml  lower right triangle are not referenced.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas saxpy ,sscal
c     fortran max,min
c
c     internal variables
c
c
c
      m = ml + mu + 1
      info = 0
c
c     zero initial fill-in columns
c
      j0 = mu + 2
      j1 = MIN(n,m) - 1
      DO jz = j0, j1
         i0 = m + 1 - jz
         abd(i0:ml,jz) = zero
      END DO
      jz = j1
      ju = 0
c
c     gaussian elimination with partial pivoting
c
      nm1 = n - 1
      IF (nm1 .ge. 1) THEN
         DO k = 1, nm1
            kp1 = k + 1
c
c        zero next fill-in column
c
            jz = jz + 1
            IF (jz .le. n) THEN
               abd(:ml,jz) = zero
            END IF
c
c        find l = pivot INDEX
c
            lm = MIN(ml,n - k)
            isamax1 = MAXLOC(ABS(abd(m:m+lm,k)))
            l = isamax1(1) + m - 1
!           l = isamax(lm + 1,abd(m,k),1) + m - 1
            ipvt(k) = l + k - m
c
c        zero pivot implies this column already triangularized
c
            IF (abd(l,k) .ne. zero) THEN
c
c           interchange if necessary
c
               IF (l .ne. m) THEN
                  t = abd(l,k)
                  abd(l,k) = abd(m,k)
                  abd(m,k) = t
               END IF
c
c           compute multipliers
c
               t = -one/abd(m,k)
               CALL dscal (lm, t, abd(m+1,k), 1)
c
c           row elimination with column indexing
c
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
!DEC$ ENDIF
      SUBROUTINE sgbfa1 (ABD, LDA, N, ML, MU, IPVT, INFO)
      USE LIPREC, ONLY: WP => SP
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER LDA, N, ML, MU, INFO
      INTEGER, DIMENSION(N) :: IPVT
      REAL(WP), DIMENSION(*) :: ABD

      CALL sgbfa (ABD, LDA, N, ML, MU, IPVT, INFO)

      END SUBROUTINE sgbfa1 


      SUBROUTINE dgbfa1 (ABD, LDA, N, ML, MU, IPVT, INFO)
      USE LIPREC, ONLY: WP => DP
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER LDA, N, ML, MU, INFO
      INTEGER, DIMENSION(N) :: IPVT
      REAL(WP), DIMENSION(*) :: ABD

      CALL dgbfa (ABD, LDA, N, ML, MU, IPVT, INFO)

      END SUBROUTINE dgbfa1 

