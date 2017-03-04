!DEC$ IF .NOT.DEFINED (CRAY)
      SUBROUTINE sgefa (A, LDA, N, IPVT, INFO)
      USE LIPREC, ONLY: WP => SP
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER LDA, N, INFO
      INTEGER, DIMENSION(N) :: IPVT
      REAL(WP), DIMENSION(LDA,N) :: A
C-----------------------------------------------
C   L O C A L   V A R I A B L E S
C-----------------------------------------------
      REAL(WP), PARAMETER :: ZERO = 0, ONE = 1
      INTEGER :: J, K, KP1, L, NM1
      INTEGER, DIMENSION(1) :: ISAMAX
      REAL(WP) :: ELEMENT
C-----------------------------------------------
c
c     sgefa factors a real matrix by gaussian elimination.
c
c     sgefa is usually called by sgeco, but it can be called
c     directly with a saving in time if rcond  is not needed.
c     (time for sgeco) = (1 + 9/n)*(time for sgefa) .
c
c     on ENTRY
c
c        a       REAL(lda, n)
c                the matrix to be factored.
c
c        lda     INTEGER
c                the leading dimension of the array  a .
c
c        n       INTEGER
c                the order of the matrix  a .
c
c     on RETURN
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    INTEGER(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k).eq.0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that sgesl or sgedi will divide by zero
c                     if called.  use  rcond  in sgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      IF (nm1 .ge. 1) THEN
         DO k = 1, nm1
            kp1 = k + 1
c
c         find l = pivot index
c
            isamax = MAXLOC(ABS(a(k:n,k)))
            l = isamax(1) + k - 1
C           l = IDAMAX(n-k+1,a(k,k),1) + k - 1
            ipvt(k) = l
c

c         zero pivot implies this column already triangularized
c
            IF (a(l,k) .ne. zero) THEN
c
c           interchange if necessary
c
               IF (l .ne. k) THEN
                  element = a(l,k)
                  a(l,k) = a(k,k)
                  a(k,k) = element
               END IF
c
c         compute multipliers
c
               element = -one/a(k,k)
!              CALL bla_scal(n-k,element,a(k+1,k),1)
               a(k+1:n,k) = element*a(k+1:n,k)
c
c           row elimination with column indexing
c
               DO j = kp1, n
                 element = a(l,j)
                 IF (l .ne. k) THEN
                    a(l,j) = a(k,j)
                    a(k,j) = element
                 END IF
!                CALL bla_axpy(n-k,element,a(k+1,k),1,a(k+1,j),1)
                 a(k+1:n,j) = a(k+1:n,j) + element*a(k+1:n,k)
               END DO

            ELSE
               info = k
            END IF
         END DO
      END IF
      ipvt(n) = n
      IF (a(n,n) .eq. zero) info = n

      END SUBROUTINE sgefa

      SUBROUTINE dgefa (A, LDA, N, IPVT, INFO)
      USE LIPREC, ONLY: WP => DP
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER LDA, N, INFO
      INTEGER, DIMENSION(N) :: IPVT
      REAL(WP), DIMENSION(LDA,N) :: A
C-----------------------------------------------
C   L O C A L   V A R I A B L E S
C-----------------------------------------------
      REAL(WP), PARAMETER :: ZERO = 0, ONE = 1
      INTEGER :: J, K, KP1, L, NM1
      INTEGER, DIMENSION(1) :: ISAMAX
      REAL(WP) :: ELEMENT
C-----------------------------------------------
c
c     sgefa factors a REAL matrix by gaussian elimination.
c
c     sgefa is usually called by sgeco, but it can be called
c     directly with a saving in time IF  rcond  is not needed.
c     (time for sgeco) = (1 + 9/n)*(time for sgefa) .
c
c     on ENTRY
c
c        a       REAL(lda, n)
c                the matrix to be factored.
c
c        lda     INTEGER
c                the leading DIMENSION of the array  a .
c
c        n       INTEGER
c                the order of the matrix  a .
c
c     on RETURN
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  WHERE
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    INTEGER(n)
c                an INTEGER vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k).eq.0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that sgesl or sgedi will divide by zero
c                     if called.  use  rcond  in sgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      IF (nm1 .ge. 1) THEN
         DO k = 1, nm1
            kp1 = k + 1
c
c         find l = pivot INDEX
c
            isamax = MAXLOC(ABS(a(k:n,k)))
            l = isamax(1) + k - 1
C           l = IDAMAX(n-k+1,a(k,k),1) + k - 1
            ipvt(k) = l
c

c         zero pivot implies this column already triangularized
c
            IF (a(l,k) .ne. zero) THEN
c
c           interchange if necessary
c
               IF (l .ne. k) THEN
                  element = a(l,k)
                  a(l,k) = a(k,k)
                  a(k,k) = element
               END IF
c
c         compute multipliers
c
               element = -one/a(k,k)
!              CALL bla_scal(n-k,element,a(k+1,k),1)
               a(k+1:n,k) = element*a(k+1:n,k)
c
c           row elimination with column indexing
c
               DO j = kp1, n
                 element = a(l,j)
                 IF (l .ne. k) THEN
                    a(l,j) = a(k,j)
                    a(k,j) = element
                 END IF
!                CALL bla_axpy(n-k,element,a(k+1,k),1,a(k+1,j),1)
                 a(k+1:n,j) = a(k+1:n,j) + element*a(k+1:n,k)
               END DO

            ELSE
               info = k
            END IF
         END DO
      END IF
      ipvt(n) = n
      IF (a(n,n) .eq. zero) info = n

      END SUBROUTINE dgefa
!DEC$ ENDIF
      SUBROUTINE sgefa1 (A, LDA, N, IPVT, INFO)
      USE LIPREC, ONLY: WP => SP
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER LDA, N, INFO
      INTEGER, DIMENSION(N) :: IPVT
      REAL(WP), DIMENSION(LDA,N) :: A

      CALL sgefa (A, LDA, N, IPVT, INFO)

      END SUBROUTINE sgefa1 

      SUBROUTINE dgefa1 (A, LDA, N, IPVT, INFO)
      USE LIPREC, ONLY: WP => DP
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER LDA, N, INFO
      INTEGER, DIMENSION(N) :: IPVT
      REAL(WP), DIMENSION(LDA,N) :: A

      CALL dgefa (A, LDA, N, IPVT, INFO)

      END SUBROUTINE dgefa1
