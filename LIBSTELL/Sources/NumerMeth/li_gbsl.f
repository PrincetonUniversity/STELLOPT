      SUBROUTINE sgbsl1 (ABD, LDA, N, ML, MU, IPVT, B, JOB)
      USE LIPREC, ONLY: WP => SP
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER LDA, N, ML, MU, JOB
      INTEGER, DIMENSION(N) :: IPVT
      REAL(WP), DIMENSION(*) :: ABD
      REAL(WP), DIMENSION(N) :: B

      CALL sgbsl (ABD, LDA, N, ML, MU, IPVT, B, JOB)

      END SUBROUTINE sgbsl1

      SUBROUTINE dgbsl1 (ABD, LDA, N, ML, MU, IPVT, B, JOB)
      USE LIPREC, ONLY: WP => DP
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER LDA, N, ML, MU, JOB
      INTEGER, DIMENSION(N) :: IPVT
      REAL(WP), DIMENSION(*) :: ABD
      REAL(WP), DIMENSION(N) :: B

      CALL dgbsl (ABD, LDA, N, ML, MU, IPVT, B, JOB)

      END SUBROUTINE dgbsl1

#ifndef CRAY
      SUBROUTINE sgbsl (ABD, LDA, N, ML, MU, IPVT, B, JOB)
      USE LIPREC, ONLY: WP => SP
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER LDA, N, ML, MU, JOB
      INTEGER, DIMENSION(N) :: IPVT
      REAL(WP), DIMENSION(LDA,N) :: ABD
      REAL(WP), DIMENSION(N) :: B
C-----------------------------------------------
C   L O C A L   V A R I A B L E S
C-----------------------------------------------
      INTEGER :: K, KB, L, LA, LB, LM, M, NM1
      REAL(WP) :: T
C-----------------------------------------------
C   E X T E R N A L   F U N C T I O N S
C-----------------------------------------------
      REAL(WP) , EXTERNAL :: sdot
      EXTERNAL saxpy
C-----------------------------------------------
c
c     sgbsl solves the real band system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by sgbco or sgbfa.
c
c     on ENTRY
c
c        abd     REAL(lda, n)
c                the output from sgbco or sgbfa.
c
c        lda     INTEGER
c                the leading DIMENSION of the array  abd .
c
c        n       INTEGER
c                the order of the original matrix.
c
c        ml      INTEGER
c                number of diagonals below the main diagonal.
c
c        mu      INTEGER
c                number of diagonals above the main diagonal.
c
c        ipvt    INTEGER(n)
c                the pivot vector from sgbco or sgbfa.
c
c        b       REAL(n)
c                the right hand side vector.
c
c        job     INTEGER
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b , WHERE
c                            trans(a)  is the TRANSPOSE.
c
c     on RETURN
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if sgbco has set rcond .gt. 0.0
c        or sgbfa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call sgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call sgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas saxpy ,sdot
c     fortran min
c
c     internal variables
c
c
      m = mu + ml + 1
      nm1 = n - 1
      IF (job .eq. 0) THEN
c
c        job = 0 , solve  a * x = b
c        first solve l*y = b
c
         IF (ml .ne. 0) THEN
            DO k = 1, nm1
               lm = MIN(ml,n - k)
               l = ipvt(k)
               t = b(l)
               IF (l .ne. k) THEN
                  b(l) = b(k)
                  b(k) = t
               END IF
               CALL saxpy (lm, t, abd(m+1,k), 1, b(k+1), 1)
            END DO
         END IF
c
c        now solve  u*x = y
c
         DO kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/abd(m,k)
            lm = MIN(k,m) - 1
            la = m - lm
            lb = k - lm
            t = -b(k)
            CALL saxpy (lm, t, abd(la,k), 1, b(lb), 1)
         END DO
      ELSE
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         DO k = 1, n
            lm = MIN(k,m) - 1
            la = m - lm
            lb = k - lm
            t = sdot(lm,abd(la,k),1,b(lb),1)
            b(k) = (b(k)-t)/abd(m,k)
         END DO
c
c        now solve trans(l)*x = y
c
         IF (ml .ne. 0) THEN
            DO kb = 1, nm1
               k = n - kb
               lm = MIN(ml,n - k)
               b(k) = b(k) + sdot(lm,abd(m+1,k),1,b(k+1),1)
               l = ipvt(k)
               IF (l .ne. k) THEN
                  t = b(l)
                  b(l) = b(k)
                  b(k) = t
               END IF
            END DO
         END IF
      END IF

      END SUBROUTINE sgbsl

      SUBROUTINE dgbsl (ABD, LDA, N, ML, MU, IPVT, B, JOB)
      USE LIPREC, ONLY: WP => DP
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER LDA, N, ML, MU, JOB
      INTEGER, DIMENSION(N) :: IPVT
      REAL(WP), DIMENSION(LDA,N) :: ABD
      REAL(WP), DIMENSION(N) :: B
C-----------------------------------------------
C   L O C A L   V A R I A B L E S
C-----------------------------------------------
      INTEGER :: K, KB, L, LA, LB, LM, M, NM1
      REAL(WP) :: T
C-----------------------------------------------
C   E X T E R N A L   F U N C T I O N S
C-----------------------------------------------
      REAL(WP) , EXTERNAL :: ddot
      EXTERNAL daxpy
C-----------------------------------------------
c
c     sgbsl solves the REAL band system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by sgbco or sgbfa.
c
c     on ENTRY
c
c        abd     REAL(lda, n)
c                the output from sgbco or sgbfa.
c
c        lda     INTEGER
c                the leading DIMENSION of the array  abd .
c
c        n       INTEGER
c                the order of the original matrix.
c
c        ml      INTEGER
c                number of diagonals below the main diagonal.
c
c        mu      INTEGER
c                number of diagonals above the main diagonal.
c
c        ipvt    INTEGER(n)
c                the pivot vector from sgbco or sgbfa.
c
c        b       REAL(n)
c                the right hand side vector.
c
c        job     INTEGER
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b , WHERE
c                            trans(a)  is the TRANSPOSE.
c
c     on RETURN
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if sgbco has set rcond .gt. 0.0
c        or sgbfa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call sgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call sgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy ,ddot
c     fortran min
c
c     internal variables
c
c
      m = mu + ml + 1
      nm1 = n - 1
      IF (job .eq. 0) THEN
c
c        job = 0 , solve  a * x = b
c        first solve l*y = b
c
         IF (ml .ne. 0) THEN
            DO k = 1, nm1
               lm = MIN(ml,n - k)
               l = ipvt(k)
               t = b(l)
               IF (l .ne. k) THEN
                  b(l) = b(k)
                  b(k) = t
               END IF
               CALL daxpy (lm, t, abd(m+1,k), 1, b(k+1), 1)
            END DO
         END IF
c
c        now solve  u*x = y
c
         DO kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/abd(m,k)
            lm = MIN(k,m) - 1
            la = m - lm
            lb = k - lm
            t = -b(k)
            CALL daxpy (lm, t, abd(la,k), 1, b(lb), 1)
         END DO
      ELSE
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         DO k = 1, n
            lm = MIN(k,m) - 1
            la = m - lm
            lb = k - lm
            t = ddot(lm,abd(la,k),1,b(lb),1)
            b(k) = (b(k)-t)/abd(m,k)
         END DO
c
c        now solve trans(l)*x = y
c
         IF (ml .ne. 0) THEN
            DO kb = 1, nm1
               k = n - kb
               lm = MIN(ml,n - k)
               b(k) = b(k) + ddot(lm,abd(m+1,k),1,b(k+1),1)
               l = ipvt(k)
               IF (l .ne. k) THEN
                  t = b(l)
                  b(l) = b(k)
                  b(k) = t
               END IF
            END DO
         END IF
      END IF

      END SUBROUTINE dgbsl
#endif
