      SUBROUTINE qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: m, n, lda, lipvt
      LOGICAL, INTENT(in) :: pivot
      INTEGER, DIMENSION(lipvt), INTENT(out) :: ipvt
      REAL(rprec), DIMENSION(lda,n), INTENT(inout) :: a
      REAL(rprec), DIMENSION(n), INTENT(out) :: rdiag, acnorm, wa
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: one = 1, p05 = 0.05_dp,
     1    zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: j, jp1, k, kmax, MINmn
      REAL(rprec) :: ajnorm, epsmch, SUM0, temp, enorm, dpmpar
      REAL(rprec) , ALLOCATABLE :: temp1u(:)
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL dpmpar,enorm
C-----------------------------------------------
c
c     SUBROUTINE qrfac
c
c     this SUBROUTINE uses householder transformations with column
c     pivoting (optional) to compute a qr factorization of the
c     m by n matrix a. that is, qrfac determines an orthogonal
c     matrix q, a permutation matrix p, and an upper trapezoidal
c     matrix r with diagonal elements of nonincreasing magnitude,
c     such that a*p = q*r. the householder transformation for
c     column k, k = 1,2,...,MIN(m,n), is of the form
c
c                           t
c           i - (1/u(k))*u*u
c
c     WHERE u has zeros in the first k-1 positions. the form of
c     this transformation and the method of pivoting first
c     appeared in the corresponding linpack SUBROUTINE.
c
c     the SUBROUTINE statement is
c
c       SUBROUTINE qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
c
c     WHERE
c
c       m is a positive INTEGER input variable set to the number
c         of rows of a.
c
c       n is a positive INTEGER input variable set to the number
c         of columns of a.
c
c       a is an m by n array. on input a CONTAINS the matrix for
c         which the qr factorization is to be computed. on output
c         the strict upper trapezoidal part of a CONTAINS the strict
c         upper trapezoidal part of r, and the lower trapezoidal
c         part of a CONTAINS a factored form of q (the non-trivial
c         elements of the u vectors described above).
c
c       lda is a positive INTEGER input variable not less than m
c         which specifies the leading DIMENSION of the array a.
c
c       pivot is a LOGICAL input variable. IF pivot is set true,
c         THEN column pivoting is enforced. IF pivot is set false,
c         THEN no column pivoting is done.
c
c       ipvt is an INTEGER output array of length lipvt. ipvt
c         defines the permutation matrix p such that a*p = q*r.
c         column j of p is column ipvt(j) of the identity matrix.
c         IF pivot is false, ipvt is not referenced.
c
c       lipvt is a positive integer input variable. if pivot is false,
c         then lipvt may be as small as 1. if pivot is true, then
c         lipvt must be at least n.
c
c       rdiag is an output array of length n which CONTAINS the
c         diagonal elements of r.
c
c       acnorm is an output array of length n which CONTAINS the
c         norms of the corresponding columns of the input matrix a.
c         IF this information is not needed, THEN acnorm can coincide
c         with rdiag.
c
c       wa is a work array of length n. IF pivot is false, THEN wa
c         can coincide with rdiag.
c
c     subprograms called
c
c       MINpack-supplied ... dpmpar,enorm
c
c       fortran-supplied ... MAX,sqrt,min
c
c     argonne national laboratory. MINpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
c
c     epsmch is the machine precision.
c
      epsmch = dpmpar(1)
c
c     compute the initial column norms and initialize several arrays.
c
      DO j = 1, n
         acnorm(j) = enorm(m,a(1,j))
         rdiag(j) = acnorm(j)
         wa(j) = rdiag(j)
         IF (pivot) ipvt(j) = j
      END DO
c
c     reduce a to r with householder transformations.
c
      MINmn = MIN(m,n)
      ALLOCATE (temp1u(m))

      DO j = 1, MINmn
         IF (pivot) THEN
c
c        bring the column of largest norm into the pivot position.
c
            kMAX = j
            DO k = j,n
               IF (rdiag(k) .gt. rdiag(kmax)) kMAX = k
            END DO
            IF (kMAX .ne. j) THEN
               temp1u(:m) = a(:m,j)
               a(:m,j) = a(:m,kmax)
               a(:m,kmax) = temp1u(:m)
               rdiag(kmax) = rdiag(j)
               wa(kmax) = wa(j)
               k = ipvt(j)
               ipvt(j) = ipvt(kmax)
               ipvt(kmax) = k
            END IF
         END IF
c
c        compute the householder transformation to reduce the
c        j-th column of a to a multiple of the j-th unit vector.
c
         ajnorm = enorm(m-j+1,a(j,j))
         IF (ajnorm .ne. zero) THEN
            IF (a(j,j) .lt. zero) ajnorm = -ajnorm
            a(j:m,j) = a(j:m,j)/ajnorm
            a(j,j) = a(j,j) + one
c
c        apply the transformation to the remaining columns
c        and update the norms.
c
            jp1 = j + 1
            IF (n .ge. jp1) THEN
               DO k = jp1, n
                  SUM0 = SUM(a(j:m,j)*a(j:m,k))
                  temp = SUM0/a(j,j)
                  a(j:m,k) = a(j:m,k) - temp*a(j:m,j)
                  IF (pivot .and. rdiag(k).ne.zero) THEN
                     temp = a(j,k)/rdiag(k)
                     rdiag(k) = rdiag(k)*SQRT(MAX(zero,one-temp**2))
                     IF (p05*(rdiag(k)/wa(k))**2 .le. epsmch) THEN
                        rdiag(k) = enorm(m - j,a(jp1,k))
                        wa(k) = rdiag(k)
                     END IF
                  END IF
               END DO
            END IF
         END IF
         rdiag(j) = -ajnorm
      END DO

      DEALLOCATE (temp1u)

      END SUBROUTINE qrfac
