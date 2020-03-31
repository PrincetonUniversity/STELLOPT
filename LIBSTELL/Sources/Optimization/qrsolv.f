      SUBROUTINE qrsolv(n, r, ldr, ipvt, diag, qtb, x, sdiag, wa)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: n, ldr
      INTEGER, DIMENSION(n), INTENT(in) :: ipvt
      REAL(rprec), DIMENSION(ldr,n), INTENT(inout) :: r
      REAL(rprec), DIMENSION(n) :: diag, qtb, x, sdiag, wa
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: one = 1, zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: j, jp1, k, kp1, l, nsing, l1
      REAL(rprec) :: COS, cotan, qtbpj, SIN, SUM0, TAN, temp
      REAL(rprec) , ALLOCATABLE :: temp1u(:)
C-----------------------------------------------
c
c     SUBROUTINE qrsolv
c
c     given an m by n matrix a, an n by n diagonal matrix d,
c     and an m-vector b, the problem is to determine an x which
c     solves the system
c
c           a*x = b ,     d*x = 0 ,
c
c     in the least squares sense.
c
c     this SUBROUTINE completes the solution of the problem
c     IF it is provided with the necessary information from the
c     qr factorization, with column pivoting, of a. that is, IF
c     a*p = q*r, WHERE p is a permutation matrix, q has orthogonal
c     columns, and r is an upper triangular matrix with diagonal
c     elements of nonincreasing magnitude, THEN qrsolv EXPects
c     the full upper triangle of r, the permutation matrix p,
c     and the first n components of (q TRANSPOSE)*b. the system
c     a*x = b, d*x = 0, is THEN equivalent to
c
c                  t       t
c           r*z = q *b ,  p *d*p*z = 0 ,
c
c     where x = p*z. If this system does not have full rank,
c     then a least squares solution is obtained. on output qrsolv
c     also provides an upper triangular matrix s such that
c
c            t   t               t
c           p *(a *a + d*d)*p = s *s .
c
c     s is computed within qrsolv and may be of separate interest.
c
c     the SUBROUTINE statement is
c
c       SUBROUTINE qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)
c
c     WHERE
c
c       n is a positive INTEGER input variable set to the order of r.
c
c       r is an n by n array. on input the full upper triangle
c         must contain the full upper triangle of the matrix r.
c         on output the full upper triangle is unaltered, and the
c         strict lower triangle CONTAINS the strict upper triangle
c         (transposed) of the upper triangular matrix s.
c
c       ldr is a positive INTEGER input variable not less than n
c         which specifies the leading DIMENSION of the array r.
c
c       ipvt is an INTEGER input array of length n which defines the
c         permutation matrix p such that a*p = q*r. column j of p
c         is column ipvt(j) of the identity matrix.
c
c       diag is an input array of length n which must contain the
c         diagonal elements of the matrix d.
c
c       qtb is an input array of length n which must contain the first
c         n elements of the vector (q TRANSPOSE)*b.
c
c       x is an output array of length n which CONTAINS the least
c         squares solution of the system a*x = b, d*x = 0.
c
c       sdiag is an output array of length n which CONTAINS the
c         diagonal elements of the upper triangular matrix s.
c
c       wa is a work array of length n.
c
c     subprograms called
c
c       fortran-supplied ... ABS,sqrt
c
c     argonne national laboratory. MINpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
c
c     copy r and (q TRANSPOSE)*b to preserve input and initialize s.
c     in particular, SAVE the diagonal elements of r in x.
c
      DO j = 1, n
         r(j:n,j) = r(j,j:n)
         x(j) = r(j,j)
      END DO
      wa(1:n) = qtb(1:n)
c
c     eliminate the diagonal matrix d using a givens rotation.
c
      DO j = 1, n
c
c        prepare the row of d to be eliminated, locating the
c        diagonal element using p from the qr factorization.
c
         l = ipvt(j)
         IF (diag(l) .ne. zero) THEN
            sdiag(j+1:n) = zero
            sdiag(j) = diag(l)
c
c        the transformations to eliminate the row of d
c        modify ONLY a single element of (q TRANSPOSE)*b
c        beyond the first n, which is initially zero.
c
            qtbpj = zero
            DO k = j, n
c
c           determine a givens rotation which eliminates the
c           appropriate element in the current row of d.
c
               IF (sdiag(k) .ne. zero) THEN
                  IF (ABS(r(k,k)) .lt. ABS(sdiag(k))) THEN
                     coTAN = r(k,k)/sdiag(k)
                     SIN = one/SQRT(one + cotan*cotan)
                     COS = SIN*cotan
                  ELSE
                     TAN = sdiag(k)/r(k,k)
                     COS = one/SQRT(one + TAN*tan)
                     SIN = COS*tan
                  END IF
c
c           compute the modified diagonal element of r and
c           the modified element of ((q TRANSPOSE)*b,0).
c
                  r(k,k) = COS*r(k,k) + SIN*sdiag(k)
                  temp = COS*wa(k) + SIN*qtbpj
                  qtbpj = (-sin*wa(k)) + COS*qtbpj
                  wa(k) = temp
c
c           accumulate the tranformation in the row of s.
c
                  kp1 = k + 1
                  IF (n .ge. kp1) THEN
                     l1 = n-kp1+1
                     ALLOCATE (temp1u(l1))
                     temp1u(:l1) = COS*r(kp1:n,k) + SIN*sdiag(kp1:n)
                     sdiag(kp1:n) = (-sin*r(kp1:n,k)) + COS*sdiag(kp1:n)
                     r(kp1:n,k) = temp1u(:l1)
                     DEALLOCATE (temp1u)
                  END IF
               END IF
            END DO
         END IF
c
c        store the diagonal element of s and restore
c        the corresponding diagonal element of r.
c
         sdiag(j) = r(j,j)
         r(j,j) = x(j)
      END DO
c
c     solve the triangular system for z. IF the system is
c     singular, then obtain a least squares solution.
c
      nsing = n
      DO j = 1, n
         IF (sdiag(j).eq.zero .and. nsing.eq.n) nsing = j - 1
         IF (nsing .lt. n) wa(j) = zero
      END DO
      IF (nsing .ge. 1) THEN
         DO k = 1, nsing
            j = nsing - k + 1
            jp1 = j + 1
            IF (nsing .lt. jp1) THEN
               SUM0 = zero
            ELSE
               SUM0 = SUM(r(jp1:nsing,j)*wa(jp1:nsing))
            END IF
            wa(j) = (wa(j)-sum0)/sdiag(j)
         END DO
      END IF
c
c     permute the components of z back to components of x.
c
      DO j = 1, n
         l = ipvt(j)
         x(l) = wa(j)
      END DO

      END SUBROUTINE qrsolv
