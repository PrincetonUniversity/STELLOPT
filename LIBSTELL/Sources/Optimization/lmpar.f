      SUBROUTINE lmpar (n, r, ldr, ipvt, diag, qtb, delta, par, x,
     1                 sdiag, wa1, wa2)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: n, ldr, ipvt(n)
      REAL(rprec), INTENT(in) :: delta
      REAL(rprec), INTENT(inout) :: par
      REAL(rprec), DIMENSION(ldr,n), INTENT(inout) :: r
      REAL(rprec), DIMENSION(n) :: wa1, wa2
      REAL(rprec), DIMENSION(n), INTENT(in) :: diag, qtb
      REAL(rprec), DIMENSION(n), INTENT(out) :: x, sdiag
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: p1=0.1_dp, zero = 0,
     1     p001 = 0.001_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: iter, j, jm1, jp1, k, l, nsing
      REAL(rprec) :: dxnorm, dwarf, fp, gnorm, parc, parl, paru,
     1   sum0, enorm, dpmpar, temp, epsmch
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL dpmpar, enorm
C-----------------------------------------------
c
c     SUBROUTINE lmpar
c
c     given an m by n matrix a, an n by n nonsingular diagonal
c     matrix d, an m-vector b, and a positive number delta,
c     the problem is to determine a value for the parameter
c     par such that IF x solves the system
c
c           a*x = b ,     SQRT(par)*d*x = 0 ,
c
c     in the least squares sense, and dxnorm is the euclidean
c     norm of d*x, THEN either par is zero and
c
c           (dxnorm-delta) .le. 0.1*delta ,
c
c     or par is positive and
c
c           ABS(dxnorm-delta) .le. 0.1*delta .
c
c     this SUBROUTINE completes the solution of the problem
c     IF it is provided with the necessary information from the
c     qr factorization, with column pivoting, of a. that is, IF
c     a*p = q*r, WHERE p is a permutation matrix, q has orthogonal
c     columns, and r is an upper triangular matrix with diagonal
c     elements of nonincreasing magnitude, THEN lmpar EXPects
c     the full upper triangle of r, the permutation matrix p,
c     and the first n components of (q TRANSPOSE)*b. on output
c     lmpar also provides an upper triangular matrix s such that
c
c            t   t                   t
c           p *(a *a + par*d*d)*p = s *s .
c
c     s is employed within lmpar and may be of separate interest.
c
c     ONLY a few iterations are generally needed for convergence
c     of the algorithm. IF, however, the limit of 10 iterations
c     is reached, THEN the output par will contain the best
c     value obtained so far.
c
c     the SUBROUTINE statement is
c
c       SUBROUTINE lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,
c                        wa1,wa2)
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
c       delta is a positive input variable which specifies an upper
c         bound on the euclidean norm of d*x.
c
c       par is a nonnegative variable. on input par CONTAINS an
c         initial estimate of the levenberg-marquardt PARAMETER.
c         on output par CONTAINS the final estimate.
c
c       x is an output array of length n which CONTAINS the least
c         squares solution of the system a*x = b, SQRT(par)*d*x = 0,
c         for the output par.
c
c       sdiag is an output array of length n which CONTAINS the
c         diagonal elements of the upper triangular matrix s.
c
c       wa1 and wa2 are work arrays of length n.
c
c     subprograms called
c
c       MINpack-supplied ... dpmpar,enorm,qrsolv
c
c       fortran-supplied ... ABS,max,min,sqrt
c
c     argonne national laboratory. MINpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********

c
c     dwarf is the smallest positive magnitude.
c
      dwarf = dpmpar(2)
      epsmch = SQRT(dpmpar(1))       !SPH: Gives better lower bound for r(j,j)
c
c     compute and store in x the gauss-newton direction. IF the
c     jacobian is rank-deficient, obtain a least squares solution.
c
      nsing = n
      DO j = 1, n
         wa1(j) = qtb(j)
!SPH     IF (r(j,j).eq.zero .and. nsing.eq.n) nsing = j - 1
         IF (ABS(r(j,j)).le.epsmch .and. nsing.eq.n) nsing = j - 1
         IF (nsing .lt. n) wa1(j) = zero
      END DO
      IF (nsing .ge. 1) THEN
         DO k = 1, nsing
            j = nsing - k + 1
            wa1(j) = wa1(j)/r(j,j)
            temp = wa1(j)
            jm1 = j - 1
            IF (jm1 .ge. 1) THEN
               wa1(:jm1) = wa1(:jm1) - r(:jm1,j)*temp
            END IF
         END DO
      END IF
      DO j = 1, n
         l = ipvt(j)
         x(l) = wa1(j)
      END DO
c
c     initialize the iteration counter.
c     evaluate the FUNCTION at the origin, and test
c     for acceptance of the gauss-newton direction.
c
      wa2 = diag * x
      dxnorm = enorm(n,wa2)
      fp = dxnorm - delta

      IF (fp .le. p1*delta) THEN
         par = zero
         RETURN
      END IF

c
c     BEGIN GAUSS-NEWTON STEP
c
c     IF the jacobian is not rank deficient, the newton
c     step provides a lower bound, parl, for the zero of
c     the FUNCTION. otherwise set this bound to zero.
c
      parl = zero
      IF (nsing .ge. n) THEN
         wa1 = diag(ipvt)*(wa2(ipvt)/dxnorm)
         DO j = 1, n
            SUM0 = zero
            jm1 = j - 1
            IF (jm1 .ge. 1) SUM0 = SUM(r(:jm1,j)*wa1(:jm1))
            wa1(j) = (wa1(j)-sum0)/r(j,j)
         END DO
         temp = enorm(n,wa1)
         parl = ((fp/delta)/temp)/temp
      END IF
c
c     calculate an upper bound, paru, for the zero of the FUNCTION.
c
      DO j = 1, n
         SUM0 = SUM(r(:j,j)*qtb(:j))
         l = ipvt(j)
         wa1(j) = SUM0/diag(l)
      END DO
      gnorm = enorm(n,wa1)
      paru = gnorm/delta
      IF (paru .eq. zero) paru = dwarf/MIN(delta,p1)
c
c     IF the input par lies outside of the interval (parl,paru),
c     set par to the CLOSEr ENDpoint.
c
      par = MAX(par,parl)
      par = MIN(par,paru)
      IF (par .eq. zero) par = gnorm/dxnorm

c
c     beginning of an iteration loop.
c
      DO iter = 1, 10
c
c        evaluate the FUNCTION at the current value of par.
c
         IF (par .le. zero) par = MAX(dwarf,p001*paru)
         temp = SQRT(par)

         wa1 = temp*diag
         CALL qrsolv (n, r, ldr, ipvt, wa1, qtb, x, sdiag, wa2)
         wa2 = diag*x
         dxnorm = enorm(n,wa2)
         temp = fp
         fp = dxnorm - delta
c
c        If the function is small enough, accept the current value
c        of par. also test for the exceptional cases where parl
c        is zero or the number of iterations has reached 10.
c
         IF (ABS(fp).le.p1*delta .or. parl.eq.zero .and.
     1       fp.le.temp .and. temp.lt.zero .or. iter.eq.10) EXIT
c
c        compute the newton correction.
c
         wa1 = diag(ipvt)*(wa2(ipvt)/dxnorm)
         DO j = 1, n
            wa1(j) = wa1(j)/sdiag(j)
            temp = wa1(j)
            jp1 = j + 1
            IF (n .ge. jp1) wa1(jp1:n) = wa1(jp1:n) - r(jp1:n,j)*temp
         END DO
         temp = enorm(n,wa1)
         parc = ((fp/delta)/temp)/temp
c
c        depending on the sign of the FUNCTION, update parl or paru.
c
         IF (fp .gt. zero) parl = MAX(parl,par)
         IF (fp .lt. zero) paru = MIN(paru,par)
c
c        compute an improved estimate for par.
c
         par = MAX(parl,par + parc)

      END DO                        !!end of an iteration.

      END SUBROUTINE lmpar
