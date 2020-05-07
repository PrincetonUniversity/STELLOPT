      SUBROUTINE svd_solve(m, n, mp, np, a, b, v, w, nw, small)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER m, n, mp, np, nw
      REAL(rprec), DIMENSION(mp,np) :: a
      REAL(rprec), DIMENSION(n,n) :: v
      REAL(rprec), DIMENSION(n) :: w
      REAL(rprec), DIMENSION(m) :: b
      REAL(rprec) :: small
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
!      REAL(rprec) :: small = 1.0e-10_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: istat, i, j
      INTEGER :: INFO, LWORK, LDA, LDU, LDVT
      REAL(rprec), ALLOCATABLE :: ap(:,:), WORK(:), u(:,:), VT(:,:)
      INTEGER, ALLOCATABLE :: IWORK(:)
      CHARACTER :: JOBZ
C-----------------------------------------------

c       Solves Matrix equation A * V = B for V using SVD method
c       Uses svd routines from numerical recipes
c       Prashant Valanju (Dec 1998) pvalanju@mail.utexas.edu

c       Almost same as SvdSolveB except this one finds
c       nw solutions V(:,i) by keeping i = 1 to nw vectors, where
c       nw is the number of vectors with nonzero weights > small

c       It also uses V internally for SVD decomposition to save space, so

c       Inputs:
c       A(M,N)   - Matrix with physical dimensions (Mp,Np)
c       B(M)     - R.H.S. of A X = B
c       M        - Number of rows of A to use
c       N        - Number of columns of A to use
c       Mp,Np    - Storage dimensions of A(Mp,Np)
c
c       Output:
c       nw      - number of vectors with normalized weights above small
c          Eventually, this will be the number of optimum weights
c          after we decide on a criterion for optimization
c       V(N,nw) - nw Solutions of A V = B, each a vector of length N
c          V(:N,iw) = solution with top iw weights kept
c          Physical dimensions of V are (N,N)
c       w(n)    - Weights in decreasing order
c

c  Allocate local arrays.
c  Note u(m,n) is enough because needed part of a(mp,np) is copied to
c  u(m,n), and a(mp,np) is never used directly. This saves space.
c  It is essential to USE a local u since svdcmp changes it

      ALLOCATE (ap(m,n),  stat=istat)
      IF(istat.ne.0) STOP 'Stop: No memory in svd_nesc'
      ALLOCATE (u(m,n),  stat=istat)
      IF(istat.ne.0) STOP 'Stop: No memory in svd_nesc'

c.......................................
c  Initialize ALL to zero to wipe out effects of old CALL
      w(:n) = 0                                !Zero ALL weights
      DO j = 1, n
         ap(:m,j) = a(:m,j)                  !Because U will be changed by svdcmp
         v(:n,j) = 0
      END DO

c  Do the SVD decomposition of a, i.e, of u into u, v and w
c      CALL svdcmp (u, m, n, m, n, w, v)
c  Sort weights and matrices with DECREASING weights so w(1) is biggest
c  Permute weight w(i) AND column vectors U(*,i), V(*,i) at the same time
c      CALL sortsvd (m, n, m, n, w, u, v)
      LDA = M
      LDU = M
      LDVT = N
      LWORK = max( 3*min(M,N) + max(max(M,N),7*min(M,N)), 
     1     3*min(M,N) + max(max(M,N),5*min(M,N)*min(M,N)+4*min(M,N)), 
     2     min(M,N)*(6+4*min(M,N))+max(M,N))
      JOBZ = 'S' ! min(M,N)
      ALLOCATE (WORK(LWORK),  stat=istat)
      IF(istat.ne.0) STOP 'Stop: No memory in svd_nesc'
      ALLOCATE (IWORK(8*min(M,N)),  stat=istat)
      IF(istat.ne.0) STOP 'Stop: No memory in svd_nesc'
      CALL DGESDD( JOBZ, M, N, AP, LDA, W, U, LDU, V, LDVT, WORK,
     1     LWORK, IWORK, INFO )
      IF(info.ne.0) print *, "Error in SVD (DGESDD): info=",INFO

c  Find nw = number of large weights (dcreasing ordered by sortsvd)
         DO nw = n, 1, -1          !Find first large weight and get out
            IF ( ABS(w(nw)/w(1)) .gt. small) EXIT
         END DO

c  Find nw solutions by successively adding up eigenvectors V
c     with correct weight given by (U(i) dot b)/w(i) in eq 14.3.17 (NR).
c  The coeff of ith vector V(i) is the dot product of the i-th column
c     of U with the rhs vector B, divided by w(i)
c  This does the svdbksb directly, faster than the NR svdbksub routine
c     and uses less memory due to the dual role of V
c  Note: any optimization scheme to find the 'best' nw will

c      First set the 1st vector with largest weight w(1)
       v(:n,1) = SUM(u(:m,1)*b(:m)) *v(:n,1) /w(1)
c      Next add the vectors with successive weights (in decreasing order)
         DO i = 2, nw
            j = i - 1
            v(:n,i) = v(:n,j) + SUM(u(:m,i)*b(:m)) *v(:n,i) / w(i)
         END DO
c................................................

      DEALLOCATE (u, stat=istat)
      DEALLOCATE (ap, stat=istat)
      DEALLOCATE (WORK, stat=istat)
      DEALLOCATE (IWORK, stat=istat)

      END SUBROUTINE svd_solve
