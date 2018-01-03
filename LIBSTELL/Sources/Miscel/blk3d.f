      SUBROUTINE blk3d(a, bm1, bp1, srces, mblk, nblocks)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE stel_kinds
      USE safe_open_mod
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: bytes_per_rprec = 8
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nblocks, mblk
      REAL(rprec), DIMENSION(mblk,mblk,nblocks), INTENT(in) :: 
     1                       a, bp1, bm1
      REAL(rprec), DIMENSION(mblk,nblocks), INTENT(inout) :: srces
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
!      INTEGER :: ibuph, incnow, irecl, incbu, iunit=102, ndisk
      INTEGER :: k, k1, ier
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ipiv
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: ainv
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: ql
C-----------------------------------------------
C   E x t e r n a l  S u b r o u t i n e s
C-----------------------------------------------
      EXTERNAL la_getrf, la_getrs
C-----------------------------------------------
c  modified (June, 2003, ORNL):         S. P. Hirshman
c-----------------------------------------------------------------------
c
c  this subroutine solves a block-tridiagonal system of equations.
c
c-----------------------------------------------------------------------
c  INPUT
c  mblk                : block size
c  nblocks             : number of blocks
c  a                   : diagonal block
c  bm1, bp1            : lower, upper block (see equation below)
c  srces(mblk,nblocks) : right-hand size source
c
c  LOCAL VARIABLES
c  iunit               : unit number for block-tridiagonal solution disk file.
c
c  solutions are indexed in m-n fourier-space, legendre-space. the tri-diagonal
c  equation is:
c
c           bm1 * f(l-1) + a * f(l) + bp1 * f(l+1) = source(l)
c
c     GENERAL SOLUTION SCHEME APPLIED TO EACH BLOCK ROW (INDEX L)
c
c     1. Start from row N and solve for x(N) in terms of x(N-1):
c
c        x(N) = -q(N)*x(N-1) + r(N)
c
c        q(N) =  a(N)[-1] * bm1;    r(N) = a(N)[-1] * s(N)
c
c        where a(N)[-1] is the inverse of a(N)
c
c     2. Substitute for lth row to get recursion equation fo q(l) and r(l):
c
c        x(l) = -q(l)*x(l-1) + r(l), in general, where:
c
c        q(l) = (a(l) - bp1(l)*q(l+1))[-1] * bm1(l)
c
c        r(l) = (a(l) - bp1(l)*q(l+1))[-1] * (s(l) - bp1(l)*r(l+1))
c
c     3. At row l = 1, bm1(1) = 0 and get an equation for x(1) corresponding to q(1) = 0:
c
c        x(1) = r(1)
c
c     4. Finally, can back-solve for x(N) in terms of x(N-1) from eqn.(1) above
c
c
c     NUMERICAL IMPLEMENTATION (USING LAPACK, la_... ROUTINES)
c
c     1. CALL la_getrf:   Perform LU factorization of diagonal block (A)
c     2. CALL la_getrs:   With multiple (mblk) right-hand sides, to do block inversion
c                         operation, A X = B  (store result in X; here B is a matrix)
c     3. CALL la_getrs:   With single right hand side (source) to solve A x = b (b a vector)
c
!      ndisk = mblk

      ALLOCATE (ql(mblk,mblk, nblocks), ainv(mblk,mblk), ipiv(mblk),
     1          stat=ier)
      IF (ier .ne. 0) STOP 'Allocation error in blk3d'

c  create disk file for doing direct access i/o.

!      incnow = ndisk
!      irecl  = bytes_per_rprec*incnow
!      incbu = 1 + (ndisk - 1)/incnow
!      ibuph = 0

!      iunit = 10
!      CALL safe_open(iunit, ier, 'NULL', 'scratch', 'unformatted',
!     1     irecl, 'DIRECT')
!      IF (ier .ne. 0) STOP 'Error opening scratch file in blk3d'

      ainv = a(:,:,nblocks)

c  main loop. load and process (backwards) block-rows nblocks to 1. 

      BLOCKS: DO k = nblocks, 1, -1
!
!     Compute (and save) ql(nblocks) = a(nblocks)[-1] * ql, and source terms A-1 * srces [A-1 == inv(A)]
!
         CALL la_getrf (mblk, mblk, ainv, mblk, ipiv, ier)
         IF (ier .ne. 0) GOTO 200

         IF (k .gt. 1) THEN
            ql(:,:,k) = bm1(:,:,k)
            CALL la_getrs('n', mblk, mblk, ainv,
     1                    mblk, ipiv, ql(1,1,k), mblk, ier)
            CALL la_getrs('n', mblk, 1,    ainv,
     1                    mblk,ipiv,srces(1,k),mblk,ier)

!            CALL wrdisk(iunit, ql, ndisk, incnow, ibuph, incbu, ier)
!            IF (ier .ne. 0) GOTO 302
         ELSE
            CALL la_getrs('n', mblk, 1,   ainv,
     1                    mblk,ipiv,srces(1,k),mblk,ier)
         END IF

         k1 = k - 1
         IF (k1 .eq. 0) EXIT

!
!      Update effective diagonal "a" matrix and source terms and store ql 2 l-steps back
!
         ainv  = a(:,:,k1) - MATMUL(bp1(:,:,k1), ql(:,:,k))
         srces(:,k1) = srces(:,k1) - MATMUL(bp1(:,:,k1),srces(:,k))

      END DO BLOCKS

c  backward solution sweep for block-rows k = 2 to nblocks
c  read blocks ql from disk.

      DO k = 2, nblocks
!         CALL rddisk (iunit, ql, ndisk, incnow, ibuph, ier)
!         IF (ier .ne. 0) GOTO 303
!         ibuph = ibuph - incbu
         srces(:,k) = srces(:,k) - MATMUL(ql(:,:,k),srces(:,k-1))
      END DO

      GOTO 400

c  error returns. ------------------------------------------------------

  200 CONTINUE
      WRITE (6, '(a,i4)') ' Error factoring matrix in blk3d: block = '
     1                     , k,' error id = ', ier
      GOTO 400
  301 CONTINUE
!      WRITE (6, '(a,i8)') ' BLK3D:   error in opening file:  ',
!     1   'RECL = ', irecl
  302 CONTINUE
      WRITE (6, '(a)') ' BLK3D:   error in I/O routine WRDISK'
  303 CONTINUE
      WRITE (6, '(a)') ' BLK3D:   error in I/O routine RDDISK'
      ier = -2
  305 CONTINUE
      WRITE (6, '(2/a,i4,2/)') ' BLK3D:   error detected:   ier =',
     1   ier
      STOP

c  destroy disk file and return. ---------------------------------------

  400 CONTINUE

!      CLOSE (iunit)

      DEALLOCATE (ainv, ql, ipiv)

      END SUBROUTINE blk3d
