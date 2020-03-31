      SUBROUTINE tridslv(a, d, b, c, jmin, jmax, mnd1, ns, nrhs)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: jmax, mnd1, ns, nrhs
      INTEGER, DIMENSION(0:mnd1), INTENT(in) :: jmin
      REAL(rprec), DIMENSION(ns,0:mnd1) :: a, d, b
      REAL(rprec), DIMENSION(ns,0:mnd1, nrhs), INTENT(inout) :: c
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: mn, in, i0, in1, jrhs
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: alf
      REAL(rprec), DIMENSION(0:mnd1) :: psi0
C-----------------------------------------------
!
!     SOLVES B(I)*X(I-1)+D(I)*X(I)+A(I)*X(I+1)=C(I), I=IN,JMAX
!     AND RETURNS ANSWER IN C(I)
!     ADDED VECTORIZATION ON FOURIER MODE ARGUMENT (01-2000)
!     AND NEW ARGUMENT (NRHS) TO DO MULTIPLE RIGHT SIDES SIMULTANEOUSLY
!
      IF (jmax .gt. ns) STOP 'jmax>ns in tridslv'

      ALLOCATE (alf(ns,0:mnd1), stat = in)
      IF (in .ne. 0) STOP 'Allocation error in tridslv'

      in = MINVAL(jmin)
!
!      FILL IN MN BELOW MAX(JMIN) WITH DUMMY VALUES
!      TO ALLOW VECTORIZATION ON MN INDEX
!
      DO mn = 0, mnd1
         in1 = jmin(mn)-1
         IF (in1 .ge. in) THEN
            d(in:in1, mn) = 1
            c(in:in1, mn, 1:nrhs) = 0
            b(in:in1, mn) = 0
            a(in:in1, mn) = 0
         END IF
      END DO

      in1 = in + 1

      psi0(:)= d(in,:)
      IF (ANY(psi0 .eq. zero)) STOP 'psi0 == 0 error in tridslv'
      psi0 = one/psi0
      DO jrhs = 1, nrhs
         c(in,:,jrhs) = c(in,:,jrhs)*psi0(:)
      END DO

      DO i0 = in1,jmax
         alf(i0-1,:) = a(i0-1,:)*psi0(:)
         psi0 = d(i0,:) - b(i0,:)*alf(i0-1,:)
         IF (ANY(ABS(psi0) .le. 1.E-8_dp*ABS(d(i0,:)))) 
     1       STOP 'psi0/d(i0) < 1.E-8: possible singularity in tridslv'
         psi0  = one/psi0
         DO jrhs = 1, nrhs
            c(i0,:,jrhs) = (c(i0,:,jrhs) - b(i0,:)*c(i0-1,:,jrhs))
     1                   * psi0
         END DO
      END DO

      DO i0 = jmax - 1, in, -1
         DO jrhs = 1,nrhs
            c(i0,:,jrhs) = c(i0,:,jrhs) - alf(i0,:)*c(i0+1,:,jrhs)
         END DO
      END DO

      DEALLOCATE (alf)

      END SUBROUTINE tridslv

!
!     ADDED 10/3/14 TO HANDLE DATA-LAYOUT IN PARALLEL VERSION OF VMEC2000
      SUBROUTINE tridslv_par(a, d, b, c, jmin, jmax, mnd1, ns, nrhs)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: jmax, mnd1, ns, nrhs
      INTEGER, DIMENSION(0:mnd1), INTENT(in) :: jmin
      REAL(rprec), DIMENSION(0:mnd1,ns) :: a, d, b
      REAL(rprec), DIMENSION(0:mnd1,ns,nrhs), INTENT(inout) :: c
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: mn, in, i0, in1, jrhs
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: alf
      REAL(rprec), DIMENSION(0:mnd1) :: psi0
C-----------------------------------------------
!
!     SOLVES B(I)*X(I-1)+D(I)*X(I)+A(I)*X(I+1)=C(I), I=IN,JMAX
!     AND RETURNS ANSWER IN C(I)
!     ADDED VECTORIZATION ON FOURIER MODE ARGUMENT (01-2000)
!     AND NEW ARGUMENT (NRHS) TO DO MULTIPLE RIGHT SIDES SIMULTANEOUSLY
!
      IF (jmax .gt. ns) STOP 'jmax>ns in tridslv_par'

      ALLOCATE (alf(0:mnd1,ns), stat=in)
      IF (in .ne. 0) STOP 'Allocation error in tridslv_par'

      in = MINVAL(jmin)
!
!      FILL IN MN BELOW MAX(JMIN) WITH DUMMY VALUES
!      TO ALLOW VECTORIZATION ON MN INDEX
!
      DO mn = 0, mnd1
         in1 = jmin(mn)-1
         IF (in1 .ge. in) THEN
            d(mn, in:in1) = 1
            c(mn, in:in1, 1:nrhs) = 0
            b(mn, in:in1) = 0
            a(mn, in:in1) = 0
         END IF
      END DO

      in1 = in + 1

      psi0(:)= d(:,in)
      IF (ANY(psi0 .eq. zero)) STOP 'psi0 == 0 error in tridslv_par'
      psi0 = one/psi0
      DO jrhs = 1, nrhs
         c(:,in,jrhs) = c(:,in,jrhs)*psi0(:)
      END DO

      DO i0 = in1,jmax
         alf(:,i0-1) = a(:,i0-1)*psi0(:)
         psi0 = d(:,i0) - b(:,i0)*alf(:,i0-1)
         IF (ANY(ABS(psi0) .le. 1.E-8_dp*ABS(d(:,i0)))) 
     1       STOP 'psi0/d(i0) < 1.E-8: possible singularity in tridslv'
         psi0  = one/psi0
         DO jrhs = 1, nrhs
            c(:,i0,jrhs) = (c(:,i0,jrhs) - b(:,i0)*c(:,i0-1,jrhs))
     1                   * psi0
         END DO
      END DO

      DO i0 = jmax - 1, in, -1
         DO jrhs = 1,nrhs
            c(:,i0,jrhs) = c(:,i0,jrhs) - alf(:,i0)*c(:,i0+1,jrhs)
         END DO
      END DO

      DEALLOCATE (alf)

      END SUBROUTINE tridslv_par

