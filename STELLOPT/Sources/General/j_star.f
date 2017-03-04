      SUBROUTINE j_star(modb, bmin, bmax, ep_mu, jstar, nzeta, ntheta)
      USE stel_kinds
      USE mpi_params                                                     !MPI
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER nzeta, ntheta
      REAL(rprec), INTENT(in) :: ep_mu
      REAL(rprec), DIMENSION(nzeta,ntheta), INTENT(in) :: modb
      REAL(rprec), DIMENSION(ntheta), INTENT(in) :: bmin, bmax
      REAL(rprec), DIMENSION(ntheta), INTENT(out) :: jstar
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, one = 1
      INTEGER :: ku_lw, ku_up, ku, istat
      REAL(rprec) :: dzeta
      REAL(rprec), DIMENSION(ntheta) :: test_bmin, test_bmax
      REAL(rprec) , DIMENSION(ntheta) :: test1, test2
      REAL(rprec) , ALLOCATABLE,  DIMENSION(:,:) :: vpar1
      LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: l3v
C-----------------------------------------------
!
!     computes the trapped branch of jstar on a single flux surface
!     and for a single value of ep/mu.  jstar is only non-zero for
!     values of theta where trapped particle orbits can exist.  at
!     theta values which are in the passing particle regime (ep/mu > bmax)
!     or the forbidden regime (ep/mu < bmin) jstar is set to 0.  a jstar
!     topology is assumed here such that as theta runs from 0 to pi,
!     one first encounters the ep/mu = bmax point and then
!     the ep/mu = bmin point.
!
      ku_lw = 1
      ku_up = ntheta
      jstar = zero
      dzeta = one
      IF (nzeta .gt. 1) dzeta = one/REAL((nzeta-1),rprec)
      test_bmax = ep_mu - bmax(:ntheta)
      test_bmin = ep_mu - bmin(:ntheta)

      test1(2:ntheta) = test_bmax(2:ntheta)*test_bmax(:ntheta-1)
      test2(2:ntheta) = test_bmin(2:ntheta)*test_bmin(:ntheta-1)
      DO ku = 2,ntheta
          IF(test1(ku) .le. zero) ku_lw = ku
          IF(test2(ku) .le. zero) ku_up = ku
      END DO

      IF (ku_lw .ge. ku_up) RETURN

      ALLOCATE (vpar1(nzeta,ku_up-ku_lw+1), l3v(nzeta,ku_up-ku_lw+1),
     1   stat=istat)
      IF (istat .ne. 0) THEN
         IF (myid .eq. master) PRINT *,' Allocation error in J_STAR!'    !MPI
         RETURN
      END IF

      vpar1 = one - modb(:,ku_lw:ku_up)/ep_mu
      l3v = vpar1 .gt. zero
      WHERE (l3v) vpar1 = SQRT(vpar1)/modb(:,ku_lw:ku_up)
      jstar(ku_lw:ku_up) = dzeta * SUM(vpar1, mask=l3v, dim=1)

      DEALLOCATE (vpar1, l3v)

      END SUBROUTINE j_star
