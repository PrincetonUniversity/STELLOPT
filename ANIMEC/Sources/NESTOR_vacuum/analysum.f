      SUBROUTINE analysum(grpmn, bvec, sl, tl, m, n, l, ivacskip, ndim)
      USE vacmod
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: m, n, l, ivacskip, ndim
      REAL(rprec), INTENT(inout) :: grpmn(0:mf,-nf:nf,nuv2,ndim)
      REAL(rprec), INTENT(inout) :: bvec(0:mf,-nf:nf,ndim)
      REAL(rprec), DIMENSION(nuv2), INTENT(in) :: sl, tl
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: istat
      REAL(rprec), ALLOCATABLE, TARGET :: trigp(:)
      REAL(rprec), POINTER :: sinp(:), cosp(:)
C-----------------------------------------------
      ALLOCATE (trigp(nuv2), stat=istat)
      IF (istat .ne. 0) STOP 'Allocation error in analysum'

      IF (n .lt. 0) STOP 'error calling analysum!'

      sinp => trigp
      sinp = (sinu1(:,m)*cosv1(:,n) - sinv1(:,n)*cosu1(:,m))
     1     *  cmns(l,m,n)                            !SIN(mu - |n|v)*cmns
      bvec(m,n,1) = bvec(m,n,1) + SUM(tl*bexni*sinp)

      IF (ivacskip .eq. 0) grpmn(m,n,:,1) = grpmn(m,n,:,1) + sl*sinp

      IF (lasym) THEN

         cosp => trigp
         cosp = (cosu1(:,m)*cosv1(:,n) + sinv1(:,n)*sinu1(:,m))
     1        *  cmns(l,m,n)                            !COS(mu - |n|v)*cmns
         bvec(m,n,2) = bvec(m,n,2) + SUM(tl*bexni*cosp)

         IF (ivacskip .eq. 0) grpmn(m,n,:,2) = grpmn(m,n,:,2) + sl*cosp

      END IF

      DEALLOCATE (trigp)

      END SUBROUTINE analysum
