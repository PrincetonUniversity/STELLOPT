      SUBROUTINE analysum2(grpmn, bvec, slp, tlp, slm, tlm,
     1    m, n, l, ivacskip, ndim)
      USE vacmod
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: m, n, l, ivacskip, ndim
      REAL(rprec), INTENT(inout) :: grpmn(0:mf,-nf:nf,nuv2,ndim)
      REAL(rprec), INTENT(inout) :: bvec(0:mf,-nf:nf,ndim)
      REAL(rprec), DIMENSION(nuv2), INTENT(in) ::
     1   slp, tlp, slm, tlm
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: istat
      REAL(rprec), ALLOCATABLE :: temp(:)
      REAL(rprec), POINTER     :: sinp(:), sinm(:)
      REAL(rprec), POINTER     :: cosp(:), cosm(:)
C-----------------------------------------------
      ALLOCATE (sinp(nuv2), sinm(nuv2), temp(nuv2), stat=istat)
      IF (istat .ne. 0) STOP 'Allocation error in analysum2'
      
      IF (n .lt. 0) STOP 'error calling analysum2!'
      
      sinp =  sinu1(:,m)*cosv1(:,n)*cmns(l,m,n)
      temp = -cosu1(:,m)*sinv1(:,n)*cmns(l,m,n)
      sinm = sinp - temp                 !SIN(mu + |n|v) * cmns (l,m,|n|)
      sinp = sinp + temp                 !SIN(mu - |n|v) * cmns (l,m,|n|)
      bvec(m,n,1)  = bvec(m,n,1)  + SUM(tlp*bexni*sinp)
      bvec(m,-n,1) = bvec(m,-n,1) + SUM(tlm*bexni*sinm)

      IF (ivacskip .eq. 0) THEN
         grpmn(m,n,:,1)  = grpmn(m,n,:,1)  + slp*sinp
         grpmn(m,-n,:,1) = grpmn(m,-n,:,1) + slm*sinm
      END IF


      IF (lasym) THEN
         cosp => sinp;  cosm => sinm
         cosp = cosu1(:,m)*cosv1(:,n)*cmns(l,m,n)
         temp = sinu1(:,m)*sinv1(:,n)*cmns(l,m,n)
         cosm = cosp - temp                !COS(mu + |n|v) * cmns (l,m,|n|)
         cosp = cosp + temp                !COS(mu - |n|v) * cmns (l,m,|n|)
         bvec(m,n,2)  = bvec(m,n,2)  + SUM(tlp*bexni*cosp)
         bvec(m,-n,2) = bvec(m,-n,2) + SUM(tlm*bexni*cosm)

         IF (ivacskip .eq. 0) THEN
            grpmn(m,n,:,2)  = grpmn(m,n,:,2)  + slp*cosp
            grpmn(m,-n,:,2) = grpmn(m,-n,:,2) + slm*cosm
         END IF

      END IF

      DEALLOCATE (sinp, sinm, temp)

      END SUBROUTINE analysum2
