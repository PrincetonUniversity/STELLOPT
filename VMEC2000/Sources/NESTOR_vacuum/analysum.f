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
      INTEGER :: i
      REAL(rprec) :: sinp, cosp
C-----------------------------------------------
      IF (n .LT. 0) STOP 'error calling analysum!'

      DO i = 1,nuv2
         sinp = (sinu1(i,m)*cosv1(i,n) - sinv1(i,n)*cosu1(i,m))
     1            *  cmns(l,m,n)                            !SIN(mu - |n|v)*cmns
         IF (ivacskip .EQ. 0) grpmn(m,n,i,1) = grpmn(m,n,i,1) 
     1                                       + sl(i)*sinp
         bvec(m,n,1) = bvec(m,n,1) + tl(i)*bexni(i)*sinp

         IF (lasym) THEN
            cosp = (cosu1(i,m)*cosv1(i,n) + sinv1(i,n)*sinu1(i,m))
     1           *  cmns(l,m,n)                                !COS(mu - |n|v)*cmns
            IF (ivacskip .EQ. 0) grpmn(m,n,i,2) = grpmn(m,n,i,2) 
     1                                          + sl(i)*cosp
            bvec(m,n,2) = bvec(m,n,2) + tl(i)*bexni(i)*cosp
         END IF
      END DO

      END SUBROUTINE analysum
