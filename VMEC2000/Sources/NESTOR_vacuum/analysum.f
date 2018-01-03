      SUBROUTINE analysum (grpmn, bvec, sl, tl, m, n, l, ivacskip, 
     1                     ndim)
      USE vacmod
      USE parallel_include_module
      USE timer_sub
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(IN) :: m, n, l, ivacskip, ndim
      REAL(dp), INTENT(INOUT) :: grpmn(0:mf,-nf:nf,ndim,nuv3)
      REAL(dp), INTENT(INOUT) :: bvec(0:mf,-nf:nf,ndim)
      REAL(dp), DIMENSION(nuv3), INTENT(IN) :: sl, tl
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER  :: i
      REAL(dp) :: sinp, cosp, ton, toff
C-----------------------------------------------
      CALL second0(ton)

      DO i = nuv3min, nuv3max
         sinp = (sinu1(i,m)*cosv1(i,n) - sinv1(i,n)*cosu1(i,m))
     1        *  cmns(l,m,n)                                   !SIN(mu - |n|v)*cmns
         IF (ivacskip .EQ. 0) grpmn(m,n,1,i) = grpmn(m,n,1,i) 
     1                                       + sl(i)*sinp
         bvec(m,n,1) = bvec(m,n,1) + tl(i)*bexni(i)*sinp

         IF (lasym) THEN
            cosp = (cosu1(i,m)*cosv1(i,n) + sinv1(i,n)*sinu1(i,m))
     1           *  cmns(l,m,n)                                !COS(mu - |n|v)*cmns

            IF (ivacskip .EQ. 0) grpmn(m,n,2,i) = grpmn(m,n,2,i) 
     1                                          + sl(i)*cosp
            bvec(m,n,2) = bvec(m,n,2) + tl(i)*bexni(i)*cosp
         END IF
      END DO
      CALL second0(toff)
      timer_vac(tasum) = timer_vac(tasum) + (toff-ton)

      END SUBROUTINE analysum
