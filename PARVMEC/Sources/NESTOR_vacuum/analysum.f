#if defined(SKS)      
      SUBROUTINE analysum_par(grpmn, bvec, sl, tl, m, n, l, ivacskip, 
     1                       ndim)
      USE vacmod
      USE parallel_include_module
      USE timer_sub
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: m, n, l, ivacskip, ndim
      REAL(rprec), INTENT(inout) :: grpmn(0:mf,-nf:nf,ndim,nuv3)
      REAL(rprec), INTENT(inout) :: bvec(0:mf,-nf:nf,ndim)
      REAL(rprec), DIMENSION(nuv3), INTENT(in) :: sl, tl
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i
      REAL(rprec) :: sinp, cosp
      REAL(rprec) :: skston, skstoff
C-----------------------------------------------

      CALL second0(skston)

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
      CALL second0(skstoff)
      timer_vac(tasum) = timer_vac(tasum) + (skstoff-skston)

      END SUBROUTINE analysum_par
#endif

      SUBROUTINE analysum(grpmn, bvec, sl, tl, m, n, l, ivacskip, ndim)
      USE vacmod
      USE parallel_include_module
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: m, n, l, ivacskip, ndim
      REAL(rprec), INTENT(inout) :: grpmn(0:mf,-nf:nf,nuv3,ndim)
      REAL(rprec), INTENT(inout) :: bvec(0:mf,-nf:nf,ndim)
      REAL(rprec), DIMENSION(nuv3), INTENT(in) :: sl, tl
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: istat, i, j, k, ll
      REAL(rprec), ALLOCATABLE, TARGET :: trigp(:)
      REAL(rprec), POINTER :: sinp(:), cosp(:)
C-----------------------------------------------

      ALLOCATE (trigp(nuv3), stat=istat)
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

!      DO k=nuv3min, nuv3max
!        DO ll=1, ndim
!          DO j=-nf, nf
!            DO i=0, mf
!              WRITE(2000+vrank,"(F20.10)") grpmn(i,j,ll,k)
!            END DO
!          END DO
!        END DO
!      END DO

      END SUBROUTINE analysum
