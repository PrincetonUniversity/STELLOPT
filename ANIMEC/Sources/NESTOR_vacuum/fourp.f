      SUBROUTINE fourp (grpmn, grp, istore, istart, iend, ndim)
      USE vacmod
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(inout) :: istart
      INTEGER, INTENT(in)    :: iend, istore, ndim
      REAL(rprec), INTENT(in)    :: grp(nuv,istore)
      REAL(rprec), INTENT(inout) :: grpmn(0:mf,-nf:nf,nuv2,ndim)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: n, kv, ku, ip, iuv, m, ireflect, isym
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:,:) :: g1, g2
      REAL(rprec), ALLOCATABLE :: kernel(:), gcos(:),gsin(:)
      REAL(rprec) :: cosm, sinm, cosn, sinn
C-----------------------------------------------
!
!     PERFORM KV (TOROIDAL ANGLE) TRANSFORM (OVER UNPRIMED MESH IN EQ. 2.14)
!     THUS, THE (m,n) INDEX HERE CORRESPONDS TO THE FIRST INDEX OF AMATRIX
!     NOTE: THE .5 FACTOR (IN COSN,SINN) ACCOUNTS FOR THE SUM IN KERNELM
!     ON ENTRY THE FIRST TIME, GRPMN IS SIN,COS * Kmn(analytic)
!
!     THE 3rd INDEX OF GRPMN IS THE PRIMED U,V MESH COORDINATE     
!
      ALLOCATE (g1(istore,nu2,0:nf,ndim), g2(istore,nu2,0:nf,ndim),
     1   kernel(istore), gcos(istore), gsin(istore), stat=m)
      IF (m .ne. 0) STOP 'Allocation error in fourp'

      g1 = 0
      g2 = 0

      DO n = 0,nf
         DO kv = 1,nv
            cosn = p5*onp*cosv(n,kv)
            sinn = p5*onp*sinv(n,kv)
            iuv = kv
            DO ku = 1,nu2
               ireflect = imirr(iuv)
               DO isym = 1, ndim
                  IF (isym .eq. 1) THEN
                     kernel(:istore) =
     1                  grp(iuv,:istore) - grp(ireflect,:istore)           !anti-symmetric part (u,v -> -u,-v)
                  ELSE IF (isym .eq. 2) THEN
                     kernel(:istore) = 
     1                  grp(iuv,:istore) + grp(ireflect,:istore)           !symmetric part
                  END IF
                  DO ip = 1,istore
                     g1(ip,ku,n,isym)=g1(ip,ku,n,isym) + cosn*kernel(ip)
                     g2(ip,ku,n,isym)=g2(ip,ku,n,isym) + sinn*kernel(ip)
                  END DO
               END DO
               iuv = iuv + nv
            END DO
         END DO
      END DO

!
!     PERFORM KU (POLOIDAL ANGLE) TRANFORM [COMPLETE SIN(mu-nv) / COS(mu-nv) TRANSFORM]
!
      
      DO m = 0,mf
         DO ku = 1,nu2
            DO isym = 1, ndim
               IF (isym .eq. 1) THEN
                  cosm = -cosui(m,ku)
                  sinm =  sinui(m,ku)
               ELSE IF (isym .eq. 2) THEN
                  sinm = cosui(m,ku)
                  cosm = sinui(m,ku)
               END IF
                  DO n= 0,nf
                     DO ip = 1,istore
                        gcos(ip) = g1(ip,ku,n,isym)*sinm
                        gsin(ip) = g2(ip,ku,n,isym)*cosm
                        grpmn(m,n,ip+istart,isym) = 
     1                  grpmn(m,n,ip+istart,isym) + gcos(ip) + gsin(ip)
                     END DO

                     IF (n .ne. 0) THEN
                        DO ip = 1,istore
                           grpmn(m,-n,ip+istart,isym) = 
     1                     grpmn(m,-n,ip+istart,isym)
     2                         + gcos(ip) - gsin(ip)
                        END DO
                     ENDIF
                  END DO
              END DO
         END DO
      END DO

      istart = iend

      DEALLOCATE (g1, g2, kernel, gcos, gsin, stat=m)

      END SUBROUTINE fourp
