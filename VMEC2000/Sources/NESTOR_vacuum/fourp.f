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
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: g1, g2
      REAL(rprec) :: cosm, sinm, cosn, sinn, kernel, gcos, gsin
C-----------------------------------------------
!
!     PERFORM KV (TOROIDAL ANGLE) TRANSFORM (OVER UNPRIMED MESH IN EQ. 2.14)
!     THUS, THE (m,n) INDEX HERE CORRESPONDS TO THE FIRST INDEX OF AMATRIX
!     NOTE: THE .5 FACTOR (IN COSN,SINN) ACCOUNTS FOR THE SUM IN KERNELM
!     ON ENTRY THE FIRST TIME, GRPMN IS SIN,COS * Kmn(analytic)
!
!     THE 3rd INDEX OF GRPMN IS THE PRIMED U,V MESH COORDINATE     
!
      IF (ndim .GT. 2) STOP 'NDIM > 2'

      ALLOCATE (g1(istore,0:nf,ndim), g2(istore,0:nf,ndim),
     1          stat=m)
      IF (m .NE. 0) STOP 'Allocation error in fourp'

      DO ku = 1,nu2
         g1 = 0
         g2 = 0
         DO kv = 1,nv
            DO n = 0,nf
               cosn = p5*onp*cosv(n,kv)
               sinn = p5*onp*sinv(n,kv)
               iuv = kv+nv*(ku-1)
               ireflect = imirr(iuv)
               DO isym = 1, ndim
                  DO ip = 1,istore
                     IF (isym .eq. 1) THEN
                        kernel =
     1                  grp(iuv,ip) - grp(ireflect,ip)           !anti-symmetric part (u,v -> -u,-v)
                     ELSE IF (isym .eq. 2) THEN
                        kernel = 
     1                  grp(iuv,ip) + grp(ireflect,ip)           !symmetric part
                     END IF
                     g1(ip,n,isym)=g1(ip,n,isym) + cosn*kernel
                     g2(ip,n,isym)=g2(ip,n,isym) + sinn*kernel
                  END DO
               END DO
            END DO
         END DO
!
!     PERFORM KU (POLOIDAL ANGLE) TRANFORM [COMPLETE SIN(mu-nv) / COS(mu-nv) TRANSFORM]
!
         DO m = 0,mf
            DO isym = 1, ndim
               IF (isym .EQ. 1) THEN
                  cosm = -cosui(m,ku)
                  sinm =  sinui(m,ku)
               ELSE IF (isym .EQ. 2) THEN
                  sinm = cosui(m,ku)
                  cosm = sinui(m,ku)
               END IF
               DO n= 0,nf
                 DO ip = 1,istore
                     gcos = g1(ip,n,isym)*sinm
                     gsin = g2(ip,n,isym)*cosm
                     grpmn(m,n,ip+istart,isym) = 
     1               grpmn(m,n,ip+istart,isym) + gcos + gsin
                     IF (n .NE. 0) THEN
                        grpmn(m,-n,ip+istart,isym) = 
     1                  grpmn(m,-n,ip+istart,isym) + gcos - gsin
                     ENDIF
                 END DO
               END DO
            END DO
         END DO
      END DO

      istart = iend

      DEALLOCATE (g1, g2, stat=m)

      END SUBROUTINE fourp
