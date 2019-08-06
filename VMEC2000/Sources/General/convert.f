      SUBROUTINE convert_par(rmnc,zmns,lmns,rmns,zmnc,lmnc,rzl_array)
      USE vmec_main
      USE vmec_params
      USE parallel_include_module
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(dp), DIMENSION(mnmax), INTENT(OUT) ::
     &    rmnc, zmns, lmns, rmns, zmnc, lmnc
      REAL(dp), DIMENSION(0:ntor,0:mpol1,ns,3*ntmax),
     &    INTENT(INOUT) :: rzl_array
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(dp), PARAMETER :: p5 = 0.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: rmncc, rmnss, rmncs, rmnsc, zmncs, zmnsc, 
     &           zmncc, zmnss, lmncs, lmnsc, lmncc, lmnss
      INTEGER :: mn, m, n, n1, bufsize, js
      REAL(dp) :: t1, sign0, mul1, tbroadon, tbroadoff
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: bcastbuf
C-----------------------------------------------
!
!     FOR EDGE (js=ns) ONLY:
!     CONVERTS INTERNAL MODE REPRESENTATION TO STANDARD
!     FORM FOR OUTPUT (COEFFICIENTS OF COS(mu-nv), SIN(mu-nv) WITHOUT mscale,nscale norms)
!
      js = ns
#if defined(MPI_OPT)
      bufsize = (ntor+1)*(mpol1+1)*3*ntmax
      ALLOCATE(bcastbuf(bufsize))
      mn=0
      DO n1 = 1, 3*ntmax
         DO m = 0, mpol1
            DO n = 0, ntor
               mn = mn + 1
               bcastbuf(mn) = rzl_array(n,m,js,n1)
            END DO
         END DO
      END DO
      CALL second0(tbroadon)
      CALL MPI_Bcast(bcastbuf, bufsize, MPI_REAL8, nranks - 1,
     &               NS_COMM, MPI_ERR)
      IF(vlactive) THEN
        CALL MPI_Bcast(bcastbuf, bufsize, MPI_REAL8, 0,
     &                 VAC_COMM, MPI_ERR)
      END IF
      CALL second0(tbroadoff)
      broadcast_time = broadcast_time + (tbroadoff -tbroadon)

      mn=0
      DO n1 = 1, 3*ntmax
         DO m = 0, mpol1
            DO n = 0, ntor
               mn = mn + 1
               rzl_array(n,m,js,n1) = bcastbuf(mn)
            END DO
         END DO
      END DO
      DEALLOCATE(bcastbuf)
#endif

      rmncc = rcc
      rmnss = rss
      rmnsc = rsc
      rmncs = rcs
      zmnsc = zsc + ntmax
      zmncc = zcc + ntmax
      zmncs = zcs + ntmax
      zmnss = zss + ntmax
      lmnsc = zsc + 2*ntmax
      lmncc = zcc + 2*ntmax
      lmncs = zcs + 2*ntmax
      lmnss = zss + 2*ntmax

!
!     DO M = 0 MODES SEPARATELY (ONLY KEEP N >= 0 HERE: COS(-NV), SIN(-NV))
!
      mn = 0; m = 0
      zmns(1:ntor+1) = 0;  lmns(1:ntor+1) = 0
      DO n = 0, ntor
         t1 = mscale(m)*nscale(n)
         mn = mn + 1
         rmnc(mn) = t1*rzl_array(n,m,js,rmncc)
         IF (.not. lthreed) CYCLE
         zmns(mn) =-t1*rzl_array(n,m,js,zmncs)
         lmns(mn) =-t1*rzl_array(n,m,js,lmncs)
      END DO

      lmns(1) = 0                         !may have been used for storing iota variation...

      DO m = 1, mpol1
         DO n = -ntor, ntor
            n1 = ABS(n)
            t1 = mscale(m)*nscale(n1)
            mn = mn + 1
            IF (n .eq. 0) THEN
               rmnc(mn) = t1*rzl_array(n,m,js,rmncc)
               zmns(mn) = t1*rzl_array(n,m,js,zmnsc)
               lmns(mn) = t1*rzl_array(n,m,js,lmnsc)
            ELSE IF (js .gt. 1) THEN
               sign0 = n/n1
               IF (.not.lthreed) sign0 = 0
               rmnc(mn) = p5*t1*(rzl_array(n1,m,js,rmncc) +
     &                           sign0*rzl_array(n1,m,js,rmnss))
               zmns(mn) = p5*t1*(rzl_array(n1,m,js,zmnsc) -
     &                           sign0*rzl_array(n1,m,js,zmncs))
               lmns(mn) = p5*t1*(rzl_array(n1,m,js,lmnsc) -
     &                           sign0*rzl_array(n1,m,js,lmncs))
            ELSE IF (js .eq. 1) THEN
               rmnc(mn) = 0
               zmns(mn) = 0
               lmns(mn) = 0
            END IF
         END DO
      END DO

      IF (mn .ne. mnmax) STOP 'Error in Convert!'

      IF (.not.lasym) THEN
         rmns = 0
         zmnc = 0
         lmnc = 0
         RETURN
      END IF

      mn = 0; m = 0
      rmns(1:ntor+1) = 0
      DO n = 0, ntor
         t1 = mscale(m)*nscale(n)
         mn = mn + 1
         zmnc(mn) = t1*rzl_array(n,m,js,zmncc)
         lmnc(mn) = t1*rzl_array(n,m,js,lmncc)
         IF (.not.lthreed) CYCLE
         rmns(mn) =-t1*rzl_array(n,m,js,rmncs)                           !ers-fixed sign
      END DO

      mul1 = 1
      IF (.not.lthreed) mul1 = 0
      DO m = 1, mpol1
         DO n = -ntor, ntor
            n1 = ABS(n)
            t1 = mscale(m)*nscale(n1)
            mn = mn + 1
            IF (n .eq. 0) THEN
               rmns(mn) = t1*rzl_array(n,m,js,rmnsc)
               zmnc(mn) = t1*rzl_array(n,m,js,zmncc)
               lmnc(mn) = t1*rzl_array(n,m,js,lmncc)
            ELSE IF (js .gt. 1) THEN
               sign0 = n/n1
               rmns(mn) = p5*t1*(mul1*rzl_array(n1,m,js,rmnsc) -
     &                           sign0*rzl_array(n1,m,js,rmncs))
               zmnc(mn) = p5*t1*(mul1*rzl_array(n1,m,js,zmncc) +
     &                           sign0*rzl_array(n1,m,js,zmnss))
               lmnc(mn) = p5*t1*(mul1*rzl_array(n1,m,js,lmncc) +
     &                           sign0*rzl_array(n1,m,js,lmnss))
            ELSE IF (js .eq. 1) THEN
               rmns(mn) = 0
               zmnc(mn) = 0
               lmnc(mn) = 0
            END IF
         END DO
      END DO

      END SUBROUTINE convert_par
      
      SUBROUTINE convert(rmnc,zmns,lmns,rmns,zmnc,lmnc,rzl_array,js)
      USE vmec_main
      USE vmec_params
      USE parallel_include_module
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(IN) :: js
      REAL(dp), DIMENSION(mnmax), INTENT(out) ::
     &    rmnc, zmns, lmns, rmns, zmnc, lmnc
      REAL(dp), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax),
     &    INTENT(in) :: rzl_array
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(dp), PARAMETER :: p5 = 0.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: rmncc, rmnss, rmncs, rmnsc, zmncs, zmnsc, 
     &           zmncc, zmnss, lmncs, lmnsc, lmncc, lmnss
      INTEGER :: mn, m, n, n1
      REAL(dp) :: t1, sign0, mul1
C-----------------------------------------------
!
!     CONVERTS INTERNAL MODE REPRESENTATION TO STANDARD
!     FORM FOR OUTPUT (COEFFICIENTS OF COS(mu-nv), SIN(mu-nv) WITHOUT mscale,nscale norms)
!
      rmncc = rcc
      rmnss = rss
      rmnsc = rsc
      rmncs = rcs
      zmnsc = zsc + ntmax
      zmncc = zcc + ntmax
      zmncs = zcs + ntmax
      zmnss = zss + ntmax
      lmnsc = zsc + 2*ntmax
      lmncc = zcc + 2*ntmax
      lmncs = zcs + 2*ntmax
      lmnss = zss + 2*ntmax

!
!     DO M = 0 MODES SEPARATELY (ONLY KEEP N >= 0 HERE: COS(-NV), SIN(-NV))
!
      mn = 0; m = 0
      zmns(1:ntor+1) = 0;  lmns(1:ntor+1) = 0
      DO n = 0, ntor
         t1 = mscale(m)*nscale(n)
         mn = mn + 1
         rmnc(mn) = t1*rzl_array(js,n,m,rmncc)
         IF (.not. lthreed) CYCLE
         zmns(mn) =-t1*rzl_array(js,n,m,zmncs)
         lmns(mn) =-t1*rzl_array(js,n,m,lmncs)
      END DO

      IF (lthreed .and. js.eq.1) THEN
         mn = 0
         DO n = 0, ntor
            t1 = mscale(m)*nscale(n)
            mn = mn + 1
            lmns(mn) =-t1*(2*rzl_array(2,n,m,lmncs)
     &               -       rzl_array(3,n,m,lmncs))
         END DO
      END IF

      lmns(1) = 0                         !may have been used for storing iota variation...

      DO m = 1, mpol1
         DO n = -ntor, ntor
            n1 = ABS(n)
            t1 = mscale(m)*nscale(n1)
            mn = mn + 1
            IF (n .eq. 0) THEN
               rmnc(mn) = t1*rzl_array(js,n,m,rmncc)
               zmns(mn) = t1*rzl_array(js,n,m,zmnsc)
               lmns(mn) = t1*rzl_array(js,n,m,lmnsc)
            ELSE IF (js .gt. 1) THEN
               sign0 = n/n1
               IF (.not.lthreed) sign0 = 0
               rmnc(mn) = p5*t1*(rzl_array(js,n1,m,rmncc) +
     &                           sign0*rzl_array(js,n1,m,rmnss))
               zmns(mn) = p5*t1*(rzl_array(js,n1,m,zmnsc) -
     &                           sign0*rzl_array(js,n1,m,zmncs))
               lmns(mn) = p5*t1*(rzl_array(js,n1,m,lmnsc) -
     &                           sign0*rzl_array(js,n1,m,lmncs))
            ELSE IF (js .eq. 1) THEN
               rmnc(mn) = 0
               zmns(mn) = 0
               lmns(mn) = 0
            END IF
         END DO
      END DO

      IF (mn .ne. mnmax) STOP 'Error in Convert!'

      IF (.not.lasym) THEN
         rmns = 0
         zmnc = 0
         lmnc = 0
         RETURN
      END IF

      mn = 0; m = 0
      rmns(1:ntor+1) = 0
      DO n = 0, ntor
         t1 = mscale(m)*nscale(n)
         mn = mn + 1
         zmnc(mn) = t1*rzl_array(js,n,m,zmncc)
         lmnc(mn) = t1*rzl_array(js,n,m,lmncc)
         IF (.not.lthreed) CYCLE
         rmns(mn) =-t1*rzl_array(js,n,m,rmncs)                           !ers-fixed sign
      END DO

      mul1 = 1
      IF (.not.lthreed) mul1 = 0
      DO m = 1, mpol1
         DO n = -ntor, ntor
            n1 = ABS(n)
            t1 = mscale(m)*nscale(n1)
            mn = mn + 1
            IF (n .eq. 0) THEN
               rmns(mn) = t1*rzl_array(js,n,m,rmnsc)
               zmnc(mn) = t1*rzl_array(js,n,m,zmncc)
               lmnc(mn) = t1*rzl_array(js,n,m,lmncc)
            ELSE IF (js .gt. 1) THEN
               sign0 = n/n1
               rmns(mn) = p5*t1*(mul1*rzl_array(js,n1,m,rmnsc) -
     &                           sign0*rzl_array(js,n1,m,rmncs))
               zmnc(mn) = p5*t1*(mul1*rzl_array(js,n1,m,zmncc) +
     &                           sign0*rzl_array(js,n1,m,zmnss))
               lmnc(mn) = p5*t1*(mul1*rzl_array(js,n1,m,lmncc) +
     &                           sign0*rzl_array(js,n1,m,lmnss))
            ELSE IF (js .eq. 1) THEN
               rmns(mn) = 0
               zmnc(mn) = 0
               lmnc(mn) = 0
            END IF
         END DO
      END DO

      END SUBROUTINE convert
