      SUBROUTINE convert(rmnc,zmns,lmns,rmns,zmnc,lmnc,rzl_array,js)
      USE vmec_main
      USE vmec_params
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER js
      REAL(rprec), DIMENSION(mnmax), INTENT(out) ::
     1    rmnc, zmns, lmns, rmns, zmnc, lmnc
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax),
     1    INTENT(in) :: rzl_array
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: p5 = 0.5_dp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: rmncc, rmnss, rmncs, rmnsc, zmncs, zmnsc, 
     1           zmncc, zmnss, lmncs, lmnsc, lmncc, lmnss
      INTEGER :: mn, m, n, n1
      REAL(rprec) :: t1, sign0 
!-----------------------------------------------
!
!     CONVERTS INTERNAL MODE REPRESENTATION TO STANDARD
!     FORM FOR OUTPUT (COEFFICIENTS OF COS(mu-nv), SIN(mu-nv) WITHOUT internal mscale,nscale norms)
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
     1               -       rzl_array(3,n,m,lmncs))
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
               rmnc(mn) = p5*t1*(rzl_array(js,n1,m,rmncc)+sign0*
     1            rzl_array(js,n1,m,rmnss))
               zmns(mn) = p5*t1*(rzl_array(js,n1,m,zmnsc)-sign0*
     1            rzl_array(js,n1,m,zmncs))
               lmns(mn) = p5*t1*(rzl_array(js,n1,m,lmnsc)-sign0*
     1            rzl_array(js,n1,m,lmncs))
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
               rmns(mn) = p5*t1*(rzl_array(js,n1,m,rmnsc)-sign0*    !ers-corrected rmnsc <-> rmncs 
     1            rzl_array(js,n1,m,rmncs))
               zmnc(mn) = p5*t1*(rzl_array(js,n1,m,zmncc)+sign0*
     1            rzl_array(js,n1,m,zmnss))
               lmnc(mn) = p5*t1*(rzl_array(js,n1,m,lmncc)+sign0*
     1            rzl_array(js,n1,m,lmnss))
            ELSE IF (js .eq. 1) THEN
               rmns(mn) = 0
               zmnc(mn) = 0
               lmnc(mn) = 0
            END IF
         END DO
      END DO

      END SUBROUTINE convert
