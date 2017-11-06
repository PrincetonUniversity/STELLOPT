      MODULE INIT_GEOMETRY

      LOGICAL     :: lflip

      CONTAINS

      SUBROUTINE flip_theta(rmn, zmn, lmn)
      USE vmec_main
      USE vmec_params, ONLY: ntmax, rcc, rss, zsc, zcs,                 &
                                    zcc, zss, rsc, rcs
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(0:ntor,0:mpol1,ntmax),                     &
         INTENT(inout) :: rmn, zmn
      REAL(rprec), DIMENSION(0:ntor,0:mpol1,ntmax),                     &
        INTENT(inout), OPTIONAL :: lmn
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: n, m
      REAL(rprec) :: mul1
      LOGICAL :: l_lmn
!-----------------------------------------------
!
!     FLIP THETA -> PI - THETA (INITIALLY, TO MAKE JACOBIAN < 0)
!
      mul1=-1
      l_lmn = PRESENT(lmn)
      DO m=1,mpol1
         DO n=0,ntor
            rmn(n,m,rcc) = mul1*rmn(n,m,rcc)
            zmn(n,m,zsc) =-mul1*zmn(n,m,zsc)
            IF (l_lmn) lmn(n,m,zsc) =-mul1*lmn(n,m,zsc)
            IF (lthreed) THEN
               rmn(n,m,rss) =-mul1*rmn(n,m,rss)
               zmn(n,m,zcs) = mul1*zmn(n,m,zcs)
               IF (l_lmn) lmn(n,m,zcs) = mul1*lmn(n,m,zcs)
            END IF
            IF (lasym) THEN
               rmn(n,m,rsc) =-mul1*rmn(n,m,rsc)
               zmn(n,m,zcc) = mul1*zmn(n,m,zcc)
               IF (l_lmn) lmn(n,m,zcc) = mul1*lmn(n,m,zcc)
               IF (lthreed) THEN
                  rmn(n,m,rcs) = mul1*rmn(n,m,rcs)
                  zmn(n,m,zss) =-mul1*zmn(n,m,zss)
                  IF (l_lmn) lmn(n,m,zss) =-mul1*lmn(n,m,zss)
               END IF
            END IF
         END DO

         mul1 = -mul1
 
      END DO

      END SUBROUTINE flip_theta

      END MODULE INIT_GEOMETRY
