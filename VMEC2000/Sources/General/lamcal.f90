      SUBROUTINE lamcal_par(overg, guu, guv, gvv)
      USE vmec_main
      USE vmec_params, ONLY: ntmax, jlam, lamscale
      USE realspace, ONLY: psqrts
      USE parallel_include_module
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(nznt,ns), INTENT(in) ::                    &
         overg, guu, guv, gvv
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: damping_fac = 2
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: m, n, js, nsmin, nsmax, numjs, i, j, k, l 
      REAL(rprec) :: tnn, tnm, tmm, power, pfactor0, pfactor

      INTEGER :: blksize
      INTEGER, ALLOCATABLE, DIMENSION(:) :: counts, disps
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:,:) :: send_buf
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: recv_buf
      REAL(rprec) :: allgvton, allgvtoff
!-----------------------------------------------

      nsmin=tlglob; nsmax=t1rglob

      blam(nsmin:nsmax) = SUM(guu(:,nsmin:nsmax)*overg(:,nsmin:nsmax), dim=1)
      clam(nsmin:nsmax) = SUM(gvv(:,nsmin:nsmax)*overg(:,nsmin:nsmax), dim=1)
      dlam(nsmin:nsmax) = SUM(guv(:,nsmin:nsmax)*overg(:,nsmin:nsmax), dim=1)
      blam(1) = blam(2)
      clam(1) = clam(2)
      dlam(1) = dlam(2)
      blam(ns+1) =  0
      clam(ns+1) =  0
      dlam(ns+1) =  0

      nsmin=MAX(2,tlglob); nsmax=trglob
      DO js = nsmin, nsmax
         blam(js) = cp5*(blam(js) + blam(js+1))
         clam(js) = cp5*(clam(js) + clam(js+1))
         dlam(js) = cp5*(dlam(js) + dlam(js+1))
      END DO

      pfaclam = 0
!      pfactor0 = 2*damping_fac/(2*r0scale*lamscale)**2
      pfactor0 = damping_fac/(2*r0scale*lamscale)**2

      DO m = 0, mpol1
         tmm = m*m
         power = MIN(tmm/256, 8._dp)
         pfactor = pfactor0
         DO n = 0, ntor
            IF (m.eq.0 .and. n.eq.0) CYCLE
            tnn = (n*nfp)**2
            tnm = 2*m*n*nfp

            DO js = MAX(jlam(m),tlglob), trglob
               pfaclam(n,m,js,1) = blam(js)*tnn + clam(js)*tmm                 &
                                 + SIGN(dlam(js),blam(js))*tnm
               IF (pfaclam(n,m,js,1) .eq. zero) THEN
                  pfaclam(n,m,js,1) = -1.E-10_dp
               END IF
               pfaclam(n,m,js,1) = (pfactor/pfaclam(n,m,js,1))                 &
                                 * psqrts(1,js)**power !Damps m > 16 modes
            END DO
         END DO
      END DO

      nsmin=tlglob; nsmax=trglob
      DO n = 2, ntmax
         pfaclam(0:ntor,0:mpol1,nsmin:nsmax,n) =                               &
            pfaclam(0:ntor,0:mpol1,nsmin:nsmax,1)
      END DO

!
!     ADD NORM FOR CHIP (PREVIOUSLY IOTA) FORCE, STORED IN lmnsc(m=0,n=0) COMPONENT
!
      nsmin=tlglob; nsmax=trglob
      DO js = nsmin, nsmax
         pfaclam(0,0,js,1) = (pfactor0*lamscale**2)/blam(js)
      END DO

      END SUBROUTINE lamcal_par

      SUBROUTINE lamcal(overg, guu, guv, gvv)
      USE vmec_main
      USE vmec_params, ONLY: ntmax, jlam, lamscale
      USE realspace, ONLY: sqrts

      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(ns,nznt), INTENT(in) ::                    &
         overg, guu, guv, gvv
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: damping_fac=2
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: m,n,js
      REAL(rprec) :: tnn, tnm, tmm, power, pfactor0, pfactor

!-----------------------------------------------


      blam(:ns) = SUM(guu*overg, dim=2)
      clam(:ns) = SUM(gvv*overg, dim=2)
      dlam(:ns) = SUM(guv*overg, dim=2)
      blam(1) = blam(2)
      clam(1) = clam(2)
      dlam(1) = dlam(2)
      blam(ns+1) =  0
      clam(ns+1) =  0
      dlam(ns+1) =  0
      DO js = 2, ns
         blam(js) = cp5*(blam(js) + blam(js+1))
         clam(js) = cp5*(clam(js) + clam(js+1))
         dlam(js) = cp5*(dlam(js) + dlam(js+1))
      END DO

      faclam = 0
      pfactor0 = damping_fac/(2*r0scale*lamscale)**2

      DO m = 0, mpol1
         tmm = m*m
         power = MIN(tmm/256, 8._dp)
         pfactor = pfactor0
         DO n = 0, ntor
            IF (m.eq.0 .and. n.eq.0) CYCLE
!            IF (n .gt. 1) pfactor = pfactor0/4                          !sometimes helps convergence
            tnn = (n*nfp)**2
            tnm = 2*m*n*nfp
            DO js = jlam(m), ns
               faclam(js,n,m,1) = blam(js)*tnn + clam(js)*tmm                  &
                                + SIGN(dlam(js),blam(js))*tnm
               IF (faclam(js,n,m,1) .eq. zero) THEN
                  faclam(js,n,m,1) = -1.E-10_dp
               END IF
               faclam(js,n,m,1) = (pfactor/faclam(js,n,m,1))                   &
                               * sqrts(js)**power !Damps m > 16 modes
            END DO
         END DO
      END DO

      DO n = 2, ntmax
         faclam(:ns,0:ntor,0:mpol1,n) = faclam(:ns,0:ntor,0:mpol1,1)
      END DO
!
!     ADD NORM FOR CHIP (PREVIOUSLY IOTA) FORCE, STORED IN lmnsc(m=0,n=0) COMPONENT
!
      DO js = 1, ns
         faclam(js,0,0,1) = (pfactor0*lamscale**2)/blam(js)
      END DO

      END SUBROUTINE lamcal
