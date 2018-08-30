      SUBROUTINE spectrum(rmn, zmn)
      USE vmec_main
      USE vmec_params, ONLY: mscale, nscale, ntmax, rss, zcs, rsc, zcc
      USE totzsp_mod, ONLY:  convert_sym, convert_asym
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,ntmax), 
     1             INTENT(inout) :: rmn, zmn
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: js, ntype, n, m
      REAL(rprec), DIMENSION(ns) :: t1, dnumer, denom
      REAL(rprec) :: scale
C-----------------------------------------------
!
!     CONVERT FROM INTERNAL XC REPRESENTATION FOR m=1 MODES, R+(at rsc) = .5(rsc + zcc),
!     R-(at zcc) = .5(rsc - zcc), TO REQUIRED rsc, zcc FORMS
!
#ifndef _HBANGLE
      IF (lthreed) CALL convert_sym (rmn(:,:,:,rss), zmn(:,:,:,zcs))
      IF (lasym)   CALL convert_asym (rmn(:,:,:,rsc), zmn(:,:,:,zcc))
#endif
      dnumer(2:ns) = zero
      denom(2:ns) = zero
      DO ntype = 1,ntmax
        DO n = 0,ntor
          DO m = 1,mpol1
             scale = (mscale(m)*nscale(n))**2
             DO js = 2,ns
                t1(js) =(rmn(js,n,m,ntype)**2 + zmn(js,n,m,ntype)**2)
     2               *scale
             END DO
             dnumer(2:ns) = dnumer(2:ns) + t1(2:ns)*xmpq(m,3)
             denom (2:ns) = denom (2:ns) + t1(2:ns)*xmpq(m,2)
          END DO
        END DO
      ENDDO

      specw(2:ns) = dnumer(2:ns)/denom(2:ns)

      END SUBROUTINE spectrum

