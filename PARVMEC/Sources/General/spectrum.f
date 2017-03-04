      SUBROUTINE spectrum_par(rmn, zmn)
      USE parallel_include_module
      USE vmec_main
      USE vmec_params, ONLY: mscale, nscale, ntmax, rss, zcs, rsc, zcc
      USE totzsp_mod, ONLY:  convert_sym, convert_asym
#if defined(SKS)
      USE totzsp_mod, ONLY:  convert_sym_par, convert_asym_par
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(0:ntor,0:mpol1,ns,ntmax), 
     1             INTENT(inout) :: rmn, zmn
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: js, ntype, n, m, nsmin, nsmax
      REAL(rprec), DIMENSION(ns) :: t1, dnumer, denom
      REAL(rprec) :: scale
C-----------------------------------------------
!
!     CONVERT FROM INTERNAL XC REPRESENTATION FOR m=1 MODES, R+(at rsc) = .5(rsc + zcc),
!     R-(at zcc) = .5(rsc - zcc), TO REQUIRED rsc, zcc FORMS
!
      nsmin=MAX(2,tlglob); nsmax=MIN(t1rglob,ns)
#ifndef _HBANGLE
      IF (lthreed) CALL convert_sym_par (rmn(:,:,:,rss),
     1           zmn(:,:,:,zcs),nsmin, nsmax)
      IF (lasym)   CALL convert_asym_par (rmn(:,:,:,rsc), 
     1           zmn(:,:,:,zcc), nsmin, nsmax)
#endif
      dnumer(nsmin:nsmax) = zero
      denom(nsmin:nsmax) = zero
      DO ntype = 1,ntmax
        DO n = 0,ntor
          DO m = 1,mpol1
             scale = (mscale(m)*nscale(n))**2
             DO js = nsmin,nsmax
                t1(js) =(rmn(n,m,js,ntype)**2 + zmn(n,m,js,ntype)**2)
     2               *scale
             END DO
             dnumer(nsmin:nsmax) = dnumer(nsmin:nsmax) + 
     1                             t1(nsmin:nsmax)*xmpq(m,3)
             denom (nsmin:nsmax) = denom (nsmin:nsmax) +
     1                             t1(nsmin:nsmax)*xmpq(m,2)        
          END DO
        END DO
      ENDDO

      specw(nsmin:nsmax) = dnumer(nsmin:nsmax)/denom(nsmin:nsmax)
#endif
      END SUBROUTINE spectrum_par

      SUBROUTINE spectrum(rmn, zmn)
      USE vmec_main
      USE vmec_params, ONLY: mscale, nscale, ntmax, rss, zcs, rsc, zcc
      USE totzsp_mod, ONLY:  convert_sym, convert_asym
      USE parallel_include_module
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,ntmax), 
     1             INTENT(inout) :: rmn, zmn
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: js, ntype, n, m, nsmin, nsmax
      REAL(rprec), DIMENSION(ns) :: t1, dnumer, denom
      REAL(rprec) :: scale
C-----------------------------------------------
!
!     CONVERT FROM INTERNAL XC REPRESENTATION FOR m=1 MODES, R+(at rsc) = .5(rsc + zcc),
!     R-(at zcc) = .5(rsc - zcc), TO REQUIRED rsc, zcc FORMS
!
      nsmin=MAX(2,tlglob); nsmax=MIN(t1rglob,ns)
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

