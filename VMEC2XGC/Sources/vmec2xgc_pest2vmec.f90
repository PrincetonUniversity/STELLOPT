!-----------------------------------------------------------------------
!     Subroutine:    vmec2xgc_pest2vmec
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          01/17/2017
!     Description:   This subroutine returns the VMEC poloidal angle
!                    given a PEST poloidal angle.  We use Fourier
!                    transforms as spline quantities (PSPLINE) was
!                    found to be insufficient.
!-----------------------------------------------------------------------
      SUBROUTINE vmec2xgc_pest2vmec(s,ustar,v,u)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE read_wout_mod
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(in) :: s, ustar, v
      REAL(rprec), INTENT(out) :: u
      INTEGER     :: mn,k1,k2,n1
      REAL(rprec) :: frac, ds, dth, lam, dlam, sinu, cosu
      REAL(rprec), ALLOCATABLE :: cnp(:),snp(:),fmns(:)
      DOUBLE PRECISION, PARAMETER      :: one = 1.0D+00
      DOUBLE PRECISION, PARAMETER      :: search_tol = 1.0D-12
      DOUBLE PRECISION, PARAMETER      :: pi2 = 6.283185482025146D+00
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ALLOCATE(cnp(mnmax),snp(mnmax),fmns(mnmax))
      FORALL(mn = 1:mnmax) cnp(mn) = cos(xn(mn)*v/nfp)
      FORALL(mn = 1:mnmax) snp(mn) = sin(xn(mn)*v/nfp)
      ! Determine s and interpolate (linear)
      ds = one/REAL(ns)
      k1 = FLOOR(s*ns)
      k2 = CEILING(s*ns)
      frac = (s - REAL(k1)/REAL(ns))*ns
      FORALL(mn = 1:mnmax) fmns(mn) = (1-frac)*lmns(mn,k1)+frac*lmns(mn,k2)
      dth = one
      n1 = 0
      u = ustar
      DO WHILE(ABS(dth) >= search_tol .and. n1 < 5000)
         lam = 0
         dlam = 0
         DO mn = 1, mnmax
            cosu = fmns(mn)*cos(xm(mn)*u)
            sinu = fmns(mn)*sin(xm(mn)*u)
            lam = lam + sinu*cnp(mn)-cosu*snp(mn)
            dlam = dlam + xm(mn)*(cosu*cnp(mn)+sinu*snp(mn))
         END DO
         dth = -(u + lam - ustar)/(one+dlam)
         n1 = n1 + 1
         u = u + 0.5*dth
      END DO
      DEALLOCATE(cnp,snp,fmns)
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE vmec2xgc_pest2vmec
