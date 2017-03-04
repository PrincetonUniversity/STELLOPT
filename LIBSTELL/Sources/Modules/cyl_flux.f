      MODULE cyl_flux
      USE stel_kinds
      USE stel_constants, ONLY: twopi, one, zero
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Module for coordinate conversion, cylindrical <- -> VMEC flux
!  CONTAINS two subroutines:
!     flx2cyl - convert from VMEC flux to cylindrical - inverse Fourier Transform
!     cyl2flx - convert from cylindrical to  VMEC flux - Root find using flx2cyl
!
!  cyl2flx CONTAINS two subroutines:
!     newt2d - implements a 2 dimensional Newton rootfind
!     get_flxcoord - interface with fewer arguments to flx2cyl
!  Note - since newt2d and get_flxcoord are CONTAINed in cyl2flx, they have 
!  access to all the variables in cyl2flx.
!
!  2011-09-06 JDH
!    Split this module from vmec_utils, to clarify interfaces and dependencies
!-------------------------------------------------------------------------------
      
      CONTAINS

!*******************************************************************************
!-------------------------------------------------------------------------------
!  Convert from flux to cylindrical cooridnates
!-------------------------------------------------------------------------------
      SUBROUTINE flx2cyl(rzl_array, c_flux, r_cyl, ns, ntor, 
     1                   mpol, ntmax, lthreed, lasym, iflag,
     2                   mscale, nscale, Ru, Rv, Zu, Zv)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(out) :: iflag
      INTEGER, INTENT(in)  :: ns, ntmax, mpol, ntor
      LOGICAL :: lthreed, lasym
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol-1,3*ntmax),
     1   INTENT(in) :: rzl_array
      REAL(rprec), INTENT(in) :: c_flux(3)
      REAL(rprec), INTENT(out) :: r_cyl(3)
      REAL(rprec), INTENT(in), OPTIONAL :: 
     1                            mscale(0:mpol-1), nscale(0:ntor)
      REAL(rprec), INTENT(out), OPTIONAL :: Ru, Rv, Zu, Zv
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: rcc = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: rss, rsc, rcs, zsc, zcs, zcc, zss
      INTEGER :: istat, jslo, jshi, mpol1, m, n
      REAL(rprec), DIMENSION(0:ntor,0:mpol-1) ::
     1           rmncc, rmnss, zmncs, zmnsc,
     2           rmncs, rmnsc, zmncc, zmnss
      REAL(rprec) :: wlo, whi, wlo_odd, whi_odd, hs1, 
     1               si, ui, vi, r11, z11
      REAL(rprec) :: wmins(0:ntor,0:mpol-1),
     1               wplus(0:ntor,0:mpol-1)
      REAL(rprec) :: cosu, sinu, cosv, sinv,
     1               cosmu(0:mpol-1), sinmu(0:mpol-1),
     2               cosnv(0:ntor),  sinnv(0:ntor), 
     3               cosnvn(0:ntor), sinnvn(0:ntor)
      REAL(rprec) :: work1(0:mpol-1,12)
      LOGICAL :: lru, lrv, lzu, lzv
!-------------------------------------------------------------------------------!
!     COMPUTES THE CYLINDRICAL COORDINATES R11 and Z11
!     AT A FIXED FLUX COORDINATE POINT si, ui(theta), vi(N*phi)
!     WHERE phi = geometric toroidal angle (0 to 2pi), N = no. field periods
!
!     INPUT:
!     c_flux:          array of (si,ui,vi) values to convert to cylindrical coordinates
!     rzl_array:       array of (r, z, lambda) Fourier coefficients for all radial, m, n values
!     ns, mpol,ntor, ntmax:  radial, poloidal, toroidal, type (r,z,l) dimensions of rzl_array
!     mscale, nscale:  option scale factors. Use ONLY if rzl_array comes from within VMEC.
!                      If arising from WOUT file, mscale, nscale => 1 are not passed
!
!     OUTPUT:
!     r_cyl    :       R = r_cyl(1);  N*PHI = r_cyl(2);   Z = r_cyl(3)
!
!     OPTIONAL OUTPUT
!                      Ru = dR/du;    Rv = dR/dv = dR/dphi / N
!                      Zu = dZ/du;    Zv = dZ/dv = dZ/dphi / N
!
!     NOTE:            User is responsible for multiplying Rv, Zv by N to get phi derivatives
!
!-------------------------------------------------------------------------------
      
      iflag = -1
      si = c_flux(1);  ui = c_flux(2);  vi = c_flux(3)
      r_cyl(2) = vi

!      IF (si.lt.zero .or. si.gt.one) THEN
      IF (si .lt. zero) THEN
         WRITE(6, *)' In flx2cyl, s(flux) must be > 0'
         RETURN
      END IF

      lru = PRESENT(ru); lrv = PRESENT(rv)
      lzu = PRESENT(zu); lzv = PRESENT(zv)
      IF (lrv .and. .not. lthreed) rv = 0
      IF (lzv .and. .not. lthreed) zv = 0

!
!     FIND INTERPOLATED s VALUE AND COMPUTE INTERPOLATION WEIGHTS wlo, whi (for even m modes)
!     AND wlo_odd, whi_odd (for odd m modes).
!     FOR si <= 1, POINT IS INSIDE PLASMA;
!     FOR si > 1, TRY EXTRAPOLATION (WITH CONTINUOUS s, u DERIVATIVES INTO "vacuum" REGION
!
      hs1 = one/(ns-1)
      IF (si .le. one) THEN
 
         jslo = 1 + si/hs1
         jshi = jslo+1
         wlo = (hs1*(jshi-1) - si)/hs1
         whi = 1 - wlo
         IF (jslo .eq. ns) jshi = jslo

!
!     USE Rmn, Zmn ~ SQRT(s) FOR ODD-m MODES, SO INTERPOLATE Xmn/SQRT(s)
! 
         whi_odd = whi*SQRT(si/(hs1*(jshi-1)))
         IF (jslo .ne. 1) THEN
            wlo_odd = wlo*SQRT(si/(hs1*(jslo-1)))
         ELSE
            wlo_odd = 0
            whi_odd = SQRT(si/(hs1*(jshi-1)))
         END IF

      ELSE

         jshi = ns
         jslo = ns-1
         wlo  = -(si - 1)/hs1;    wlo_odd = wlo
         whi  = 1 - wlo;          whi_odd = whi

      ENDIF


      mpol1 = mpol-1

      wmins(:,0:mpol1:2) = wlo
      wplus(:,0:mpol1:2) = whi
      wmins(:,1:mpol1:2) = wlo_odd
      wplus(:,1:mpol1:2) = whi_odd



      IF (.not.lasym) THEN
         IF (lthreed) THEN
            IF (ntmax .ne. 2) STOP 'ntmax != 2 in flx2cyl!'
            rss = 2;  zcs = 2
         ELSE
            IF (ntmax .ne. 1) STOP 'ntmax != 1 in flx2cyl!'
         END IF
      ELSE
         IF (lthreed) THEN
             IF (ntmax .ne. 4) STOP 'ntmax != 4 in flx2cyl!'
             rss = 2;  rsc = 3;  rcs = 4
             zcs = 2;  zcc = 3;  zss = 4
         ELSE
             IF (ntmax .ne. 2) STOP 'ntmax != 2 in flx2cyl!'
             rsc = 2;  zcc = 2
         END IF
      END IF

      zsc = 1+ntmax; zcs = zcs+ntmax; zcc = zcc+ntmax; zss = zss+ntmax

      rmncc = wmins*rzl_array(jslo,:,:,rcc) 
     1      + wplus*rzl_array(jshi,:,:,rcc)        !!COS(mu) COS(nv)
      zmnsc = wmins*rzl_array(jslo,:,:,zsc) 
     1      + wplus*rzl_array(jshi,:,:,zsc)        !!SIN(mu) COS(nv)

      IF (lthreed) THEN
         rmnss = wmins*rzl_array(jslo,:,:,rss) 
     1         + wplus*rzl_array(jshi,:,:,rss)     !!SIN(mu) SIN(nv)
         zmncs = wmins*rzl_array(jslo,:,:,zcs)
     1         + wplus*rzl_array(jshi,:,:,zcs)     !!COS(mu) SIN(nv)
      END IF

!
!     SETUP TRIG ARRAYS
!
      cosu = COS(ui);   sinu = SIN(ui)
      cosv = COS(vi);   sinv = SIN(vi)

      cosmu(0) = 1;    sinmu(0) = 0
      cosnv(0) = 1;    sinnv(0) = 0
      DO m = 1, mpol1
         cosmu(m) = cosmu(m-1)*cosu - sinmu(m-1)*sinu
         sinmu(m) = sinmu(m-1)*cosu + cosmu(m-1)*sinu
      END DO

      IF (PRESENT(mscale)) THEN
         cosmu = cosmu*mscale;  sinmu = sinmu*mscale
      END IF

      DO n = 1, ntor
         cosnv(n) = cosnv(n-1)*cosv - sinnv(n-1)*sinv
         sinnv(n) = sinnv(n-1)*cosv + cosnv(n-1)*sinv
      END DO

      IF (PRESENT(nscale)) THEN
         cosnv = cosnv*nscale;  sinnv = sinnv*nscale
      END IF

      IF (lrv .or. lzv) THEN
         DO n = 0, ntor
            cosnvn(n) = n*cosnv(n)
            sinnvn(n) =-n*sinnv(n)
         END DO
      END IF

      iflag = 0

!
!     COMPUTE R11, Z11 IN REAL SPACE
!
!     FIRST, INVERSE TRANSFORM IN N-V SPACE, FOR FIXED M
!
      DO m = 0, mpol1
 
         work1(m,1) = SUM(rmncc(:,m)*cosnv(:))
         work1(m,2) = SUM(zmnsc(:,m)*cosnv(:))
         IF (lru) work1(m,3) =-m*work1(m,1)
         IF (lzu) work1(m,4) = m*work1(m,2)
         IF (lthreed) THEN
            IF (lrv) work1(m,5) = SUM(rmncc(:,m)*sinnvn(:))
            IF (lzv) work1(m,6) = SUM(zmnsc(:,m)*sinnvn(:))
            work1(m,7) = SUM(rmnss(:,m)*sinnv(:))
            work1(m,8) = SUM(zmncs(:,m)*sinnv(:))
            IF (lru) work1(m,9) = m*work1(m,7)
            IF (lzu) work1(m,10) =-m*work1(m,8)
            IF (lrv) work1(m,11) = SUM(rmnss(:,m)*cosnvn(:))
            IF (lzv) work1(m,12) = SUM(zmncs(:,m)*cosnvn(:))
         END IF

      END DO

!
!     NEXT, INVERSE TRANSFORM IN M-U SPACE
!
      IF (lthreed) THEN
         r11 = SUM(work1(:,1)*cosmu(:) + work1(:,7)*sinmu(:))
         z11 = SUM(work1(:,2)*sinmu(:) + work1(:,8)*cosmu(:))
         IF (lru) ru = SUM(work1(:,3)*sinmu(:) + work1(:,9)*cosmu(:))
         IF (lzu) zu = SUM(work1(:,4)*cosmu(:) + work1(:,10)*sinmu(:))
         IF (lrv) rv = SUM(work1(:,5)*cosmu(:) + work1(:,11)*sinmu(:))
         IF (lzv) zv = SUM(work1(:,6)*sinmu(:) + work1(:,12)*cosmu(:))
      ELSE          
         r11 = SUM(work1(:,1)*cosmu(:))
         z11 = SUM(work1(:,2)*sinmu(:))
         IF (lru) ru = SUM(work1(:,3)*sinmu(:))
         IF (lzu) zu = SUM(work1(:,4)*cosmu(:))
      END IF


      IF (.not.lasym) GOTO 1000

      rmnsc = wmins*rzl_array(jslo,:,:,rsc) 
     1      + wplus*rzl_array(jshi,:,:,rsc)        !!SIN(mu) COS(nv)
      zmncc = wmins*rzl_array(jslo,:,:,zcc) 
     1      + wplus*rzl_array(jshi,:,:,zcc)        !!COS(mu) COS(nv)

      IF (lthreed) THEN
         rmncs = wmins*rzl_array(jslo,:,:,rcs) 
     1         + wplus*rzl_array(jshi,:,:,rcs)     !!COS(mu) SIN(nv)
         zmnss = wmins*rzl_array(jslo,:,:,zss)
     1         + wplus*rzl_array(jshi,:,:,zss)     !!SIN(mu) SIN(nv)
      END IF

!
!     COMPUTE R11, Z11 IN REAL SPACE
!
!     FIRST, INVERSE TRANSFORM IN N-V SPACE, FOR FIXED M
!
      DO m = 0, mpol1
 
         work1(m,1) = SUM(rmnsc(:,m)*cosnv(:))
         work1(m,2) = SUM(zmncc(:,m)*cosnv(:))
         IF (lru) work1(m,3) = m*work1(m,1)
         IF (lzu) work1(m,4) =-m*work1(m,2)

         IF (lthreed) THEN
            IF (lrv) work1(m,5) = SUM(rmnsc(:,m)*sinnvn(:))
            IF (lzv) work1(m,6) = SUM(zmncc(:,m)*sinnvn(:))
            work1(m,7) = SUM(rmncs(:,m)*sinnv(:))
            work1(m,8) = SUM(zmnss(:,m)*sinnv(:))
            IF (lru) work1(m,9) =-m*work1(m,7)
            IF (lzu) work1(m,10) = m*work1(m,8)
            IF (lrv) work1(m,11) = SUM(rmncs(:,m)*cosnvn(:))
            IF (lzv) work1(m,12) = SUM(zmnss(:,m)*cosnvn(:))
         END IF

      END DO

!
!     NEXT, INVERSE TRANSFORM IN M-U SPACE
!
      IF (lthreed) THEN
         r11 = r11 + SUM(work1(:,1)*sinmu(:) + work1(:,7)*cosmu(:))
         z11 = z11 + SUM(work1(:,2)*cosmu(:) + work1(:,8)*sinmu(:))
         IF (lru) ru = ru + 
     1                 SUM(work1(:,3)*cosmu(:) + work1(:,9)*sinmu(:))
         IF (lzu) zu = zu + 
     1                 SUM(work1(:,4)*sinmu(:) + work1(:,10)*cosmu(:))
         IF (lrv) rv = rv + 
     1                 SUM(work1(:,5)*sinmu(:) + work1(:,11)*cosmu(:))
         IF (lzv) zv = zv +
     1                 SUM(work1(:,6)*cosmu(:) + work1(:,12)*sinmu(:))
      ELSE          
         r11 = r11 + SUM(work1(:,1)*sinmu(:))
         z11 = z11 + SUM(work1(:,2)*cosmu(:))
         IF (lru) ru = ru + SUM(work1(:,3)*cosmu(:))
         IF (lzu) zu = zu + SUM(work1(:,4)*sinmu(:))
      END IF

 1000 CONTINUE

      r_cyl(1) = r11;  r_cyl(3) = z11

      END SUBROUTINE flx2cyl

!*******************************************************************************
!-------------------------------------------------------------------------------
!  convert from cylindrical to flux coordinates
!-------------------------------------------------------------------------------
      SUBROUTINE cyl2flx(rzl_in, r_cyl, c_flx, ns_in, ntor_in, mpol_in, 
     1      ntmax_in, lthreed_in, lasym_in, info, nfe, fmin, 
     1      mscale, nscale, ru, zu, rv, zv)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(out)       :: info, nfe
      INTEGER, INTENT(in)        :: ns_in, ntor_in, mpol_in, ntmax_in
      REAL(rprec), INTENT(in)    :: r_cyl(3)
      REAL(rprec), INTENT(inout) :: c_flx(3)
      REAL(rprec), INTENT(in), TARGET :: 
     1                 rzl_in(ns_in,0:ntor_in,0:mpol_in-1,3*ntmax_in)
      REAL(rprec), TARGET, OPTIONAL :: 
     1                 mscale(0:mpol_in-1), nscale(0:ntor_in)
      REAL(rprec), INTENT(out), OPTIONAL :: ru, zu, rv, zv
      REAL(rprec), INTENT(out)   :: fmin
      LOGICAL, INTENT(in) :: lthreed_in, lasym_in
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: nvar = 2
      REAL(rprec), PARAMETER :: ftol = 1.e-16_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec) :: xc_opt(nvar), r_cyl_out(3), fmin0
      INTEGER     :: iflag, itry, nfe_out

C-----------------------------------------------
C   Variables to communicate with internal subroutines newt2d and get_flxcoord
C-----------------------------------------------
      INTEGER     :: ns_loc, ntmax_loc, mpol_loc, ntor_loc
      REAL(rprec) :: r_target, phi_target, z_target, fnorm
      REAL(rprec), POINTER :: rzl_array(:,:,:,:)
      REAL(rprec), POINTER :: mscale_loc(:), nscale_loc(:)
      LOGICAL :: lthreed_loc, lasym_loc, lscale

C-----------------------------------------------
!     LOCAL PARAMETERS:
!     ftol    :   nominally, set to 1.E-16. Gives a maximum (relative)
!                 error in matching R and Z of sqrt(ftol), or 1.E-8. 
!                 To increase accuracy, ftol should be lowered, but this 
!                 may require more Newton iterations (slows code).
!       
!     INPUT:
!     rzl_in  :   4D array with r,z (lambda) Fourier coefficients vs. radius
!                 
!     r_cyl   :   vector specifying cylindrical point to match, R = r_cyl(1), 
!                 N*phi = r_cyl(2), Z = r_cyl(3)
!                 NOTE: N*phi (N=no. field periods) is input, NOT phi!
!     ns_in   :   number of radial nodes in input array rzl_in
!     ntmax_in:   number of different types of modes (cos-cos, sin-sin, etc)
!     mpol_in :   number of poloidal modes (0:mpol_in-1)
!     ntor_in :   number of toroidal modes = ntor_in+1 (0:ntor)
!     lthreed_in :true if this is a 3D plasma
!     lasym_in:   true if this is an asymmetric plasma
!     mscale  (nscale) : 
!                 optional scaling arrays for cos, sin arrays. Used
!                 only if this routine is called from within VMEC.
!
!     OUTPUT:
!     nfe     :   number of function evaluations
!     info    :   = 0, no errors, -1, fmin > ftol (tolerance exceeded on output)
!                 = -3, s > 1 outside plasma, probably
!     fmin    :   minimum value of f = (r - Rin)**2 + (z - Zin)**2 at c_flx

!     INPUT/OUTPUT:
!     c_flx   :   array of flux coordinates (s = c_flx(1), u=theta= c_flx(2), 
!                 v = N*phi= c_flx(3))
!                 on input, initial guess (recommend magnetic axis if "cold" start)
!                 on output, s, u values corresponding to r_cyl
C-----------------------------------------------
!     Initialize global variables
      rzl_array => rzl_in
      lthreed_loc = lthreed_in;  lasym_loc = lasym_in
      mpol_loc = mpol_in;  ntor_loc = ntor_in
      ns_loc = ns_in; ntmax_loc = ntmax_in
      lscale = PRESENT(mscale)
      IF (lscale) THEN
         mscale_loc => mscale;  nscale_loc => nscale
      END IF
      r_target = r_cyl(1); phi_target = r_cyl(2);  z_target = r_cyl(3)

!     Initialize local variables
      xc_opt(1) = c_flx(1); xc_opt(2) = c_flx(2)

!     Avoid exact magnetic axis, which is singular point
      IF (c_flx(1) .eq. zero) xc_opt(1) = one/(ns_loc-1)   

      fnorm = r_target**2 + z_target**2
      IF (fnorm .lt. EPSILON(fnorm)) fnorm = 1
      fnorm = one/fnorm


      nfe = 0
      fmin0 = 1

      DO itry = 1, 4

         CALL newt2d(xc_opt, fmin, ftol, nfe_out, nvar, info)
         nfe = nfe + nfe_out

         IF (fmin.le.ftol .or. info.eq.-3) EXIT
!
!        JOG POINT (BY ROTATING ANGLE) TO IMPROVE CONVERGENCE
!
         IF (fmin .gt. 1.E-3*fmin0) THEN
            xc_opt(2) = xc_opt(2) + twopi/20
         ELSE 
            xc_opt(2) = xc_opt(2) + twopi/40
         END IF

         fmin0 = MIN(fmin, fmin0)
!        PRINT *,' ITRY = ', itry+1,' FMIN = ', fmin
            
      END DO
         
      c_flx(1) = xc_opt(1); c_flx(2) = xc_opt(2); c_flx(3) = phi_target
!SPH      IF (info.eq.0 .and. c_flx(1).gt.one) c_flx(1) = one

      c_flx(2) = MOD(c_flx(2), twopi)
      DO WHILE (c_flx(2) .lt. zero)
         c_flx(2) = c_flx(2) + twopi
      END DO

!
!     COMPUTE Ru, Zu, Rv, Zv IF REQUIRED
!
      IF ((PRESENT(ru) .or. PRESENT(zu) .or. 
     1     PRESENT(rv) .or. PRESENT(zv)) .and. info.eq.0)
     2    CALL flx2cyl(rzl_in, c_flx, r_cyl_out, ns_loc, ntor_loc, 
     3         mpol_loc, ntmax_loc, lthreed_loc, lasym_loc, 
     4         iflag, MSCALE=mscale, NSCALE=nscale, 
     5         RU=ru, ZU=zu, RV=rv, ZV=zv)

      CONTAINS ! internal subprograms newt2d and get_flxcoord
    
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      SUBROUTINE newt2d(xc_opt, fmin, ftol, nfe, nvar, iflag)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nvar
      INTEGER, INTENT(out) :: nfe, iflag
      REAL(rprec), INTENT(inout) :: xc_opt(nvar)
      REAL(rprec), INTENT(in)    :: ftol
      REAL(rprec), INTENT(out)   :: fmin
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, PARAMETER :: niter = 50
      INTEGER     :: ieval, isgt1
      REAL(rprec) :: c_flx(3), r_cyl_out(3), fvec(nvar), sflux, 
     1               uflux, eps0, eps, epu, xc_min(2), factor
      REAL(rprec) :: x0(3), xs(3), xu(3), dels, delu, tau, fmin0,
     1               ru1, zu1, edge_value, snew
!-------------------------------------------------------------------------------
!     INPUT/OUTPUT:
!     xc_opt:   sflux = xc_opt(1), uflux = xc_opt(2) are the toroidal flux
!               coordinate and poloidal angle coordinate, respectively
!     iflag:    = 0, successfully find s,u point
!               =-1, did not converge
!               =-3, sflux > 1, probably
!
!     LOCAL VARIABLES:
!     tau:      d(R,Z)/d(s,u) (Jacobian)
!     isgt1:    counter for number of times s>1

!     FIND FLUX COORDINATES (s,u) WHICH CORRESPOND TO ZERO OF TARGET FUNCTION
!
!            F == (R - R_TARGET)**2 + (Z - Z_TARGET)**2
!
!     FOR A GIVEN CYLINDRICAL COORDINATE POINT (R_TARGET, N*PHI=PHI_TARGET, Z_TARGET)
!
!     Reference:  S.E.Attenberger, W.A.Houlberg, S.P.Hirshman, J. Comp Phys 72 (1987) 435.
!
!     The algorithm used here modifies this slightly to improve "faltering" convergence
!     by choosing a steepest-descent path when the step size has been decreased sufficiently
!     without yielding a lower value of F.
!
!-------------------------------------------------------------------------------

      iflag = -1      
      eps0 = SQRT(EPSILON(eps))
      xc_min = xc_opt

      c_flx(3) = phi_target
      fmin0 = 1.E10_dp
      factor = 1
      nfe = 0
      edge_value = one + one/(ns_loc-1)
      isgt1 = 0

      DO ieval = 1, niter
         nfe = nfe + 1

         sflux = MAX(xc_opt(1), zero)
!        sflux = MIN(MAX(xc_opt(1), zero), one)
         uflux = xc_opt(2)
         c_flx(1) = sflux;  c_flx(2) = uflux

!        COMPUTE R,Z, Ru, Zu
         CALL get_flxcoord(x0, c_flx, ru=ru1, zu=zu1)
         xu(1) = ru1; xu(3) = zu1

!        MAKE SURE sflux IS LARGE ENOUGH
!        TO COMPUTE d(sqrt(s))/ds ACCURATELY NEAR ORIGIN
         IF (sflux .ge. 1000*eps0) THEN
            eps = eps0
         ELSE
            eps = eps0*sflux
         END IF

!        COMPUTE Rs, Zs NUMERICALLY
         eps = ABS(eps)
         IF (sflux .ge. 1-eps) eps = -eps
         c_flx(1) = sflux + eps
         CALL get_flxcoord(r_cyl_out, c_flx)
         xs = (r_cyl_out - x0)/eps
         c_flx(1) = sflux

         x0(1) = x0(1) - r_target
         x0(3) = x0(3) - z_target
         fmin = (x0(1)**2 + x0(3)**2)*fnorm

         IF (fmin .gt. fmin0) THEN
            factor = (2*factor)/3
            xc_opt = xc_min
!           REDIRECT ALONG STEEPEST-DESCENT PATH
            IF (6*factor .lt. one) THEN
               dels =-(x0(1)*xs(1) + x0(3)*xs(3))/(xs(1)**2 + xs(3)**2)
               delu =-(x0(1)*xu(1) + x0(3)*xu(3))/(xu(1)**2 + xu(3)**2)
            END IF
         ELSE
            fmin0 = fmin
            factor = 1
            xc_min = xc_opt

!           NEWTON STEP
            tau = xu(1)*xs(3) - xu(3)*xs(1)
            IF (ABS(tau) .le. ABS(eps)*r_target**2) THEN
               iflag = -2
               EXIT
            END IF
            dels = ( x0(1)*xu(3) - x0(3)*xu(1))/tau
            delu = (-x0(1)*xs(3) + x0(3)*xs(1))/tau
            IF (fmin .gt. 1.E-3_dp) THEN
               dels = dels/2; delu = delu/2
            END IF
 
         END IF
 
         IF (fmin .le. ftol) EXIT

         IF (ABS(dels) .gt. one)   dels = SIGN(one, dels)
!         IF (ABS(delu) .gt. twopi/2) delu = SIGN(twopi/2, delu)

         snew = xc_opt(1) + dels*factor
         IF (snew .lt. zero) THEN
            xc_opt(1) = -snew/2               !Prevents oscillations around origin s=0
            xc_opt(2) = xc_opt(2) + twopi/2 
            delu = -delu
!             factor = (-snew/2-xc_opt(1))/dels
!             xc_opt(1) = -snew/2
         ELSE
            xc_opt(1) = snew
         END IF
         xc_opt(2) = xc_opt(2) + delu*factor

         IF (xc_opt(1) .gt. edge_value) THEN
            isgt1 = isgt1+1
            IF (xc_opt(1) .gt. 2._dp) isgt1 = isgt1+1
            IF (isgt1 .gt. 5) EXIT
         END IF

      END DO

      IF (isgt1.gt.5) THEN
         iflag = -3
         xc_min = xc_opt
      ELSE IF (xc_min(1) .gt. edge_value) THEN
         iflag = -3
      ELSE IF (fmin0 .le. ftol) THEN
         iflag = 0
      ELSE
         iflag = -1
      END IF
      
      fmin = fmin0
      xc_opt = xc_min
      xc_opt(2) = MOD(xc_opt(2), twopi)

      END SUBROUTINE newt2d    ! Internal subroutine to cyl2flx

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      SUBROUTINE get_flxcoord(x1, c_flx, ru, zu)  ! Internal subroutine to cyl2flx
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), INTENT(out) :: x1(3)
      REAL(rprec), INTENT(in)  :: c_flx(3)
      REAL(rprec), INTENT(out), OPTIONAL :: ru, zu
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: iflag
C-----------------------------------------------
      IF (lscale) THEN
         CALL flx2cyl(rzl_array, c_flx, x1, ns_loc, ntor_loc, mpol_loc, 
     1              ntmax_loc, lthreed_loc, lasym_loc, iflag, 
     2              MSCALE=mscale_loc, NSCALE=nscale_loc, RU=ru, ZU=zu)
      ELSE
         CALL flx2cyl(rzl_array, c_flx, x1, ns_loc, ntor_loc, mpol_loc, 
     1              ntmax_loc, lthreed_loc, lasym_loc, iflag, 
     2              RU=ru, ZU=zu)
      END IF

      END SUBROUTINE get_flxcoord  ! Internal subroutine to cyl2flx

      END SUBROUTINE cyl2flx

      END MODULE cyl_flux
