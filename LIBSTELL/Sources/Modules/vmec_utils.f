      MODULE vmec_utils
      USE stel_kinds
      USE stel_constants, ONLY: twopi, one, zero
      IMPLICIT NONE

      INTEGER, PRIVATE     :: ns_loc, ntmax_loc, mpol_loc, ntor_loc
      REAL(rprec), PRIVATE :: r_target, phi_target, z_target, fnorm
      REAL(rprec), POINTER, PRIVATE :: rzl_array(:,:,:,:)
      REAL(rprec), POINTER, PRIVATE :: mscale_loc(:), nscale_loc(:)
      LOGICAL, PRIVATE :: lthreed_loc, lasym_loc, lscale
      PRIVATE :: newt2d, get_flxcoord
!
!     THIS MODULE CONTAINS USEFUL UTILITIES FOR PROCESSING VMEC 
!     DATA. MOST FUNCTIONS ARE OVERLOADED TO BE ABLE TO USE EITHER 
!     INTERNALLY DATA (LOCAL FROM WITHIN VMEC) OR DATA FROM WOUT FILE
!

!
!     OVERLOADED FUNCTIONS
!
      INTERFACE GetBcyl
          MODULE PROCEDURE GetBcyl_WOUT, GetBcyl_VMEC
      END INTERFACE

      INTERFACE GetAcyl
          MODULE PROCEDURE GetAcyl_WOUT
      END INTERFACE

      INTERFACE GetJcyl
          MODULE PROCEDURE GetJcyl_WOUT
      END INTERFACE

      INTERFACE MSE_pitch
          MODULE PROCEDURE MSE_pitch_WOUT
          MODULE PROCEDURE MSE_pitch_VMEC
      END INTERFACE

      CONTAINS

      SUBROUTINE GetBcyl_WOUT(R1, Phi, Z1, Br, Bphi, Bz, 
     1                        sflx, uflx, info)
      USE read_wout_mod, phi_wout=>phi, ns_w=>ns, ntor_w=>ntor,
     1     mpol_w=>mpol, ntmax_w=>ntmax, lthreed_w=>lthreed,
     2     lasym_w=>lasym
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, OPTIONAL, INTENT(out) :: info
      REAL(rprec), INTENT(in)  :: R1, Z1, Phi
      REAL(rprec), INTENT(out) :: Br, Bphi, Bz
      REAL(rprec), INTENT(out), OPTIONAL :: sflx, uflx
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: fmin_acceptable = 1.E-12_dp
      INTEGER     :: nfe, info_loc
      REAL(rprec) :: r_cyl(3), c_flx(3), fmin
      REAL(rprec) :: Ru1, Zu1, Rv1, Zv1
      REAL(rprec) :: bsupu1, bsupv1
C-----------------------------------------------
      IF (.not.lwout_opened) THEN
         WRITE(6, '(2a,/,a)')
     1   ' This form of GetBcyl can only be called if WOUT has been',
     2   ' previously opened!',' Try GetBcyl_VMEC form instead.'
         RETURN
      END IF

      CALL LoadRZL

!     Computes cylindrical components of the magnetic field, Br, Bphi, Bz,
!     at the specified cylindrical coordinate point (R1, Phi, Z1), where
!     Phi is the true geometric toroidal angle (NOT N*Phi)
!
!     INPUT
!     R1, Phi, Z1  : cylindrical coordinates at which evaluation is to take place
!     
!     OUTPUT
!     Br, Bphi, Bz : computed cylindrical components of B at input point
!     sflx, uflx   : computed flux and theta angle at the cylindrical point
!
!     1. Convert to point in flux-coordinates: cflux = (s, u, v=N*phi)
!        and evaluate Ru, Zu, Rv, Zv at that point
!
      r_cyl(1) = R1;  r_cyl(2) = nfp*Phi;  r_cyl(3) = Z1
      c_flx(1) = 0;   c_flx(2) = 0;        c_flx(3) = r_cyl(2)
      IF (PRESENT(sflx)) c_flx(1) = sflx
      IF (PRESENT(uflx)) c_flx(2) = uflx
      CALL cyl2flx(rzl_local, r_cyl, c_flx, ns_w, ntor_w, mpol_w, 
     1     ntmax_w, lthreed_w, lasym_w, info_loc, nfe, fmin, 
     2     RU=Ru1, ZU=Zu1, RV=Rv1, ZV=Zv1)
      Rv1 = nfp*Rv1;  Zv1 = nfp*Zv1

      IF (info_loc.eq.-1 .and. (fmin .le. fmin_acceptable)) info_loc = 0

      IF (PRESENT(info)) info = info_loc
      IF (info_loc .ne. 0) RETURN

      IF (PRESENT(sflx)) sflx = c_flx(1)  
      IF (PRESENT(uflx)) uflx = c_flx(2)

      IF (c_flx(1) .gt. 2) THEN
         Br = 0;  Bphi = 0;  Bz = 0
         RETURN
      ELSE IF (c_flx(1) .gt. one) THEN
         c_flx(1) = one
      END IF
!
!     2. Evaluate Bsupu, Bsupv at this point
!
      CALL tosuvspace (c_flx(1), c_flx(2), c_flx(3), 
     1                 BSUPU=bsupu1, BSUPV=bsupv1)
!
!     3. Form Br, Bphi, Bz
!
      Br   = Ru1*bsupu1 + Rv1*bsupv1
      Bphi = R1 *bsupv1
      Bz   = Zu1*bsupu1 + Zv1*bsupv1
      
      END SUBROUTINE GetBcyl_WOUT

      SUBROUTINE GetAcyl_WOUT(R1, Phi, Z1, Ar, Aphi, Az, 
     1                        sflx, uflx, info)
      USE read_wout_mod, phi_wout=>phi, ns_w=>ns, ntor_w=>ntor,
     1     mpol_w=>mpol, ntmax_w=>ntmax, lthreed_w=>lthreed,
     2     lasym_w=>lasym, chi_wout=>chi, phipf_wout => phipf,
     3     isigng_w=>isigng
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, OPTIONAL, INTENT(out) :: info
      REAL(rprec), INTENT(in)  :: R1, Z1, Phi
      REAL(rprec), INTENT(out) :: Ar, Aphi, Az
      REAL(rprec), INTENT(out), OPTIONAL :: sflx, uflx
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: fmin_acceptable = 1.E-12_dp
      INTEGER     :: nfe, info_loc, js_lo, js_hi
      REAL(rprec) :: r_cyl(3), c_flx(3), fmin
      REAL(rprec) :: Ru1, Zu1, Rv1, Zv1
      REAL(rprec) :: Rs1, Zs1, c_flx2(3), r_cyl2(3)
      REAL(rprec) :: lam1, ds
      REAL(rprec) :: g11,g12,g13,g22,g23,g33
      REAL(rprec) :: g11i,g12i,g13i,g22i,g23i,g33i
      REAL(rprec) :: phi_flux, chi_flux,wegt, phip_flux
      REAL(rprec) :: asubs1, asubu1, asubv1, gdet
      REAL(rprec) :: asups1, asupu1, asupv1
C-----------------------------------------------
      IF (.not.lwout_opened) THEN
         WRITE(6, '(2a,/,a)')
     1   ' This form of GetBcyl can only be called if WOUT has been',
     2   ' previously opened!',' Try GetBcyl_VMEC form instead.'
         RETURN
      END IF

      CALL LoadRZL

!     Computes cylindrical components of the vector potential field, Ar, Aphi, Az,
!     at the specified cylindrical coordinate point (R1, Phi, Z1), where
!     Phi is the true geometric toroidal angle (NOT N*Phi)
!
!     INPUT
!     R1, Phi, Z1  : cylindrical coordinates at which evaluation is to take place
!     
!     OUTPUT
!     Ar, Aphi, Az : computed cylindrical components of A at input point
!     sflx, uflx   : computed flux and theta angle at the cylindrical point
!
!     1. Convert to point in flux-coordinates: cflux = (s, u, v=N*phi)
!        and evaluate Ru, Zu, Rv, Zv at that point
!
      r_cyl(1) = R1;  r_cyl(2) = nfp*Phi;  r_cyl(3) = Z1
      c_flx(1) = 0;   c_flx(2) = 0;        c_flx(3) = r_cyl(2)
      IF (PRESENT(sflx)) c_flx(1) = sflx
      IF (PRESENT(uflx)) c_flx(2) = uflx
      CALL cyl2flx(rzl_local, r_cyl, c_flx, ns_w, ntor_w, mpol_w, 
     1     ntmax_w, lthreed_w, lasym_w, info_loc, nfe, fmin, 
     2     RU=Ru1, ZU=Zu1, RV=Rv1, ZV=Zv1)
      Rv1 = nfp*Rv1;  Zv1 = nfp*Zv1

      IF (info_loc.eq.-1 .and. (fmin .le. fmin_acceptable)) info_loc = 0

      IF (PRESENT(info)) info = info_loc
      IF (info_loc .ne. 0) RETURN

      IF (PRESENT(sflx)) sflx = c_flx(1)  
      IF (PRESENT(uflx)) uflx = c_flx(2)

      IF (c_flx(1) .gt. one) THEN
         Ar = 0;  Aphi = 0;  Az = 0
         RETURN
      END IF
!
!     2. Interpolate Chi and phi
!        Note that arrays are indexed from 1:ns
!
      js_lo = FLOOR(c_flx(1)*(ns_w-1))
      js_hi = js_lo+1
      IF (js_hi > ns_w) THEN
        js_lo = ns_w -1
        js_hi = ns_w
      ENDIF
      wegt  = c_flx(1)*ns_w - js_lo
      phi_flux =   (1.0-wegt)*phi_wout(js_lo)
     1            + wegt*phi_wout(js_hi)
      chi_flux =   (1.0-wegt)*chi_wout(js_lo)
     1            + wegt*chi_wout(js_hi)
      phip_flux =  (1.0-wegt)*phipf_wout(js_lo)
     1            + wegt*phipf_wout(js_hi)

!
!     3. Get lambda and radial derivatives
!
      CALL tosuvspace(c_flx(1),c_flx(2),c_flx(3),LAM=lam1)
      ds = 0.25/ns_w
      c_flx2 = c_flx
      c_flx2(1) = c_flx2(1) + ds
      CALL flx2cyl(rzl_local,c_flx2,r_cyl2,ns_w,ntor_w, mpol_w,
     1             ntmax_w, lthreed_w, lasym_w, info_loc)
      Rs1 = (r_cyl2(1)-r_cyl(1))/ds
      Zs1 = (r_cyl2(3)-r_cyl(3))/ds
!      Rs1 = Rs1/phip_flux
!      Zs1 = Zs1/phip_flux

!
!     2. Evaluate Metric Elements
!       g_ij = e_i * e_j = e^i * e^j
!           e^i = e_j x e_k / sqrt(g)
!           e_i = d X / di
!     g21 = g12
!     g31 = g13
!     g32 = g23
!
      g11 = Rs1*Rs1         + Zs1*Zs1
      g12 = Rs1*Ru1         + Zs1*Zu1
      g13 = Rs1*Rv1         + Zs1*Zv1
      g22 = Ru1*Ru1         + Zu1*Zu1
      g23 = Ru1*Rv1         + Zu1*Zv1
      g33 = Rv1*Rv1 + R1*R1 + Zv1*Zv1
      gdet =  g11*(g22*g33-g23*g23)
     1       -g12*(g12*g33-g23*g13)
     2       +g13*(g12*g23-g22*g13)
      g11i = g22*g33 - g23*g23
      g12i = g13*g23 - g12*g33
      g13i = g12*g23 - g13*g22
      g22i = g11*g33 - g13*g13
      g23i = g13*g12 - g11*g23
      g33i = g11*g22 - g12*g12
!
!     3. Calculate the Fields
!        A = A_i * e^i
!           e^i = e_j x e_k / sqrt(g)
!           e_i = d X / di
!             g = det(g_ij)
!
      asubs1 = -lam1 * phip_flux*isigng_w    
      asubu1 =  phi_flux*isigng_w  
      asubv1 = -chi_flux*isigng_w 
      asups1 = (asubs1*g11i+asubu1*g12i+asubv1*g13i)/gdet
      asupu1 = (asubs1*g12i+asubu1*g22i+asubv1*g23i)/gdet
      asupv1 = (asubs1*g13i+asubu1*g23i+asubv1*g33i)/gdet
      Ar   = asups1*Rs1 + asupu1*Ru1 + asupv1*Rv1
      Aphi = asupv1*R1
      Az   = asups1*Zs1 + asupu1*Zu1 + asupv1*Zv1

!      Ar   = -asubs1*(R1*Zu1)          + asubu1*(R1*Zs1)
!      Az   =  asubs1*(R1*Ru1)          - asubu1*(R1*Rs1)
!      Aphi =  asubs1*(Ru1*Zv1-Rv1*Zu1) + asubu1*(Rs1*Zv1-Rv1*Zs1)
!     1         + asubv1*(Ru1*Zs1-Rs1*Zu1)
!      Ar   = Ar/(gdet)
!      Az   = Az/(gdet)
!      Aphi = Aphi/(gdet)
!      Aphi = gdet
!      Ar   = Ar/(gdet*phi_wout(ns_w))
!      Az   = Az/(gdet*phi_wout(ns_w))
!      Aphi   = Aphi/(gdet*phi_wout(ns_w))
!      Ar   = Ar/sqrt(gdet)
!      Az   = Az/sqrt(gdet)
!      Aphi = Aphi/sqrt(gdet)
      
      END SUBROUTINE GetAcyl_WOUT


      SUBROUTINE GetBcyl_VMEC(R1, Phi, Z1, Br, Bphi, Bz, sflx, uflx, 
     1     bsupu, bsupv, rzl_array, ns_in, ntor_in, mpol_in, ntmax_in, 
     2     nzeta, ntheta3, nper, mscale, nscale, lthreed_in, lasym_in,  
     3     info)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: ns_in, ntor_in, mpol_in, ntmax_in, 
     1                       nzeta, ntheta3, nper
      INTEGER, OPTIONAL, INTENT(out) :: info
      LOGICAL, INTENT(in) :: lthreed_in, lasym_in
      REAL(rprec), INTENT(in)  :: R1, Z1, Phi
      REAL(rprec), INTENT(in)  :: 
     1             rzl_array(ns_in,0:ntor_in,0:mpol_in-1,2*ntmax_in),
     2             mscale(0:mpol_in-1), nscale(0:ntor_in)
      REAL(rprec), DIMENSION(ns_in,nzeta,ntheta3), INTENT(in) 
     1                         :: bsupu, bsupv
      REAL(rprec), INTENT(out) :: Br, Bphi, Bz, sflx, uflx
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: c1p5 = 1.5_dp
      REAL(rprec), PARAMETER :: fmin_acceptable = 1.E-12_dp
      INTEGER     :: nfe, info_loc, jslo, jshi, julo, juhi, 
     1               kvlo, kvhi, ntheta1
      REAL(rprec) :: r_cyl(3), c_flx(3), vflx, vflx_norm, 
     1               uflx_norm, fmin
      REAL(rprec) :: wgt_s, wgt_u, wgt_v, hs1, hu1, hv1
      REAL(rprec) :: Ru1, Zu1, Rv1, Zv1
      REAL(rprec) :: bsupu1, bsupv1, bsupu2, bsupv2
C-----------------------------------------------

!     Computes cylindrical components of the magnetic field, Br, Bphi, Bz,
!     at the specified cylindrical coordinate point (R1, Phi, Z1), where
!     Phi is the true geometric toroidal angle (NOT NPER*Phi)
!     Also, sflx, uflx are the computed flux and theta angle at the point
!
!     This routine is callable from within the VMEC code (in contrast to
!     the GetBcyl routine, which requires WOUT output file).
!
!     1. Convert to point in flux-coordinates: cflux = (s, u, v=N*phi)
!        and evaluate Ru, Zu, Rv, Zv at that point
!
      r_cyl(1) = R1;  r_cyl(2) = nper*Phi;  r_cyl(3) = Z1
      c_flx(1) = 0;   c_flx(2) = 0;         c_flx(3) = r_cyl(2)
      CALL cyl2flx(rzl_array, r_cyl, c_flx, ns_in, ntor_in, mpol_in, 
     1     ntmax_in, lthreed_in, lasym_in, info_loc, nfe, fmin, 
     2     mscale, nscale, RU=Ru1, ZU=Zu1, RV=Rv1, ZV=Zv1)
      Rv1 = nper*Rv1;  Zv1 = nper*Zv1

      IF (info_loc.eq.-1 .and. (fmin .le. fmin_acceptable)) info_loc = 0

      IF (PRESENT(info)) info = info_loc
      IF (info_loc .ne. 0) RETURN

      sflx = c_flx(1);  uflx = c_flx(2);  vflx = c_flx(3)
      IF (c_flx(1) .gt. one) THEN
         Br = 0;  Bphi = 0;  Bz = 0
         RETURN
      END IF

!
!     2. Evaluate Bsupu, Bsupv at this flux coordinate point by 2D interpolation in s, u space
!        This is not quite as accurate as the 1D (s) interpolation based on the Fourier coefficients
!        of bsupu, bsupv...done in GetBcyl...
!        Formula 25.2.66 (Bivariate, 4pt Formula) in Abramowitz and Stegun
!
      hs1 = one/(ns_in - 1)
      jslo = INT(c1p5 + sflx/hs1)
      jshi = jslo+1
      wgt_s = (sflx - hs1*(jslo-c1p5))/hs1
      IF (jslo .eq. ns_in) THEN
!        USE Xhalf(ns+1) = 2*Xhalf(ns) - Xhalf(ns-1) FOR "GHOST" POINT VALUE hs/2 OUTSIDE EDGE
!        THEN, X = wlo*Xhalf(ns) + whi*Xhalf(ns+1) == Xhalf(ns) + whi*(Xhalf(ns) - Xhalf(ns-1)) 
!        WHERE wlo = 1 - wgt_s, whi = wgt_s
         jshi = jslo-1
         wgt_s = 1+wgt_s
      ELSE IF (jslo .eq. 1) THEN
         jslo = 2
      END IF

      IF (lasym_in) THEN
         ntheta1 = ntheta3
      ELSE
         ntheta1 = 2*(ntheta3 - 1)
      END IF
      
      uflx = MOD(uflx, twopi)
      DO WHILE (uflx .lt. zero) 
         uflx = uflx+twopi
      END DO

      hu1 = one/ntheta1
      uflx_norm = uflx/twopi
      julo = INT(1 + uflx_norm/hu1)
      IF (julo .gt. ntheta3) THEN
         IF (ABS(uflx_norm - 1) .lt. 1.E-2*hu1) THEN
            julo = 1
            uflx_norm = 0
         ELSE IF (ABS(uflx_norm - .5_dp) .lt. 1.E-2*hu1) THEN
            julo = ntheta3
            uflx_norm = .5_dp
         ELSE
            PRINT *, 'julo=', julo,' > ntheta3=', ntheta3,
     1      ' uflx_norm=', uflx_norm, ' in GetBcyl!'
            IF (PRESENT(info)) info = -10
            RETURN
         END IF
      END IF
      juhi = julo + 1
      IF (julo .eq. ntheta3) juhi = 1         !Periodic point at u = 0
      wgt_u = (uflx_norm - hu1*(julo-1))/hu1

      
      DO WHILE (vflx .lt. zero) 
         vflx = vflx+twopi
      END DO
      vflx = MOD(vflx, twopi)
      hv1 = one/nzeta
      vflx_norm = vflx/twopi
      kvlo = INT(1 + vflx_norm/hv1)
      kvhi = kvlo+1
      IF (kvlo .eq. nzeta) kvhi = 1
      wgt_v = (vflx_norm - hv1*(kvlo-1))/hv1

!
!     BIVARIATE INTERPOLATION IN S, U AT 2 kv PLANES
!
      bsupu1 = (1-wgt_s)*((1-wgt_u)*bsupu(jslo,kvlo,julo)
     2       +               wgt_u *bsupu(jslo,kvlo,juhi))
     1       + wgt_s*    ((1-wgt_u)*bsupu(jshi,kvlo,julo)
     3       +               wgt_u *bsupu(jshi,kvlo,juhi))

      bsupv1 = (1-wgt_s)*((1-wgt_u)*bsupv(jslo,kvlo,julo)
     2       +               wgt_u *bsupv(jslo,kvlo,juhi))
     1       + wgt_s*    ((1-wgt_u)*bsupv(jshi,kvlo,julo)
     3       +               wgt_u *bsupv(jshi,kvlo,juhi))

      bsupu2 = (1-wgt_s)*((1-wgt_u)*bsupu(jslo,kvhi,julo)
     2       +               wgt_u *bsupu(jslo,kvhi,juhi))
     1       + wgt_s*    ((1-wgt_u)*bsupu(jshi,kvhi,julo)
     3       +               wgt_u *bsupu(jshi,kvhi,juhi))

      bsupv2 = (1-wgt_s)*((1-wgt_u)*bsupv(jslo,kvhi,julo)
     2       +               wgt_u *bsupv(jslo,kvhi,juhi))
     1       + wgt_s*    ((1-wgt_u)*bsupv(jshi,kvhi,julo)
     3       +               wgt_u *bsupv(jshi,kvhi,juhi))

!
!     LINEAR INTERPOLATION IN V
!
      bsupu1 = (1-wgt_v)*bsupu1 + wgt_v*bsupu2
      bsupv1 = (1-wgt_v)*bsupv1 + wgt_v*bsupv2

!
!     3. Form Br, Bphi, Bz
!
      Br   = Ru1*bsupu1 + Rv1*bsupv1
      Bphi = R1 *bsupv1
      Bz   = Zu1*bsupu1 + Zv1*bsupv1
      
      END SUBROUTINE GetBcyl_VMEC


      SUBROUTINE GetJcyl_WOUT(R1, Phi, Z1, JR, JPHI, JZ, 
     1                        sflx, uflx, info)
      USE read_wout_mod, phi_wout1=>phi, ns_w1=>ns, ntor_w1=>ntor,
     1     mpol_w1=>mpol, ntmax_w1=>ntmax, lthreed_w1=>lthreed,
     2     lasym_w1=>lasym
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, OPTIONAL, INTENT(out) :: info
      REAL(rprec), INTENT(in)  :: R1, Z1, Phi
      REAL(rprec), INTENT(out) :: JR, JPHI, JZ
      REAL(rprec), INTENT(out), OPTIONAL :: sflx, uflx
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: fmin_acceptable = 1.E-6_dp
      INTEGER     :: nfe, info_loc
      REAL(rprec) :: r_cyl(3), c_flx(3), fmin
      REAL(rprec) :: Ru1, Zu1, Rv1, Zv1
      REAL(rprec) :: jsupu1, jsupv1, gsqrt1
C-----------------------------------------------
      IF (.not.lwout_opened) THEN
         WRITE(6, '(2a,/,a)')
     1   ' This form of GetBcyl can only be called if WOUT has been',
     2   ' previously opened!'
         RETURN
      END IF

      CALL LoadRZL

!     Computes cylindrical components of the current, Jr, Jphi, Jz,
!     at the specified cylindrical coordinate point (R1, Phi, Z1), where
!     Phi is the true geometric toroidal angle (NOT N*Phi)
!
!     INPUT
!     R1, Phi, Z1  : cylindrical coordinates at which evaluation is to take place
!     
!     OUTPUT
!     Br, Bphi, Bz : computed cylindrical components of B at input point
!     sflx, uflx   : computed flux and theta angle at the cylindrical point
!
!     1. Convert to point in flux-coordinates: cflux = (s, u, v=N*phi)
!        and evaluate Ru, Zu, Rv, Zv at that point
!
      r_cyl(1) = R1;  r_cyl(2) = nfp*Phi;  r_cyl(3) = Z1
      c_flx(1) = 0;   c_flx(2) = 0;        c_flx(3) = r_cyl(2)
      CALL cyl2flx(rzl_local, r_cyl, c_flx, ns_w1, ntor_w1, mpol_w1, 
     1     ntmax_w1, lthreed_w1, lasym_w1, info_loc, nfe, fmin, 
     2     RU=Ru1, ZU=Zu1, RV=Rv1, ZV=Zv1)
      Rv1 = nfp*Rv1;  Zv1 = nfp*Zv1

      IF (info_loc.eq.-1 .and. (fmin .le. fmin_acceptable)) info_loc = 0

      IF (PRESENT(info)) info = info_loc
      IF (info_loc .ne. 0) RETURN

      IF (PRESENT(sflx)) sflx = c_flx(1)  
      IF (PRESENT(uflx)) uflx = c_flx(2)

      IF (c_flx(1) .gt. one) THEN
         Jr = 0;  Jphi = 0;  Jz = 0
         RETURN
      END IF

!     3. Evaluate d(Bsubs)/du and d(Bsubs)/dv, d(Bsubu)/ds, d(Bsubv)/ds at this point
      CALL tosuvspace (c_flx(1), c_flx(2), c_flx(3), 
     1                 GSQRT=gsqrt1, JSUPU=jsupu1, JSUPV=jsupv1)

!      WRITE (36, '(1p4e12.4)') R1*jsupv1, dbsubuds1, dbsubsdu1, gsqrt1
!
!     4. Return Jr, Jphi, Jz
!
      Jr   = Ru1*jsupu1 + Rv1*jsupv1
      Jphi =              R1 *jsupv1
      Jz   = Zu1*jsupu1 + Zv1*jsupv1
      
      END SUBROUTINE GetJcyl_WOUT


      FUNCTION MSE_pitch_WOUT(r1, phi1, z1, acoef, efield, info)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(out) :: info
      REAL(rprec), INTENT(in) :: r1, phi1, z1, acoef(6)
      REAL(rprec), INTENT(in), OPTIONAL :: efield(2)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec) :: MSE_pitch_WOUT, Er, Ez, Br, Bphi, Bz
C-----------------------------------------------
!
!     Computes the Motional Stark Effect tan(pitch angle) == mse_pitch
!     at a given cylindrical point (R, phi, Z) inside the plasma
!
!     INPUT
!     Acoef : array of constants defined by the viewing geometry and beam velocity
!     r1, f1, z1: cylindrical coordinate of measurement point
!     Efield: (optional) electric field components (Er, Ez) in rest frame
!
!     OUTPUT
!     MSE_pitch  pitch angle at the input point
!     info       info = 0, calculation is valid
!
      IF (PRESENT(efield)) THEN
         Er = efield(1);  Ez = efield(2)
      ELSE
         Er = 0; Ez = 0
      END IF
      info = -1
!
!     Compute cylindrical components of B-field at given point R1, phi=f1, Z1
!
      CALL GetBcyl_WOUT(r1, phi1, z1, br, bphi, bz, INFO=info)

      MSE_pitch_WOUT = (acoef(1)*Bz   + acoef(5)*Er)/
     1                 (acoef(2)*Bphi + acoef(3)*Br 
     2               + acoef(4)*Bz   + acoef(6)*Ez)

      END FUNCTION MSE_pitch_WOUT

      FUNCTION MSE_pitch_VMEC(r1, phi1, z1, acoef, efield, sflx, uflx, 
     1     bsupu, bsupv, rzl_array, ns_in, ntor_in, mpol_in, ntmax_in, 
     2     nzeta, ntheta3, nper, mscale, nscale, lthreed_in, lasym_in,  
     3     info)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), INTENT(in) :: r1, phi1, z1, acoef(6), efield(2)
      INTEGER, INTENT(in) :: ns_in, ntor_in, mpol_in, ntmax_in, 
     1                       nzeta, ntheta3, nper
      LOGICAL, INTENT(in) :: lthreed_in, lasym_in
      REAL(rprec), INTENT(in)  :: 
     1             rzl_array(ns_in,0:ntor_in,0:mpol_in-1,2*ntmax_in),
     2             mscale(0:mpol_in-1), nscale(0:ntor_in)
      REAL(rprec), DIMENSION(ns_in,nzeta,ntheta3), INTENT(in) 
     1                         :: bsupu, bsupv
      REAL(rprec), INTENT(out) :: sflx, uflx
      INTEGER, INTENT(out) :: info
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec) :: MSE_pitch_VMEC, Er, Ez, Br, Bphi, Bz
C-----------------------------------------------
!
!     Computes the Motional Stark Effect tan(pitch angle) == mse_pitch
!     at a given cylindrical point (R, phi, Z) inside the plasma
!
!     INPUT
!     Acoef : array of constants defined by the viewing geometry and beam velocity
!     r1, f1, z1: cylindrical coordinate of measurement point
!     Efield: (optional) electric field components (Er, Ez) in rest frame
!
!     OUTPUT
!     MSE_pitch  pitch angle at the input point
!     info       info = 0, calculation is valid
!
      Er = efield(1);  Ez = efield(2)
      info = -1
!
!     Compute cylindrical components of B-field at given point R1, phi=f1, Z1
!
      CALL GetBcyl_VMEC(r1, phi1, z1, br, bphi, bz, sflx, uflx, 
     1     bsupu, bsupv, rzl_array, ns_in, ntor_in, mpol_in, ntmax_in, 
     2     nzeta, ntheta3, nper, mscale, nscale, lthreed_in, lasym_in,  
     3     info)

      MSE_pitch_VMEC = (acoef(1)*Bz   + acoef(5)*Er)/
     1                 (acoef(2)*Bphi + acoef(3)*Br 
     2               +  acoef(4)*Bz   + acoef(6)*Ez)

      END FUNCTION MSE_pitch_VMEC


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
C-----------------------------------------------
!
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
     1                 rzl_in(ns_in,0:ntor_in,0:mpol_in-1,2*ntmax_in)
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
     1     PRESENT(rv) .or. PRESENT(zv)) .and. info.eq.0) THEN
         IF (lscale) THEN
            CALL flx2cyl(rzl_in, c_flx, r_cyl_out, ns_loc, ntor_loc, 
     1         mpol_loc, ntmax_loc, lthreed_loc, lasym_loc, 
     2         iflag, MSCALE=mscale_loc, NSCALE=nscale_loc, 
     3         RU=ru, ZU=zu, RV=rv, ZV=zv)
         ELSE
            CALL flx2cyl(rzl_in, c_flx, r_cyl_out, ns_loc, ntor_loc, 
     1         mpol_loc, ntmax_loc, lthreed_loc, lasym_loc, 
     2         iflag, 
     3         RU=ru, ZU=zu, RV=rv, ZV=zv)
         END IF
      END IF

    
      END SUBROUTINE cyl2flx

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
      !INTEGER, PARAMETER :: niter = 50
      INTEGER, PARAMETER :: niter = 500
      INTEGER     :: ieval, isgt1
      REAL(rprec) :: c_flx(3), r_cyl_out(3), fvec(nvar), sflux, 
     1               uflux, eps0, eps, epu, xc_min(2), factor
      REAL(rprec) :: x0(3), xs(3), xu(3), dels, delu, tau, fmin0,
     1               ru1, zu1, edge_value, snew
C-----------------------------------------------
!
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
      iflag = -1      
      eps0 = SQRT(EPSILON(eps))
      xc_min = xc_opt

      c_flx(3) = phi_target
      fmin0 = 1.E10_dp
      factor = 1
      nfe = 0
      edge_value = one + one/(ns_loc-1)
      edge_value = 2*one
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
            !xc_opt(2) = xc_opt(2) + twopi/2 
            !delu = -delu
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

      END SUBROUTINE newt2d

      SUBROUTINE get_flxcoord(x1, c_flx, ru, zu)
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

      END SUBROUTINE get_flxcoord

      END MODULE vmec_utils
