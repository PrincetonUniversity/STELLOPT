!-----------------------------------------------------------------------
!     Module:        booz_ad
!     Purpose:       Exact tangent-linear (JVP) and adjoint (VJP) of the
!                    Boozer spectrum |B|_mn computed by boozer_coords.
!
!     Differentiated inputs on the half-mesh surface jrad:
!        bmodmnc, bsubumnc, bsubvmnc (Nyquist modes), lmns (VMEC modes),
!        hiota, and for lasym_b also bmodmns, bsubumns, bsubvmns, lmnc.
!     Differentiated outputs on that surface:
!        bmncb (and bmnsb for lasym_b), iota, G = gpsi, I = ipsi.
!     R and Z do not enter |B|_mn, iota, G or I, so their derivatives
!     vanish identically and they are not arguments here.
!
!     The derivatives are those of the discrete map that boozer_coords
!     evaluates (same VMEC angle grid, same quadrature weights and the
!     same normalisation scl), so FD of boozer_coords converges to them.
!     Linear stages (transpmn, vcoords_w, vcoords_rz) are reused as their
!     own tangents, their adjoints are the explicitly transposed sums,
!     and the pointwise harfun map and the final projection in boozer
!     are differentiated by hand.
!
!     Callers must have run read_wout_booz and boozer_setup first.
!-----------------------------------------------------------------------
      MODULE booz_ad
      USE stel_kinds, ONLY: rprec
      USE booz_params, ONLY: ns, nfp, mnmax, mnmax_nyq, mboz, nboz,     &
     &   mnboz, nunv, nv_boz, hiota, gpsi, ipsi, bsubumnc, bsubvmnc,    &
     &   bsubumns, bsubvmns, bmodmnc, bmodmns, ntorsum, lasym_b
      USE booz_persistent, ONLY: nu2_b, sfull, scl, thgrd, ztgrd, xm,   &
     &   xn, xm_nyq, xn_nyq, xmb, xnb, rmnc, zmns, lmns, rmns, zmnc,    &
     &   lmnc, cosm_b, sinm_b, cosn_b, sinn_b, cosm_nyq, sinm_nyq,      &
     &   cosn_nyq, sinn_nyq
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: boozer_coords_tl, boozer_coords_ad

      TYPE booz_ad_state
         REAL(rprec), ALLOCATABLE, DIMENSION(:) :: bmod, w, wt, wz,     &
     &      lam, lt, lz, uboz, vboz, xjac, psubu, psubv, qw
         REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: cosmm, sinmm,      &
     &      cosnn, sinnn
         REAL(rprec) :: g, i, iota, jacfac
      END TYPE booz_ad_state

      CONTAINS

!     Recompute the primal intermediates of boozer_coords at jrad with
!     the original routines.
      SUBROUTINE booz_ad_primal(jrad, st)
      INTEGER, INTENT(in) :: jrad
      TYPE(booz_ad_state), INTENT(out) :: st
      REAL(rprec), DIMENSION(mnmax_nyq) :: pmns_l, pmnc_l
      REAL(rprec), DIMENSION(nunv) :: r1, z1
      REAL(rprec), DIMENSION(ns) :: g_l, i_l
      INTEGER :: nparity, i1, imax

      ALLOCATE(st%bmod(nunv), st%w(nunv), st%wt(nunv), st%wz(nunv),    &
     &   st%lam(nunv), st%lt(nunv), st%lz(nunv), st%uboz(nunv),        &
     &   st%vboz(nunv), st%xjac(nunv), st%psubu(nunv), st%psubv(nunv), &
     &   st%qw(nunv), st%cosmm(nunv,0:mboz), st%sinmm(nunv,0:mboz),    &
     &   st%cosnn(nunv,0:nboz), st%sinnn(nunv,0:nboz))

      pmns_l = 0; pmnc_l = 0; g_l = 0; i_l = 0
      CALL transpmn(pmns_l, bsubumnc(1,jrad), bsubvmnc(1,jrad), pmnc_l, &
     &   bsubumns(1,jrad), bsubvmns(1,jrad), xm_nyq, xn_nyq, g_l, i_l,  &
     &   mnmax_nyq, jrad, lasym_b)
      st%g = g_l(jrad); st%i = i_l(jrad); st%iota = hiota(jrad)

      DO nparity = 0, 1
         CALL vcoords_rz(rmnc, zmns, lmns, rmns, zmnc, lmnc, xm, xn,    &
     &      ntorsum, ns, jrad, mnmax, r1, z1, st%lt, st%lz, st%lam,     &
     &      sfull, nparity, nunv, nfp, lasym_b)
      END DO

      CALL vcoords_w(bmodmnc(1,jrad), bmodmns(1,jrad), pmns_l, pmnc_l,  &
     &   xm_nyq, xn_nyq, jrad, mnmax_nyq, st%bmod, st%wt, st%wz, st%w,  &
     &   nunv, nfp, lasym_b)

      g_l(jrad) = st%g; i_l(jrad) = st%i
      CALL harfun(st%jacfac, hiota, g_l, i_l, jrad, nunv, st%lt, st%lz, &
     &   st%lam, st%wt, st%wz, st%w, st%uboz, st%vboz, st%xjac)
      st%psubu = (st%wt - st%i*st%lt)/st%jacfac
      st%psubv = (st%wz - st%i*st%lz)/st%jacfac

      CALL trigfunc(thgrd+st%uboz, ztgrd+st%vboz, st%cosmm, st%sinmm,   &
     &   st%cosnn, st%sinnn, mboz, nboz, nunv)

!     Quadrature weights of boozer: half weight on u=0 and u=pi rows
      st%qw = 1
      IF (.not.lasym_b) THEN
         i1 = nv_boz*(nu2_b-1)+1
         imax = i1-1+nv_boz
         st%qw(1:nv_boz) = 0.5_rprec
         st%qw(i1:imax) = 0.5_rprec
      END IF
      END SUBROUTINE booz_ad_primal

!     cos and sin of (m*theta_B - n*zeta_B) on the grid for Boozer mode mn
      SUBROUTINE booz_ad_phase(st, mn, cphi, sphi)
      TYPE(booz_ad_state), INTENT(in) :: st
      INTEGER, INTENT(in) :: mn
      REAL(rprec), DIMENSION(nunv), INTENT(out) :: cphi, sphi
      INTEGER :: m, n
      REAL(rprec) :: sgn
      m = NINT(xmb(mn))
      n = NINT(ABS(xnb(mn)))/nfp
      sgn = SIGN(1.0_rprec, xnb(mn))
      cphi = st%cosmm(:,m)*st%cosnn(:,n) + st%sinmm(:,m)*st%sinnn(:,n)*sgn
      sphi = st%sinmm(:,m)*st%cosnn(:,n) - st%cosmm(:,m)*st%sinnn(:,n)*sgn
      END SUBROUTINE booz_ad_phase

!     VMEC-grid basis functions for a VMEC (lnyq=F) or Nyquist mode
      SUBROUTINE booz_ad_vbasis(lnyq, mn, tcos, tsin, xmv, xnv)
      LOGICAL, INTENT(in) :: lnyq
      INTEGER, INTENT(in) :: mn
      REAL(rprec), DIMENSION(nunv), INTENT(out) :: tcos, tsin
      REAL(rprec), INTENT(out) :: xmv, xnv
      INTEGER :: m, n
      REAL(rprec) :: sgn
      IF (lnyq) THEN
         xmv = xm_nyq(mn); xnv = xn_nyq(mn)
      ELSE
         xmv = xm(mn); xnv = xn(mn)
      END IF
      m = NINT(xmv)
      n = NINT(ABS(xnv))/nfp
      sgn = SIGN(1.0_rprec, xnv)
      IF (lnyq) THEN
         tcos = cosm_nyq(:,m)*cosn_nyq(:,n) + sinm_nyq(:,m)*sinn_nyq(:,n)*sgn
         tsin = sinm_nyq(:,m)*cosn_nyq(:,n) - cosm_nyq(:,m)*sinn_nyq(:,n)*sgn
      ELSE
         tcos = cosm_b(:,m)*cosn_b(:,n) + sinm_b(:,m)*sinn_b(:,n)*sgn
         tsin = sinm_b(:,m)*cosn_b(:,n) - cosm_b(:,m)*sinn_b(:,n)*sgn
      END IF
      END SUBROUTINE booz_ad_vbasis

!-----------------------------------------------------------------------
!     Tangent-linear model.  Inputs d_* are tangents of the jrad columns
!     of bmodmnc, bsubumnc, bsubvmnc, lmns, of hiota(jrad) and (lasym_b)
!     of bmodmns, bsubumns, bsubvmns, lmnc.  Outputs are the tangents of
!     bmncb, bmnsb, iota, G and I on that surface.
!-----------------------------------------------------------------------
      SUBROUTINE boozer_coords_tl(jrad, d_bmnc, d_bsubu, d_bsubv,       &
     &   d_lmns, d_iota, d_bmns, d_bsubus, d_bsubvs, d_lmnc,            &
     &   d_bmncb, d_bmnsb, d_iota_out, d_g, d_i)
      INTEGER, INTENT(in) :: jrad
      REAL(rprec), DIMENSION(mnmax_nyq), INTENT(in) :: d_bmnc, d_bsubu,  &
     &   d_bsubv, d_bmns, d_bsubus, d_bsubvs
      REAL(rprec), DIMENSION(mnmax), INTENT(in) :: d_lmns, d_lmnc
      REAL(rprec), INTENT(in) :: d_iota
      REAL(rprec), DIMENSION(mnboz), INTENT(out) :: d_bmncb, d_bmnsb
      REAL(rprec), INTENT(out) :: d_iota_out, d_g, d_i
      TYPE(booz_ad_state) :: st
      REAL(rprec), DIMENSION(mnmax_nyq) :: d_pmns, d_pmnc
      REAL(rprec), DIMENSION(ns) :: g_l, i_l
      REAL(rprec), DIMENSION(nunv) :: d_bmod, d_w, d_wt, d_wz, d_lam,   &
     &   d_lt, d_lz, d_v, d_u, d_pu, d_pv, d_jac, r1, z1, cphi, sphi,   &
     &   bj, dbj
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: zr1, zr2, zr3, zr4,    &
     &   dl2s, dl2c
      REAL(rprec) :: dd, dem
      INTEGER :: mn, nparity

      CALL booz_ad_primal(jrad, st)

!     transpmn and vcoords_w are linear: apply them to the tangents
      d_pmns = 0; d_pmnc = 0; g_l = 0; i_l = 0
      CALL transpmn(d_pmns, d_bsubu, d_bsubv, d_pmnc, d_bsubus,         &
     &   d_bsubvs, xm_nyq, xn_nyq, g_l, i_l, mnmax_nyq, jrad, lasym_b)
      d_g = g_l(jrad); d_i = i_l(jrad)
      CALL vcoords_w(d_bmnc, d_bmns, d_pmns, d_pmnc, xm_nyq, xn_nyq,    &
     &   jrad, mnmax_nyq, d_bmod, d_wt, d_wz, d_w, nunv, nfp, lasym_b)

!     lambda synthesis of vcoords_rz is linear in lmns, lmnc
      ALLOCATE(zr1(mnmax,ns), zr2(mnmax,ns), zr3(mnmax,ns),             &
     &   zr4(mnmax,ns), dl2s(mnmax,ns), dl2c(mnmax,ns))
      zr1 = 0; zr2 = 0; zr3 = 0; zr4 = 0; dl2s = 0; dl2c = 0
      dl2s(:,jrad) = d_lmns
      IF (lasym_b) dl2c(:,jrad) = d_lmnc
      DO nparity = 0, 1
         CALL vcoords_rz(zr1, zr2, dl2s, zr3, zr4, dl2c, xm, xn,        &
     &      ntorsum,                                                    &
     &      ns, jrad, mnmax, r1, z1, d_lt, d_lz, d_lam, sfull, nparity, &
     &      nunv, nfp, lasym_b)
      END DO
      DEALLOCATE(zr1, zr2, zr3, zr4, dl2s, dl2c)

!     harfun, differentiated pointwise
      dem = 1/st%jacfac
      dd = d_g + d_iota*st%i + st%iota*d_i
      d_v  = dem*(d_w  - d_i*st%lam - st%i*d_lam - st%vboz*dd)
      d_pu = dem*(d_wt - d_i*st%lt  - st%i*d_lt  - st%psubu*dd)
      d_pv = dem*(d_wz - d_i*st%lz  - st%i*d_lz  - st%psubv*dd)
      d_u  = d_lam + d_iota*st%vboz + st%iota*d_v
      d_jac = d_lt*(1+st%psubv) + (1+st%lt)*d_pv                        &
     &      + (d_iota - d_lz)*st%psubu + (st%iota - st%lz)*d_pu

!     boozer projection
      bj  = st%qw*st%bmod*st%xjac
      dbj = st%qw*(d_bmod*st%xjac + st%bmod*d_jac)
      d_bmnsb = 0
      DO mn = 1, mnboz
         CALL booz_ad_phase(st, mn, cphi, sphi)
         d_bmncb(mn) = scl(mn)*SUM(dbj*cphi                             &
     &      - bj*sphi*(xmb(mn)*d_u - xnb(mn)*d_v))
         IF (lasym_b) d_bmnsb(mn) = scl(mn)*SUM(dbj*sphi                &
     &      + bj*cphi*(xmb(mn)*d_u - xnb(mn)*d_v))
      END DO
      d_iota_out = d_iota
      END SUBROUTINE boozer_coords_tl

!-----------------------------------------------------------------------
!     Adjoint model.  Given cotangents of bmncb, bmnsb, iota, G and I on
!     surface jrad, return the cotangents of the jrad columns of the
!     inputs listed for boozer_coords_tl (overwritten, not accumulated).
!-----------------------------------------------------------------------
      SUBROUTINE boozer_coords_ad(jrad, a_bmncb, a_bmnsb, a_iota_out,   &
     &   a_g_out, a_i_out, a_bmnc, a_bsubu, a_bsubv, a_lmns, a_iota,    &
     &   a_bmns, a_bsubus, a_bsubvs, a_lmnc)
      INTEGER, INTENT(in) :: jrad
      REAL(rprec), DIMENSION(mnboz), INTENT(in) :: a_bmncb, a_bmnsb
      REAL(rprec), INTENT(in) :: a_iota_out, a_g_out, a_i_out
      REAL(rprec), DIMENSION(mnmax_nyq), INTENT(out) :: a_bmnc, a_bsubu, &
     &   a_bsubv, a_bmns, a_bsubus, a_bsubvs
      REAL(rprec), DIMENSION(mnmax), INTENT(out) :: a_lmns, a_lmnc
      REAL(rprec), INTENT(out) :: a_iota
      TYPE(booz_ad_state) :: st
      REAL(rprec), DIMENSION(nunv) :: s1, s2m, s2n, cphi, sphi, ac, as, &
     &   a_bmod, a_jac, a_u, a_v, a_pu, a_pv, a_w, a_wt, a_wz, a_lam,   &
     &   a_lt, a_lz, tcos, tsin
      REAL(rprec), DIMENSION(mnmax_nyq) :: a_pmns, a_pmnc
      REAL(rprec) :: a_g, a_i, a_d, dem, xmv, xnv
      INTEGER :: mn

      CALL booz_ad_primal(jrad, st)

!     transpose of the boozer projection
      s1 = 0; s2m = 0; s2n = 0
      DO mn = 1, mnboz
         CALL booz_ad_phase(st, mn, cphi, sphi)
         ac = scl(mn)*a_bmncb(mn)*cphi
         as = -scl(mn)*a_bmncb(mn)*sphi
         IF (lasym_b) THEN
            ac = ac + scl(mn)*a_bmnsb(mn)*sphi
            as = as + scl(mn)*a_bmnsb(mn)*cphi
         END IF
         s1  = s1  + ac
         s2m = s2m + xmb(mn)*as
         s2n = s2n + xnb(mn)*as
      END DO
      a_bmod = st%qw*st%xjac*s1
      a_jac  = st%qw*st%bmod*s1
      a_u    = st%qw*st%bmod*st%xjac*s2m
      a_v    = -st%qw*st%bmod*st%xjac*s2n

!     transpose of the pointwise harfun map
      dem = 1/st%jacfac
      a_iota = a_iota_out + SUM(a_jac*st%psubu)
      a_lt = a_jac*(1+st%psubv)
      a_pv = a_jac*(1+st%lt)
      a_lz = -a_jac*st%psubu
      a_pu = a_jac*(st%iota - st%lz)
!     U = lam + iota*v
      a_lam = a_u
      a_iota = a_iota + SUM(a_u*st%vboz)
      a_v = a_v + st%iota*a_u
!     v, pu, pv = (w* - I*lam*)/D
      a_w  = dem*a_v
      a_wt = dem*a_pu
      a_wz = dem*a_pv
      a_lam = a_lam - st%i*dem*a_v
      a_lt  = a_lt  - st%i*dem*a_pu
      a_lz  = a_lz  - st%i*dem*a_pv
      a_i = a_i_out - dem*SUM(a_v*st%lam + a_pu*st%lt + a_pv*st%lz)
      a_d = -dem*SUM(a_v*st%vboz + a_pu*st%psubu + a_pv*st%psubv)
!     D = G + iota*I
      a_g = a_g_out + a_d
      a_iota = a_iota + a_d*st%i
      a_i = a_i + a_d*st%iota

!     transpose of vcoords_w (Nyquist synthesis of |B| and w)
      a_bmnc = 0; a_bmns = 0; a_pmns = 0; a_pmnc = 0
      DO mn = 1, mnmax_nyq
         CALL booz_ad_vbasis(.true., mn, tcos, tsin, xmv, xnv)
         a_bmnc(mn) = SUM(tcos*a_bmod)
         a_pmns(mn) = SUM(tsin*a_w + tcos*(xmv*a_wt - xnv*a_wz))
         IF (lasym_b) THEN
            a_bmns(mn) = SUM(tsin*a_bmod)
            a_pmnc(mn) = SUM(tcos*a_w - tsin*(xmv*a_wt - xnv*a_wz))
         END IF
      END DO

!     transpose of transpmn
      a_bsubu = 0; a_bsubv = 0; a_bsubus = 0; a_bsubvs = 0
      DO mn = 1, mnmax_nyq
         IF (NINT(xm_nyq(mn)) .ne. 0) THEN
            a_bsubu(mn) = a_pmns(mn)/xm_nyq(mn)
            IF (lasym_b) a_bsubus(mn) = -a_pmnc(mn)/xm_nyq(mn)
         ELSE IF (NINT(xn_nyq(mn)) .ne. 0) THEN
            a_bsubv(mn) = -a_pmns(mn)/xn_nyq(mn)
            IF (lasym_b) a_bsubvs(mn) = a_pmnc(mn)/xn_nyq(mn)
         ELSE
            a_bsubv(mn) = a_g
            a_bsubu(mn) = a_i
         END IF
      END DO

!     transpose of the lambda synthesis in vcoords_rz
      a_lmns = 0; a_lmnc = 0
      DO mn = 1, mnmax
         CALL booz_ad_vbasis(.false., mn, tcos, tsin, xmv, xnv)
         a_lmns(mn) = SUM(tsin*a_lam + tcos*(xmv*a_lt - xnv*a_lz))
         IF (lasym_b) a_lmnc(mn) =                                      &
     &      SUM(tcos*a_lam - tsin*(xmv*a_lt - xnv*a_lz))
      END DO
      END SUBROUTINE boozer_coords_ad

      END MODULE booz_ad
