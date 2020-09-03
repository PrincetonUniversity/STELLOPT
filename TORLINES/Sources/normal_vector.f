
      subroutine normal_vector(u, v)
      use stel_kinds
      USE torlines_background
      USE torlines_runtime, ONLY: lvessel, vsurf
      USE wall_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec) :: v, u
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0, one = 1
      integer :: mn, m, n, ik
      real(rprec) :: rbn, zbn, rub, rvb, zub, zvb,
     1    cosnv, sinnv, cosmn1, sinmn1, cosmu, sinmu, snx, sny
      real(rprec) :: x2, y2, v2, norm, twopi, dist
      REAL(rprec) :: surf_norm(3)
C-----------------------------------------------

   
      twopi = 8*atan(one)

      rbn = zero;  rub = zero;   rvb = zero
      zbn = zero;  zub = zero;   zvb = zero

      do mn = 1, mnmax_m
         m = ixm_m(mn)
         n = ixn_m(mn)                          !!nfp factored out already, nescoil convention
         cosmu = cos(m*u)
         sinmu = sin(m*u)
         cosnv = cos(n*v)
         sinnv = sin(n*v)
         cosmn1 = cosmu*cosnv - sinmu*sinnv   !!cos(mu+nv), nescoil convention
         sinmn1 = sinmu*cosnv + cosmu*sinnv
         rbn = rbn + rmnc_m(mn) * cosmn1
     1             + rmns_m(mn) * sinmn1
         rub = rub - ixm_m(mn) * rmnc_m(mn) * sinmn1
     1             + ixm_m(mn) * rmns_m(mn) * cosmn1
         rvb = rvb - nfp*ixn_m(mn) * rmnc_m(mn) * sinmn1
     1             + nfp*ixn_m(mn) * rmns_m(mn) * cosmn1
         zbn = zbn + zmns_m(mn) * sinmn1
     1             + zmnc_m(mn) * cosmn1
         zub = zub + ixm_m(mn) * zmns_m(mn) * cosmn1
     1             - ixm_m(mn) * zmnc_m(mn) * sinmn1
         zvb = zvb + nfp*ixn_m(mn) * zmns_m(mn) * cosmn1
     1             - nfp*ixn_m(mn) * zmnc_m(mn) * sinmn1
      end do

!
!     un-normalized r,phi,z components of outward normal vector
!
      surf_norm(1) =  rbn * zub
      surf_norm(2) =  (rub * zvb - rvb * zub)
      surf_norm(3) = -rbn * rub

      surf_norm = isgn * surf_norm
      
      
      norm = sqrt(sum(surf_norm*surf_norm))

      surf_norm(:) = surf_norm(:) / norm

      cosnv = cos(v/nfp)
      sinnv = sin(v/nfp)

      snx = surf_norm(1)*cosnv - surf_norm(2)*sinnv
      sny = surf_norm(1)*sinnv + surf_norm(2)*cosnv

      x2    = rbn*cosnv + dr_m*snx
      y2    = rbn*sinnv + dr_m*sny
      zb_ws = zbn+dr_m*surf_norm(3)


      v2    = atan2(y2,x2)
      !if(v2 .lt. zero .and. x2 .lt. zero) v2 = v2 + twopi   !Cut at v = +- pi
      v2 = v2*nfp
      vb_ws = v2

      rb_ws = sqrt (x2*x2 + y2*y2)
      RETURN

      end subroutine normal_vector
