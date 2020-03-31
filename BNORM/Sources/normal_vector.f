
      subroutine normal_vector(u, v)
      use stel_kinds
      use neswrite
      use normal_info
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec) :: v, u
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0, one = 1
      integer :: mn, m, n
      real(rprec) :: rbn, zbn, rub, rvb, zub, zvb,
     1    cosnv, sinnv, cosmn1, sinmn1, cosmu, sinmu, snx, sny
      real(rprec) :: x2, y2, v2, norm, twopi
C-----------------------------------------------

      twopi = 8*atan(one)

      rbn = zero;  rub = zero;   rvb = zero
      zbn = zero;  zub = zero;   zvb = zero

      do mn = 1, mnmax
         m = ixm(mn)
         n = ixn(mn)                          !!nfp factored out already, nescoil convention
         cosmu = cos(m*u)
         sinmu = sin(m*u)
         cosnv = cos(n*v)
         sinnv = sin(n*v)
         cosmn1 = cosmu*cosnv - sinmu*sinnv   !!cos(mu+nv), nescoil convention
         sinmn1 = sinmu*cosnv + cosmu*sinnv
         rbn = rbn + rmnc(mn) * cosmn1
     1             + rmns(mn) * sinmn1
         rub = rub - ixm(mn) * rmnc(mn) * sinmn1
     1             + ixm(mn) * rmns(mn) * cosmn1
         rvb = rvb - nfp*ixn(mn) * rmnc(mn) * sinmn1
     1             + nfp*ixn(mn) * rmns(mn) * cosmn1
         zbn = zbn + zmns(mn) * sinmn1
     1             + zmnc(mn) * cosmn1
         zub = zub + ixm(mn) * zmns(mn) * cosmn1
     1             - ixm(mn) * zmnc(mn) * sinmn1
         zvb = zvb + nfp*ixn(mn) * zmns(mn) * cosmn1
     1             - nfp*ixn(mn) * zmnc(mn) * sinmn1
      end do

!
!     un-normalized r,phi,z components of outward normal vector
!
      surf_norm(1) = rbn * zub
      surf_norm(2) = rub * zvb - rvb * zub
      surf_norm(3) =-rbn * rub

      norm = sqrt(sum(surf_norm**2))

      surf_norm(:) = surf_norm(:) / norm

      cosnv = cos(v/nfp)
      sinnv = sin(v/nfp)

      snx = surf_norm(1)*cosnv - surf_norm(2)*sinnv
      sny = surf_norm(1)*sinnv + surf_norm(2)*cosnv

      x2    = rbn*cosnv + coil_separation*snx
      y2    = rbn*sinnv + coil_separation*sny
      v2    = atan2(y2,x2)                                  !REAL phi on winding surface
      if(v2 .lt. zero .and. x2 .lt. zero) v2 = v2 + twopi   !Cut at v = +- pi

      v2 = v2*nfp
      vb_ws = v2

      rb_ws = sqrt (x2*x2 + y2*y2)
      zb_ws = zbn + coil_separation*surf_norm(3)

      end subroutine normal_vector
