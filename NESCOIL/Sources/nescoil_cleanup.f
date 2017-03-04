
! ----------------------------------------------------------------------
      subroutine nescoil_cleanup
!................................................................
c     Purpose:
c     deallocate everything at end of all nescoil runs
c................................................................
      USE Vsurface1, ONLY: x, y, z, r, nuv
      USE Vsurface2, ONLY: xu, yu, xv, ru, rv
      USE Vsurface3, ONLY: yv, zu, zv
      USE Vsurface4, ONLY: xuu, yuu, zuu, ruu
      USE Vsurface5, ONLY: xuv, yuv, zuv, ruv
      USE Vsurface6, ONLY: xvv, yvv, zvv, rvv
      USE Vsurface7, ONLY: snx, sny, snz
      USE Vsurface8, ONLY: xcur, ycur, zcur
      USE Vsurface9, ONLY: x1, y1, z1, r1
      USE Vsurfac10, ONLY: x1u, y1u, z1u
      USE Vsurfac11, ONLY: x1v, y1v, z1v
      USE Vsurfac12, ONLY: r1uu, z1uu, r1u
      USE Vsurfac13, ONLY: snx1, sny1, snz1, dsur1
      USE Vbiopre
      USE Vculine1
      USE Vprecal3, ONLY: factor, eps
      USE Vprecal4
      use Vprecal5
      USE Vprecal6, ONLY: comu, simu
      USE Vprecal7, ONLY: conv, sinv
      USE Vprecal8, ONLY: coh1, sih1
      USE Vprecal9, ONLY: comu1, simu1
      USE Vprecal10
      USE Vdiagno2, ONLY: pot
      USE Vdiagno3
      USE Vvacuum1, ONLY: cr, cz
      USE Vvacuum2, ONLY: cr1, cz1, cl1
      USE Vvacuum4, ONLY: cr2, cz2
      USE Vvacuum5, ONLY: cr3, cz3
      USE Vvacuum6, ONLY: cf, sf
      USE Vfourier2, ONLY: trigs, trigsv
      USE Vbnorm, ONLY: bn_ext
      use NumParams, ONLY: inesc
      implicit none
c................................................................

c      from nescoil:
      deallocate (vx, vy, vz, dx, dy, dz, xw, yw, zw, curre, cok, sik)

c      from read_nescoil_input:
      deallocate (x, y, z, r, xu, yu, xv, ru, rv, yv, zu, zv,
     1   xuu, yuu, zuu, ruu, xuv, yuv, zuv, ruv, xvv, yvv, zvv, rvv,
     2   snx, sny, snz, xcur, ycur, zcur)

c      from precal:
      deallocate (trigs, trigsv,
     1   conv1, sinv1, conv, sinv, coh, sih, comu, simu,
     2   factor, coh1, sih1, comu1, simu1, eps)

c      from surface_plas:
      deallocate (x1, y1, z1, r1, x1u, y1u, z1u, x1v, y1v, z1v,
     1   r1uu, z1uu, r1u, snx1, sny1, snz1, dsur1)

c      from surface_coil:
      deallocate (cr, cz, cr1, cz1,
     4   cl1, cr2, cz2, cr3, cz3, cf, sf)

c      from various places:
      deallocate (pote, pot, bn_ext)

      write (inesc, '(a)') '---- Deallocated all nescoil arrays ----'
      close (inesc)

      end subroutine nescoil_cleanup
