
! ----------------------------------------------------------------------
      subroutine surface_coil
!................................................................
c                                                             01/01/89
c     purpose:
c     Set up coil surface from fourier coeffs
c     Calculate dsur, normals etc
c
c................................................................
      use Vmeshes
      use NumParams
      use OutCtrl, ONLY: w_csurf
      use LoopCtrl
      USE Vvacuum3
      USE Vprecal1
      USE Vprecal5
      USE Vprecal6
      USE Vprecal7, ONLY: conv, sinv, nd
      USE Vvacuum1, ONLY: cr, cz, ms, ns
      USE Vsurface1, ONLY: x, y, z, r, nuv
      USE Vsurface2, ONLY: xu, yu, xv, ru, rv
      USE Vsurface3, ONLY: yv, zu, zv
      USE Vsurface4, ONLY: xuu, yuu, zuu, ruu
      USE Vsurface5, ONLY: xuv, yuv, zuv, ruv
      USE Vsurface6, ONLY: xvv, yvv, zvv, rvv
      USE Vsurface7, ONLY: snx, sny, snz
      USE Vsurface8, ONLY: xcur, ycur, zcur
      USE Vsurface9
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, m, n, kv, ku, k
      real(rprec) :: cm, cn, cofp, cofm, sifp, sifm
c ----------------------------------------------------------------------
c   outer surface   1

      i=0
      if (.not.allocated(x))
     1 allocate (x(nuv), y(nuv), z(nuv), r(nuv), xu(nuv), yu(nuv),
     1   xv(nuv), ru(nuv), rv(nuv), yv(nuv), zu(nuv), zv(nuv),
     2   xuu(nuv), yuu(nuv), zuu(nuv), ruu(nuv), xuv(nuv), yuv(nuv),
     3   zuv(nuv), ruv(nuv), xvv(nuv), yvv(nuv), zvv(nuv), rvv(nuv),
     4   snx(nuv), sny(nuv), snz(nuv), xcur(nuv), ycur(nuv), zcur(nuv),
     5   stat = i)
      if (i .ne. 0) stop 'allocation error in NESCOIL surface_coil'
      x(:) = zero
      y(:) = zero
      z(:) = zero
      r(:) = zero

      xu(:) = zero
      yu(:) = zero
      zu(:) = zero
      ru(:) = zero

      xv(:) = zero
      yv(:) = zero
      zv(:) = zero
      rv(:) = zero

      xuu(:) = zero
      yuu(:) = zero
      zuu(:) = zero
      ruu(:) = zero

      xuv(:) = zero
      yuv(:) = zero
      zuv(:) = zero
      ruv(:) = zero

      xvv(:) = zero
      yvv(:) = zero
      zvv(:) = zero
      rvv(:) = zero

      snx(:) = zero
      sny(:) = zero
      snz(:) = zero

      xcur(:) = zero
      ycur(:) = zero
      zcur(:) = zero

      do 20 m = 0,ms
      do 20 n = 0,ns
         cm     = m*pi2
         cn     = n*pi2
         i      = 0
      do 20 kv = 1 , nv
      do 20 ku = 1 , nu
         i      = i + 1
         cofp   = comu(ku,m)*conv(kv,n)-simu(ku,m)*sinv(kv,n)
         cofm   = comu(ku,m)*conv(kv,n)+simu(ku,m)*sinv(kv,n)
         sifp   = simu(ku,m)*conv(kv,n)+comu(ku,m)*sinv(kv,n)
         sifm   = simu(ku,m)*conv(kv,n)-comu(ku,m)*sinv(kv,n)
         r(i)   = r(i)    +           cr(m,n)*cofp + cr(m,-n)*cofm
         z(i)   = z(i)    +           cz(m,n)*sifp + cz(m,-n)*sifm
         ru(i)  = ru(i)   - cm *    ( cr(m,n)*sifp + cr(m,-n)*sifm )
         rv(i)  = rv(i)   - cn *    ( cr(m,n)*sifp - cr(m,-n)*sifm )
         zu(i)  = zu(i)   + cm *    ( cz(m,n)*cofp + cz(m,-n)*cofm )
         zv(i)  = zv(i)   + cn *    ( cz(m,n)*cofp - cz(m,-n)*cofm )
         ruu(i) = ruu(i)  - cm * cm*( cr(m,n)*cofp + cr(m,-n)*cofm )
         ruv(i) = ruv(i)  - cm * cn*( cr(m,n)*cofp - cr(m,-n)*cofm )
         rvv(i) = rvv(i)  - cn * cn*( cr(m,n)*cofp + cr(m,-n)*cofm )
         zuu(i) = zuu(i)  - cm * cm*( cz(m,n)*sifp + cz(m,-n)*sifm )
         zuv(i) = zuv(i)  - cm * cn*( cz(m,n)*sifp - cz(m,-n)*sifm )
         zvv(i) = zvv(i)  - cn * cn*( cz(m,n)*sifp + cz(m,-n)*sifm )
   20 continue
      do 50 kv = 1 , nv
      do 50 ku = 1,nu
         i      = nu*(kv-1)+ku
         x(i)   =   coh(kv) * r(i)
         y(i)   =   sih(kv) * r(i)
         xu(i)  =   coh(kv) * ru(i)
         yu(i)  =   sih(kv) * ru(i)
         xv(i)  =   coh(kv) * rv(i) - alp*y(i)
         yv(i)  =   sih(kv) * rv(i) + alp*x(i)
         yuu(i) =   sih(kv) * ruu(i)
         yuv(i) =   sih(kv) * ruv(i) + alp*xu(i)
         yvv(i) =   sih(kv) * rvv(i) + alp*(2.*xv(i)+alp*y(i))
         xuu(i) =   coh(kv) * ruu(i)
         xuv(i) =   coh(kv) * ruv(i) - alp*yu(i)
         xvv(i) =   coh(kv) * rvv(i) - alp*(2.*yv(i)-alp*x(i))
   50 continue
      do 60  i = 1 , nuv
         snx(i) = yu(i)*zv(i)-zu(i)*yv(i)
         sny(i) = zu(i)*xv(i)-xu(i)*zv(i)
         snz(i) = xu(i)*yv(i)-yu(i)*xv(i)
         xcur(i) = cut*xv(i)-cup*xu(i)
         ycur(i) = cut*yv(i)-cup*yu(i)
         zcur(i) = cut*zv(i)-cup*zu(i)
   60 continue

c.....Write coilsurface info to output if asked to do so
      if(w_csurf>0 .or. w_csurf==-1) then
         write (inesc, 320)
         write (inesc, 330) (i,i=0,ms)
         write (inesc, 340)
         do k = -ns, -1
            write (inesc, 360) k, (cr(i,0),i=0,ms)
         end do
         write (inesc, 360) k, (2*cr(i,k),i=0,ms)
         do k = 1, ns
            write (inesc, 360) k, (cr(i,k),i=0,ms)
         end do
         write (inesc, 370)
         write (inesc, 330) (i,i=0,ms)
         write (inesc, 340)
         do k = -ns, -1
            write (inesc, 360) k, (cz(i,0),i=0,ms)
         end do
         write (inesc, 360) k, (2*cz(i,k),i=0,ms)
         do k = 1, ns
            write (inesc, 360) k, (cz(i,k),i=0,ms)
         end do
 320     format(/1x,'     Coil Surface cr(m,n)'/)
 330     format(1x,'      m ',50(i5,:,5x))
 340     format(1x,'    n   ')
 360     format(4x,i2,4x,50(f9.6,:,1x))
 370     format(/1x,'     Coil Surface cz(m,n)'/)
      endif

      if(w_csurf>1 .or. w_csurf==-2) then
         write (inesc, '(a,i4,a)') '---- coil x,y,z,r(i) for i = 1,',
     1       nuvh,' ----'
         do i = 1, nuvh
            write (inesc,"(4g16.6)") x(i),y(i),z(i),r(i)
         enddo
         write (inesc, '(a)') '---- end coil x,y,z,r ----'
      endif

      if(w_csurf>2 .or. w_csurf==-3) then
         write (inesc, '(a,i4,a)')
     1      '---- coil snx,sny,snz(i) for i = 1,',nuvh,' ----'
         do i = 1, nuvh
            write (inesc,"(3g16.6)") snx(i),sny(i),snz(i)
         enddo
         write (inesc, '(a)') '---- end coil snx,sny,snz ----'
      endif

      if(w_csurf>3 .or. w_csurf==-4) then
         write (inesc, '(a,i4,a)')
     1     '---- coil xcur,ycur,zcur(i) for i = 1,',nuvh,' ----'
         do i = 1, nuvh
            write (inesc,"(3g16.6)") xcur(i),ycur(i),zcur(i)
         enddo
         write (inesc, '(a)') '---- end coil xcur,ycur,zcur ----'
      endif

      end subroutine surface_coil
