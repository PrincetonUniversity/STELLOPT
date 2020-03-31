
! -------------------------------------------------------------
      subroutine surfcur_diag
!................................................................
c     Purpose:
c................................................................

      USE Vmeshes
      USE NumParams
      use OutCtrl, ONLY: w_jsurf
      USE LoopCtrl
      USE Vvacuum3
      USE Vprecal1
      USE Vprecal2
      USE Vprecal3
      USE Vprecal6
      USE Vprecal7
      USE Vsurface1
      USE Vsurface2, ONLY: xu, yu, xv
      USE Vsurface3, ONLY: yv, zu, zv
      USE Vsurface4, ONLY: xuu, yuu, zuu
      USE Vsurface5, ONLY: xuv, yuv, zuv
      USE Vsurface6, ONLY: xvv, yvv, zvv
      USE Vsurface7, ONLY: snx, sny, snz
      USE Vdiagno2
      USE Vdiagno3
      use Vsurfcurdiag
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, m, n, ku, kv
      real(rprec), dimension(:), allocatable ::  wx, wy, wz,
     1     curvat, cude, potu, potv, potuu, potuv, potvv, vx, vy, vz
      real(rprec) :: cw,cwm,sw,swm,cm,cn,v2,w2,vw
c ----------------------------------------------------------------------
      i=0
      if (.not.allocated(wx))
     1 allocate (wx(nuv), wy(nuv), wz(nuv), curvat(nuv), cude(nuv),
     1   potuu(nuv), potuv(nuv), potvv(nuv), vx(nuv), vy(nuv),
     2   vz(nuv), potu(nuv), potv(nuv), pote(nuv), stat = i)
      if (i .ne. 0) return
c
      vx= zero
      vy= zero
      vz= zero
      curvat= zero
      cude= zero
      wx= zero
      wy= zero
      wz= zero

      pote = zero
      potu = cut
      potv = cup
      potuu = zero
      potuv = zero
      potvv = zero

      do 30 n= 0,nf
         cn      = pi2*n
      do 30 m= 0,mf
         cm      = pi2*m
         i       = 0
      do 30  kv = 1,nv
      do 30  ku = 1,nu
         i      = i+1
         cw     = comu(ku,m)*conv(kv,n)-simu(ku,m)*sinv(kv,n)
         cwm    = comu(ku,m)*conv(kv,n)+simu(ku,m)*sinv(kv,n)
         sw     = simu(ku,m)*conv(kv,n)+comu(ku,m)*sinv(kv,n)
         swm    = simu(ku,m)*conv(kv,n)-comu(ku,m)*sinv(kv,n)
         pote(i)  = pote(i)+ eps(m)*eps(n)*(pot(m,n)*sw + pot(m,-n)*swm)
         potu(i)  = potu(i) + cm   *eps(n)*(pot(m,n)*cw + pot(m,-n)*cwm)
         potv(i)  = potv(i) + cn   *eps(m)*(pot(m,n)*cw - pot(m,-n)*cwm)
         potuu(i) = potuu(i)- cm*cm*eps(n)*(pot(m,n)*sw + pot(m,-n)*swm)
         potuv(i) = potuv(i)- cm      *cn* (pot(m,n)*sw - pot(m,-n)*swm)
         potvv(i) = potvv(i)- cn*cn*eps(m)*(pot(m,n)*sw + pot(m,-n)*swm)
   30 continue

c     curvature  of the surface current lines
c
      do i=1,nuv
         vx(i)  = xu(i)*potv(i)-xv(i)*potu(i)
         vy(i)  = yu(i)*potv(i)-yv(i)*potu(i)
         vz(i)  = zu(i)*potv(i)-zv(i)*potu(i)
         wx(i)  = xu(i)  *(potuv(i)*potv(i)-potvv(i)*potu(i))
     $          + xv(i)  *(potuv(i)*potu(i)-potuu(i)*potv(i))
     $          - potu(i)*(  xuv(i)*potv(i)-  xvv(i)*potu(i))
     $          - potv(i)*(  xuv(i)*potu(i)-  xuu(i)*potv(i))
c
         wy(i)  = yu(i)  *(potuv(i)*potv(i)-potvv(i)*potu(i))
     $          + yv(i)  *(potuv(i)*potu(i)-potuu(i)*potv(i))
     $          - potu(i)*(  yuv(i)*potv(i)-  yvv(i)*potu(i))
     $          - potv(i)*(  yuv(i)*potu(i)-  yuu(i)*potv(i))
c
         wz(i)  = zu(i)  *(potuv(i)*potv(i)-potvv(i)*potu(i))
     $          + zv(i)  *(potuv(i)*potu(i)-potuu(i)*potv(i))
     $          - potu(i)*(  zuv(i)*potv(i)-  zvv(i)*potu(i))
     $          - potv(i)*(  zuv(i)*potu(i)-  zuu(i)*potv(i))
      end do

      do i=1,nuvh
         v2     = vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i)
         w2     = wx(i)*wx(i)+wy(i)*wy(i)+wz(i)*wz(i)
         vw     = vx(i)*wx(i)+vy(i)*wy(i)+vz(i)*wz(i)
         curvat(i) = sqrt(w2-vw*vw/v2)/v2
         cude(i) = sqrt(v2/(snx(i)*snx(i)+sny(i)*sny(i)+snz(i)*snz(i)))
      end do

      curvra = maxval(curvat(1:nuvh))
      if (curvra .ne. zero) curvra = 1./curvra
      curvri = minval(curvat(1:nuvh))
      if (curvri .ne. zero ) curvri = 1./curvri
      cudema = maxval(cude(1:nuvh))
      cudemi = minval(cude(1:nuvh))
      avcude = (sum(cude(1:nu)) + 2.*sum(cude(1+nu:nuvh-nu))
     1       +  sum(cude(nuvh-nu+1:nuvh)))/float(nuv)

c.....Always Write this surface current info to output
      write (inesc, 20) 'J Surf max, min, ave = ',cudema, cudemi, avcude
      write (inesc, 20) 'Curv_R of J min, max = ',curvra, curvri

c.....Write extra surface current info to output if asked to do so
      if(w_jsurf>0 .or. w_jsurf==-1) then
         write (inesc, '(a,i4,a)') '---- potential phi(i) for i = 1,',
     1        nuv,' ----'
         write (inesc,*) (pote(i), i = 1, nuv)
         write (inesc, '(a)') '---- end  potential phi(i) ----'
      endif

      if(w_jsurf>1 .or. w_jsurf==-2) then
         write (inesc, '(a,i4,a)')
     1      '---- current jsurf(i) for i = 1,',nuvh,' ----'
         write (inesc,*) (cude(i), i = 1, nuvh)
         write (inesc, '(a)') '---- end  current jsurf(i) ----'
      endif

      if (iloop .le. 0) deallocate (wx, wy, wz, curvat, cude,
     1   potuu, potuv, potvv, vx, vy, vz, potu, potv)

 20   format(a,1p3e16.8)

      end subroutine surfcur_diag
