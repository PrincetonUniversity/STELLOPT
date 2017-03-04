

! ----------------------------------------------------------------------
      subroutine bn_surface
! ----------------------------------------------------------------------
c                                                            11.08.99
c     purpose:
c
c
c ----------------------------------------------------------------------
      use meshes
      use bnvariables
      implicit none
c ----------------------------------------------------------------------
      integer :: i, m, n, ku, kv, k, np2
      real(rprec), dimension(:), allocatable :: r, ru, rv
      real(rprec) :: snx, sny, snz, coh, sih, co, si, cofp, sifp,
     1    cofm, sifm, cm, cn
c ----------------------------------------------------------------------
      allocate (r(nuv), ru(nuv), rv(nuv), stat=i)

      do  i = 1 , nuv
         r(i)    = 0._dp
         z(i)    = 0._dp
         ru(i)   = 0._dp
         rv(i)   = 0._dp
         zu(i)   = 0._dp
         zv(i)   = 0._dp
      enddo
      do  m = 0,mb
         do  n = -nb,nb
            cm     = m*pi2
            cn     = n*pi2
            do  i  = 1 , nuv                                             
               cofp   = conu(i,m)*conv(i,n)-sinu(i,m)*sinv(i,n)         !cos(u+v)
               cofm   = conu(i,m)*conv(i,n)+sinu(i,m)*sinv(i,n)         !cos(u-v)
               sifp   = sinu(i,m)*conv(i,n)+conu(i,m)*sinv(i,n)         !sin(u+v)
               sifm   = sinu(i,m)*conv(i,n)-conu(i,m)*sinv(i,n)         !sin(u-v)
               r(i)   = r(i)    +       cr(m,n)*cofp + crs(m,n)*sifp
               ru(i)  = ru(i)   - cm *( cr(m,n)*sifp - crs(m,n)*cofp)
               rv(i)  = rv(i)   - cn *( cr(m,n)*sifp - crs(m,n)*cofp)
               z(i)   = z(i)    +       cz(m,n)*sifp + czc(m,n)*cofp
               zu(i)  = zu(i)   + cm *( cz(m,n)*cofp - czc(m,n)*sifp)
               zv(i)  = zv(i)   + cn *( cz(m,n)*cofp - czc(m,n)*sifp)
            enddo
         enddo
      enddo
c----------------------------------------------------------
      do  kv = 1, nv
         coh    = cos(alvp*(kv-1))
         sih    = sin(alvp*(kv-1))
         do  ku = 1, nu
            i      = nu*(kv-1)+ku
            x(i)   = coh * r(i)
            y(i)   = sih * r(i)
            xu(i)  = coh * ru(i)
            yu(i)  = sih * ru(i)
            xv(i)  = coh * rv(i) - alp*y(i)
            yv(i)  = sih * rv(i) + alp*x(i)
         enddo
      enddo
      
      do   i = 1 , nuv
         snx    = yu(i)*zv(i)-zu(i)*yv(i)
         sny    = zu(i)*xv(i)-xu(i)*zv(i)
         snz    = xu(i)*yv(i)-yu(i)*xv(i)
         sqf(i) = sqrt(snx*snx+sny*sny+snz*snz)
         guu(i) = xu(i)*xu(i)+yu(i)*yu(i)+zu(i)*zu(i)
         guv(i) = xu(i)*xv(i)+yu(i)*yv(i)+zu(i)*zv(i)
         gvv(i) = xv(i)*xv(i)+yv(i)*yv(i)+zv(i)*zv(i)
         dju(i) = bu(i)*guv(i)- bv(i)*guu(i)
         djv(i) = bu(i)*gvv(i)- bv(i)*guv(i)
      enddo
      do i=1,nuv
         djx(i)  = bu(i)*xv(i)-bv(i)*xu(i)
         djy(i)  = bu(i)*yv(i)-bv(i)*yu(i)
         djz(i)  = bu(i)*zv(i)-bv(i)*zu(i)
      enddo

      do k=1,np-1
        co  =cos(alp*k)
        si  =sin(alp*k)
       do i=1,nuv
           x(i+k*nuv) = x(i)*co - y(i)*si
           y(i+k*nuv) = y(i)*co + x(i)*si
           z(i+k*nuv) = z(i)
         djx(i+k*nuv) = djx(i)*co-djy(i)*si
         djy(i+k*nuv) = djy(i)*co+djx(i)*si
         djz(i+k*nuv) = djz(i)
       enddo
      enddo
      np2 =  np*np
      do i=1,nuv
         guv(i) = np*guv(i)
         gvv(i) = np2*gvv(i)
      enddo

      deallocate (r, ru, rv, stat=i)

      end subroutine bn_surface
