
! ----------------------------------------------------------------------
      subroutine surface_plas
!................................................................
c     purpose:                                                01/01/89
c     Set up plasma surface from fourier coeffs
c     Calculate dsur, normals etc
c
c................................................................
      use Vmeshes
      use NumParams
      use OutCtrl, ONLY: w_psurf
      use LoopCtrl
      USE Vvacuum1
      USE Vvacuum2, ONLY: cr1, cz1, ms1, ns1
      USE Vvacuum3
      USE Vprecal1
      USE Vprecal2
c   following line removed: these are all from Vmeshes, some compilers can't cope
c     USE Vprecal3, ONLY: nuvh1, nuv1, nu1, nv1
      USE Vprecal4
      USE Vprecal5
c     USE Vprecal6, ONLY: nu       !  same here
      USE Vprecal7
      USE Vprecal8, ONLY: coh1, sih1
      USE Vprecal9, ONLY: comu1, simu1
      USE Vprecal10, ONLY: conv1, sinv1
      USE Vsurface9, ONLY: x1, y1, z1, r1
      USE Vsurfac10, ONLY: x1u, y1u, z1u
      USE Vsurfac11, ONLY: x1v, y1v, z1v
      USE Vsurfac12, ONLY: r1uu, z1uu, r1u
      USE Vsurfac13, ONLY: snx1, sny1, snz1, dsur1, fla
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, m, n, kv, ku, k
      real(rprec) :: cm, cn, cofp, cofm, sifp, sifm, u, v
c ----------------------------------------------------------------------
c   plasma surface
c
      i=0
      if (.not.allocated(x1))
     1 allocate (x1(nuvh1), y1(nuvh1), z1(nuvh1), r1(nuvh1),
     1          x1u(nuvh1), y1u(nuvh1), z1u(nuvh1), x1v(nuvh1),
     2          y1v(nuvh1), z1v(nuvh1), r1uu(nuvh1), z1uu(nuvh1),
     3          r1u(nuvh1), snx1(nuvh1), sny1(nuvh1), snz1(nuvh1),
     4          dsur1(nuvh1), stat = i)
      if (i .ne. 0) stop 'allocation error in NESCOIL surface_plas'

         r1(:)   = zero
         z1(:)   = zero
         r1u(:)  = zero
         y1v(:)  = zero
         z1u(:)  = zero
         z1v(:)  = zero
         r1uu(:)  = zero
         z1uu(:)  = zero

      do 20 m = 0,ms1
      do 20 n = 0,ns1
         cm     = m*pi2
         cn     = n*pi2
         i      = 0
      do 20 kv = 1 , nv1/2 + 1
      do 20 ku = 1 , nu1
         i      = i + 1
         cofp   = comu1(ku,m)*conv1(kv,n)-simu1(ku,m)*sinv1(kv,n)
         cofm   = comu1(ku,m)*conv1(kv,n)+simu1(ku,m)*sinv1(kv,n)
         sifp   = simu1(ku,m)*conv1(kv,n)+comu1(ku,m)*sinv1(kv,n)
         sifm   = simu1(ku,m)*conv1(kv,n)-comu1(ku,m)*sinv1(kv,n)
         r1(i)   = r1(i)    +       cr1(m,n)*cofp + cr1(m,-n)*cofm
         z1(i)   = z1(i)    +       cz1(m,n)*sifp + cz1(m,-n)*sifm
         r1u(i)  = r1u(i)   - cm *( cr1(m,n)*sifp + cr1(m,-n)*sifm )
         y1v(i)  = y1v(i)   - cn *( cr1(m,n)*sifp - cr1(m,-n)*sifm )
         z1u(i)  = z1u(i)   + cm *( cz1(m,n)*cofp + cz1(m,-n)*cofm )
         z1v(i)  = z1v(i)   + cn *( cz1(m,n)*cofp - cz1(m,-n)*cofm )
         r1uu(i) = r1uu(i) -cm*cm*( cr1(m,n)*cofp + cr1(m,-n)*cofm )
         z1uu(i) = z1uu(i) -cm*cm*( cz1(m,n)*sifp + cz1(m,-n)*sifm )
   20 continue
c
      do kv = 1, nv1/2 + 1
         do ku = 1, nu1
            i = nu1*(kv - 1) + ku
            x1(i) = coh1(kv)*r1(i)
            y1(i) = sih1(kv)*r1(i)
            x1u(i) = coh1(kv)*r1u(i)
            y1u(i) = sih1(kv)*r1u(i)
            x1v(i) = coh1(kv)*y1v(i) - alp*y1(i)
            y1v(i) = sih1(kv)*y1v(i) + alp*x1(i)
            snx1(i) = y1u(i)*z1v(i) - z1u(i)*y1v(i)
            sny1(i) = z1u(i)*x1v(i) - x1u(i)*z1v(i)
            snz1(i) = x1u(i)*y1v(i) - y1u(i)*x1v(i)
            dsur1(i) = sqrt(snx1(i)*snx1(i)+sny1(i)*sny1(i)+
     1         snz1(i)*snz1(i))
            snx1(i) = snx1(i)/dsur1(i)
            sny1(i) = sny1(i)/dsur1(i)
            snz1(i) = snz1(i)/dsur1(i)
         end do
      end do

c     Compute total plasma surface area
      fla = sum(dsur1(1:nu1)) + 2.*sum(dsur1(nu1+1:nuvh1-nu1)) +
     1   sum(dsur1(nuvh1-nu1+1:nuvh1))
      write (inesc, '(a,1pe10.3)') 'Total plasma surface area = ', fla

c.....Write plasma surface info to output if asked to do so
      if(w_psurf>0 .or. w_psurf==-1) then
         write (inesc, 320)
         write (inesc, 330) (i,i=0,ms1)
         write (inesc, 340)
         do k = -ns1, -1
            write (inesc, 360) k, (cr1(i,0),i=0,ms1)
         end do
         write (inesc, 360) k, (2*cr1(i,k),i=0,ms1)
         do k = 1, ns1
            write (inesc, 360) k, (cr1(i,k),i=0,ms1)
         end do
         write (inesc, 370)
         write (inesc, 330) (i,i=0,ms1)
         write (inesc, 340)
         do k = -ns1, -1
            write (inesc, 360) k, (cz1(i,0),i=0,ms1)
         end do
         write (inesc, 360) k, (2*cz1(i,k),i=0,ms1)
         do k = 1, ns1
            write (inesc, 360) k, (cz1(i,k),i=0,ms1)
         end do
 320     format(/1x,'     Plasma Surface cr1(m,n)'/)
 330     format(1x,'      m ',50(i5,:,5x))
 340     format(1x,'    n   ')
 360     format(4x,i2,4x,50(f9.6,:,1x))
 370     format(/1x,'     Plasma Surface cz1(m,n)'/)
      endif

      if(w_psurf>1 .or. w_psurf==-2) then
         write (inesc, '(a,i4,a)') '---- plasma x,y,z,r(i) for i = 1,',
     1       nuvh1,' ----'
         do i = 1, nuvh1
            write (inesc,"(4g16.6)") x1(i),y1(i),z1(i),r1(i)
         enddo
         write (inesc, '(a)') '---- end plasma x,y,z,r(i) ----'
      endif

      if(w_psurf>2 .or. w_psurf==-3) then
         write (inesc, '(a,i4,a)')
     1       '---- plasma dsur,snx,sny,snz for i = 1,',nuvh1,' ----'
         do i = 1, nuvh1
            write (inesc,"(4g16.6)") dsur1(i),snx1(i),sny1(i),snz1(i)
         enddo
         write (inesc,'(a)') '---- end plasma dsur,snx,sny,snz ----'
      endif

      if(w_psurf>3 .or. w_psurf==-4) then
         write (inesc, '(a,i4,a)')
     1     '---- plasma xu,yu,xv,yv(i) for i = 1,',nuvh1,' ----'
         do i = 1, nuvh1
            write (inesc,"(4g16.6)") x1u(i),y1u(i),x1v(i),y1v(i)
         enddo
         write (inesc, '(a)') '---- end plasma xu,yu,xv,yv(i) ----'
      endif

c     Compute bn_ext (external normal field)
      call bnfld

      end subroutine surface_plas
