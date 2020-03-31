
!----------------------------------------------------------------------
      subroutine precal
!................................................................
c                                                             01/01/89
c     purpose:
c     Calculate everything that does not need pot.coeffs, e.g.
c     surface and fft stuff
c................................................................
      use Vmeshes
      use NumParams
      use LoopCtrl
      USE Vprecal1
      USE Vprecal2
      USE Vprecal3
      USE Vprecal4
      USE Vprecal5
      USE Vprecal6, ONLY: comu, simu
      USE Vprecal7, ONLY: conv, sinv
      USE Vprecal8, ONLY: coh1, sih1
      USE Vprecal9, ONLY: comu1, simu1
      USE Vprecal10, ONLY: conv1, sinv1
      USE Vfourier2
      implicit none
c................................................................
C   L o c a l   V a r i a b l e s
c................................................................
      integer :: k, m, ku, n, kv, i
c................................................................
c
      i=0
      if (.not.allocated(trigs))
     1 allocate (trigs(2*nu), trigsv(2*nv), conv1(nv1,0:nd),
     1   sinv1(nv1,0:nd), comu(nu,0:md), simu(nu,0:md), factor(nuvh1),
     2   coh(nv), sih(nv), conv(nv,0:nd), sinv(nv,0:nd), coh1(nv1),
     3   sih1(nv1), comu1(nu1,0:md), simu1(nu1,0:md), eps(0:md + nd),
     4   stat = i)
      if (i .ne. 0) stop 'allocation error in NESCOIL precal'

      alp  = pi2/np
      alu  = pi2/nu
      alv  = pi2/nv
      alvp = pi2/(nv*np)
      fnuv = one/nuv
      fnv  = one/nv
      alu1 = pi2/nu1
      alv1 = pi2/nv1
      alvp1 = pi2/(nv1*np)
      fnuv1 = one/nuv1
      fnv1  = one/nv1
      call fftfax_g (nu, ifax, trigs)
      call cftfax_g (nv, ifaxv, trigsv)
      do 10 k=0,np-1
         cok(k) = cos(alp* k)
         sik(k) = sin(alp* k)
   10 continue
      do 20 m=0,md
      do 20 ku=1,nu
         comu(ku,m)   = cos(m*alu*(ku-1))
         simu(ku,m)   = sin(m*alu*(ku-1))
   20 continue
      do 30 n=0,nd
      do 30 kv=1,nv
         conv(kv,n)   = cos(n*alv*(kv-1))
         sinv(kv,n)   = sin(n*alv*(kv-1))
   30 continue
      do 40 kv=1,nv
         coh(kv)= cos(alvp*(kv-1))
         sih(kv)= sin(alvp*(kv-1))
   40 continue
      do 50 m=0,md
      do 50 ku=1,nu1
         comu1(ku,m)   = cos(m*alu1*(ku-1))
         simu1(ku,m)   = sin(m*alu1*(ku-1))
   50 continue
      do 60 n=0,nd
      do 60 kv=1,nv1
         conv1(kv,n)   = cos(n*alv1*(kv-1))
         sinv1(kv,n)   = sin(n*alv1*(kv-1))
   60 continue
      do 70 kv=1,nv1
         coh1(kv)= cos(alvp1*(kv-1))
         sih1(kv)= sin(alvp1*(kv-1))
   70 continue
      do 80 i=1,nu1
         factor(i) = one
         factor(i+nuv1/2) = one
   80 continue
      do 90 i=nu1+1,nuv1/2
         factor(i) = two
   90 continue
         eps(0) = one/two
      do 100 i=1,md+nd
         eps(i) = one
  100 continue

      end subroutine precal
