
      subroutine bn_fouri(bnfou,bnfou_c)
c ----------------------------------------------------------------------
c     purpose:                                               06.04.00
c
c
c ---------------------------------------------------------------------
      use meshes
      use bnvariables
      implicit none
      integer :: i, m, n, nlo
c ---------------------------------------------------------------------
      real(rprec), dimension(0:mf,-nf:nf)   :: bnfou,bnfou_c
      real(rprec), dimension(0:md,-nd:nd)   :: aufou,avfou
      real(rprec), dimension(0:md,-nd:nd)   :: aufou_s,avfou_s
      real(rprec), dimension(nuv)           :: bn
      real(rprec)                           :: faz,cofp,sifp
c ----------------------------------------------------------------------
      do  m  = 0,mf
         faz   = 2*fnuv
         if(m.eq.0) faz = fnuv
         do  n  = -nf,nf
            aufou(m,n) = 0
            avfou(m,n) = 0
            aufou_s(m,n) = 0
            avfou_s(m,n) = 0
            do  i = 1,nuv
               cofp   = conu(i,m)*conv(i,n)-sinu(i,m)*sinv(i,n)         !cos(u+v)
               sifp   = sinu(i,m)*conv(i,n)+conu(i,m)*sinv(i,n)         !sin(u+v)
               aufou(m,n) = aufou(m,n) +au(i)*cofp*faz
               avfou(m,n) = avfou(m,n) +av(i)*cofp*faz
               !aufou_s(m,n) = aufou_s(m,n) +au(i)*sifp*fnuv
               !avfou_s(m,n) = avfou_s(m,n) +av(i)*sifp*fnuv
            enddo
         enddo
      enddo
      bn(:nuv) = 0
      do  m  = 0,mf
         do  n  = -nf,nf
            do  i = 1,nuv
               cofp   = conu(i,m)*conv(i,n)-sinu(i,m)*sinv(i,n)         !cos(u+v)
               sifp   = sinu(i,m)*conv(i,n)+conu(i,m)*sinv(i,n)         !sin(u+v)
               bn(i) = bn(i) +pi2*(m*avfou(m,n)-n*aufou(m,n))*sifp
!     1                       +pi2*(m*avfou_s(m,n)+n*aufou_s(m,n))*cofp
            enddo
         enddo
      enddo
      do  m  = 0,mf
         faz   = 2*fnuv
         if(m.eq.0) faz = fnuv
         do  n  = -nf,nf
            bnfou(m,n) = 0
            bnfou_c(m,n) = 0
            do  i = 1,nuv
               cofp   = conu(i,m)*conv(i,n)-sinu(i,m)*sinv(i,n)         !cos(u+v)
               sifp   = sinu(i,m)*conv(i,n)+conu(i,m)*sinv(i,n)         !sin(u+v)
               bnfou(m,n) = bnfou(m,n) +bn(i)/sqf(i)*sifp*faz
!               bnfou_c(m,n)=bnfou_c(m,n) +bn(i)/sqf(i)*cofp*faz
            enddo
         enddo
      enddo
!     end  subroutine bn_fouri
      end
