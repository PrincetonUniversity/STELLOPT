
!----------------------------------------------------------------------
      subroutine bn_precal
! ----------------------------------------------------------------------
c                                                          11.09.99
c     purpose:
c
c
c ----------------------------------------------------------------------
      use meshes
      use bnvariables
      implicit none
      integer :: i, m, n, ku, kv
c ----------------------------------------------------------------------
         pi     = 2*asin(1._dp)
         pi2    = 2 * pi
         alp    = pi2 / (np)
         alu    = pi2 / (nu)
         alv    = pi2 / (nv)
         alvp   = pi2 / (nvp)
         fnuv   = 1._dp/(nuv)
c
      do  m=0,md
         i      = 0
       do  kv = 1,nv
        do  ku = 1,nu
         i      = i + 1
         conu(i,m) = cos(alu*m*(ku-1))
         sinu(i,m) = sin(alu*m*(ku-1))
        enddo
       enddo
      enddo
      do  n=-nd,nd
         i      = 0
       do  kv = 1,nv
        do  ku = 1,nu
         i      = i + 1
         sinv(i,n) = sin(alv*n*(kv-1))
         conv(i,n) = cos(alv*n*(kv-1))
        enddo

       enddo
       do  i=1,nuv
         sinv(i,-n) = -sinv(i,n)
         conv(i,-n) =  conv(i,n)
       enddo
      enddo
      do kv=1,nvp
       do ku=1,nu
         i    = ku + nu*(kv-1)
         indv(i)  = kv
         indu(i)  = ku
       enddo
      enddo
      do  ku = -nu+1,nu-1
       if(iabs(ku).ne.nu/2) then
         tanu(ku) = tan(.5_dp*alu*ku)/pi
       else
         tanu(ku)=1.e+20_dp
       endif
      enddo
      do  kv = -nvp+1,nvp-1
       if(iabs(kv).ne.nvp/2) then
         tanv(kv) = tan(.5_dp*alvp*kv)/pi
        else
         tanv(kv) = 1.e+20_dp
        endif
      enddo

      end subroutine bn_precal
