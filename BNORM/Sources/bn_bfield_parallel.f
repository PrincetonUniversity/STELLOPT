
! ----------------------------------------------------------------------
      subroutine bn_bfield_parallel
! ----------------------------------------------------------------------
c                                                            11.09.99
c     purpose:
c
c
c ----------------------------------------------------------------------
      use meshes
      use bnvariables
      implicit none
      integer :: i, m, n
      real(rprec) :: cop, sip
c ----------------------------------------------------------------------
      curpol = alp*bsubv(0,0)
c     write(6,3000) curpol
 3000 format('   curpol= ',1pe12.4)
      do m=  0,md
         do n=-nd,nd
            bsubu(m,n) = pi2*bsubu(m,n)/curpol
            bsubv(m,n) = alp*bsubv(m,n)/curpol
            bsubus(m,n) = pi2*bsubus(m,n)/curpol
            bsubvs(m,n) = alp*bsubvs(m,n)/curpol
         enddo
      enddo

      do  i = 1 , nuv
         bu (i)   = 0._dp
         bv (i)   = 0._dp
      enddo

      do  m = 0,mb
         do  n = -nb,nb
            do  i  = 1,nuv
               cop      = conu(i,m)*conv(i,n)-sinu(i,m)*sinv(i,n)       !cos(u+v)
               sip      = sinu(i,m)*conv(i,n)+conu(i,m)*sinv(i,n)       !sin(u+v)
               bu (i)   = bu (i) + bsubu(m,n)*cop
     1                           + bsubus(m,n)*sip
               bv (i)   = bv (i) + bsubv(m,n)*cop
     1                           + bsubvs(m,n)*sip
            enddo
         enddo
      enddo

      end subroutine   bn_bfield_parallel
