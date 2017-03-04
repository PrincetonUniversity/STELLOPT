      subroutine bn_vecpot
c ----------------------------------------------------------------------
c     purpose:                                               11.08.99   

c
c
c ---------------------------------------------------------------------
      use meshes
      use bnvariables
      implicit none
c ----------------------------------------------------------------------
      integer :: i, ip
      real(rprec), dimension(:), allocatable :: ax, ay, az, analyt
      real(rprec) :: pi41, dintu, dintv, pi21      
c ----------------------------------------------------------------------
      allocate (ax(nuv), ay(nuv), az(nuv), analyt(nuv),stat=i)

      pi41  = .5_dp/pi2
      pi21  = 1._dp/pi2
      ax   = 0;      ay   = 0;      az   = 0 

      do i = 1,nuv
         call vecpot (ax(i), ay(i), az(i), i, 1, i-1)
         call vecpot (ax(i), ay(i), az(i), i, i+1, nuvp)
      enddo
      
!      do i=1,nuv
!      print *,ax(i),ay(i),az(i)
!      enddo

!     FREE UP BIG ARRAYS
      deallocate (x, y, z, indu, indv, djx, djy, djz)
      
      call regint(analyt, guu, guv, gvv, nuv) ! compute  singular  integral                

      do  i  = 1,nuv
         dintu = (ax(i)*xu(i)+ay(i)*yu(i)+az(i)*zu(i))*fnuv
         dintv = (ax(i)*xv(i)+ay(i)*yv(i)+az(i)*zv(i))*fnuv
         au(i) = pi41*(dintu +dju(i)*analyt(i)*np)                      
         av(i) = pi41*(dintv +djv(i)*analyt(i)*np)                      
      enddo

      deallocate (ax, ay, az, analyt)

!     end  subroutine bn_vecpot
      end
