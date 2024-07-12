      function ddsum(n, din)
      use stel_kinds
      implicit none
      integer :: n
      real(rprec), intent(in) :: din(n)
      real(rprec) :: ddsum
!       
!     accurate summation by sorting and
!     adding in increasing order
!
      logical is_all_negative, is_all_positive
      integer istart
      real(rprec) :: possum,negsum
      real(rprec), parameter :: zero = 0


      integer i, info
      integer, allocatable :: iperm(:)
      real(rprec), allocatable :: dinabs(:)
      real(rprec) dsum

!SIMPLE FUNCTIONALITY      ddsum = SUM(din)
!      return

      allocate (iperm(n), dinabs(n))

!     ------------------------
!     use shell sort to sort increasing abs value
!     ------------------------
      dinabs = abs(din)
      call dshell(n, dinabs, iperm )

      dsum = zero
      do i=1,n
         dsum = dsum + din(iperm(i))
      enddo

      ddsum = dsum

      deallocate (iperm, dinabs)

      return
      end function ddsum


      subroutine dshell(n,dx,iperm)
      use stel_kinds
      integer, intent(in) :: n
      integer, intent(out) :: iperm(n)
      real(rprec), intent(in) :: dx(n)
      integer :: i
      real(rprec) :: dtemp,dxert

      integer j,k,m,maxmn
      integer itemp
      
      do i=1,n
        iperm(i) = i
      enddo

      m = n/2
      do while (m.gt.0)
         maxmn = n - m


         do 23004 j = 1,maxmn
            do 23006 k = j,1,-m
            dxert = dx(iperm(k+m))
            if (dxert .ge. dx(iperm(k))) exit

            itemp = iperm(k+m)
            iperm(k+m) = iperm(k)
            iperm(k) = itemp


c      dtemp = dx(k+m)
c      dx(k+m) = dx(k)
c      dx(k) = dtemp


23006       continue
C end do k

23004    continue
C end do j

      m = m/2

      end do
C end while

      end subroutine dshell

      subroutine dshell2(n,arrin,iperm)
      use stel_kinds
      implicit none
      integer, intent(in) :: n
      integer, intent(out) :: iperm(n)
      real(rprec), intent(in) :: arrin(n)
      integer :: i, j, ir, m, indx
      real(rprec) :: dtemp
      
      do i=1,n
        iperm(i) = i
      enddo

      m = n/2+1
      ir = n

      do while (m .gt. 0)
         if (m .gt. 1) then
            m = m-1
            indx = iperm(m)
            dtemp = arrin(indx)
         else
            indx = iperm(ir)
            dtemp = arrin(indx)
            iperm(ir) = iperm(1)
            ir = ir-1
            if (ir .le. 1) then
!            if (ir .eq. 1) then
               iperm(1) = indx
               exit
            end if
         end if
    
         i = m
         j = m+m
         do while (j .le. ir) 
            if (j .lt. ir) then
               if (arrin(iperm(j)) .lt. arrin(iperm(j+1))) j=j+1
            endif
            if (dtemp .lt. arrin(iperm(j))) then
               iperm(i) = iperm(j)
               i = j
               j = j+j
            else
               j = ir+1
            endif
         end do

         iperm(i) = indx

      end do

      end subroutine dshell2
