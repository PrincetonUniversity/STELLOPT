
! ----------------------------------------------------------------------
      subroutine regint(f,a,b,c,n)
! ----------------------------------------------------------------------
      use stel_kinds
      implicit none
      integer :: i, n
      real(rprec), dimension(n) :: f, a, b, c
      real(rprec) :: sqp, sqm, sqa, sqc, top, tom
c ----------------------------------------------------------------------
      do  i=1,n
         sqp    = sqrt(a(i)+2._dp*b(i)+c(i))
         sqm    = sqrt(a(i)-2._dp*b(i)+c(i))
         sqa    = sqrt(a(i))
         sqc    = sqrt(c(i))
         top    = log((sqc*sqp+c(i)+b(i))/(sqa*sqp-a(i)-b(i)))/sqp
         tom    = log((sqc*sqm+c(i)-b(i))/(sqa*sqm-a(i)+b(i)))/sqm
         f(i)   = top + tom
      enddo

      end subroutine regint
