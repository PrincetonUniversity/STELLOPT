
      subroutine positiv(ya, n, ivar)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      use parambs
      use vmec0
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer n, ivar
      real(rprec), dimension(*) :: ya
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
       real(rprec), parameter :: zero = 0, D36 = 1.E-36_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i
C-----------------------------------------------
c--
c     make corrections to ensure a positive array
c--
c
      if (ivar .eq. 1) then
         ya(:n) = abs(ya(:n)) + D36
      else
         where (ya(:n) .le. zero) ya(:n) = D36
      endif

      end subroutine positiv
