
      subroutine smooth1(ya, n1, n2, wk, frac)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      use parambs
      use vmec0
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer n1, n2
      real(rprec) frac
      real(rprec), dimension(*) :: ya, wk
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
       real(rprec), parameter :: zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: n11, n21, n, i
      real(rprec) :: as, a, a1, a2, a3
C-----------------------------------------------
c--
c       smoothing real array ya(*) for i: n1.le.i.le.n2
c     frac - defines how strong is smoothing;
c     if frac=0 then only smoothes the peaks
c--
      n11 = n1 + 1
      n21 = n2 - 1
      n = n21 - n11
      if (n .le. 0) return                    !too little number of points
c
c- check the edges
c
      as = sum(abs(ya(n11:n21-1)-ya(n11+1:n21)))
      as = as/n
c
      if (as .le. zero) then
         ya(n1) = ya(n11)
         ya(n2) = ya(n21)
         return
      else
         a = 3*(as + abs(ya(n11)-ya(n11+1)))
         a1 = abs(ya(n1)-ya(n11))
         if (a1 > a) ya(n1) = 2*ya(n11) - ya(n11+1)
         a = 3*(as + abs(ya(n21)-ya(n21-1)))
         a1 = abs(ya(n2)-ya(n21))
         if (a1 > a) ya(n2) = 2*ya(n21) - ya(n21-1)
      endif
c
c- work array
c
      wk(n1:n2) = ya(n1:n2)
c
c- check for strong peaks and remove
c
      do i = n11, n21
         a1 = .5_dp*(ya(i+1)+ya(i-1))
         a2 = 3*abs(ya(i+1)-ya(i-1))
         a3 = abs(ya(i)-a1)
         if (a3 > a2) ya(i) = a1
      end do
c
c- smoothing with a factor frac
c
      if (frac .le. zero) return
      ya(n11:n21) =(ya(n11:n21)+.5_dp*(wk(n11+1:n21+1)+wk(n11-1:n21-1))*
     1   frac)/(one + frac)

      end subroutine smooth1
