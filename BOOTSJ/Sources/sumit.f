
      function sumit (f, mbuse, nbuse)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      use stel_kinds
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer mbuse, nbuse
      real(rprec), dimension(-mbuse:mbuse,0:nbuse) :: f
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: n, m
      real(rprec) :: sumit
C-----------------------------------------------
c-
c specific sum of terms in f(i,j) table
c-
c
      sumit = 0
      do m = -mbuse, mbuse
         sumit = sumit + sum(f(m,1:nbuse))
      end do
      sumit = 2*sumit
      sumit = sumit + sum(f(-mbuse:mbuse,0))
      sumit = sumit - f(0,0)

      end function sumit
