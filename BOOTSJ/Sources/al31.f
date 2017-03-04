
      function al31 (x, zeff, alphae, alphai)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      use vmec0
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec) :: x, zeff, alphae, alphai
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec) :: z, z2, d, a, al31
C-----------------------------------------------
c-
c calculates L31 transport coefficient, according to
c S. Hirshman, Phys. Fluids 31, 3150 (1988)
c x - ratio of trapped to circulated particles
c-
      z = zeff                              !effective ion charge number
      z2 = z**2
      d=1.414_dp*z+z2+x*(0.754_dp+2.657_dp*z+2.0_dp*z2)+
     1   x*x*(0.348_dp+1.243_dp*z+z2)
      a = 0.754_dp + 2.21_dp*z + z2 + x*(0.348_dp + 1.243_dp*z + z2)
      al31 = x*a/d
      alphae = 1 - (0.884_dp + 2.074_dp*z)/a
      alphai = 1 - 1.172_dp/(1 + 0.462_dp*x)

      end function al31
