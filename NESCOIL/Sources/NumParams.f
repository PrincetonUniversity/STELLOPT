      module NumParams
      use stel_kinds
      integer, parameter :: inescoil = 6
      integer :: inesc = inescoil
      real(rprec), parameter :: zero = 0,
     1    one = 1, two = 2, three = 3,
     2    pi = 3.14159265358979312_dp, pi2 = two*pi, mu0 = 4.0e-7_dp*pi

      contains
        logical function is_zero(x)
        use stel_kinds
        real(rprec) :: x
        is_zero = abs(x) < tiny(x)
        end function is_zero

      end module NumParams
