      MODULE legendre_params                                                 !! LEGENDRE
      USE stel_kinds                                                          !! LEGENDRE
      INTEGER :: n_leg                                                       !! LEGENDRE
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: a_leg, b_leg,              !! LEGENDRE
     1   a_leg_inv, b_leg_inv                                                !! LEGENDRE
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: tc, ti, tm                   !! LEGENDRE
      END MODULE legendre_params                                             !! LEGENDRE
