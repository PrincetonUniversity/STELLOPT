subroutine quasisymmetry_Jacobian

  use quasisymmetry_variables

  implicit none

  integer :: j

  ! d (Riccati equation) / d sigma:
  ! For convenience we will fill all the columns now, and re-write the first column in a moment.
  Jacobian = d_d_zeta
  do j = 1,N_phi
     Jacobian(j,j) = Jacobian(j,j) + iota * 2 * sigma(j)
  end do

  ! d (Riccati equation) / d iota:
  Jacobian(:, 1) = B1Squared_over_curvatureSquared * B1Squared_over_curvatureSquared + 1 + sigma * sigma

end subroutine quasisymmetry_Jacobian
