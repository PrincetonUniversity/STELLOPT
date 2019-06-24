subroutine quasisymmetry_residual

  use quasisymmetry_variables

  implicit none

  integer :: j

  residual = matmul(d_d_zeta, sigma) ! Could use BLAS for speed if needed

  residual = residual &
       + (iota + axis_helicity*nfp) * ( B1Squared_over_curvatureSquared * B1Squared_over_curvatureSquared + 1 + sigma * sigma) &
       - 2 * B1Squared_over_curvatureSquared * (-sign_psi * torsion + I2_over_B0) / B0_over_abs_G0 * sign_G

end subroutine quasisymmetry_residual
