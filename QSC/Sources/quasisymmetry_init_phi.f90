subroutine quasisymmetry_init_phi

  use quasisymmetry_variables

  implicit none

  real(dp), dimension(:,:), allocatable :: temp_matrix
  real(dp), dimension(:), allocatable :: temp_vector
  integer :: option, quadrature_option, i

  if (allocated(phi)) deallocate(phi)
  if (allocated(phi_extended)) deallocate(phi_extended)
  if (allocated(d_d_phi)) deallocate(d_d_phi)

  if (allocated(X1s)) deallocate(X1s)
  if (allocated(X1c)) deallocate(X1c)
  if (allocated(Y1s)) deallocate(Y1s)
  if (allocated(Y1c)) deallocate(Y1c)
  if (allocated(R1s)) deallocate(R1s)
  if (allocated(R1c)) deallocate(R1c)
  if (allocated(Z1s)) deallocate(Z1s)
  if (allocated(Z1c)) deallocate(Z1c)
  if (allocated(sigma)) deallocate(sigma)
  if (allocated(elongation)) deallocate(elongation)

  if (allocated(Jacobian)) deallocate(Jacobian)
  if (allocated(residual)) deallocate(residual)
  if (allocated(step_direction)) deallocate(step_direction)
  
  allocate(phi(N_phi))
  allocate(d_d_phi(N_phi,N_phi))
  allocate(temp_matrix(N_phi,N_phi))
  allocate(temp_vector(N_phi))

  allocate(X1s(N_phi))
  allocate(X1c(N_phi))
  allocate(Y1s(N_phi))
  allocate(Y1c(N_phi))
  allocate(R1s(N_phi))
  allocate(R1c(N_phi))
  allocate(Z1s(N_phi))
  allocate(Z1c(N_phi))
  allocate(sigma(N_phi))
  allocate(elongation(N_phi))

  matrix_size = N_phi
  allocate(Jacobian(matrix_size, matrix_size))
  allocate(residual(matrix_size))
  allocate(step_direction(matrix_size))

  option = 20
  quadrature_option = 0
  call quasisymmetry_differentiation_matrix(N_phi,0_dp, 2*pi/nfp, option, quadrature_option, phi, temp_vector, d_d_phi, temp_matrix)

  phi_extended = [( 2*pi*i/(N_phi*nfp), i=0,N_phi*nfp-1 )]

  dimension_Fourier = max(max_axis_nmax+1, (N_phi+1) / 2)
  if (allocated(sin_n_phi)) deallocate(sin_n_phi)
  if (allocated(cos_n_phi)) deallocate(cos_n_phi)
  allocate(sin_n_phi(N_phi,dimension_Fourier))
  allocate(cos_n_phi(N_phi,dimension_Fourier))
  sin_n_phi(:,1) = 0
  cos_n_phi(:,1) = 1
  do i = 2, dimension_Fourier
     sin_n_phi(:,i) = sin((i-1) * nfp * phi)
     cos_n_phi(:,i) = cos((i-1) * nfp * phi)
  end do

  d_phi = phi(2) - phi(1)

  deallocate(temp_matrix, temp_vector)

end subroutine quasisymmetry_init_phi
