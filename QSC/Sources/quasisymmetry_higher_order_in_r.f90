subroutine quasisymmetry_higher_order_in_r

  use quasisymmetry_variables

  implicit none

  real(dp), dimension(:), allocatable :: V1, V2, V3, qs, qc, rs, rc
  real(dp), dimension(:), allocatable :: Y2s_from_X20, Y2s_inhomogeneous, Y2c_from_X20, Y2c_inhomogeneous
  real(dp), dimension(:), allocatable :: fX0_from_X20, fX0_from_Y20, fX0_inhomogeneous
  real(dp), dimension(:), allocatable :: fXs_from_X20, fXs_from_Y20, fXs_inhomogeneous
  real(dp), dimension(:), allocatable :: fXc_from_X20, fXc_from_Y20, fXc_inhomogeneous
  real(dp), dimension(:), allocatable :: fY0_from_X20, fY0_from_Y20, fY0_inhomogeneous
  real(dp), dimension(:), allocatable :: fYs_from_X20, fYs_from_Y20, fYs_inhomogeneous
  real(dp), dimension(:), allocatable :: fYc_from_X20, fYc_from_Y20, fYc_inhomogeneous
  real(dp) :: factor, iota_N, beta_1s, beta_2c, beta_2s, normalizer, max_eq1residual, max_eq2residual
  real(dp), dimension(:,:), allocatable :: matrix
  real(dp), dimension(:), allocatable :: right_hand_side
  integer :: j, iunit=20
  real(dp), dimension(:), allocatable :: fX0, fXs, fXc, fY0, fYs, fYc, eq1residual, eq2residual
  real(dp), dimension(:), allocatable :: d_X20_d_zeta, d_X2s_d_zeta, d_X2c_d_zeta
  real(dp), dimension(:), allocatable :: d_Y20_d_zeta, d_Y2s_d_zeta, d_Y2c_d_zeta
  real(dp), dimension(:), allocatable :: d_Z20_d_zeta, d_Z2s_d_zeta, d_Z2c_d_zeta, d_X3c3_d_zeta, d_X3s3_d_zeta
  real(dp), dimension(:), allocatable :: d_X3c1_d_zeta, d_X3s1_d_zeta, d_Y3c1_d_zeta, d_Y3s1_d_zeta, d_Y3c3_d_zeta, d_Y3s3_d_zeta
  real(dp), dimension(:), allocatable :: d_Z3c1_d_zeta, d_Z3s1_d_zeta, d_Z3c3_d_zeta, d_Z3s3_d_zeta
  real(dp), dimension(:), allocatable :: flux_constraint_coefficient
  real(dp), dimension(:), allocatable :: Q, predicted_flux_constraint_coefficient
  integer :: N_helicity
  integer :: vector_size, index_mixedPartialsEquation_0, index_mixedPartialsEquation_s, index_mixedPartialsEquation_c
  integer :: index_XYEquation_0, index_XYEquation_s, index_XYEquation_c, index_initialCondition
  integer :: index_X3c1, index_X3s1, index_Y3c1, index_Y3s1, index_Y3c3, index_Y3s3, index_iota2
  real(dp) :: I2, G0, B1c, Bbar, I4, G2
  ! Variables needed by LAPACK:
  integer :: INFO
  integer, dimension(:), allocatable :: IPIV

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  iota_N = iota + axis_helicity*nfp
  abs_G0_over_B0 = 1 / B0_over_abs_G0

  if ((verbose) .and. abs(iota_N) < 1e-8) print "(a,es21.14)","Warning: |iota_N| is very small so O(r^2) solve will be poorly conditioned. iota_N=",iota_N

  if (allocated(X20)) deallocate(X20)
  if (allocated(Y20)) deallocate(Y20)
  if (allocated(Z20)) deallocate(Z20)
  if (allocated(X2s)) deallocate(X2s)
  if (allocated(Y2s)) deallocate(Y2s)
  if (allocated(Z2s)) deallocate(Z2s)
  if (allocated(X2c)) deallocate(X2c)
  if (allocated(Y2c)) deallocate(Y2c)
  if (allocated(Z2c)) deallocate(Z2c)

  if (allocated(X20_untwisted)) deallocate(X20_untwisted)
  if (allocated(Y20_untwisted)) deallocate(Y20_untwisted)
  if (allocated(Z20_untwisted)) deallocate(Z20_untwisted)
  if (allocated(X2s_untwisted)) deallocate(X2s_untwisted)
  if (allocated(Y2s_untwisted)) deallocate(Y2s_untwisted)
  if (allocated(Z2s_untwisted)) deallocate(Z2s_untwisted)
  if (allocated(X2c_untwisted)) deallocate(X2c_untwisted)
  if (allocated(Y2c_untwisted)) deallocate(Y2c_untwisted)
  if (allocated(Z2c_untwisted)) deallocate(Z2c_untwisted)

  if (allocated(R20)) deallocate(R20)
  if (allocated(R2s)) deallocate(R2s)
  if (allocated(R2c)) deallocate(R2c)
  if (allocated(z20_cylindrical)) deallocate(z20_cylindrical)
  if (allocated(z2s_cylindrical)) deallocate(z2s_cylindrical)
  if (allocated(z2c_cylindrical)) deallocate(z2c_cylindrical)

  allocate(X20(N_phi))
  allocate(X2s(N_phi))
  allocate(X2c(N_phi))
  allocate(Y20(N_phi))
  allocate(Y2s(N_phi))
  allocate(Y2c(N_phi))
  allocate(Z20(N_phi))
  allocate(Z2s(N_phi))
  allocate(Z2c(N_phi))

  allocate(X20_untwisted(N_phi))
  allocate(X2s_untwisted(N_phi))
  allocate(X2c_untwisted(N_phi))
  allocate(Y20_untwisted(N_phi))
  allocate(Y2s_untwisted(N_phi))
  allocate(Y2c_untwisted(N_phi))
  allocate(Z20_untwisted(N_phi))
  allocate(Z2s_untwisted(N_phi))
  allocate(Z2c_untwisted(N_phi))

  allocate(R20(N_phi))
  allocate(R2s(N_phi))
  allocate(R2c(N_phi))
  allocate(z20_cylindrical(N_phi))
  allocate(z2s_cylindrical(N_phi))
  allocate(z2c_cylindrical(N_phi))

  allocate(V1(N_phi))
  allocate(V2(N_phi))
  allocate(V3(N_phi))
  allocate(qs(N_phi))
  allocate(qc(N_phi))
  allocate(rs(N_phi))
  allocate(rc(N_phi))

  V1 = X1c * X1c + Y1c * Y1c + Y1s * Y1s
  V2 = 2 * Y1s * Y1c
  V3 = X1c * X1c + Y1c * Y1c - Y1s * Y1s

  ! The "matmul"s that follow could be sped up using BLAS2
  factor = - B0_over_abs_G0 / 8;
  Z20 = factor*matmul(d_d_zeta,V1)
  Z2s = factor*(matmul(d_d_zeta,V2) - 2 * iota_N * V3)
  Z2c = factor*(matmul(d_d_zeta,V3) + 2 * iota_N * V2)

  qs = -iota_N * X1c - Y1s * torsion * abs_G0_over_B0
  qc = matmul(d_d_zeta,X1c) - Y1c * torsion * abs_G0_over_B0
  rs = matmul(d_d_zeta,Y1s) - iota_N * Y1c
  rc = matmul(d_d_zeta,Y1c) + iota_N * Y1s + X1c * torsion * abs_G0_over_B0

  X2s = B0_over_abs_G0 * (matmul(d_d_zeta,Z2s) - 2*iota_N*Z2c + B0_over_abs_G0 * ( abs_G0_over_B0*abs_G0_over_B0*B2s/B0 + (qc * qs + rc * rs)/2)) / curvature

  X2c = B0_over_abs_G0 * (matmul(d_d_zeta,Z2c) + 2*iota_N*Z2s - B0_over_abs_G0 * (-abs_G0_over_B0*abs_G0_over_B0*B2c/B0 &
       + abs_G0_over_B0*abs_G0_over_B0*eta_bar*eta_bar/2 - (qc * qc - qs * qs + rc * rc - rs * rs)/4)) / curvature

  beta_1s = -4 * sign_psi * sign_G * mu0 * p2 * eta_bar * abs_G0_over_B0 / (iota_N * B0 * B0)

  allocate(Y2s_from_X20(N_phi))
  allocate(Y2s_inhomogeneous(N_phi))
  allocate(Y2c_from_X20(N_phi))
  allocate(Y2c_inhomogeneous(N_phi))

  allocate(fX0_from_X20(N_phi))
  allocate(fX0_from_Y20(N_phi))
  allocate(fX0_inhomogeneous(N_phi))

  allocate(fXs_from_X20(N_phi))
  allocate(fXs_from_Y20(N_phi))
  allocate(fXs_inhomogeneous(N_phi))

  allocate(fXc_from_X20(N_phi))
  allocate(fXc_from_Y20(N_phi))
  allocate(fXc_inhomogeneous(N_phi))

  allocate(fY0_from_X20(N_phi))
  allocate(fY0_from_Y20(N_phi))
  allocate(fY0_inhomogeneous(N_phi))

  allocate(fYs_from_X20(N_phi))
  allocate(fYs_from_Y20(N_phi))
  allocate(fYs_inhomogeneous(N_phi))

  allocate(fYc_from_X20(N_phi))
  allocate(fYc_from_Y20(N_phi))
  allocate(fYc_inhomogeneous(N_phi))


  Y2s_from_X20 = -sign_G * sign_psi * curvature * curvature / (eta_bar * eta_bar)
  Y2s_inhomogeneous = sign_G * sign_psi * (-curvature/2 + curvature*curvature/(eta_bar*eta_bar)*(-X2c + X2s * sigma))

  Y2c_from_X20 = -sign_G * sign_psi * curvature * curvature * sigma / (eta_bar * eta_bar)
  Y2c_inhomogeneous = sign_G * sign_psi * curvature * curvature / (eta_bar * eta_bar) * (X2s + X2c * sigma)

  ! Note: in the fX* and fY* quantities below, I've omitted the contributions from X20 and Y20 to the d/dzeta terms. These contributions are
  ! handled later when we assemble the large matrix.

  fX0_from_X20 = -4 * sign_G * sign_psi * abs_G0_over_B0 * (Y2c_from_X20 * Z2s - Y2s_from_X20 * Z2c)
  fX0_from_Y20 = -torsion * abs_G0_over_B0 - 4 * sign_G * sign_psi * abs_G0_over_B0 * (Z2s) &
       - sign_psi * I2_over_B0 * (-2) * abs_G0_over_B0
  fX0_inhomogeneous = curvature * abs_G0_over_B0 * Z20 - 4 * sign_G * sign_psi * abs_G0_over_B0 * (Y2c_inhomogeneous * Z2s - Y2s_inhomogeneous * Z2c) &
       - sign_psi * I2_over_B0 * (0.5d+0 * curvature * sign_G * sign_psi) * abs_G0_over_B0 + beta_1s * abs_G0_over_B0 / 2 * Y1c

  fXs_from_X20 = -torsion * abs_G0_over_B0 * Y2s_from_X20 - 4 * sign_psi * sign_G * abs_G0_over_B0 * (Y2c_from_X20 * Z20) &
       - sign_psi * I2_over_B0 * (- 2 * Y2s_from_X20) * abs_G0_over_B0
  fXs_from_Y20 = - 4 * sign_psi * sign_G * abs_G0_over_B0 * (-Z2c + Z20)
  fXs_inhomogeneous = matmul(d_d_zeta,X2s) - 2 * iota_N * X2c - torsion * abs_G0_over_B0 * Y2s_inhomogeneous + curvature * abs_G0_over_B0 * Z2s &
       - 4 * sign_psi * sign_G * abs_G0_over_B0 * (Y2c_inhomogeneous * Z20) &
       - sign_psi * I2_over_B0 * (0.5d+0 * curvature * sign_psi * sign_G - 2 * Y2s_inhomogeneous) * abs_G0_over_B0 &
       - (0.5d+0) * abs_G0_over_B0 * beta_1s * Y1s

  fXc_from_X20 = - torsion * abs_G0_over_B0 * Y2c_from_X20 - 4 * sign_psi * sign_G * abs_G0_over_B0 * (-Y2s_from_X20 * Z20) &
       - sign_psi * I2_over_B0 * (- 2 * Y2c_from_X20) * abs_G0_over_B0
  fXc_from_Y20 = - torsion * abs_G0_over_B0 - 4 * sign_psi * sign_G * abs_G0_over_B0 * (Z2s) &
       - sign_psi * I2_over_B0 * (-2) * abs_G0_over_B0
  fXc_inhomogeneous = matmul(d_d_zeta,X2c) + 2 * iota_N * X2s - torsion * abs_G0_over_B0 * Y2c_inhomogeneous + curvature * abs_G0_over_B0 * Z2c &
       - 4 * sign_psi * sign_G * abs_G0_over_B0 * (-Y2s_inhomogeneous * Z20) &
       - sign_psi * I2_over_B0 * (0.5d+0 * curvature * sign_G * sign_psi - 2 * Y2c_inhomogeneous) * abs_G0_over_B0 &
       - (0.5d+0) * abs_G0_over_B0 * beta_1s * Y1c

  fY0_from_X20 = torsion * abs_G0_over_B0 - sign_psi * I2_over_B0 * (2) * abs_G0_over_B0
  fY0_from_Y20 = 0
  fY0_inhomogeneous = -4 * sign_psi * sign_G * abs_G0_over_B0 * (X2s * Z2c - X2c * Z2s) &
       - sign_psi * I2_over_B0 * (-0.5d+0 * curvature * X1c * X1c) * abs_G0_over_B0 - (0.5d+0) * abs_G0_over_B0 * beta_1s * X1c

  fYs_from_X20 = -2 * iota_N * Y2c_from_X20 - 4 * sign_psi * sign_G * abs_G0_over_B0 * (Z2c)
  fYs_from_Y20 = -2 * iota_N
  fYs_inhomogeneous = matmul(d_d_zeta,Y2s_inhomogeneous) - 2 * iota_N * Y2c_inhomogeneous + torsion * abs_G0_over_B0 * X2s &
       - 4 * sign_psi * sign_G * abs_G0_over_B0 * (-X2c * Z20) - 2 * sign_psi * I2_over_B0 * X2s * abs_G0_over_B0

  fYc_from_X20 = 2 * iota_N * Y2s_from_X20 - 4 * sign_psi * sign_G * abs_G0_over_B0 * (-Z2s)
  fYc_from_Y20 = 0
  fYc_inhomogeneous = matmul(d_d_zeta,Y2c_inhomogeneous) + 2 * iota_N * Y2s_inhomogeneous + torsion * abs_G0_over_B0 * X2c &
       - 4 * sign_psi * sign_G * abs_G0_over_B0 * (X2s * Z20) &
       - sign_psi * I2_over_B0 * (-0.5d+0 * curvature * X1c * X1c + 2 * X2c) * abs_G0_over_B0 + 0.5d+0 * abs_G0_over_B0 * beta_1s * X1c

  allocate(matrix(2*N_phi,2*N_phi))
  allocate(right_hand_side(2*N_phi))

  matrix = 0
  do j = 1, N_phi
     ! Handle the terms involving d X_0 / d zeta and d Y_0 / d zeta:
     ! ----------------------------------------------------------------

     ! Equation 1, terms involving X0:
     ! Contributions arise from Y1c * fYs - Y1s * fYc.
     matrix(j,1:N_phi) = Y1c(j) * d_d_zeta(j,:) * Y2s_from_X20 - Y1s(j) * d_d_zeta(j,:) * Y2c_from_X20

     ! Equation 1, terms involving Y0:
     ! Contributions arise from -Y1s * fY0 - Y1s * fYc, and they happen to be equal.
     matrix(j,(N_phi+1):(2*N_phi)) = -2 * Y1s(j) * d_d_zeta(j,:)

     ! Equation 2, terms involving X0:
     ! Contributions arise from -X1c * fX0 + Y1s * fYs + Y1c * fYc
     matrix(j+N_phi,1:N_phi) = -X1c(j) * d_d_zeta(j,:) + Y1s(j) * d_d_zeta(j,:) * Y2s_from_X20 + Y1c(j) * d_d_zeta(j,:) * Y2c_from_X20

     ! Equation 2, terms involving Y0:
     ! Contributions arise from -Y1c * fY0 + Y1c * fYc, but they happen to cancel.

     ! Now handle the terms involving X_0 and Y_0 without d/dzeta derivatives:
     ! ----------------------------------------------------------------

     matrix(j,j      ) = matrix(j,j      ) + X1c(j) * fXs_from_X20(j) - Y1s(j) * fY0_from_X20(j) + Y1c(j) * fYs_from_X20(j) - Y1s(j) * fYc_from_X20(j)
     matrix(j,j+N_phi) = matrix(j,j+N_phi) + X1c(j) * fXs_from_Y20(j) - Y1s(j) * fY0_from_Y20(j) + Y1c(j) * fYs_from_Y20(j) - Y1s(j) * fYc_from_Y20(j)

     matrix(j+N_phi,j      ) = matrix(j+N_phi,j      ) - X1c(j) * fX0_from_X20(j) + X1c(j) * fXc_from_X20(j) - Y1c(j) * fY0_from_X20(j) + Y1s(j) * fYs_from_X20(j) + Y1c(j) * fYc_from_X20(j)
     matrix(j+N_phi,j+N_phi) = matrix(j+N_phi,j+N_phi) - X1c(j) * fX0_from_Y20(j) + X1c(j) * fXc_from_Y20(j) - Y1c(j) * fY0_from_Y20(j) + Y1s(j) * fYs_from_Y20(j) + Y1c(j) * fYc_from_Y20(j)
  end do

  right_hand_side(1:N_phi) = -(X1c * fXs_inhomogeneous - Y1s * fY0_inhomogeneous + Y1c * fYs_inhomogeneous - Y1s * fYc_inhomogeneous)
  right_hand_side((N_phi+1):(2*N_phi)) = -(- X1c * fX0_inhomogeneous + X1c * fXc_inhomogeneous - Y1c * fY0_inhomogeneous + Y1s * fYs_inhomogeneous + Y1c * fYc_inhomogeneous)

!!$  print *,"Here comes abs_G0_over_B0:",abs_G0_over_B0
!!$  print *,"Here comes X1c:"
!!$  print *,X1c
!!$  print *,"Here comes X1s:"
!!$  print *,X1s
!!$  print *,"Here comes Y1c:"
!!$  print *,Y1c
!!$  print *,"Here comes Y1s:"
!!$  print *,Y1s
!!$  print *,"Here comes torsion:"
!!$  print *,torsion
!!$
!!$  print *,"Here comes Z20:"
!!$  print *,Z20
!!$  print *,"Here comes Z2s:"
!!$  print *,Z2s
!!$  print *,"Here comes Z2c:"
!!$  print *,Z2c
!!$
!!$  print *,"Here comes qs:"
!!$  print *,qs
!!$  print *,"Here comes qc:"
!!$  print *,qc
!!$  print *,"Here comes rs:"
!!$  print *,rs
!!$  print *,"Here comes rc:"
!!$  print *,rc
!!$
!!$  print *,"Here comes X2c:"
!!$  print *,X2c
!!$  print *,"Here comes X2s:"
!!$  print *,X2s
!!$  print *,"Here comes Y2c_inhomogeneous:"
!!$  print *,Y2c_inhomogeneous
!!$  print *,"Here comes Y2s_inhomogeneous:"
!!$  print *,Y2s_inhomogeneous
!!$  print *," "
!!$  print *,"Here comes fX0_inhomogeneous:"
!!$  print *,fX0_inhomogeneous
!!$  print *,"Here comes fXs_inhomogeneous:"
!!$  print *,fXs_inhomogeneous
!!$  print *,"Here comes fXc_inhomogeneous:"
!!$  print *,fXc_inhomogeneous
!!$  print *,"Here comes fY0_inhomogeneous:"
!!$  print *,fY0_inhomogeneous
!!$  print *,"Here comes fYs_inhomogeneous:"
!!$  print *,fYs_inhomogeneous
!!$  print *,"Here comes fYc_inhomogeneous:"
!!$  print *,fYc_inhomogeneous
!!$
!!$  print *,"Here comes right_hand_side:"
!!$  print *,right_hand_side
!!$
!!$  open(unit=iunit,file="matrix.dat")
!!$  do j = 1, N_phi*2
!!$     write(iunit, "(*(es24.15))") matrix(j,:)
!!$  end do
!!$  close(iunit)

  ! We will use the LAPACK subroutine DGESV to solve a general (asymmetric) linear system
  ! solution = matrix \ right_hand_side
  ! Note that LAPACK will over-write "right_hand_side" with the solution, and over-write "matrix" with the LU factorization.
  allocate(IPIV(2*N_phi))
  call DGESV(2*N_phi, 1, matrix, 2*N_phi, IPIV, right_hand_side, 2*N_phi, INFO)
  deallocate(IPIV)
  if (INFO /= 0) then
     print *, "Error in LAPACK call DGESV: info = ", INFO
     stop
  end if

  X20 = right_hand_side(1:N_phi)
  Y20 = right_hand_side((N_phi+1):(2*N_phi))

  ! Now that we have X20 and Y20 explicitly, we can reconstruct Y2s, Y2c, and B20:
  Y2s = Y2s_inhomogeneous + Y2s_from_X20 * X20
  Y2c = Y2c_inhomogeneous + Y2c_from_X20 * X20 + Y20

  B20 = B0 * (curvature * X20 - B0_over_abs_G0 * matmul(d_d_zeta,Z20) + (0.5d+0) * eta_bar * eta_bar - mu0 * p2 / (B0 * B0) &
       - (0.25d+0) * B0_over_abs_G0 * B0_over_abs_G0 * (qc * qc + qs * qs + rc * rc + rs * rs))

!!$  print *,"Here comes X20:"
!!$  print *,X20
!!$  print *,"Here comes Y20:"
!!$  print *,Y20
!!$  print *,"Here comes Y2s:"
!!$  print *,Y2s
!!$  print *,"Here comes Y2c:"
!!$  print *,Y2c
!!$  print *,"Here comes B20:"
!!$  print *,B20

!  if (.true.) then
  if (.false.) then
     ! For a sanity test, compute the residuals of two equations in a more direct way (now that X20 and Y20 are available explicitly) to make sure we get 0.
     allocate(fX0(N_phi))
     allocate(fXs(N_phi))
     allocate(fXc(N_phi))
     allocate(fY0(N_phi))
     allocate(fYs(N_phi))
     allocate(fYc(N_phi))
     allocate(eq1residual(N_phi))
     allocate(eq2residual(N_phi))

     fX0 = matmul(d_d_zeta,X20) - torsion * abs_G0_over_B0 * Y20 + curvature * abs_G0_over_B0 * Z20 &
          -4*sign_G*sign_psi*abs_G0_over_B0*(Y2c * Z2s - Y2s * Z2c) &
          - sign_psi * I2_over_B0 * (curvature/2 * X1c * Y1c - 2 * Y20) * abs_G0_over_B0 + abs_G0_over_B0 * beta_1s * Y1c / 2

     fXs = matmul(d_d_zeta,X2s) - 2 * iota_N * X2c - torsion * abs_G0_over_B0 * Y2s + curvature * abs_G0_over_B0 * Z2s &
          -4*sign_G*sign_psi*abs_G0_over_B0*(-Y20 * Z2c + Y2c * Z20) &
          -sign_psi * I2_over_B0 * (curvature/2 * X1c * Y1s - 2 * Y2s) * abs_G0_over_B0 - abs_G0_over_B0 * beta_1s * Y1s / 2

     fXc = matmul(d_d_zeta,X2c) + 2 * iota_N * X2s - torsion * abs_G0_over_B0 * Y2c + curvature * abs_G0_over_B0 * Z2c &
          -4*sign_G*sign_psi*abs_G0_over_B0*(Y20 * Z2s - Y2s * Z20) &
          -sign_psi * I2_over_B0 * (curvature/2 * X1c * Y1c - 2 * Y2c) * abs_G0_over_B0 - abs_G0_over_B0 * beta_1s * Y1c / 2

     fY0 = matmul(d_d_zeta,Y20) + torsion * abs_G0_over_B0 * X20 - 4*sign_G*sign_psi*abs_G0_over_B0*(X2s * Z2c - X2c * Z2s) &
          -sign_psi * I2_over_B0 * (-curvature/2*X1c*X1c + 2*X20) * abs_G0_over_B0 - abs_G0_over_B0 * beta_1s * X1c / 2

     fYs = matmul(d_d_zeta,Y2s) - 2 * iota_N * Y2c + torsion * abs_G0_over_B0 * X2s &
          -4*sign_G*sign_psi*abs_G0_over_B0*(X20 * Z2c - X2c * Z20) - 2*sign_psi* I2_over_B0 * X2s * abs_G0_over_B0

     fYc = matmul(d_d_zeta,Y2c) + 2 * iota_N * Y2s + torsion * abs_G0_over_B0 * X2c &
          -4*sign_G*sign_psi*abs_G0_over_B0*(X2s * Z20 - X20 * Z2s) &
          -sign_psi * I2_over_B0 * (-curvature/2 * X1c * X1c + 2 * X2c) * abs_G0_over_B0 + abs_G0_over_B0 * beta_1s * X1c / 2


     eq1residual = X1c * fXs - Y1s * fY0 + Y1c * fYs - Y1s * fYc

     eq2residual = -X1c * fX0 + X1c * fXc - Y1c * fY0 + Y1s * fYs + Y1c * fYc

     max_eq1residual = maxval(abs(eq1residual))
     max_eq2residual = maxval(abs(eq2residual))
     print *,"max(abs(eq1residual)):",max_eq1residual
     print *,"max(abs(eq2residual)):",max_eq2residual

     if (max_eq1residual > 1e-8) stop "Equation 1 residual is large !!!"
     if (max_eq2residual > 1e-8) stop "Equation 2 residual is large !!!"

     ! Now check the two equations that were used to determine Y2s and Y2c:

     eq1residual = -X1c * Y2c + X1c * Y20 + X2s * Y1s + X2c * Y1c - X20 * Y1c

     eq2residual = X1c * Y2s + X2c * Y1s - X2s * Y1c + X20 * Y1s + sign_G * sign_psi * X1c * curvature / 2

     max_eq1residual = maxval(abs(eq1residual))
     max_eq2residual = maxval(abs(eq2residual))
     print *,"max(abs(Y2c eq residual)):",max_eq1residual
     print *,"max(abs(Y2s eq residual)):",max_eq2residual

     if (max_eq1residual > 1e-8) stop "Y2c equation residual is large !!!"
     if (max_eq2residual > 1e-8) stop "Y2s equation residual is large !!!"

     deallocate(fX0, fXs, fXc, fY0, fYs, fYc, eq1residual, eq2residual)
  end if

  normalizer = 1 / sum(d_l_d_phi)
  B20_mean = sum(B20 * d_l_d_phi) * normalizer
  B20_residual = sqrt(sum((B20 - B20_mean) * (B20 - B20_mean) * d_l_d_phi) * normalizer) / B0

  deallocate(V1, V2, V3, qs, qc, rs, rc)
  deallocate(matrix, right_hand_side)
  deallocate(Y2s_from_X20, Y2s_inhomogeneous, Y2c_from_X20, Y2c_inhomogeneous)
  deallocate(fX0_from_X20, fX0_from_Y20, fX0_inhomogeneous)
  deallocate(fXs_from_X20, fXs_from_Y20, fXs_inhomogeneous)
  deallocate(fXc_from_X20, fXc_from_Y20, fXc_inhomogeneous)
  deallocate(fY0_from_X20, fY0_from_Y20, fY0_inhomogeneous)
  deallocate(fYs_from_X20, fYs_from_Y20, fYs_inhomogeneous)
  deallocate(fYc_from_X20, fYc_from_Y20, fYc_inhomogeneous)

  N_helicity = - axis_helicity*nfp
  I2 = I2_over_B0 * B0
  G0 = sign_G * abs_G0_over_B0 * B0
  Bbar = sign_psi * B0
  G2 = -mu0 * p2 * G0 / (B0 * B0) - iota * I2
  allocate(d_Z20_d_zeta(N_phi))
  d_Z20_d_zeta = matmul(d_d_zeta,Z20)
  if (allocated(B0_order_a_squared_to_cancel)) deallocate(B0_order_a_squared_to_cancel)
  allocate(B0_order_a_squared_to_cancel(N_phi))
  B0_order_a_squared_to_cancel = -sign_G * B0 * B0 * (G2 + I2 * N_helicity) * abs_G0_over_B0 / (2*G0*G0) &
       -sign_G * sign_psi * B0 * 2 * (X2c * Y2s - X2s * Y2c) &
       -sign_G * B0 * B0 / (2*G0) * (abs_G0_over_B0 * X20 * curvature - d_Z20_d_zeta) &
       -sign_G * sign_psi * B0 * I2 / (4*G0) * (-abs_G0_over_B0 * torsion * (X1c*X1c + Y1c*Y1c + Y1s*Y1s) + Y1c * d_X1c_d_zeta - X1c * d_Y1c_d_zeta)


!  print *,"AAA"
  if (trim(order_r_option) == order_r_option_r2) then
     deallocate(d_Z20_d_zeta)
     return
  end if
!  print *,"BBB"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Beginning of O(r^3) calculation
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if( trim(order_r_option) == order_r_option_r3_X3s3_Y3s3 &
       .or. trim(order_r_option) == order_r_option_r3_X3c3_Y3c3 &
       .or. trim(order_r_option) == order_r_option_r3_Y3s3_Y3c3) then
     stop "This order_r_option is not valid."
  end if

  if (allocated(X3s1)) deallocate(X3s1)
  if (allocated(X3s3)) deallocate(X3s3)
  if (allocated(X3c1)) deallocate(X3c1)
  if (allocated(X3c3)) deallocate(X3c3)
  if (allocated(Y3s1)) deallocate(Y3s1)
  if (allocated(Y3s3)) deallocate(Y3s3)
  if (allocated(Y3c1)) deallocate(Y3c1)
  if (allocated(Y3c3)) deallocate(Y3c3)
  if (allocated(Z3s1)) deallocate(Z3s1)
  if (allocated(Z3s3)) deallocate(Z3s3)
  if (allocated(Z3c1)) deallocate(Z3c1)
  if (allocated(Z3c3)) deallocate(Z3c3)
  if (allocated(R3s1)) deallocate(R3s1)
  if (allocated(R3s3)) deallocate(R3s3)
  if (allocated(R3c1)) deallocate(R3c1)
  if (allocated(R3c3)) deallocate(R3c3)
  if (allocated(z3s1_cylindrical)) deallocate(z3s1_cylindrical)
  if (allocated(z3s3_cylindrical)) deallocate(z3s3_cylindrical)
  if (allocated(z3c1_cylindrical)) deallocate(z3c1_cylindrical)
  if (allocated(z3c3_cylindrical)) deallocate(z3c3_cylindrical)

  if (allocated(X3s1_untwisted)) deallocate(X3s1_untwisted)
  if (allocated(X3s3_untwisted)) deallocate(X3s3_untwisted)
  if (allocated(X3c1_untwisted)) deallocate(X3c1_untwisted)
  if (allocated(X3c3_untwisted)) deallocate(X3c3_untwisted)
  if (allocated(Y3s1_untwisted)) deallocate(Y3s1_untwisted)
  if (allocated(Y3s3_untwisted)) deallocate(Y3s3_untwisted)
  if (allocated(Y3c1_untwisted)) deallocate(Y3c1_untwisted)
  if (allocated(Y3c3_untwisted)) deallocate(Y3c3_untwisted)
  if (allocated(Z3s1_untwisted)) deallocate(Z3s1_untwisted)
  if (allocated(Z3s3_untwisted)) deallocate(Z3s3_untwisted)
  if (allocated(Z3c1_untwisted)) deallocate(Z3c1_untwisted)
  if (allocated(Z3c3_untwisted)) deallocate(Z3c3_untwisted)

  allocate(X3s1(N_phi))
  allocate(X3s3(N_phi))
  allocate(X3c1(N_phi))
  allocate(X3c3(N_phi))
  allocate(Y3s1(N_phi))
  allocate(Y3s3(N_phi))
  allocate(Y3c1(N_phi))
  allocate(Y3c3(N_phi))
  allocate(Z3s1(N_phi))
  allocate(Z3s3(N_phi))
  allocate(Z3c1(N_phi))
  allocate(Z3c3(N_phi))
  allocate(R3s1(N_phi))
  allocate(R3s3(N_phi))
  allocate(R3c1(N_phi))
  allocate(R3c3(N_phi))
  allocate(z3s1_cylindrical(N_phi))
  allocate(z3s3_cylindrical(N_phi))
  allocate(z3c1_cylindrical(N_phi))
  allocate(z3c3_cylindrical(N_phi))
  allocate(X3s1_untwisted(N_phi))
  allocate(X3s3_untwisted(N_phi))
  allocate(X3c1_untwisted(N_phi))
  allocate(X3c3_untwisted(N_phi))
  allocate(Y3s1_untwisted(N_phi))
  allocate(Y3s3_untwisted(N_phi))
  allocate(Y3c1_untwisted(N_phi))
  allocate(Y3c3_untwisted(N_phi))
  allocate(Z3s1_untwisted(N_phi))
  allocate(Z3s3_untwisted(N_phi))
  allocate(Z3c1_untwisted(N_phi))
  allocate(Z3c3_untwisted(N_phi))

  ! Derivatives of Z20, Z2c, and Z2s are needed to compute (X3,Y3).
  allocate(d_Z2c_d_zeta(N_phi))
  allocate(d_Z2s_d_zeta(N_phi))
  d_Z2c_d_zeta = matmul(d_d_zeta,Z2c)
  d_Z2s_d_zeta = matmul(d_d_zeta,Z2s)

 ! print *,"DDD"
  if (trim(order_r_option) == order_r_option_r3_flux_constraint .or. trim(order_r_option) == order_r_option_r3_flux_constraint_const_B20) then
     X3s1 = 0
     X3s3 = 0
     X3c3 = 0

     Y3c3 = 0
     Y3s3 = 0

     Z3s1 = 0
     Z3s3 = 0
     Z3c1 = 0
     Z3c3 = 0

     allocate(flux_constraint_coefficient(N_phi))

     ! The formula below is copied from "20190305-01 GarrenBoozer r2 corrected radius.nb"
     B1c = B0 * eta_bar
     flux_constraint_coefficient = (-4*B0**2*G0*X20**2*Y1c**2 + 8*B0**2*G0*X20*X2c*Y1c**2 - 4*B0**2*G0*X2c**2*Y1c**2 - &
           4*B0**2*G0*X2s**2*Y1c**2 + 8*B0*G0*B1c*X1c*X2s*Y1c*Y1s + 16*B0**2*G0*X20*X2s*Y1c*Y1s + &
           2*B0**2*I2*iota*X1c**2*Y1s**2 - G0*B1c**2*X1c**2*Y1s**2 - 4*B0*G0*B20*X1c**2*Y1s**2 - &
           8*B0*G0*B1c*X1c*X20*Y1s**2 - 4*B0**2*G0*X20**2*Y1s**2 - 8*B0*G0*B1c*X1c*X2c*Y1s**2 - &
           8*B0**2*G0*X20*X2c*Y1s**2 - 4*B0**2*G0*X2c**2*Y1s**2 - 4*B0**2*G0*X2s**2*Y1s**2 + &
           8*B0**2*G0*X1c*X20*Y1c*Y20 - 8*B0**2*G0*X1c*X2c*Y1c*Y20 - 8*B0**2*G0*X1c*X2s*Y1s*Y20 - &
           4*B0**2*G0*X1c**2*Y20**2 - 8*B0**2*G0*X1c*X20*Y1c*Y2c + 8*B0**2*G0*X1c*X2c*Y1c*Y2c + &
           24*B0**2*G0*X1c*X2s*Y1s*Y2c + 8*B0**2*G0*X1c**2*Y20*Y2c - 4*B0**2*G0*X1c**2*Y2c**2 + &
           8*B0**2*G0*X1c*X2s*Y1c*Y2s - 8*B0*G0*B1c*X1c**2*Y1s*Y2s - 8*B0**2*G0*X1c*X20*Y1s*Y2s - &
           24*B0**2*G0*X1c*X2c*Y1s*Y2s - 4*B0**2*G0*X1c**2*Y2s**2 - 4*B0**2*G0*X1c**2*Z20**2 - &
           4*B0**2*G0*Y1c**2*Z20**2 - 4*B0**2*G0*Y1s**2*Z20**2 - 4*B0**2*abs_G0_over_B0*I2*Y1c*Y1s*Z2c + &
           8*B0**2*G0*X1c**2*Z20*Z2c + 8*B0**2*G0*Y1c**2*Z20*Z2c - 8*B0**2*G0*Y1s**2*Z20*Z2c - &
           4*B0**2*G0*X1c**2*Z2c**2 - 4*B0**2*G0*Y1c**2*Z2c**2 - 4*B0**2*G0*Y1s**2*Z2c**2 + &
           2*B0**2*abs_G0_over_B0*I2*X1c**2*Z2s + 2*B0**2*abs_G0_over_B0*I2*Y1c**2*Z2s - 2*B0**2*abs_G0_over_B0*I2*Y1s**2*Z2s + &
           16*B0**2*G0*Y1c*Y1s*Z20*Z2s - 4*B0**2*G0*X1c**2*Z2s**2 - 4*B0**2*G0*Y1c**2*Z2s**2 - &
           4*B0**2*G0*Y1s**2*Z2s**2 + B0**2*abs_G0_over_B0*I2*X1c**3*Y1s*torsion + B0**2*abs_G0_over_B0*I2*X1c*Y1c**2*Y1s*torsion + &
           B0**2*abs_G0_over_B0*I2*X1c*Y1s**3*torsion - B0**2*I2*X1c*Y1c*Y1s*d_X1c_d_zeta + &
           B0**2*I2*X1c**2*Y1s*d_Y1c_d_zeta)/(16*B0**2*G0*X1c**2*Y1s**2)

     if (.true.) then
        ! Check equations from paper:
        allocate(Q(N_phi))
        allocate(predicted_flux_constraint_coefficient(N_phi))
        
        Q = sign_psi * B0 * abs_G0_over_B0 / (2*G0*G0) * (iota_N * I2 + mu0 * p2 * G0 / (B0 * B0)) + 2 * (X2c * Y2s - X2s * Y2c) &
             + sign_psi * B0 / (2*G0) * (abs_G0_over_B0 * X20 * curvature - d_Z20_d_zeta) &
             + I2 / (4 * G0) * (-abs_G0_over_B0 * torsion * (X1c*X1c + Y1s*Y1s + Y1c*Y1c) + Y1c * d_X1c_d_zeta - X1c * d_Y1c_d_zeta)
        predicted_flux_constraint_coefficient = - Q / (2 * sign_G * sign_psi)
!!$     print *,"flux_constraint_coefficient:"
!!$     print *,flux_constraint_coefficient
        print *,"flux_constraint_coefficient - predicted_flux_constraint_coefficient:"
        print *,flux_constraint_coefficient - predicted_flux_constraint_coefficient
        
!!$     print *,flux_constraint_coefficient - B0_order_a_squared_to_cancel/(2*B0)

        deallocate(Q,predicted_flux_constraint_coefficient)
     end if

     if (trim(order_r_option) == order_r_option_r3_flux_constraint_const_B20) then
        flux_constraint_coefficient = flux_constraint_coefficient + (B20 - B20_mean) / (2 * B0)
     end if

     X3c1 = X1c * flux_constraint_coefficient
     Y3c1 = Y1c * flux_constraint_coefficient
     Y3s1 = Y1s * flux_constraint_coefficient

     deallocate(flux_constraint_coefficient)
  end if
!  print *,"FFF"
  if (trim(order_r_option) == order_r_option_r3_simplified .or. trim(order_r_option) == order_r_option_r3_simplified_with_Z3) then
     X3s1 = 0
     X3s3 = 0
     X3c3 = 0

     Y3c3 = 0
     Y3s3 = 0

     Z3s1 = 0
     Z3s3 = 0
     Z3c1 = 0
     Z3c3 = 0

     X3c1 = -(1/(8 * G0 * Y1s)) * (2 * I2 * N_helicity * X1c * Y1s + 2 * (-I2 * iota - (G0 * p2 * mu0)/ (B0*B0)) * X1c * Y1s - &
          16 * G0 * X2s * Y20 - 8 * G0 * X2s * Y2c + 16 * G0 * X20 * Y2s + 8 * G0 * X2c * Y2s - 8 * B0 * sign_psi * (iota - N_helicity) * Z2s + &
          2 * B0 * abs_G0_over_B0 * sign_psi * X20 * curvature + 4 * B0 * abs_G0_over_B0 * sign_psi * X2c * curvature - 3 * abs_G0_over_B0 * I2 * X1c**2 * torsion - &
          3 * abs_G0_over_B0 * I2 * Y1c**2 * torsion + abs_G0_over_B0 * I2 * Y1s**2 * torsion + 3 * I2 * Y1c * d_X1c_d_zeta - 3 * I2 * X1c * d_Y1c_d_zeta - &
          2 * B0 * sign_psi * d_Z20_d_zeta - 4 * B0 * sign_psi * d_Z2c_d_zeta)
     
     
     Y3c1 =-(1/(8 * G0 * X1c * Y1s)) * (2 * I2 * N_helicity * X1c * Y1c * Y1s + 2 * (-I2 * iota - (G0 * p2 * mu0)/(B0**2)) * X1c * Y1c * Y1s - &
          16 * G0 * X2s * Y1c * Y20 + 32 * G0 * X2c * Y1s * Y20 - 8 * G0 * X2s * Y1c * Y2c - 32 * G0 * X20 * Y1s * Y2c + 16 * G0 * X20 * Y1c * Y2s + &
          8 * G0 * X2c * Y1c * Y2s + 16 * B0 * sign_psi * (iota - N_helicity) * Y1s * Z2c - 8 * B0 * sign_psi * (iota - N_helicity) * Y1c * Z2s + &
          2 * B0 * abs_G0_over_B0 * sign_psi * X20 * Y1c * curvature + 4 * B0 * abs_G0_over_B0 * sign_psi * X2c * Y1c * curvature + &
          8 * B0 * abs_G0_over_B0 * sign_psi * X2s * Y1s * curvature - 3 * abs_G0_over_B0 * I2 * X1c**2 * Y1c * torsion - 3 * abs_G0_over_B0 * I2 * Y1c**3 * torsion - &
          7 * abs_G0_over_B0 * I2 * Y1c * Y1s**2 * torsion + 3 * I2 * Y1c**2 * d_X1c_d_zeta + 4 * I2 * Y1s**2 * d_X1c_d_zeta - 3 * I2 * X1c * Y1c * d_Y1c_d_zeta - &
          4 * I2 * X1c * Y1s * d_Y1s_d_zeta - 2 * B0 * sign_psi * Y1c * d_Z20_d_zeta - 4 * B0 * sign_psi * Y1c * d_Z2c_d_zeta - 8 * B0 * sign_psi * Y1s * d_Z2s_d_zeta)
     
     
     Y3s1 = -(1/(8 * G0 * X1c)) * (2 * I2 * N_helicity * X1c * Y1s + 2 * (-I2 * iota - (G0 * p2 * mu0)/(B0**2)) * X1c * Y1s + 16 * G0 * X2s * Y20 - &
          8 * G0 * X2s * Y2c - 16 * G0 * X20 * Y2s + 8 * G0 * X2c * Y2s + 8 * B0 * sign_psi * (iota - N_helicity) * Z2s + 2 * B0 * abs_G0_over_B0 * sign_psi * X20 * curvature &
          - 4 * B0 * abs_G0_over_B0 * sign_psi * X2c * curvature + abs_G0_over_B0 * I2 * X1c**2 * torsion + abs_G0_over_B0 * I2 * Y1c**2 * torsion - &
          3 * abs_G0_over_B0 * I2 * Y1s**2 * torsion - I2 * Y1c * d_X1c_d_zeta + I2 * X1c * d_Y1c_d_zeta - 2 * B0 * sign_psi * d_Z20_d_zeta + 4 * B0 * sign_psi * d_Z2c_d_zeta)
  end if
 ! print *,"HHH"
  if (trim(order_r_option) == order_r_option_r3_simplified_with_Z3 &
       .or. trim(order_r_option) == order_r_option_r3_B3 &
       .or. trim(order_r_option) == order_r_option_r3_X3s3_X3c3 &
       .or. trim(order_r_option) == order_r_option_r3_X3s3_Y3s3 &
       .or. trim(order_r_option) == order_r_option_r3_X3c3_Y3c3 &
       .or. trim(order_r_option) == order_r_option_r3_Y3s3_Y3c3) then
     ! Compute Z3 terms
     ! These are worked out in "20190318-01 Wrick's streamlined Garren-Boozer method, MHD.nb"

     allocate(d_X20_d_zeta(N_phi))
     allocate(d_X2c_d_zeta(N_phi))
     allocate(d_X2s_d_zeta(N_phi))
     allocate(d_Y20_d_zeta(N_phi))
     allocate(d_Y2c_d_zeta(N_phi))
     allocate(d_Y2s_d_zeta(N_phi))
     d_X20_d_zeta = matmul(d_d_zeta,X20)
     d_X2c_d_zeta = matmul(d_d_zeta,X2c)
     d_X2s_d_zeta = matmul(d_d_zeta,X2s)
     d_Y20_d_zeta = matmul(d_d_zeta,Y20)
     d_Y2c_d_zeta = matmul(d_d_zeta,Y2c)
     d_Y2s_d_zeta = matmul(d_d_zeta,Y2s)

     Z3s1 = (8*iota_N*X1c*X20 + 8*iota_N*Y1c*Y20 + 4*beta_1s*abs_G0_over_B0*X1c*Y1s + &
          iota_N*X1c**3*curvature + iota_N*X1c*Y1c**2*curvature - iota_N*X1c*Y1s**2*curvature - &
          2*abs_G0_over_B0*X1c*Z2s*curvature + 2*abs_G0_over_B0*X2s*Y1c*torsion + 4*abs_G0_over_B0*X20*Y1s*torsion - &
          2*abs_G0_over_B0*X2c*Y1s*torsion - 2*abs_G0_over_B0*X1c*Y2s*torsion - 4*X2s*d_X1c_d_zeta - &
          2*X1c*d_X2s_d_zeta - 4*Y2s*d_Y1c_d_zeta - &
          X1c*Y1s*curvature*d_Y1c_d_zeta - 8*Y20*d_Y1s_d_zeta + &
          4*Y2c*d_Y1s_d_zeta - X1c*Y1c*curvature*d_Y1s_d_zeta - &
          4*Y1s*d_Y20_d_zeta + 2*Y1s*d_Y2c_d_zeta - 2*Y1c*d_Y2s_d_zeta)/ &
        (12*abs_G0_over_B0)

     Z3c1 = (-8*iota_N*Y1s*Y20 - 2*iota_N*X1c*Y1c*Y1s*curvature - 4*abs_G0_over_B0*X1c*Z20*curvature - &
          2*abs_G0_over_B0*X1c*Z2c*curvature + 4*abs_G0_over_B0*X20*Y1c*torsion + 2*abs_G0_over_B0*X2c*Y1c*torsion + &
          2*abs_G0_over_B0*X2s*Y1s*torsion - 4*abs_G0_over_B0*X1c*Y20*torsion - 2*abs_G0_over_B0*X1c*Y2c*torsion - &
          8*X20*d_X1c_d_zeta - 4*X2c*d_X1c_d_zeta - &
          3*X1c**2*curvature*d_X1c_d_zeta - 4*X1c*d_X20_d_zeta - &
          2*X1c*d_X2c_d_zeta - 8*Y20*d_Y1c_d_zeta - 4*Y2c*d_Y1c_d_zeta - &
          3*X1c*Y1c*curvature*d_Y1c_d_zeta - 4*Y2s*d_Y1s_d_zeta - &
          X1c*Y1s*curvature*d_Y1s_d_zeta - 4*Y1c*d_Y20_d_zeta - &
          2*Y1c*d_Y2c_d_zeta - 2*Y1s*d_Y2s_d_zeta)/(12*abs_G0_over_B0)

     Z3s3 = (8*iota_N*X1c*X2c + 8*iota_N*Y1c*Y2c - 8*iota_N*Y1s*Y2s + iota_N*X1c**3*curvature + &
          iota_N*X1c*Y1c**2*curvature - iota_N*X1c*Y1s**2*curvature - 2*abs_G0_over_B0*X1c*Z2s*curvature + &
          2*abs_G0_over_B0*X2s*Y1c*torsion + 2*abs_G0_over_B0*X2c*Y1s*torsion - 2*abs_G0_over_B0*X1c*Y2s*torsion - &
          4*X2s*d_X1c_d_zeta - 2*X1c*d_X2s_d_zeta - 4*Y2s*d_Y1c_d_zeta - &
          X1c*Y1s*curvature*d_Y1c_d_zeta - 4*Y2c*d_Y1s_d_zeta - &
          X1c*Y1c*curvature*d_Y1s_d_zeta - 2*Y1s*d_Y2c_d_zeta - &
          2*Y1c*d_Y2s_d_zeta)/(12*abs_G0_over_B0)

     Z3c3 = (-8*iota_N*X1c*X2s - 8*iota_N*Y1s*Y2c - 8*iota_N*Y1c*Y2s - &
          2*iota_N*X1c*Y1c*Y1s*curvature - 2*abs_G0_over_B0*X1c*Z2c*curvature + 2*abs_G0_over_B0*X2c*Y1c*torsion - &
          2*abs_G0_over_B0*X2s*Y1s*torsion - 2*abs_G0_over_B0*X1c*Y2c*torsion - 4*X2c*d_X1c_d_zeta - &
          X1c**2*curvature*d_X1c_d_zeta - 2*X1c*d_X2c_d_zeta - &
          4*Y2c*d_Y1c_d_zeta - X1c*Y1c*curvature*d_Y1c_d_zeta + &
          4*Y2s*d_Y1s_d_zeta + X1c*Y1s*curvature*d_Y1s_d_zeta - &
          2*Y1c*d_Y2c_d_zeta + 2*Y1s*d_Y2s_d_zeta)/(12*abs_G0_over_B0)

  end if
  !print *,"JJJ"
  if (trim(order_r_option) == order_r_option_r3_B3 &
       .or. trim(order_r_option) == order_r_option_r3_X3s3_X3c3 &
       .or. trim(order_r_option) == order_r_option_r3_X3s3_Y3s3 &
       .or. trim(order_r_option) == order_r_option_r3_X3c3_Y3c3 &
       .or. trim(order_r_option) == order_r_option_r3_Y3s3_Y3c3) then

     allocate(d_X3c3_d_zeta(N_phi))
     allocate(d_X3s3_d_zeta(N_phi))
     allocate(d_Z3c1_d_zeta(N_phi))
     allocate(d_Z3s1_d_zeta(N_phi))
     allocate(d_Z3c3_d_zeta(N_phi))
     allocate(d_Z3s3_d_zeta(N_phi))
     d_Z3c1_d_zeta = matmul(d_d_zeta,Z3c1)
     d_Z3s1_d_zeta = matmul(d_d_zeta,Z3s1)
     d_Z3c3_d_zeta = matmul(d_d_zeta,Z3c3)
     d_Z3s3_d_zeta = matmul(d_d_zeta,Z3s3)
  end if
  !print *,"LLL"
  if (trim(order_r_option) == order_r_option_r3_B3) then
     ! Compute X3s3 and X3c3 from B3s3 and B3c3
 
     ! Equations copied from "20190318-01 Wrick's streamlined Garren-Boozer method, MHD.nb"
     X3c3 = ((G0**2*B1c**3)/B0**5 - (3*G0**2*B1c*B2c)/B0**4 + (2*G0**2*B3c3_input)/B0**3 - &
         2*iota_N**2*X1c*X2c - 2*iota_N**2*Y1c*Y2c + 2*iota_N**2*Y1s*Y2s + 6*abs_G0_over_B0*iota_N*Z3s3 - &
         abs_G0_over_B0*iota_N*X1c*Z2s*curvature + abs_G0_over_B0**2*X1c*X2c*curvature**2 - abs_G0_over_B0*iota_N*X2s*Y1c*torsion - &
         abs_G0_over_B0*iota_N*X2c*Y1s*torsion + abs_G0_over_B0*iota_N*X1c*Y2s*torsion - abs_G0_over_B0**2*Y1c*Z2c*curvature*torsion + &
         abs_G0_over_B0**2*Y1s*Z2s*curvature*torsion + abs_G0_over_B0**2*X1c*X2c*torsion**2 + abs_G0_over_B0**2*Y1c*Y2c*torsion**2 - &
         abs_G0_over_B0**2*Y1s*Y2s*torsion**2 + 2*iota_N*X2s*d_X1c_d_zeta + &
         abs_G0_over_B0*Z2c*curvature*d_X1c_d_zeta - abs_G0_over_B0*Y2c*torsion*d_X1c_d_zeta - &
         abs_G0_over_B0*Y1c*torsion*d_X2c_d_zeta + d_X1c_d_zeta*d_X2c_d_zeta + &
         iota_N*X1c*d_X2s_d_zeta + abs_G0_over_B0*Y1s*torsion*d_X2s_d_zeta + &
         2*iota_N*Y2s*d_Y1c_d_zeta + abs_G0_over_B0*X2c*torsion*d_Y1c_d_zeta + &
         2*iota_N*Y2c*d_Y1s_d_zeta - abs_G0_over_B0*X2s*torsion*d_Y1s_d_zeta + &
         iota_N*Y1s*d_Y2c_d_zeta + abs_G0_over_B0*X1c*torsion*d_Y2c_d_zeta + &
         d_Y1c_d_zeta*d_Y2c_d_zeta + iota_N*Y1c*d_Y2s_d_zeta - &
         d_Y1s_d_zeta*d_Y2s_d_zeta - abs_G0_over_B0*X1c*curvature*d_Z2c_d_zeta + &
         2*abs_G0_over_B0*d_Z3c3_d_zeta)/(2*abs_G0_over_B0**2*curvature)

     X3s3 = ((-3*G0**2*B1c*B2s)/B0**4 + (2*G0**2*B3s3_input)/B0**3 - 2*iota_N**2*X1c*X2s - &
         2*iota_N**2*Y1s*Y2c - 2*iota_N**2*Y1c*Y2s - 6*abs_G0_over_B0*iota_N*Z3c3 + abs_G0_over_B0*iota_N*X1c*Z2c*curvature + &
         abs_G0_over_B0**2*X1c*X2s*curvature**2 + abs_G0_over_B0*iota_N*X2c*Y1c*torsion - abs_G0_over_B0*iota_N*X2s*Y1s*torsion - &
         abs_G0_over_B0*iota_N*X1c*Y2c*torsion - abs_G0_over_B0**2*Y1s*Z2c*curvature*torsion - abs_G0_over_B0**2*Y1c*Z2s*curvature*torsion + &
         abs_G0_over_B0**2*X1c*X2s*torsion**2 + abs_G0_over_B0**2*Y1s*Y2c*torsion**2 + abs_G0_over_B0**2*Y1c*Y2s*torsion**2 - &
         2*iota_N*X2c*d_X1c_d_zeta + abs_G0_over_B0*Z2s*curvature*d_X1c_d_zeta - &
         abs_G0_over_B0*Y2s*torsion*d_X1c_d_zeta - iota_N*X1c*d_X2c_d_zeta - &
         abs_G0_over_B0*Y1s*torsion*d_X2c_d_zeta - abs_G0_over_B0*Y1c*torsion*d_X2s_d_zeta + &
         d_X1c_d_zeta*d_X2s_d_zeta - 2*iota_N*Y2c*d_Y1c_d_zeta + &
         abs_G0_over_B0*X2s*torsion*d_Y1c_d_zeta + 2*iota_N*Y2s*d_Y1s_d_zeta + &
         abs_G0_over_B0*X2c*torsion*d_Y1s_d_zeta - iota_N*Y1c*d_Y2c_d_zeta + &
         d_Y1s_d_zeta*d_Y2c_d_zeta + iota_N*Y1s*d_Y2s_d_zeta + &
         abs_G0_over_B0*X1c*torsion*d_Y2s_d_zeta + d_Y1c_d_zeta*d_Y2s_d_zeta - &
         abs_G0_over_B0*X1c*curvature*d_Z2s_d_zeta + 2*abs_G0_over_B0*d_Z3s3_d_zeta)/(2*abs_G0_over_B0**2*curvature)

     d_X3s3_d_zeta = matmul(d_d_zeta, X3s3)
     d_X3c3_d_zeta = matmul(d_d_zeta, X3c3)
  elseif (trim(order_r_option) == order_r_option_r3_X3s3_X3c3) then
     X3c3 = 0
     X3s3 = 0
     d_X3s3_d_zeta = 0
     d_X3c3_d_zeta = 0
  end if
  !print *,"NNN"
  if (trim(order_r_option) == order_r_option_r3_B3 &
       .or. trim(order_r_option) == order_r_option_r3_X3s3_X3c3 &
       .or. trim(order_r_option) == order_r_option_r3_X3s3_Y3s3 &
       .or. trim(order_r_option) == order_r_option_r3_X3c3_Y3c3 &
       .or. trim(order_r_option) == order_r_option_r3_Y3s3_Y3c3) then
     ! Full O(r^3) solution.

     beta_2c = 0
     beta_2s = 0
     I4 = 0
     
     ! Set up and solve the O(r^3) linear system
     vector_size = 6 * N_phi + 1
     allocate(matrix(vector_size, vector_size))
     allocate(right_hand_side(vector_size))
     matrix = 0

     ! Rows of the vector of unknowns, and columns of the matrix:
     index_X3c1      = 1 + 0 * N_phi
     index_X3s1      = 1 + 1 * N_phi
     index_Y3c1      = 1 + 2 * N_phi
     index_Y3s1      = 1 + 3 * N_phi
     index_Y3c3      = 1 + 4 * N_phi
     index_Y3s3      = 1 + 5 * N_phi
     index_iota2     = 1 + 6 * N_phi
     !index_unknown1  = 1 + 4 * N_phi
     !index_unknown2  = 1 + 5 * N_phi

     ! Rows of the right-hand-side vector, and rows of the matrix:
     index_XYEquation_0            = 1 + 0 * N_phi
     index_XYEquation_s            = 1 + 1 * N_phi
     index_XYEquation_c            = 1 + 2 * N_phi
     index_mixedPartialsEquation_0 = 1 + 3 * N_phi
     index_mixedPartialsEquation_s = 1 + 4 * N_phi
     index_mixedPartialsEquation_c = 1 + 5 * N_phi
     index_initialCondition        = 1 + 6 * N_phi

!!$     if(trim(order_r_option)==order_r_option_r3_B3 .or. trim(order_r_option)==order_r_option_r3_X3s3_X3c3) then
!!$        index_Y3s3 = index_unknown1
!!$        index_Y3c3 = index_unknown2
!!$     elseif(trim(order_r_option)==order_r_option_r3_X3s3_Y3s3) then
!!$        index_X3c3 = index_unknown1
!!$        index_Y3c3 = index_unknown2
!!$     elseif(trim(order_r_option)==order_r_option_r3_X3c3_Y3c3) then
!!$        index_X3s3 = index_unknown1
!!$        index_Y3s3 = index_unknown2
!!$     elseif(trim(order_r_option)==order_r_option_r3_Y3s3_Y3c3) then
!!$        index_X3s3 = index_unknown1
!!$        index_X3c3 = index_unknown2
!!$     else
!!$        stop "Should not get here."
!!$     end if

     ! Impose initial condition for Y3c1
     matrix(index_initialCondition, index_Y3c1) = 1 
     right_hand_side(index_initialCondition) = Y3c1_initial

     ! These next equations come from "20190318-01 Wrick's streamlined Garren-Boozer method, MHD.nb"
        
     ! sin(2 theta) part of the XY equation:
     right_hand_side(index_XYEquation_s:(index_XYEquation_s+N_phi-1)) = -( 4*X2c*Y20 - 4*X20*Y2c + (2*Bbar*iota_N*Z2c)/G0 + &
            (Bbar*abs_G0_over_B0*X2s*curvature)/G0 - (abs_G0_over_B0*I2*Y1c*Y1s*torsion)/G0 + &
            (I2*Y1s*d_X1c_d_zeta)/(2*G0) - (I2*X1c*d_Y1s_d_zeta)/(2*G0) - (Bbar*d_Z2s_d_zeta)/G0) &
            - 3*Y1c*X3c3 - 3*Y1s*X3s3
     do j = 1, N_phi
        matrix(index_XYEquation_s+j-1,index_X3c1+j-1) = -Y1c(j)
        matrix(index_XYEquation_s+j-1,index_X3s1+j-1) = Y1s(j)
        matrix(index_XYEquation_s+j-1,index_Y3c1+j-1) = X1c(j)
        matrix(index_XYEquation_s+j-1,index_Y3c3+j-1) = -3*X1c(j)
     end do

     ! cos(2 theta) part of the XY equation:
     right_hand_side(index_XYEquation_c:(index_XYEquation_c+N_phi-1)) = -(- 4*X2s*Y20 + &
          4*X20*Y2s - (2*Bbar*iota_N*Z2s)/G0 + &
          (Bbar*abs_G0_over_B0*X2c*curvature)/G0 - (abs_G0_over_B0*I2*(X1c**2)*torsion)/(2*G0) - (abs_G0_over_B0*I2*(Y1c**2)*torsion)/(2*G0) + &
          (abs_G0_over_B0*I2*(Y1s**2)*torsion)/(2*G0) + (I2*Y1c*d_X1c_d_zeta)/(2*G0) - &
          (I2*X1c*d_Y1c_d_zeta)/(2*G0) - (Bbar*d_Z2c_d_zeta)/G0) &
          - 3*Y1s*X3c3 -(-3*Y1c)*X3s3
     do j = 1, N_phi
        matrix(index_XYEquation_c+j-1,index_X3c1+j-1) = Y1s(j)
        !matrix(index_XYEquation_c+j-1,index_X3c3+j-1) = 3*Y1s(j)
        matrix(index_XYEquation_c+j-1,index_X3s1+j-1) = Y1c(j)
        !matrix(index_XYEquation_c+j-1,index_X3s3+j-1) = -3*Y1c(j)
        matrix(index_XYEquation_c+j-1,index_Y3s1+j-1) = -X1c(j)
        matrix(index_XYEquation_c+j-1,index_Y3s3+j-1) = 3*X1c(j)
     end do

     ! cos(0 theta) part of the XY equation:
     right_hand_side(index_XYEquation_0:(index_XYEquation_0+N_phi-1)) = -((Bbar*abs_G0_over_B0*G2)/(G0*G0) + (Bbar*abs_G0_over_B0*I2*N_helicity)/(G0*G0) - &
          4*X2s*Y2c + 4*X2c*Y2s + (Bbar*abs_G0_over_B0*X20*curvature)/G0 - &
          (abs_G0_over_B0*I2*X1c*X1c*torsion)/(2*G0) - (abs_G0_over_B0*I2*Y1c*Y1c*torsion)/(2*G0) - &
          (abs_G0_over_B0*I2*Y1s*Y1s*torsion)/(2*G0) + (I2*Y1c*d_X1c_d_zeta)/(2*G0) - &
          (I2*X1c*d_Y1c_d_zeta)/(2*G0) - (Bbar*d_Z20_d_zeta)/G0)
     do j = 1, N_phi
        matrix(index_XYEquation_0+j-1,index_X3c1+j-1) = 2*Y1s(j)
        matrix(index_XYEquation_0+j-1,index_X3s1+j-1) = -2*Y1c(j)
        matrix(index_XYEquation_0+j-1,index_Y3s1+j-1) = 2*X1c(j)
     end do

     ! sin(2 theta) part of the mixed-partials equation:
     right_hand_side(index_mixedPartialsEquation_s:(index_mixedPartialsEquation_s+N_phi-1)) = -((-8*iota_N*X20*X2s)/abs_G0_over_B0   + &
          + (16*I2*X2c*Y20)/Bbar - (16*I2*X20*Y2c)/Bbar - &
          (8*iota_N*Y20*Y2s)/abs_G0_over_B0  - (4*iota_N**2*X1c**2*Z2c)/abs_G0_over_B0**2 - &
          (4*iota_N**2*Y1c**2*Z2c)/abs_G0_over_B0**2 + (8*I2*iota_N*X1c*Y1s*Z2c)/(Bbar*abs_G0_over_B0) - &
          (4*iota_N**2*Y1s**2*Z2c)/abs_G0_over_B0**2 + 2*X20*Y1c*beta_1s - 2*X2c*Y1c*beta_1s - &
          2*X2s*Y1s*beta_1s - 2*X1c*Y20*beta_1s + 2*X1c*Y2c*beta_1s + 2*X1c*Y1s*beta_2c - &
          (4*iota_N*X1c**2*X2s*curvature)/abs_G0_over_B0 - (2*iota_N*X2s*Y1c**2*curvature)/abs_G0_over_B0 + &
          (2*iota_N*X20*Y1c*Y1s*curvature)/abs_G0_over_B0 - (2*iota_N*X2s*Y1s**2*curvature)/abs_G0_over_B0 - &
          (2*iota_N*X1c*Y1s*Y20*curvature)/abs_G0_over_B0 + (6*iota_N*X1c*Y1s*Y2c*curvature)/abs_G0_over_B0 - &
          (2*iota_N*X1c*Y1c*Y2s*curvature)/abs_G0_over_B0 + 4*X2c*Z20*curvature - 4*X20*Z2c*curvature + &
          X1c*Z3c1*curvature - 3*X1c*Z3c3*curvature + (iota_N*X1c**2*Y1c*Y1s*curvature**2)/abs_G0_over_B0 + &
          X1c**2*Z20*curvature**2 - X1c**2*Z2c*curvature**2     + (I2*iota_N*X1c**2*Y1c*Y1s*torsion)/(Bbar*abs_G0_over_B0) + &
          (I2*iota_N*Y1c**3*Y1s*torsion)/(Bbar*abs_G0_over_B0) - (4*I2**2*X1c*Y1c*Y1s**2*torsion)/Bbar**2 + &
          (I2*iota_N*Y1c*Y1s**3*torsion)/(Bbar*abs_G0_over_B0) - 8*X2c*Y20*torsion + 8*X20*Y2c*torsion  - (2*I2*X1c**2*Z20*torsion)/Bbar - &
          (2*I2*Y1c**2*Z20*torsion)/Bbar + (2*I2*Y1s**2*Z20*torsion)/Bbar + &
          (2*I2*X1c**2*Z2c*torsion)/Bbar + (2*I2*Y1c**2*Z2c*torsion)/Bbar - &
          (8*iota_N*X1c*Y1s*Z2c*torsion)/abs_G0_over_B0 + (2*I2*Y1s**2*Z2c*torsion)/Bbar + &
          3*X1c*X20*Y1c*curvature*torsion - 3*X1c*X2c*Y1c*curvature*torsion - &
          6*X1c*X2s*Y1s*curvature*torsion - 3*X1c**2*Y20*curvature*torsion + 3*X1c**2*Y2c*curvature*torsion + &
          (4*I2*X1c*Y1c*Y1s**2*torsion**2)/Bbar   - &
          (I2*iota_N*X1c**2*Y1s*d_X1c_d_zeta)/(Bbar*abs_G0_over_B0**2) + &
          (2*I2**2*X1c*Y1s**2*d_X1c_d_zeta)/(Bbar**2*abs_G0_over_B0) - &
          (I2*iota_N*Y1s**3*d_X1c_d_zeta)/(Bbar*abs_G0_over_B0**2) + &
          (2*I2*Y1c*Z20*d_X1c_d_zeta)/(Bbar*abs_G0_over_B0) - &
          (2*I2*Y1c*Z2c*d_X1c_d_zeta)/(Bbar*abs_G0_over_B0) + &
          (2*iota_N*X1c*Z2s*d_X1c_d_zeta)/abs_G0_over_B0**2 - &
          (4*I2*Y1s*Z2s*d_X1c_d_zeta)/(Bbar*abs_G0_over_B0) - &
          (X1c*X20*curvature*d_X1c_d_zeta)/abs_G0_over_B0 + (X1c*X2c*curvature*d_X1c_d_zeta)/abs_G0_over_B0 - &
          (3*I2*X1c*Y1s**2*torsion*d_X1c_d_zeta)/(Bbar*abs_G0_over_B0) + &
          (4*X2c*d_X20_d_zeta)/abs_G0_over_B0 + (X1c**2*curvature*d_X20_d_zeta)/abs_G0_over_B0 - &
          (4*X20*d_X2c_d_zeta)/abs_G0_over_B0 - (X1c**2*curvature*d_X2c_d_zeta)/abs_G0_over_B0   - &
          (I2*iota_N*X1c*Y1c*Y1s*d_Y1c_d_zeta)/(Bbar*abs_G0_over_B0**2)  - &
          (2*I2*X1c*Z20*d_Y1c_d_zeta)/(Bbar*abs_G0_over_B0) + &
          (2*I2*X1c*Z2c*d_Y1c_d_zeta)/(Bbar*abs_G0_over_B0) - &
          (4*iota_N*Y1s*Z2c*d_Y1c_d_zeta)/abs_G0_over_B0**2 + &
          (2*iota_N*Y1c*Z2s*d_Y1c_d_zeta)/abs_G0_over_B0**2 + &
          (X20*Y1c*curvature*d_Y1c_d_zeta)/abs_G0_over_B0 - (X2c*Y1c*curvature*d_Y1c_d_zeta)/abs_G0_over_B0 - &
          (2*X2s*Y1s*curvature*d_Y1c_d_zeta)/abs_G0_over_B0 - &
          (2*X1c*Y20*curvature*d_Y1c_d_zeta)/abs_G0_over_B0 + &
          (2*X1c*Y2c*curvature*d_Y1c_d_zeta)/abs_G0_over_B0 + &
          (I2*Y1c*Y1s**2*torsion*d_Y1c_d_zeta)/(Bbar*abs_G0_over_B0) - &
          (I2*Y1s**2*d_X1c_d_zeta*d_Y1c_d_zeta)/(Bbar*abs_G0_over_B0**2) + &
          (I2*iota_N*X1c**3*d_Y1s_d_zeta)/(Bbar*abs_G0_over_B0**2) + &
          (I2*iota_N*X1c*Y1c**2*d_Y1s_d_zeta)/(Bbar*abs_G0_over_B0**2) -&
          (2*I2**2*X1c**2*Y1s*d_Y1s_d_zeta)/(Bbar**2*abs_G0_over_B0) + &
          (I2*iota_N*X1c*Y1s**2*d_Y1s_d_zeta)/(Bbar*abs_G0_over_B0**2)  + &
          (4*iota_N*Y1c*Z2c*d_Y1s_d_zeta)/abs_G0_over_B0**2 - &
          (4*I2*X1c*Z2s*d_Y1s_d_zeta)/(Bbar*abs_G0_over_B0) + &
          (2*iota_N*Y1s*Z2s*d_Y1s_d_zeta)/abs_G0_over_B0**2 + &
          (2*X2s*Y1c*curvature*d_Y1s_d_zeta)/abs_G0_over_B0 - (X20*Y1s*curvature*d_Y1s_d_zeta)/abs_G0_over_B0 - &
          (X2c*Y1s*curvature*d_Y1s_d_zeta)/abs_G0_over_B0 - &
          (X1c**2*Y1s*curvature**2*d_Y1s_d_zeta)/abs_G0_over_B0 + &
          (3*I2*X1c**2*Y1s*torsion*d_Y1s_d_zeta)/(Bbar*abs_G0_over_B0) - &
          (I2*Y1c**2*Y1s*torsion*d_Y1s_d_zeta)/(Bbar*abs_G0_over_B0) + &
          (2*I2*X1c*Y1s*d_Y1c_d_zeta*d_Y1s_d_zeta)/(Bbar*abs_G0_over_B0**2) - &
          (I2*X1c*Y1c*d_Y1s_d_zeta**2)/(Bbar*abs_G0_over_B0**2) + (4*Y2c*d_Y20_d_zeta)/abs_G0_over_B0 + &
          (X1c*Y1c*curvature*d_Y20_d_zeta)/abs_G0_over_B0 - (4*Y20*d_Y2c_d_zeta)/abs_G0_over_B0 - &
          (X1c*Y1c*curvature*d_Y2c_d_zeta)/abs_G0_over_B0 - (2*X1c*Y1s*curvature*d_Y2s_d_zeta)/abs_G0_over_B0    - &
          (2*iota_N*Y1c*Y1s*d_Z20_d_zeta)/abs_G0_over_B0**2 - &
          (X1c*d_X1c_d_zeta*d_Z20_d_zeta)/abs_G0_over_B0**2 - &
          (Y1c*d_Y1c_d_zeta*d_Z20_d_zeta)/abs_G0_over_B0**2 + &
          (Y1s*d_Y1s_d_zeta*d_Z20_d_zeta)/abs_G0_over_B0**2 + &
          (X1c*d_X1c_d_zeta*d_Z2c_d_zeta)/abs_G0_over_B0**2 + &
          (Y1c*d_Y1c_d_zeta*d_Z2c_d_zeta)/abs_G0_over_B0**2 + &
          (Y1s*d_Y1s_d_zeta*d_Z2c_d_zeta)/abs_G0_over_B0**2 + &
          (2*iota_N*X1c**2*d_Z2s_d_zeta)/abs_G0_over_B0**2 + (2*iota_N*Y1c**2*d_Z2s_d_zeta)/abs_G0_over_B0**2 + &
          (2*iota_N*Y1s**2*d_Z2s_d_zeta)/abs_G0_over_B0**2 + (4*X1c*Y1s*torsion*d_Z2s_d_zeta)/abs_G0_over_B0 + &
          (2*Y1s*d_Y1c_d_zeta*d_Z2s_d_zeta)/abs_G0_over_B0**2 - &
          (2*Y1c*d_Y1s_d_zeta*d_Z2s_d_zeta)/abs_G0_over_B0**2) &
          -(- (3*X1c)/abs_G0_over_B0 * d_X3c3_d_zeta + ((12*I2*Y1c)/Bbar - 6*Y1c*torsion + (3*d_X1c_d_zeta)/abs_G0_over_B0) * X3c3 &
          + (- (12*iota_N*X1c)/abs_G0_over_B0 + (12*I2*Y1s)/Bbar - 6*Y1s*torsion) * X3s3)

     matrix(index_mixedPartialsEquation_s:(index_mixedPartialsEquation_s+N_phi-1),index_iota2) = (2*Y1c*Y1s)/abs_G0_over_B0
     !matrix(index_mixedPartialsEquation_s,index_X3c1) = diag(X1c/abs_G0_over_B0) * d_d_zeta + diag(- (4*I2*Y1c)/Bbar + 2*Y1c*torsion - (d_X1c_d_zeta)/abs_G0_over_B0)
     !matrix(index_mixedPartialsEquation_s,index_X3c3) = diag(- (3*X1c)/abs_G0_over_B0) * d_d_zeta + diag((12*I2*Y1c)/Bbar - 6*Y1c*torsion + (3*d_X1c_d_zeta)/abs_G0_over_B0)
     !matrix(index_mixedPartialsEquation_s,index_X3s1) = diag((4*I2*Y1s)/Bbar - 2*Y1s*torsion)
     !matrix(index_mixedPartialsEquation_s,index_X3s3) = diag(- (12*iota_N*X1c)/abs_G0_over_B0 + (12*I2*Y1s)/Bbar - 6*Y1s*torsion)
     !matrix(index_mixedPartialsEquation_s,index_Y3c1) = diag((Y1c)/abs_G0_over_B0) * d_d_zeta + diag((4*I2*X1c)/Bbar - 2*X1c*torsion - (d_Y1c_d_zeta)/abs_G0_over_B0)
     !matrix(index_mixedPartialsEquation_s,index_Y3c3) = diag(- (3*Y1c)/abs_G0_over_B0) * d_d_zeta + diag(- (12*I2*X1c)/Bbar + (12*iota_N*Y1s)/abs_G0_over_B0 + 6*X1c*torsion  + (3*d_Y1c_d_zeta)/abs_G0_over_B0)
     !matrix(index_mixedPartialsEquation_s,index_Y3s1) = diag(- (Y1s)/abs_G0_over_B0) * d_d_zeta + diag((d_Y1s_d_zeta)/abs_G0_over_B0)
     !matrix(index_mixedPartialsEquation_s,index_Y3s3) = diag(- (3*Y1s)/abs_G0_over_B0) * d_d_zeta + diag(- (12*iota_N*Y1c)/abs_G0_over_B0  + (3*d_Y1s_d_zeta)/abs_G0_over_B0)
     !matrix(index_mixedPartialsEquation_s,index_iota2) = (2*Y1c*Y1s)/abs_G0_over_B0
     do j = 1, N_phi
        matrix(index_mixedPartialsEquation_s+j-1,index_X3c1:(index_X3c1+N_phi-1)) = X1c(j)/abs_G0_over_B0 * d_d_zeta(j,:)
        matrix(index_mixedPartialsEquation_s+j-1,index_Y3c1:(index_Y3c1+N_phi-1)) = (Y1c(j))/abs_G0_over_B0 * d_d_zeta(j,:)
        matrix(index_mixedPartialsEquation_s+j-1,index_Y3c3:(index_Y3c3+N_phi-1)) = - (3*Y1c(j))/abs_G0_over_B0 * d_d_zeta(j,:)
        matrix(index_mixedPartialsEquation_s+j-1,index_Y3s1:(index_Y3s1+N_phi-1)) = - (Y1s(j))/abs_G0_over_B0 * d_d_zeta(j,:)
        matrix(index_mixedPartialsEquation_s+j-1,index_Y3s3:(index_Y3s3+N_phi-1)) = - (3*Y1s(j))/abs_G0_over_B0 * d_d_zeta(j,:)
     end do
     do j = 1, N_phi
        ! Since d_d_zeta has 0's on the diagonal, we can over-write the diagonal elements here, rather than add to the diagonal elements.
        matrix(index_mixedPartialsEquation_s+j-1,index_X3c1+j-1) = - (4*I2*Y1c(j))/Bbar + 2*Y1c(j)*torsion(j) - (d_X1c_d_zeta(j))/abs_G0_over_B0
        matrix(index_mixedPartialsEquation_s+j-1,index_X3s1+j-1) = (4*I2*Y1s(j))/Bbar - 2*Y1s(j)*torsion(j)
        matrix(index_mixedPartialsEquation_s+j-1,index_Y3c1+j-1) = (4*I2*X1c(j))/Bbar - 2*X1c(j)*torsion(j) - (d_Y1c_d_zeta(j))/abs_G0_over_B0
        matrix(index_mixedPartialsEquation_s+j-1,index_Y3c3+j-1) = - (12*I2*X1c(j))/Bbar + (12*iota_N*Y1s(j))/abs_G0_over_B0 + 6*X1c(j)*torsion(j)  + (3*d_Y1c_d_zeta(j))/abs_G0_over_B0
        matrix(index_mixedPartialsEquation_s+j-1,index_Y3s1+j-1) = (d_Y1s_d_zeta(j))/abs_G0_over_B0
        matrix(index_mixedPartialsEquation_s+j-1,index_Y3s3+j-1) = - (12*iota_N*Y1c(j))/abs_G0_over_B0  + (3*d_Y1s_d_zeta(j))/abs_G0_over_B0
     end do

     ! cos(2 theta) part of the mixed-partials equation:
     right_hand_side(index_mixedPartialsEquation_c:(index_mixedPartialsEquation_c+N_phi-1)) = -((-8*iota_N*X20*X2c)/abs_G0_over_B0      - (16*I2*X2s*Y20)/Bbar - &
          (8*iota_N*Y20*Y2c)/abs_G0_over_B0 + (16*I2*X20*Y2s)/Bbar    + &
          (4*iota_N**2*X1c**2*Z2s)/abs_G0_over_B0**2 + (4*iota_N**2*Y1c**2*Z2s)/abs_G0_over_B0**2 - &
          (8*I2*iota_N*X1c*Y1s*Z2s)/(Bbar*abs_G0_over_B0) + (4*iota_N**2*Y1s**2*Z2s)/abs_G0_over_B0**2 + &
          2*X2s*Y1c*beta_1s - 2*X20*Y1s*beta_1s - 2*X2c*Y1s*beta_1s - 2*X1c*Y2s*beta_1s - &
          2*X1c*Y1s*beta_2s - (iota_N*X1c**2*X20*curvature)/abs_G0_over_B0 - (4*iota_N*X1c**2*X2c*curvature)/abs_G0_over_B0 + &
          (iota_N*X20*Y1c**2*curvature)/abs_G0_over_B0 - (2*iota_N*X2c*Y1c**2*curvature)/abs_G0_over_B0 - &
          (iota_N*X20*Y1s**2*curvature)/abs_G0_over_B0 - (2*iota_N*X2c*Y1s**2*curvature)/abs_G0_over_B0 - &
          (2*iota_N*X1c*Y1c*Y20*curvature)/abs_G0_over_B0 - (2*iota_N*X1c*Y1c*Y2c*curvature)/abs_G0_over_B0 - &
          (6*iota_N*X1c*Y1s*Y2s*curvature)/abs_G0_over_B0 - 4*X2s*Z20*curvature + 4*X20*Z2s*curvature - &
          X1c*Z3s1*curvature + 3*X1c*Z3s3*curvature - (iota_N*X1c**4*curvature**2)/(2*abs_G0_over_B0) - &
          (iota_N*X1c**2*Y1c**2*curvature**2)/(2*abs_G0_over_B0) - (3*iota_N*X1c**2*Y1s**2*curvature**2)/(2*abs_G0_over_B0) + &
          X1c**2*Z2s*curvature**2 + (I2*iota_N*X1c**4*torsion)/(2*Bbar*abs_G0_over_B0)   + (I2*iota_N*X1c**2*Y1c**2*torsion)/(Bbar*abs_G0_over_B0) + &
          (I2*iota_N*Y1c**4*torsion)/(2*Bbar*abs_G0_over_B0) - (2*I2**2*X1c**3*Y1s*torsion)/Bbar**2  - (2*I2**2*X1c*Y1c**2*Y1s*torsion)/Bbar**2 + &
          (2*I2**2*X1c*Y1s**3*torsion)/Bbar**2 - (I2*iota_N*Y1s**4*torsion)/(2*Bbar*abs_G0_over_B0) + &
          8*X2s*Y20*torsion - 8*X20*Y2s*torsion   + &
          (4*I2*Y1c*Y1s*Z20*torsion)/Bbar - (2*I2*X1c**2*Z2s*torsion)/Bbar - &
          (2*I2*Y1c**2*Z2s*torsion)/Bbar + (8*iota_N*X1c*Y1s*Z2s*torsion)/abs_G0_over_B0 - &
          (2*I2*Y1s**2*Z2s*torsion)/Bbar + 3*X1c*X2s*Y1c*curvature*torsion - &
          3*X1c*X20*Y1s*curvature*torsion - 6*X1c*X2c*Y1s*curvature*torsion - 3*X1c**2*Y2s*curvature*torsion - &
          2*X1c**3*Y1s*curvature**2*torsion + (2*I2*X1c**3*Y1s*torsion**2)/Bbar + &
          (2*I2*X1c*Y1c**2*Y1s*torsion**2)/Bbar - (2*I2*X1c*Y1s**3*torsion**2)/Bbar   - &
          (I2*iota_N*X1c**2*Y1c*d_X1c_d_zeta)/(2*Bbar*abs_G0_over_B0**2) - &
          (I2*iota_N*Y1c**3*d_X1c_d_zeta)/(2*Bbar*abs_G0_over_B0**2) + &
          (2*I2**2*X1c*Y1c*Y1s*d_X1c_d_zeta)/(Bbar**2*abs_G0_over_B0) - &
          (3*I2*iota_N*Y1c*Y1s**2*d_X1c_d_zeta)/(2*Bbar*abs_G0_over_B0**2) - &
          (2*I2*Y1s*Z20*d_X1c_d_zeta)/(Bbar*abs_G0_over_B0) + &
          (2*iota_N*X1c*Z2c*d_X1c_d_zeta)/abs_G0_over_B0**2 - &
          (4*I2*Y1s*Z2c*d_X1c_d_zeta)/(Bbar*abs_G0_over_B0) + &
          (2*I2*Y1c*Z2s*d_X1c_d_zeta)/(Bbar*abs_G0_over_B0) - &
          (X1c*X2s*curvature*d_X1c_d_zeta)/abs_G0_over_B0 - &
          (3*I2*X1c*Y1c*Y1s*torsion*d_X1c_d_zeta)/(Bbar*abs_G0_over_B0) + &
          (I2*X1c*Y1s*d_X1c_d_zeta**2)/(2*Bbar*abs_G0_over_B0**2) - &
          (4*X2s*d_X20_d_zeta)/abs_G0_over_B0 + (4*X20*d_X2s_d_zeta)/abs_G0_over_B0 + &
          (X1c**2*curvature*d_X2s_d_zeta)/abs_G0_over_B0   + (I2*iota_N*X1c**3*d_Y1c_d_zeta)/(2*Bbar*abs_G0_over_B0**2) + &
          (I2*iota_N*X1c*Y1c**2*d_Y1c_d_zeta)/(2*Bbar*abs_G0_over_B0**2) - &
          (2*I2**2*X1c**2*Y1s*d_Y1c_d_zeta)/(Bbar**2*abs_G0_over_B0) + &
          (3*I2*iota_N*X1c*Y1s**2*d_Y1c_d_zeta)/(2*Bbar*abs_G0_over_B0**2)   + &
          (2*iota_N*Y1c*Z2c*d_Y1c_d_zeta)/abs_G0_over_B0**2 - &
          (2*I2*X1c*Z2s*d_Y1c_d_zeta)/(Bbar*abs_G0_over_B0) + &
          (4*iota_N*Y1s*Z2s*d_Y1c_d_zeta)/abs_G0_over_B0**2 + &
          (X2s*Y1c*curvature*d_Y1c_d_zeta)/abs_G0_over_B0 - (X20*Y1s*curvature*d_Y1c_d_zeta)/abs_G0_over_B0 - &
          (2*X2c*Y1s*curvature*d_Y1c_d_zeta)/abs_G0_over_B0 - &
          (2*X1c*Y2s*curvature*d_Y1c_d_zeta)/abs_G0_over_B0 - &
          (3*X1c**2*Y1s*curvature**2*d_Y1c_d_zeta)/(2*abs_G0_over_B0) + &
          (7*I2*X1c**2*Y1s*torsion*d_Y1c_d_zeta)/(2*Bbar*abs_G0_over_B0) + &
          (I2*Y1c**2*Y1s*torsion*d_Y1c_d_zeta)/(2*Bbar*abs_G0_over_B0) - &
          (I2*Y1s**3*torsion*d_Y1c_d_zeta)/(2*Bbar*abs_G0_over_B0) - &
          (I2*Y1c*Y1s*d_X1c_d_zeta*d_Y1c_d_zeta)/(Bbar*abs_G0_over_B0**2) + &
          (3*I2*X1c*Y1s*d_Y1c_d_zeta**2)/(2*Bbar*abs_G0_over_B0**2)   + &
          (2*I2*X1c*Z20*d_Y1s_d_zeta)/(Bbar*abs_G0_over_B0) - &
          (4*I2*X1c*Z2c*d_Y1s_d_zeta)/(Bbar*abs_G0_over_B0) + &
          (2*iota_N*Y1s*Z2c*d_Y1s_d_zeta)/abs_G0_over_B0**2 - &
          (4*iota_N*Y1c*Z2s*d_Y1s_d_zeta)/abs_G0_over_B0**2 - &
          (X20*Y1c*curvature*d_Y1s_d_zeta)/abs_G0_over_B0 + (2*X2c*Y1c*curvature*d_Y1s_d_zeta)/abs_G0_over_B0 + &
          (X2s*Y1s*curvature*d_Y1s_d_zeta)/abs_G0_over_B0 + (2*X1c*Y20*curvature*d_Y1s_d_zeta)/abs_G0_over_B0 + &
          (X1c**2*Y1c*curvature**2*d_Y1s_d_zeta)/(2*abs_G0_over_B0) - &
          (I2*X1c**2*Y1c*torsion*d_Y1s_d_zeta)/(2*Bbar*abs_G0_over_B0) - &
          (I2*Y1c**3*torsion*d_Y1s_d_zeta)/(2*Bbar*abs_G0_over_B0) + &
          (I2*Y1c*Y1s**2*torsion*d_Y1s_d_zeta)/(2*Bbar*abs_G0_over_B0) - &
          (I2*X1c**2*d_X1c_d_zeta*d_Y1s_d_zeta)/(2*Bbar*abs_G0_over_B0**2) + &
          (I2*Y1c**2*d_X1c_d_zeta*d_Y1s_d_zeta)/(2*Bbar*abs_G0_over_B0**2) + &
          (I2*Y1s**2*d_X1c_d_zeta*d_Y1s_d_zeta)/(2*Bbar*abs_G0_over_B0**2) - &
          (I2*X1c*Y1c*d_Y1c_d_zeta*d_Y1s_d_zeta)/(Bbar*abs_G0_over_B0**2) - &
          (I2*X1c*Y1s*d_Y1s_d_zeta**2)/(2*Bbar*abs_G0_over_B0**2) - &
          (4*Y2s*d_Y20_d_zeta)/abs_G0_over_B0 - (X1c*Y1s*curvature*d_Y20_d_zeta)/abs_G0_over_B0 - &
          (2*X1c*Y1s*curvature*d_Y2c_d_zeta)/abs_G0_over_B0 + (4*Y20*d_Y2s_d_zeta)/abs_G0_over_B0 + &
          (X1c*Y1c*curvature*d_Y2s_d_zeta)/abs_G0_over_B0     - (iota_N*X1c**2*d_Z20_d_zeta)/abs_G0_over_B0**2 - &
          (iota_N*Y1c**2*d_Z20_d_zeta)/abs_G0_over_B0**2 + (iota_N*Y1s**2*d_Z20_d_zeta)/abs_G0_over_B0**2 + &
          (Y1s*d_Y1c_d_zeta*d_Z20_d_zeta)/abs_G0_over_B0**2 + &
          (Y1c*d_Y1s_d_zeta*d_Z20_d_zeta)/abs_G0_over_B0**2 + &
          (2*iota_N*X1c**2*d_Z2c_d_zeta)/abs_G0_over_B0**2 + (2*iota_N*Y1c**2*d_Z2c_d_zeta)/abs_G0_over_B0**2 + &
          (2*iota_N*Y1s**2*d_Z2c_d_zeta)/abs_G0_over_B0**2 + (4*X1c*Y1s*torsion*d_Z2c_d_zeta)/abs_G0_over_B0 + &
          (2*Y1s*d_Y1c_d_zeta*d_Z2c_d_zeta)/abs_G0_over_B0**2 - &
          (2*Y1c*d_Y1s_d_zeta*d_Z2c_d_zeta)/abs_G0_over_B0**2 - &
          (X1c*d_X1c_d_zeta*d_Z2s_d_zeta)/abs_G0_over_B0**2 - &
          (Y1c*d_Y1c_d_zeta*d_Z2s_d_zeta)/abs_G0_over_B0**2 - &
          (Y1s*d_Y1s_d_zeta*d_Z2s_d_zeta)/abs_G0_over_B0**2) &
          -(   (- (12*iota_N*X1c)/abs_G0_over_B0 + (12*I2*Y1s)/Bbar  - 6*Y1s*torsion) * X3c3 &
          + (3*X1c)/abs_G0_over_B0 * d_X3s3_d_zeta + (-  (12*I2*Y1c)/Bbar + 6*Y1c*torsion - (3*d_X1c_d_zeta)/abs_G0_over_B0) * X3s3)

     matrix(index_mixedPartialsEquation_c:(index_mixedPartialsEquation_c+N_phi-1),index_iota2) = (X1c**2/abs_G0_over_B0 + Y1c**2/abs_G0_over_B0 - Y1s**2/abs_G0_over_B0)
     do j = 1, N_phi
        matrix(index_mixedPartialsEquation_c+j-1,index_X3s1:(index_X3s1+N_phi-1)) = - (X1c(j))/abs_G0_over_B0 * d_d_zeta(j,:)
        matrix(index_mixedPartialsEquation_c+j-1,index_Y3c1:(index_Y3c1+N_phi-1)) = - (Y1s(j))/abs_G0_over_B0 * d_d_zeta(j,:)
        matrix(index_mixedPartialsEquation_c+j-1,index_Y3c3:(index_Y3c3+N_phi-1)) = - (3*Y1s(j))/abs_G0_over_B0 * d_d_zeta(j,:)
        matrix(index_mixedPartialsEquation_c+j-1,index_Y3s1:(index_Y3s1+N_phi-1)) = - (Y1c(j))/abs_G0_over_B0 * d_d_zeta(j,:)
        matrix(index_mixedPartialsEquation_c+j-1,index_Y3s3:(index_Y3s3+N_phi-1)) = (3*Y1c(j))/abs_G0_over_B0 * d_d_zeta(j,:)
     end do
     do j = 1, N_phi
        ! Since d_d_zeta has 0's on the diagonal, we can over-write the diagonal elements here, rather than add to the diagonal elements.
        matrix(index_mixedPartialsEquation_c+j-1,index_X3c1+j-1) = (4*I2*Y1s(j))/Bbar - 2*Y1s(j)*torsion(j)
        matrix(index_mixedPartialsEquation_c+j-1,index_X3s1+j-1) = (4*I2*Y1c(j))/Bbar - 2*Y1c(j)*torsion(j) + d_X1c_d_zeta(j)/abs_G0_over_B0
        matrix(index_mixedPartialsEquation_c+j-1,index_Y3c1+j-1) = (d_Y1s_d_zeta(j))/abs_G0_over_B0
        matrix(index_mixedPartialsEquation_c+j-1,index_Y3c3+j-1) = - (12*iota_N*Y1c(j))/abs_G0_over_B0 + (3*d_Y1s_d_zeta(j))/abs_G0_over_B0
        matrix(index_mixedPartialsEquation_c+j-1,index_Y3s1+j-1) = - (4*I2*X1c(j))/Bbar + 2*X1c(j)*torsion(j) + (d_Y1c_d_zeta(j))/abs_G0_over_B0
        matrix(index_mixedPartialsEquation_c+j-1,index_Y3s3+j-1) = (12*I2*X1c(j))/Bbar - (12*iota_N*Y1s(j))/abs_G0_over_B0 - 6*X1c(j)*torsion(j) - (3*d_Y1c_d_zeta(j))/abs_G0_over_B0
     end do

     !matrix(index_mixedPartialsEquation_c,index_X3c1) = diag((4*I2*Y1s)/Bbar - 2*Y1s*torsion)
     !matrix(index_mixedPartialsEquation_c,index_X3c3) = diag(- (12*iota_N*X1c)/abs_G0_over_B0 + (12*I2*Y1s)/Bbar  - 6*Y1s*torsion)
     !matrix(index_mixedPartialsEquation_c,index_X3s1) = diag(- (X1c)/abs_G0_over_B0) * d_d_zeta + diag((4*I2*Y1c)/Bbar - 2*Y1c*torsion + (d_X1c_d_zeta)/abs_G0_over_B0)
     !matrix(index_mixedPartialsEquation_c,index_X3s3) = diag((3*X1c)/abs_G0_over_B0) * d_d_zeta + diag(-  (12*I2*Y1c)/Bbar + 6*Y1c*torsion - (3*d_X1c_d_zeta)/abs_G0_over_B0)
     !matrix(index_mixedPartialsEquation_c,index_Y3c1) = diag(- (Y1s)/abs_G0_over_B0) * d_d_zeta + diag((d_Y1s_d_zeta)/abs_G0_over_B0)
     !matrix(index_mixedPartialsEquation_c,index_Y3c3) = diag(- (3*Y1s)/abs_G0_over_B0) * d_d_zeta + diag(- (12*iota_N*Y1c)/abs_G0_over_B0 + (3*d_Y1s_d_zeta)/abs_G0_over_B0)
     !matrix(index_mixedPartialsEquation_c,index_Y3s1) = diag(- (Y1c)/abs_G0_over_B0) * d_d_zeta + diag(- (4*I2*X1c)/Bbar + 2*X1c*torsion + (d_Y1c_d_zeta)/abs_G0_over_B0)
     !matrix(index_mixedPartialsEquation_c,index_Y3s3) = diag((3*Y1c)/abs_G0_over_B0) * d_d_zeta + diag((12*I2*X1c)/Bbar - (12*iota_N*Y1s)/abs_G0_over_B0 - 6*X1c*torsion - (3*d_Y1c_d_zeta)/abs_G0_over_B0)
     !matrix(index_mixedPartialsEquation_c,index_iota2) = (X1c**2/abs_G0_over_B0 + Y1c**2/abs_G0_over_B0 - Y1s**2/abs_G0_over_B0)
     
     ! cos(0 theta) part of the mixed-partials equation:
     right_hand_side(index_mixedPartialsEquation_0:(index_mixedPartialsEquation_0+N_phi-1)) = -((-8*iota_N*X2c**2)/abs_G0_over_B0 - (8*iota_N*X2s**2)/abs_G0_over_B0   + (4*I4*X1c*Y1s)/Bbar    - (16*I2*X2s*Y2c)/Bbar - &
          (8*iota_N*Y2c**2)/abs_G0_over_B0 + (16*I2*X2c*Y2s)/Bbar - (8*iota_N*Y2s**2)/abs_G0_over_B0    + &
          (4*iota_N**2*Y1c*Y1s*Z2c)/abs_G0_over_B0**2 - (2*iota_N**2*X1c**2*Z2s)/abs_G0_over_B0**2 - &
          (2*iota_N**2*Y1c**2*Z2s)/abs_G0_over_B0**2 + (2*iota_N**2*Y1s**2*Z2s)/abs_G0_over_B0**2 - &
          (2*iota_N*X1c**2*X20*curvature)/abs_G0_over_B0 - (3*iota_N*X1c**2*X2c*curvature)/abs_G0_over_B0 - &
          (2*iota_N*X20*Y1c**2*curvature)/abs_G0_over_B0 + (iota_N*X2c*Y1c**2*curvature)/abs_G0_over_B0 + &
          (2*iota_N*X2s*Y1c*Y1s*curvature)/abs_G0_over_B0 - (2*iota_N*X20*Y1s**2*curvature)/abs_G0_over_B0 - &
          (iota_N*X2c*Y1s**2*curvature)/abs_G0_over_B0 - (4*iota_N*X1c*Y1c*Y2c*curvature)/abs_G0_over_B0 - &
          (4*iota_N*X1c*Y1s*Y2s*curvature)/abs_G0_over_B0 - 4*X2s*Z2c*curvature + 4*X2c*Z2s*curvature + &
          2*X1c*Z3s1*curvature - (iota_N*X1c**4*curvature**2)/(2*abs_G0_over_B0) - &
          (iota_N*X1c**2*Y1c**2*curvature**2)/(2*abs_G0_over_B0) - (3*iota_N*X1c**2*Y1s**2*curvature**2)/(2*abs_G0_over_B0) + &
          X1c**2*Z2s*curvature**2 + (I2*iota_N*X1c**4*torsion)/(2*Bbar*abs_G0_over_B0)  + &
          (I2*iota_N*X1c**2*Y1c**2*torsion)/(Bbar*abs_G0_over_B0) + (I2*iota_N*Y1c**4*torsion)/(2*Bbar*abs_G0_over_B0) - &
          (2*I2**2*X1c**3*Y1s*torsion)/Bbar**2  - &
          (2*I2**2*X1c*Y1c**2*Y1s*torsion)/Bbar**2 + (3*I2*iota_N*X1c**2*Y1s**2*torsion)/(Bbar*abs_G0_over_B0) + &
          (I2*iota_N*Y1c**2*Y1s**2*torsion)/(Bbar*abs_G0_over_B0) - (2*I2**2*X1c*Y1s**3*torsion)/Bbar**2 + &
          (I2*iota_N*Y1s**4*torsion)/(2*Bbar*abs_G0_over_B0) + 8*X2s*Y2c*torsion - 8*X2c*Y2s*torsion  + (4*I2*Y1c*Y1s*Z2c*torsion)/Bbar - (2*I2*X1c**2*Z2s*torsion)/Bbar - &
          (2*I2*Y1c**2*Z2s*torsion)/Bbar + (2*I2*Y1s**2*Z2s*torsion)/Bbar + &
          3*X1c*X2s*Y1c*curvature*torsion - 6*X1c*X20*Y1s*curvature*torsion - &
          3*X1c*X2c*Y1s*curvature*torsion - 3*X1c**2*Y2s*curvature*torsion - 2*X1c**3*Y1s*curvature**2*torsion + &
          (2*I2*X1c**3*Y1s*torsion**2)/Bbar + (2*I2*X1c*Y1c**2*Y1s*torsion**2)/Bbar + &
          (2*I2*X1c*Y1s**3*torsion**2)/Bbar  - &
          (I2*iota_N*X1c**2*Y1c*d_X1c_d_zeta)/(2*Bbar*abs_G0_over_B0**2) - &
          (I2*iota_N*Y1c**3*d_X1c_d_zeta)/(2*Bbar*abs_G0_over_B0**2) + &
          (2*I2**2*X1c*Y1c*Y1s*d_X1c_d_zeta)/(Bbar**2*abs_G0_over_B0) - &
          (I2*iota_N*Y1c*Y1s**2*d_X1c_d_zeta)/(2*Bbar*abs_G0_over_B0**2) - &
          (4*I2*Y1s*Z20*d_X1c_d_zeta)/(Bbar*abs_G0_over_B0) + &
          (2*iota_N*X1c*Z2c*d_X1c_d_zeta)/abs_G0_over_B0**2 - &
          (2*I2*Y1s*Z2c*d_X1c_d_zeta)/(Bbar*abs_G0_over_B0) + &
          (2*I2*Y1c*Z2s*d_X1c_d_zeta)/(Bbar*abs_G0_over_B0) - &
          (X1c*X2s*curvature*d_X1c_d_zeta)/abs_G0_over_B0 - &
          (3*I2*X1c*Y1c*Y1s*torsion*d_X1c_d_zeta)/(Bbar*abs_G0_over_B0) + &
          (I2*X1c*Y1s*d_X1c_d_zeta**2)/(2*Bbar*abs_G0_over_B0**2) - &
          (4*X2s*d_X2c_d_zeta)/abs_G0_over_B0 + (4*X2c*d_X2s_d_zeta)/abs_G0_over_B0 + &
          (X1c**2*curvature*d_X2s_d_zeta)/abs_G0_over_B0  + &
          (I2*iota_N*X1c**3*d_Y1c_d_zeta)/(2*Bbar*abs_G0_over_B0**2) + &
          (I2*iota_N*X1c*Y1c**2*d_Y1c_d_zeta)/(2*Bbar*abs_G0_over_B0**2) - &
          (2*I2**2*X1c**2*Y1s*d_Y1c_d_zeta)/(Bbar**2*abs_G0_over_B0) + &
          (3*I2*iota_N*X1c*Y1s**2*d_Y1c_d_zeta)/(2*Bbar*abs_G0_over_B0**2)  + (2*iota_N*Y1c*Z2c*d_Y1c_d_zeta)/abs_G0_over_B0**2 - &
          (2*I2*X1c*Z2s*d_Y1c_d_zeta)/(Bbar*abs_G0_over_B0) + &
          (2*iota_N*Y1s*Z2s*d_Y1c_d_zeta)/abs_G0_over_B0**2 + &
          (X2s*Y1c*curvature*d_Y1c_d_zeta)/abs_G0_over_B0 - (2*X20*Y1s*curvature*d_Y1c_d_zeta)/abs_G0_over_B0 - &
          (X2c*Y1s*curvature*d_Y1c_d_zeta)/abs_G0_over_B0 - (2*X1c*Y2s*curvature*d_Y1c_d_zeta)/abs_G0_over_B0 - &
          (3*X1c**2*Y1s*curvature**2*d_Y1c_d_zeta)/(2*abs_G0_over_B0) + &
          (7*I2*X1c**2*Y1s*torsion*d_Y1c_d_zeta)/(2*Bbar*abs_G0_over_B0) + &
          (I2*Y1c**2*Y1s*torsion*d_Y1c_d_zeta)/(2*Bbar*abs_G0_over_B0) + &
          (I2*Y1s**3*torsion*d_Y1c_d_zeta)/(2*Bbar*abs_G0_over_B0) - &
          (I2*Y1c*Y1s*d_X1c_d_zeta*d_Y1c_d_zeta)/(Bbar*abs_G0_over_B0**2) + &
          (3*I2*X1c*Y1s*d_Y1c_d_zeta**2)/(2*Bbar*abs_G0_over_B0**2) - &
          (I2*iota_N*X1c*Y1c*Y1s*d_Y1s_d_zeta)/(Bbar*abs_G0_over_B0**2)  - (4*I2*X1c*Z20*d_Y1s_d_zeta)/(Bbar*abs_G0_over_B0) + &
          (2*I2*X1c*Z2c*d_Y1s_d_zeta)/(Bbar*abs_G0_over_B0) - &
          (2*iota_N*Y1s*Z2c*d_Y1s_d_zeta)/abs_G0_over_B0**2 + &
          (2*iota_N*Y1c*Z2s*d_Y1s_d_zeta)/abs_G0_over_B0**2 + &
          (2*X20*Y1c*curvature*d_Y1s_d_zeta)/abs_G0_over_B0 - (X2c*Y1c*curvature*d_Y1s_d_zeta)/abs_G0_over_B0 - &
          (X2s*Y1s*curvature*d_Y1s_d_zeta)/abs_G0_over_B0 + (2*X1c*Y2c*curvature*d_Y1s_d_zeta)/abs_G0_over_B0 + &
          (X1c**2*Y1c*curvature**2*d_Y1s_d_zeta)/(2*abs_G0_over_B0) - &
          (I2*X1c**2*Y1c*torsion*d_Y1s_d_zeta)/(2*Bbar*abs_G0_over_B0) - &
          (I2*Y1c**3*torsion*d_Y1s_d_zeta)/(2*Bbar*abs_G0_over_B0) - &
          (I2*Y1c*Y1s**2*torsion*d_Y1s_d_zeta)/(2*Bbar*abs_G0_over_B0) - &
          (I2*X1c**2*d_X1c_d_zeta*d_Y1s_d_zeta)/(2*Bbar*abs_G0_over_B0**2) + &
          (I2*Y1c**2*d_X1c_d_zeta*d_Y1s_d_zeta)/(2*Bbar*abs_G0_over_B0**2) - &
          (I2*Y1s**2*d_X1c_d_zeta*d_Y1s_d_zeta)/(2*Bbar*abs_G0_over_B0**2) - &
          (I2*X1c*Y1c*d_Y1c_d_zeta*d_Y1s_d_zeta)/(Bbar*abs_G0_over_B0**2) + &
          (I2*X1c*Y1s*d_Y1s_d_zeta**2)/(2*Bbar*abs_G0_over_B0**2) - &
          (2*X1c*Y1s*curvature*d_Y20_d_zeta)/abs_G0_over_B0 - (4*Y2s*d_Y2c_d_zeta)/abs_G0_over_B0 - &
          (X1c*Y1s*curvature*d_Y2c_d_zeta)/abs_G0_over_B0 + (4*Y2c*d_Y2s_d_zeta)/abs_G0_over_B0 + &
          (X1c*Y1c*curvature*d_Y2s_d_zeta)/abs_G0_over_B0   + (2*iota_N*X1c**2*d_Z20_d_zeta)/abs_G0_over_B0**2 + &
          (2*iota_N*Y1c**2*d_Z20_d_zeta)/abs_G0_over_B0**2 + (2*iota_N*Y1s**2*d_Z20_d_zeta)/abs_G0_over_B0**2 + &
          (4*X1c*Y1s*torsion*d_Z20_d_zeta)/abs_G0_over_B0 + &
          (2*Y1s*d_Y1c_d_zeta*d_Z20_d_zeta)/abs_G0_over_B0**2 - &
          (2*Y1c*d_Y1s_d_zeta*d_Z20_d_zeta)/abs_G0_over_B0**2 - &
          (iota_N*X1c**2*d_Z2c_d_zeta)/abs_G0_over_B0**2 - (iota_N*Y1c**2*d_Z2c_d_zeta)/abs_G0_over_B0**2 + &
          (iota_N*Y1s**2*d_Z2c_d_zeta)/abs_G0_over_B0**2 + &
          (Y1s*d_Y1c_d_zeta*d_Z2c_d_zeta)/abs_G0_over_B0**2 + &
          (Y1c*d_Y1s_d_zeta*d_Z2c_d_zeta)/abs_G0_over_B0**2 - &
          (2*iota_N*Y1c*Y1s*d_Z2s_d_zeta)/abs_G0_over_B0**2 - &
          (X1c*d_X1c_d_zeta*d_Z2s_d_zeta)/abs_G0_over_B0**2 - &
          (Y1c*d_Y1c_d_zeta*d_Z2s_d_zeta)/abs_G0_over_B0**2 + &
          (Y1s*d_Y1s_d_zeta*d_Z2s_d_zeta)/abs_G0_over_B0**2)

     matrix(index_mixedPartialsEquation_0:(index_mixedPartialsEquation_0+N_phi-1),index_iota2) = ((-2*X1c**2)/abs_G0_over_B0 - (2*Y1c**2)/abs_G0_over_B0 - (2*Y1s**2)/abs_G0_over_B0)
     do j = 1, N_phi
        matrix(index_mixedPartialsEquation_0+j-1,index_X3s1:(index_X3s1+N_phi-1)) = (2*X1c(j))/abs_G0_over_B0 * d_d_zeta(j,:)
        matrix(index_mixedPartialsEquation_0+j-1,index_Y3c1:(index_Y3c1+N_phi-1)) = - (2*Y1s(j))/abs_G0_over_B0 * d_d_zeta(j,:)
        matrix(index_mixedPartialsEquation_0+j-1,index_Y3s1:(index_Y3s1+N_phi-1)) = (2*Y1c(j))/abs_G0_over_B0 * d_d_zeta(j,:)
     end do
     do j = 1, N_phi
        matrix(index_mixedPartialsEquation_0+j-1,index_X3c1+j-1) = - (4*iota_N*X1c(j))/abs_G0_over_B0 - 4*Y1s(j)*torsion(j) + (8*I2*Y1s(j))/Bbar
        matrix(index_mixedPartialsEquation_0+j-1,index_X3s1+j-1) = - (8*I2*Y1c(j))/Bbar + 4*Y1c(j)*torsion(j) - (2*d_X1c_d_zeta(j))/abs_G0_over_B0
        matrix(index_mixedPartialsEquation_0+j-1,index_Y3c1+j-1) = - (4*iota_N*Y1c(j))/abs_G0_over_B0 + (2*d_Y1s_d_zeta(j))/abs_G0_over_B0
        matrix(index_mixedPartialsEquation_0+j-1,index_Y3s1+j-1) = (8*I2*X1c(j))/Bbar - (4*iota_N*Y1s(j))/abs_G0_over_B0 - 4*X1c(j)*torsion(j) - (2*d_Y1c_d_zeta(j))/abs_G0_over_B0
     end do
!!$     matrix(index_mixedPartialsEquation_0,index_X3c1) = diag(- (4*iota_N*X1c)/abs_G0_over_B0 - 4*Y1s*torsion + (8*I2*Y1s)/Bbar)
!!$     matrix(index_mixedPartialsEquation_0,index_X3s1) = diag((2*X1c)/abs_G0_over_B0) * d_d_zeta + diag(- (8*I2*Y1c)/Bbar + 4*Y1c*torsion - (2*d_X1c_d_zeta)/abs_G0_over_B0)
!!$     matrix(index_mixedPartialsEquation_0,index_Y3c1) = diag(- (2*Y1s)/abs_G0_over_B0) * d_d_zeta + diag(- (4*iota_N*Y1c)/abs_G0_over_B0 + (2*d_Y1s_d_zeta)/abs_G0_over_B0)
!!$     matrix(index_mixedPartialsEquation_0,index_Y3s1) = diag(+ (2*Y1c)/abs_G0_over_B0) * d_d_zeta + diag( (8*I2*X1c)/Bbar - (4*iota_N*Y1s)/abs_G0_over_B0 - 4*X1c*torsion - (2*d_Y1c_d_zeta)/abs_G0_over_B0)
!!$     matrix(index_mixedPartialsEquation_0,index_iota2) = ((-2*X1c**2)/abs_G0_over_B0 - (2*Y1c**2)/abs_G0_over_B0 - (2*Y1s**2)/abs_G0_over_B0)
     
     
     ! We will use the LAPACK subroutine DGESV to solve a general (asymmetric) linear system
     ! solution = matrix \ right_hand_side
     ! Note that LAPACK will over-write "right_hand_side" with the solution, and over-write "matrix" with the LU factorization.
     print *,"Solving the FULL O(r^3) problem!"
     allocate(IPIV(vector_size))
     call DGESV(vector_size, 1, matrix, vector_size, IPIV, right_hand_side, vector_size, INFO)
     deallocate(IPIV)
     if (INFO /= 0) then
        print *, "Error in LAPACK call DGESV: info = ", INFO
        stop
     end if

     ! Extract the parts from the solution vector:
     X3s1 = right_hand_side(index_X3s1:(index_X3s1+N_phi-1))
     X3c1 = right_hand_side(index_X3c1:(index_X3c1+N_phi-1))
     Y3s1 = right_hand_side(index_Y3s1:(index_Y3s1+N_phi-1))
     Y3c1 = right_hand_side(index_Y3c1:(index_Y3c1+N_phi-1))
     Y3s3 = right_hand_side(index_Y3s3:(index_Y3s3+N_phi-1))
     Y3c3 = right_hand_side(index_Y3c3:(index_Y3c3+N_phi-1))
     iota2 = right_hand_side(index_iota2)

     deallocate(matrix, right_hand_side)

     if (.true.) then
        ! Verify that we have actually solved the equations

        allocate(d_X3c1_d_zeta(N_phi))
        allocate(d_X3s1_d_zeta(N_phi))
        allocate(d_Y3c1_d_zeta(N_phi))
        allocate(d_Y3s1_d_zeta(N_phi))
        allocate(d_Y3c3_d_zeta(N_phi))
        allocate(d_Y3s3_d_zeta(N_phi))
        d_X3c1_d_zeta = matmul(d_d_zeta,X3c1)
        d_X3s1_d_zeta = matmul(d_d_zeta,X3s1)
        d_Y3c1_d_zeta = matmul(d_d_zeta,Y3c1)
        d_Y3s1_d_zeta = matmul(d_d_zeta,Y3s1)
        d_Y3c3_d_zeta = matmul(d_d_zeta,Y3c3)
        d_Y3s3_d_zeta = matmul(d_d_zeta,Y3s3)

        ! sin(2*theta) component of the XY equation:
        eq1residual = -(X3c1*Y1c) + 3*X3c3*Y1c + X3s1*Y1s + 3*X3s3*Y1s + 4*X2c*Y20 - &
             4*X20*Y2c + X1c*Y3c1 - 3*X1c*Y3c3 + (2*Bbar*iota_N*Z2c)/G0 + &
             (Bbar*abs_G0_over_B0*X2s*curvature)/G0 - (abs_G0_over_B0*I2*Y1c*Y1s*torsion)/G0 + &
             (I2*Y1s*d_X1c_d_zeta)/(2*G0) - (I2*X1c*d_Y1s_d_zeta)/(2*G0) - &
             (Bbar*d_Z2s_d_zeta)/G0

        max_eq1residual = maxval(abs(eq1residual))
        print *,"max(abs(residual of sin(2*theta) XY equation)):",max_eq1residual
        if (max_eq1residual > 1e-8) stop "Residual is large !!!"

        ! cos(2*theta) component of the XY equation: 
        eq1residual = X3s1*Y1c - 3*X3s3*Y1c + X3c1*Y1s + 3*X3c3*Y1s - 4*X2s*Y20 + &
             4*X20*Y2s - X1c*Y3s1 + 3*X1c*Y3s3 - (2*Bbar*iota_N*Z2s)/G0 + &
             (Bbar*abs_G0_over_B0*X2c*curvature)/G0 - (abs_G0_over_B0*I2*X1c**2*torsion)/(2*G0) - (abs_G0_over_B0*I2*Y1c**2*torsion)/(2*G0) + &
             (abs_G0_over_B0*I2*Y1s**2*torsion)/(2*G0) + (I2*Y1c*d_X1c_d_zeta)/(2*G0) - &
             (I2*X1c*d_Y1c_d_zeta)/(2*G0) - (Bbar*d_Z2c_d_zeta)/G0

        max_eq1residual = maxval(abs(eq1residual))
        print *,"max(abs(residual of cos(2*theta) XY equation)):",max_eq1residual
        if (max_eq1residual > 1e-8) stop "Residual is large !!!"

        ! cos(0*theta) component of the XY equation: 
        eq1residual = (Bbar*abs_G0_over_B0*G2)/G0**2 + (Bbar*abs_G0_over_B0*I2*N_helicity)/G0**2 - 2*X3s1*Y1c + 2*X3c1*Y1s - &
             4*X2s*Y2c + 4*X2c*Y2s + 2*X1c*Y3s1 + (Bbar*abs_G0_over_B0*X20*curvature)/G0 - &
             (abs_G0_over_B0*I2*X1c**2*torsion)/(2*G0) - (abs_G0_over_B0*I2*Y1c**2*torsion)/(2*G0) - &
             (abs_G0_over_B0*I2*Y1s**2*torsion)/(2*G0) + (I2*Y1c*d_X1c_d_zeta)/(2*G0) - &
             (I2*X1c*d_Y1c_d_zeta)/(2*G0) - (Bbar*d_Z20_d_zeta)/G0

        max_eq1residual = maxval(abs(eq1residual))
        print *,"max(abs(residual of cos(0*theta) XY equation)):",max_eq1residual
        if (max_eq1residual > 1e-8) stop "Residual is large !!!"

        ! sin(2*theta) component of the mixed-partials equation:
        eq1residual = (-8*iota_N*X20*X2s)/abs_G0_over_B0 - (12*iota_N*X1c*X3s3)/abs_G0_over_B0 - (4*I2*X3c1*Y1c)/Bbar + &
             (12*I2*X3c3*Y1c)/Bbar + (4*I2*X3s1*Y1s)/Bbar + (12*I2*X3s3*Y1s)/Bbar + &
             (2*iota2*Y1c*Y1s)/abs_G0_over_B0 + (16*I2*X2c*Y20)/Bbar - (16*I2*X20*Y2c)/Bbar - &
             (8*iota_N*Y20*Y2s)/abs_G0_over_B0 + (4*I2*X1c*Y3c1)/Bbar - (12*I2*X1c*Y3c3)/Bbar + &
             (12*iota_N*Y1s*Y3c3)/abs_G0_over_B0 - (12*iota_N*Y1c*Y3s3)/abs_G0_over_B0 - (4*iota_N**2*X1c**2*Z2c)/abs_G0_over_B0**2 - &
             (4*iota_N**2*Y1c**2*Z2c)/abs_G0_over_B0**2 + (8*I2*iota_N*X1c*Y1s*Z2c)/(Bbar*abs_G0_over_B0) - &
             (4*iota_N**2*Y1s**2*Z2c)/abs_G0_over_B0**2 + 2*X20*Y1c*beta_1s - 2*X2c*Y1c*beta_1s - &
             2*X2s*Y1s*beta_1s - 2*X1c*Y20*beta_1s + 2*X1c*Y2c*beta_1s + 2*X1c*Y1s*beta_2c - &
             (4*iota_N*X1c**2*X2s*curvature)/abs_G0_over_B0 - (2*iota_N*X2s*Y1c**2*curvature)/abs_G0_over_B0 + &
             (2*iota_N*X20*Y1c*Y1s*curvature)/abs_G0_over_B0 - (2*iota_N*X2s*Y1s**2*curvature)/abs_G0_over_B0 - &
             (2*iota_N*X1c*Y1s*Y20*curvature)/abs_G0_over_B0 + (6*iota_N*X1c*Y1s*Y2c*curvature)/abs_G0_over_B0 - &
             (2*iota_N*X1c*Y1c*Y2s*curvature)/abs_G0_over_B0 + 4*X2c*Z20*curvature - 4*X20*Z2c*curvature + &
             X1c*Z3c1*curvature - 3*X1c*Z3c3*curvature + (iota_N*X1c**2*Y1c*Y1s*curvature**2)/abs_G0_over_B0 + &
             X1c**2*Z20*curvature**2 - X1c**2*Z2c*curvature**2 + 2*X3c1*Y1c*torsion - 6*X3c3*Y1c*torsion - &
             2*X3s1*Y1s*torsion - 6*X3s3*Y1s*torsion + (I2*iota_N*X1c**2*Y1c*Y1s*torsion)/(Bbar*abs_G0_over_B0) + &
             (I2*iota_N*Y1c**3*Y1s*torsion)/(Bbar*abs_G0_over_B0) - (4*I2**2*X1c*Y1c*Y1s**2*torsion)/Bbar**2 + &
             (I2*iota_N*Y1c*Y1s**3*torsion)/(Bbar*abs_G0_over_B0) - 8*X2c*Y20*torsion + 8*X20*Y2c*torsion - &
             2*X1c*Y3c1*torsion + 6*X1c*Y3c3*torsion - (2*I2*X1c**2*Z20*torsion)/Bbar - &
             (2*I2*Y1c**2*Z20*torsion)/Bbar + (2*I2*Y1s**2*Z20*torsion)/Bbar + &
             (2*I2*X1c**2*Z2c*torsion)/Bbar + (2*I2*Y1c**2*Z2c*torsion)/Bbar - &
             (8*iota_N*X1c*Y1s*Z2c*torsion)/abs_G0_over_B0 + (2*I2*Y1s**2*Z2c*torsion)/Bbar + &
             3*X1c*X20*Y1c*curvature*torsion - 3*X1c*X2c*Y1c*curvature*torsion - &
             6*X1c*X2s*Y1s*curvature*torsion - 3*X1c**2*Y20*curvature*torsion + 3*X1c**2*Y2c*curvature*torsion + &
             (4*I2*X1c*Y1c*Y1s**2*torsion**2)/Bbar - (X3c1*d_X1c_d_zeta)/abs_G0_over_B0 + &
             (3*X3c3*d_X1c_d_zeta)/abs_G0_over_B0 - &
             (I2*iota_N*X1c**2*Y1s*d_X1c_d_zeta)/(Bbar*abs_G0_over_B0**2) + &
             (2*I2**2*X1c*Y1s**2*d_X1c_d_zeta)/(Bbar**2*abs_G0_over_B0) - &
             (I2*iota_N*Y1s**3*d_X1c_d_zeta)/(Bbar*abs_G0_over_B0**2) + &
             (2*I2*Y1c*Z20*d_X1c_d_zeta)/(Bbar*abs_G0_over_B0) - &
             (2*I2*Y1c*Z2c*d_X1c_d_zeta)/(Bbar*abs_G0_over_B0) + &
             (2*iota_N*X1c*Z2s*d_X1c_d_zeta)/abs_G0_over_B0**2 - &
             (4*I2*Y1s*Z2s*d_X1c_d_zeta)/(Bbar*abs_G0_over_B0) - &
             (X1c*X20*curvature*d_X1c_d_zeta)/abs_G0_over_B0 + (X1c*X2c*curvature*d_X1c_d_zeta)/abs_G0_over_B0 - &
             (3*I2*X1c*Y1s**2*torsion*d_X1c_d_zeta)/(Bbar*abs_G0_over_B0) + &
             (4*X2c*d_X20_d_zeta)/abs_G0_over_B0 + (X1c**2*curvature*d_X20_d_zeta)/abs_G0_over_B0 - &
             (4*X20*d_X2c_d_zeta)/abs_G0_over_B0 - (X1c**2*curvature*d_X2c_d_zeta)/abs_G0_over_B0 + &
             (X1c*d_X3c1_d_zeta)/abs_G0_over_B0 - (3*X1c*d_X3c3_d_zeta)/abs_G0_over_B0 - &
             (I2*iota_N*X1c*Y1c*Y1s*d_Y1c_d_zeta)/(Bbar*abs_G0_over_B0**2) - &
             (Y3c1*d_Y1c_d_zeta)/abs_G0_over_B0 + (3*Y3c3*d_Y1c_d_zeta)/abs_G0_over_B0 - &
             (2*I2*X1c*Z20*d_Y1c_d_zeta)/(Bbar*abs_G0_over_B0) + &
             (2*I2*X1c*Z2c*d_Y1c_d_zeta)/(Bbar*abs_G0_over_B0) - &
             (4*iota_N*Y1s*Z2c*d_Y1c_d_zeta)/abs_G0_over_B0**2 + &
             (2*iota_N*Y1c*Z2s*d_Y1c_d_zeta)/abs_G0_over_B0**2 + &
             (X20*Y1c*curvature*d_Y1c_d_zeta)/abs_G0_over_B0 - (X2c*Y1c*curvature*d_Y1c_d_zeta)/abs_G0_over_B0 - &
             (2*X2s*Y1s*curvature*d_Y1c_d_zeta)/abs_G0_over_B0 - &
             (2*X1c*Y20*curvature*d_Y1c_d_zeta)/abs_G0_over_B0 + &
             (2*X1c*Y2c*curvature*d_Y1c_d_zeta)/abs_G0_over_B0 + &
             (I2*Y1c*Y1s**2*torsion*d_Y1c_d_zeta)/(Bbar*abs_G0_over_B0) - &
             (I2*Y1s**2*d_X1c_d_zeta*d_Y1c_d_zeta)/(Bbar*abs_G0_over_B0**2) + &
             (I2*iota_N*X1c**3*d_Y1s_d_zeta)/(Bbar*abs_G0_over_B0**2) + &
             (I2*iota_N*X1c*Y1c**2*d_Y1s_d_zeta)/(Bbar*abs_G0_over_B0**2) - &
             (2*I2**2*X1c**2*Y1s*d_Y1s_d_zeta)/(Bbar**2*abs_G0_over_B0) + &
             (I2*iota_N*X1c*Y1s**2*d_Y1s_d_zeta)/(Bbar*abs_G0_over_B0**2) + &
             (Y3s1*d_Y1s_d_zeta)/abs_G0_over_B0 + (3*Y3s3*d_Y1s_d_zeta)/abs_G0_over_B0 + &
             (4*iota_N*Y1c*Z2c*d_Y1s_d_zeta)/abs_G0_over_B0**2 - &
             (4*I2*X1c*Z2s*d_Y1s_d_zeta)/(Bbar*abs_G0_over_B0) + &
             (2*iota_N*Y1s*Z2s*d_Y1s_d_zeta)/abs_G0_over_B0**2 + &
             (2*X2s*Y1c*curvature*d_Y1s_d_zeta)/abs_G0_over_B0 - (X20*Y1s*curvature*d_Y1s_d_zeta)/abs_G0_over_B0 - &
             (X2c*Y1s*curvature*d_Y1s_d_zeta)/abs_G0_over_B0 - &
             (X1c**2*Y1s*curvature**2*d_Y1s_d_zeta)/abs_G0_over_B0 + &
             (3*I2*X1c**2*Y1s*torsion*d_Y1s_d_zeta)/(Bbar*abs_G0_over_B0) - &
             (I2*Y1c**2*Y1s*torsion*d_Y1s_d_zeta)/(Bbar*abs_G0_over_B0) + &
             (2*I2*X1c*Y1s*d_Y1c_d_zeta*d_Y1s_d_zeta)/(Bbar*abs_G0_over_B0**2) - &
             (I2*X1c*Y1c*d_Y1s_d_zeta**2)/(Bbar*abs_G0_over_B0**2) + (4*Y2c*d_Y20_d_zeta)/abs_G0_over_B0 + &
             (X1c*Y1c*curvature*d_Y20_d_zeta)/abs_G0_over_B0 - (4*Y20*d_Y2c_d_zeta)/abs_G0_over_B0 - &
             (X1c*Y1c*curvature*d_Y2c_d_zeta)/abs_G0_over_B0 - (2*X1c*Y1s*curvature*d_Y2s_d_zeta)/abs_G0_over_B0 + &
             (Y1c*d_Y3c1_d_zeta)/abs_G0_over_B0 - (3*Y1c*d_Y3c3_d_zeta)/abs_G0_over_B0 - &
             (Y1s*d_Y3s1_d_zeta)/abs_G0_over_B0 - (3*Y1s*d_Y3s3_d_zeta)/abs_G0_over_B0 - &
             (2*iota_N*Y1c*Y1s*d_Z20_d_zeta)/abs_G0_over_B0**2 - &
             (X1c*d_X1c_d_zeta*d_Z20_d_zeta)/abs_G0_over_B0**2 - &
             (Y1c*d_Y1c_d_zeta*d_Z20_d_zeta)/abs_G0_over_B0**2 + &
             (Y1s*d_Y1s_d_zeta*d_Z20_d_zeta)/abs_G0_over_B0**2 + &
             (X1c*d_X1c_d_zeta*d_Z2c_d_zeta)/abs_G0_over_B0**2 + &
             (Y1c*d_Y1c_d_zeta*d_Z2c_d_zeta)/abs_G0_over_B0**2 + &
             (Y1s*d_Y1s_d_zeta*d_Z2c_d_zeta)/abs_G0_over_B0**2 + &
             (2*iota_N*X1c**2*d_Z2s_d_zeta)/abs_G0_over_B0**2 + (2*iota_N*Y1c**2*d_Z2s_d_zeta)/abs_G0_over_B0**2 + &
             (2*iota_N*Y1s**2*d_Z2s_d_zeta)/abs_G0_over_B0**2 + (4*X1c*Y1s*torsion*d_Z2s_d_zeta)/abs_G0_over_B0 + &
             (2*Y1s*d_Y1c_d_zeta*d_Z2s_d_zeta)/abs_G0_over_B0**2 - &
             (2*Y1c*d_Y1s_d_zeta*d_Z2s_d_zeta)/abs_G0_over_B0**2
        
        max_eq1residual = maxval(abs(eq1residual))
        print *,"max(abs(residual of sin(2*theta) mixed-partials equation)):",max_eq1residual
        if (max_eq1residual > 1e-8) stop "Residual is large !!!"
        
        ! cos(2*theta) component of the mixed-partials equation:
        eq1residual = (-8*iota_N*X20*X2c)/abs_G0_over_B0 - (12*iota_N*X1c*X3c3)/abs_G0_over_B0 + (4*I2*X3s1*Y1c)/Bbar - &
             (12*I2*X3s3*Y1c)/Bbar + (4*I2*X3c1*Y1s)/Bbar + (12*I2*X3c3*Y1s)/Bbar + &
             iota2*(X1c**2/abs_G0_over_B0 + Y1c**2/abs_G0_over_B0 - Y1s**2/abs_G0_over_B0) - (16*I2*X2s*Y20)/Bbar - &
             (8*iota_N*Y20*Y2c)/abs_G0_over_B0 + (16*I2*X20*Y2s)/Bbar - (12*iota_N*Y1c*Y3c3)/abs_G0_over_B0 - &
             (4*I2*X1c*Y3s1)/Bbar + (12*I2*X1c*Y3s3)/Bbar - (12*iota_N*Y1s*Y3s3)/abs_G0_over_B0 + &
             (4*iota_N**2*X1c**2*Z2s)/abs_G0_over_B0**2 + (4*iota_N**2*Y1c**2*Z2s)/abs_G0_over_B0**2 - &
             (8*I2*iota_N*X1c*Y1s*Z2s)/(Bbar*abs_G0_over_B0) + (4*iota_N**2*Y1s**2*Z2s)/abs_G0_over_B0**2 + &
             2*X2s*Y1c*beta_1s - 2*X20*Y1s*beta_1s - 2*X2c*Y1s*beta_1s - 2*X1c*Y2s*beta_1s - &
             2*X1c*Y1s*beta_2s - (iota_N*X1c**2*X20*curvature)/abs_G0_over_B0 - (4*iota_N*X1c**2*X2c*curvature)/abs_G0_over_B0 + &
             (iota_N*X20*Y1c**2*curvature)/abs_G0_over_B0 - (2*iota_N*X2c*Y1c**2*curvature)/abs_G0_over_B0 - &
             (iota_N*X20*Y1s**2*curvature)/abs_G0_over_B0 - (2*iota_N*X2c*Y1s**2*curvature)/abs_G0_over_B0 - &
             (2*iota_N*X1c*Y1c*Y20*curvature)/abs_G0_over_B0 - (2*iota_N*X1c*Y1c*Y2c*curvature)/abs_G0_over_B0 - &
             (6*iota_N*X1c*Y1s*Y2s*curvature)/abs_G0_over_B0 - 4*X2s*Z20*curvature + 4*X20*Z2s*curvature - &
             X1c*Z3s1*curvature + 3*X1c*Z3s3*curvature - (iota_N*X1c**4*curvature**2)/(2*abs_G0_over_B0) - &
             (iota_N*X1c**2*Y1c**2*curvature**2)/(2*abs_G0_over_B0) - (3*iota_N*X1c**2*Y1s**2*curvature**2)/(2*abs_G0_over_B0) + &
             X1c**2*Z2s*curvature**2 + (I2*iota_N*X1c**4*torsion)/(2*Bbar*abs_G0_over_B0) - 2*X3s1*Y1c*torsion + &
             6*X3s3*Y1c*torsion + (I2*iota_N*X1c**2*Y1c**2*torsion)/(Bbar*abs_G0_over_B0) + &
             (I2*iota_N*Y1c**4*torsion)/(2*Bbar*abs_G0_over_B0) - (2*I2**2*X1c**3*Y1s*torsion)/Bbar**2 - &
             2*X3c1*Y1s*torsion - 6*X3c3*Y1s*torsion - (2*I2**2*X1c*Y1c**2*Y1s*torsion)/Bbar**2 + &
             (2*I2**2*X1c*Y1s**3*torsion)/Bbar**2 - (I2*iota_N*Y1s**4*torsion)/(2*Bbar*abs_G0_over_B0) + &
             8*X2s*Y20*torsion - 8*X20*Y2s*torsion + 2*X1c*Y3s1*torsion - 6*X1c*Y3s3*torsion + &
             (4*I2*Y1c*Y1s*Z20*torsion)/Bbar - (2*I2*X1c**2*Z2s*torsion)/Bbar - &
             (2*I2*Y1c**2*Z2s*torsion)/Bbar + (8*iota_N*X1c*Y1s*Z2s*torsion)/abs_G0_over_B0 - &
             (2*I2*Y1s**2*Z2s*torsion)/Bbar + 3*X1c*X2s*Y1c*curvature*torsion - &
             3*X1c*X20*Y1s*curvature*torsion - 6*X1c*X2c*Y1s*curvature*torsion - 3*X1c**2*Y2s*curvature*torsion - &
             2*X1c**3*Y1s*curvature**2*torsion + (2*I2*X1c**3*Y1s*torsion**2)/Bbar + &
             (2*I2*X1c*Y1c**2*Y1s*torsion**2)/Bbar - (2*I2*X1c*Y1s**3*torsion**2)/Bbar + &
             (X3s1*d_X1c_d_zeta)/abs_G0_over_B0 - (3*X3s3*d_X1c_d_zeta)/abs_G0_over_B0 - &
             (I2*iota_N*X1c**2*Y1c*d_X1c_d_zeta)/(2*Bbar*abs_G0_over_B0**2) - &
             (I2*iota_N*Y1c**3*d_X1c_d_zeta)/(2*Bbar*abs_G0_over_B0**2) + &
             (2*I2**2*X1c*Y1c*Y1s*d_X1c_d_zeta)/(Bbar**2*abs_G0_over_B0) - &
             (3*I2*iota_N*Y1c*Y1s**2*d_X1c_d_zeta)/(2*Bbar*abs_G0_over_B0**2) - &
             (2*I2*Y1s*Z20*d_X1c_d_zeta)/(Bbar*abs_G0_over_B0) + &
             (2*iota_N*X1c*Z2c*d_X1c_d_zeta)/abs_G0_over_B0**2 - &
             (4*I2*Y1s*Z2c*d_X1c_d_zeta)/(Bbar*abs_G0_over_B0) + &
             (2*I2*Y1c*Z2s*d_X1c_d_zeta)/(Bbar*abs_G0_over_B0) - &
             (X1c*X2s*curvature*d_X1c_d_zeta)/abs_G0_over_B0 - &
             (3*I2*X1c*Y1c*Y1s*torsion*d_X1c_d_zeta)/(Bbar*abs_G0_over_B0) + &
             (I2*X1c*Y1s*d_X1c_d_zeta**2)/(2*Bbar*abs_G0_over_B0**2) - &
             (4*X2s*d_X20_d_zeta)/abs_G0_over_B0 + (4*X20*d_X2s_d_zeta)/abs_G0_over_B0 + &
             (X1c**2*curvature*d_X2s_d_zeta)/abs_G0_over_B0 - (X1c*d_X3s1_d_zeta)/abs_G0_over_B0 + &
             (3*X1c*d_X3s3_d_zeta)/abs_G0_over_B0 + (I2*iota_N*X1c**3*d_Y1c_d_zeta)/(2*Bbar*abs_G0_over_B0**2) + &
             (I2*iota_N*X1c*Y1c**2*d_Y1c_d_zeta)/(2*Bbar*abs_G0_over_B0**2) - &
             (2*I2**2*X1c**2*Y1s*d_Y1c_d_zeta)/(Bbar**2*abs_G0_over_B0) + &
             (3*I2*iota_N*X1c*Y1s**2*d_Y1c_d_zeta)/(2*Bbar*abs_G0_over_B0**2) + &
             (Y3s1*d_Y1c_d_zeta)/abs_G0_over_B0 - (3*Y3s3*d_Y1c_d_zeta)/abs_G0_over_B0 + &
             (2*iota_N*Y1c*Z2c*d_Y1c_d_zeta)/abs_G0_over_B0**2 - &
             (2*I2*X1c*Z2s*d_Y1c_d_zeta)/(Bbar*abs_G0_over_B0) + &
             (4*iota_N*Y1s*Z2s*d_Y1c_d_zeta)/abs_G0_over_B0**2 + &
             (X2s*Y1c*curvature*d_Y1c_d_zeta)/abs_G0_over_B0 - (X20*Y1s*curvature*d_Y1c_d_zeta)/abs_G0_over_B0 - &
             (2*X2c*Y1s*curvature*d_Y1c_d_zeta)/abs_G0_over_B0 - &
             (2*X1c*Y2s*curvature*d_Y1c_d_zeta)/abs_G0_over_B0 - &
             (3*X1c**2*Y1s*curvature**2*d_Y1c_d_zeta)/(2*abs_G0_over_B0) + &
             (7*I2*X1c**2*Y1s*torsion*d_Y1c_d_zeta)/(2*Bbar*abs_G0_over_B0) + &
             (I2*Y1c**2*Y1s*torsion*d_Y1c_d_zeta)/(2*Bbar*abs_G0_over_B0) - &
             (I2*Y1s**3*torsion*d_Y1c_d_zeta)/(2*Bbar*abs_G0_over_B0) - &
             (I2*Y1c*Y1s*d_X1c_d_zeta*d_Y1c_d_zeta)/(Bbar*abs_G0_over_B0**2) + &
             (3*I2*X1c*Y1s*d_Y1c_d_zeta**2)/(2*Bbar*abs_G0_over_B0**2) + &
             (Y3c1*d_Y1s_d_zeta)/abs_G0_over_B0 + (3*Y3c3*d_Y1s_d_zeta)/abs_G0_over_B0 + &
             (2*I2*X1c*Z20*d_Y1s_d_zeta)/(Bbar*abs_G0_over_B0) - &
             (4*I2*X1c*Z2c*d_Y1s_d_zeta)/(Bbar*abs_G0_over_B0) + &
             (2*iota_N*Y1s*Z2c*d_Y1s_d_zeta)/abs_G0_over_B0**2 - &
             (4*iota_N*Y1c*Z2s*d_Y1s_d_zeta)/abs_G0_over_B0**2 - &
             (X20*Y1c*curvature*d_Y1s_d_zeta)/abs_G0_over_B0 + (2*X2c*Y1c*curvature*d_Y1s_d_zeta)/abs_G0_over_B0 + &
             (X2s*Y1s*curvature*d_Y1s_d_zeta)/abs_G0_over_B0 + (2*X1c*Y20*curvature*d_Y1s_d_zeta)/abs_G0_over_B0 + &
             (X1c**2*Y1c*curvature**2*d_Y1s_d_zeta)/(2*abs_G0_over_B0) - &
             (I2*X1c**2*Y1c*torsion*d_Y1s_d_zeta)/(2*Bbar*abs_G0_over_B0) - &
             (I2*Y1c**3*torsion*d_Y1s_d_zeta)/(2*Bbar*abs_G0_over_B0) + &
             (I2*Y1c*Y1s**2*torsion*d_Y1s_d_zeta)/(2*Bbar*abs_G0_over_B0) - &
             (I2*X1c**2*d_X1c_d_zeta*d_Y1s_d_zeta)/(2*Bbar*abs_G0_over_B0**2) + &
             (I2*Y1c**2*d_X1c_d_zeta*d_Y1s_d_zeta)/(2*Bbar*abs_G0_over_B0**2) + &
             (I2*Y1s**2*d_X1c_d_zeta*d_Y1s_d_zeta)/(2*Bbar*abs_G0_over_B0**2) - &
             (I2*X1c*Y1c*d_Y1c_d_zeta*d_Y1s_d_zeta)/(Bbar*abs_G0_over_B0**2) - &
             (I2*X1c*Y1s*d_Y1s_d_zeta**2)/(2*Bbar*abs_G0_over_B0**2) - &
             (4*Y2s*d_Y20_d_zeta)/abs_G0_over_B0 - (X1c*Y1s*curvature*d_Y20_d_zeta)/abs_G0_over_B0 - &
             (2*X1c*Y1s*curvature*d_Y2c_d_zeta)/abs_G0_over_B0 + (4*Y20*d_Y2s_d_zeta)/abs_G0_over_B0 + &
             (X1c*Y1c*curvature*d_Y2s_d_zeta)/abs_G0_over_B0 - (Y1s*d_Y3c1_d_zeta)/abs_G0_over_B0 - &
             (3*Y1s*d_Y3c3_d_zeta)/abs_G0_over_B0 - (Y1c*d_Y3s1_d_zeta)/abs_G0_over_B0 + &
             (3*Y1c*d_Y3s3_d_zeta)/abs_G0_over_B0 - (iota_N*X1c**2*d_Z20_d_zeta)/abs_G0_over_B0**2 - &
             (iota_N*Y1c**2*d_Z20_d_zeta)/abs_G0_over_B0**2 + (iota_N*Y1s**2*d_Z20_d_zeta)/abs_G0_over_B0**2 + &
             (Y1s*d_Y1c_d_zeta*d_Z20_d_zeta)/abs_G0_over_B0**2 + &
             (Y1c*d_Y1s_d_zeta*d_Z20_d_zeta)/abs_G0_over_B0**2 + &
             (2*iota_N*X1c**2*d_Z2c_d_zeta)/abs_G0_over_B0**2 + (2*iota_N*Y1c**2*d_Z2c_d_zeta)/abs_G0_over_B0**2 + &
             (2*iota_N*Y1s**2*d_Z2c_d_zeta)/abs_G0_over_B0**2 + (4*X1c*Y1s*torsion*d_Z2c_d_zeta)/abs_G0_over_B0 + &
             (2*Y1s*d_Y1c_d_zeta*d_Z2c_d_zeta)/abs_G0_over_B0**2 - &
             (2*Y1c*d_Y1s_d_zeta*d_Z2c_d_zeta)/abs_G0_over_B0**2 - &
             (X1c*d_X1c_d_zeta*d_Z2s_d_zeta)/abs_G0_over_B0**2 - &
             (Y1c*d_Y1c_d_zeta*d_Z2s_d_zeta)/abs_G0_over_B0**2 - &
             (Y1s*d_Y1s_d_zeta*d_Z2s_d_zeta)/abs_G0_over_B0**2
        
        max_eq1residual = maxval(abs(eq1residual))
        print *,"max(abs(residual of cos(2*theta) mixed-partials equation)):",max_eq1residual
        if (max_eq1residual > 1e-8) stop "Residual is large !!!"
        
        ! cos(0*theta) component of the mixed-partials equation:
        eq1residual = (-8*iota_N*X2c**2)/abs_G0_over_B0 - (8*iota_N*X2s**2)/abs_G0_over_B0 - (4*iota_N*X1c*X3c1)/abs_G0_over_B0 - &
             (8*I2*X3s1*Y1c)/Bbar + (4*I4*X1c*Y1s)/Bbar + (8*I2*X3c1*Y1s)/Bbar + &
             iota2*((-2*X1c**2)/abs_G0_over_B0 - (2*Y1c**2)/abs_G0_over_B0 - (2*Y1s**2)/abs_G0_over_B0) - (16*I2*X2s*Y2c)/Bbar - &
             (8*iota_N*Y2c**2)/abs_G0_over_B0 + (16*I2*X2c*Y2s)/Bbar - (8*iota_N*Y2s**2)/abs_G0_over_B0 - &
             (4*iota_N*Y1c*Y3c1)/abs_G0_over_B0 + (8*I2*X1c*Y3s1)/Bbar - (4*iota_N*Y1s*Y3s1)/abs_G0_over_B0 + &
             (4*iota_N**2*Y1c*Y1s*Z2c)/abs_G0_over_B0**2 - (2*iota_N**2*X1c**2*Z2s)/abs_G0_over_B0**2 - &
             (2*iota_N**2*Y1c**2*Z2s)/abs_G0_over_B0**2 + (2*iota_N**2*Y1s**2*Z2s)/abs_G0_over_B0**2 - &
             (2*iota_N*X1c**2*X20*curvature)/abs_G0_over_B0 - (3*iota_N*X1c**2*X2c*curvature)/abs_G0_over_B0 - &
             (2*iota_N*X20*Y1c**2*curvature)/abs_G0_over_B0 + (iota_N*X2c*Y1c**2*curvature)/abs_G0_over_B0 + &
             (2*iota_N*X2s*Y1c*Y1s*curvature)/abs_G0_over_B0 - (2*iota_N*X20*Y1s**2*curvature)/abs_G0_over_B0 - &
             (iota_N*X2c*Y1s**2*curvature)/abs_G0_over_B0 - (4*iota_N*X1c*Y1c*Y2c*curvature)/abs_G0_over_B0 - &
             (4*iota_N*X1c*Y1s*Y2s*curvature)/abs_G0_over_B0 - 4*X2s*Z2c*curvature + 4*X2c*Z2s*curvature + &
             2*X1c*Z3s1*curvature - (iota_N*X1c**4*curvature**2)/(2*abs_G0_over_B0) - &
             (iota_N*X1c**2*Y1c**2*curvature**2)/(2*abs_G0_over_B0) - (3*iota_N*X1c**2*Y1s**2*curvature**2)/(2*abs_G0_over_B0) + &
             X1c**2*Z2s*curvature**2 + (I2*iota_N*X1c**4*torsion)/(2*Bbar*abs_G0_over_B0) + 4*X3s1*Y1c*torsion + &
             (I2*iota_N*X1c**2*Y1c**2*torsion)/(Bbar*abs_G0_over_B0) + (I2*iota_N*Y1c**4*torsion)/(2*Bbar*abs_G0_over_B0) - &
             (2*I2**2*X1c**3*Y1s*torsion)/Bbar**2 - 4*X3c1*Y1s*torsion - &
             (2*I2**2*X1c*Y1c**2*Y1s*torsion)/Bbar**2 + (3*I2*iota_N*X1c**2*Y1s**2*torsion)/(Bbar*abs_G0_over_B0) + &
             (I2*iota_N*Y1c**2*Y1s**2*torsion)/(Bbar*abs_G0_over_B0) - (2*I2**2*X1c*Y1s**3*torsion)/Bbar**2 + &
             (I2*iota_N*Y1s**4*torsion)/(2*Bbar*abs_G0_over_B0) + 8*X2s*Y2c*torsion - 8*X2c*Y2s*torsion - &
             4*X1c*Y3s1*torsion + (4*I2*Y1c*Y1s*Z2c*torsion)/Bbar - (2*I2*X1c**2*Z2s*torsion)/Bbar - &
             (2*I2*Y1c**2*Z2s*torsion)/Bbar + (2*I2*Y1s**2*Z2s*torsion)/Bbar + &
             3*X1c*X2s*Y1c*curvature*torsion - 6*X1c*X20*Y1s*curvature*torsion - &
             3*X1c*X2c*Y1s*curvature*torsion - 3*X1c**2*Y2s*curvature*torsion - 2*X1c**3*Y1s*curvature**2*torsion + &
             (2*I2*X1c**3*Y1s*torsion**2)/Bbar + (2*I2*X1c*Y1c**2*Y1s*torsion**2)/Bbar + &
             (2*I2*X1c*Y1s**3*torsion**2)/Bbar - (2*X3s1*d_X1c_d_zeta)/abs_G0_over_B0 - &
             (I2*iota_N*X1c**2*Y1c*d_X1c_d_zeta)/(2*Bbar*abs_G0_over_B0**2) - &
             (I2*iota_N*Y1c**3*d_X1c_d_zeta)/(2*Bbar*abs_G0_over_B0**2) + &
             (2*I2**2*X1c*Y1c*Y1s*d_X1c_d_zeta)/(Bbar**2*abs_G0_over_B0) - &
             (I2*iota_N*Y1c*Y1s**2*d_X1c_d_zeta)/(2*Bbar*abs_G0_over_B0**2) - &
             (4*I2*Y1s*Z20*d_X1c_d_zeta)/(Bbar*abs_G0_over_B0) + &
             (2*iota_N*X1c*Z2c*d_X1c_d_zeta)/abs_G0_over_B0**2 - &
             (2*I2*Y1s*Z2c*d_X1c_d_zeta)/(Bbar*abs_G0_over_B0) + &
             (2*I2*Y1c*Z2s*d_X1c_d_zeta)/(Bbar*abs_G0_over_B0) - &
             (X1c*X2s*curvature*d_X1c_d_zeta)/abs_G0_over_B0 - &
             (3*I2*X1c*Y1c*Y1s*torsion*d_X1c_d_zeta)/(Bbar*abs_G0_over_B0) + &
             (I2*X1c*Y1s*d_X1c_d_zeta**2)/(2*Bbar*abs_G0_over_B0**2) - &
             (4*X2s*d_X2c_d_zeta)/abs_G0_over_B0 + (4*X2c*d_X2s_d_zeta)/abs_G0_over_B0 + &
             (X1c**2*curvature*d_X2s_d_zeta)/abs_G0_over_B0 + (2*X1c*d_X3s1_d_zeta)/abs_G0_over_B0 + &
             (I2*iota_N*X1c**3*d_Y1c_d_zeta)/(2*Bbar*abs_G0_over_B0**2) + &
             (I2*iota_N*X1c*Y1c**2*d_Y1c_d_zeta)/(2*Bbar*abs_G0_over_B0**2) - &
             (2*I2**2*X1c**2*Y1s*d_Y1c_d_zeta)/(Bbar**2*abs_G0_over_B0) + &
             (3*I2*iota_N*X1c*Y1s**2*d_Y1c_d_zeta)/(2*Bbar*abs_G0_over_B0**2) - &
             (2*Y3s1*d_Y1c_d_zeta)/abs_G0_over_B0 + (2*iota_N*Y1c*Z2c*d_Y1c_d_zeta)/abs_G0_over_B0**2 - &
             (2*I2*X1c*Z2s*d_Y1c_d_zeta)/(Bbar*abs_G0_over_B0) + &
             (2*iota_N*Y1s*Z2s*d_Y1c_d_zeta)/abs_G0_over_B0**2 + &
             (X2s*Y1c*curvature*d_Y1c_d_zeta)/abs_G0_over_B0 - (2*X20*Y1s*curvature*d_Y1c_d_zeta)/abs_G0_over_B0 - &
             (X2c*Y1s*curvature*d_Y1c_d_zeta)/abs_G0_over_B0 - (2*X1c*Y2s*curvature*d_Y1c_d_zeta)/abs_G0_over_B0 - &
             (3*X1c**2*Y1s*curvature**2*d_Y1c_d_zeta)/(2*abs_G0_over_B0) + &
             (7*I2*X1c**2*Y1s*torsion*d_Y1c_d_zeta)/(2*Bbar*abs_G0_over_B0) + &
             (I2*Y1c**2*Y1s*torsion*d_Y1c_d_zeta)/(2*Bbar*abs_G0_over_B0) + &
             (I2*Y1s**3*torsion*d_Y1c_d_zeta)/(2*Bbar*abs_G0_over_B0) - &
             (I2*Y1c*Y1s*d_X1c_d_zeta*d_Y1c_d_zeta)/(Bbar*abs_G0_over_B0**2) + &
             (3*I2*X1c*Y1s*d_Y1c_d_zeta**2)/(2*Bbar*abs_G0_over_B0**2) - &
             (I2*iota_N*X1c*Y1c*Y1s*d_Y1s_d_zeta)/(Bbar*abs_G0_over_B0**2) + &
             (2*Y3c1*d_Y1s_d_zeta)/abs_G0_over_B0 - (4*I2*X1c*Z20*d_Y1s_d_zeta)/(Bbar*abs_G0_over_B0) + &
             (2*I2*X1c*Z2c*d_Y1s_d_zeta)/(Bbar*abs_G0_over_B0) - &
             (2*iota_N*Y1s*Z2c*d_Y1s_d_zeta)/abs_G0_over_B0**2 + &
             (2*iota_N*Y1c*Z2s*d_Y1s_d_zeta)/abs_G0_over_B0**2 + &
             (2*X20*Y1c*curvature*d_Y1s_d_zeta)/abs_G0_over_B0 - (X2c*Y1c*curvature*d_Y1s_d_zeta)/abs_G0_over_B0 - &
             (X2s*Y1s*curvature*d_Y1s_d_zeta)/abs_G0_over_B0 + (2*X1c*Y2c*curvature*d_Y1s_d_zeta)/abs_G0_over_B0 + &
             (X1c**2*Y1c*curvature**2*d_Y1s_d_zeta)/(2*abs_G0_over_B0) - &
             (I2*X1c**2*Y1c*torsion*d_Y1s_d_zeta)/(2*Bbar*abs_G0_over_B0) - &
             (I2*Y1c**3*torsion*d_Y1s_d_zeta)/(2*Bbar*abs_G0_over_B0) - &
             (I2*Y1c*Y1s**2*torsion*d_Y1s_d_zeta)/(2*Bbar*abs_G0_over_B0) - &
             (I2*X1c**2*d_X1c_d_zeta*d_Y1s_d_zeta)/(2*Bbar*abs_G0_over_B0**2) + &
             (I2*Y1c**2*d_X1c_d_zeta*d_Y1s_d_zeta)/(2*Bbar*abs_G0_over_B0**2) - &
             (I2*Y1s**2*d_X1c_d_zeta*d_Y1s_d_zeta)/(2*Bbar*abs_G0_over_B0**2) - &
             (I2*X1c*Y1c*d_Y1c_d_zeta*d_Y1s_d_zeta)/(Bbar*abs_G0_over_B0**2) + &
             (I2*X1c*Y1s*d_Y1s_d_zeta**2)/(2*Bbar*abs_G0_over_B0**2) - &
             (2*X1c*Y1s*curvature*d_Y20_d_zeta)/abs_G0_over_B0 - (4*Y2s*d_Y2c_d_zeta)/abs_G0_over_B0 - &
             (X1c*Y1s*curvature*d_Y2c_d_zeta)/abs_G0_over_B0 + (4*Y2c*d_Y2s_d_zeta)/abs_G0_over_B0 + &
             (X1c*Y1c*curvature*d_Y2s_d_zeta)/abs_G0_over_B0 - (2*Y1s*d_Y3c1_d_zeta)/abs_G0_over_B0 + &
             (2*Y1c*d_Y3s1_d_zeta)/abs_G0_over_B0 + (2*iota_N*X1c**2*d_Z20_d_zeta)/abs_G0_over_B0**2 + &
             (2*iota_N*Y1c**2*d_Z20_d_zeta)/abs_G0_over_B0**2 + (2*iota_N*Y1s**2*d_Z20_d_zeta)/abs_G0_over_B0**2 + &
             (4*X1c*Y1s*torsion*d_Z20_d_zeta)/abs_G0_over_B0 + &
             (2*Y1s*d_Y1c_d_zeta*d_Z20_d_zeta)/abs_G0_over_B0**2 - &
             (2*Y1c*d_Y1s_d_zeta*d_Z20_d_zeta)/abs_G0_over_B0**2 - &
             (iota_N*X1c**2*d_Z2c_d_zeta)/abs_G0_over_B0**2 - (iota_N*Y1c**2*d_Z2c_d_zeta)/abs_G0_over_B0**2 + &
             (iota_N*Y1s**2*d_Z2c_d_zeta)/abs_G0_over_B0**2 + &
             (Y1s*d_Y1c_d_zeta*d_Z2c_d_zeta)/abs_G0_over_B0**2 + &
             (Y1c*d_Y1s_d_zeta*d_Z2c_d_zeta)/abs_G0_over_B0**2 - &
             (2*iota_N*Y1c*Y1s*d_Z2s_d_zeta)/abs_G0_over_B0**2 - &
             (X1c*d_X1c_d_zeta*d_Z2s_d_zeta)/abs_G0_over_B0**2 - &
             (Y1c*d_Y1c_d_zeta*d_Z2s_d_zeta)/abs_G0_over_B0**2 + &
             (Y1s*d_Y1s_d_zeta*d_Z2s_d_zeta)/abs_G0_over_B0**2
        
        max_eq1residual = maxval(abs(eq1residual))
        print *,"max(abs(residual of cos(0*theta) mixed-partials equation)):",max_eq1residual
        if (max_eq1residual > 1e-8) stop "Residual is large !!!"

        deallocate(d_X3c1_d_zeta, d_X3s1_d_zeta, d_Y3c1_d_zeta, d_Y3s1_d_zeta, d_Y3c3_d_zeta, d_Y3s3_d_zeta)
     end if

     ! Now that we have X3, compute B3.

     B3s1 = -(B0**3*((-3*G0**2*B1c*B2s)/B0**4 + 2*iota_N**2*X1c*X2s - 2*iota_N**2*Y1s*Y2c + &
            2*iota_N**2*Y1c*Y2s - 2*abs_G0_over_B0*iota_N*Z3c1 - 2*abs_G0_over_B0**2*X3s1*curvature - &
            2*abs_G0_over_B0*iota_N*X1c*Z20*curvature + 3*abs_G0_over_B0*iota_N*X1c*Z2c*curvature + abs_G0_over_B0**2*X1c*X2s*curvature**2 - &
            2*abs_G0_over_B0*iota_N*X20*Y1c*torsion + 3*abs_G0_over_B0*iota_N*X2c*Y1c*torsion + 3*abs_G0_over_B0*iota_N*X2s*Y1s*torsion + &
            2*abs_G0_over_B0*iota_N*X1c*Y20*torsion - 3*abs_G0_over_B0*iota_N*X1c*Y2c*torsion - &
            2*abs_G0_over_B0**2*Y1s*Z20*curvature*torsion + abs_G0_over_B0**2*Y1s*Z2c*curvature*torsion - &
            abs_G0_over_B0**2*Y1c*Z2s*curvature*torsion + abs_G0_over_B0**2*X1c*X2s*torsion**2 + &
            2*abs_G0_over_B0**2*Y1s*Y20*torsion**2 - abs_G0_over_B0**2*Y1s*Y2c*torsion**2 + abs_G0_over_B0**2*Y1c*Y2s*torsion**2 - &
            2*iota_N*X2c*d_X1c_d_zeta + abs_G0_over_B0*Z2s*curvature*d_X1c_d_zeta - &
            abs_G0_over_B0*Y2s*torsion*d_X1c_d_zeta - 2*iota_N*X1c*d_X20_d_zeta - &
            2*abs_G0_over_B0*Y1s*torsion*d_X20_d_zeta + iota_N*X1c*d_X2c_d_zeta + &
            abs_G0_over_B0*Y1s*torsion*d_X2c_d_zeta - abs_G0_over_B0*Y1c*torsion*d_X2s_d_zeta + &
            d_X1c_d_zeta*d_X2s_d_zeta - 2*iota_N*Y2c*d_Y1c_d_zeta + &
            abs_G0_over_B0*X2s*torsion*d_Y1c_d_zeta - 2*iota_N*Y2s*d_Y1s_d_zeta + &
            2*abs_G0_over_B0*X20*torsion*d_Y1s_d_zeta - abs_G0_over_B0*X2c*torsion*d_Y1s_d_zeta - &
            2*iota_N*Y1c*d_Y20_d_zeta + 2*d_Y1s_d_zeta*d_Y20_d_zeta + &
            iota_N*Y1c*d_Y2c_d_zeta - d_Y1s_d_zeta*d_Y2c_d_zeta + &
            iota_N*Y1s*d_Y2s_d_zeta + abs_G0_over_B0*X1c*torsion*d_Y2s_d_zeta + &
            d_Y1c_d_zeta*d_Y2s_d_zeta - abs_G0_over_B0*X1c*curvature*d_Z2s_d_zeta + &
            2*abs_G0_over_B0*d_Z3s1_d_zeta))/(2*G0**2)

     B3c1 = -(B0**3*((4*G0*G2*B1c)/B0**3 + (4*G0*I2*N_helicity*B1c)/B0**3 + (4*G0*I2*iota_N*B1c)/B0**3 + &
            (3*G0**2*B1c**3)/B0**5 - (6*G0**2*B1c*B20)/B0**4 - (3*G0**2*B1c*B2c)/B0**4 + &
            2*iota_N**2*X1c*X2c + 2*iota_N**2*Y1c*Y2c + 2*iota_N**2*Y1s*Y2s + 2*abs_G0_over_B0*iota_N*Z3s1 - &
            2*abs_G0_over_B0**2*X3c1*curvature - 3*abs_G0_over_B0*iota_N*X1c*Z2s*curvature + 2*abs_G0_over_B0**2*X1c*X20*curvature**2 + &
            abs_G0_over_B0**2*X1c*X2c*curvature**2 - 3*abs_G0_over_B0*iota_N*X2s*Y1c*torsion + 2*abs_G0_over_B0*iota_N*X20*Y1s*torsion + &
            3*abs_G0_over_B0*iota_N*X2c*Y1s*torsion + 3*abs_G0_over_B0*iota_N*X1c*Y2s*torsion - &
            2*abs_G0_over_B0**2*Y1c*Z20*curvature*torsion - abs_G0_over_B0**2*Y1c*Z2c*curvature*torsion - &
            abs_G0_over_B0**2*Y1s*Z2s*curvature*torsion + 2*abs_G0_over_B0**2*X1c*X20*torsion**2 + &
            abs_G0_over_B0**2*X1c*X2c*torsion**2 + 2*abs_G0_over_B0**2*Y1c*Y20*torsion**2 + abs_G0_over_B0**2*Y1c*Y2c*torsion**2 + &
            abs_G0_over_B0**2*Y1s*Y2s*torsion**2 + 2*iota_N*X2s*d_X1c_d_zeta + &
            2*abs_G0_over_B0*Z20*curvature*d_X1c_d_zeta + abs_G0_over_B0*Z2c*curvature*d_X1c_d_zeta - &
            2*abs_G0_over_B0*Y20*torsion*d_X1c_d_zeta - abs_G0_over_B0*Y2c*torsion*d_X1c_d_zeta - &
            2*abs_G0_over_B0*Y1c*torsion*d_X20_d_zeta + 2*d_X1c_d_zeta*d_X20_d_zeta - &
            abs_G0_over_B0*Y1c*torsion*d_X2c_d_zeta + d_X1c_d_zeta*d_X2c_d_zeta - &
            iota_N*X1c*d_X2s_d_zeta - abs_G0_over_B0*Y1s*torsion*d_X2s_d_zeta + &
            2*iota_N*Y2s*d_Y1c_d_zeta + 2*abs_G0_over_B0*X20*torsion*d_Y1c_d_zeta + &
            abs_G0_over_B0*X2c*torsion*d_Y1c_d_zeta - 2*iota_N*Y2c*d_Y1s_d_zeta + &
            abs_G0_over_B0*X2s*torsion*d_Y1s_d_zeta + 2*iota_N*Y1s*d_Y20_d_zeta + &
            2*abs_G0_over_B0*X1c*torsion*d_Y20_d_zeta + 2*d_Y1c_d_zeta*d_Y20_d_zeta + &
            iota_N*Y1s*d_Y2c_d_zeta + abs_G0_over_B0*X1c*torsion*d_Y2c_d_zeta + &
            d_Y1c_d_zeta*d_Y2c_d_zeta - iota_N*Y1c*d_Y2s_d_zeta + &
            d_Y1s_d_zeta*d_Y2s_d_zeta - 2*abs_G0_over_B0*X1c*curvature*d_Z20_d_zeta - &
            abs_G0_over_B0*X1c*curvature*d_Z2c_d_zeta + 2*abs_G0_over_B0*d_Z3c1_d_zeta))/(2*G0**2)

     B3s3 = -(B0**3*((-3*G0**2*B1c*B2s)/B0**4 - 2*iota_N**2*X1c*X2s - 2*iota_N**2*Y1s*Y2c - &
            2*iota_N**2*Y1c*Y2s - 6*abs_G0_over_B0*iota_N*Z3c3 - 2*abs_G0_over_B0**2*X3s3*curvature + &
            abs_G0_over_B0*iota_N*X1c*Z2c*curvature + abs_G0_over_B0**2*X1c*X2s*curvature**2 + abs_G0_over_B0*iota_N*X2c*Y1c*torsion - &
            abs_G0_over_B0*iota_N*X2s*Y1s*torsion - abs_G0_over_B0*iota_N*X1c*Y2c*torsion - abs_G0_over_B0**2*Y1s*Z2c*curvature*torsion - &
            abs_G0_over_B0**2*Y1c*Z2s*curvature*torsion + abs_G0_over_B0**2*X1c*X2s*torsion**2 + abs_G0_over_B0**2*Y1s*Y2c*torsion**2 + &
            abs_G0_over_B0**2*Y1c*Y2s*torsion**2 - 2*iota_N*X2c*d_X1c_d_zeta + &
            abs_G0_over_B0*Z2s*curvature*d_X1c_d_zeta - abs_G0_over_B0*Y2s*torsion*d_X1c_d_zeta - &
            iota_N*X1c*d_X2c_d_zeta - abs_G0_over_B0*Y1s*torsion*d_X2c_d_zeta - &
            abs_G0_over_B0*Y1c*torsion*d_X2s_d_zeta + d_X1c_d_zeta*d_X2s_d_zeta - &
            2*iota_N*Y2c*d_Y1c_d_zeta + abs_G0_over_B0*X2s*torsion*d_Y1c_d_zeta + &
            2*iota_N*Y2s*d_Y1s_d_zeta + abs_G0_over_B0*X2c*torsion*d_Y1s_d_zeta - &
            iota_N*Y1c*d_Y2c_d_zeta + d_Y1s_d_zeta*d_Y2c_d_zeta + &
            iota_N*Y1s*d_Y2s_d_zeta + abs_G0_over_B0*X1c*torsion*d_Y2s_d_zeta + &
            d_Y1c_d_zeta*d_Y2s_d_zeta - abs_G0_over_B0*X1c*curvature*d_Z2s_d_zeta + &
            2*abs_G0_over_B0*d_Z3s3_d_zeta))/(2*G0**2)

     B3c3 = -(B0**3*((G0**2*B1c**3)/B0**5 - (3*G0**2*B1c*B2c)/B0**4 - 2*iota_N**2*X1c*X2c - &
            2*iota_N**2*Y1c*Y2c + 2*iota_N**2*Y1s*Y2s + 6*abs_G0_over_B0*iota_N*Z3s3 - 2*abs_G0_over_B0**2*X3c3*curvature - &
            abs_G0_over_B0*iota_N*X1c*Z2s*curvature + abs_G0_over_B0**2*X1c*X2c*curvature**2 - abs_G0_over_B0*iota_N*X2s*Y1c*torsion - &
            abs_G0_over_B0*iota_N*X2c*Y1s*torsion + abs_G0_over_B0*iota_N*X1c*Y2s*torsion - abs_G0_over_B0**2*Y1c*Z2c*curvature*torsion + &
            abs_G0_over_B0**2*Y1s*Z2s*curvature*torsion + abs_G0_over_B0**2*X1c*X2c*torsion**2 + abs_G0_over_B0**2*Y1c*Y2c*torsion**2 - &
            abs_G0_over_B0**2*Y1s*Y2s*torsion**2 + 2*iota_N*X2s*d_X1c_d_zeta + &
            abs_G0_over_B0*Z2c*curvature*d_X1c_d_zeta - abs_G0_over_B0*Y2c*torsion*d_X1c_d_zeta - &
            abs_G0_over_B0*Y1c*torsion*d_X2c_d_zeta + d_X1c_d_zeta*d_X2c_d_zeta + &
            iota_N*X1c*d_X2s_d_zeta + abs_G0_over_B0*Y1s*torsion*d_X2s_d_zeta + &
            2*iota_N*Y2s*d_Y1c_d_zeta + abs_G0_over_B0*X2c*torsion*d_Y1c_d_zeta + &
            2*iota_N*Y2c*d_Y1s_d_zeta - abs_G0_over_B0*X2s*torsion*d_Y1s_d_zeta + &
            iota_N*Y1s*d_Y2c_d_zeta + abs_G0_over_B0*X1c*torsion*d_Y2c_d_zeta + &
            d_Y1c_d_zeta*d_Y2c_d_zeta + iota_N*Y1c*d_Y2s_d_zeta - &
            d_Y1s_d_zeta*d_Y2s_d_zeta - abs_G0_over_B0*X1c*curvature*d_Z2c_d_zeta + &
            2*abs_G0_over_B0*d_Z3c3_d_zeta))/(2*G0**2)

     deallocate(d_Z3s1_d_zeta, d_Z3c1_d_zeta, d_Z3s3_d_zeta, d_Z3c3_d_zeta)
     deallocate(d_X3s3_d_zeta, d_X3c3_d_zeta)
  end if
  !print *,"QQQ"
  if(trim(order_r_option)==order_r_option_r3_B3) then
     ! Sanity test: we got back the B3 we requested:

     max_eq1residual = maxval(abs(B3c3 - B3c3_input))
     print *,"max(abs(B3c3-B3c3_input)):",max_eq1residual
     if (max_eq1residual > 1e-13) stop "B3c3 difference is large !!!"
     
     max_eq1residual = maxval(abs(B3s3 - B3s3_input))
     print *,"max(abs(B3s3-B3s3_input)):",max_eq1residual
     if (max_eq1residual > 1e-13) stop "B3s3 difference is large !!!"
  end if

  if(allocated(d_X20_d_zeta)) deallocate(d_X20_d_zeta)
  if(allocated(d_X2s_d_zeta)) deallocate(d_X2s_d_zeta)
  if(allocated(d_X2c_d_zeta)) deallocate(d_X2c_d_zeta)
  if(allocated(d_Y20_d_zeta)) deallocate(d_Y20_d_zeta)
  if(allocated(d_Y2s_d_zeta)) deallocate(d_Y2s_d_zeta)
  if(allocated(d_Y2c_d_zeta)) deallocate(d_Y2c_d_zeta)
  if(allocated(d_Z20_d_zeta)) deallocate(d_Z20_d_zeta)
  if(allocated(d_Z2s_d_zeta)) deallocate(d_Z2s_d_zeta)
  if(allocated(d_Z2c_d_zeta)) deallocate(d_Z2c_d_zeta)

end subroutine quasisymmetry_higher_order_in_r
