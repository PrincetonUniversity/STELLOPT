subroutine quasisymmetry_init_axis

  use quasisymmetry_variables

  implicit none

  integer :: n, j
  real(dp), dimension(:), allocatable :: sinangle, cosangle
  real(dp), dimension(:), allocatable :: d2_l_d_phi2, torsion_numerator, torsion_denominator
  real(dp), dimension(:,:), allocatable :: d_r_d_phi_cylindrical, d2_r_d_phi2_cylindrical, d3_r_d_phi3_cylindrical
  real(dp), dimension(:,:), allocatable :: d_tangent_d_l_cylindrical

  if (allocated(R0)) deallocate(R0)
  if (allocated(Z0)) deallocate(Z0)
  if (allocated(R0_extended)) deallocate(R0_extended)
  if (allocated(Z0_extended)) deallocate(Z0_extended)
  if (allocated(R0p)) deallocate(R0p)
  if (allocated(Z0p)) deallocate(Z0p)
  if (allocated(R0pp)) deallocate(R0pp)
  if (allocated(Z0pp)) deallocate(Z0pp)
  if (allocated(R0ppp)) deallocate(R0ppp)
  if (allocated(Z0ppp)) deallocate(Z0ppp)
  if (allocated(d_l_d_phi)) deallocate(d_l_d_phi)
  if (allocated(curvature)) deallocate(curvature)
  if (allocated(torsion)) deallocate(torsion)
  if (allocated(tangent_cylindrical)) deallocate(tangent_cylindrical)
  if (allocated(normal_cylindrical)) deallocate(normal_cylindrical)
  if (allocated(binormal_cylindrical)) deallocate(binormal_cylindrical)
  if (allocated(tangent_Cartesian)) deallocate(tangent_Cartesian)
  if (allocated(normal_Cartesian)) deallocate(normal_Cartesian)
  if (allocated(binormal_Cartesian)) deallocate(binormal_Cartesian)
  if (allocated(B1Squared_over_curvatureSquared)) deallocate(B1Squared_over_curvatureSquared)

  allocate(sinangle(N_phi))
  allocate(cosangle(N_phi))
  allocate(R0(N_phi))
  allocate(Z0(N_phi))
  allocate(R0_extended(N_phi*nfp))
  allocate(Z0_extended(N_phi*nfp))
  allocate(R0p(N_phi))
  allocate(Z0p(N_phi))
  allocate(R0pp(N_phi))
  allocate(Z0pp(N_phi))
  allocate(R0ppp(N_phi))
  allocate(Z0ppp(N_phi))

  allocate(d_l_d_phi(N_phi))
  allocate(d2_l_d_phi2(N_phi))
  allocate(curvature(N_phi))
  allocate(torsion(N_phi))
  allocate(torsion_numerator(N_phi))
  allocate(torsion_denominator(N_phi))

  allocate(d_r_d_phi_cylindrical(N_phi,3))
  allocate(d2_r_d_phi2_cylindrical(N_phi,3))
  allocate(d3_r_d_phi3_cylindrical(N_phi,3))
  allocate(tangent_cylindrical(N_phi,3))
  allocate(d_tangent_d_l_cylindrical(N_phi,3))
  allocate(normal_cylindrical(N_phi,3))
  allocate(binormal_cylindrical(N_phi,3))
  allocate(tangent_Cartesian(N_phi,3))
  allocate(normal_Cartesian(N_phi,3))
  allocate(binormal_Cartesian(N_phi,3))

  R0 = R0c(1)
  Z0 = Z0c(1)
  R0_extended = R0c(1)
  Z0_extended = Z0c(1)
  R0p = 0
  Z0p = 0
  R0pp = 0
  Z0pp = 0
  R0ppp = 0
  Z0ppp = 0

  do n = 1, axis_nmax
     !sinangle = sin(n*nfp*phi)
     !cosangle = cos(n*nfp*phi)
     sinangle = sin_n_phi(:,n+1)
     cosangle = cos_n_phi(:,n+1)

     R0 = R0 + R0c(n+1) * cosangle + R0s(n+1) * sinangle
     Z0 = Z0 + Z0c(n+1) * cosangle + Z0s(n+1) * sinangle
     R0_extended = R0_extended + R0c(n+1) * cosangle + R0s(n+1) * sinangle
     Z0_extended = Z0_extended + R0c(n+1) * cosangle + R0s(n+1) * sinangle
     R0p = R0p + R0c(n+1) * (-n*nfp)*sinangle + R0s(n+1) * (n*nfp)*cosangle
     Z0p = Z0p + Z0c(n+1) * (-n*nfp)*sinangle + Z0s(n+1) * (n*nfp)*cosangle
     R0pp = R0pp + R0c(n+1) * (-n*nfp*n*nfp)*cosangle + R0s(n+1) * (-n*nfp*n*nfp)*sinangle
     Z0pp = Z0pp + Z0c(n+1) * (-n*nfp*n*nfp)*cosangle + Z0s(n+1) * (-n*nfp*n*nfp)*sinangle
     R0ppp = R0ppp + R0c(n+1) * (n*nfp*n*nfp*n*nfp)*sinangle + R0s(n+1) * (-n*nfp*n*nfp*n*nfp)*cosangle
     Z0ppp = Z0ppp + Z0c(n+1) * (n*nfp*n*nfp*n*nfp)*sinangle + Z0s(n+1) * (-n*nfp*n*nfp*n*nfp)*cosangle

  end do

  d_l_d_phi = sqrt(R0 * R0 + R0p * R0p + Z0p * Z0p)
  d2_l_d_phi2 = (R0 * R0p + R0p * R0pp + Z0p * Z0pp) / d_l_d_phi
  B0_over_abs_G0 = N_phi / sum(d_l_d_phi)

  d_r_d_phi_cylindrical(:,1) = R0p
  d_r_d_phi_cylindrical(:,2) = R0
  d_r_d_phi_cylindrical(:,3) = Z0p

  d2_r_d_phi2_cylindrical(:,1) = R0pp - R0
  d2_r_d_phi2_cylindrical(:,2) = 2 * R0p
  d2_r_d_phi2_cylindrical(:,3) = Z0pp

  d3_r_d_phi3_cylindrical(:,1) = R0ppp - 3*R0p
  d3_r_d_phi3_cylindrical(:,2) = 3 * R0pp - R0
  d3_r_d_phi3_cylindrical(:,3) = Z0ppp

  ! Compute the Frenet-Serret frame:

  do j = 1,3
     tangent_cylindrical(:,j) = d_r_d_phi_cylindrical(:,j) / d_l_d_phi

     d_tangent_d_l_cylindrical(:,j) = (-d_r_d_phi_cylindrical(:,j) * d2_l_d_phi2 / d_l_d_phi &
          + d2_r_d_phi2_cylindrical(:,j)) / (d_l_d_phi * d_l_d_phi)
  end do

  curvature = sqrt( &
       d_tangent_d_l_cylindrical(:,1) * d_tangent_d_l_cylindrical(:,1) + &
       d_tangent_d_l_cylindrical(:,2) * d_tangent_d_l_cylindrical(:,2) + &
       d_tangent_d_l_cylindrical(:,3) * d_tangent_d_l_cylindrical(:,3))

  do j = 1,3
     normal_cylindrical(:,j) = d_tangent_d_l_cylindrical(:,j) / curvature
  end do
  
  ! b = t x n
  binormal_cylindrical(:,1) = tangent_cylindrical(:,2) * normal_cylindrical(:,3) - tangent_cylindrical(:,3) * normal_cylindrical(:,2)
  binormal_cylindrical(:,2) = tangent_cylindrical(:,3) * normal_cylindrical(:,1) - tangent_cylindrical(:,1) * normal_cylindrical(:,3)
  binormal_cylindrical(:,3) = tangent_cylindrical(:,1) * normal_cylindrical(:,2) - tangent_cylindrical(:,2) * normal_cylindrical(:,1)

  ! The minus sign in the next line is absent in wikipedia and mathworld.wolfram.com/Torsion.html but
  ! present in Garren & Boozer's sign convention.
  torsion_numerator = -(0 &
       + d_r_d_phi_cylindrical(:,1) * (d2_r_d_phi2_cylindrical(:,2) * d3_r_d_phi3_cylindrical(:,3) - d2_r_d_phi2_cylindrical(:,3) * d3_r_d_phi3_cylindrical(:,2)) &
       + d_r_d_phi_cylindrical(:,2) * (d2_r_d_phi2_cylindrical(:,3) * d3_r_d_phi3_cylindrical(:,1) - d2_r_d_phi2_cylindrical(:,1) * d3_r_d_phi3_cylindrical(:,3)) &
       + d_r_d_phi_cylindrical(:,3) * (d2_r_d_phi2_cylindrical(:,1) * d3_r_d_phi3_cylindrical(:,2) - d2_r_d_phi2_cylindrical(:,2) * d3_r_d_phi3_cylindrical(:,1)))
    
  torsion_denominator = 0 &
       + (d_r_d_phi_cylindrical(:,2) * d2_r_d_phi2_cylindrical(:,3) - d_r_d_phi_cylindrical(:,3) * d2_r_d_phi2_cylindrical(:,2)) ** 2 &
       + (d_r_d_phi_cylindrical(:,3) * d2_r_d_phi2_cylindrical(:,1) - d_r_d_phi_cylindrical(:,1) * d2_r_d_phi2_cylindrical(:,3)) ** 2 &
       + (d_r_d_phi_cylindrical(:,1) * d2_r_d_phi2_cylindrical(:,2) - d_r_d_phi_cylindrical(:,2) * d2_r_d_phi2_cylindrical(:,1)) ** 2
  
  torsion = torsion_numerator / torsion_denominator

  sinangle = sin(phi)
  cosangle = cos(phi)
    
  ! Compute the Cartesian components of the Frenet unit vectors:
  tangent_Cartesian(:,1) = tangent_cylindrical(:,1) * cosangle - tangent_cylindrical(:,2) * sinangle
  tangent_Cartesian(:,2) = tangent_cylindrical(:,1) * sinangle + tangent_cylindrical(:,2) * cosangle
  tangent_Cartesian(:,3) = tangent_cylindrical(:,3)
    
  normal_Cartesian(:,1) = normal_cylindrical(:,1) * cosangle - normal_cylindrical(:,2) * sinangle
  normal_Cartesian(:,2) = normal_cylindrical(:,1) * sinangle + normal_cylindrical(:,2) * cosangle
  normal_Cartesian(:,3) = normal_cylindrical(:,3)
    
  binormal_Cartesian(:,1) = binormal_cylindrical(:,1) * cosangle - binormal_cylindrical(:,2) * sinangle
  binormal_Cartesian(:,2) = binormal_cylindrical(:,1) * sinangle + binormal_cylindrical(:,2) * cosangle
  binormal_Cartesian(:,3) = binormal_cylindrical(:,3)

  B1Squared_over_curvatureSquared = (B1s_over_B0*B1s_over_B0 + B1c_over_B0*B1c_over_B0) / (curvature * curvature)

  if (allocated(d_d_zeta)) deallocate(d_d_zeta)
  allocate(d_d_zeta(N_phi, N_phi))
  do j=1,N_phi
     d_d_zeta(j,:) = d_d_phi(j,:) / (B0_over_abs_G0 * d_l_d_phi(j))
  end do

  X1s = B1s_over_B0 / curvature
  X1c = B1c_over_B0 / curvature

  deallocate(sinangle, cosangle)
  deallocate(d2_l_d_phi2, torsion_numerator, torsion_denominator)
  deallocate(d_r_d_phi_cylindrical, d2_r_d_phi2_cylindrical, d3_r_d_phi3_cylindrical)


end subroutine quasisymmetry_init_axis
