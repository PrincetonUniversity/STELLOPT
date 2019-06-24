subroutine quasisymmetry_max_curvature

  use quasisymmetry_variables

  implicit none

  integer :: index_of_max
  real(dp) :: phi_of_max_curvature, fmin_tolerance, quasisymmetry_fmin, maxval_curvature

  already_found_max_curvature = .true.

  index_of_max = maxloc(curvature,1)
  maxval_curvature = curvature(index_of_max)

  fmin_tolerance = 1.0d-10
  phi_of_max_curvature =  quasisymmetry_fmin((index_of_max-2)*d_phi, (index_of_max)*d_phi, &
       minus_curvature, fmin_tolerance)

  max_curvature = -minus_curvature(phi_of_max_curvature)
  if (verbose) then
     print *,"maxval(curvature):      ",maxval_curvature
     print *,"max curvature from fmin:",max_curvature
  end if
  if (maxval_curvature > max_curvature * (1 + 1.0d-10)) then
     print *,"Error! Something went wrong with the max_curvature search on proc",mpi_rank,". maxval_curvature=",maxval_curvature,", max_curvature=",max_curvature,", curvature=",curvature
     max_curvature = maxval_curvature ! Use the bigger number, going forward.
     !if (maxval_curvature < 1000) stop ! If the curvature is larger than this, we don't care about the solution much, so don't bother aborting.
  end if

contains

  real(dp) function minus_curvature(this_phi)
    
    implicit none
    
    real(dp) :: this_phi, cosangle, sinangle
    real(dp) :: local_d_tangent_d_l_cylindrical(3), local_d_r_d_phi_cylindrical(3), local_d2_r_d_phi2_cylindrical(3)
    real(dp) :: local_R0, local_R0p, local_Z0p, local_R0pp, local_Z0pp, local_d_l_d_phi, local_d2_l_d_phi2
    integer :: n, j

    local_R0 = R0c(1)
    local_R0p = 0
    local_Z0p = 0
    local_R0pp = 0
    local_Z0pp = 0

    do n = 1, axis_nmax
       sinangle = sin(n*nfp*this_phi)
       cosangle = cos(n*nfp*this_phi)
       
       local_R0 = local_R0 + R0c(n+1) * cosangle + R0s(n+1) * sinangle
       local_R0p = local_R0p + R0c(n+1) * (-n*nfp)*sinangle + R0s(n+1) * (n*nfp)*cosangle
       local_Z0p = local_Z0p + Z0c(n+1) * (-n*nfp)*sinangle + Z0s(n+1) * (n*nfp)*cosangle
       local_R0pp = local_R0pp + R0c(n+1) * (-n*nfp*n*nfp)*cosangle + R0s(n+1) * (-n*nfp*n*nfp)*sinangle
       local_Z0pp = local_Z0pp + Z0c(n+1) * (-n*nfp*n*nfp)*cosangle + Z0s(n+1) * (-n*nfp*n*nfp)*sinangle

    end do

    local_d_l_d_phi = sqrt(local_R0 * local_R0 + local_R0p * local_R0p + local_Z0p * local_Z0p)
    local_d2_l_d_phi2 = (local_R0 * local_R0p + local_R0p * local_R0pp + local_Z0p * local_Z0pp) / local_d_l_d_phi

    local_d_r_d_phi_cylindrical(1) = local_R0p
    local_d_r_d_phi_cylindrical(2) = local_R0
    local_d_r_d_phi_cylindrical(3) = local_Z0p

    local_d2_r_d_phi2_cylindrical(1) = local_R0pp - local_R0
    local_d2_r_d_phi2_cylindrical(2) = 2 * local_R0p
    local_d2_r_d_phi2_cylindrical(3) = local_Z0pp

    do j = 1,3
       local_d_tangent_d_l_cylindrical(j) = (-local_d_r_d_phi_cylindrical(j) * local_d2_l_d_phi2 / local_d_l_d_phi &
            + local_d2_r_d_phi2_cylindrical(j)) / (local_d_l_d_phi * local_d_l_d_phi)
    end do
    
    minus_curvature = -sqrt( &
         local_d_tangent_d_l_cylindrical(1) * local_d_tangent_d_l_cylindrical(1) + &
         local_d_tangent_d_l_cylindrical(2) * local_d_tangent_d_l_cylindrical(2) + &
         local_d_tangent_d_l_cylindrical(3) * local_d_tangent_d_l_cylindrical(3))

    return
    
  end function minus_curvature

end subroutine quasisymmetry_max_curvature


