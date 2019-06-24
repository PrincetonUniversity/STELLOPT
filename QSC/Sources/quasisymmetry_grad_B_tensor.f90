subroutine quasisymmetry_grad_B_tensor

  use quasisymmetry_variables

  implicit none

  real(dp), dimension(:,:,:), allocatable :: grad_B_tensor
  real(dp), dimension(:), allocatable :: div_B, should_be_curvature
  real(dp), dimension(:,:), allocatable :: should_be_curvature_times_normal, should_be_normal_Cartesian
  integer :: j, k
  real(dp) :: iota_N

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  allocate(grad_B_tensor(N_phi,3,3))

  if (allocated(modBinv_sqrt_half_grad_B_colon_grad_B)) deallocate(modBinv_sqrt_half_grad_B_colon_grad_B)
  allocate(modBinv_sqrt_half_grad_B_colon_grad_B(N_phi))

  iota_N = iota + axis_helicity*nfp

  ! The formula below is derived in "20181018-02 Gradient of B from near-axis expansion- v2.docx"

  modBinv_sqrt_half_grad_B_colon_grad_B = 0
  do j = 1, 3
     do k = 1, 3
        grad_B_tensor(:,k,j) = sign_psi * B0 * B0_over_abs_G0 * ( &
             (sign_G*sign_psi*abs_G0_over_B0*curvature*tangent_Cartesian(:,j) &
             + (d_X1c_d_zeta*Y1s + iota_N*X1c*Y1c) * normal_Cartesian(:,j) &
             + (d_Y1c_d_zeta*Y1s - d_Y1s_d_zeta*Y1c + sign_G*sign_psi*abs_G0_over_B0*torsion + iota_N*(Y1s*Y1s + Y1c*Y1c)) * binormal_Cartesian(:,j) &
             ) * normal_Cartesian(:,k) &
             + ((-sign_G*sign_psi*abs_G0_over_B0*torsion - iota_N*X1c*X1c) * normal_Cartesian(:,j) &
             + (X1c*d_Y1s_d_zeta - iota_N*X1c*Y1c) * binormal_Cartesian(:,j) &
             ) * binormal_Cartesian(:,k)) &
             + sign_G*curvature*B0*normal_Cartesian(:,j)*tangent_Cartesian(:,k)

        modBinv_sqrt_half_grad_B_colon_grad_B = modBinv_sqrt_half_grad_B_colon_grad_B + grad_B_tensor(:,k,j) ** 2
     end do
  end do

  modBinv_sqrt_half_grad_B_colon_grad_B = sqrt((0.5d+0) * modBinv_sqrt_half_grad_B_colon_grad_B) / B0

  !if (.true.) then
  if (.false.) then
     ! Sanity tests

     print *,"modBinv_sqrt_half_grad_B_colon_grad_B:",modBinv_sqrt_half_grad_B_colon_grad_B

     allocate(div_B(N_phi))
     allocate(should_be_curvature(N_phi))
     allocate(should_be_curvature_times_normal(N_phi,3))
     allocate(should_be_normal_Cartesian(N_phi,3))
     div_B = grad_B_tensor(:,1,1) + grad_B_tensor(:,2,2) + grad_B_tensor(:,3,3)
     print *,"div_B:",div_B
     print *,"max(abs(div_B)):",maxval(abs(div_B))

     print *," "
     print *,"grad_B_tensor(:,1,2):",grad_B_tensor(:,1,2)
     print *,"grad_B_tensor(:,2,1):",grad_B_tensor(:,2,1)
     print *,"max(abs(grad_B_tensor(:,1,2) - grad_B_tensor(:,2,1))):",maxval(abs(grad_B_tensor(:,1,2)-grad_B_tensor(:,2,1)))
     print *," "
     print *,"grad_B_tensor(:,1,3):",grad_B_tensor(:,1,3)
     print *,"grad_B_tensor(:,3,1):",grad_B_tensor(:,3,1)
     print *,"max(abs(grad_B_tensor(:,1,3) - grad_B_tensor(:,3,1))):",maxval(abs(grad_B_tensor(:,1,3)-grad_B_tensor(:,3,1)))

     print *," "
     print *,"grad_B_tensor(:,2,3):",grad_B_tensor(:,2,3)
     print *,"grad_B_tensor(:,3,2):",grad_B_tensor(:,3,2)
     print *,"max(abs(grad_B_tensor(:,2,3) - grad_B_tensor(:,3,2))):",maxval(abs(grad_B_tensor(:,2,3)-grad_B_tensor(:,3,2)))

     ! kappa \vec{N} = (1/|B|) grad |B| = (1/B^2) (grad \vec{B}) dot \vec{B} = (sign_G/|B|) (grad \vec{B}) dot \vec{t}
     should_be_curvature_times_normal = 0
     do j = 1, 3
        do k = 1, 3
           should_be_curvature_times_normal(:,k) = should_be_curvature_times_normal(:,k) + (1/B0) * grad_B_tensor(:,k,j) * tangent_Cartesian(:,j)
        end do
     end do
     should_be_curvature = sqrt(should_be_curvature_times_normal(:,1) ** 2 + should_be_curvature_times_normal(:,2) ** 2 + should_be_curvature_times_normal(:,3) ** 2)
     do k = 1, 3
        should_be_normal_Cartesian(:,k) = should_be_curvature_times_normal(:,k) / should_be_curvature
     end do

     print *," "
     print *,"should_be_curvature:",should_be_curvature
     print *,"          curvature:",curvature
     print *,"max(abs(curvature - should_be_curvature)):",maxval(abs(curvature - should_be_curvature))
     print *," "
     print *,"max(abs(normal_Cartesian - should_be_normal_Cartesian)):",maxval(abs(normal_Cartesian - should_be_normal_Cartesian))
     print *," "

     deallocate(div_B, should_be_curvature, should_be_curvature_times_normal, should_be_normal_Cartesian)

  end if

  call quasisymmetry_max_modBinv_sqrt_half_grad_B_colon_grad_B()

  deallocate(grad_B_tensor)

end subroutine quasisymmetry_grad_B_tensor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

subroutine quasisymmetry_max_modBinv_sqrt_half_grad_B_colon_grad_B

  use quasisymmetry_variables

  implicit none

  real(dp), dimension(:), allocatable :: modBinv_sqrt_half_grad_B_colon_grad_B_sin, modBinv_sqrt_half_grad_B_colon_grad_B_cos, modBinv_sqrt_half_grad_B_colon_grad_B_reconstructed
  integer :: n, index_of_max
  real(dp) :: phi_of_max_modBinv_sqrt_half_grad_B_colon_grad_B, fmin_tolerance, quasisymmetry_fmin, maxval_modBinv_sqrt_half_grad_B_colon_grad_B


  ! Search for maximum using Fourier interpolation...
  index_of_max = maxloc(modBinv_sqrt_half_grad_B_colon_grad_B,1)
  maxval_modBinv_sqrt_half_grad_B_colon_grad_B = modBinv_sqrt_half_grad_B_colon_grad_B(index_of_max)
!!$  if (maxval_modBinv_sqrt_half_grad_B_colon_grad_B > max_precise_modBinv_sqrt_half_grad_B_colon_grad_B) then
!!$     max_modBinv_sqrt_half_grad_B_colon_grad_B = maxval_modBinv_sqrt_half_grad_B_colon_grad_B
!!$     if (verbose) then
!!$        print *,"maxval(modBinv_sqrt_half_grad_B_colon_grad_B):      ",maxval_modBinv_sqrt_half_grad_B_colon_grad_B
!!$        print "(a)"," modBinv_sqrt_half_grad_B_colon_grad_B > max_precise_modBinv_sqrt_half_grad_B_colon_grad_B so skipping precise solve."
!!$     end if
!!$     return
!!$  end if


  ! In preparation for searching for max modBinv_sqrt_half_grad_B_colon_grad_B, Fourier transform the modBinv_sqrt_half_grad_B_colon_grad_B.
  allocate(modBinv_sqrt_half_grad_B_colon_grad_B_sin((N_phi+1)/2))
  allocate(modBinv_sqrt_half_grad_B_colon_grad_B_cos((N_phi+1)/2))
  allocate(modBinv_sqrt_half_grad_B_colon_grad_B_reconstructed(N_phi))
  modBinv_sqrt_half_grad_B_colon_grad_B_sin(1) = 0
  modBinv_sqrt_half_grad_B_colon_grad_B_cos(1) = sum(modBinv_sqrt_half_grad_B_colon_grad_B) / N_phi
  do n=1,((N_phi-1)/2)
     modBinv_sqrt_half_grad_B_colon_grad_B_sin(n+1) = sum(modBinv_sqrt_half_grad_B_colon_grad_B * sin_n_phi(:,n+1)) * 2 / (N_phi)
     modBinv_sqrt_half_grad_B_colon_grad_B_cos(n+1) = sum(modBinv_sqrt_half_grad_B_colon_grad_B * cos_n_phi(:,n+1)) * 2 / (N_phi)
  end do

!!$  modBinv_sqrt_half_grad_B_colon_grad_B_reconstructed = modBinv_sqrt_half_grad_B_colon_grad_B_cos(1)
!!$  do n = 1,((N_phi-1)/2)
!!$     modBinv_sqrt_half_grad_B_colon_grad_B_reconstructed = modBinv_sqrt_half_grad_B_colon_grad_B_reconstructed + modBinv_sqrt_half_grad_B_colon_grad_B_sin(n+1) * sin_n_phi(:,n+1) + modBinv_sqrt_half_grad_B_colon_grad_B_cos(n+1) * cos_n_phi(:,n+1)
!!$  end do
!!$
!!$  print *,"modBinv_sqrt_half_grad_B_colon_grad_B:"
!!$  print *,modBinv_sqrt_half_grad_B_colon_grad_B
!!$  print *,"modBinv_sqrt_half_grad_B_colon_grad_B_reconstructed:"
!!$  print *,modBinv_sqrt_half_grad_B_colon_grad_B_reconstructed
!!$
!!$  print *,"modBinv_sqrt_half_grad_B_colon_grad_B_sin:",modBinv_sqrt_half_grad_B_colon_grad_B_sin
!!$  print *,"modBinv_sqrt_half_grad_B_colon_grad_B_cos:",modBinv_sqrt_half_grad_B_colon_grad_B_cos

!!$  print *,"index_of_max:",index_of_max
!!$  print *,"d_phi:",d_phi,"(index_of_max-1)*d_phi:",(index_of_max-1)*d_phi
!!$  print *,"Interpolated modBinv_sqrt_half_grad_B_colon_grad_B at index_of_max:",-minus_modBinv_sqrt_half_grad_B_colon_grad_B((index_of_max-1)*d_phi)

  fmin_tolerance = 0
  phi_of_max_modBinv_sqrt_half_grad_B_colon_grad_B =  quasisymmetry_fmin((index_of_max-2)*d_phi, (index_of_max)*d_phi, &
       minus_modBinv_sqrt_half_grad_B_colon_grad_B, fmin_tolerance)

  max_modBinv_sqrt_half_grad_B_colon_grad_B = -minus_modBinv_sqrt_half_grad_B_colon_grad_B(phi_of_max_modBinv_sqrt_half_grad_B_colon_grad_B)
  if (verbose) then
     print *,"maxval(modBinv_sqrt_half_grad_B_colon_grad_B):      ",maxval_modBinv_sqrt_half_grad_B_colon_grad_B
     print *,"max modBinv_sqrt_half_grad_B_colon_grad_B from fmin:",max_modBinv_sqrt_half_grad_B_colon_grad_B
  end if
  if (maxval_modBinv_sqrt_half_grad_B_colon_grad_B > max_modBinv_sqrt_half_grad_B_colon_grad_B * (1 + 1.0d-10)) then
     print *,"Error! Something went wrong with the max_modBinv_sqrt_half_grad_B_colon_grad_B search on proc",mpi_rank,". maxval_modBinv_sqrt_half_grad_B_colon_grad_B=",maxval_modBinv_sqrt_half_grad_B_colon_grad_B,", max_modBinv_sqrt_half_grad_B_colon_grad_B=",max_modBinv_sqrt_half_grad_B_colon_grad_B,", modBinv_sqrt_half_grad_B_colon_grad_B=",modBinv_sqrt_half_grad_B_colon_grad_B
     !if (maxval_modBinv_sqrt_half_grad_B_colon_grad_B < max_precise_modBinv_sqrt_half_grad_B_colon_grad_B) stop ! If the modBinv_sqrt_half_grad_B_colon_grad_B is larger than this, we don't care about the solution much, so don't bother aborting.
  end if

  deallocate(modBinv_sqrt_half_grad_B_colon_grad_B_sin, modBinv_sqrt_half_grad_B_colon_grad_B_cos, modBinv_sqrt_half_grad_B_colon_grad_B_reconstructed)

contains

  real(dp) function minus_modBinv_sqrt_half_grad_B_colon_grad_B(this_phi)
    
    implicit none
    
    real(dp) :: f, this_phi
    integer :: nn
    
    f = modBinv_sqrt_half_grad_B_colon_grad_B_cos(1)

    do nn = 1,((N_phi-1)/2)
       f = f + modBinv_sqrt_half_grad_B_colon_grad_B_sin(nn+1) * sin(nn*nfp*this_phi) + modBinv_sqrt_half_grad_B_colon_grad_B_cos(nn+1) * cos(nn*nfp*this_phi)
    end do
    minus_modBinv_sqrt_half_grad_B_colon_grad_B = -f
    return
    
  end function minus_modBinv_sqrt_half_grad_B_colon_grad_B

end subroutine quasisymmetry_max_modBinv_sqrt_half_grad_B_colon_grad_B


