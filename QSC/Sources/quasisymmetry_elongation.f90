subroutine quasisymmetry_elongation

  use quasisymmetry_variables

  implicit none

  real(dp), dimension(:), allocatable :: p, q
  real(dp), dimension(:), allocatable :: elongation_sin, elongation_cos, elongation_reconstructed
  integer :: n, index_of_max
  real(dp) :: phi_of_max_elongation, fmin_tolerance, quasisymmetry_fmin, maxval_elongation

  allocate(q(N_phi))
  allocate(p(N_phi))

  ! See my note 20180329-03 for derivation of the formula below for elongation:
  p = R1s*R1s + R1c*R1c + Z1s*Z1s + Z1c*Z1c
  q = R1s*Z1c - R1c*Z1s
  !elongation = 2*abs(q) / (p - sqrt(p*p-4*q*q)) ! This version suffers precision loss for large elongation.
  elongation = (p + sqrt(p*p-4*q*q))/(2*abs(q)) ! This version is more stable numerically.

  ! Search for maximum using Fourier interpolation...
  index_of_max = maxloc(elongation,1)
  maxval_elongation = elongation(index_of_max)
  if (maxval_elongation > max_precise_elongation) then
     max_elongation = maxval_elongation
     if (verbose) then
        print *,"maxval(elongation):      ",maxval_elongation
        print "(a)"," elongation > max_precise_elongation so skipping precise solve."
     end if
     return
  end if


  ! In preparation for searching for max elongation, Fourier transform the elongation.
  allocate(elongation_sin((N_phi+1)/2))
  allocate(elongation_cos((N_phi+1)/2))
  allocate(elongation_reconstructed(N_phi))
  elongation_sin(1) = 0
  elongation_cos(1) = sum(elongation) / N_phi
  do n=1,((N_phi-1)/2)
     elongation_sin(n+1) = sum(elongation * sin_n_phi(:,n+1)) * 2 / (N_phi)
     elongation_cos(n+1) = sum(elongation * cos_n_phi(:,n+1)) * 2 / (N_phi)
  end do

!!$  elongation_reconstructed = elongation_cos(1)
!!$  do n = 1,((N_phi-1)/2)
!!$     elongation_reconstructed = elongation_reconstructed + elongation_sin(n+1) * sin_n_phi(:,n+1) + elongation_cos(n+1) * cos_n_phi(:,n+1)
!!$  end do
!!$
!!$  print *,"elongation:"
!!$  print *,elongation
!!$  print *,"elongation_reconstructed:"
!!$  print *,elongation_reconstructed
!!$
!!$  print *,"elongation_sin:",elongation_sin
!!$  print *,"elongation_cos:",elongation_cos

!!$  print *,"index_of_max:",index_of_max
!!$  print *,"d_phi:",d_phi,"(index_of_max-1)*d_phi:",(index_of_max-1)*d_phi
!!$  print *,"Interpolated elongation at index_of_max:",-minus_elongation((index_of_max-1)*d_phi)

  fmin_tolerance = 0
  phi_of_max_elongation =  quasisymmetry_fmin((index_of_max-2)*d_phi, (index_of_max)*d_phi, &
       minus_elongation, fmin_tolerance)

  max_elongation = -minus_elongation(phi_of_max_elongation)
  if (verbose) then
     print *,"maxval(elongation):      ",maxval_elongation
     print *,"max elongation from fmin:",max_elongation
  end if
  if (maxval_elongation > max_elongation * (1 + 1.0d-10)) then
     print *,"Error! Something went wrong with the max_elongation search."
     if (maxval_elongation < max_precise_elongation) stop ! If the elongation is larger than this, we don't care about the solution much, so don't bother aborting.
  end if

  deallocate(p,q)
  deallocate(elongation_sin, elongation_cos, elongation_reconstructed)

contains

  real(dp) function minus_elongation(this_phi)
    
    implicit none
    
    real(dp) :: f, this_phi
    integer :: nn
    
    f = elongation_cos(1)

    do nn = 1,((N_phi-1)/2)
       f = f + elongation_sin(nn+1) * sin(nn*nfp*this_phi) + elongation_cos(nn+1) * cos(nn*nfp*this_phi)
    end do
    minus_elongation = -f
    return
    
  end function minus_elongation

end subroutine quasisymmetry_elongation


