subroutine quasisymmetry_solve

  use quasisymmetry_variables

  implicit none

  integer :: iteration, j_line_search
  real(dp) :: residual_norm, last_residual_norm, initial_residual_norm, step_scale
  real(dp), dimension(:), allocatable :: state, state0

  ! Variables needed by LAPACK:                                                                                            
  integer :: INFO
  integer, dimension(:), allocatable :: IPIV

  allocate(state(matrix_size))
  allocate(state0(matrix_size))
  allocate(IPIV(matrix_size))

  ! Initialize state
  sigma = sigma_initial ! sigma_initial only needs to fix the first element of sigma, but let's set all the other elements of sigma to this same value as an initial guess.
  iota = 0

  state(1) = iota
  state(2:N_phi) = sigma(2:N_phi)

  call quasisymmetry_residual()
  initial_residual_norm = sqrt(sum(residual * residual))
  residual_norm = initial_residual_norm
  !print "(a,es10.3)","                 Initial relatiresidual norm:"

  ! Here is the main Newton iteration:
  Newton_tolerance_achieved = .false.
  Newton: do iteration = 1, N_iterations
     last_residual_norm = residual_norm
     !if (residual_norm / initial_residual_norm < Newton_tolerance) then
     if (residual_norm < Newton_tolerance) then
        Newton_tolerance_achieved = .true.
        exit Newton
     end if

     call quasisymmetry_Jacobian()

     state0 = state
     if (verbose) print "(a,i3)","  Newton iteration ",iteration
     ! We will use the LAPACK subroutine DGESV to solve a general (asymmetric) linear system
     ! step_direction = - matrix \ residual
     step_direction = -residual ! Note that LAPACK will over-write step_direction with the solution, and over-write Jacobian with the LU factorization.
     call DGESV(matrix_size, 1, Jacobian, matrix_size, IPIV, step_direction, matrix_size, INFO)
     if (INFO /= 0) then
        print *, "Error in LAPACK call DGESV: info = ", INFO
        stop
     end if

     step_scale = 1
     line_search: do j_line_search = 1, N_line_search
        state = state0 + step_scale * step_direction

        sigma = state ! This over-writes sigma(1) with the wrong value. Now fix it:
        sigma(1) = sigma_initial
        iota = state(1)
        call quasisymmetry_residual()
        residual_norm = sqrt(sum(residual * residual))
        !if (verbose) print "(a,i3,a,es10.3,a,es23.15)","    Line search step",j_line_search,"  Relative residual L2 norm:",residual_norm / initial_residual_norm,"  iota:",iota
        if (verbose) print "(a,i3,a,es10.3,a,es23.15)","    Line search step",j_line_search,"  Residual L2 norm:",residual_norm,"  iota:",iota
        if (residual_norm < last_residual_norm) exit line_search

        step_scale = step_scale / 2
     end do line_search

     if (residual_norm > last_residual_norm) then
        if (verbose) print *,"Line search failed to reduce residual."
        exit Newton
     end if
  end do Newton
  ! End of Newton solve.
  ! Now compute quantities that are derived from the solution:


  !Y1s = sign_G * curvature * ( B1c_over_B0 + B1s_over_B0 * sigma) / (B1c_over_B0*B1c_over_B0 + B1s_over_B0*B1s_over_B0)
  !Y1c = sign_G * curvature * (-B1s_over_B0 + B1c_over_B0 * sigma) / (B1c_over_B0*B1c_over_B0 + B1s_over_B0*B1s_over_B0)
  Y1s = sign_G * sign_psi * curvature / eta_bar
  Y1c = sign_G * sign_psi * curvature * sigma / eta_bar

  call quasisymmetry_elongation()

  !print *,"elongation:"
  !print *,elongation
  if (verbose) then
     print "(a,es23.15,a,es23.15)", " Final iota:",iota,"  max elongation:",max_elongation
     print "(a,es23.15)", " Final sigma(0): ",sigma(1)
  end if

  deallocate(state, state0, IPIV)

  ! Print R1c in MATLAB format, for comparing with the matlab version:
!!$  print *,"R1c_fortran=[",R1c(1),";"
!!$  do iteration = 2,N_phi-1
!!$     print *,R1c(iteration),";"
!!$  end do
!!$  print *,R1c(N_phi),"];"

end subroutine quasisymmetry_solve
