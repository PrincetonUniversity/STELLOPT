! Link the unmodified serial symforce routine with only its required globals.
module vmec_main
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    integer, parameter :: ns = 3, nzeta = 1, ntheta1 = 8
    integer, parameter :: ntheta2 = 5, ntheta3 = 8
    real(dp), parameter :: cp5 = 0.5_dp
    logical :: lthreed = .false.
    integer :: ireflect(3) = [1, 2, 3]
end module

module realspace
    implicit none
    integer :: ireflect(3) = [1, 2, 3]
end module

module parallel_include_module
    implicit none
end module

module timer_sub
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp) :: tforon = 0, tforoff = 0, s_symforces_time = 0, timer(1) = 0
    integer, parameter :: tfor = 1
end module

subroutine second0(t)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(out) :: t
    call cpu_time(t)
end subroutine

program serial_axis_oracle
    use vmec_main, only: dp, ns, ntheta1, ntheta2
    implicit none
    real(dp) :: force(3, 8, 0:1, 10), scratch(3, 8, 0:1, 10)
    real(dp) :: r1(24, 0:1), z1(24, 0:1), initial_r(24, 0:1), initial_z(24, 0:1)
    real(dp) :: theta, error
    real(dp) :: force_after(3, 8, 0:1, 10)
    real(dp), allocatable :: axis_r_save(:,:), axis_z_save(:,:)
    integer :: iter2 = 1
    logical :: lmove_axis = .true.
    character(len=32) :: mode
    integer :: i, j, flat
    call get_command_argument(1, mode)
    if (mode == "later_iteration") iter2 = 2
    force = 0
    scratch = 0
    r1 = 0
    z1 = 0
    do i = 1, ntheta1
        theta = 2*acos(-1.0_dp)*real(i - 1, dp)/ntheta1
        force(:, i, 0, 1) = -7000.0_dp*sin(theta)
        force(:, i, 0, 4) = -4000.0_dp*cos(theta)
        do j = 1, ns
            flat = ns*(i - 1) + j
            r1(flat, 0) = 6.2_dp + 0.62_dp*cos(theta)
            z1(flat, 0) = 0.3_dp + 0.62_dp*sin(theta)
        end do
    end do
    initial_r = r1
    initial_z = z1
    include "serial_axis_save.inc"
    call symforce(force(:, :, :, 1), force(:, :, :, 2), force(:, :, :, 3), &
        force(:, :, :, 4), force(:, :, :, 5), force(:, :, :, 6), &
        force(:, :, :, 7), force(:, :, :, 8), force(:, :, :, 9), &
        force(:, :, :, 10), r1, scratch(:, :, :, 2), scratch(:, :, :, 3), &
        z1, scratch(:, :, :, 5), scratch(:, :, :, 6), scratch(:, :, :, 7), &
        scratch(:, :, :, 8), scratch(:, :, :, 9), scratch(:, :, :, 10))
    error = 0
    do i = 1, ntheta2
        theta = 2*acos(-1.0_dp)*real(i - 1, dp)/ntheta1
        do j = 1, ns
            flat = ns*(i - 1) + j
            error = max(error, abs(r1(flat, 0) + 7000.0_dp*sin(theta)))
            error = max(error, abs(z1(flat, 0) + 4000.0_dp*cos(theta)))
        end do
    end do
    if (error > 1e-10_dp) error stop 'Serial native force parity differs from oracle'
    if (minval(r1(:ns*ntheta2, 0)) >= 0) error stop 'Serial overwrite not reproduced'
    if (minval(initial_r(:, 0)) <= 5) error stop 'Invalid initial geometry fixture'
    print *, 'Serial native force parity maximum error:', error
    print *, 'Initial physical R range:', minval(initial_r(:, 0)), maxval(initial_r(:, 0))
    print *, 'Written force range:', minval(r1(:ns*ntheta2, 0))
    force_after = force
    include 'serial_axis_restore.inc'
    if (iter2 == 1) then
        if (any(r1 /= initial_r)) error stop 'Axis scan would receive corrupted geometry'
        if (any(z1 /= initial_z)) error stop 'Axis scan would receive corrupted geometry'
    else
        if (allocated(axis_r_save)) error stop 'Later iteration saved axis geometry'
        if (minval(r1(:ns*ntheta2, 0)) >= 0) error stop 'Later force scratch changed'
    end if
    if (any(force /= force_after)) error stop 'Restoration changed force arrays'
    print *, 'Serial geometry branch and unchanged native force arrays verified.'
end program
