! Compile with the unmodified production symforce_par and extracted repair blocks.
module vmec_main
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    integer, parameter :: nzeta = 1, ntheta1 = 8, ntheta2 = 5, ntheta3 = 8
    integer, parameter :: ns = 3, nznt = 8
    real(dp), parameter :: cp5 = 0.5_dp
    logical :: lthreed = .false.
end module

module realspace
    implicit none
    integer :: ireflect_par(1) = [1]
end module

module parallel_include_module
    implicit none
    integer :: t1lglob = 1, t1rglob = 3
end module

module timer_sub
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp) :: tforon = 0, tforoff = 0, symforces_time = 0, timer(1) = 0
    integer, parameter :: tfor = 1
end module

subroutine second0(t)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(out) :: t
    call cpu_time(t)
end subroutine

program axis_storage_oracle
    use vmec_main, only: dp, nznt, ns, ntheta1, ntheta2
    implicit none
    real(dp) :: force(1, 8, 3, 0:1, 10), scratch(1, 8, 3, 0:1, 10)
    real(dp) :: pr1(8, 3, 0:1), pz1(8, 3, 0:1)
    real(dp) :: original_r(8, 3, 0:1), original_z(8, 3, 0:1)
    real(dp) :: force_after(1, 8, 3, 0:1, 10)
    real(dp), allocatable :: axis_r_save(:, :, :), axis_z_save(:, :, :)
    real(dp) :: theta, error
    integer :: i, iter2 = 1
    logical :: lmove_axis = .true.
    character(len=16) :: shape
    call get_command_argument(1, shape)
    force = 0
    scratch = 0
    pr1 = 0
    pz1 = 0
    do i = 1, ntheta1
        theta = 2*acos(-1.0_dp)*real(i - 1, dp)/ntheta1
        force(1, i, :, 0, 1) = -7000.0_dp*sin(theta)
        force(1, i, :, 0, 4) = -4000.0_dp*cos(theta)
        pr1(i, :, 0) = 6.2_dp + 0.62_dp*cos(theta)
        pz1(i, :, 0) = 0.3_dp + 0.62_dp*sin(theta)
        if (shape == "shaped") then
            pr1(i, :, 0) = pr1(i, :, 0) + 0.09_dp*cos(2*theta)
            pz1(i, :, 0) = 0.3_dp + 1.05_dp*sin(theta) + 0.06_dp*cos(2*theta)
        end if
    end do
    original_r = pr1
    original_z = pz1
    include 'official_axis_save.inc'
    ! Actual funct3d_par alias: pr1 and pz1 become ara and aza force outputs.
    call symforce_par(force(:, :, :, :, 1), force(:, :, :, :, 2), &
        force(:, :, :, :, 3), force(:, :, :, :, 4), force(:, :, :, :, 5), &
        force(:, :, :, :, 6), force(:, :, :, :, 7), force(:, :, :, :, 8), &
        force(:, :, :, :, 9), force(:, :, :, :, 10), pr1, &
        scratch(:, :, :, :, 2), scratch(:, :, :, :, 3), pz1, &
        scratch(:, :, :, :, 5), scratch(:, :, :, :, 6), &
        scratch(:, :, :, :, 7), scratch(:, :, :, :, 8), &
        scratch(:, :, :, :, 9), scratch(:, :, :, :, 10))
    error = 0
    do i = 1, ntheta2
        theta = 2*acos(-1.0_dp)*real(i - 1, dp)/ntheta1
        error = max(error, abs(pr1(i, 3, 0) + 7000.0_dp*sin(theta)))
        error = max(error, abs(pz1(i, 3, 0) + 4000.0_dp*cos(theta)))
    end do
    if (error > 1e-10_dp) error stop 'Native sine/cosine force parity oracle failed'
    if (minval(pr1(:ntheta2, 3, 0)) >= 0) error stop 'Force scratch overwrite absent'
    print *, 'Native force parity maximum error:', error
    print *, 'Physical R before scratch:', &
        minval(original_r(:, 3, 0)), maxval(original_r(:, 3, 0))
    print *, 'Force values after scratch:', minval(pr1(:ntheta2, 3, 0))
    force_after = force
    include 'official_axis_restore.inc'
    if (any(pr1 /= original_r)) error stop 'Axis scan would receive corrupted geometry'
    if (any(pz1 /= original_z)) error stop 'Axis scan would receive corrupted geometry'
    if (any(force /= force_after)) error stop 'Restoration changed native force arrays'
    print *, 'Native R/Z restored exactly; force arrays preserved.'
end program
