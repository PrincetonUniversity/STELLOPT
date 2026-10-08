program test_mpi_params_serial
    use mpi_params, only: mpi_calc_myrange, mpi_stel_abort
    implicit none
    integer :: first, last, communicator

    communicator = 0
    call mpi_calc_myrange(communicator, -2, 7, first, last)
    if (first /= -2 .or. last /= 7) error stop 'Incorrect serial work interval'
    call mpi_calc_myrange(communicator, 3, 2, first, last)
    if (first /= 3 .or. last /= 2) error stop 'Empty serial interval changed'
    call mpi_stel_abort(42)
    print '(a)', 'PASS: serial interval and non-MPI diagnostic'
end program test_mpi_params_serial
