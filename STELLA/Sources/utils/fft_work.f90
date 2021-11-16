module fft_work

  use constants, only: kind_id

  implicit none

  include "fftw3.f"

  public :: fft_type, delete_fft
  public :: init_ccfftw, init_crfftw, init_rcfftw
  public :: FFT_BACKWARD, FFT_FORWARD

  private

  integer, parameter :: FFT_BACKWARD = FFTW_BACKWARD
  integer, parameter :: FFT_FORWARD  = FFTW_FORWARD

  type :: fft_type
     integer :: n, is, type
     integer (kind_id) :: plan
     real :: scale
  end type fft_type

contains

  subroutine init_ccfftw (fft, is, n, data_in, data_out)

    use mp, only: mp_abort

    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, n
    complex, dimension(:), intent(inout) :: data_in, data_out

    integer :: j
    
    fft%n = n
    fft%is = is
    fft%scale = 1./real(n)
    if (is > 0) fft%scale = 1.
    fft%type = 1

    j = FFTW_PATIENT
    call dfftw_plan_dft_1d(fft%plan, n, data_in, data_out, is, j)

  end subroutine init_ccfftw

  subroutine init_crfftw (fft, is, n, data_in, data_out)

    use mp, only: mp_abort

    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, n
    complex, dimension(:), intent(in out) :: data_in
    real, dimension (:), intent (in out) :: data_out

    integer :: j
    
    fft%n = n
    fft%is = is
    fft%scale = 1./real(n)
    if (is > 0) fft%scale = 1.
    fft%type = 1
    
    j = FFTW_PATIENT
    call dfftw_plan_dft_c2r_1d(fft%plan, n, data_in, data_out, j)

  end subroutine init_crfftw

  subroutine init_rcfftw (fft, is, n, data_in, data_out)

    use mp, only: mp_abort

    type (fft_type), intent (out) :: fft
    integer, intent (in) :: is, n
    real, dimension(:), intent(in out) :: data_in
    complex, dimension (:), intent (in out) :: data_out

    integer :: j
    
    fft%n = n
    fft%is = is
    fft%scale = 1./real(n)
    if (is > 0) fft%scale = 1.
    fft%type = 1
    
    j = FFTW_PATIENT
    call dfftw_plan_dft_r2c_1d(fft%plan, n, data_in, data_out, j)

  end subroutine init_rcfftw

  subroutine delete_fft(fft)
    
    type (fft_type), intent (in out) :: fft

    call dfftw_destroy_plan(fft%plan)

  end subroutine delete_fft

end module fft_work
