module stella_time

  implicit none

  private

  real :: code_dt 
  real :: cfl_dt = -1.
  real :: code_dt_min, code_dt_max

  ! added May 18, 2009 to take care of problems
  ! in exb_shear calculation after change in time step size
  real :: code_dt_old = 0.
  real :: code_time = 0.

  public :: code_dt, update_time, code_dt_old
  public :: code_time
  public :: save_dt_min, save_dt, save_dt_cfl, write_dt
  public :: init_tstart, init_delt
  public :: cfl_dt
  public :: code_dt_min, code_dt_max

contains

  subroutine init_tstart (tstart)

    real, intent (in) :: tstart

    code_time = tstart

  end subroutine init_tstart

  subroutine init_delt (delt)
    real, intent (in) :: delt

    code_dt = delt
    ! do not allow code_dt to increase beyond input value
    code_dt_max = code_dt

  end subroutine init_delt

  subroutine update_time
! MAB+CMR, 21/5/09: set code_dt_old to code_dt BEFORE any changes in timestep
    code_dt_old = code_dt
    code_time = code_time + code_dt

  end subroutine update_time

  subroutine save_dt_cfl (delt_cfl)

    real, intent (in) :: delt_cfl

    cfl_dt = delt_cfl
       
  end subroutine save_dt_cfl

  subroutine save_dt_min (dt_min)

    real, intent (in) :: dt_min

    code_dt_min = dt_min

  end subroutine save_dt_min

  subroutine save_dt(delt)
    
    real, intent (in) :: delt

    code_dt = delt
    
  end subroutine save_dt

  subroutine write_dt

    if (cfl_dt > 0. .and. cfl_dt < 1.e7) &
         write (*,*) 'TIME STEP:'
         write (*,'(A12, ES10.2E2)') "   cfl_dt:"//REPEAT(' ',50),cfl_dt
         write (*,'(A12, ES10.2E2)') "   code_dt:"//REPEAT(' ',50),code_dt 
       
  end subroutine write_dt

end module stella_time
