program fluxes

  ! takes argument 1 as .fluxes output file path
  ! argument 2 as number of time steps in .fluxes file
  ! argument 3 as number of species in simulation
  ! argument 4 as starting time for time average
  ! writes the time average of the fluxes to screen

  implicit none

  integer :: iargc
  integer :: flxunit = 101
  integer :: it, nstep, nspec
  integer :: target_it
  real :: tstart = 0.0
  logical :: tstart_flag
  character (500) :: line
  character (500) :: flxfile
  real, dimension (:), allocatable :: time
  real, dimension (:,:), allocatable :: pflx, vflx, qflx
  real, dimension (:), allocatable :: pflxavg, vflxavg, qflxavg, pi_over_q
  
  call getarg (1,flxfile)
  call getarg (2,line)
  read (line,*) nstep
  call getarg (3,line)
  read (line,*) nspec
  if (iargc() > 3) then
     call getarg (4,line)
     read (line,*) tstart
  end if

  write (*,*) nstep, tstart, nspec, trim(flxfile)

  allocate (time(nstep))
  allocate (pflx(nspec,nstep))
  allocate (vflx(nspec,nstep))
  allocate (qflx(nspec,nstep))
  allocate (pflxavg(nspec))
  allocate (vflxavg(nspec))
  allocate (qflxavg(nspec))
  allocate (pi_over_q(nspec))

  target_it = 1
  tstart_flag = .true.

  open (unit=flxunit, file=trim(flxfile)//".fluxes")
  read (flxunit,*) line
  do it = 1, nstep
     read (flxunit,*) time(it), pflx(:,it), vflx(:,it), qflx(:,it)
     ! find the time index corresponding to the requested start time
     if (tstart_flag .and. time(it) > tstart) then
        target_it = it
        tstart_flag = .false.
     end if
  end do

  call time_average (time, target_it, pflx, pflxavg)
  call time_average (time, target_it, vflx, vflxavg)
  call time_average (time, target_it, qflx, qflxavg)
  call time_average (time, target_it, vflx/qflx, pi_over_q)

  close (flxunit)

  open (flxunit, file=trim(flxfile)//".fluxes_tavg")
  write (flxunit,*) 'pflx: ', pflxavg, 'vflx: ', vflxavg, 'qflx: ', qflxavg, 'vflx/qflx: ', pi_over_q, vflxavg/qflxavg
  close (flxunit)

  deallocate (time)
  deallocate (pflx, vflx, qflx)
  deallocate (pflxavg, vflxavg, qflxavg, pi_over_q)

contains

  subroutine time_average (t, it, flx, flxavg)

    implicit none

    real, dimension (:), intent (in) :: t
    integer, intent (in) :: it
    real, dimension (:,:), intent (in) :: flx
    real, dimension (:), intent (out) :: flxavg

    integer :: i, nt

    nt = size(t)

    flxavg = 0.5*(t(it+1)-t(it))*flx(:,it)
    do i = it+1, nt-1
       flxavg = flxavg + (t(i+1)-t(i))*flx(:,i)
    end do
    flxavg = flxavg + 0.5*(t(nt)-t(nt-1))*flx(:,nt)

    flxavg = flxavg/(t(nt)-t(it))

  end subroutine time_average

end program fluxes
