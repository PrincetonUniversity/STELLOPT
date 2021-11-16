module physics_parameters

  implicit none

  public :: init_physics_parameters
  public :: finish_physics_parameters
  public :: beta, zeff, tite, nine, rhostar, vnew_ref
  public :: g_exb, g_exbfac, omprimfac

  private

  real :: beta, zeff, tite, nine, rhostar, irhostar, vnew_ref
  real :: g_exb, g_exbfac, omprimfac
  logical :: initialized = .false.

contains

  subroutine init_physics_parameters

    implicit none

    if (initialized) return
    initialized = .true.

    call read_parameters

  end subroutine init_physics_parameters

  subroutine read_parameters

    use file_utils, only: input_unit_exist
    use mp, only: proc0, broadcast

    implicit none

    integer :: in_file
    logical :: rpexist
    
    namelist /parameters/ beta, zeff, tite, nine, rhostar, vnew_ref &
                          ,g_exb, g_exbfac, omprimfac, irhostar

    if (proc0) then
       beta = 0.0
       vnew_ref = -1.0 ! various input options will override this value if it is negative
       rhostar  = -1.0 ! = m_ref * vt_ref / (e * B_ref * a_ref), with refs in SI
       irhostar = -1.0
       zeff = 1.0
       tite = 1.0
       nine = 1.0

       g_exb    = 0.0
       g_exbfac = 1.0
       omprimfac = 1.0

       in_file = input_unit_exist("parameters", rpexist)
       if (rpexist) read (unit=in_file,nml=parameters)

       if (irhostar.gt.0) rhostar = 1./irhostar
    end if

    call broadcast (beta)
    call broadcast (vnew_ref)
    call broadcast (zeff)
    call broadcast (rhostar)
    call broadcast (tite)
    call broadcast (nine)
    call broadcast (g_exb)
    call broadcast (g_exbfac)
    call broadcast (omprimfac)

  end subroutine read_parameters

  subroutine finish_physics_parameters

    implicit none

    initialized = .false.

  end subroutine finish_physics_parameters

end module physics_parameters
