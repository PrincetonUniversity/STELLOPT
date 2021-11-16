module text_options
  implicit none
  private
  public :: text_option
  public :: get_option_value

  integer, parameter :: maxlen = 30

  type :: text_option
     character(maxlen) :: name
     integer :: value
  end type text_option

contains

  subroutine get_option_value (selection, options, value, &
       error_unit, selection_name, stop_on_error)

    use mp, only: mp_abort

    implicit none
    character(*), intent (in) :: selection
    type (text_option), dimension (:), intent (in) :: options
    integer, intent (in out) :: value
    integer, intent (in), optional :: error_unit
    character(*), intent (in), optional :: selection_name
    logical, intent (in), optional :: stop_on_error
    integer :: i, l, n_partial_matches, v_partial_match
    integer :: err
    logical :: local_stop

    local_stop = .false.
    if (present(stop_on_error)) local_stop = stop_on_error

    do i = 1, size(options)
       if (trim(selection) == trim(options(i)%name)) then
          value = options(i)%value
          return
       end if
    end do

    ! look for partial matches
    l = len_trim(selection)
    n_partial_matches = 0
    do i = 1, size(options)
       if (l < len_trim(options(i)%name)) then
          if (trim(selection) == options(i)%name(1:l)) then
             n_partial_matches = n_partial_matches + 1
             v_partial_match = options(i)%value
          end if
       end if
    end do

    if (n_partial_matches == 1) then
       value = v_partial_match
       return
    end if

    if (present(error_unit)) then
       err = error_unit
    else
       err = 6
    end if

    if (n_partial_matches == 0) then
       if (present(selection_name)) then
          write (unit=err, fmt="('Invalid selection for ', a, ': ', a)") &
               trim(selection_name), trim(selection)
       else
          write (unit=err, fmt="('Invalid selection: ',a)") trim(selection)
       end if
       write (unit=err, fmt="('Valid selections are:')")
       do i = 1, size(options)
          write (unit=err, fmt="(3x,a)") trim(options(i)%name)
       end do
    else
       if (present(selection_name)) then
          write (unit=err, fmt="('Ambiguous selection for ', a, ': ', a)") &
               trim(selection_name), trim(selection)
       else
          write (unit=err, fmt="('Ambiguous selection: ',a)") trim(selection)
       end if
       write (unit=err, fmt="('Matching selections are:')")
       do i = 1, size(options)
          if (l < len_trim(options(i)%name)) then
             if (trim(selection) == options(i)%name(1:l)) then
                write (unit=err, fmt="(3x,a)") trim(options(i)%name)
             end if
          end if
       end do
    end if

    if (local_stop) then
      call mp_abort ('STOP error in get_option_value')
    endif

    write (unit=err, fmt="('Continuing with default selection...')")
  end subroutine get_option_value

end module text_options
