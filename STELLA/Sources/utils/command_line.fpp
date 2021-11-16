# include "define.inc"

module command_line

  ! (c) Copyright 1991 to 1998 by Michael A. Beer, William D. Dorland, 
  ! P. B. Snyder, Q. P. Liu, and Gregory W. Hammett. ALL RIGHTS RESERVED.
  !
  ! <doc>
  !  A wrapper module for handling command line arguments.
  !  This module provides subroutine cl_getarg and integer function cl_iargc.
  !  Most of the compilers have getarg and iargc as their extensions.
  !  If not, one can use POSIX pxfgetarg and ipxfargc.
  !
  !  Note that Fortran 2003 includes get_command_argument and 
  !  command_argument_count which will replace getarg and iargc.
  ! </doc>

  implicit none

  private

  public :: cl_getarg, cl_iargc

contains

  function cl_iargc()

    ! <doc>
    !  returns the number of arguments
    !  using intrinsic iargc or POSIX ipxfargc
    ! </doc>

# ifdef POSIX
# if FCOMPILER == _INTEL_
    use ifposix, only: iargc => ipxfargc
# endif
# else
# if FCOMPILER == _NAG_
    use f90_unix, only: iargc
# endif
# endif

    implicit none
    integer :: cl_iargc
# if ( POSIX == _NONE_ && FCOMPILER != _NAG_ )
    integer :: iargc
# endif

    cl_iargc = iargc()

  end function cl_iargc

  subroutine cl_getarg (k, arg, len, ierr)
    
    ! <doc>
    !  gets k-th argument string and its length
    !  using intrinsic getarg or POSIX pxfgetarg
    ! </doc>

# ifdef POSIX
    use ifposix, only: pxfgetarg
# else
# if FCOMPILER == _NAG_
    use f90_unix, only: getarg
# endif
# endif

    implicit none
    integer,           intent (in)  :: k
    character (len=*), intent (out) :: arg
    integer,           intent (out) :: len
    integer,           intent (out) :: ierr

# ifdef POSIX
    call pxfgetarg (k, arg, len, ierr)
# else
    call getarg (k, arg)
    len = len_trim(arg)
    ierr = 0
# endif

  end subroutine cl_getarg

end module command_line
