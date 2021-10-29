# include "define.inc"

module system_fortran
  implicit none

  private

  public :: systemf

#ifdef ISO_C_BINDING
    interface
      subroutine call_system_c (command) bind(C, name='system')
        use, intrinsic :: iso_c_binding, only: c_char, c_int
        character(kind=c_char), intent(in) :: command(*)
     end subroutine call_system_c
   end interface  
#endif

  contains

  subroutine systemf (command)
#ifdef ISO_C_BINDING
    use, intrinsic :: iso_c_binding, only: c_null_char
    implicit none

    character (*), intent (in) :: command

    call call_system_c(command//c_null_char)

#else
# if FCOMPILER == _INTEL_
!   for system call with intel compiler
    use ifport, only: system
# endif
    
    implicit none

    character (*), intent (in) :: command
# if FCOMPILER == _CRAY_
    integer :: ierr
    ierr=system(command)
#endif
#endif /* ISO_C_BINDING */

  end subroutine systemf

end module system_fortran
