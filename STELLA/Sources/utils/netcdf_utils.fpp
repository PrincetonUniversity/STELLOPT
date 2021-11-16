# include "define.inc"

module netcdf_utils

# ifdef NETCDF
  use netcdf, only: NF90_FLOAT, NF90_DOUBLE
  use netcdf, only: NF90_NOWRITE, NF90_CLOBBER, NF90_NOERR
  use netcdf, only: nf90_strerror
  use netcdf, only: nf90_inquire_variable
  use netcdf, only: nf90_inquire_dimension
  use netcdf, only: nf90_open, nf90_close
  use netcdf, only: nf90_inq_varid
  use netcdf, only: NF90_INT
# endif
  implicit none

  public :: netcdf_error
  public :: get_netcdf_code_precision
  public :: check_netcdf_file_precision
  public :: netcdf_real, kind_nf, netcdf_int
  private

# ifdef NETCDF
  integer, parameter :: kind_nf = kind (NF90_NOERR)
# else
  integer, parameter :: kind_nf = kind (1) 
# endif
  integer (kind_nf) :: netcdf_real=0, netcdf_int=0
  logical :: test = .false.

contains

  function get_netcdf_code_precision () result (code_real)
    use constants, only: pi, kind_rs, kind_rd
    use file_utils, only: error_unit
    integer :: code_real
# ifdef NETCDF

    ! second condition for Cray
    if ( (kind(pi)==kind_rs) .or. (kind_rs==kind_rd) ) then
       code_real = NF90_FLOAT
    else if (kind(pi)==kind_rd) then
       code_real = NF90_DOUBLE
    else
       write (error_unit(),*) &
            'ERROR: precision mismatch in get_netcdf_code_precision'
    end if
# endif
  end function get_netcdf_code_precision

  subroutine check_netcdf_file_precision (ncid, filename)

    use file_utils, only: error_unit

    integer (kind_nf), intent (in), optional :: ncid
    character (*), intent (in), optional :: filename
# ifdef NETCDF
    integer (kind_nf) :: file_real
    integer (kind_nf) :: ist, ncid_private, tid
    integer :: ierr

    !SET integer type: NOTE: This is not checked for compatability!!!! GGH 20 JAN 2012
    netcdf_int=NF90_INT

    ist = NF90_NOERR
    file_real = -1

    if (present(ncid)) then
       if (present(filename)) then
          ierr = error_unit()
          write (ierr,*) 'WARNING: in calling check_netcdf_file_precision'
          write (ierr,*) &
               'WARNING: both filename and ncid given -- filename ignored'
       end if
       ncid_private = ncid
    else
       if (present(filename)) then
          ist = nf90_open (filename, NF90_NOWRITE, ncid_private)
          if (test) write (error_unit(),*) &
               'opened netcdf file ', trim(filename), ' with ncid: ', &
               ncid_private, ' in check_netcdf_file_precision'
          if (ist /= NF90_NOERR) then
             call netcdf_error (ist, file=filename)
             return
          end if
       else
          ierr = error_unit()
          write (ierr,*) 'ERROR: in calling check_netcdf_file_precision'
          write (ierr,*) 'ERROR: either filename or ncid should be given'
          return
       end if
    end if

    ist = nf90_inq_varid (ncid_private, 't0', tid)
    if (ist /= NF90_NOERR) call netcdf_error (ist, var='t0')

    ! get file_real
    if (ist == NF90_NOERR) then
       ist = nf90_inquire_variable (ncid_private, tid, xtype=file_real)
       if (ist /= NF90_NOERR) call netcdf_error (ist, ncid_private, tid)
    end if

    if (.not.present(ncid)) then
       ist = nf90_close (ncid_private)
       if (ist /= NF90_NOERR) call netcdf_error (ist, file=filename)
    end if

    ! check if file_real == code_real
    if (file_real /= netcdf_real) then
       ierr = error_unit()
       write (ierr,*) 'WARNING: precision mismatch in input netcdf file and running code'
       if (file_real == NF90_FLOAT) then
          write (ierr,*) 'WARNING: file_real = NF90_FLOAT'
       else if (file_real == NF90_DOUBLE) then
          write (ierr,*) 'WARNING: file_real = NF90_DOUBLE'
       else
          write (ierr,*) 'WARNING: unknown file_real', file_real
       end if
       if (netcdf_real == NF90_FLOAT) then
          write (ierr,*) 'WARNING: code_real = NF90_FLOAT'
       else if (netcdf_real == NF90_DOUBLE) then
          write (ierr,*) 'WARNING: code_real = NF90_DOUBLE'
       else
          write (ierr,*) 'WARNING: unknown code_real'
       end if
    end if
# endif
  end subroutine check_netcdf_file_precision

  subroutine netcdf_error &
       (istatus, ncid, varid, dimid, file, dim, var, att, message, abort)

    use file_utils, only: error_unit
    use mp, only: proc0, finish_mp
# ifdef NETCDF
    use netcdf, only: NF90_GLOBAL
# endif

    integer (kind_nf), intent (in) :: istatus
    integer (kind_nf), intent (in), optional :: ncid
    integer (kind_nf), intent (in), optional :: varid
    integer (kind_nf), intent (in), optional :: dimid
    character (*), intent (in), optional :: file
    character (*), intent (in), optional :: dim
    character (*), intent (in), optional :: var
    character (*), intent (in), optional :: att
    character (*), intent (in), optional :: message
    logical, intent (in), optional :: abort
# ifdef NETCDF
    integer (kind_nf) :: ist
    integer :: ierr
    character (20) :: varname, dimname

    ierr = error_unit()

    write (ierr, '(2a)', advance='no') 'ERROR: ', trim (nf90_strerror (istatus))
    ! TT: If $ control fails, there is an alternative advance='no' specifier

    if (present(file)) &
         write (ierr, '(2a)', advance='no') ' in file: ', trim (file)

    if (present(dim)) &
         write (ierr, '(2a)', advance='no') ' in dimension: ', trim (dim)

    if (present(var)) &
         write (ierr, '(2a)', advance='no') ' in variable: ', trim (var)

    if (present(varid)) then
       if (present(ncid)) then
          if ( (varid == NF90_GLOBAL) .and. present(att) ) then
             write (ierr, '(2a)') ' in global attribute: ', trim(att)
             return
          else
             ist = nf90_inquire_variable (ncid, varid, varname)
             if (ist == NF90_NOERR) then
                write (ierr, '(a,i8,2a)', advance='no') ' in varid: ', varid, &
                     & ' variable name: ', trim (varname)
             else
                write (ierr, *)
                write (ierr, '(3a,i8,a,i8)', advance='no') 'ERROR in netcdf_error: ', &
                     trim (nf90_strerror(ist)), ' in varid: ', varid, &
                     ', ncid: ', ncid
             end if
          end if
          if (present(att)) &
               write (ierr, '(2a)') ' with the attribute: ', trim(att)
       else
          write (ierr, *)
          write (ierr, '(2a)', advance='no') 'ERROR in netcdf_error: ', &
               & 'ncid missing while varid present in the argument'
       end if
    end if

    if (present(dimid)) then
       if (present(ncid)) then
          ist = nf90_inquire_dimension (ncid, dimid, dimname)
          if (ist == NF90_NOERR) then
             write (ierr, '(a,i8,2a)', advance='no') ' in dimid: ', dimid, &
                  & ' dimension name: ', trim (dimname)
          else
             write (ierr, *)
             write (ierr, '(3a,i8,a,i8)', advance='no') 'ERROR in netcdf_error: ', &
                  trim (nf90_strerror(ist)), ' in dimid: ', dimid, &
                  ', ncid: ', ncid
          end if
       else
          write (ierr, *)
          write (ierr, '(2a)', advance='no') 'ERROR in netcdf_error: ', &
               & 'ncid missing while dimid present in the argument'
       end if
    end if

    if (present(message)) write (ierr, '(a)', advance='no') trim(message)
    
    ! append line-break
    write(ierr,*)

    ! if error is detected, the program should abort immediately
    if(present(abort)) then
       if(abort) then
           call finish_mp
           if(proc0) stop 'Aborted by netcdf_error'
          stop
       endif
    endif
# endif
  end subroutine netcdf_error

end module netcdf_utils
