subroutine handle_err(status,nam3,nam1,nam2)
!DEC$ IF DEFINED (NETCDF)
  include "netcdf.inc"
  INTEGER, intent(in) :: status
  character*(*), intent(in) :: nam1, nam2, nam3
  if (status .ne. nf_noerr) then
     WRITE(*,10) nam1,nam2,nam3
10   format('% ',a,'--E-- A netCDF error has occurred in: ' ,a/      &
          &          'while processing: ',a)
     print *, nf_strerror(status)
  endif
!DEC$ ELSE
return
!DEC$ ENDIF
end subroutine handle_err
