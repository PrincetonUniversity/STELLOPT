MODULE ezcdf_attrib
   USE ezcdf_inqvar
#ifdef NETCDF
   include "netcdf.inc"
   INTEGER, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
   INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
   CHARACTER(len=nf_max_name) :: varnam_noalpha
   PRIVATE :: r4, r8, varnam_noalpha
   PUBLIC :: cdf_setatt, cdf_getatt
 
   INTERFACE cdf_setatt
      MODULE PROCEDURE cdfSetatt, cdfsa_i, cdfsa_d, cdfsa_f
   END INTERFACE
 
   INTERFACE cdf_getatt
      MODULE PROCEDURE cdfGetatt, cdfga_i, cdfga_d, cdfga_f
   END INTERFACE

CONTAINS
 
SUBROUTINE cdfSetatt(ncid,varnam,long_name,units,ier,varid)
  implicit none
!  include "netcdf.inc"
  ! Input
  integer,       intent(in)           :: ncid
  character*(*), intent(in)           :: varnam
  character*(*), intent(in), optional :: long_name, units
  ! Output
  integer, optional,     intent(out) :: ier, varid
  ! Local
  integer                 :: varid0
  integer                 :: status
  character*(*), parameter :: cmplx_name = '__CmPlx_Re_Im'
  character*(*), parameter :: logical_name = '__logical__'
 
  if (PRESENT (ier)) ier = 1
 
  varnam_noalpha = varnam
  CALL alpha_numeric(varnam_noalpha)

!Find varid
  status = nf_inq_varid(ncid,varnam_noalpha,varid0)
  if (status .ne. 0) then
     ! perhaps varnam is complex, try...
     status = nf_inq_varid(ncid,trim(varnam_noalpha)//cmplx_name,varid0)
     if(status .ne. 0) then
        status = nf_inq_varid(ncid,trim(varnam_noalpha)//logical_name,varid0)
     endif
  endif
 
!  call handle_err(status,varnam,'cdf_setatt','nf_inq_varid')
  if(status .ne. 0) GO TO 100
 
  if (PRESENT(long_name)) then
     status = nf_put_att_text (ncid, varid0, 'long_name', LEN_TRIM(long_name),  &
                          TRIM(long_name))
     CALL handle_err(status,long_name,'cdf_setatt','nf_put_att_text')
     if (status .ne. 0) GO TO 100
  end if
 
  if (PRESENT(units)) then
     status = nf_put_att_text (ncid, varid0, 'units', LEN_TRIM(units),          &
                          TRIM(units))
     CALL handle_err(status,units,'cdf_setatt','nf_put_att_text')
     if (status .ne. 0) GO TO 100
  end if
 
100 CONTINUE
 
  if (PRESENT (ier)) ier = status
  if (PRESENT (varid)) varid = varid0
 
END SUBROUTINE cdfSetatt
 
SUBROUTINE cdf_SetTitle(ncid,title,ier)
  implicit none
!  include "netcdf.inc"
  ! Input
  integer,       intent(in)           :: ncid
  character*(*), intent(in)           :: title
  ! Output
  integer, optional,     intent(out) :: ier
  ! Local
  integer                 :: status
 
  if (PRESENT (ier)) ier = 1
 
  status = nf_put_att_text (ncid, nf_global, 'title', LEN_TRIM(title),  &
                          TRIM(title))
  CALL handle_err(status,title,'cdf_settitle','nf_put_att_text')
 
  if (PRESENT (ier)) ier = status
 
END SUBROUTINE cdf_SetTitle
 
SUBROUTINE cdfsa_i(ncid,varnam,valid_range,long_name,units,ier)
  implicit none
!  include "netcdf.inc"
  ! Input
  integer,       intent(in)           :: ncid
  character*(*), intent(in)           :: varnam
  character*(*), intent(in), optional :: long_name, units
  integer, intent(in), dimension(2)   :: valid_range
  ! Output
  integer, optional,     intent(out) :: ier
  ! Local
  integer                 :: varid
  integer                 :: status
 
  call cdfSetatt(ncid,varnam,long_name,units,ier,varid)
 
  status = nf_put_att_int(ncid, varid, 'valid_range', nf_double, &
                          2, valid_range)
  CALL handle_err(status,'valid_range','cdf_setatt','nf_put_att_double')
 
END SUBROUTINE cdfsa_i
 
 
SUBROUTINE cdfsa_f(ncid,varnam,valid_range,long_name,units,ier)
  implicit none
!  include "netcdf.inc"
  ! Input
  integer,       intent(in)               :: ncid
  character*(*), intent(in)               :: varnam
  character*(*), intent(in), optional     :: long_name, units
  real(kind=r4), intent(in), dimension(2) :: valid_range
  ! Output
  integer, optional,     intent(out) :: ier
  ! Local
  integer                 :: varid
  integer                 :: status
 
  call cdfSetatt(ncid,varnam,long_name,units,ier,varid)
 
  status = nf_put_att_real(ncid, varid, 'valid_range', nf_double, &
                           2, valid_range)
  CALL handle_err(status,'valid_range','cdf_setatt','nf_put_att_double')
 
END SUBROUTINE cdfsa_f
 
SUBROUTINE cdfsa_d(ncid,varnam,valid_range,long_name,units,ier)
  implicit none
!  include "netcdf.inc"
  ! Input
  integer,       intent(in)               :: ncid
  character*(*), intent(in)               :: varnam
  character*(*), intent(in), optional     :: long_name, units
  real(kind=r8), intent(in), dimension(2) :: valid_range
  ! Output
  integer, optional,     intent(out) :: ier
  ! Local
  integer                 :: varid
  integer                 :: status
 
  call cdfSetatt(ncid,varnam,long_name,units,ier,varid)
 
  status = nf_put_att_double(ncid, varid, 'valid_range', nf_double, &
                             2, valid_range)
  CALL handle_err(status,'valid_range','cdf_setatt','nf_put_att_double')
 
END SUBROUTINE cdfsa_d
 
 
SUBROUTINE cdfGetatt(ncid,varnam,long_name,units,ier,varid)
  implicit none
!  include "netcdf.inc"
  ! Input
  integer,       intent(in)           :: ncid
  character*(*), intent(in)           :: varnam
  ! Output
  character*(*), intent(out), optional :: long_name, units
  integer, optional,     intent(out) :: ier, varid
  ! Local
  integer                 :: varid0
  integer                 :: status
  character*(*), parameter :: cmplx_name = '__CmPlx_Re_Im'
  character*(*), parameter :: logical_name = '__logical__'
  character*(nf_max_name) :: name
 
  if (PRESENT (ier)) ier = 1

  varnam_noalpha = varnam
  CALL alpha_numeric(varnam_noalpha)
 
!Find varid0
  status = nf_inq_varid(ncid,varnam_noalpha,varid0)
  if (status .ne. 0) then
     ! perhaps varnam is complex, try...
     status = nf_inq_varid(ncid,trim(varnam_noalpha)//cmplx_name,varid0)
     if(status .ne. 0) then
        status = nf_inq_varid(ncid,trim(varnam_noalpha)//logical_name,varid0)
     endif
  endif
 
!  call handle_err(status,varnam,'cdf_setatt','nf_inq_varid')
  if(status .ne. 0) GO TO 100
 
  if (PRESENT(long_name)) then
     name = ""
     status = nf_get_att_text (ncid, varid0, 'long_name', name)
     if (status .eq. nf_noerr) then
        long_name = name(1:len(name))
     else
        long_name = ""
     end if
  end if
 
  if (PRESENT(units)) then
     name = ""
     status = nf_get_att_text (ncid, varid0, 'units', name)
     if (status .eq. nf_noerr) then
        units = name(1:len(name))
     else
        units = ""
     end if
  end if
 
100 CONTINUE
 
  if (PRESENT (ier)) ier = status
  if (PRESENT (varid)) varid = varid0
 
END SUBROUTINE cdfGetatt
 
SUBROUTINE cdf_GetTitle(ncid,title,ier)
  implicit none
!  include "netcdf.inc"
  ! Input
  integer,       intent(in)           :: ncid
  ! Output
  character*(*), intent(out)         :: title
  integer, optional,     intent(out) :: ier
  ! Local
  integer                 :: status
  character*(nf_max_name) :: name
 
  if (PRESENT (ier)) ier = 1
 
  status = nf_get_att_text (ncid, nf_global, 'title', name)
  if (status .eq. nf_noerr) then
     title = name(1:len(name))
  else
     title = ""
  end if
 
  if (PRESENT (ier)) ier = status
 
END SUBROUTINE cdf_GetTitle
 
SUBROUTINE cdfga_i(ncid,varnam,valid_range,long_name,units,ier)
  implicit none
!  include "netcdf.inc"
  ! Input
  integer,       intent(in)            :: ncid
  character*(*), intent(in)            :: varnam
  ! Output
  character*(*), intent(out), optional :: long_name, units
  integer, intent(out), dimension(2)   :: valid_range
  integer, optional,     intent(out)   :: ier
  ! Local
  integer                 :: varid
  integer                 :: status
 
  call cdfGetatt(ncid,varnam,long_name,units,ier,varid)
 
  status = nf_get_att_int(ncid, varid, 'valid_range', valid_range)
  CALL handle_err(status,'valid_range','cdf_getatt','nf_get_att_int')
 
END SUBROUTINE cdfga_i
 
SUBROUTINE cdfga_f(ncid,varnam,valid_range,long_name,units,ier)
  implicit none
!  include "netcdf.inc"
  ! Input
  integer,       intent(in)           :: ncid
  character*(*), intent(in)           :: varnam
  ! Output
  character*(*), intent(out), optional     :: long_name, units
  real(kind=r4), intent(out), dimension(2) :: valid_range
  integer, optional,     intent(out) :: ier
  ! Local
  integer                 :: varid
  integer                 :: status
 
  call cdfGetatt(ncid,varnam,long_name,units,ier,varid)
 
  status = nf_get_att_real(ncid, varid, 'valid_range', valid_range)
  CALL handle_err(status,'valid_range','cdf_getatt','nf_get_att_real')
 
END SUBROUTINE cdfga_f
 
SUBROUTINE cdfga_d(ncid,varnam,valid_range,long_name,units,ier)
  implicit none
!  include "netcdf.inc"
  ! Input
  integer,       intent(in)           :: ncid
  character*(*), intent(in)           :: varnam
  ! Output
  character*(*), intent(out), optional     :: long_name, units
  real(kind=r8), intent(out), dimension(2) :: valid_range
  integer, optional,     intent(out) :: ier
  ! Local
  integer                 :: varid
  integer                 :: status
 
  call cdfGetatt(ncid,varnam,long_name,units,ier,varid)
 
  status = nf_get_att_double(ncid, varid, 'valid_range', valid_range)
  CALL handle_err(status,'valid_range','cdf_getatt','nf_get_att_double')
 
END SUBROUTINE cdfga_d
#endif
END MODULE ezcdf_attrib
 
