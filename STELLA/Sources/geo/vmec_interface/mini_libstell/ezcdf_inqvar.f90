MODULE ezcdf_inqvar
 
  IMPLICIT NONE
!DEC$ IF DEFINED (NETCDF)
  PUBLIC :: cdfInqVar, cdfgv, cdfInqV, cdf_inquire, alpha_numeric

  INTERFACE cdfInqVarDim
      MODULE PROCEDURE cdfInqV, cdfgv, cdf_inquire
  END INTERFACE

  PRIVATE
  INCLUDE "netcdf.inc"

  CHARACTER*(nf_max_name) :: varnam_noalpha
  PRIVATE varnam_noalpha
 
CONTAINS
 
subroutine cdfInqVar(ncid,varnam,dimlens,eztype,ier)
  ! Inquire a Variable and its dimensions
  ! 03/08/99 C. Ludescher
  ! C. Ludescher/A. Pletzer Tue Apr  4 10:11:33 EDT 2000
  ! + support for complex type (ap) Wed May 16 15:18:05 EDT 2001
  implicit none
  ! Input
  integer,       intent(in)          :: ncid
  character*(*), intent(in)          :: varnam
  ! Output
  integer, dimension(:), intent(out) :: dimlens
  character*(*),         intent(out) :: eztype
  integer, optional,     intent(out) :: ier
  ! Local
  integer                 :: ndims, varid, natts, xtype
  integer                 :: status, i
  integer, dimension(3)   :: dimids
  character*(nf_max_name) :: name
  integer, parameter      :: cmplx_len = 13
  character(cmplx_len), parameter :: cmplx_name = '__CmPlx_Re_Im'
  character*11, parameter :: logical_name = '__logical__'
  logical :: is_complex, is_logical
 
  if (PRESENT (ier)) ier = 1
  is_complex = .false.
  is_logical = .false.
 
! replace any non-alphanumeric characters with underscores
  varnam_noalpha = varnam
  CALL alpha_numeric(varnam_noalpha)

  status = nf_inq_varid(ncid,varnam_noalpha,varid)
  if (status .ne. 0) then
     ! perhaps varnam is complex, try...
     status = nf_inq_varid(ncid,trim(varnam_noalpha)//cmplx_name,varid)
     if(status .eq. 0) then
        is_complex = .true.
     else
        status = nf_inq_varid(ncid,trim(varnam_noalpha)//logical_name,varid)
        if(status .eq. 0) is_logical = .true.
     endif
  endif

  if(is_complex)  varnam_noalpha = trim(varnam_noalpha)//cmplx_name
  if(is_logical)  varnam_noalpha = trim(varnam_noalpha)//logical_name
! call handle_err(status,varnam_noalpha,'cdfInqVar','nf_inq_varid')
  if(status .ne. 0) return
 
  status = nf_inq_var(ncid,varid,name,xtype,ndims,dimids,natts)
  call handle_err(status,varnam_noalpha,'cdfInqVar','nf_inq_var')
  if (status .ne. 0) return
 
  if (size(dimlens) .lt. ndims) return
  dimlens = 0
 
  select case (xtype)
  case (nf_double)
     eztype = 'R8'
     if(is_complex) eztype = 'C16'
  case (nf_int)
     eztype = 'INT'
     if(is_logical) eztype = 'LOG'
!  case (nf_byte)
!     eztype = 'LOG'
  case (nf_float)
     eztype = 'R4'
     if(is_complex) eztype = 'C8'
  case (nf_char)
     eztype = 'CHAR'
  end select
 
  do i=1,ndims
     status = nf_inq_dim(ncid,dimids(i),name,dimlens(i))
     call handle_err(status,varnam,'cdfInqVar','nf_inq_dim')
  end do
 
  if(is_complex) then
     dimlens(1) = dimlens(1)/2
  endif
  if (PRESENT (ier)) ier = status
 
end subroutine cdfInqVar
 
  ! automatic conversion to free f90 compatible form
  ! free.pl cdfgv.for
  ! linewidth: 72
  ! file names: cdfgv.for
  !
  SUBROUTINE cdfgv(ncid,varnam,varid,dimlens,sizes,xtype,status)
    !
    !     Get Variable id, etc.
    !     02/11/99 C.Ludescher
    ! C. Ludescher/A. Pletzer Tue Apr  4 10:11:33 EDT 2000
    !
    implicit none
!!$      include "netcdf.inc"
    ! Input
    integer                 :: ncid
    character*(*)           :: varnam
    character*1             :: xtype
    integer, dimension(:)   :: sizes
    ! Output
    integer, dimension(:)   :: dimlens
    integer                 :: varid, status
    ! Local
    integer                 :: i, vartyp, ndims, atts, rank
    integer, dimension(3)   :: dimids
    character*(nf_max_name) :: name

!   replace any non-alphanumeric characters with underscores
    varnam_noalpha = varnam
    CALL alpha_numeric(varnam_noalpha)
    
    status = nf_inq_varid(ncid,varnam_noalpha,varid)
!   call handle_err(status,varnam,'cdfgv','nf_inq_varid')
    if (status .ne. 0) return
    status = nf_inq_var(ncid,varid,name,vartyp,ndims,dimids,atts)
    call handle_err(status,varnam,'cdfgv','nf_inq_var')
    if (status .ne. 0) return
    ! Verify input dimension is correct
    rank = size(sizes)
    status = 1
    if (ndims .eq. 3 .and. rank .ne. 3) then
       print "('% cdfgv: --E-- The variable ',a,                      &
            &         ' is 3 dimensional')",varnam
       return
    else if (ndims .eq. 2 .and. rank .ne. 2) then
       print "('% cdfgv: --E-- The variable ',a,                      &
            &         ' is 2 dimensional')",varnam
       return
    endif
    if (ndims .eq. 0 .and. sizes(1) .ne. 0) then
       print "('% cdfgv: --E-- The variable ',a,                      &
            &           ' is a Scalar')",varnam
       return
    else if (ndims .eq. 1 .and. rank .ne. 1 ) then
       print "('% cdfgv: --E-- The variable ',a,                      &
            &           ' is 1 dimensional')",                                 &
            &         varnam
       return
    endif
    ! Verify data type is matching
    select case (xtype)
    case ('i')
       if(vartyp .ne. nf_int) then
          print "('% cdfgv: --E-- ',a,' is of type Integer !')",      &
               &             varnam
          return
       end if
    case ('l')
       if(vartyp .ne. nf_byte) then
          print "('% cdfgv: --E-- ',a,' is of type logical !')",         &
               &            varnam
          return
       end if
    case ('d')
       if(vartyp .ne. nf_double) then
          print "('% cdfgv: --E-- ',a,' is of type REAL*8 !')",       &
               &            varnam
          return
       end if
    case ('r')
       if(vartyp .ne. nf_real) then
          print "('% cdfgv: --E-- ',a,' is of type default REAL !')",       &
               &            varnam
          return
       end if
    case ('c')
       if(vartyp .ne. nf_char) then
          print "('% cdfgv: --E-- ',a,' is of type Character !')",         &
               &            varnam
          return
       end if
    end select
    status = 0
    do i=1,ndims
       dimlens(i) = 0
       status = nf_inq_dim(ncid,dimids(i),name,dimlens(i))
       call handle_err(status,varnam,'cdfgv','nf_inq_dim')
    end do
    ! Check array size is big enough
    select case (ndims)
    case (1)
       if(dimlens(1) .gt. sizes(1) ) then
          print "('% cdfgv: --W-- Output array size =',I6,/           &
               &           '                is smaller than ',                    &
               &           a,' size =',I6/,                                       &
               &           '                output will be truncated !')",        &
               &         sizes(1),varnam,dimlens(1)
       endif
    case (2)
       if(dimlens(1) .gt. sizes(1) .or. dimlens(2) .gt. sizes(2)) then
          print "('% cdfgv: --W-- Output array size =',I6,' *',I6,/      &
               &           '                is smaller than ',                    &
               &           a,' size =',I6,' *',I6/,                               &
               &           '                output will be truncated !')",        &
               &         sizes(1),sizes(2),varnam,dimlens(1),dimlens(2)
       endif
    case (3)
       if(dimlens(1) .gt. sizes(1) .or. dimlens(2) .gt. sizes(2)         &
            &   .or. dimlens(3) .gt. sizes(3)) then
          print "('% cdfgv: --W-- Output array size =',                  &
               &          I5,' *',I5,' *',I5/,                                    &
               &           '                is smaller than ',                    &
               &           a,' size =',I5,' *',I5,' *',I5/,                       &
               &           '                output will be truncated !')",        &
               &         sizes(1),sizes(2),sizes(3),varnam,                       &
               &         dimlens(1),dimlens(2),dimlens(3)
       endif
    end select

  end SUBROUTINE cdfgv
 
  SUBROUTINE cdfInqV(ncid,varnam,varid,dimlens,ndims,status)
    ! Inquire variable-id and dimlens
    ! 03/09/99 C.Ludescher
    !
    implicit none
    !
!!$      include "netcdf.inc"
    ! Input
    integer,       intent(in)  :: ncid
    character*(*), intent(in)  :: varnam
    ! Returns
    !      integer, dimension(3)  ::  dimlens
    !      integer ::  ndims, varid
    !      integer ::  status
    integer, dimension(:), intent(out)   ::  dimlens
    integer,               intent(out)   ::  ndims, varid
    integer,               intent(out)   ::  status
    ! Local
    integer                 :: natts, xtype, i
    integer, dimension(3)   :: dimids
    character*(nf_max_name) :: name
    !---------------------------------------------------------------------------
!   Initialize values SAL 07012014
    varid = 0 ! SAL 07012014
    xtype=0; ndims=0; dimids=0; natts=0; dimlens=0 ! SAL 07012014

!   replace any non-alphanumeric characters with underscores
    varnam_noalpha = varnam
    CALL alpha_numeric(varnam_noalpha)

    status = nf_inq_varid(ncid,varnam_noalpha,varid)
!    call handle_err(status,varnam,'cdfInqV','nf_inq_varid')
    if (status .ne. 0) return
    status = nf_inq_var(ncid,varid,name,xtype,ndims,dimids,natts)
    call handle_err(status,varnam,'cdfInqV','nf_inq_var')

    if (ndims .gt. size(dimlens)) stop 'dimlens too small in cdfInqV'

    do i=1,ndims
       status = nf_inq_dimlen(ncid,dimids(i),dimlens(i))
       call handle_err(status,varnam,'cdfInqV','nf_inq_dimlen')
    end do
  END SUBROUTINE cdfInqV
 
  SUBROUTINE cdf_inquire(ncid,varnam,dimlens,xtype,ier)
       implicit none
       ! Input
       integer,       intent(in)            :: ncid
       character*(*), intent(in)            :: varnam
       ! Output
       integer, dimension(:),   intent(out) :: dimlens
       character*(*), optional, intent(out) :: xtype
       integer, optional,       intent(out) :: ier
       integer :: ezerror
       character*4 :: eztype
 
       CALL cdfInqVar(ncid, varnam, dimlens, eztype, ezerror)
 
       IF (PRESENT(xtype)) xtype = eztype
       IF (PRESENT(ier)) ier = ezerror
 
  END SUBROUTINE cdf_inquire
 
  SUBROUTINE alpha_numeric(string)
     IMPLICIT NONE
     CHARACTER*(*), INTENT(INOUT)     :: string
     CHARACTER*(LEN_TRIM(string))     :: temp
     INTEGER                          :: iascii, i

  ! 04/03/03 S. Hirshman
  ! replaces any non-alphanumeric characters with underscores
     temp = adjustl(string)
     string = trim(temp)
     do i = 2, len_trim(string)
       iascii = iachar(string(i:i))
       if ( ((iascii >= iachar('0')) .and. (iascii <= iachar('9')))           &
      .or.  ((iascii >= iachar('A')) .and. (iascii <= iachar('z')))) cycle
       string(i:i) = '_'
     end do

  END SUBROUTINE alpha_numeric
!DEC$ ENDIF
END MODULE ezcdf_inqvar
