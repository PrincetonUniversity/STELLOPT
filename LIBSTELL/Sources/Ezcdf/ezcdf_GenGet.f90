MODULE ezcdf_GenGet
  USE ezcdf_opncls
  USE ezcdf_inqvar
#ifdef NETCDF
  EXTERNAL handle_err
 
  PRIVATE
 
  include "netcdf.inc"

  INTEGER, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
  INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
  CHARACTER*(*), PARAMETER :: cmplx_name = '__CmPlx_Re_Im'
  PRIVATE :: r4, r8, cmplx_name
 
  PUBLIC :: cdfr_3i, cdfr_3l, cdfr_3d, cdfr_3c16, cdfr_3f, cdfr_3c8,           &
       & cdfr_2i, cdfr_2l, cdfr_2d, cdfr_2c16, cdfr_2f, cdfr_2c8, cdfr_2c,     &
       & cdfr_1i, cdfr_1l, cdfr_1d, cdfr_1c16, cdfr_1f, cdfr_1c8, cdfr_1c,     &
       & cdfr_0i, cdfr_0l, cdfr_0d, cdfr_0c16, cdfr_0f, cdfr_0c8, cdfGetVar,   &
       & cdf_read
 

  ! Generic Interface to Read netcdf data Variables
  ! 03/10/99 C. Ludescher
  ! C. Ludescher/A. Pletzer Tue Apr  4 10:11:33 EDT 2000
  ! + support for complex types (ap) Wed May 16 15:18:05 EDT 2001
  !====================================================================
  ! Generic Read Routines: cdfGetVar

  ! alias to cdfGetVar (ifc does not like => )
  INTERFACE cdf_read
     MODULE PROCEDURE cdfr_3i, cdfr_3l, cdfr_3d, cdfr_3c16, cdfr_3f, cdfr_3c8,           &
                      cdfr_2i, cdfr_2l, cdfr_2d, cdfr_2c16, cdfr_2f, cdfr_2c8, cdfr_2c,  &
                      cdfr_1i, cdfr_1l, cdfr_1d, cdfr_1c16, cdfr_1f, cdfr_1c8, cdfr_1c,  &
                      cdfr_0i, cdfr_0l, cdfr_0d, cdfr_0c16, cdfr_0f, cdfr_0c8
  END INTERFACE
  ! same as above (Intel compiler does not handle well aliases)
  INTERFACE cdfGetVar
     MODULE PROCEDURE cdfr_3i, cdfr_3l, cdfr_3d, cdfr_3c16, cdfr_3f, cdfr_3c8,           &
                      cdfr_2i, cdfr_2l, cdfr_2d, cdfr_2c16, cdfr_2f, cdfr_2c8, cdfr_2c,  &
                      cdfr_1i, cdfr_1l, cdfr_1d, cdfr_1c16, cdfr_1f, cdfr_1c8, cdfr_1c,  &
                      cdfr_0i, cdfr_0l, cdfr_0d, cdfr_0c16, cdfr_0f, cdfr_0c8
  END INTERFACE

CONTAINS
!---------------------------------------------
!cdfGetVar implementation routines 

SUBROUTINE cdfr_3i(ncid,varnam,varval,ier)
  ! Read 3 dimensional Integer array
  !
  implicit none
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  integer, dimension(:,:,:), intent(inout) :: varval
  integer, optional,         intent(out) :: ier
  ! Local
  integer, dimension(3)   :: st,     cnt,    ldim
  integer                 :: varid,  status, j,  k
  integer, dimension(3)   :: dimlens

  integer, dimension(:,:,:), allocatable :: temp
  integer :: ndim1, ndim2, ndim3 

  if (PRESENT (ier)) ier = 1
  ldim = shape(varval)
  call cdfgv(ncid,varnam,varid,dimlens,ldim,'i',status)
  if (status .ne. 0) return
  !
  if(dimlens(1)==ldim(1).and.dimlens(2)==ldim(2).and.dimlens(3)==ldim(3))then
     status = nf_get_var_int(ncid, varid, varval)
     if(status/=NF_NOERR) then
        call handle_err(status,varnam,'cdfr_3i', 'nf_get_var_int')
        return
     endif
  else
     allocate(temp(dimlens(1), dimlens(2), dimlens(3)))
     status = nf_get_var_int(ncid, varid, temp)
     if(status/=NF_NOERR) then
        call handle_err(status,varnam,'cdfr_3i', 'nf_get_var_int')
        deallocate(temp)
        return
     endif
     ndim1 = min(dimlens(1), ldim(1))
     ndim2 = min(dimlens(2), ldim(2))
     ndim3 = min(dimlens(3), ldim(3))
     varval(1:ndim1, 1:ndim2, 1:ndim3) = temp(1:ndim1, 1:ndim2, 1:ndim3)
     deallocate(temp)
  endif

!!$  st(1) = 1
!!$  cnt(1) = min(dimlens(1),ldim(1)) ! x count
!!$  cnt(2) = 1
!!$  cnt(3) = 1
!!$  do k = 1,min(dimlens(2),ldim(2))
!!$     st(2) = k
!!$     do j = 1,min(dimlens(3),ldim(3)) ! For each Z : read slab in varval
!!$        st(3) = j           ! Start of slab
!!$        status = nf_get_vara_int(ncid,varid,st,cnt,varval(1,k,j))
!!$        if (status .ne. NF_NOERR) then
!!$           call handle_err(status,varnam,'cdfr_3i',                 &
!!$                &              'nf_get_vara_int')
!!$           return
!!$        end if
!!$     end do
!!$  end do
  if (PRESENT (ier)) ier = 0
END SUBROUTINE cdfr_3i
 
SUBROUTINE cdfr_3l(ncid,varnam,varval,ier)
  ! Read 3 dimensional logical array
  !
  implicit none
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  logical, dimension(:,:,:), intent(inout) :: varval
  integer, optional,         intent(out) :: ier
  ! Local
  integer                                :: status
  integer, dimension(:,:,:), allocatable :: varval_i
  character*11, parameter :: logical_name = '__logical__'
 
  ALLOCATE (varval_i(size(varval,1), size(varval,2), size(varval,3)), stat=status)
  if (status .ne. 0) STOP 'Allocation error in cdf_getvar'
 
  call cdfr_3i(ncid,trim(varnam)//logical_name,varval_i,ier)
 
  WHERE (varval_i == 0)
     varval = .false.
  ELSEWHERE
     varval = .true.
  END WHERE
  DEALLOCATE (varval_i)
END SUBROUTINE cdfr_3l
 
SUBROUTINE cdfr_3d(ncid,varnam,varval,ier)
  implicit none
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  REAL(KIND=r8), dimension(:,:,:), intent(inout) ::  varval
  integer, optional,               intent(out) :: ier
  ! Local
  integer, dimension(3)   :: st,     cnt,    ldim
  integer                 :: varid,  status, j,  k
  integer, dimension(3)   :: dimlens

  real(r8), dimension(:,:,:), allocatable :: temp
  integer ndim1, ndim2, ndim3

  if (PRESENT (ier)) ier = 1
  ldim = shape(varval)
  call cdfgv(ncid,varnam,varid,dimlens,ldim,'d',status)
  if (status .ne. 0) return
  !

  if(dimlens(1)==ldim(1).and.dimlens(2)==ldim(2).and.dimlens(3)==ldim(3))then
     status = nf_get_var_double(ncid, varid, varval)
     if(status/=NF_NOERR) then
        call handle_err(status,varnam,'cdfr_3d', 'nf_get_var_double')
        return
     endif
  else
     allocate(temp(dimlens(1), dimlens(2), dimlens(3)))
     status = nf_get_var_double(ncid, varid, temp)
     if(status/=NF_NOERR) then
        call handle_err(status,varnam,'cdfr_3d', 'nf_get_var_double')
        deallocate(temp)
        return
     endif
     ndim1 = min(dimlens(1), ldim(1))
     ndim2 = min(dimlens(2), ldim(2))
     ndim3 = min(dimlens(3), ldim(3))
     varval(1:ndim1, 1:ndim2, 1:ndim3) = temp(1:ndim1, 1:ndim2, 1:ndim3)
     deallocate(temp)
  endif

!!$  st(1) = 1
!!$  cnt(1) = min(dimlens(1),ldim(1)) ! x count
!!$  cnt(2) = 1
!!$  cnt(3) = 1
!!$  do k = 1,min(dimlens(2),ldim(2))
!!$     st(2) = k
!!$     do j = 1,min(dimlens(3),ldim(3)) ! For each Z : read slab into varval
!!$        st(3) = j                      ! Start of slab
!!$        status = nf_get_vara_double(ncid,varid,st,cnt,varval(1,k,j))
!!$        if (status .ne. NF_NOERR) then
!!$           call handle_err(status,varnam,'cdfr_3d',                    &
!!$                &           'nf_get_vara_double')
!!$           return
!!$        end if
!!$     end do
!!$  end do

  if (PRESENT (ier)) ier = 0
END SUBROUTINE cdfr_3d

SUBROUTINE cdfr_3c16(ncid,varnam,varval,ier)
  implicit none
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  COMPLEX(KIND=r8), dimension(:,:,:), intent(inout) ::  varval
  integer, optional,               intent(out) :: ier
  ! Local
  integer, dimension(3)   :: st,     cnt,    ldim
  integer                 :: varid,  status, j,  k, i
  integer, dimension(3)   :: dimlens

  real(r8), dimension(:,:,:), allocatable :: temp
  integer ndim1, ndim2, ndim3
 
  if (PRESENT (ier)) ier = 1
  ldim = shape(varval)
  ldim(1) = 2*ldim(1) ! Re/Im pairs
  call cdfgv(ncid,trim(varnam)//cmplx_name,varid,dimlens,ldim,'d',status)
  if (status .ne. 0) return
  !
  if(dimlens(1)==ldim(1).and.dimlens(2)==ldim(2).and.dimlens(3)==ldim(3))then
     status = nf_get_var_double(ncid, varid, varval)
     if(status/=NF_NOERR) then
        call handle_err(status,varnam,'cdfr_3c16', 'nf_get_var_double')
        return
     endif
  else
     allocate(temp(dimlens(1), dimlens(2), dimlens(3)))
     status = nf_get_var_double(ncid, varid, temp)
     if(status/=NF_NOERR) then
        call handle_err(status,varnam,'cdfr_3c16', 'nf_get_var_double')
        deallocate(temp)
        return
     endif
     ndim1 = min(dimlens(1), ldim(1))
     ndim2 = min(dimlens(2), ldim(2))
     ndim3 = min(dimlens(3), ldim(3))
     do i = 1, ndim1/2
        varval(i, 1:ndim2, 1:ndim3) = temp(2*(i-1)+1, 1:ndim2, 1:ndim3) +  &
             & (0._r8,1._r8)*temp(2*(i-1)+2, 1:ndim2, 1:ndim3)
     enddo 
     deallocate(temp)
  endif
!!$
!!$  st(1) = 1
!!$  cnt(1) = min(dimlens(1),ldim(1)) ! x count
!!$  cnt(2) = 1
!!$  cnt(3) = 1
!!$  do k = 1,min(dimlens(2),ldim(2))
!!$     st(2) = k
!!$     do j = 1,min(dimlens(3),ldim(3)) ! For each Z : read slab into varval
!!$        st(3) = j                      ! Start of slab
!!$        status = nf_get_vara_double(ncid,varid,st,cnt,varval(1,k,j))
!!$        if (status .ne. NF_NOERR) then
!!$           call handle_err(status,varnam,'cdfr_3c16',                    &
!!$                &           'nf_get_vara_double')
!!$           return
!!$        end if
!!$     end do
!!$  end do
  if (PRESENT (ier)) ier = 0
END SUBROUTINE cdfr_3c16

SUBROUTINE cdfr_3f(ncid,varnam,varval,ier)
  implicit none
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  REAL(KIND=r4), dimension(:,:,:), intent(inout) ::  varval
  integer, optional,               intent(out) :: ier
  ! Local
  integer, dimension(3)   :: st,     cnt,    ldim
  integer                 :: varid,  status, j,  k
  integer, dimension(3)   :: dimlens

  real, dimension(:,:,:), allocatable :: temp
  integer ndim1, ndim2, ndim3

  if (PRESENT (ier)) ier = 1
  ldim = shape(varval)
  call cdfgv(ncid,varnam,varid,dimlens,ldim,'r',status)
  if (status .ne. 0) return
  !
  if(dimlens(1)==ldim(1).and.dimlens(2)==ldim(2).and.dimlens(3)==ldim(3))then
     status = nf_get_var_real(ncid, varid, varval)
     if(status/=NF_NOERR) then
        call handle_err(status,varnam,'cdfr_3f', 'nf_get_var_real')
        return
     endif
  else
     allocate(temp(dimlens(1), dimlens(2), dimlens(3)))
     status = nf_get_var_real(ncid, varid, temp)
     if(status/=NF_NOERR) then
        call handle_err(status,varnam,'cdfr_3f', 'nf_get_var_real')
        deallocate(temp)
        return
     endif
     ndim1 = min(dimlens(1), ldim(1))
     ndim2 = min(dimlens(2), ldim(2))
     ndim3 = min(dimlens(3), ldim(3))
     varval(1:ndim1, 1:ndim2, 1:ndim3) = temp(1:ndim1, 1:ndim2, 1:ndim3)
     deallocate(temp)
  endif

!!$  st(1) = 1
!!$  cnt(1) = min(dimlens(1),ldim(1)) ! x count
!!$  cnt(2) = 1
!!$  cnt(3) = 1
!!$  do k = 1,min(dimlens(2),ldim(2))
!!$     st(2) = k
!!$     do j = 1,min(dimlens(3),ldim(3)) ! For each Z : read slab into varval
!!$        st(3) = j                      ! Start of slab
!!$        status = nf_get_vara_real(ncid,varid,st,cnt,varval(1,k,j))
!!$        if (status .ne. NF_NOERR) then
!!$           call handle_err(status,varnam,'cdfr_3f',                    &
!!$                &           'nf_get_vara_real')
!!$           return
!!$        end if
!!$     end do
!!$  end do
  if (PRESENT (ier)) ier = 0
END SUBROUTINE cdfr_3f

SUBROUTINE cdfr_3c8(ncid,varnam,varval,ier)
  implicit none
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  COMPLEX(KIND=r4), dimension(:,:,:), intent(inout) ::  varval
  integer, optional,               intent(out) :: ier
  ! Local
  integer, dimension(3)   :: st,     cnt,    ldim
  integer                 :: varid,  status, j,  k
  integer, dimension(3)   :: dimlens

  real, dimension(:,:,:), allocatable :: temp
  integer ndim1, ndim2, ndim3, i
 
  if (PRESENT (ier)) ier = 1
  ldim = shape(varval)
  ldim(1) = 2*ldim(1) ! Re/Im pairs
  call cdfgv(ncid,trim(varnam)//cmplx_name,varid,dimlens,ldim,'r',status)
  if (status .ne. 0) return
  !
  if(dimlens(1)==ldim(1).and.dimlens(2)==ldim(2).and.dimlens(3)==ldim(3))then
     status = nf_get_var_real(ncid, varid, varval)
     if(status/=NF_NOERR) then
        call handle_err(status,varnam,'cdfr_3c8', 'nf_get_var_real')
        return
     endif
  else
     allocate(temp(dimlens(1), dimlens(2), dimlens(3)))
     status = nf_get_var_real(ncid, varid, temp)
     if(status/=NF_NOERR) then
        call handle_err(status,varnam,'cdfr_3c8', 'nf_get_var_real')
        deallocate(temp)
        return
     endif
     ndim1 = min(dimlens(1), ldim(1))
     ndim2 = min(dimlens(2), ldim(2))
     ndim3 = min(dimlens(3), ldim(3))
     do i = 1, ndim1/2
        varval(i, 1:ndim2, 1:ndim3) = temp(2*(i-1)+1, 1:ndim2, 1:ndim3) +  &
             & (0.,1.)*temp(2*(i-1)+2, 1:ndim2, 1:ndim3)
     enddo 
     deallocate(temp)
  endif

!!$  st(1) = 1
!!$  cnt(1) = min(dimlens(1),ldim(1)) ! x count
!!$  cnt(2) = 1
!!$  cnt(3) = 1
!!$  do k = 1,min(dimlens(2),ldim(2))
!!$     st(2) = k
!!$     do j = 1,min(dimlens(3),ldim(3)) ! For each Z : read slab into varval
!!$        st(3) = j                      ! Start of slab
!!$        status = nf_get_vara_real(ncid,varid,st,cnt,varval(1,k,j))
!!$        if (status .ne. NF_NOERR) then
!!$           call handle_err(status,varnam,'cdfr_3c8',                    &
!!$                &           'nf_get_vara_real')
!!$           return
!!$        end if
!!$     end do
!!$  end do
  if (PRESENT (ier)) ier = 0
END SUBROUTINE cdfr_3c8
 
SUBROUTINE cdfr_2i(ncid,varnam,varval,ier)
  ! Read 2 dimensional Integer array
  !
  implicit none
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  integer, dimension(:,:), intent(inout) :: varval
  integer, optional,       intent(out) :: ier
  ! Local
  integer, dimension(2)   :: st,     cnt,    ldim
  integer                 :: varid,  status, j
  integer, dimension(2)   :: dimlens

  integer, dimension(:,:), allocatable :: temp
  integer ndim1, ndim2

  if (PRESENT (ier)) ier = 1
  ldim = shape(varval)
  call cdfgv(ncid,varnam,varid,dimlens,ldim,'i',status)
  if (status .ne. 0) return
  !
  if(dimlens(1)==ldim(1).and.dimlens(2)==ldim(2))then
     status = nf_get_var_int(ncid, varid, varval)
     if(status/=NF_NOERR) then
        call handle_err(status,varnam,'cdfr_2i', 'nf_get_var_int')
        return
     endif
  else
     allocate(temp(dimlens(1), dimlens(2)))
     status = nf_get_var_int(ncid, varid, temp)
     if(status/=NF_NOERR) then
        call handle_err(status,varnam,'cdfr_2i', 'nf_get_var_int')
        deallocate(temp)
        return
     endif
     ndim1 = min(dimlens(1), ldim(1))
     ndim2 = min(dimlens(2), ldim(2))
     varval(1:ndim1, 1:ndim2) = temp(1:ndim1, 1:ndim2)
     deallocate(temp)
  endif

!!$  st(1) = 1
!!$  cnt(1) = min(dimlens(1),ldim(1)) ! x count
!!$  cnt(2) = 1
!!$  do j = 1,min(dimlens(2),ldim(2)) ! For each Y : read slab into varval
!!$     st(2) = j              ! Start of slab
!!$     status = nf_get_vara_int(ncid,varid,st,cnt,varval(1,j))
!!$     if (status .ne. NF_NOERR) then
!!$        call handle_err(status,varnam,'cdfr_2i',                    &
!!$             &           'nf_get_vara_int')
!!$        return
!!$     end if
!!$  end do
  if (PRESENT (ier)) ier = 0
END SUBROUTINE cdfr_2i
 
SUBROUTINE cdfr_2l(ncid,varnam,varval,ier)
  ! Read 2 dimensional logical array
  !
  implicit none
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  logical, dimension(:,:), intent(inout) :: varval
  integer, optional,       intent(out) :: ier
  ! Local
  integer                              :: status
  integer, dimension(:,:), allocatable :: varval_i
  character*11, parameter :: logical_name = '__logical__'
 
  ALLOCATE (varval_i(size(varval,1), size(varval,2)), stat=status)
  if (status .ne. 0) STOP 'Allocation error in cdf_getvar'
 
  call cdfr_2i(ncid,trim(varnam)//logical_name,varval_i,ier)
 
  WHERE (varval_i == 0)
     varval = .false.
  ELSEWHERE
     varval = .true.
  END WHERE
  DEALLOCATE (varval_i)
END SUBROUTINE cdfr_2l
 
SUBROUTINE cdfr_2d(ncid,varnam,varval,ier)
  implicit none
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  REAL(KIND=r8), dimension(:,:), intent(inout) ::  varval
  integer, optional,             intent(out) :: ier
  ! Local
  integer, dimension(2)   :: st,     cnt,    ldim
  integer                 :: varid,  status, j
  integer, dimension(2)   :: dimlens

  real(r8), dimension(:,:), allocatable :: temp
  integer ndim1, ndim2

  if (PRESENT (ier)) ier = 1
  ldim = shape(varval)
  call cdfgv(ncid,varnam,varid,dimlens,ldim,'d',status)
  if (status .ne. 0) return
  !
  if(dimlens(1)==ldim(1).and.dimlens(2)==ldim(2))then
     status = nf_get_var_double(ncid, varid, varval)
     if(status/=NF_NOERR) then
        call handle_err(status,varnam,'cdfr_2d', 'nf_get_var_double')
        return
     endif
  else
     allocate(temp(dimlens(1), dimlens(2)))
     status = nf_get_var_double(ncid, varid, temp)
     if(status/=NF_NOERR) then
        call handle_err(status,varnam,'cdfr_2d', 'nf_get_var_double')
        deallocate(temp)
        return
     endif
     ndim1 = min(dimlens(1), ldim(1))
     ndim2 = min(dimlens(2), ldim(2))
     varval(1:ndim1, 1:ndim2) = temp(1:ndim1, 1:ndim2)
     deallocate(temp)
  endif

!!$  st(1) = 1
!!$  cnt(1) = min(dimlens(1),ldim(1)) ! x count
!!$  cnt(2) = 1
!!$  do j = 1,min(dimlens(2),ldim(2)) ! For each Y : read slab into varval
!!$     st(2) = j              ! Start of slab
!!$     status = nf_get_vara_double(ncid,varid,st,cnt,varval(1,j))
!!$     if (status .ne. NF_NOERR) then
!!$        call handle_err(status,varnam,'cdfr_2d',                    &
!!$             &           'nf_get_vara_double')
!!$        return
!!$     end if
!!$  end do
  if (PRESENT (ier)) ier = 0
END SUBROUTINE cdfr_2d

SUBROUTINE cdfr_2c16(ncid,varnam,varval,ier)
  implicit none
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  COMPLEX(KIND=r8), dimension(:,:), intent(inout) ::  varval
  integer, optional,             intent(out) :: ier
  ! Local
  integer, dimension(2)   :: st,     cnt,    ldim
  integer                 :: varid,  status, j, i
  integer, dimension(2)   :: dimlens

  real(r8), dimension(:,:), allocatable :: temp
  integer ndim1, ndim2
 
  if (PRESENT (ier)) ier = 1
  ldim = shape(varval)
  ldim(1) = 2*ldim(1) ! Re/Im pairs
  call cdfgv(ncid,trim(varnam)//cmplx_name,varid,dimlens,ldim,'d',status)
  if (status .ne. 0) return
  !
  if(dimlens(1)==ldim(1).and.dimlens(2)==ldim(2))then
     status = nf_get_var_double(ncid, varid, varval)
     if(status/=NF_NOERR) then
        call handle_err(status,varnam,'cdfr_2c16', 'nf_get_var_double')
        return
     endif
  else
     allocate(temp(dimlens(1), dimlens(2)))
     status = nf_get_var_double(ncid, varid, temp)
     if(status/=NF_NOERR) then
        call handle_err(status,varnam,'cdfr_2c16', 'nf_get_var_double')
        deallocate(temp)
        return
     endif
     ndim1 = min(dimlens(1), ldim(1))
     ndim2 = min(dimlens(2), ldim(2))
     do i = 1, ndim1/2
        varval(i, 1:ndim2) = temp(2*(i-1)+1, 1:ndim2) + &
             & (0._r8, 1._r8)*temp(2*(i-1)+2, 1:ndim2)
     enddo
     deallocate(temp)
  endif
!!$  st(1) = 1
!!$  cnt(1) = min(dimlens(1),ldim(1)) ! x count
!!$  cnt(2) = 1
!!$  do j = 1,min(dimlens(2),ldim(2)) ! For each Y : read slab into varval
!!$     st(2) = j              ! Start of slab
!!$     status = nf_get_vara_double(ncid,varid,st,cnt,varval(1,j))
!!$     if (status .ne. NF_NOERR) then
!!$        call handle_err(status,varnam,'cdfr_2c16',                    &
!!$             &           'nf_get_vara_double')
!!$        return
!!$     end if
!!$  end do
  if (PRESENT (ier)) ier = 0
END SUBROUTINE cdfr_2c16

SUBROUTINE cdfr_2f(ncid,varnam,varval,ier)
  implicit none
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  REAL(KIND=r4), dimension(:,:), intent(inout) ::  varval
  integer, optional,             intent(out) :: ier
  ! Local
  integer, dimension(2)   :: st,     cnt,    ldim
  integer                 :: varid,  status, j
  integer, dimension(2)   :: dimlens

  real, dimension(:,:), allocatable :: temp
  integer ndim1, ndim2

  if (PRESENT (ier)) ier = 1
  ldim = shape(varval)
  call cdfgv(ncid,varnam,varid,dimlens,ldim,'r',status)
  if (status .ne. 0) return
  !
  if(dimlens(1)==ldim(1).and.dimlens(2)==ldim(2))then
     status = nf_get_var_real(ncid, varid, varval)
     if(status/=NF_NOERR) then
        call handle_err(status,varnam,'cdfr_2f', 'nf_get_var_real')
        return
     endif
  else
     allocate(temp(dimlens(1), dimlens(2)))
     status = nf_get_var_real(ncid, varid, temp)
     if(status/=NF_NOERR) then
        call handle_err(status,varnam,'cdfr_2f', 'nf_get_var_real')
        deallocate(temp)
        return
     endif
     ndim1 = min(dimlens(1), ldim(1))
     ndim2 = min(dimlens(2), ldim(2))
     varval(1:ndim1, 1:ndim2) = temp(1:ndim1, 1:ndim2)
     deallocate(temp)
  endif

!!$  st(1) = 1
!!$  cnt(1) = min(dimlens(1),ldim(1)) ! x count
!!$  cnt(2) = 1
!!$  do j = 1,min(dimlens(2),ldim(2)) ! For each Y : read slab into varval
!!$     st(2) = j              ! Start of slab
!!$     status = nf_get_vara_real(ncid,varid,st,cnt,varval(1,j))
!!$     if (status .ne. NF_NOERR) then
!!$        call handle_err(status,varnam,'cdfr_2f',                    &
!!$             &           'nf_get_vara_real')
!!$        return
!!$     end if
!!$  end do
  if (PRESENT (ier)) ier = 0
END SUBROUTINE cdfr_2f

SUBROUTINE cdfr_2c8(ncid,varnam,varval,ier)
  implicit none
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  COMPLEX(KIND=r4), dimension(:,:), intent(inout) ::  varval
  integer, optional,             intent(out) :: ier
  ! Local
  integer, dimension(2)   :: st,     cnt,    ldim
  integer                 :: varid,  status, j, i
  integer, dimension(2)   :: dimlens

  real(r8), dimension(:,:), allocatable :: temp
  integer ndim1, ndim2
 
  if (PRESENT (ier)) ier = 1
  ldim = shape(varval)
  ldim(1) = 2*ldim(1) ! Re/Pairs
  call cdfgv(ncid,trim(varnam)//cmplx_name,varid,dimlens,ldim,'r',status)
  if (status .ne. 0) return
  !
  if(dimlens(1)==ldim(1).and.dimlens(2)==ldim(2))then
     status = nf_get_var_real(ncid, varid, varval)
     if(status/=NF_NOERR) then
        call handle_err(status,varnam,'cdfr_2c8', 'nf_get_var_real')
        return
     endif
  else
     allocate(temp(dimlens(1), dimlens(2)))
     status = nf_get_var_real(ncid, varid, temp)
     if(status/=NF_NOERR) then
        call handle_err(status,varnam,'cdfr_2c8', 'nf_get_var_real')
        deallocate(temp)
        return
     endif
     ndim1 = min(dimlens(1), ldim(1))
     ndim2 = min(dimlens(2), ldim(2))
     do i = 1, ndim1/2
        varval(i, 1:ndim2) = temp(2*(i-1)+1, 1:ndim2) + &
             & (0.,1.)*temp(2*(i-1)+2, 1:ndim2)
     enddo
     deallocate(temp)
  endif
!!$  st(1) = 1
!!$  cnt(1) = min(dimlens(1),ldim(1)) ! x count
!!$  cnt(2) = 1
!!$  do j = 1,min(dimlens(2),ldim(2)) ! For each Y : read slab into varval
!!$     st(2) = j              ! Start of slab
!!$     status = nf_get_vara_real(ncid,varid,st,cnt,varval(1,j))
!!$     if (status .ne. NF_NOERR) then
!!$        call handle_err(status,varnam,'cdfr_2c8',                    &
!!$             &           'nf_get_vara_real')
!!$        return
!!$     end if
!!$  end do
  if (PRESENT (ier)) ier = 0
END SUBROUTINE cdfr_2c8

SUBROUTINE cdfr_2c(ncid,varnam,varval,ier)
  implicit none
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  character*(*), dimension(:),intent(inout) :: varval
  integer, optional,          intent(out) :: ier
  ! Local
  integer, dimension(2)   :: st,     cnt,    ldim
  integer                 :: varid,  status, charlen, j
  integer, dimension(2)   :: dimlens
  if (PRESENT (ier)) ier = 1
  ldim(1) = len(varval)
  ldim(2)=size(varval)
  call cdfgv(ncid,varnam,varid,dimlens,ldim,'c',status)
  if (status .ne. 0) return
  !
  st(1) = 1
  cnt(1) = min(dimlens(1),ldim(1)) ! x count
  cnt(2) = 1
  do j = 1,min(dimlens(2),ldim(2)) ! For each Y : read slab into varval
     st(2) = j              ! Start of slab
     status = nf_get_vara_text(ncid,varid,st,cnt,varval(j))
     if (status .ne. NF_NOERR) then
        call handle_err(status,varnam,'cdfr_2c',                       &
             &        'nf_get_var_text')
        return
     end if
  end do
  if (PRESENT (ier)) ier = 0
END SUBROUTINE cdfr_2c
 
SUBROUTINE cdfr_1i(ncid,varnam,varval,ier)
  implicit none
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  integer, dimension(:), intent(inout) :: varval
  integer, optional,     intent(out) :: ier
  ! Local
  integer                 :: varid,  status
  integer, dimension(1)   :: dimlens,ldim
  if (PRESENT (ier)) ier = 1
  ldim = shape(varval)
  call cdfgv(ncid,varnam,varid,dimlens,ldim,'i',status)
  if (status .ne. 0) return
  !
  status = nf_get_var_int(ncid,varid,varval)
  call handle_err(status,varnam,'cdfr_1i','nf_get_var_int')
  if (PRESENT (ier)) ier = status
END SUBROUTINE cdfr_1i
 
SUBROUTINE cdfr_1l(ncid,varnam,varval,ier)
  ! Read 1 dimensional logical array
  !
  implicit none
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  logical, dimension(:), intent(inout) :: varval
  integer, optional,     intent(out) :: ier
  ! Local
  integer                            :: status
  integer, dimension(:), allocatable :: varval_i
  character*11, parameter :: logical_name = '__logical__'
 
  ALLOCATE (varval_i(size(varval,1)), stat=status)
  if (status .ne. 0) STOP 'Allocation error in cdf_getvar'
 
  call cdfr_1i(ncid,trim(varnam)//logical_name,varval_i,ier)
 
  WHERE (varval_i == 0)
     varval = .false.
  ELSEWHERE
     varval = .true.
  END WHERE
  DEALLOCATE (varval_i)
END SUBROUTINE cdfr_1l
 
SUBROUTINE cdfr_1d(ncid,varnam,varval,ier)
  implicit none
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  REAL(KIND=r8), dimension(:), intent(inout) :: varval
  integer, optional,           intent(out) :: ier
  ! Local
  integer                 :: varid,  status
  integer, dimension(1)   :: dimlens,ldim
  if (PRESENT (ier)) ier = 1
  ldim = shape(varval)
  call cdfgv(ncid,varnam,varid,dimlens,ldim,'d',status)
  if (status .ne. 0) return
  !
  status = nf_get_var_double(ncid,varid,varval)
  call handle_err(status,varnam,'cdfr_1d','nf_get_var_double')
  if (PRESENT (ier)) ier = status
END SUBROUTINE cdfr_1d

SUBROUTINE cdfr_1c16(ncid,varnam,varval,ier)
  implicit none
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  COMPLEX(KIND=r8), dimension(:), intent(inout) :: varval
  integer, optional,           intent(out) :: ier
 
  ! Local
  integer                 :: varid,  status
  integer, dimension(1)   :: dimlens,ldim
  if (PRESENT (ier)) ier = 1
  ldim = shape(varval)
  ldim(1) = 2*ldim(1) ! Re/Im pairs
  call cdfgv(ncid,trim(varnam)//cmplx_name,varid,dimlens,ldim,'d',status)
  if (status .ne. 0) return
  !
  status = nf_get_var_double(ncid,varid,varval)
  call handle_err(status,varnam,'cdfr_1c16','nf_get_var_double')
  if (PRESENT (ier)) ier = status
END SUBROUTINE cdfr_1c16

SUBROUTINE cdfr_1f(ncid,varnam,varval,ier)
  implicit none
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  REAL(KIND=r4), dimension(:), intent(inout) :: varval
  integer, optional,           intent(out) :: ier
  ! Local
  integer                 :: varid,  status
  integer, dimension(1)   :: dimlens,ldim
  if (PRESENT (ier)) ier = 1
  ldim = shape(varval)
  call cdfgv(ncid,varnam,varid,dimlens,ldim,'r',status)
  if (status .ne. 0) return
  !
  status = nf_get_var_real(ncid,varid,varval)
  call handle_err(status,varnam,'cdfr_1f','nf_get_var_real')
  if (PRESENT (ier)) ier = status
END SUBROUTINE cdfr_1f

SUBROUTINE cdfr_1c8(ncid,varnam,varval,ier)
  implicit none
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  COMPLEX(KIND=r4), dimension(:), intent(inout) :: varval
  integer, optional,           intent(out) :: ier
  ! Local
  integer                 :: varid,  status
  integer, dimension(1)   :: dimlens,ldim
 
  if (PRESENT (ier)) ier = 1
  ldim = shape(varval)
  ldim(1) = 2*ldim(1)
  call cdfgv(ncid,trim(varnam)//cmplx_name,varid,dimlens,ldim,'r',status)
  if (status .ne. 0) return
  !
  status = nf_get_var_real(ncid,varid,varval)
  call handle_err(status,varnam,'cdfr_1c8','nf_get_var_real')
  if (PRESENT (ier)) ier = status
END SUBROUTINE cdfr_1c8

SUBROUTINE cdfr_1c(ncid,varnam,varval,ier)
  implicit none
  ! Input
  integer,       intent(in)  :: ncid
  character*(*), intent(in)  :: varnam
  ! Output
  character*(*),     intent(inout) :: varval
  integer, optional, intent(out) :: ier
  ! Local
  integer                 :: varid,  status
  integer, dimension(1)   :: dimlens,ldim
  if (PRESENT (ier)) ier = 1
  ldim(1) = len(varval)
  call cdfgv(ncid,varnam,varid,dimlens,ldim,'c',status)
  if (status .ne. 0) return
  !
  status = nf_get_var_text(ncid,varid,varval)
  call handle_err(status,varnam,'cdfr_1c','nf_get_var_text')
  if (PRESENT (ier)) ier = status
END SUBROUTINE cdfr_1c
 
SUBROUTINE cdfr_0i(ncid,varnam,varval,ier)
  implicit none
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  integer,           intent(out) ::  varval
  integer, optional, intent(out) :: ier
  ! Local
  integer                 :: varid,  status
  integer, dimension(1)   :: dimlens,ldim
  if (PRESENT (ier)) ier = 1
  ldim(1) = 0
  call cdfgv(ncid,varnam,varid,dimlens,ldim,'i',status)
  if (status .ne. 0) return
  !
  status = nf_get_var_int(ncid,varid,varval)
  call handle_err(status,varnam,'cdfr_0i','nf_get_var_int')
  if (PRESENT (ier)) ier = status
END SUBROUTINE cdfr_0i
 
SUBROUTINE cdfr_0l(ncid,varnam,varval,ier)
  ! Read scalar logical array
  !
  implicit none
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  logical,           intent(out) :: varval
  integer, optional, intent(out) :: ier
  ! Local
  integer :: varval_i
  character*11, parameter :: logical_name = '__logical__'
 
  call cdfr_0i(ncid,trim(varnam)//logical_name,varval_i,ier)
 
  IF (varval_i == 0) THEN
     varval = .false.
  ELSE
     varval = .true.
  END IF
END SUBROUTINE cdfr_0l
 
SUBROUTINE cdfr_0d(ncid,varnam,varval,ier)
  implicit none
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  REAL(KIND=r8),     intent(out) :: varval
  integer, optional, intent(out) :: ier
  ! Local
  integer                 :: varid,  status
  integer, dimension(1)   :: dimlens,ldim
  if (PRESENT (ier)) ier = 1
  ldim(1) = 0
  call cdfgv(ncid,varnam,varid,dimlens,ldim,'d',status)
  if (status .ne. 0) return
  !
  status = nf_get_var_double(ncid,varid,varval)
  call handle_err(status,varnam,'cdfr_0d','nf_get_var_double')
  if (PRESENT (ier)) ier = status
END SUBROUTINE cdfr_0d

SUBROUTINE cdfr_0c16(ncid,varnam,varval,ier)
  implicit none
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  COMPLEX(KIND=r8),     intent(out) :: varval
  integer, optional, intent(out) :: ier
  ! Local
  integer                 :: varid,  status
  integer, dimension(1)   :: dimlens,ldim
 
  if (PRESENT (ier)) ier = 1
  ldim(1) = 2 ! Re/Im pair
  call cdfgv(ncid,trim(varnam)//cmplx_name,varid,dimlens,ldim,'d',status)
  if (status .ne. 0) return
  !
  status = nf_get_var_double(ncid,varid,varval)
  call handle_err(status,varnam,'cdfr_0d','nf_get_var_double')
  if (PRESENT (ier)) ier = status
END SUBROUTINE cdfr_0c16

SUBROUTINE cdfr_0f(ncid,varnam,varval,ier)
  implicit none
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  REAL(KIND=r4),     intent(out) :: varval
  integer, optional, intent(out) :: ier
  ! Local
  integer                 :: varid,  status
  integer, dimension(1)   :: dimlens,ldim
  if (PRESENT (ier)) ier = 1
  ldim(1) = 0
  call cdfgv(ncid,varnam,varid,dimlens,ldim,'r',status)
  if (status .ne. 0) return
  !
  status = nf_get_var_real(ncid,varid,varval)
  call handle_err(status,varnam,'cdfr_0f','nf_get_var_real')
  if (PRESENT (ier)) ier = status
END SUBROUTINE cdfr_0f

SUBROUTINE cdfr_0c8(ncid,varnam,varval,ier)
  implicit none
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  COMPLEX(KIND=r4),     intent(out) :: varval
  integer, optional, intent(out) :: ier
  ! Local
  integer                 :: varid,  status
  integer, dimension(1)   :: dimlens,ldim
 
  if (PRESENT (ier)) ier = 1
  ldim(1) = 2 ! Re/Im pair
  call cdfgv(ncid,trim(varnam)//cmplx_name,varid,dimlens,ldim,'r',status)
  if (status .ne. 0) return
  !
  status = nf_get_var_real(ncid,varid,varval)
  call handle_err(status,varnam,'cdfr_0f','nf_get_var_real')
  if (PRESENT (ier)) ier = status
END SUBROUTINE cdfr_0c8
#endif
END MODULE ezcdf_GenGet
