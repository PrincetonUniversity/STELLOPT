MODULE ezcdf_GenPut
  USE ezcdf_opncls
  USE ezcdf_inqvar

  ! Generic Interface to Write netcdf data Variables
  ! 03/10/99 C. Ludescher
  ! C. Ludescher/A. Pletzer Tue Apr  4 10:11:33 EDT 2000
  ! + support for complex numbers (ap) Wed May 16 15:18:05 EDT 2001
  ! added support for logicals (sph) Oct 1, 2002
  IMPLICIT NONE
!DEC$ IF DEFINED (NETCDF)

  !
  ! C. Ludescher/A. Pletzer Tue Apr  4 10:11:33 EDT 2000
  !
  ! added support (ap) Wed May 16 15:06:34 EDT 2001
  PUBLIC :: cdfw_3i, cdfw_3l, cdfw_3d, cdfw_3c16, cdfw_3f, cdfw_3c8,  &
       & cdfw_2i, cdfw_2l, cdfw_2d, cdfw_2c16, cdfw_2f, cdfw_2c8, cdfw_2c, &
       & cdfw_1i, cdfw_1l, cdfw_1d, cdfw_1c16, cdfw_1f, cdfw_1c8, cdfw_1c, &
       & cdfw_0i, cdfw_0l, cdfw_0d, cdfw_0c16, cdfw_0f, cdfw_0c8,          &
       & cdfPutVar, cdf_write

  ! cdf_write is an alias of cdfPutVar (required by ifc compiler)
  INTERFACE cdf_write
     MODULE PROCEDURE cdfw_3i, cdfw_3l, cdfw_3d, cdfw_3c16, cdfw_3f, cdfw_3c8, &
          cdfw_2i, cdfw_2l, cdfw_2d, cdfw_2c16, cdfw_2f, cdfw_2c8, cdfw_2c, &
          cdfw_1i, cdfw_1l, cdfw_1d, cdfw_1c16, cdfw_1f, cdfw_1c8, cdfw_1c, &
          cdfw_0i, cdfw_0l, cdfw_0d, cdfw_0c16, cdfw_0f, cdfw_0c8
  END INTERFACE

  ! same as above (Intel compiler does not handle well aliases)
  INTERFACE cdfPutVar
     MODULE PROCEDURE cdfw_3i, cdfw_3l, cdfw_3d, cdfw_3c16, cdfw_3f, cdfw_3c8, &
          cdfw_2i, cdfw_2l, cdfw_2d, cdfw_2c16, cdfw_2f, cdfw_2c8, cdfw_2c, &
          cdfw_1i, cdfw_1l, cdfw_1d, cdfw_1c16, cdfw_1f, cdfw_1c8, cdfw_1c, &
          cdfw_0i, cdfw_0l, cdfw_0d, cdfw_0c16, cdfw_0f, cdfw_0c8
  END INTERFACE

  PUBLIC :: cdfDefVar, cdf_define,                                          &
       & cdfd_3i, cdfd_3l, cdfd_3d, cdfd_3c16, cdfd_3f, cdfd_3c8,          &
       & cdfd_2i, cdfd_2l, cdfd_2d, cdfd_2c16, cdfd_2f, cdfd_2c8, cdfd_2c, &
       & cdfd_1i, cdfd_1l, cdfd_1d, cdfd_1c16, cdfd_1f, cdfd_1c8, cdfd_1c, &
       & cdfd_0i, cdfd_0l, cdfd_0d, cdfd_0c16, cdfd_0f, cdfd_0c8

  INTERFACE cdf_define
     MODULE PROCEDURE cdfDefVar,                                            &
          cdfd_3i, cdfd_3l, cdfd_3d, cdfd_3c16, cdfd_3f, cdfd_3c8,          &
          cdfd_2i, cdfd_2l, cdfd_2d, cdfd_2c16, cdfd_2f, cdfd_2c8, cdfd_2c, &
          cdfd_1i, cdfd_1l, cdfd_1d, cdfd_1c16, cdfd_1f, cdfd_1c8, cdfd_1c, &
          cdfd_0i, cdfd_0l, cdfd_0d, cdfd_0c16, cdfd_0f, cdfd_0c8
  END INTERFACE

  PRIVATE
  include "netcdf.inc"

  INTEGER, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
  INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
  CHARACTER*13, parameter :: cmplx_name = '__CmPlx_Re_Im'

  CHARACTER(len=nf_max_name)  :: varnam_noalpha
  PRIVATE :: r4, r8, cmplx_name, varnam_noalpha

  EXTERNAL handle_err


CONTAINS
!-----------------------------------------------------------------
!cdfPutVar implementation routines

  SUBROUTINE cdfw_3i(ncid,varnam,varval,ier)
    !     Write 3 dimensional Integer data array
    !     Use cdfDefVar to define the Variable
    implicit none
    ! Input
    integer,                   intent(in) :: ncid
    character*(*),             intent(in) :: varnam
    integer, dimension(:,:,:), intent(in) :: varval
    ! Output
    integer, optional,        intent(out) :: ier
    ! Local
    integer, dimension(3) :: st = (/1,1,1/)
    integer, dimension(3) :: dimlens
    integer               :: varid, ndims, sts
    if (PRESENT (ier)) ier = 1
    sts = nf_enddef(ncid)
    call cdfInqV(ncid,varnam,varid,dimlens,ndims,sts)
    if (sts .ne. 0) return
    if (ndims .ne. 3) then
       print "('% cdfPutVar_3i: --E-- The variable ',a,               &
            &         ' was defined as',i2,' dimensional')",varnam,ndims
       return
    end if
    sts = nf_put_vara_int(ncid,varid,st,dimlens,varval)
    call handle_err(sts,varnam,'cdfPutVar_3i','nf_put_vara_int')
    if (PRESENT (ier)) ier = sts
  END SUBROUTINE cdfw_3i

  SUBROUTINE cdfw_3l(ncid,varnam,varval,ier)
    !     Write 3 dimensional logical array
    !     Use cdfDefVar to define the Variable
    implicit none
    ! Input
    integer,                   intent(in) :: ncid
    character*(*),             intent(in) :: varnam
    logical, dimension(:,:,:), intent(in) :: varval
    ! Output
    integer, optional,         intent(out) :: ier
    ! Local
    integer                                :: sts
    integer, allocatable, dimension(:,:,:) :: varval_i
    character*11, parameter :: logical_name = '__logical__'


    ALLOCATE (varval_i(size(varval,1), size(varval,2), size(varval,3)), stat=sts)
    if (sts .ne. 0) stop 'Allocation error in cdf_putvar'

    WHERE (varval) 
       varval_i = 1
    ELSEWHERE 
       varval_i = 0
    END WHERE

    call cdfw_3i(ncid,trim(varnam)//logical_name,varval_i,ier)

    DEALLOCATE (varval_i)
  END SUBROUTINE cdfw_3l

  SUBROUTINE cdfw_3d(ncid,varnam,varval,ier)
    !     Write 3 dimensional 64-bit Real data array
    !     Use cdfDefVar to define the Variable
    implicit none
    ! Input
    integer,                         intent(in) :: ncid
    character*(*),                   intent(in) :: varnam
    REAL(KIND=r8), dimension(:,:,:), intent(in) :: varval
    ! Output
    integer, optional,      intent(out) :: ier
    ! Local
    integer, dimension(3) :: st = (/1,1,1/)
    integer, dimension(3) :: dimlens
    integer               :: varid, ndims, sts
    if (PRESENT (ier)) ier = 1
    sts= nf_enddef(ncid)
    call cdfInqV(ncid,varnam,varid,dimlens,ndims,sts)
    if (sts .ne. 0)  return
    if (ndims .ne. 3) then
       print "('% cdfPutVar_3d: --E-- The variable ',a,               &
            &         ' was defined as',i2,' dimensional')",varnam,ndims
       return
    end if
    sts = nf_put_vara_double(ncid,varid,st,dimlens,varval)
    call handle_err(sts,varnam,'cdf_3d','nf_put_vara_double')
    if (PRESENT (ier)) ier = sts
  END SUBROUTINE cdfw_3d

  SUBROUTINE cdfw_3c16(ncid,varnam,varval,ier)
    !     Write 3 dimensional 128-bit complex data array
    !     Use cdfDefVar to define the Variable
    implicit none
    ! Input
    integer,                         intent(in) :: ncid
    character*(*),                   intent(in) :: varnam
    COMPLEX(KIND=r8), dimension(:,:,:), intent(in) :: varval
    ! Output
    integer, optional,      intent(out) :: ier
    ! Local
    integer, dimension(3) :: st = (/1,1,1/)
    integer, dimension(3) :: dimlens
    integer               :: varid, ndims, sts

    if (PRESENT (ier)) ier = 1
    sts= nf_enddef(ncid)
    call cdfInqV(ncid,trim(varnam)//cmplx_name,varid,dimlens,ndims,sts)
    if (sts .ne. 0)  return
    if (ndims .ne. 3) then
       print "('% cdfPutVar_3c16: --E-- The variable ',a,               &
            &         ' was defined as',i2,' dimensional')",varnam,ndims
       return
    end if
    sts = nf_put_vara_double(ncid,varid,st,dimlens,varval)
    call handle_err(sts,varnam,'cdf_3c16','nf_put_vara_double')
    if (PRESENT (ier)) ier = sts
  END SUBROUTINE cdfw_3c16

  SUBROUTINE cdfw_3f(ncid,varnam,varval,ier)
    !     Write 3 dimensional default Real data array
    !     Use cdfDefVar to define the Variable
    implicit none
    ! Input
    integer,                         intent(in) :: ncid
    character*(*),                   intent(in) :: varnam
    REAL(KIND=r4), dimension(:,:,:), intent(in) :: varval
    ! Output
    integer, optional,      intent(out) :: ier
    ! Local
    integer, dimension(3) :: st = (/1,1,1/)
    integer, dimension(3) :: dimlens
    integer               :: varid, ndims, sts
    if (PRESENT (ier)) ier = 1
    sts= nf_enddef(ncid)
    call cdfInqV(ncid,varnam,varid,dimlens,ndims,sts)
    if (sts .ne. 0)  return
    if (ndims .ne. 3) then
       print "('% cdfPutVar_3f: --E-- The variable ',a,               &
            &         ' was defined as',i2,' dimensional')",varnam,ndims
       return
    end if
    sts = nf_put_vara_real(ncid,varid,st,dimlens,varval)
    call handle_err(sts,varnam,'cdf_3f','nf_put_vara_real')
    if (PRESENT (ier)) ier = sts
  END SUBROUTINE cdfw_3f

  SUBROUTINE cdfw_3c8(ncid,varnam,varval,ier)
    !     Write 3 dimensional complex*8 data array
    !     Use cdfDefVar to define the Variable
    implicit none
    ! Input
    integer,                         intent(in) :: ncid
    character*(*),                   intent(in) :: varnam
    COMPLEX(KIND=r4), dimension(:,:,:), intent(in) :: varval
    ! Output
    integer, optional,      intent(out) :: ier
    ! Local
    integer, dimension(3) :: st = (/1,1,1/)
    integer, dimension(3) :: dimlens
    integer               :: varid, ndims, sts

    if (PRESENT (ier)) ier = 1
    sts= nf_enddef(ncid)
    call cdfInqV(ncid,trim(varnam)//cmplx_name,varid,dimlens,ndims,sts)
    if (sts .ne. 0)  return
    if (ndims .ne. 3) then
       print "('% cdfPutVar_3c8: --E-- The variable ',a,               &
            &         ' was defined as',i2,' dimensional')",varnam,ndims
       return
    end if
    sts = nf_put_vara_real(ncid,varid,st,dimlens,varval)
    call handle_err(sts,varnam,'cdf_3c8','nf_put_vara_real')
    if (PRESENT (ier)) ier = sts
  END SUBROUTINE cdfw_3c8

  SUBROUTINE cdfw_2i(ncid,varnam,varval,ier)
    !     Write 2 dimensional Integer data array
    !     Use cdfDefVar to define the Variable
    implicit none
    ! Input
    integer,                 intent(in) :: ncid
    character*(*),           intent(in) :: varnam
    integer, dimension(:,:), intent(in) :: varval
    ! Output
    integer, optional,      intent(out) :: ier
    ! Local
    integer, dimension(2) :: st = (/1,1/)
    integer, dimension(2) :: dimlens
    integer               :: varid, ndims, sts
    if (PRESENT (ier)) ier = 1
    sts = nf_enddef(ncid)
    call cdfInqV(ncid,varnam,varid,dimlens,ndims,sts)
    if (sts .ne. 0) return
    if (ndims .ne. 2) then
       print "('% cdfPutVar_2i: --E-- The variable ',a,               &
            &         ' was defined as',i2,' dimensional')",varnam,ndims
       return
    end if
    sts = nf_put_vara_int(ncid,varid,st,dimlens,varval)
    call handle_err(sts,varnam,'cdfPutVar_2i','nf_put_vara_int')
    if (PRESENT (ier)) ier = sts
  END SUBROUTINE cdfw_2i

  SUBROUTINE cdfw_2l(ncid,varnam,varval,ier)
    !     Write 2 dimensional logical array
    !     Use cdfDefVar to define the Variable
    implicit none
    ! Input
    integer,                 intent(in) :: ncid
    character*(*),           intent(in) :: varnam
    logical, dimension(:,:), intent(in) :: varval
    ! Output
    integer, optional,       intent(out) :: ier
    ! Local
    integer                              :: sts
    integer, allocatable, dimension(:,:) :: varval_i
    character*11, parameter :: logical_name = '__logical__'

    ALLOCATE (varval_i(size(varval,1), size(varval,2)), stat=sts)
    if (sts .ne. 0) stop 'Allocation error in cdf_putvar'

    WHERE (varval) 
       varval_i = 1
    ELSEWHERE 
       varval_i = 0
    END WHERE

    call cdfw_2i(ncid,trim(varnam)//logical_name,varval_i,ier)

    DEALLOCATE (varval_i)
  END SUBROUTINE cdfw_2l

  SUBROUTINE cdfw_2d(ncid,varnam,varval,ier)
    !     Write 2 dimensional 64-bit Real data array
    !     Use cdfDefVar to define the Variable
    implicit none
    ! Input
    integer,                       intent(in) :: ncid
    character*(*),                 intent(in) :: varnam
    REAL(KIND=r8), dimension(:,:), intent(in) :: varval
    ! Output
    integer, optional,      intent(out) :: ier
    ! Local
    integer, dimension(2) :: st = (/1,1/)
    integer, dimension(2) :: dimlens
    integer               :: varid, ndims, sts
    if (PRESENT (ier)) ier = 1
    sts= nf_enddef(ncid)
    call cdfInqV(ncid,varnam,varid,dimlens,ndims,sts)
    if (sts .ne. 0)  return
    if (ndims .ne. 2) then
       print "('% cdfPutVar_2d: --E-- The variable ',a,               &
            &         ' was defined as',i2,' dimensional')",varnam,ndims
       return
    end if
    sts = nf_put_vara_double(ncid,varid,st,dimlens,varval)
    call handle_err(sts,varnam,'cdf_2d','nf_put_vara_double')
    if (PRESENT (ier)) ier = sts
  END SUBROUTINE cdfw_2d

  SUBROUTINE cdfw_2c16(ncid,varnam,varval,ier)
    !     Write 2 dimensional 128-bit complex data array
    !     Use cdfDefVar to define the Variable
    implicit none
    ! Input
    integer,                       intent(in) :: ncid
    character*(*),                 intent(in) :: varnam
    COMPLEX(KIND=r8), dimension(:,:), intent(in) :: varval
    ! Output
    integer, optional,      intent(out) :: ier
    ! Local
    integer, dimension(2) :: st = (/1,1/)
    integer, dimension(2) :: dimlens
    integer               :: varid, ndims, sts

    if (PRESENT (ier)) ier = 1
    sts= nf_enddef(ncid)
    call cdfInqV(ncid,trim(varnam)//cmplx_name,varid,dimlens,ndims,sts)
    if (sts .ne. 0)  return
    if (ndims .ne. 2) then
       print "('% cdfPutVar_2c16: --E-- The variable ',a,               &
            &         ' was defined as',i2,' dimensional')",varnam,ndims
       return
    end if
    sts = nf_put_vara_double(ncid,varid,st,dimlens,varval)
    call handle_err(sts,varnam,'cdf_2c16','nf_put_vara_double')
    if (PRESENT (ier)) ier = sts
  END SUBROUTINE cdfw_2c16

  SUBROUTINE cdfw_2f(ncid,varnam,varval,ier)
    !     Write 2 dimensional default Real data array
    !     Use cdfDefVar to define the Variable
    implicit none
    ! Input
    integer,                       intent(in) :: ncid
    character*(*),                 intent(in) :: varnam
    REAL(KIND=r4), dimension(:,:), intent(in) :: varval
    ! Output
    integer, optional,      intent(out) :: ier
    ! Local
    integer, dimension(2) :: st = (/1,1/)
    integer, dimension(2) :: dimlens
    integer               :: varid, ndims, sts
    if (PRESENT (ier)) ier = 1
    sts= nf_enddef(ncid)
    call cdfInqV(ncid,varnam,varid,dimlens,ndims,sts)
    if (sts .ne. 0)  return
    if (ndims .ne. 2) then
       print "('% cdfPutVar_2f: --E-- The variable ',a,               &
            &         ' was defined as',i2,' dimensional')",varnam,ndims
       return
    end if
    sts = nf_put_vara_real(ncid,varid,st,dimlens,varval)
    call handle_err(sts,varnam,'cdf_2f','nf_put_vara_real')
    if (PRESENT (ier)) ier = sts
  END SUBROUTINE cdfw_2f

  SUBROUTINE cdfw_2c8(ncid,varnam,varval,ier)
    !     Write 2 dimensional complex*8 data array
    !     Use cdfDefVar to define the Variable
    implicit none
    ! Input
    integer,                       intent(in) :: ncid
    character*(*),                 intent(in) :: varnam
    COMPLEX(KIND=r4), dimension(:,:), intent(in) :: varval
    ! Output
    integer, optional,      intent(out) :: ier
    ! Local
    integer, dimension(2) :: st = (/1,1/)
    integer, dimension(2) :: dimlens
    integer               :: varid, ndims, sts

    if (PRESENT (ier)) ier = 1
    sts= nf_enddef(ncid)
    call cdfInqV(ncid,trim(varnam)//cmplx_name,varid,dimlens,ndims,sts)
    if (sts .ne. 0)  return
    if (ndims .ne. 2) then
       print "('% cdfPutVar_2c8: --E-- The variable ',a,               &
            &         ' was defined as',i2,' dimensional')",varnam,ndims
       return
    end if
    sts = nf_put_vara_real(ncid,varid,st,dimlens,varval)
    call handle_err(sts,varnam,'cdf_2c8','nf_put_vara_real')
    if (PRESENT (ier)) ier = sts
  END SUBROUTINE cdfw_2c8

  SUBROUTINE cdfw_2c(ncid,varnam,varval,ier)
    !     Write 2 dimensional character array
    !     Use cdfDefVar to define the Variable
    implicit none
    ! Input
    integer,                     intent(in) :: ncid
    character*(*),               intent(in) :: varnam
    character*(*), dimension(:), intent(in) :: varval
    ! Output
    integer, optional,      intent(out) :: ier
    ! Local
    integer, dimension(2) :: st = (/1,1/)
    integer, dimension(2) :: dimlens
    integer               :: varid, ndims, sts
    if (PRESENT (ier)) ier = 1
    sts = nf_enddef(ncid)
    call cdfInqV(ncid,varnam,varid,dimlens,ndims,sts)
    if (sts .ne. 0) return
    if (ndims .ne. 2) then
       print "('% cdfPutVar_2c: --E-- The variable ',a,               &
            &         ' was defined as',i2,' dimensional')",varnam,ndims
       return
    end if
    sts = nf_put_vara_text(ncid,varid,st,dimlens,varval)
    call handle_err(sts,varnam,'cdfPutVar_2c','nf_put_var_text')
    if (PRESENT (ier)) ier = sts
  END SUBROUTINE cdfw_2c

  SUBROUTINE cdfw_1i(ncid,varnam,varval,ier)
    !     write 1 dimensional Integer array
    implicit none
    ! Input
    integer,                 intent(in) :: ncid
    character*(*),           intent(in) :: varnam
    integer, dimension(:),   intent(in) :: varval
    ! Output
    integer, optional,      intent(out) :: ier
    ! Local
    integer, dimension(1) :: dimlens
    integer               :: varid, ndims, sts
    if (PRESENT (ier)) ier = 1
    sts = nf_enddef(ncid)
    call cdfInqV(ncid,varnam,varid,dimlens,ndims,sts)
    if (sts .ne. 0) return
    if (ndims .ne. 1) then
       print "('% cdfPutVar_1i: --E-- The variable ',a,               &
            &         ' was defined as',i2,' dimensional')",varnam,ndims
       return
    end if
    sts = nf_put_var_int(ncid,varid,varval)
    call handle_err(sts,varnam,'cdfPutVar_1i','nf_put_var_int')
    if (PRESENT (ier)) ier = sts
  END SUBROUTINE cdfw_1i

  SUBROUTINE cdfw_1l(ncid,varnam,varval,ier)
    !     Write 1 dimensional logical array
    implicit none
    ! Input
    integer,               intent(in) :: ncid
    character*(*),         intent(in) :: varnam
    logical, dimension(:), intent(in) :: varval
    ! Output
    integer, optional,     intent(out) :: ier
    ! Local
    integer                            :: sts
    integer, allocatable, dimension(:) :: varval_i
    character*11, parameter :: logical_name = '__logical__'

    ALLOCATE (varval_i(size(varval,1)), stat=sts)
    if (sts .ne. 0) stop 'Allocation error in cdf_putvar'

    WHERE (varval) 
       varval_i = 1
    ELSEWHERE 
       varval_i = 0
    END WHERE

    call cdfw_1i(ncid,trim(varnam)//logical_name,varval_i,ier)

    DEALLOCATE (varval_i)
  END SUBROUTINE cdfw_1l

  SUBROUTINE cdfw_1d(ncid,varnam,varval,ier)
    !     Write 1 dimensional 64-bit data array
    !     Use cdfDefVar to define the Variable
    !
    implicit none
    ! Input
    integer,                     intent(in) :: ncid
    character*(*),               intent(in) :: varnam
    REAL(KIND=r8), dimension(:), intent(in) :: varval
    ! Output
    integer, optional,      intent(out) :: ier
    ! Local
    integer, dimension(1) :: dimlens
    integer               :: varid, ndims, sts
    if (PRESENT (ier)) ier = 1
    sts = nf_enddef(ncid)
    call cdfInqV(ncid,varnam,varid,dimlens,ndims,sts)
    if (sts .ne. 0) return
    if (ndims .ne. 1) then
       print "('% cdfPutVar_1d: --E-- The variable ',a,               &
            &         ' was defined as',i2,' dimensional')",varnam,ndims
       return
    end if
    sts = nf_put_var_double(ncid,varid,varval)
    call handle_err(sts,varnam,'cdfPurVar_1d','nf_put_var_double')
    if (PRESENT (ier)) ier = sts
  END SUBROUTINE cdfw_1d

  SUBROUTINE cdfw_1c16(ncid,varnam,varval,ier)
    !     Write 1 dimensional 128-bit complex data array
    !     Use cdfDefVar to define the Variable
    !
    implicit none
    ! Input
    integer,                     intent(in) :: ncid
    character*(*),               intent(in) :: varnam
    COMPLEX(KIND=r8), dimension(:), intent(in) :: varval
    ! Output
    integer, optional,      intent(out) :: ier
    ! Local
    integer, dimension(1) :: dimlens
    integer               :: varid, ndims, sts

    if (PRESENT (ier)) ier = 1
    sts = nf_enddef(ncid)
    call cdfInqV(ncid,trim(varnam)//cmplx_name,varid,dimlens,ndims,sts)
    if (sts .ne. 0) return
    if (ndims .ne. 1) then
       print "('% cdfPutVar_1c16: --E-- The variable ',a,               &
            &         ' was defined as',i2,' dimensional')",varnam,ndims
       return
    end if
    sts = nf_put_var_double(ncid,varid,varval)
    call handle_err(sts,varnam,'cdfPurVar_1c16','nf_put_var_double')
    if (PRESENT (ier)) ier = sts
  END SUBROUTINE cdfw_1c16

  SUBROUTINE cdfw_1f(ncid,varnam,varval,ier)
    !     Write 1 dimensional default real array
    !     Use cdfDefVar to define the Variable
    !
    implicit none
    ! Input
    integer,                     intent(in) :: ncid
    character*(*),               intent(in) :: varnam
    REAL(KIND=r4), dimension(:), intent(in) :: varval
    ! Output
    integer, optional,      intent(out) :: ier
    ! Local
    integer, dimension(1) :: dimlens
    integer               :: varid, ndims, sts
    if (PRESENT (ier)) ier = 1
    sts = nf_enddef(ncid)
    call cdfInqV(ncid,varnam,varid,dimlens,ndims,sts)
    if (sts .ne. 0) return
    if (ndims .ne. 1) then
       print "('% cdfPutVar_1f: --E-- The variable ',a,               &
            &         ' was defined as',i2,' dimensional')",varnam,ndims
       return
    end if
    sts = nf_put_var_real(ncid,varid,varval)
    call handle_err(sts,varnam,'cdfPurVar_1f','nf_put_var_real')
    if (PRESENT (ier)) ier = sts
  END SUBROUTINE cdfw_1f

  SUBROUTINE cdfw_1c8(ncid,varnam,varval,ier)
    !     Write 1 dimensional complex*8 array
    !     Use cdfDefVar to define the Variable
    !
    implicit none
    ! Input
    integer,                     intent(in) :: ncid
    character*(*),               intent(in) :: varnam
    COMPLEX(KIND=r4), dimension(:), intent(in) :: varval
    ! Output
    integer, optional,      intent(out) :: ier
    ! Local
    integer, dimension(1) :: dimlens
    integer               :: varid, ndims, sts

    if (PRESENT (ier)) ier = 1
    sts = nf_enddef(ncid)
    call cdfInqV(ncid,trim(varnam)//cmplx_name,varid,dimlens,ndims,sts)
    if (sts .ne. 0) return
    if (ndims .ne. 1) then
       print "('% cdfPutVar_1c8: --E-- The variable ',a,               &
            &         ' was defined as',i2,' dimensional')",varnam,ndims
       return
    end if
    sts = nf_put_var_real(ncid,varid,varval)
    call handle_err(sts,varnam,'cdfPurVar_1c8','nf_put_var_real')
    if (PRESENT (ier)) ier = sts
  END SUBROUTINE cdfw_1c8

  SUBROUTINE cdfw_1c(ncid,varnam,varval,ier)
    !     Write 1 dimensional character data array
    !     Use cdfDefVar to define the Variable
    !
    implicit none
    ! Input
    integer,       intent(in) :: ncid
    character*(*), intent(in) :: varnam
    character*(*), intent(in) :: varval
    ! Output
    integer, optional,      intent(out) :: ier
    ! Local
    integer, dimension(1) :: dimlens
    integer               :: varid, ndims, sts
    if (PRESENT (ier)) ier = 1
    sts = nf_enddef(ncid)
    call cdfInqV(ncid,varnam,varid,dimlens,ndims,sts)
    if (sts .ne. 0) return
    if (ndims .ne. 1) then
       print "('% cdfPutVar_1c: --E-- The variable ',a,               &
            &         ' was defined as',i2,' dimensional')",varnam,ndims
       return
    end if
    sts = nf_put_var_text(ncid,varid,varval)
    call handle_err(sts,varnam,'cdfPutVar_1c','nf_put_var_text')
    if (PRESENT (ier)) ier = sts
  END SUBROUTINE cdfw_1c

  SUBROUTINE cdfw_0i(ncid,varnam,varval,ier)
    !     write integer scalar
    implicit none
    ! Input
    integer,       intent(in) :: ncid
    character*(*), intent(in) :: varnam
    integer,       intent(in) :: varval
    ! Output
    integer, optional,      intent(out) :: ier
    ! Local
    integer, dimension(1) :: dimlens
    integer               :: varid, ndims, sts
    if (PRESENT (ier)) ier = 1
    sts = nf_enddef(ncid)
    call cdfInqV(ncid,varnam,varid,dimlens,ndims,sts)
    if (sts .ne. 0) return
    if (ndims .ne. 0) then
       print "('% cdfPutVar_0i: --E-- The variable ',a,               &
            &         ' was defined as',i2,' dimensional')",varnam,ndims
       return
    end if
    sts = nf_put_var_int(ncid,varid,varval)
    call handle_err(sts,varnam,'cdfPutVar_0i','nf_put_var_int')
    if (PRESENT (ier)) ier = sts
  END SUBROUTINE cdfw_0i

  SUBROUTINE cdfw_0l(ncid,varnam,varval,ier)
    !     Write logical scalar
    implicit none
    ! Input
    integer,            intent(in) :: ncid
    character*(*),      intent(in) :: varnam
    logical,            intent(in) :: varval
    ! Output
    integer, optional,  intent(out) :: ier
    ! Local
    integer                         :: varval_i
    character*11, parameter :: logical_name = '__logical__'

    IF (varval) THEN
       varval_i = 1
    ELSE
       varval_i = 0
    END IF

    call cdfw_0i(ncid,trim(varnam)//logical_name,varval_i,ier)

  END SUBROUTINE cdfw_0l

  SUBROUTINE cdfw_0d(ncid,varnam,varval,ier)
    !     Write Real*8 Scalar
    !
    implicit none
    ! Input
    integer,       intent(in) :: ncid
    character*(*), intent(in) :: varnam
    REAL(KIND=r8), intent(in) :: varval
    ! Output
    integer, optional,      intent(out) :: ier
    ! Local
    integer, dimension(1) :: dimlens
    integer               :: varid, ndims, sts
    if (PRESENT (ier)) ier = 1
    sts = nf_enddef(ncid)
    call cdfInqV(ncid,varnam,varid,dimlens,ndims,sts)
    if (sts .ne. 0) return
    if (ndims .ne. 0) then
       print "('% cdfPutVar_0d: --E-- The variable ',a,               &
            &         ' was defined as',i2,' dimensional')",varnam,ndims
       return
    end if
    sts = nf_put_var_double(ncid,varid,varval)
    call handle_err(sts,varnam,'cdfPutVar_0d','nf_put_var_double')
    if (PRESENT (ier)) ier = sts
  END SUBROUTINE cdfw_0d

  SUBROUTINE cdfw_0c16(ncid,varnam,varval,ier)
    !     Write Complex*16 Scalar
    !
    implicit none
    ! Input
    integer,       intent(in) :: ncid
    character*(*), intent(in) :: varnam
    COMPLEX(KIND=r8), intent(in) :: varval
    ! Output
    integer, optional,      intent(out) :: ier
    ! Local
    integer, dimension(1) :: dimlens
    integer               :: varid, ndims, sts

    if (PRESENT (ier)) ier = 1
    sts = nf_enddef(ncid)
    call cdfInqV(ncid,trim(varnam)//cmplx_name,varid,dimlens,ndims,sts)
    if (sts .ne. 0) return
    if (ndims .ne. 1) then ! scalar complex stored as 2 element real array
       print "('% cdfPutVar_0c16: --E-- The variable ',a,               &
            &         ' was defined as',i2,' dimensional')",varnam,ndims
       return
    end if
    sts = nf_put_var_double(ncid,varid,varval)
    call handle_err(sts,varnam,'cdfPutVar_0c16','nf_put_var_double')
    if (PRESENT (ier)) ier = sts
  END SUBROUTINE cdfw_0c16

  SUBROUTINE cdfw_0f(ncid,varnam,varval,ier)
    !     Write default real
    !
    implicit none
    ! Input
    integer,       intent(in) :: ncid
    character*(*), intent(in) :: varnam
    REAL(KIND=r4), intent(in) :: varval
    ! Output
    integer, optional,      intent(out) :: ier
    ! Local
    integer, dimension(1) :: dimlens
    integer               :: varid, ndims, sts
    if (PRESENT (ier)) ier = 1
    sts = nf_enddef(ncid)
    call cdfInqV(ncid,varnam,varid,dimlens,ndims,sts)
    if (sts .ne. 0) return
    if (ndims .ne. 0) then
       print "('% cdfPutVar_0f: --E-- The variable ',a,               &
            &         ' was defined as',i2,' dimensional')",varnam,ndims
       return
    end if
    sts = nf_put_var_real(ncid,varid,varval)
    call handle_err(sts,varnam,'cdfPutVar_0f','nf_put_var_real')
    if (PRESENT (ier)) ier = sts
  END SUBROUTINE cdfw_0f

  SUBROUTINE cdfw_0c8(ncid,varnam,varval,ier)
    !     Write complex*8
    !
    implicit none
    ! Input
    integer,       intent(in) :: ncid
    character*(*), intent(in) :: varnam
    complex(KIND=r4), intent(in) :: varval
    ! Output
    integer, optional,      intent(out) :: ier
    ! Local
    integer, dimension(1) :: dimlens
    integer               :: varid, ndims, sts

    if (PRESENT (ier)) ier = 1
    sts = nf_enddef(ncid)
    call cdfInqV(ncid,trim(varnam)//cmplx_name,varid,dimlens,ndims,sts)
    if (sts .ne. 0) return
    if (ndims .ne. 1) then
       print "('% cdfPutVar_0c8: --E-- The variable ',a,               &
            &         ' was defined as',i2,' dimensional')",varnam,ndims
       return
    end if
    sts = nf_put_var_real(ncid,varid,varval)
    call handle_err(sts,varnam,'cdfPutVar_0c8','nf_put_var_real')
    if (PRESENT (ier)) ier = sts
  END SUBROUTINE cdfw_0c8

!------------------------------------------------------------------
! cdfDefVar implementation routines

  subroutine cdfDefVar(ncid,varnam,dimlens,xtype,ier,dimname)
    ! Define a Variable and its dimensions
    ! 03/08/99 C. Ludescher
    ! C. Ludescher/A. Pletzer Tue Apr  4 10:11:33 EDT 2000
    ! added support for complex types (ap) Wed May 16 15:06:34 EDT 2001
    ! added support for logical type Oct 1 2002 (SPH)

    implicit none
!!$  include "netcdf.inc"
    ! Input
    integer,              intent(in)  :: ncid
    character*(*),        intent(in)  :: varnam
    integer, dimension(:),intent(in)  :: dimlens
    character*(*),        intent(in)  :: xtype
    character*(*), optional, dimension(:), intent(in) :: dimname
    ! Output
    integer, optional,      intent(out) :: ier
    ! Local
    integer                 :: ndims, varid, vartype, maxdims
    integer                 :: i,     n,     status
    integer, dimension(3)   :: dimids
    character*(nf_max_name) :: dimnams(3)
    character*(5)           :: zdim
    character*(9)           :: zdim1
    character*(*), parameter :: cmplx_name = '__CmPlx_Re_Im'
    character*(*), parameter :: logical_name = '__logical__'
    logical :: is_complex, is_logical

    integer dims(3)
    integer flag

    is_complex=.FALSE.
    is_logical=.FALSE.

    if (PRESENT (ier)) ier = 1
    status = 0

    dims = dimlens
    if (xtype .eq. 'R8' .or. xtype .eq. 'r8') then
       vartype=nf_double
    else if (xtype .eq. 'C16' .or. xtype .eq. 'c16') then
       ! complex array = real array with double leading length
       vartype=nf_double
       is_complex = .TRUE.
       dims(1) = 2*abs(dimlens(1))
    else if (xtype .eq. 'INT' .or. xtype .eq. 'int') then
       vartype=nf_int
    else if (xtype .eq. 'LOG' .or. xtype .eq. 'log' .or. xtype .eq. 'BOOL' .or. xtype .eq.'bool') then
       vartype=nf_int !if nf_byte, MUST use integer*1
       is_logical = .TRUE.
    else if (xtype .eq. 'R4' .or. xtype .eq. 'r4') then
       vartype=nf_float
    else if (xtype .eq. 'C8' .or. xtype .eq. 'c8') then
       ! complex array = real array with double leading length
       vartype=nf_float
       is_complex = .TRUE.
       dims(1) = 2*abs(dimlens(1))
    else if (xtype .eq. 'CHAR' .or. xtype .eq. 'char') then
       vartype=nf_char
    end if

    ! The normal behavior is to decrement the rank if all
    ! sizes to the right are 1. However, in some cases it
    ! may be desirable to keep this dimension, in which the
    ! user needs to pass -1 as dimension.

    ndims = size(dims)
    n=ndims
    flag = 1
    do i=n,1,-1
       if ((dims(i) == 1 .or. dims(i) == 0) .and. flag==1 ) then
          ndims = i-1
       else
          flag=0
       endif
       dims(i) = abs(dims(i))
    end do

    ! Check ndims <= maxdims
    if (xtype .eq. 'CHAR' .or. xtype .eq. 'char') then
       maxdims = 2
    else
       maxdims = 3
    endif
    if (ndims .gt. maxdims) then
       WRITE(*,10)ndims,xtype
10     format('% cdfDefVar --E--  Rank',i3,' not supported for ',a)
       return
    endif

    do i=1,ndims
       if (dims(i).gt. 999999999) then
          print *,'% cdfDefVar --E-- dimension >= 1.e9 not supported'
          return
       else if (dims(i) .lt. 100000) then
          write(zdim,'(i5.5)') dims(i)
          dimnams(i) = 'dim_'//zdim
       else
          write(zdim1,'(i9.9)') dims(i)
          dimnams(i) = 'dim_'//zdim1
       endif
       if (PRESENT(dimname)) then 
          if (SIZE(dimname).ge.i) then
             if (LEN_TRIM(dimname(i)).gt.0) dimnams(i) = dimname(i)
          end if
       end if
       status = nf_inq_dimid(ncid,dimnams(i),dimids(i))
       if (status .ne. nf_noerr) then
          status = nf_def_dim(ncid,dimnams(i),dims(i),dimids(i))
          call handle_err(status,dimnams(i),'cdfdef','nf_def_dim')
       endif
    end do

    if (status .ne. 0) return

  ! replace any embedded blanks or non-alphanumeric characters with underscores
    varnam_noalpha = varnam
    CALL alpha_numeric(varnam_noalpha)
    ! do some name mangling to keep track of complex and 
    ! logical types
    if(is_complex) varnam_noalpha = trim(varnam_noalpha)//cmplx_name
    if(is_logical) varnam_noalpha = trim(varnam_noalpha)//logical_name

    status = nf_def_var(ncid,varnam_noalpha,vartype,ndims,dimids,varid)

    call handle_err(status,varnam,'cdfDefVar','nf_def_var')
    if (PRESENT (ier)) ier = status
  end subroutine cdfDefVar

  SUBROUTINE cdfd_3i(ncid,varnam,varval,ier,dimname)
    !     define 3 dimensional Integer data array
    implicit none
!!$  include "netcdf.inc"
    ! Input
    integer,                   intent(in) :: ncid
    character*(*),             intent(in) :: varnam
    integer, dimension(:,:,:), intent(in) :: varval
    character*(*), optional, dimension(:), intent(in) :: dimname
    ! Output
    integer, optional,        intent(out) :: ier
    ! Local
    character*3             :: vartype
    integer                 :: dims(3)

    vartype='INT'
    dims = -shape(varval)
    call cdfDefVar(ncid,varnam,dims,vartype,ier,dimname)

  END SUBROUTINE cdfd_3i

  SUBROUTINE cdfd_2i(ncid,varnam,varval,ier,dimname)
    !     define 2 dimensional Integer data array
    implicit none
!!$  include "netcdf.inc"
    ! Input
    integer,                 intent(in) :: ncid
    character*(*),           intent(in) :: varnam
    integer, dimension(:,:), intent(in) :: varval
    character*(*), optional, dimension(:), intent(in) :: dimname
    ! Output
    integer, optional,       intent(out) :: ier
    ! Local
    character*3             :: vartype
    integer                 :: dims(3)

    vartype='INT'
    dims = (/ -shape(varval), 1 /)
    call cdfDefVar(ncid,varnam,dims,vartype,ier,dimname)

  END SUBROUTINE cdfd_2i

  SUBROUTINE cdfd_1i(ncid,varnam,varval,ier,dimname)
    !     define 1 dimensional Integer data array
    implicit none
!!$  include "netcdf.inc"
    ! Input
    integer,               intent(in) :: ncid
    character*(*),         intent(in) :: varnam
    integer, dimension(:), intent(in) :: varval
    character*(*), optional, dimension(:), intent(in) :: dimname
    ! Output
    integer, optional,     intent(out) :: ier
    ! Local
    character*3             :: vartype
    integer                 :: dims(3)

    vartype='INT'
    dims = (/ -shape(varval), 1, 1 /)
    call cdfDefVar(ncid,varnam,dims,vartype,ier,dimname)

  END SUBROUTINE cdfd_1i

  SUBROUTINE cdfd_0i(ncid,varnam,varval,ier,dimname)
    !     define scalar Integer
    implicit none
!!$  include "netcdf.inc"
    ! Input
    integer,            intent(in) :: ncid
    character*(*),      intent(in) :: varnam
    integer,            intent(in) :: varval
    character*(*), optional, dimension(:), intent(in) :: dimname
    ! Output
    integer, optional,  intent(out) :: ier
    ! Local
    character*3             :: vartype
    integer                 :: dims(3)

    vartype='INT'
    dims = 0
    call cdfDefVar(ncid,varnam,dims,vartype,ier,dimname)

  END SUBROUTINE cdfd_0i

  SUBROUTINE cdfd_3l(ncid,varnam,varval,ier,dimname)
    !     define 3 dimensional logical array
    implicit none
!!$  include "netcdf.inc"
    ! Input
    integer,                   intent(in) :: ncid
    character*(*),             intent(in) :: varnam
    logical, dimension(:,:,:), intent(in) :: varval
    character*(*), optional, dimension(:), intent(in) :: dimname
    ! Output
    integer, optional,        intent(out) :: ier
    ! Local
    character*4             :: vartype
    integer                 :: dims(3)

    vartype='LOG'
    dims = -shape(varval)
    call cdfDefVar(ncid,varnam,dims,vartype,ier,dimname)

  END SUBROUTINE cdfd_3l

  SUBROUTINE cdfd_2l(ncid,varnam,varval,ier,dimname)
    !     define 2 dimensional logical array
    implicit none
!!$  include "netcdf.inc"
    ! Input
    integer,                 intent(in) :: ncid
    character*(*),           intent(in) :: varnam
    logical, dimension(:,:), intent(in) :: varval
    character*(*), optional, dimension(:), intent(in) :: dimname
    ! Output
    integer, optional,       intent(out) :: ier
    ! Local
    character*3            :: vartype
    integer                 :: dims(3)

    vartype='LOG'
    dims = (/ -shape(varval), 1 /)
    call cdfDefVar(ncid,varnam,dims,vartype,ier,dimname)

  END SUBROUTINE cdfd_2l

  SUBROUTINE cdfd_1l(ncid,varnam,varval,ier,dimname)
    !     define 1 dimensional logical array
    implicit none
!!$  include "netcdf.inc"
    ! Input
    integer,               intent(in) :: ncid
    character*(*),         intent(in) :: varnam
    logical, dimension(:), intent(in) :: varval
    character*(*), optional, dimension(:), intent(in) :: dimname
    ! Output
    integer, optional,     intent(out) :: ier
    ! Local
    character*3            :: vartype
    integer                 :: dims(3)

    vartype='LOG'
    dims = (/ -shape(varval), 1, 1 /)
    call cdfDefVar(ncid,varnam,dims,vartype,ier,dimname)

  END SUBROUTINE cdfd_1l

  SUBROUTINE cdfd_0l(ncid,varnam,varval,ier,dimname)
    !     define scalar logical
    implicit none
!!$  include "netcdf.inc"
    ! Input
    integer,            intent(in) :: ncid
    character*(*),      intent(in) :: varnam
    logical,            intent(in) :: varval
    character*(*), optional, dimension(:), intent(in) :: dimname
    ! Output
    integer, optional,  intent(out) :: ier
    ! Local
    character*3            :: vartype
    integer                 :: dims(3)

    vartype='LOG'
    dims = 0
    call cdfDefVar(ncid,varnam,dims,vartype,ier,dimname)

  END SUBROUTINE cdfd_0l

  SUBROUTINE cdfd_3d(ncid,varnam,varval,ier,dimname)
    !     define 3 dimensional double data array
    implicit none
!!$  include "netcdf.inc"
    ! Input
    integer,                   intent(in) :: ncid
    character*(*),             intent(in) :: varnam
    REAL(kind=r8), dimension(:,:,:), intent(in) :: varval
    character*(*), optional, dimension(:), intent(in) :: dimname
    ! Output
    integer, optional,        intent(out) :: ier
    ! Local
    character*2             :: vartype
    integer                 :: dims(3)

    vartype='R8'
    dims = -shape(varval)
    call cdfDefVar(ncid,varnam,dims,vartype,ier,dimname)

  END SUBROUTINE cdfd_3d

  SUBROUTINE cdfd_2d(ncid,varnam,varval,ier,dimname)
    !     define 2 dimensional double data array
    implicit none
!!$  include "netcdf.inc"
    ! Input
    integer,                 intent(in) :: ncid
    character*(*),           intent(in) :: varnam
    REAL(kind=r8), dimension(:,:), intent(in) :: varval
    character*(*), optional, dimension(:), intent(in) :: dimname
    ! Output
    integer, optional,       intent(out) :: ier
    ! Local
    character*2             :: vartype
    integer                 :: dims(3)

    vartype='R8'
    dims = (/ -shape(varval), 1 /)
    call cdfDefVar(ncid,varnam,dims,vartype,ier,dimname)

  END SUBROUTINE cdfd_2d

  SUBROUTINE cdfd_1d(ncid,varnam,varval,ier,dimname)
    !     define 1 dimensional double data array
    implicit none
!!$  include "netcdf.inc"
    ! Input
    integer,               intent(in) :: ncid
    character*(*),         intent(in) :: varnam
    REAL(kind=r8), dimension(:), intent(in) :: varval
    character*(*), optional, dimension(:), intent(in) :: dimname
    ! Output
    integer, optional,     intent(out) :: ier
    ! Local
    character*2             :: vartype
    integer                 :: dims(3)

    vartype='R8'
    dims = (/ -shape(varval), 1, 1 /)
    call cdfDefVar(ncid,varnam,dims,vartype,ier,dimname)

  END SUBROUTINE cdfd_1d

  SUBROUTINE cdfd_0d(ncid,varnam,varval,ier,dimname)
    !     define scalar double
    implicit none
!!$  include "netcdf.inc"
    ! Input
    integer,            intent(in) :: ncid
    character*(*),      intent(in) :: varnam
    REAL(kind=r8),      intent(in) :: varval
    character*(*), optional, dimension(:), intent(in) :: dimname
    ! Output
    integer, optional,  intent(out) :: ier
    ! Local
    character*2             :: vartype
    integer                 :: dims(3)

    vartype='R8'
    dims = 0
    call cdfDefVar(ncid,varnam,dims,vartype,ier,dimname)

  END SUBROUTINE cdfd_0d

  SUBROUTINE cdfd_3f(ncid,varnam,varval,ier,dimname)
    !     define 3 dimensional real data array
    implicit none
!!$  include "netcdf.inc"
    ! Input
    integer,                   intent(in) :: ncid
    character*(*),             intent(in) :: varnam
    REAL(kind=r4), dimension(:,:,:), intent(in) :: varval
    character*(*), optional, dimension(:), intent(in) :: dimname
    ! Output
    integer, optional,        intent(out) :: ier
    ! Local
    character*2             :: vartype
    integer                 :: dims(3)

    vartype='R4'
    dims = -shape(varval)
    call cdfDefVar(ncid,varnam,dims,vartype,ier,dimname)

  END SUBROUTINE cdfd_3f

  SUBROUTINE cdfd_2f(ncid,varnam,varval,ier,dimname)
    !     define 2 dimensional real data array
    implicit none
!!$  include "netcdf.inc"
    ! Input
    integer,                 intent(in) :: ncid
    character*(*),           intent(in) :: varnam
    REAL(kind=r4), dimension(:,:), intent(in) :: varval
    character*(*), optional, dimension(:), intent(in) :: dimname
    ! Output
    integer, optional,       intent(out) :: ier
    ! Local
    character*2             :: vartype
    integer                 :: dims(3)

    vartype='R4'
    dims = (/ -shape(varval), 1 /)
    call cdfDefVar(ncid,varnam,dims,vartype,ier,dimname)

  END SUBROUTINE cdfd_2f

  SUBROUTINE cdfd_1f(ncid,varnam,varval,ier,dimname)
    !     define 1 dimensional real data array
    implicit none
!!$  include "netcdf.inc"
    ! Input
    integer,               intent(in) :: ncid
    character*(*),         intent(in) :: varnam
    REAL(kind=r4), dimension(:), intent(in) :: varval
    character*(*), optional, dimension(:), intent(in) :: dimname
    ! Output
    integer, optional,     intent(out) :: ier
    ! Local
    character*2             :: vartype
    integer                 :: dims(3)

    vartype='R4'
    dims = (/ -shape(varval), 1, 1 /)
    call cdfDefVar(ncid,varnam,dims,vartype,ier,dimname)

  END SUBROUTINE cdfd_1f

  SUBROUTINE cdfd_0f(ncid,varnam,varval,ier,dimname)
    !     define scalar real
    implicit none
!!$  include "netcdf.inc"
    ! Input
    integer,            intent(in) :: ncid
    character*(*),      intent(in) :: varnam
    REAL(kind=r4),      intent(in) :: varval
    character*(*), optional, dimension(:), intent(in) :: dimname
    ! Output
    integer, optional,  intent(out) :: ier
    ! Local
    character*2             :: vartype
    integer                 :: dims(3)

    vartype='R4'
    dims = 0
    call cdfDefVar(ncid,varnam,dims,vartype,ier,dimname)

  END SUBROUTINE cdfd_0f

  SUBROUTINE cdfd_3c16(ncid,varnam,varval,ier,dimname)
    !     define 3 dimensional 128-bit complex data array
    implicit none
!!$  include "netcdf.inc"
    ! Input
    integer,                   intent(in) :: ncid
    character*(*),             intent(in) :: varnam
    COMPLEX(kind=r8), dimension(:,:,:), intent(in) :: varval
    character*(*), optional, dimension(:), intent(in) :: dimname
    ! Output
    integer, optional,        intent(out) :: ier
    ! Local
    character*3             :: vartype
    integer                 :: dims(3)

    vartype='C16'
    dims = -shape(varval)
    call cdfDefVar(ncid,varnam,dims,vartype,ier,dimname)

  END SUBROUTINE cdfd_3c16

  SUBROUTINE cdfd_2c16(ncid,varnam,varval,ier,dimname)
    !     define 2 dimensional 128-bit complex data array
    implicit none
!!$  include "netcdf.inc"
    ! Input
    integer,                 intent(in) :: ncid
    character*(*),           intent(in) :: varnam
    COMPLEX(kind=r8), dimension(:,:), intent(in) :: varval
    character*(*), optional, dimension(:), intent(in) :: dimname
    ! Output
    integer, optional,       intent(out) :: ier
    ! Local
    character*3             :: vartype
    integer                 :: dims(3)

    vartype='C16'
    dims = (/ -shape(varval), 1 /)
    call cdfDefVar(ncid,varnam,dims,vartype,ier,dimname)

  END SUBROUTINE cdfd_2c16

  SUBROUTINE cdfd_1c16(ncid,varnam,varval,ier,dimname)
    !     define 1 dimensional 128-bit complex data array
    implicit none
!!$  include "netcdf.inc"
    ! Input
    integer,               intent(in) :: ncid
    character*(*),         intent(in) :: varnam
    COMPLEX(kind=r8), dimension(:), intent(in) :: varval
    character*(*), optional, dimension(:), intent(in) :: dimname
    ! Output
    integer, optional,     intent(out) :: ier
    ! Local
    character*3             :: vartype
    integer                 :: dims(3)

    vartype='C16'
    dims = (/ -shape(varval), 1, 1 /)
    call cdfDefVar(ncid,varnam,dims,vartype,ier,dimname)

  END SUBROUTINE cdfd_1c16

  SUBROUTINE cdfd_0c16(ncid,varnam,varval,ier,dimname)
    !     define scalar 128-bit complex
    implicit none
!!$  include "netcdf.inc"
    ! Input
    integer,            intent(in) :: ncid
    character*(*),      intent(in) :: varnam
    COMPLEX(kind=r8),      intent(in) :: varval
    character*(*), optional, dimension(:), intent(in) :: dimname
    ! Output
    integer, optional,  intent(out) :: ier
    ! Local
    character*3             :: vartype
    integer                 :: dims(3)

    vartype='C16'
    dims = 0
    call cdfDefVar(ncid,varnam,dims,vartype,ier,dimname)

  END SUBROUTINE cdfd_0c16

  SUBROUTINE cdfd_3c8(ncid,varnam,varval,ier,dimname)
    !     define 3 dimensional 64-bit complex data array
    implicit none
!!$  include "netcdf.inc"
    ! Input
    integer,                   intent(in) :: ncid
    character*(*),             intent(in) :: varnam
    COMPLEX(kind=r4), dimension(:,:,:), intent(in) :: varval
    character*(*), optional, dimension(:), intent(in) :: dimname
    ! Output
    integer, optional,        intent(out) :: ier
    ! Local
    character*2             :: vartype
    integer                 :: dims(3)

    vartype='C8'
    dims = -shape(varval)
    call cdfDefVar(ncid,varnam,dims,vartype,ier,dimname)

  END SUBROUTINE cdfd_3c8

  SUBROUTINE cdfd_2c8(ncid,varnam,varval,ier,dimname)
    !     define 2 dimensional 64-bit complex data array
    implicit none
!!$  include "netcdf.inc"
    ! Input
    integer,                 intent(in) :: ncid
    character*(*),           intent(in) :: varnam
    COMPLEX(kind=r4), dimension(:,:), intent(in) :: varval
    character*(*), optional, dimension(:), intent(in) :: dimname
    ! Output
    integer, optional,       intent(out) :: ier
    ! Local
    character*2             :: vartype
    integer                 :: dims(3)

    vartype='C8'
    dims = (/ -shape(varval), 1 /)
    call cdfDefVar(ncid,varnam,dims,vartype,ier,dimname)

  END SUBROUTINE cdfd_2c8

  SUBROUTINE cdfd_1c8(ncid,varnam,varval,ier,dimname)
    !     define 1 dimensional 64-bit complex data array
    implicit none
!!$  include "netcdf.inc"
    ! Input
    integer,               intent(in) :: ncid
    character*(*),         intent(in) :: varnam
    COMPLEX(kind=r4), dimension(:), intent(in) :: varval
    character*(*), optional, dimension(:), intent(in) :: dimname
    ! Output
    integer, optional,     intent(out) :: ier
    ! Local
    character*2             :: vartype
    integer                 :: dims(3)

    vartype='C8'
    dims = (/ -shape(varval), 1, 1 /)
    call cdfDefVar(ncid,varnam,dims,vartype,ier,dimname)

  END SUBROUTINE cdfd_1c8

  SUBROUTINE cdfd_0c8(ncid,varnam,varval,ier,dimname)
    !     define scalar 64-bit complex
    implicit none
!!$  include "netcdf.inc"
    ! Input
    integer,            intent(in) :: ncid
    character*(*),      intent(in) :: varnam
    COMPLEX(kind=r4),      intent(in) :: varval
    character*(*), optional, dimension(:), intent(in) :: dimname
    ! Output
    integer, optional,  intent(out) :: ier
    ! Local
    character*2             :: vartype
    integer                 :: dims(3)

    vartype='C8'
    dims = 0
    call cdfDefVar(ncid,varnam,dims,vartype,ier,dimname)

  END SUBROUTINE cdfd_0c8

  SUBROUTINE cdfd_2c(ncid,varnam,varval,ier,dimname)
    !     define 2 dimensional character array
    implicit none
!!$  include "netcdf.inc"
    ! Input
    integer,                     intent(in) :: ncid
    character*(*),               intent(in) :: varnam
    character*(*), dimension(:), intent(in) :: varval
    character*(*), optional, dimension(:), intent(in) :: dimname
    ! Output
    integer, optional,       intent(out) :: ier
    ! Local
    character*4             :: vartype
    integer                 :: dims(3)

    vartype='CHAR'
    dims(1) = len(varval(1))
    dims(2) = size(varval);   dims(3) = 1
    call cdfDefVar(ncid,varnam,dims,vartype,ier,dimname)

  END SUBROUTINE cdfd_2c

  SUBROUTINE cdfd_1c(ncid,varnam,varval,ier,dimname)
    !     define 1 dimensional character array
    implicit none
!!$  include "netcdf.inc"
    ! Input
    integer,         intent(in) :: ncid
    character*(*),   intent(in) :: varnam
    character*(*),   intent(in) :: varval
    character*(*), optional, dimension(:), intent(in) :: dimname
    ! Output
    integer, optional,     intent(out) :: ier
    ! Local
    character*4             :: vartype
    integer                 :: dims(3)

    vartype='CHAR'
    dims(1) = len(varval);  dims(2:3) = 1
    if (dims(1) .le. 1) dims(1) = -1
    call cdfDefVar(ncid,varnam,dims,vartype,ier,dimname)

  END SUBROUTINE cdfd_1c
!DEC$ ENDIF
END MODULE ezcdf_GenPut
