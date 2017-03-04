subroutine ezspline_cdfget3(ncid,zname,fspl,idim1,idim2,idim3,ifail)

  use EZcdf
  implicit NONE

  ! (due to ezspline rank limitation) read 4d object as 3d object
  ! use f77 interface style, should prevent unnecessary array copy.

  integer, intent(in) :: ncid             ! opened NetCDF file
  character*(*),intent(in) :: zname       ! name to us for writing
  integer, intent(in) :: idim1,idim2,idim3   ! spline data dimensions
  real*8, intent(out) :: fspl(idim1,idim2,idim3)  ! spline data & coefficients
  integer, intent(out) :: ifail           ! status retruned from NetCDF

  call cdfGetVar(ncid, trim(zname), fspl, ifail)

end subroutine ezspline_cdfget3
