      subroutine gridintrp2d(
     >   x_newgrid,nx_new,y_newgrid,ny_new,f_new,if1,
     >   nx,xpkg,ny,ypkg,jspline,fspl,icoeff,ixdim,iydim,iwarn,ier)
c
c  regrid a spline function f defined vs. x,y as in xpkg,ypkg
c  to a new grid, given by x_newgrid, y_newgrid
c
c  set warning flag if the range x_newgrid, y_newgrid
c  exceeds the range of the original xpkg,ypkg.
c
c  (xpkg,ypkg -- axis data, see genxpkg subroutine)
c
c  input:
c
      real x_newgrid(nx_new)            ! new x grid
      real y_newgrid(ny_new)            ! new y grid
c
c  output:
c
      integer if1                       ! 1st dimension of f_new
      real f_new(if1,ny_new)            ! f evaluated on this grid
c
c  input:
c
      integer nx                        ! size of old grid
      real xpkg(nx,4)                   ! old grid "package"
      integer ny                        ! size of old grid
      real ypkg(ny,4)                   ! old grid "package"
c
      integer jspline(2)                ! interpolation type by dimension
      integer :: icoeff                 ! coefficients per data point
      integer :: ixdim,iydim            ! fspl dimensions 
      !  =nx,ny unless zonal step function interpolation is used

      real fspl(icoeff,ixdim,iydim)     ! spline coefficients of f
c
c  output:
c  condition codes, =0 for normal exit
c
      integer iwarn                     ! =1 if new grid points out of range
      integer ier                       ! =1 if there is an argument error
c
c--------------------------------------------
c  local
c
      real ytmp(nx_new)
      integer ict(6)
c
      data ict/1,0,0,0,0,0/
c
c--------------------------------------------
c
      do iy=1,ny_new
         ytmp=y_newgrid(iy)
         call vecintrp2d(ict,nx_new,x_newgrid,ytmp,nx_new,f_new(1,iy),
     >      nx,xpkg,ny,ypkg,jspline,fspl,icoeff,ixdim,iydim,iwarn,ier)
      enddo
c
      return
      end
