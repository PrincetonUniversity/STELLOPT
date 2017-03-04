      subroutine r8spgrid(x_newgrid,nx_new,f_new,nx,xpkg,fspl,iwarn,ier)
c
c  regrid a spline function f defined vs. x as in xpkg
c  to a new grid, given by x_newgrid.
c
c  set warning flag if the range x_newgrid exceeds the range of the
c  original xpkg.
c
c  (xpkg -- see genxpkg subroutine)
c
c  input:
c
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER nx_new
!============
      REAL*8 x_newgrid(nx_new)            ! new grid
c
c  output:
c
      REAL*8 f_new(nx_new)                ! f evaluated on this grid
c
c  input:
c
      integer nx                        ! size of old grid
      REAL*8 xpkg(nx,4)                   ! old grid "package"
      REAL*8 fspl(4,nx)                   ! spline coefficients of f
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
      integer ict(3)
c
      data ict/1,0,0/
c
c--------------------------------------------
c
      call r8spvec(ict,nx_new,x_newgrid,nx_new,f_new,nx,xpkg,fspl,
     >   iwarn,ier)
c
      return
      end
