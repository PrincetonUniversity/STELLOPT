      subroutine r8gridherm3(
     >   x_newgrid,nx_new,y_newgrid,ny_new,z_newgrid,nz_new,
     >   f_new,if1,if2,
     >   nx,xpkg,ny,ypkg,nz,zpkg,fspl,inf2,inf3,iwarn,ier)
c
c  regrid a hermite function f defined vs. x,y,z as in xpkg, etc.
c  to a new grid, given by x_newgrid, y_newgrid, z_newgrid
c
c  set warning flag if the range exceeds the range of the
c  original x/y/zpkg's.
c
c  (xpkg/ypkg/zpkg -- axis data, see genxpkg subroutine)
c
c  input:
c
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ny_new,nz_new,nx_new,iz,iy
!============
      REAL*8 x_newgrid(nx_new)            ! new x grid
      REAL*8 y_newgrid(ny_new)            ! new y grid
      REAL*8 z_newgrid(nz_new)            ! new z grid
c
c  output:
c
      integer if1,if2                   ! 1st dimensions of f_new
      REAL*8 f_new(if1,if2,nz_new)        ! f evaluated on this grid
c
c  input:
c
      integer nx                        ! size of old grid
      REAL*8 xpkg(nx,4)                   ! old grid "package"
      integer ny                        ! size of old grid
      REAL*8 ypkg(ny,4)                   ! old grid "package"
      integer nz                        ! size of old grid
      REAL*8 zpkg(nz,4)                   ! old grid "package"
c
      integer inf2,inf3                 ! array dimensions
      REAL*8 fspl(0:7,inf2,inf3,nz)       ! hermite coefficients of f
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
      REAL*8 ytmp(nx_new)
      REAL*8 ztmp(nx_new)
      integer ict(10)
c
      data ict/1,0,0,0,0,0,0,0,0,0/
c
c--------------------------------------------
c
      do iz=1,nz_new
         ztmp=z_newgrid(iz)
         do iy=1,ny_new
            ytmp=y_newgrid(iy)
            call r8vecherm3(ict,nx_new,x_newgrid,ytmp,ztmp,
     >         nx_new,f_new(1,iy,iz),
     >         nx,xpkg,ny,ypkg,nz,zpkg,fspl,inf2,inf3,iwarn,ier)
         enddo
      enddo
c
      return
      end
