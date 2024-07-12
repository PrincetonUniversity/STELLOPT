      subroutine mkherm2(fun,x,nx,y,ny,fherm)
C
C  create a data set for Hermite interpolation, from evaluation of
C  function and derivatives.  function of 2 indep. coordinates.
C
C  input:
C
      external fun                      ! passed real function(x,y)
      real x(nx)                        ! x coordinate array
      real y(ny)                        ! y coordinate array
C
C  the passed function fun must have the interface:
C
C        real function <name>(x,y,dfdx,dfdy,d2fdxdy)
C  where x,y are input, the function returns the function value,
C  and the arguments dfdx and dfdy return as output the function
C  derivative at the point (x,y).
C
C  output:
C
      real fherm(0:3,nx,ny)             ! function data & derivatives
C
C  fherm(0,i,j) = function value f at x(i),y(j)
C  fherm(1,i,j) = derivative df/dx at x(i),y(j)
C  fherm(2,i,j) = derivative df/dy at x(i),y(j)
C  fherm(3,i,j) = derivative d2f/dxdy at x(i),y(j)
C
C----------------------------
C
      do ix=1,nx
         do iy=1,ny
            fherm(0,ix,iy)=fun(x(ix),y(iy),dfdx,dfdy,d2fdxdy)
            fherm(1,ix,iy)=dfdx
            fherm(2,ix,iy)=dfdy
            fherm(3,ix,iy)=d2fdxdy
         enddo
      enddo
C
      return
      end
