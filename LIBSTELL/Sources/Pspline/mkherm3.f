      subroutine mkherm3(fun,x,nx,y,ny,z,nz,fherm)
C
C  create a data set for Hermite interpolation, from evaluation of
C  function and derivatives.  Function of 3 indep. coordinates.
C
C  input:
C
      external fun                      ! passed real function(x,y,z)
      real x(nx)                        ! x coordinate array  (1st dim)
      real y(ny)                        ! y coordinate array  (2nd dim)
      real z(nz)                        ! z coordinate array  (3rd dim)
C
C  the passed function fun must have the interface:
C
C        real function <name>(x,y,z,dfdx,dfdy,dfdz,
C               d2fdxdy,d2fdxdz,d2fdydz,d3fdxdydz)
C  where x,y,z are input, the function returns the function value,
C  and the arguments dfdx, dfdy and dfdz return as output the function
C  derivative at the point (x,y,z).
C
C  output:
C
      real fherm(0:7,nx,ny,nz)          ! function data & derivatives
C
C  fherm(0,i,j,k) = function value f at x(i),y(j),z(k)
C  fherm(1,i,j,k) = derivative df/dx at x(i),y(j),z(k)
C  fherm(2,i,j,k) = derivative df/dy at x(i),y(j),z(k)
C  fherm(3,i,j,k) = derivative df/dz at x(i),y(j),z(k)
C  fherm(4,i,j,k) = derivative d2f/dxdy at x(i),y(j),z(k)
C  fherm(5,i,j,k) = derivative d2f/dxdz at x(i),y(j),z(k)
C  fherm(6,i,j,k) = derivative d2f/dydz at x(i),y(j),z(k)
C  fherm(7,i,j,k) = derivative d3f/dxdydz at x(i),y(j),z(k)
C
C----------------------------
C
      do iz=1,nz
         do iy=1,ny
            do ix=1,nx
               fherm(0,ix,iy,iz)=
     >            fun(x(ix),y(iy),z(iz),dfdx,dfdy,dfdz,
     >                d2fdxdy,d2fdxdz,d2fdydz,d3fdxyz)
               fherm(1,ix,iy,iz)=dfdx
               fherm(2,ix,iy,iz)=dfdy
               fherm(3,ix,iy,iz)=dfdz
               fherm(4,ix,iy,iz)=d2fdxdy
               fherm(5,ix,iy,iz)=d2fdxdz
               fherm(6,ix,iy,iz)=d2fdydz
               fherm(7,ix,iy,iz)=d3fdxyz
            enddo
         enddo
      enddo
C
      return
      end
