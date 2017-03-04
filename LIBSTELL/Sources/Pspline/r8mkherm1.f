      subroutine r8mkherm1(fun,x,nx,fherm)
C
C  create a data set for Hermite interpolation, from evaluation of
C  function and derivatives.  function of 1 indep. coordinate.
C
C  input:
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER nx,ix
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 fun,dfdx
!============
      external fun                      ! passed real function(x,dfdx)
      REAL*8 x(nx)                        ! x coordinate array
C
C  the passed function fun must have the interface:
C
C        real function <name>(x,dfdx)
C  where x is input, the function returns the function value and the
C  argument dfdx returns as output the function derivative at x.
C
C  output:
C
      REAL*8 fherm(0:1,nx)                ! function data & derivative
C
C  fherm(0,j) = function value f at x(j)
C  fherm(1,j) = derivative df/dx at x(j)
C
C----------------------------
C
      do ix=1,nx
         fherm(0,ix)=fun(x(ix),dfdx)
         fherm(1,ix)=dfdx
      enddo
C
      return
      end
