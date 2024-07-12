      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      USE stel_kinds, ONLY: dp
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER n
      REAL(dp) :: yp1,ypn,x(n),y(n),y2(n)
!Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., yi = f(xi), with
!x1 < x2 < .. . < xN, and given values yp1 and ypn for the first derivative of the interpolating
!function at points 1 and n, respectively, this routine returns an array y2(1:n) of
!length n which contains the second derivatives of the interpolating function at the tabulated
!points xi. If yp1 and/or ypn are equal to 1 × 10**30 or larger, the routine is signaled to set
!the corresponding boundary condition for a natural spline, with zero second derivative on
!that boundary. If yp1/ypn are -1 X 10**30 or smaller, set the bdy conditions so y2(1) = y2(2) or
!y2(n) = y2(n-1) (added by SPH)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER i,k
      REAL(dp) :: p,qn,sig,un,u(n)
!-----------------------------------------------
      IF (yp1.gt..99e30_dp) THEN  !The lower boundary condition is set either to be "natural"
         y2(1)=0
         u(1)=0
      ELSE IF (yp1.lt.-0.99e30_dp) THEN  !continuous 2nd derivative at boundary, y2(1) = y2(2)
         y2(1) = 1;   u(1) = 0
      ELSE                        !or else to have a specified first derivative.
         y2(1)=-0.5_dp
         u(1)=(3.0_dp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      ENDIF
      
      DO i=2,n-1               !This is the decomposition loop of the tridiagonal
                                  !algorithm. y2 and u are used for temporary
                                  !storage of the decomposed factors.
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2
         y2(i)=(sig-1._dp)/p
         u(i)=(6*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     1           /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      ENDDO 
      IF (ypn.gt..99e30_dp) THEN     !The upper boundary condition is set either to be "natural"
         qn=0
         un=0
      ELSE IF (ypn.lt.-0.99e30_dp) THEN  !continuous 2nd derivative at bdy
         qn=-1;  un=0
      ELSE                        !or else to have a specified first derivative.
         qn=0.5_dp
         un=(3._dp/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      ENDIF
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1)
      DO k=n-1,1,-1            !This is the backsubstitution loop of the tridiagonal algorithm. 
         y2(k)=y2(k)*y2(k+1)+u(k)
      ENDDO 

      END SUBROUTINE spline
