      subroutine spline_cubic_int(x,y,xx,yy,n,iflag)
      USE stel_kinds
      implicit none
!----------------------------------------------------------------------
! iflag = 0   normal return
! iflag =-1   x-request outside of bounds
! iflag =-2   xx arrays with two equal entries or not in correct order
!----------------------------------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)  :: n
      integer, intent(inout)  :: iflag
      real(rprec), intent(in)  :: x
      real(rprec), dimension(n), intent(in)  :: xx,yy
      real(rprec), intent(out) :: y
      real(rprec), dimension(n)  :: y2, dxx
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(rprec) :: yp1, ypn
      real(rprec) :: c
!-----------------------------------------------

      iflag = 0  !initialization
      if(x < xx(1) .or. x > xx(n)) then
        iflag = -1
        y=0.d+0
        return
      endif
      dxx(1:n-1)=xx(2:n)-xx(1:n-1)
      if(minval(dxx(1:n-1)) <= 0.d+0) then
        iflag=-2
        return
      endif

! fix boundary derivatives by quadratic fit
! left part
      c=((yy(3)-yy(1))/(xx(3)-xx(1))-(yy(2)-yy(1))/(xx(2)-xx(1)))
     >  /(xx(3)-xx(2))
      yp1=(yy(2)-yy(1))/(xx(2)-xx(1))-c*(xx(2)-xx(1))
! right part
      c=((yy(n-2)-yy(n))/(xx(n-2)-xx(n)) 
     >  -(yy(n-1)-yy(n))/(xx(n-1)-xx(n)))
     >  /(xx(n-2)-xx(n-1))
      ypn=(yy(n-1)-yy(n))/(xx(n-1)-xx(n))-c*(xx(n-1)-xx(n))

      call spline_int(xx,yy,n,yp1,ypn,y2)
      call splint_int(xx,yy,y2,n,x,y)

      return

      contains

      subroutine spline_int(x,y,n,yp1,ypn,y2)
! taken from numerical recipes f77 and recoded.
! Given the arrays x(1:n) and y(1:n) containing the tabulated function
! with x(1) < x(2) <...< x(n) and given values yp1 and ypn for the first
! derivative of the interpolating function at points q and n, respectively,
! this routine returns an array y2(1:n) of length n which contains the
! second derivatives of the interpolating function at the tabulated points x(i).
! If yp1 and/or ypn are equatl to 1ed+30 or larger, the routine is signaled
! to set the correspoinding boundary condition for a natural spline with zero
! derivative on that boundary.
! nmax is the largest anticipated value of n.
      integer, intent(in) :: n
      integer, parameter :: nmax=500
      real(rprec), intent(in) :: yp1, ypn
      real(rprec), dimension(n), intent(in) :: x, y
      real(rprec), dimension(n), intent(out) :: y2
      integer :: i,k
      real(rprec) :: p, qn, sig, un
      real(rprec), dimension(n) :: u

      if(yp1 > .99d+30) then
        y2(1)=0.d+0
        u(1) =0.d+0
      else
        y2(1)=-0.5d+0
        u(1)=(3.d+0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif

      do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.d+0
        y2(i)=(sig-1.d+0)/p
        u(i)=(6.d+0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     >       /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo
      if (ypn > .99d+30)then
        qn=0.d+0
        un=0.d+0
      else
        qn=0.5d+0
        un=(3.d+0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d+0)
      do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      enddo
      return
      end subroutine spline_int

      subroutine splint_int(xa,ya,y2a,n,x,y)
! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function
! with the xa(i)'s in order), and given the array y2a(1:n), which is the
! output from spline above, and given a value of x, this routine returns
! a cubic-spline interpolated value y.
      implicit none
      integer, intent(in) :: n
      real(rprec), intent(in) :: x
      real(rprec), intent(out) :: y
      real(rprec), dimension(n), intent(in) :: xa, ya, y2a
!- local -------------------
      real(rprec) :: h, c, dx
      integer :: k, khi, klo
      real(rprec), dimension(n) :: dxa,dxa3, dya, dy2a, my2a, mya

      klo=1
      khi=n
      do
        if (khi-klo <= 1) exit  !inverted num.rec. condition for endless do-loop exit
        k=(khi+klo)/2
        if(xa(k) > x) then
          khi=k
        else
          klo=k
        endif
      enddo

      dxa(1:n-1)=xa(2:n)-xa(1:n-1)
      dxa3=dxa**3
      dya(1:n-1)=ya(2:n)-ya(1:n-1)
      mya(1:n-1)=ya(2:n)+ya(1:n-1)
      dy2a(1:n-1)=y2a(2:n)-y2a(1:n-1)
      my2a(1:n-1)=y2a(2:n)+y2a(1:n-1)
! integral up to specific interval [klo:khi] in which x is.
      c = .5d+0*dot_product(dxa(1:klo-1),mya(1:klo-1))
     >    -dot_product(dxa3(1:klo-1),my2a(1:klo-1))/24.0d+0
      h=xa(khi)-xa(klo)
      if(h ==0.d+0)
     >    stop "spline_cubic: bad xa input! xa(i) have to be distinct!"
      dx=x-xa(klo)
      y = ya(klo)*dx +.5d+0*dya(klo)*dx**2/h +
     >    y2a(klo)*dx**2*(2*dx-3*h)/12.d+0 +
     >    dy2a(klo)*dx**2*(dx**2-2*h**2)/(h*24.d+0)
      y= y+c
      return
      end subroutine splint_int

      end subroutine spline_cubic_int

