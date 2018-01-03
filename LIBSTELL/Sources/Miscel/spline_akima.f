      SUBROUTINE spline_akima(x,y,xx,yy,npts,iflag)
! 2012-07-27: Typo correction, J Geiger
      USE stel_kinds
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: npts        ! number of active points
      integer, intent(inout) :: iflag
      real(rprec), intent(in)  :: x
      real(rprec), intent(out) :: y
      real(rprec), dimension(npts), intent(in) :: xx, yy
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(rprec), dimension(-1:size(xx)+2) :: a, b, c, d
      real(rprec), dimension(-1:size(xx)+2) :: xloc, yloc
      real(rprec), dimension(-1:size(xx)+2) :: m, t, dm, p, q
      integer :: i, ix, iv, ivm1
      real(rprec) :: cl,bl,cr,br,dx
!-----------------------------------------------

!     print * , x, y
!     print * , npts
!     print * , xx
!     print * , yy
!     stop 'Test!'

      iflag = 0  !initialization
      ix = size(xx)
!     print *, size(xx)
!     print *, size(yy)
      if(npts > ix)
     >   stop'spline_akima: more active points requested than available'
      if(ix /= size(yy)) stop 'size mismatch of xx and yy!'
      iv=npts
      ivm1  = iv-1
! initialize local variables
      a = 0._DP ; b = 0._DP ; c = 0._DP ; d = 0._DP
      xloc = 0._DP ; yloc = 0._DP
      m = 0._DP ; t = 0._DP ; dm = 0._DP
      p = 0._DP ; q = 0._DP

      xloc(1:iv)=xx
      xloc(-1)= 2*xloc(1)-xloc(3)
      xloc( 0)= xloc(1)+xloc(2)-xloc(3)
      xloc(iv+2)= 2*xloc(iv)-xloc(iv-2)
      xloc(iv+1)= xloc(iv)+xloc(iv-1)-xloc(iv-2)
      yloc(1:iv)=yy
! calculate linear derivatives as far as existent
      m(1:iv-1) = (yloc(2:iv)-yloc(1:iv-1))/
     >            (xloc(2:iv)-xloc(1:iv-1))
! values for i=0, -1 and iv, iv+1 by quadratic extrapolation:
      cl = (m(2)-m(1))/(xloc(3)-xloc(1))
      bl = m(1) - cl*(xloc(2)-xloc(1))
      cr = (m(iv-2)-m(iv-1))/(xloc(iv)-xloc(iv-2))
      br = m(iv-2) - cr*(xloc(iv-1)-xloc(iv-2))
      yloc( 0)=yloc(1)+bl*(xloc( 0)-xloc(1))+
     >               cl*(xloc( 0)-xloc(1))**2
      yloc(-1)=yloc(1)+bl*(xloc(-1)-xloc(1))+
     >               cl*(xloc(-1)-xloc(1))**2
      yloc(iv+1)=yloc(iv)+br*(xloc(iv+1)-xloc(iv))+
     >               cr*(xloc(iv+1)-xloc(iv))**2
      yloc(iv+2)=yloc(iv)+br*(xloc(iv+2)-xloc(iv))+
     >               cr*(xloc(iv+2)-xloc(iv))**2
! rest of linear derivatives
      m(-1) = (yloc(0)-yloc(-1))/(xloc(0)-xloc(-1))
      m( 0) = (yloc(1)-yloc( 0))/(xloc(1)-xloc( 0))
      m(iv  ) = (yloc(iv+1)-yloc(iv  ))/(xloc(iv+1)-xloc(iv  ))
      m(iv+1) = (yloc(iv+2)-yloc(iv+1))/(xloc(iv+2)-xloc(iv+1))
! calculate weights for derivatives
      dm(-1:iv)= abs(m(0:iv+1)-m(-1:iv))
      where (dm /= 0._DP) !exclude division by zero
        p(1:iv) = dm(1:iv)/(dm(1:iv)+dm(-1:iv-2))
      end where
      where (dm /= 0._DP) !exclude division by zero
        q(1:iv) = dm(-1:iv-2)/(dm(1:iv)+dm(-1:iv-2))
      end where
      t(1:iv) = p(1:iv)*m(0:iv-1)+q(1:iv)*m(1:iv)
      where ( p(1:iv)+q(1:iv) < TINY(1._DP)) ! in case of two zeros give equal weight
        t(1:iv) = 0.5_DP*m(0:iv-1)+0.5_DP*m(1:iv)
      end where
! fix coefficients
      a = yloc
      b = t
      c(1:iv-1) = (3*m(1:iv-1)-t(2:iv)-2*t(1:iv-1))/
     >               (xloc(2:iv)-xloc(1:iv-1))
      d(1:iv-1) = (t(2:iv)+t(1:iv-1)-2*m(1:iv-1))/
     >               (xloc(2:iv)-xloc(1:iv-1))**2

! calculation
      if(x<xloc(1) .or. x>xloc(iv))then
        y=0.0
        iflag=-1
        return
      endif
      if(x == xloc(iv)) then
        y = yy(iv) 
        iflag = 0
        return
      endif
      do i=1,iv-1
        if(x >= xloc(i) .and. x < xloc(i+1))then
          dx=x-xloc(i)
          y=a(i)+dx*(b(i)+dx*(c(i)+d(i)*dx))
          iflag = 0
          return
        endif
      enddo

      END SUBROUTINE spline_akima

