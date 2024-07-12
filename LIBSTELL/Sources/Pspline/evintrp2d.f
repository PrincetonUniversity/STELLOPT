      subroutine evintrp2d(xget,yget,x,nx,y,ny,jspline,
     >                   f,icoeff,ixdim,iydim,ict,fval,ier)
C
      implicit NONE
C
C  evaluate a 2d Hybrid interpolant on a rectilinear grid
C
C  this subroutine calls three subroutines:
C     vecin2d_argchk -- error check
C     herm2xy  -- find cell containing (xget,yget)
C     fvintrp2d  -- evaluate interpolant function and (optionally) derivatives
C
C  input arguments:
C  ================
C
      real xget,yget                    ! target of this interpolation
      integer nx,ny                     ! grid sizes
      real x(nx)                        ! ordered x grid
      real y(ny)                        ! ordered y grid
C
      integer :: jspline(2)             ! interpolation method for each
C           dimension: jspline(1) for x; jspline(2) for y
C           =-1: zonal step function; =0: pc lin; =1: Hermite; =2: Spline
C
      integer :: icoeff                 ! no. of coefficients per data point

      integer :: ixdim,iydim            ! dimensioning:
      !  ixdim=nx-1 for zonal step in x; otherwise nx
      !  iydim=ny-1 for zonal step in y; otherwise ny

      real f(icoeff,ixdim,iydim)        ! function data w/coefficients
C
      integer ict(6)                    ! code specifying output desired
C
C  Note on derivatives: for dimensions along which zonal step function
C    interpolation is done, ANY derivative request returns ZERO.
C    For dimensions along which piecewise linear or Hermite interpolation
C    are done, more than one differentiation returns ZERO!
C
C  Derivative controls are the same as for the compact 2d spline evaluation
C  routine (evbicub):
C
C  ict(1)=1 -- return f  (0, don't)
C  ict(2)=1 -- return df/dx  (0, don't)
C  ict(3)=1 -- return df/dy  (0, don't)
C  ict(4)=1 -- return d2f/dx2  (0, don't)
C  ict(5)=1 -- return d2f/dy2  (0, don't)
C  ict(6)=1 -- return d2f/dxdy (0, don't)
c                   the number of non zero values ict(1:6)
c                   determines the number of outputs...
c
c  new dmc December 2005 -- access to higher derivatives (even if not
c  continuous-- but can only go up to 3rd derivatives on any one coordinate.
c     if ict(1)=3 -- want 3rd derivatives
c          ict(2)=1 for d3f/dx3
c          ict(3)=1 for d3f/dx2dy
c          ict(4)=1 for d3f/dxdy2
c          ict(5)=1 for d3f/dy3
c               number of non-zero values ict(2:5) gives no. of outputs
c     if ict(1)=4 -- want 4th derivatives
c          ict(2)=1 for d4f/dx3dy
c          ict(3)=1 for d4f/dx2dy2
c          ict(4)=1 for d4f/dxdy3
c               number of non-zero values ict(2:4) gives no. of outputs
c     if ict(1)=5 -- want 5th derivatives
c          ict(2)=1 for d5f/dx3dy2
c          ict(3)=1 for d5f/dx2dy3
c               number of non-zero values ict(2:3) gives no. of outputs
c     if ict(1)=6 -- want 6th derivatives
c          d6f/dx3dy3 -- one value is returned.
C
C output arguments:
C =================
C
      real fval(6)                      ! output data
      integer ier                       ! error code =0 ==> no error
C
C  fval(1) receives the first output (depends on ict(...) spec)
C  fval(2) receives the second output (depends on ict(...) spec)
C  fval(3) receives the third output (depends on ict(...) spec)
C  fval(4) receives the fourth output (depends on ict(...) spec)
C  fval(5) receives the fifth output (depends on ict(...) spec)
C  fval(6) receives the sixth output (depends on ict(...) spec)
C
C  examples:
C    on input ict = [1,1,1,0,0,1]
C   on output fval= [f,df/dx,df/dy,d2f/dxdy], elements 5 & 6 not referenced.
C
C    on input ict = [1,0,0,0,0,0]
C   on output fval= [f] ... elements 2 -- 6 never referenced.
C
C    on input ict = [0,0,0,1,1,0]
C   on output fval= [d2f/dx2,d2f/dy2] ... elements 3 -- 6 never referenced.
C
C    on input ict = [0,0,1,0,0,0]
C   on output fval= [df/dy] ... elements 2 -- 6 never referenced.
C
C  ier -- completion code:  0 means OK
C-------------------
C  local:
C
      integer i,j                       ! cell indices
C
C  normalized displacement from (x(i),y(j)) corner of cell.
C    xparam=0 @x(i)  xparam=1 @x(i+1)
C    yparam=0 @y(j)  yparam=1 @y(j+1)
C
      real xparam,yparam
C
C  cell dimensions and
C  inverse cell dimensions hxi = 1/(x(i+1)-x(i)), hyi = 1/(y(j+1)-y(j))
C
      real hx,hy
      real hxi,hyi
C
C  0 .le. xparam .le. 1
C  0 .le. yparam .le. 1
C
C  ** the interface is very similar to herm2ev.for; can use herm2xy **
C---------------------------------------------------------------------
C
      call vecin2d_argchk('evintrp2d',jspline,
     >     icoeff,nx,ny,ixdim,iydim,ier)
      if(ier.ne.0) return
c
      call herm2xy(xget,yget,x,nx,y,ny,0,0,
     >     i,j,xparam,yparam,hx,hxi,hy,hyi,ier)
      if(ier.ne.0) return
c
      call fvintrp2d(ict,1,1,
     >     fval,i,j,xparam,yparam,hx,hxi,hy,hyi,
     >     jspline,f,icoeff,ixdim,iydim)
C
      return
      end
C---------------------------------------------------------------------
C  evaluate Hybrid function interpolation -- 2d fcn
C   NO ERROR CHECKING
C
      subroutine fvintrp2d(ict,ivec,ivecd,
     >     fval,ii,jj,xparam,yparam,hx,hxi,hy,hyi,
     >     jspline,fin,icoeff,ixdim,iydim)
C
      implicit NONE
C
      integer ict(6)                    ! requested output control
      integer ivec                      ! vector length
      integer ivecd                     ! vector dimension (1st dim of fval)
C
      integer ii(ivec),jj(ivec)         ! target cells (i,j)
      real xparam(ivec),yparam(ivec)
                          ! normalized displacements from (i,j) corners
C
      real hx(ivec),hy(ivec)            ! grid spacing, and
      real hxi(ivec),hyi(ivec)          ! inverse grid spacing 1/(x(i+1)-x(i))
                                        ! & 1/(y(j+1)-y(j))
C
      integer :: jspline(2)             ! interpolation type control
      integer :: icoeff                 ! coefficient dimension
      integer :: ixdim,iydim            ! x & y dimensions
      real fin(icoeff,ixdim,iydim)      ! data on which to interpolate...
C
      real fval(ivecd,6)                ! output returned
C
C  for detailed description of fin, ict and fval see subroutine
C  evintrp2d comments.  Note ict is not vectorized; the same output
C  is expected to be returned for all input vector data points.
C
C  note that the index inputs ii,jj and parameter inputs
C     xparam,yparam,hx,hxi,hy,hyi are vectorized, and the
C     output array fval has a vector ** 1st dimension ** whose
C     size must be given as a separate argument
C
C  to use this routine in scalar mode, pass in ivec=ivecd=1
C
C---------------
C  local:
C
      integer :: i,j,i1,i2,zonrank,linrank,cubrank
      integer :: imaxx,imaxy,imaxdlin,imaxdcub
      logical :: splin_flag
C
      integer :: inum,indx
      integer :: iderivs(2,6),ict4(4),maxd_lin(4),maxd_cub(4)
C
      integer :: idcub,idlin
      real :: xp,dxlin,dxlini,h,hi,ans1,ans2
      real :: f22(2,2)
C
      real, parameter :: ONE = 1.0d0
C---------------
C
      linrank=0
      cubrank=0
      zonrank=0

      do i=1,2
         if(jspline(i).eq.-1) then
            zonrank = zonrank + 1

         else if(jspline(i).eq.0) then
            linrank = linrank + 1

         else
            splin_flag = jspline(i).eq.2
            cubrank = cubrank + 1
         endif
      enddo
C
      if(cubrank.eq.2) then
         if(splin_flag) then
            !  bicubic spline
            call fvbicub(ict,ivec,ivecd,
     >           fval,ii,jj,xparam,yparam,hx,hxi,hy,hyi,
     >           fin,ixdim,iydim)

         else
            !  bicubic Hermite -- translate derivative code
            call dtrans
            do i=1,inum
               if(maxval(iderivs(1:2,i)).le.1) then
                  ict4=0
                  indx=2*iderivs(2,i) + iderivs(1,i) + 1
                  ict4(indx)=1
                  call herm2fcn(ict4,ivec,ivecd,
     >                 fval(1,i),ii,jj,xparam,yparam,hx,hxi,hy,hyi,
     >                 fin,ixdim,iydim)

               else
                  fval(1:ivec,i) = 0
               endif
            enddo
         endif

         return
      endif

      call dtrans

      if(cubrank.eq.0) then
         if(jspline(1).eq.0) then
            imaxx=1
         else
            imaxx=0
         endif

         if(jspline(2).eq.0) then
            imaxy=1
         else
            imaxy=0
         endif

         do i=1,inum
            if((iderivs(1,i).le.imaxx).and.(iderivs(2,i).le.imaxy)) then
               if(linrank.eq.2) then
                  !  bilinear interpolation
                  ict4=0
                  indx=2*iderivs(2,i) + iderivs(1,i) + 1
                  ict4(indx)=1
                  call pc2fcn(ict4,ivec,ivecd,
     >                 fval(1,i),ii,jj,xparam,yparam,hx,hxi,hy,hyi,
     >                 fin,ixdim,iydim)
               else if(zonrank.eq.2) then
                  !  2d step function
                  do j=1,ivec
                     fval(j,i)=fin(1,ii(j),jj(j))
                  enddo
               else
                  !  linear-zone hybrid
                  if(jspline(1).eq.0) then
                     ! linear in x, zonal in y
                     do j=1,ivec
                        if(iderivs(1,i).eq.0) then
                           fval(j,i)=(ONE-xparam(j))*fin(1,ii(j),jj(j))
     >                          +xparam(j)*fin(1,ii(j)+1,jj(j))
                        else
                           fval(j,i)=(fin(1,ii(j)+1,jj(j)) - 
     >                          fin(1,ii(j),jj(j)))*hxi(j)
                        endif
                     enddo
                  else
                     ! linear in y, zonal in x
                     do j=1,ivec
                        if(iderivs(2,i).eq.0) then
                           fval(j,i)=(ONE-yparam(j))*fin(1,ii(j),jj(j))
     >                          +yparam(j)*fin(1,ii(j),jj(j)+1)
                        else
                           fval(j,i)=(fin(1,ii(j),jj(j)+1) - 
     >                          fin(1,ii(j),jj(j)))*hyi(j)
                        endif
                     enddo
                  endif
               endif
            else
               fval(1:ivec,i) = 0
            endif
         enddo

         return
      endif

      !  OK -- Hybrid: Hermite or Spline in one dimension
      !        pclin or step in the other dimension

      if(splin_flag) then
         imaxdcub=3
      else
         imaxdcub=1
      endif

      if(linrank.eq.1) then
         imaxdlin=1
      else
         imaxdlin=0
      endif

      do i=1,inum
         if(jspline(1).ge.1) then
            idcub=iderivs(1,i)
            idlin=iderivs(2,i)
         else
            idcub=iderivs(2,i)
            idlin=iderivs(1,i)
         endif

         if((maxd_lin(i).le.imaxdlin).and. 
     >        (maxd_cub(i).le.imaxdcub)) then
            if(linrank.eq.1) then
               !  spline or Hermite in 1 dimension; pclin in other.
               if(jspline(1).ge.1) then
                  do j=1,ivec
                     xp=xparam(j)
                     dxlin=yparam(j)
                     dxlini=1.0-dxlin
                     h=hx(j)
                     hi=hxi(j)
                     i1=ii(j)
                     i2=jj(j)
                     f22(1:2,1)=fin(1:2,i1,i2)
                     f22(1:2,2)=fin(1:2,i1+1,i2)
                     if(splin_flag) then
                        call sp1d(ans1)
                     else
                        call hm1d(ans1)
                     endif
                     f22(1:2,1)=fin(1:2,i1,i2+1)
                     f22(1:2,2)=fin(1:2,i1+1,i2+1)
                     if(splin_flag) then
                        call sp1d(ans2)
                     else
                        call hm1d(ans2)
                     endif
                     if(idlin.eq.1) then
                        fval(j,i)=(ans2-ans1)*hyi(j)
                     else
                        fval(j,i)=ans2*dxlin + ans1*dxlini
                     endif
                  enddo

               else
                  do j=1,ivec
                     xp=yparam(j)
                     dxlin=xparam(j)
                     dxlini=1.0-dxlin
                     h=hy(j)
                     hi=hyi(j)
                     i1=ii(j)
                     i2=jj(j)
                     f22(1:2,1)=fin(1:2,i1,i2)
                     f22(1:2,2)=fin(1:2,i1,i2+1)
                     if(splin_flag) then
                        call sp1d(ans1)
                     else
                        call hm1d(ans1)
                     endif
                     f22(1:2,1)=fin(1:2,i1+1,i2)
                     f22(1:2,2)=fin(1:2,i1+1,i2+1)
                     if(splin_flag) then
                        call sp1d(ans2)
                     else
                        call hm1d(ans2)
                     endif
                     if(idlin.eq.1) then
                        fval(j,i)=(ans2-ans1)*hxi(j)
                     else
                        fval(j,i)=ans2*dxlin + ans1*dxlini
                     endif
                  enddo
               endif

            else
               !  spline or Hermite in 1 dimension; step fcn in other.
               if(jspline(1).ge.1) then
                  do j=1,ivec
                     xp=xparam(j)
                     h=hx(j)
                     hi=hxi(j)
                     i1=ii(j)
                     i2=jj(j)
                     f22(1:2,1)=fin(1:2,i1,i2)
                     f22(1:2,2)=fin(1:2,i1+1,i2)
                     if(splin_flag) then
                        call sp1d(ans1)
                     else
                        call hm1d(ans1)
                     endif
                     fval(j,i)=ans1
                  enddo
               else
                  do j=1,ivec
                     xp=yparam(j)
                     h=hy(j)
                     hi=hyi(j)
                     i1=ii(j)
                     i2=jj(j)
                     f22(1:2,1)=fin(1:2,i1,i2)
                     f22(1:2,2)=fin(1:2,i1,i2+1)
                     if(splin_flag) then
                        call sp1d(ans1)
                     else
                        call hm1d(ans1)
                     endif
                     fval(j,i)=ans1
                  enddo
               endif
            endif
         else
            fval(1:ivec,i) = 0
         endif
      enddo

      return

      contains
        subroutine dtrans

          ! convert ict(...) codes into intermediate coding:
          !   iderivs(i,j) = number of differentiations in coordinate i
          !                  of the j'th value evaluated

          !   maxd_lin(j) = max # of differentiations for a linear dimension
          !   maxd_cub(j) = max # of differentiations for a spline dimension

          integer :: i

          inum=0   ! actual number of vectors to be evaluated (summed here)

          if(ict(1).le.(2)) then
             if(ict(1).eq.1) then
                call add1(0,0)
             endif
             if(ict(2).eq.1) then
                call add1(1,0)
             endif
             if(ict(3).eq.1) then
                call add1(0,1)
             endif
             if(ict(4).eq.1) then
                call add1(2,0)
             endif
             if(ict(5).eq.1) then
                call add1(0,2)
             endif
             if(ict(6).eq.1) then
                call add1(1,1)
             endif

          else if(ict(1).eq.3) then
             if(ict(2).eq.1) then
                call add1(3,0)
             endif
             if(ict(3).eq.1) then
                call add1(2,1)
             endif
             if(ict(4).eq.1) then
                call add1(1,2)
             endif
             if(ict(5).eq.1) then
                call add1(0,3)
             endif

          else if(ict(1).eq.4) then
             if(ict(2).eq.1) then
                call add1(3,1)
             endif
             if(ict(3).eq.1) then
                call add1(2,2)
             endif
             if(ict(4).eq.1) then
                call add1(1,3)
             endif

          else if(ict(1).eq.5) then
             if(ict(2).eq.1) then
                call add1(3,2)
             endif
             if(ict(3).eq.1) then
                call add1(2,3)
             endif

          else if(ict(1).eq.6) then
             if(ict(2).eq.1) then
                call add1(3,3)
             endif
          endif

        end subroutine dtrans

        subroutine add1(idx,idy)
          ! insert record of derivative d[idx+idy]f/dx[idx]dy[idy]

          integer, intent(in) :: idx,idy

          inum=inum+1

          iderivs(1,inum)=idx
          if(jspline(1).le.0) then
             maxd_lin(inum)=idx
             maxd_cub(inum)=0
          else
             maxd_lin(inum)=0
             maxd_cub(inum)=idx
          endif

          iderivs(2,inum)=idy
          if(jspline(2).le.0) then
             maxd_lin(inum)=max(maxd_lin(inum),idy)
          else
             maxd_cub(inum)=max(maxd_cub(inum),idy)
          endif

        end subroutine add1

C===========================================================================
        subroutine sp1d(ans)
          real, intent(out) :: ans

C----------------------
C contained 1d spline evaluation
C----------------------

          real :: xpi,xp2,xpi2,cx,cxi,hx2,cxd,cxdi
C
          real sixth
          data sixth/0.166666666666666667/
C----------------------

          if(idcub.eq.0) then
             ! value

             xpi=1.0-xp
             xp2=xp*xp
             xpi2=xpi*xpi

             cx=xp*(xp2-1.0)
             cxi=xpi*(xpi2-1.0)
             hx2=h*h

             ans=xpi*f22(1,1) +xp*f22(1,2)
             ans=ans+sixth*hx2*(cxi*f22(2,1) +cx*f22(2,2))

          else if(idcub.eq.1) then
             ! 1st derivative

             xpi=1.0-xp
             xp2=xp*xp
             xpi2=xpi*xpi

             cxd=3.0*xp2-1.0
             cxdi=-3.0*xpi2+1.0

             ans=hi*(f22(1,2)-f22(1,1))
             ans=ans+sixth*h*(cxdi*f22(2,1) +cxd*f22(2,2))

          else if(idcub.eq.2) then
             ! 2nd derivative

             xpi=1.0-xp
             ans=xpi*f22(2,1) + xp*f22(2,2)

          else
             ! 3rd derivative

             ans = hi*(f22(2,2)-f22(2,1))

          endif

        end subroutine sp1d

C===========================================================================
        subroutine hm1d(ans)
          real, intent(out) :: ans

C----------------------
C contained 1d hermite evaluation
C----------------------

          real :: xpi,xp2,xpi2,ax,axbar,bx,bxbar
          real :: axp,axbarp,bxp,bxbarp

          xpi=1.0-xp
          xp2=xp*xp
          xpi2=xpi*xpi

          if(idcub.eq.0) then
             ! value

             ax=xp2*(3.0-2.0*xp)
             axbar=1.0-ax
             bx=-xp2*xpi
             bxbar=xpi2*xp

             ans=axbar*f22(1,1) + ax*f22(1,2)
             ans=ans + h*(bxbar*f22(2,1) + bx*f22(2,2))

          else
             ! 1st derivative

             axp=6.0*xp*xpi
             axbarp=-axp
             bxp=xp*(3.0*xp-2.0)
             bxbarp=xpi*(3.0*xpi-2.0)

             ans=hi*(axbarp*f22(1,1) +axp*f22(1,2))
             ans=ans + bxbarp*f22(2,1) + bxp*f22(2,2)

          endif

        end subroutine hm1d
                
      end
