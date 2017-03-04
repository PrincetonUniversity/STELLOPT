      subroutine evintrp3d(xget,yget,zget,x,nx,y,ny,z,nz,jspline,
     >                   f,icoeff,ixdim,iydim,izdim,ict,fval,ier)
C
      implicit NONE
C
C  use mkintrp3d to set up spline coefficients...
C
C  evaluate a 3d Hybrid interpolant on a rectilinear grid
C
C  this subroutine calls two subroutines:
C     vecin3d_argchk -- error check
C     herm3xyz  -- find cell containing (xget,yget,zget)
C     fvintrp3d  -- evaluate the spline function (w/derivatives if req.)
C
C  input arguments:
C  ================
C
      real xget,yget,zget               ! target of this interpolation
      integer nx,ny,nz                  ! grid sizes
      real x(nx)                        ! ordered x grid
      real y(ny)                        ! ordered y grid
      real z(nz)                        ! ordered z grid
C
      integer :: jspline(3)             ! interpolation method for each
C           dimension: jspline(1) for x; jspline(2) for y; jspline(3) for z
C           =-1: zonal step function; =0: pc lin; =1: Hermite; =2: Spline
C
      integer :: icoeff                 ! no. of coefficients per data point

      integer :: ixdim,iydim,izdim      ! dimensioning:
      !  ixdim=nx-1 for zonal step in x; otherwise nx
      !  iydim=ny-1 for zonal step in y; otherwise ny
      !  izdim=nz-1 for zonal step in z; otherwise nz
C
      real f(icoeff,ixdim,iydim,izdim)  ! function data
C
      integer ict(10)                   ! code specifying output desired
C
C  Note on derivatives: for dimensions along which zonal step function
C    interpolation is done, ANY derivative request returns ZERO.
C    For dimensions along which piecewise linear or Hermite interpolation
C    are done, more than one differentiation returns ZERO!
C
C  Derivative controls are the same as for the compact 3d spline evaluation
C  routine (evtricub):
C
C  ict(1)=1 -- return f  (0, don't)
C  ict(2)=1 -- return df/dx  (0, don't)
C  ict(3)=1 -- return df/dy  (0, don't)
C  ict(4)=1 -- return df/dz  (0, don't)
C  ict(5)=1 -- return d2f/dx2  (0, don't)
C  ict(6)=1 -- return d2f/dy2  (0, don't)
C  ict(7)=1 -- return d2f/dz2  (0, don't)
C  ict(8)=1 -- return d2f/dxdy (0, don't)
C  ict(9)=1 -- return d2f/dxdz (0, don't)
C  ict(10)=1-- return d2f/dydz (0, don't)
c
c  (new dmc Dec 2005 -- higher derivatives available)
c    ict(1)=3 --> 3rd derivative, .le.2 diff. in any coordinate
c      ict(2:8) select: fxxy fxxz fxyy fxyz fxzz fyyz fyzz
c      ->note ict(1)=3, ict(5)=1 gives fxyz = d3f/dxdydz
c    ict(1)=-3 --> 3rd derivative, 3 in one coordinate
c      ict(2:4) select: fxxx fyyy fzzz
c    ict(1)=4 --> 3rd derivative, .le.2 diff. in any coordinate
c      ict(2:7) select: fxxyy fxxyz fxxzz fxyyz fxyzz fyyzz
c    ict(1)=-4 --> 3rd derivative, 3 in one coordinate
c      ict(2:7) select: fxxxy fxxxz fxyyy fxzzz fyyyz fyzzz
c    ict(1)=5 --> 3rd derivative, .le.2 diff. in any coordinate
c      ict(2:4) select: fxxyyz fxxyzz fxyyzz
c    ict(1)=-5 --> 3rd derivative, 3 in one coordinate
c      ict(2:10) select:  fxxxyy fxxxyz fxxxzz fxxyyy fxxzzz
c                         fxyyyz fxyzzz fyyyzz fzzzyy
c    ict(1)=6 --> 3rd derivative, .le.2 diff. in any coordinate
c      fxxyyzz
c    ict(1)=-6 --> 3rd derivative, 3 in one coordinate
c      ict(2:10) select: fxxxyyy fxxxyyz fxxxyzz fxxxyyz
c                        fxxyyyz fxxyzzz fxyyyzz fxyyzzz fyyyzzz
c    ict(1)=-7 --> 7th derivative
c      ict(2:7) select: fxxxyyyz fxxxyyzz fxxxyzzz
c                       fxxyyyzz fxxyyzzz fxyyyzzz
c    ict(1)=-8 --> 8th derivative
c      ict(2:4) select: fxxxyyyzz fxxxyyzzz fxxyyyzzz
c    ict(1)=-9 --> 9th derivative:  fxxxyyyzzz
c
C
C output arguments:
C =================
C
      real fval(10)                     ! output data
      integer ier                       ! error code =0 ==> no error
C
C  fval(1) receives the first output (depends on ict(...) spec)
C  fval(2) receives the second output (depends on ict(...) spec)
C  fval(3) receives the third output (depends on ict(...) spec)
C  fval(4) receives the 4th output (depends on ict(...) spec)
C  fval(5-10) receive 5th thru 10th outputs (if required by ict(...) spec)
C
C  examples:
C    on input ict = [1,1,1,1,0,0,0,0,0,0,0]
C   on output fval= [f,df/dx,df/dy,df/dz]
C
C    on input ict = [1,0,0,0,0,0,0,0,0,0,0]
C   on output fval= [f] ... elements 2-10 never referenced
C
C    on input ict = [0,1,1,0,0,0,0,0,0,0,0]
C   on output fval= [df/dx,df/dy] ... elements 3-10 never referenced
C
C    on input ict = [0,0,0,0,1,0,0,0,0,0,0]
C   on output fval= [d2f/dx2] ... elements 2-10 never referenced.
C
C  ier -- completion code:  0 means OK
C-------------------
C  local:
C
      integer i,j,k                     ! cell indices
C
C  normalized displacement from (x(i),y(j),z(k)) corner of cell.
C    xparam=0 @x(i)  xparam=1 @x(i+1)
C    yparam=0 @y(j)  yparam=1 @y(j+1)
C    zparam=0 @z(k)  zparam=1 @z(k+1)
C
      real xparam,yparam,zparam
C
C  cell dimensions and
C  inverse cell dimensions hxi = 1/(x(i+1)-x(i)), hyi = 1/(y(j+1)-y(j))
C
      real hx,hy,hz
      real hxi,hyi,hzi
C
C  0 .le. xparam .le. 1
C  0 .le. yparam .le. 1
C  0 .le. zparam .le. 1
C
C---------------------------------------------------------------------
C
      call vecin3d_argchk('evintrp3d',jspline,
     >     icoeff,nx,ny,nz,ixdim,iydim,izdim,ier)
      if(ier.ne.0) return
c
C  use lookup routine as in Hermite interpolation
C
      call herm3xyz(xget,yget,zget,x,nx,y,ny,z,nz,0,0,0,
     >     i,j,k,xparam,yparam,zparam,hx,hxi,hy,hyi,hz,hzi,ier)
      if(ier.ne.0) return
c
      call fvintrp3d(ict,1,1,
     >     fval,i,j,k,xparam,yparam,zparam,
     >     hx,hxi,hy,hyi,hz,hzi,
     >     jspline,f,icoeff,ixdim,iydim,izdim)

C
      return
      end
C---------------------------------------------------------------------
C  evaluate Hybrid function interpolation -- 3d fcn
C   NO ERROR CHECKING
C
      subroutine fvintrp3d(ict,ivec,ivecd,
     >     fval,ii,jj,kk,xparam,yparam,zparam,
     >     hx,hxi,hy,hyi,hz,hzi,
     >     jspline,fin,icoeff,ixdim,iydim,izdim)
C
      implicit NONE
C
C  use mkintrp3d to set up spline coefficients...
C
      integer ict(10)                   ! requested output control
      integer ivec                      ! vector length
      integer ivecd                     ! vector dimension (1st dim of fval)
C
      integer ii(ivec),jj(ivec),kk(ivec) ! target cells (i,j,k)
      real xparam(ivec),yparam(ivec),zparam(ivec)
                          ! normalized displacements from (i,j,k) corners
C
      real hx(ivec),hy(ivec),hz(ivec)   ! grid spacing, and
      real hxi(ivec),hyi(ivec),hzi(ivec) ! inverse grid spacing
           ! 1/(x(i+1)-x(i)) & 1/(y(j+1)-y(j)) & 1/(z(k+1)-z(i))
C
      integer :: jspline(3)             ! interpolation type control
      integer :: icoeff                 ! coefficient dimension
      integer :: ixdim,iydim,izdim      ! x & y & z dimensions
      real fin(icoeff,ixdim,iydim,izdim)  ! data on which to interpolate...
C
      real fval(ivecd,10)               ! output returned
C
C  for detailed description of fin, ict and fval see subroutine evintrp3d
C  comments.  Note ict is not vectorized; the same output
C  is expected to be returned for all input vector data points.
C
C  note that the index inputs ii,jj,kk and parameter inputs
C     xparam,yparam,zparam,hx,hxi,hy,hyi,hz,hzi are vectorized, and the
C     output array fval has a vector ** 1st dimension ** whose
C     size must be given as a separate argument
C
C  to use this routine in scalar mode, pass in ivec=ivecd=1
C
C---------------
C  local:
C
      integer :: i,j,k1,k2,imaxx,imaxy,imaxz,linrank,cubrank,zonrank
      logical :: splin_flag
C
      integer :: inum
      integer :: iderivs(3,10),ict8(8),maxd_lin(8),maxd_cub(8)
      integer :: idcub,idlin,idcub2(2),idlin2(2)
C
      real :: xp,yp,hhx,hhy,hhxi,hhyi,fac1,fac2(2)
      real :: h,hi,hhi,h2(2),h2i(2)
      real :: f22(2,2),f22b(2,2),f2222(2,2,2,2),ans22(2,2),ans1,ans2
      real :: f422(4,2,2),f422b(4,2,2)
C
      real, parameter :: ONE = 1.0d0
      real, parameter :: sixth = 0.166666666666666667
      real, parameter :: z36th = sixth*sixth
C
C---------------
C
      zonrank=0
      linrank=0
      cubrank=0

      do i=1,3
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
      if(cubrank.eq.3) then
         if(splin_flag) then
            !  bicubic spline
            call fvtricub(ict,ivec,ivecd,
     >           fval,ii,jj,kk,xparam,yparam,zparam,
     >           hx,hxi,hy,hyi,hz,hzi,
     >           fin,ixdim,iydim,izdim)

         else
            !  bicubic Hermite -- translate derivative code
            call dtrans
            do i=1,inum
               if(maxval(iderivs(1:3,i)).le.1) then
                  call gen_ict8

                  call herm3fcn(ict8,ivec,ivecd,
     >                 fval(1,i),ii,jj,kk,xparam,yparam,zparam,
     >                 hx,hxi,hy,hyi,hz,hzi,
     >                 fin,ixdim,iydim,izdim)

               else
                  fval(1:ivec,i) = 0
               endif
            enddo
         endif

         return
      endif

      call dtrans

      if(cubrank.eq.0) then
         if(jspline(1).eq.-1) then
            imaxx=0
         else
            imaxx=1
         endif
         if(jspline(2).eq.-1) then
            imaxy=0
         else
            imaxy=1
         endif
         if(jspline(3).eq.-1) then
            imaxz=0
         else
            imaxz=1
         endif

         do i=1,inum
            if((iderivs(1,i).le.imaxx).and.(iderivs(2,i).le.imaxy)
     >           .and.(iderivs(3,i).le.imaxz)) then
               if(linrank.eq.3) then
                  !  trilinear interpolation
                  call gen_ict8

                  call pc3fcn(ict8,ivec,ivecd,
     >                 fval(1,i), ii,jj,kk, xparam,yparam,zparam,
     >                 hx,hxi, hy,hyi, hz,hzi,
     >                 fin, ixdim,iydim,izdim)
               else if(zonrank.eq.3) then
                  !  3d step fcn
                  do j=1,ivec
                     fval(j,i)=fin(1,ii(j),jj(j),kk(j))
                  enddo
               else if(linrank.eq.1) then
                  if(jspline(1).eq.0) then
                     ! linear in x, zonal in y,z
                     do j=1,ivec
                        if(iderivs(1,i).eq.0) then
                           fval(j,i) = 
     >                          (ONE-xparam(j))*fin(1,ii(j),jj(j),kk(j))
     >                          +xparam(j)*fin(1,ii(j)+1,jj(j),kk(j))
                        else
                           fval(j,i) = hxi(j)*
     >                          (fin(1,ii(j)+1,jj(j),kk(j))-
     >                          fin(1,ii(j),jj(j),kk(j)))
                        endif
                     enddo
                  else if(jspline(2).eq.0) then
                     ! linear in y, zonal in x,z
                     do j=1,ivec
                        if(iderivs(2,i).eq.0) then
                           fval(j,i) = 
     >                          (ONE-yparam(j))*fin(1,ii(j),jj(j),kk(j))
     >                          +yparam(j)*fin(1,ii(j),jj(j)+1,kk(j))
                        else
                           fval(j,i) = hyi(j)*
     >                          (fin(1,ii(j),jj(j)+1,kk(j))-
     >                          fin(1,ii(j),jj(j),kk(j)))
                        endif
                     enddo
                  else
                     ! linear in z, zonal in x,y
                     do j=1,ivec
                        if(iderivs(3,i).eq.0) then
                           fval(j,i) = 
     >                          (ONE-zparam(j))*fin(1,ii(j),jj(j),kk(j))
     >                          +zparam(j)*fin(1,ii(j),jj(j),kk(j)+1)
                        else
                           fval(j,i) = hzi(j)*
     >                          (fin(1,ii(j),jj(j),kk(j)+1)-
     >                          fin(1,ii(j),jj(j),kk(j)))
                        endif
                     enddo
                  endif
               else
                  if(jspline(1).eq.-1) then
                     ! zonal in x, linear in y,z
                     do j=1,ivec
                        if((iderivs(2,i).eq.0).and.
     >                       (iderivs(3,i).eq.0)) then
                           fval(j,i) = (ONE-zparam(j))*
     >                          ((ONE-yparam(j))*
     >                              fin(1,ii(j),jj(j),kk(j))
     >                            +yparam(j)*fin(1,ii(j),jj(j)+1,kk(j)))
     >                          +zparam(j)*
     >                          ((ONE-yparam(j))*
     >                              fin(1,ii(j),jj(j),kk(j)+1)
     >                          +yparam(j)*fin(1,ii(j),jj(j)+1,kk(j)+1))

                        else if((iderivs(2,i).eq.1).and.
     >                          (iderivs(3,i).eq.0)) then
                           fval(j,i) = hyi(j)*(
     >                        (ONE-zparam(j))*fin(1,ii(j),jj(j)+1,kk(j))
     >                         +zparam(j)*fin(1,ii(j),jj(j)+1,kk(j)+1)
     >                       -(ONE-zparam(j))*fin(1,ii(j),jj(j),kk(j))
     >                         -zparam(j)*fin(1,ii(j),jj(j),kk(j)+1)  )

                        else if((iderivs(2,i).eq.0).and.
     >                          (iderivs(3,i).eq.1)) then
                           fval(j,i) = hzi(j)*(
     >                        (ONE-yparam(j))*fin(1,ii(j),jj(j),kk(j)+1)
     >                         +yparam(j)*fin(1,ii(j),jj(j)+1,kk(j)+1)
     >                       -(ONE-yparam(j))*fin(1,ii(j),jj(j),kk(j))
     >                         -yparam(j)*fin(1,ii(j),jj(j)+1,kk(j))  )

                        else
                           fval(j,i) = hyi(j)*hzi(j)*(
     >                        fin(1,ii(j),jj(j),kk(j))
     >                        -fin(1,ii(j),jj(j)+1,kk(j))
     >                        -fin(1,ii(j),jj(j),kk(j)+1)
     >                        +fin(1,ii(j),jj(j)+1,kk(j)+1))

                        endif
                     enddo
                  else if(jspline(2).eq.-1) then
                     ! zonal in y, linear in x,z
                     do j=1,ivec
                        if((iderivs(1,i).eq.0).and.
     >                       (iderivs(3,i).eq.0)) then
                           fval(j,i) = (ONE-zparam(j))*
     >                          ((ONE-xparam(j))*
     >                              fin(1,ii(j),jj(j),kk(j))
     >                            +xparam(j)*fin(1,ii(j)+1,jj(j),kk(j)))
     >                          +zparam(j)*
     >                          ((ONE-xparam(j))*
     >                              fin(1,ii(j),jj(j),kk(j)+1)
     >                          +xparam(j)*fin(1,ii(j)+1,jj(j),kk(j)+1))

                        else if((iderivs(1,i).eq.1).and.
     >                          (iderivs(3,i).eq.0)) then
                           fval(j,i) = hxi(j)*(
     >                        (ONE-zparam(j))*fin(1,ii(j)+1,jj(j),kk(j))
     >                         +zparam(j)*fin(1,ii(j)+1,jj(j),kk(j)+1)
     >                       -(ONE-zparam(j))*fin(1,ii(j),jj(j),kk(j))
     >                         -zparam(j)*fin(1,ii(j),jj(j),kk(j)+1)  )

                        else if((iderivs(1,i).eq.0).and.
     >                          (iderivs(3,i).eq.1)) then
                           fval(j,i) = hzi(j)*(
     >                        (ONE-xparam(j))*fin(1,ii(j),jj(j),kk(j)+1)
     >                         +xparam(j)*fin(1,ii(j)+1,jj(j),kk(j)+1)
     >                       -(ONE-xparam(j))*fin(1,ii(j),jj(j),kk(j))
     >                         -xparam(j)*fin(1,ii(j)+1,jj(j),kk(j))  )

                        else
                           fval(j,i) = hxi(j)*hzi(j)*(
     >                        fin(1,ii(j),jj(j),kk(j))
     >                        -fin(1,ii(j)+1,jj(j),kk(j))
     >                        -fin(1,ii(j),jj(j),kk(j)+1)
     >                        +fin(1,ii(j)+1,jj(j),kk(j)+1))

                        endif
                     enddo
                  else
                     ! zonal in z, linear in x,y
                     do j=1,ivec
                        if((iderivs(1,i).eq.0).and.
     >                       (iderivs(2,i).eq.0)) then
                           fval(j,i) = (ONE-yparam(j))*
     >                          ((ONE-xparam(j))*
     >                              fin(1,ii(j),jj(j),kk(j))
     >                            +xparam(j)*fin(1,ii(j)+1,jj(j),kk(j)))
     >                          +yparam(j)*
     >                          ((ONE-xparam(j))*
     >                              fin(1,ii(j),jj(j)+1,kk(j))
     >                          +xparam(j)*fin(1,ii(j)+1,jj(j)+1,kk(j)))

                        else if((iderivs(1,i).eq.1).and.
     >                          (iderivs(2,i).eq.0)) then
                           fval(j,i) = hxi(j)*(
     >                        (ONE-yparam(j))*fin(1,ii(j)+1,jj(j),kk(j))
     >                         +yparam(j)*fin(1,ii(j)+1,jj(j)+1,kk(j))
     >                       -(ONE-yparam(j))*fin(1,ii(j),jj(j),kk(j))
     >                         -yparam(j)*fin(1,ii(j),jj(j)+1,kk(j))  )

                        else if((iderivs(1,i).eq.0).and.
     >                          (iderivs(2,i).eq.1)) then
                           fval(j,i) = hyi(j)*(
     >                        (ONE-xparam(j))*fin(1,ii(j),jj(j)+1,kk(j))
     >                         +xparam(j)*fin(1,ii(j)+1,jj(j)+1,kk(j))
     >                       -(ONE-xparam(j))*fin(1,ii(j),jj(j),kk(j))
     >                         -xparam(j)*fin(1,ii(j)+1,jj(j),kk(j))  )

                        else
                           fval(j,i) = hxi(j)*hyi(j)*(
     >                        fin(1,ii(j),jj(j),kk(j))
     >                        -fin(1,ii(j)+1,jj(j),kk(j))
     >                        -fin(1,ii(j),jj(j)+1,kk(j))
     >                        +fin(1,ii(j)+1,jj(j)+1,kk(j)))

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

      !  OK -- Hybrid: Hermite or Spline in one or two dimensions
      !        pclin or step in the other dimension(s).

      if(jspline(1).eq.2) then
         imaxx=3
      else if(jspline(1).eq.-1) then
         imaxx=0
      else
         imaxx=1
      endif

      if(jspline(2).eq.2) then
         imaxy=3
      else if(jspline(2).eq.-1) then
         imaxy=0
      else
         imaxy=1
      endif

      if(jspline(3).eq.2) then
         imaxz=3
      else if(jspline(3).eq.-1) then
         imaxz=0
      else
         imaxz=1
      endif

      if(cubrank.eq.1) then
         !  cubic in 1 dimension; bilinear or lookup in other dimensions

         do i=1,inum
            if(jspline(1).ge.1) then
               idcub=iderivs(1,i)
               idlin2(1)=iderivs(2,i)
               idlin2(2)=iderivs(3,i)
            else if(jspline(2).ge.1) then
               idcub=iderivs(2,i)
               idlin2(1)=iderivs(1,i)
               idlin2(2)=iderivs(3,i)
            else if(jspline(3).ge.1) then
               idcub=iderivs(3,i)
               idlin2(1)=iderivs(1,i)
               idlin2(2)=iderivs(2,i)
            endif

            if((iderivs(1,i).le.imaxx).and.(iderivs(2,i).le.imaxy)
     >           .and.(iderivs(3,i).le.imaxz)) then
               if(linrank.eq.2) then
                  do j=1,ivec
                     if(jspline(1).ge.1) then
                        xp=xparam(j)
                        fac2(1)=yparam(j)
                        fac2(2)=zparam(j)
                        h=hx(j)
                        hi=hxi(j)
                        h2i(1)=hyi(j)
                        h2i(2)=hzi(j)
                        f2222(1:2,1,1:2,1:2) =
     >                       fin(1:2,ii(j),jj(j):jj(j)+1,kk(j):kk(j)+1)
                        f2222(1:2,2,1:2,1:2) =
     >                     fin(1:2,ii(j)+1,jj(j):jj(j)+1,kk(j):kk(j)+1)
                     else if(jspline(2).ge.1) then
                        xp=yparam(j)
                        fac2(1)=xparam(j)
                        fac2(2)=zparam(j)
                        h=hy(j)
                        hi=hyi(j)
                        h2i(1)=hxi(j)
                        h2i(2)=hzi(j)
                        f2222(1:2,1,1:2,1:2) =
     >                       fin(1:2,ii(j):ii(j)+1,jj(j),kk(j):kk(j)+1)
                        f2222(1:2,2,1:2,1:2) =
     >                     fin(1:2,ii(j):ii(j)+1,jj(j)+1,kk(j):kk(j)+1)
                     else if(jspline(3).ge.1) then
                        xp=zparam(j)
                        fac2(1)=xparam(j)
                        fac2(2)=yparam(j)
                        h=hz(j)
                        hi=hzi(j)
                        h2i(1)=hxi(j)
                        h2i(2)=hyi(j)
                        f2222(1:2,1,1:2,1:2) =
     >                       fin(1:2,ii(j):ii(j)+1,jj(j):jj(j)+1,kk(j))
                        f2222(1:2,2,1:2,1:2) =
     >                     fin(1:2,ii(j):ii(j)+1,jj(j):jj(j)+1,kk(j)+1)
                     endif

                     !  evaluate 2x2 array of 1d splines
                     do k2=1,2
                        do k1=1,2
                           if(splin_flag) then
                              call sp1d(ans22(k1,k2),
     >                             f2222(1:2,1:2,k1,k2))
                           else
                              call hm1d(ans22(k1,k2),
     >                             f2222(1:2,1:2,k1,k2))
                           endif
                        enddo
                     enddo

                     !  bilinear interpolation
                     if((idlin2(1).eq.0).and.(idlin2(2).eq.0)) then
                        !  value
                        ans1 = (ONE-fac2(2))*
     >                     ((ONE-fac2(1))*ans22(1,1)+fac2(1)*ans22(2,1))
     >                         +fac2(2)*
     >                     ((ONE-fac2(1))*ans22(1,2)+fac2(1)*ans22(2,2))
                     else if((idlin2(1).eq.1).and.(idlin2(2).eq.0)) then
                        !  1st deriv along 1st dimension
                        ans1 = (ONE-fac2(2))*(ans22(2,1)-ans22(1,1))
     >                        +fac2(2)*(ans22(2,2)-ans22(1,2))
                        ans1 = ans1*h2i(1)
                     else if((idlin2(1).eq.0).and.(idlin2(2).eq.1)) then
                        !  1st deriv along 2nd dimension
                        ans1 = (ONE-fac2(1))*(ans22(1,2)-ans22(1,1))
     >                        +fac2(1)*(ans22(2,2)-ans22(2,1))
                        ans1 = ans1*h2i(2)
                     else
                        !  derivatives along both dimensions
                        ans1 = ans22(2,2)-ans22(1,2)
     >                        -ans22(2,1)+ans22(1,1)
                        ans1 = ans1*h2i(1)*h2i(2)
                     endif
                     fval(j,i)=ans1
                  enddo

               else if(linrank.eq.0) then
                  do j=1,ivec
                     if(jspline(1).ge.1) then
                        xp=xparam(j)
                        h=hx(j)
                        hi=hxi(j)
                        f22(1:2,1)=fin(1:2,ii(j),jj(j),kk(j))
                        f22(1:2,2)=fin(1:2,ii(j)+1,jj(j),kk(j))
                     else if(jspline(2).ge.1) then
                        xp=yparam(j)
                        h=hy(j)
                        hi=hyi(j)
                        f22(1:2,1)=fin(1:2,ii(j),jj(j),kk(j))
                        f22(1:2,2)=fin(1:2,ii(j),jj(j)+1,kk(j))
                     else if(jspline(3).ge.1) then
                        xp=zparam(j)
                        h=hz(j)
                        hi=hzi(j)
                        f22(1:2,1)=fin(1:2,ii(j),jj(j),kk(j))
                        f22(1:2,2)=fin(1:2,ii(j),jj(j),kk(j)+1)
                     endif
                     if(splin_flag) then
                        call sp1d(ans1,f22)
                     else
                        call hm1d(ans1,f22)
                     endif
                     fval(j,i)=ans1
                  enddo

               else if(linrank.eq.1) then
                  do j=1,ivec
                     if(jspline(1).ge.1) then
                        xp=xparam(j)
                        h=hx(j)
                        hi=hxi(j)
                        if(jspline(2).eq.0) then
                           hhi=hyi(j)
                           fac1=yparam(j)
                           idlin=iderivs(2,i)
                           f22(1:2,1)=fin(1:2,ii(j),jj(j),kk(j))
                           f22(1:2,2)=fin(1:2,ii(j)+1,jj(j),kk(j))
                           f22b(1:2,1)=fin(1:2,ii(j),jj(j)+1,kk(j))
                           f22b(1:2,2)=fin(1:2,ii(j)+1,jj(j)+1,kk(j))
                        else
                           hhi=hzi(j)
                           fac1=zparam(j)
                           idlin=iderivs(3,i)
                           f22(1:2,1)=fin(1:2,ii(j),jj(j),kk(j))
                           f22(1:2,2)=fin(1:2,ii(j)+1,jj(j),kk(j))
                           f22b(1:2,1)=fin(1:2,ii(j),jj(j),kk(j)+1)
                           f22b(1:2,2)=fin(1:2,ii(j)+1,jj(j),kk(j)+1)
                        endif
                     else if(jspline(2).ge.1) then
                        xp=yparam(j)
                        h=hy(j)
                        hi=hyi(j)
                        if(jspline(1).eq.0) then
                           hhi=hxi(j)
                           fac1=xparam(j)
                           idlin=iderivs(1,i)
                           f22(1:2,1)=fin(1:2,ii(j),jj(j),kk(j))
                           f22(1:2,2)=fin(1:2,ii(j),jj(j)+1,kk(j))
                           f22b(1:2,1)=fin(1:2,ii(j)+1,jj(j),kk(j))
                           f22b(1:2,2)=fin(1:2,ii(j)+1,jj(j)+1,kk(j))
                        else
                           hhi=hzi(j)
                           fac1=zparam(j)
                           idlin=iderivs(3,i)
                           f22(1:2,1)=fin(1:2,ii(j),jj(j),kk(j))
                           f22(1:2,2)=fin(1:2,ii(j),jj(j)+1,kk(j))
                           f22b(1:2,1)=fin(1:2,ii(j),jj(j),kk(j)+1)
                           f22b(1:2,2)=fin(1:2,ii(j),jj(j)+1,kk(j)+1)
                        endif
                     else if(jspline(3).ge.1) then
                        xp=zparam(j)
                        h=hz(j)
                        hi=hzi(j)
                        if(jspline(1).eq.0) then
                           hhi=hxi(j)
                           fac1=xparam(j)
                           idlin=iderivs(1,i)
                           f22(1:2,1)=fin(1:2,ii(j),jj(j),kk(j))
                           f22(1:2,2)=fin(1:2,ii(j),jj(j),kk(j)+1)
                           f22b(1:2,1)=fin(1:2,ii(j)+1,jj(j),kk(j))
                           f22b(1:2,2)=fin(1:2,ii(j)+1,jj(j),kk(j)+1)
                        else
                           hhi=hyi(j)
                           fac1=yparam(j)
                           idlin=iderivs(2,i)
                           f22(1:2,1)=fin(1:2,ii(j),jj(j),kk(j))
                           f22(1:2,2)=fin(1:2,ii(j),jj(j),kk(j)+1)
                           f22b(1:2,1)=fin(1:2,ii(j),jj(j)+1,kk(j))
                           f22b(1:2,2)=fin(1:2,ii(j),jj(j)+1,kk(j)+1)
                        endif
                     endif
                     if(splin_flag) then
                        call sp1d(ans1,f22)
                        call sp1d(ans2,f22b)
                     else
                        call hm1d(ans1,f22)
                        call hm1d(ans2,f22b)
                     endif
                     if(idlin.eq.0) then
                        fval(j,i)=(ONE-fac1)*ans1 + fac1*ans2
                     else
                        fval(j,i)=(ans2-ans1)*hhi
                     endif
                  enddo

               endif
            else
               fval(1:ivec,i) = 0
            endif
         enddo
         return
      endif

      !  OK:  cubic in 2 dimensions; linear or lookup in other dimension.

      do i=1,inum
         if(jspline(1).le.0) then
            idlin=iderivs(1,i)
            idcub2(1)=iderivs(2,i)
            idcub2(2)=iderivs(3,i)
         else if(jspline(2).le.0) then
            idlin=iderivs(2,i)
            idcub2(1)=iderivs(1,i)
            idcub2(2)=iderivs(3,i)
         else if(jspline(3).le.0) then
            idlin=iderivs(3,i)
            idcub2(1)=iderivs(1,i)
            idcub2(2)=iderivs(2,i)
         endif
         
         if((iderivs(1,i).le.imaxx).and.(iderivs(2,i).le.imaxy)
     >        .and.(iderivs(3,i).le.imaxz)) then

            if(linrank.eq.1) then
               ! linear interpolation in non-spline dimension
               do j=1,ivec
                  if(jspline(1).le.0) then
                     fac1=xparam(j)
                     hi=hxi(j)
                     xp=yparam(j)
                     yp=zparam(j)
                     hhx=hy(j)
                     hhxi=hyi(j)
                     hhy=hz(j)
                     hhyi=hzi(j)
                     f422 = fin(1:4,ii(j),jj(j):jj(j)+1,kk(j):kk(j)+1)
                     f422b= fin(1:4,ii(j)+1,jj(j):jj(j)+1,kk(j):kk(j)+1)
                  else if(jspline(2).le.0) then
                     fac1=yparam(j)
                     hi=hyi(j)
                     xp=xparam(j)
                     yp=zparam(j)
                     hhx=hx(j)
                     hhxi=hxi(j)
                     hhy=hz(j)
                     hhyi=hzi(j)
                     f422 = fin(1:4,ii(j):ii(j)+1,jj(j),kk(j):kk(j)+1)
                     f422b= fin(1:4,ii(j):ii(j)+1,jj(j)+1,kk(j):kk(j)+1)
                  else if(jspline(3).le.0) then
                     fac1=zparam(j)
                     hi=hzi(j)
                     xp=xparam(j)
                     yp=yparam(j)
                     hhx=hx(j)
                     hhxi=hxi(j)
                     hhy=hy(j)
                     hhyi=hyi(j)
                     f422 = fin(1:4,ii(j):ii(j)+1,jj(j):jj(j)+1,kk(j))
                     f422b= fin(1:4,ii(j):ii(j)+1,jj(j):jj(j)+1,kk(j)+1)
                  endif
                  if(splin_flag) then
                     call sp2d(ans1,f422)
                     call sp2d(ans2,f422b)
                  else
                     call hm2d(ans1,f422)
                     call hm2d(ans2,f422b)
                  endif
                  if(idlin.eq.1) then
                     fval(j,i) = (ans2-ans1)*hi
                  else
                     fval(j,i) = (ONE-fac1)*ans1 + fac1*ans2
                  endif
               enddo
            else
               ! lookup in the non-spline dimension
               do j=1,ivec
                  if(jspline(1).le.0) then
                     xp=yparam(j)
                     yp=zparam(j)
                     hhx=hy(j)
                     hhxi=hyi(j)
                     hhy=hz(j)
                     hhyi=hzi(j)
                     f422 = fin(1:4,ii(j),jj(j):jj(j)+1,kk(j):kk(j)+1)
                  else if(jspline(2).le.0) then
                     xp=xparam(j)
                     yp=zparam(j)
                     hhx=hx(j)
                     hhxi=hxi(j)
                     hhy=hz(j)
                     hhyi=hzi(j)
                     f422 = fin(1:4,ii(j):ii(j)+1,jj(j),kk(j):kk(j)+1)
                  else if(jspline(3).le.0) then
                     xp=xparam(j)
                     yp=yparam(j)
                     hhx=hx(j)
                     hhxi=hxi(j)
                     hhy=hy(j)
                     hhyi=hyi(j)
                     f422 = fin(1:4,ii(j):ii(j)+1,jj(j):jj(j)+1,kk(j))
                  endif
                  if(splin_flag) then
                     call sp2d(ans1,f422)
                  else
                     call hm2d(ans1,f422)
                  endif
                  fval(j,i) = ans1
               enddo
            endif
         else
            fval(1:ivec,i) = 0
         endif
      enddo
      contains
        subroutine gen_ict8
          ict8=0
          if(iderivs(1,i).eq.0) then
             if(iderivs(2,i).eq.0) then
                if(iderivs(3,i).eq.0) then
                   ict8(1)=1    ! f
                else
                   ict8(4)=1    ! df/dz
                endif
             else
                if(iderivs(3,i).eq.0) then
                   ict8(3)=1    ! df/dy
                else
                   ict8(7)=1    ! d2f/dydz
                endif
             endif
          else
             if(iderivs(2,i).eq.0) then
                if(iderivs(3,i).eq.0) then
                   ict8(2)=2    ! df/dx
                else
                   ict8(6)=6    ! d2f/dxdz
                endif
             else
                if(iderivs(3,i).eq.0) then
                   ict8(5)=1    ! d2f/dxdy
                else
                   ict8(8)=1    ! d3f/dxdydz
                endif
             endif
          endif
        end subroutine gen_ict8

        subroutine dtrans

          ! convert ict(...) codes into intermediate coding:
          !   iderivs(i,j) = number of differentiations in coordinate i
          !                  of the j'th value evaluated

          !   maxd_lin(j) = max # of differentiations for a linear dimension
          !   maxd_cub(j) = max # of differentiations for a spline dimension

          integer :: i

          inum=0   ! actual number of vectors to be evaluated (summed here)

          if(abs(ict(1)).le.(2)) then
             if(ict(1).eq.1) then
                call add1(0,0,0)
             endif
             if(ict(2).eq.1) then
                call add1(1,0,0)
             endif
             if(ict(3).eq.1) then
                call add1(0,1,0)
             endif
             if(ict(4).eq.1) then
                call add1(0,0,1)
             endif
             if(ict(5).eq.1) then
                call add1(2,0,0)
             endif
             if(ict(6).eq.1) then
                call add1(0,2,0)
             endif
             if(ict(7).eq.1) then
                call add1(0,0,2)
             endif
             if(ict(8).eq.1) then
                call add1(1,1,0)
             endif
             if(ict(9).eq.1) then
                call add1(1,0,1)
             endif
             if(ict(10).eq.1) then
                call add1(0,1,1)
             endif

          else if(ict(1).eq.3) then
             if(ict(2).eq.1) then
                call add1(2,1,0)
             endif
             if(ict(3).eq.1) then
                call add1(2,0,1)
             endif
             if(ict(4).eq.1) then
                call add1(1,2,0)
             endif
             if(ict(5).eq.1) then
                call add1(1,1,1)
             endif
             if(ict(6).eq.1) then
                call add1(1,0,2)
             endif
             if(ict(7).eq.1) then
                call add1(0,2,1)
             endif
             if(ict(8).eq.1) then
                call add1(0,1,2)
             endif

          else if(ict(1).eq.-3) then
             if(ict(2).eq.1) then
                call add1(3,0,0)
             endif
             if(ict(3).eq.1) then
                call add1(0,3,0)
             endif
             if(ict(4).eq.1) then
                call add1(0,0,3)
             endif

          else if(ict(1).eq.4) then
             if(ict(2).eq.1) then
                call add1(2,2,0)
             endif
             if(ict(3).eq.1) then
                call add1(2,1,1)
             endif
             if(ict(4).eq.1) then
                call add1(2,0,2)
             endif
             if(ict(5).eq.1) then
                call add1(1,2,1)
             endif
             if(ict(6).eq.1) then
                call add1(1,1,2)
             endif
             if(ict(7).eq.1) then
                call add1(0,2,2)
             endif

          else if(ict(1).eq.-4) then
             if(ict(2).eq.1) then
                call add1(3,1,0)
             endif
             if(ict(3).eq.1) then
                call add1(3,0,1)
             endif
             if(ict(4).eq.1) then
                call add1(1,3,0)
             endif
             if(ict(5).eq.1) then
                call add1(1,0,3)
             endif
             if(ict(6).eq.1) then
                call add1(0,3,1)
             endif
             if(ict(7).eq.1) then
                call add1(0,1,3)
             endif

          else if(ict(1).eq.5) then
             if(ict(2).eq.1) then
                call add1(2,2,1)
             endif
             if(ict(3).eq.1) then
                call add1(2,1,2)
             endif
             if(ict(4).eq.1) then
                call add1(1,2,2)
             endif

          else if(ict(1).eq.-5) then
             if(ict(2).eq.1) then
                call add1(3,2,0)
             endif
             if(ict(3).eq.1) then
                call add1(3,1,1)
             endif
             if(ict(4).eq.1) then
                call add1(3,0,2)
             endif
             if(ict(5).eq.1) then
                call add1(2,3,0)
             endif
             if(ict(6).eq.1) then
                call add1(2,0,3)
             endif
             if(ict(7).eq.1) then
                call add1(1,3,1)
             endif
             if(ict(8).eq.1) then
                call add1(1,1,3)
             endif
             if(ict(9).eq.1) then
                call add1(0,3,2)
             endif
             if(ict(10).eq.1) then
                call add1(0,2,3)
             endif

          else if(ict(1).eq.6) then
             call add1(2,2,2)

          else if(ict(1).eq.-6) then
             if(ict(2).eq.1) then
                call add1(3,3,0)
             endif
             if(ict(3).eq.1) then
                call add1(3,2,1)
             endif
             if(ict(4).eq.1) then
                call add1(3,1,2)
             endif
             if(ict(5).eq.1) then
                call add1(3,0,3)
             endif
             if(ict(6).eq.1) then
                call add1(2,3,1)
             endif
             if(ict(7).eq.1) then
                call add1(2,1,3)
             endif
             if(ict(8).eq.1) then
                call add1(1,3,2)
             endif
             if(ict(9).eq.1) then
                call add1(1,2,3)
             endif
             if(ict(10).eq.1) then
                call add1(0,3,3)
             endif

          else if(abs(ict(1)).eq.7) then
             if(ict(2).eq.1) then
                call add1(3,3,1)
             endif
             if(ict(3).eq.1) then
                call add1(3,2,2)
             endif
             if(ict(4).eq.1) then
                call add1(3,1,3)
             endif
             if(ict(5).eq.1) then
                call add1(2,3,2)
             endif
             if(ict(6).eq.1) then
                call add1(2,2,3)
             endif
             if(ict(7).eq.1) then
                call add1(1,3,3)
             endif

          else if(abs(ict(1)).eq.8) then
             if(ict(2).eq.1) then
                call add1(3,3,2)
             endif
             if(ict(3).eq.1) then
                call add1(3,2,3)
             endif
             if(ict(4).eq.1) then
                call add1(2,3,3)
             endif

          else if(abs(ict(1)).eq.9) then
             call add1(3,3,3)

          endif

        end subroutine dtrans

        subroutine add1(idx,idy,idz)
          ! insert record of derivative d[idx+idy+idz]f/dx[idx]dy[idy]dz[idz]

          integer, intent(in) :: idx,idy,idz

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


          iderivs(3,inum)=idz
          if(jspline(3).le.0) then
             maxd_lin(inum)=max(maxd_lin(inum),idz)
          else
             maxd_cub(inum)=max(maxd_cub(inum),idz)
          endif

        end subroutine add1

C===========================================================================
        subroutine sp1d(ans,f22)
          real, intent(in) :: f22(2,2)
          real, intent(out) :: ans

C----------------------
C contained 1d spline evaluation
C----------------------

          real :: xpi,xp2,xpi2,cx,cxi,hx2,cxd,cxdi
C
C----------------------

          xpi=1.0-xp
          xp2=xp*xp
          xpi2=xpi*xpi

          if(idcub.eq.0) then
             ! value

             cx=xp*(xp2-1.0)
             cxi=xpi*(xpi2-1.0)
             hx2=h*h

             ans=xpi*f22(1,1) +xp*f22(1,2)
             ans=ans+sixth*hx2*(cxi*f22(2,1) +cx*f22(2,2))

          else if(idcub.eq.1) then
             ! 1st derivative

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
        subroutine hm1d(ans,f22)
          real, intent(in) :: f22(2,2)
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

C===========================================================================
        subroutine sp2d(ans,ff)
          real, intent(in) :: ff(0:3,2,2)
          real, intent(out) :: ans

C----------------------
C contained 2d spline evaluation
C----------------------

          real xpi,xp2,xpi2,cx,cxi,hx2,cxd,cxdi
          real ypi,yp2,ypi2,cy,cyi,hy2,cyd,cydi

C----------------------

          xpi=1.0-xp
          xp2=xp*xp
          xpi2=xpi*xpi

          ypi=1.0-yp
          yp2=yp*yp
          ypi2=ypi*ypi
C
C  get desired values:
C
          if((idcub2(1).eq.0).and.(idcub2(2).eq.0)) then

C  function value:

             cy=yp*(yp2-1.0)
             cyi=ypi*(ypi2-1.0)
             hy2=hhy*hhy

             cx=xp*(xp2-1.0)
             cxi=xpi*(xpi2-1.0)
             hx2=hhx*hhx

             ans=xpi*(ypi*ff(0,1,1)  +yp*ff(0,1,2))+
     >            xp*(ypi*ff(0,2,1)+yp*ff(0,2,2))

             ans=ans+sixth*hx2*(
     >            cxi*(ypi*ff(1,1,1)  +yp*ff(1,1,2))+
     >            cx*(ypi*ff(1,2,1)+yp*ff(1,2,2)))

             ans=ans+sixth*hy2*(
     >            xpi*(cyi*ff(2,1,1)  +cy*ff(2,1,2))+
     >            xp*(cyi*ff(2,2,1)+cy*ff(2,2,2)))

             ans=ans+z36th*hx2*hy2*(
     >            cxi*(cyi*ff(3,1,1)  +cy*ff(3,1,2))+
     >            cx*(cyi*ff(3,2,1)+cy*ff(3,2,2)))


          else if((idcub2(1).eq.1).and.(idcub2(2).eq.0)) then

C  df/dx:

             cy=yp*(yp2-1.0)
             cyi=ypi*(ypi2-1.0)
             hy2=hhy*hhy

             cxd=3.0*xp2-1.0
             cxdi=-3.0*xpi2+1.0

             ans=hhxi*(
     >            -(ypi*ff(0,1,1)  +yp*ff(0,1,2))
     >            +(ypi*ff(0,2,1)+yp*ff(0,2,2)))

             ans=ans+sixth*hhx*(
     >            cxdi*(ypi*ff(1,1,1)  +yp*ff(1,1,2))+
     >            cxd*(ypi*ff(1,2,1)+yp*ff(1,2,2)))

             ans=ans+sixth*hhxi*hy2*(
     >            -(cyi*ff(2,1,1)  +cy*ff(2,1,2))
     >            +(cyi*ff(2,2,1)+cy*ff(2,2,2)))

             ans=ans+z36th*hhx*hy2*(
     >            cxdi*(cyi*ff(3,1,1)  +cy*ff(3,1,2))+
     >            cxd*(cyi*ff(3,2,1)+cy*ff(3,2,2)))


          else if((idcub2(1).eq.0).and.(idcub2(2).eq.1)) then

C  df/dy:

             cyd=3.0*yp2-1.0
             cydi=-3.0*ypi2+1.0

             cx=xp*(xp2-1.0)
             cxi=xpi*(xpi2-1.0)
             hx2=hhx*hhx

             ans=hhyi*(
     >            xpi*(-ff(0,1,1)  +ff(0,1,2))+
     >            xp*(-ff(0,2,1)+ff(0,2,2)))

             ans=ans+sixth*hx2*hhyi*(
     >            cxi*(-ff(1,1,1)  +ff(1,1,2))+
     >            cx*(-ff(1,2,1)+ff(1,2,2)))

             ans=ans+sixth*hhy*(
     >            xpi*(cydi*ff(2,1,1)  +cyd*ff(2,1,2))+
     >            xp*(cydi*ff(2,2,1)+cyd*ff(2,2,2)))

             ans=ans+z36th*hx2*hhy*(
     >            cxi*(cydi*ff(3,1,1)  +cyd*ff(3,1,2))+
     >            cx*(cydi*ff(3,2,1)+cyd*ff(3,2,2)))


          else if((idcub2(1).eq.2).and.(idcub2(2).eq.0)) then

C  d2f/dx2:

             cy=yp*(yp2-1.0)
             cyi=ypi*(ypi2-1.0)
             hy2=hhy*hhy

             ans=(
     >            xpi*(ypi*ff(1,1,1)  +yp*ff(1,1,2))+
     >            xp*(ypi*ff(1,2,1)+yp*ff(1,2,2)))

             ans=ans+sixth*hy2*(
     >            xpi*(cyi*ff(3,1,1)  +cy*ff(3,1,2))+
     >            xp*(cyi*ff(3,2,1)+cy*ff(3,2,2)))


          else if((idcub2(1).eq.0).and.(idcub2(2).eq.2)) then

C  d2f/dy2:

             cx=xp*(xp2-1.0)
             cxi=xpi*(xpi2-1.0)
             hx2=hhx*hhx

             ans=(
     >            xpi*(ypi*ff(2,1,1)  +yp*ff(2,1,2))+
     >            xp*(ypi*ff(2,2,1)+yp*ff(2,2,2)))

             ans=ans+sixth*hx2*(
     >            cxi*(ypi*ff(3,1,1)  +yp*ff(3,1,2))+
     >            cx*(ypi*ff(3,2,1)+yp*ff(3,2,2)))


          else if((idcub2(1).eq.1).and.(idcub2(2).eq.1)) then

C  d2f/dxdy:

             cyd=3.0*yp2-1.0
             cydi=-3.0*ypi2+1.0

             cxd=3.0*xp2-1.0
             cxdi=-3.0*xpi2+1.0

             ans=hhxi*hhyi*(
     >            ff(0,1,1)  -ff(0,1,2)
     >            -ff(0,2,1)+ff(0,2,2))

             ans=ans+sixth*hhx*hhyi*(
     >            cxdi*(-ff(1,1,1)  +ff(1,1,2))+
     >            cxd*(-ff(1,2,1)+ff(1,2,2)))

             ans=ans+sixth*hhxi*hhy*(
     >            -(cydi*ff(2,1,1)  +cyd*ff(2,1,2))
     >            +(cydi*ff(2,2,1)+cyd*ff(2,2,2)))

             ans=ans+z36th*hhx*hhy*(
     >            cxdi*(cydi*ff(3,1,1)  +cyd*ff(3,1,2))+
     >            cxd*(cydi*ff(3,2,1)+cyd*ff(3,2,2)))

          else if((idcub2(1).eq.3).and.(idcub2(2).eq.0)) then

c  evaluate d3f/dx3 (not continuous):

             cy=yp*(yp2-1.0)
             cyi=ypi*(ypi2-1.0)
             hy2=hhy*hhy

             ans=hhxi*(
     >            -(ypi*ff(1,1,1)  +yp*ff(1,1,2))
     >            +(ypi*ff(1,2,1)+yp*ff(1,2,2)))

             ans=ans+sixth*hy2*hhxi*(
     >            -(cyi*ff(3,1,1)  +cy*ff(3,1,2))
     >            +(cyi*ff(3,2,1)+cy*ff(3,2,2)))


          else if((idcub2(1).eq.2).and.(idcub2(2).eq.1)) then

c  evaluate d3f/dx2dy

             cyd=3.0*yp2-1.0
             cydi=-3.0*ypi2+1.0

             ans=hhyi*(
     >            xpi*(-ff(1,1,1)  +ff(1,1,2))+
     >            xp*(-ff(1,2,1) +ff(1,2,2)))

             ans=ans+sixth*hhy*(
     >            xpi*(cydi*ff(3,1,1) +cyd*ff(3,1,2))+
     >            xp*(cydi*ff(3,2,1)+cyd*ff(3,2,2)))


          else if((idcub2(1).eq.1).and.(idcub2(2).eq.2)) then

c  evaluate d3f/dxdy2

             cxd=3.0*xp2-1.0
             cxdi=-3.0*xpi2+1.0

             ans=hhxi*(
     >            -(ypi*ff(2,1,1)  +yp*ff(2,1,2))
     >            +(ypi*ff(2,2,1)+yp*ff(2,2,2)))

             ans=ans+sixth*hhx*(
     >            cxdi*(ypi*ff(3,1,1)  +yp*ff(3,1,2))+
     >            cxd*(ypi*ff(3,2,1)+yp*ff(3,2,2)))

          else if((idcub2(1).eq.0).and.(idcub2(2).eq.3)) then

c  evaluate d3f/dy3 (not continuous):

             cx=xp*(xp2-1.0)
             cxi=xpi*(xpi2-1.0)
             hx2=hhx*hhx

             ans=hhyi*(
     >            xpi*(-ff(2,1,1)  +ff(2,1,2))+
     >            xp*(-ff(2,2,1) +ff(2,2,2)))

             ans=ans+sixth*hx2*hhyi*(
     >            cxi*(-ff(3,1,1)  +ff(3,1,2))+
     >            cx*(-ff(3,2,1) +ff(3,2,2)))


          else if((idcub2(1).eq.3).and.(idcub2(2).eq.1)) then

c  evaluate d4f/dx3dy (not continuous):

             cyd=3.0*yp2-1.0
             cydi=-3.0*ypi2+1.0

             ans=hhxi*hhyi*(
     >            +( ff(1,1,1)  -ff(1,1,2))
     >            +(-ff(1,2,1)+ff(1,2,2)))

             ans=ans+sixth*hhy*hhxi*(
     >            -(cydi*ff(3,1,1)  +cyd*ff(3,1,2))
     >            +(cydi*ff(3,2,1)+cyd*ff(3,2,2)))


          else if((idcub2(1).eq.2).and.(idcub2(2).eq.2)) then

c  evaluate d4f/dx2dy2

             ans=xpi*(ypi*ff(3,1,1)  +yp*ff(3,1,2))+
     >            xp*(ypi*ff(3,2,1)  +yp*ff(3,2,2))

          else if((idcub2(1).eq.1).and.(idcub2(2).eq.3)) then

c  evaluate d4f/dxdy3 (not continuous):

             cxd=3.0*xp2-1.0
             cxdi=-3.0*xpi2+1.0

             ans=hhyi*hhxi*(
     >            +( ff(2,1,1)  -ff(2,1,2))
     >            +(-ff(2,2,1)+ff(2,2,2)))

             ans=ans+sixth*hhx*hhyi*(
     >            cxdi*(-ff(3,1,1)  +ff(3,1,2))+
     >            cxd*(-ff(3,2,1) +ff(3,2,2)))


          else if((idcub2(1).eq.3).and.(idcub2(2).eq.2)) then

c  evaluate d5f/dx3dy2 (not continuous):

             ans=hhxi*(
     >            -(ypi*ff(3,1,1)  +yp*ff(3,1,2))
     >            +(ypi*ff(3,2,1)+yp*ff(3,2,2)))

          else if((idcub2(1).eq.2).and.(idcub2(2).eq.3)) then

c  evaluate d5f/dx2dy3 (not continuous):

             ans=hhyi*(
     >            xpi*(-ff(3,1,1)  +ff(3,1,2))+
     >            xp*(-ff(3,2,1)+ff(3,2,2)))


          else if((idcub2(1).eq.3).and.(idcub2(2).eq.3)) then

c  evaluate d6f/dx3dy3 (not continuous):

             ans=hhxi*hhyi*(
     >            +( ff(3,1,1)  -ff(3,1,2))
     >            +(-ff(3,2,1)  +ff(3,2,2)))

          endif

        end subroutine sp2d

C===========================================================================
        subroutine hm2d(ans,ff)
          real, intent(in) :: ff(0:3,2,2)
          real, intent(out) :: ans

C----------------------
C contained 2d hermite evaluation
C----------------------

          real :: xpi,xp2,xpi2,ax,axbar,bx,bxbar
          real :: axp,axbarp,bxp,bxbarp

          real :: ypi,yp2,ypi2,ay,aybar,by,bybar
          real :: ayp,aybarp,byp,bybarp

          xpi=1.0-xp
          xp2=xp*xp
          xpi2=xpi*xpi

          ypi=1.0-yp
          yp2=yp*yp
          ypi2=ypi*ypi

          if((idcub2(1).eq.0).and.(idcub2(2).eq.0)) then
             !  evaluation

             ax=xp2*(3.0-2.0*xp)
             axbar=1.0-ax
             bx=-xp2*xpi
             bxbar=xpi2*xp

             ay=yp2*(3.0-2.0*yp)
             aybar=1.0-ay
             by=-yp2*ypi
             bybar=ypi2*yp

             ans=axbar*(aybar*ff(0,1,1)  +ay*ff(0,1,2))+
     >            ax*(aybar*ff(0,2,1)+ay*ff(0,2,2))

             ans=ans+hhx*(
     >            bxbar*(aybar*ff(1,1,1)  +ay*ff(1,1,2))+
     >            bx*(aybar*ff(1,2,1)+ay*ff(1,2,2)))

             ans=ans+hhy*(
     >            axbar*(bybar*ff(2,1,1)  +by*ff(2,1,2))+
     >            ax*(bybar*ff(2,2,1)+by*ff(2,2,2)))

             ans=ans+hhx*hhy*(
     >            bxbar*(bybar*ff(3,1,1)  +by*ff(3,1,2))+
     >            bx*(bybar*ff(3,2,1)+by*ff(3,2,2)))

          else if((idcub2(1).eq.1).and.(idcub2(2).eq.0)) then
             !  1st deriv, 1st dim

             axp=6.0*xp*xpi
             axbarp=-axp
             bxp=xp*(3.0*xp-2.0)
             bxbarp=xpi*(3.0*xpi-2.0)

             ay=yp2*(3.0-2.0*yp)
             aybar=1.0-ay
             by=-yp2*ypi
             bybar=ypi2*yp

             ans=hhxi*(
     >            axbarp*(aybar*ff(0,1,1)  +ay*ff(0,1,2))+
     >            axp*(aybar*ff(0,2,1)+ay*ff(0,2,2)))

             ans=ans+
     >            bxbarp*(aybar*ff(1,1,1)  +ay*ff(1,1,2))+
     >            bxp*(aybar*ff(1,2,1)+ay*ff(1,2,2))

             ans=ans+hhxi*hhy*(
     >            axbarp*(bybar*ff(2,1,1)  +by*ff(2,1,2))+
     >            axp*(bybar*ff(2,2,1)+by*ff(2,2,2)))

             ans=ans+hhy*(
     >            bxbarp*(bybar*ff(3,1,1)  +by*ff(3,1,2))+
     >            bxp*(bybar*ff(3,2,1)+by*ff(3,2,2)))

          else if((idcub2(1).eq.0).and.(idcub2(2).eq.1)) then
             !  1st deriv, 2nd dim

             ax=xp2*(3.0-2.0*xp)
             axbar=1.0-ax
             bx=-xp2*xpi
             bxbar=xpi2*xp

             ayp=6.0*yp*ypi
             aybarp=-ayp
             byp=yp*(3.0*yp-2.0)
             bybarp=ypi*(3.0*ypi-2.0)

             ans=hhyi*(
     >            axbar*(aybarp*ff(0,1,1)  +ayp*ff(0,1,2))+
     >            ax*(aybarp*ff(0,2,1)+ayp*ff(0,2,2)))

             ans=ans+hhx*hhyi*(
     >            bxbar*(aybarp*ff(1,1,1)  +ayp*ff(1,1,2))+
     >            bx*(aybarp*ff(1,2,1)+ayp*ff(1,2,2)))

             ans=ans+
     >            axbar*(bybarp*ff(2,1,1)  +byp*ff(2,1,2))+
     >            ax*(bybarp*ff(2,2,1)+byp*ff(2,2,2))

             ans=ans+hhx*(
     >            bxbar*(bybarp*ff(3,1,1)  +byp*ff(3,1,2))+
     >            bx*(bybarp*ff(3,2,1)+byp*ff(3,2,2)))

          else if((idcub2(1).eq.1).and.(idcub2(2).eq.1)) then
             !  1st deriv, in both dimensions

             axp=6.0*xp*xpi
             axbarp=-axp
             bxp=xp*(3.0*xp-2.0)
             bxbarp=xpi*(3.0*xpi-2.0)

             ayp=6.0*yp*ypi
             aybarp=-ayp
             byp=yp*(3.0*yp-2.0)
             bybarp=ypi*(3.0*ypi-2.0)

             ans=hhxi*hhyi*(
     >            axbarp*(aybarp*ff(0,1,1)  +ayp*ff(0,1,2))+
     >            axp*(aybarp*ff(0,2,1)+ayp*ff(0,2,2)))

             ans=ans+hhyi*(
     >            bxbarp*(aybarp*ff(1,1,1)  +ayp*ff(1,1,2))+
     >            bxp*(aybarp*ff(1,2,1)+ayp*ff(1,2,2)))

             ans=ans+hhxi*(
     >            axbarp*(bybarp*ff(2,1,1)  +byp*ff(2,1,2))+
     >            axp*(bybarp*ff(2,2,1)+byp*ff(2,2,2)))

             ans=ans+
     >            bxbarp*(bybarp*ff(3,1,1)  +byp*ff(3,1,2))+
     >            bxp*(bybarp*ff(3,2,1)+byp*ff(3,2,2))

          endif

        end subroutine hm2d

      end
