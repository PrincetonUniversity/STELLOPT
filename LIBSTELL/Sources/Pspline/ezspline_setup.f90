!/////
! R8 !
!/////
subroutine EZspline_setup1_r8(spline_o, f, ier, exact_dim)
  use ezspline_obj
  implicit none
  type(EZspline1_r8) spline_o
  real(ezspline_r8), dimension(:), intent(in) :: f
  ! ier:
  ! 0=ok
  ! 98=some error occurred in EZspline_setup
  integer, intent(out) :: ier

  logical, intent(in), OPTIONAL :: exact_dim  ! set T to require exact 
  !  dimensioning match between f and spline_o%fspl; default is F

  logical :: iexact
  integer ifail
 
  integer :: ipx
  integer iper, imsg, itol, inum, in1
  real(ezspline_r8) ztol, df1, df2
 
  !-------------------------
  iexact=.FALSE.
  if(present(exact_dim)) iexact = exact_dim

  if( .not.EZspline_allocated(spline_o) ) then
     ier = 98
     return
  endif
 
  in1 = size(spline_o%fspl,2)

  ier = 57
  if(size(f,1).lt.in1) return

  if(iexact) then
     ier = 58
     if(size(f,1).gt.in1) return
  endif
  ier = 0

  !
  ! recompute min/max in case user changed the grid manually
  spline_o%x1max = maxval(spline_o%x1)
  spline_o%x1min = minval(spline_o%x1)
  !
  ! set xpkg, useful for cloud and array interpolation
  imsg=0
  itol=0        ! range tolerance option
  ztol=5.e-7_ezspline_r8 ! range tolerance, if itol is set
  iper=0
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) iper=1
  call r8genxpkg(spline_o%n1, spline_o%x1(1), spline_o%x1pkg(1,1),&
       & iper,imsg,itol,ztol, spline_o%klookup1 ,ifail)
  if(ifail/=0) ier=27
 
  spline_o%isReady = 0
 
  spline_o%fspl(1, 1:in1) = &
       &  f(1:in1)
 
  if (spline_o%isHermite == 0 .and. spline_o%isLinear == 0) then
 
     call r8mkspline(             &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%fspl(1,1), &
          &   spline_o%ibctype1(1), spline_o%bcval1min, &
          &   spline_o%ibctype1(2), spline_o%bcval1max, &
          &   spline_o%ilin1, ifail)
 
     if(ifail /= 0) then
        ier = 98
     else
        spline_o%isReady = 1
     endif

  else if (spline_o%isLinear == 1) then

     spline_o%ilin1=0
     if(in1.gt.2) then
        if(spline_o%x1pkg(3,4).eq.0.0_ezspline_r8) spline_o%ilin1=1  ! evenly spaced grid
     else
        spline_o%ilin1=1
     end if

     spline_o%isReady = 1   ! no coefficient setup necessary
 
  else
     !
     ! Hermite polynomial coefficient setup based on Akima's method
     ! (df/dx..etc not required)
     !
     ipx = 0
     if (spline_o%ibctype1(1)==-1 .or. spline_o%ibctype1(2)==-1) then
        ipx=1
     else if (spline_o%ibctype1(1)<-1 .or. spline_o%ibctype1(1)>1 .or. &
          spline_o%ibctype1(2)<-1 .or. spline_o%ibctype1(2)>1 ) then
        ipx=99  ! an error...
     else if (spline_o%ibctype1(1)==1 .or. spline_o%ibctype1(2)==1) then
        ipx=2
        if(spline_o%ibctype1(1)==1) then
           spline_o%fspl(2,1)=spline_o%bcval1min
        else
           ! hand implemented default BC
           df1=(spline_o%fspl(1,2)-spline_o%fspl(1,1))/ &
                (spline_o%x1(2)-spline_o%x1(1))
           df2=(spline_o%fspl(1,3)-spline_o%fspl(1,2))/ &
                (spline_o%x1(3)-spline_o%x1(2))
           spline_o%fspl(2,1)=(3*df1-df2)/2
        endif
        inum=spline_o%n1
        if(spline_o%ibctype1(2)==1) then
           spline_o%fspl(2,inum)=spline_o%bcval1max
        else
           ! hand implemented default BC
           df1=(spline_o%fspl(1,inum)-spline_o%fspl(1,inum-1))/ &
                (spline_o%x1(inum)-spline_o%x1(inum-1))
           df2=(spline_o%fspl(1,inum-1)-spline_o%fspl(1,inum-2))/ &
                (spline_o%x1(inum-1)-spline_o%x1(inum-2))
           spline_o%fspl(2,inum)=(3*df1-df2)/2
        endif
     endif
     ifail = 0
     call r8akherm1p(spline_o%x1(1), spline_o%n1, &
          & spline_o%fspl(1,1), &
          & spline_o%ilin1, &
          & ipx, ifail)
 
     if (ifail /=0 ) then
        ier = 91
     else
        spline_o%isReady = 1
     endif
 
  endif
 
 
end subroutine EZspline_setup1_r8
 
 
subroutine EZspline_setup2_r8(spline_o, f, ier, exact_dim)
  use ezspline_obj
  implicit none
  type(EZspline2_r8) spline_o
  real(ezspline_r8), dimension(:,:), intent(in) :: f
  ! ier:
  ! 0=ok
  ! 98=some error occurred in EZspline_setup
  integer, intent(out) :: ier

  logical, intent(in), OPTIONAL :: exact_dim  ! set T to require exact
  !  dimensioning match between f and spline_o%fspl; default is F

  logical :: iexact
  integer ifail
 
  integer :: ipx, ipy
  integer iper, imsg, itol, inum, ii, jj, in0, in1, in2
  real(ezspline_r8) ztol, df1, df2
 
  !-------------------------
  iexact=.FALSE.
  if(present(exact_dim)) iexact = exact_dim
 
  if( .not.EZspline_allocated(spline_o) ) then
     ier = 98
     return
  endif
 
  in0 = size(spline_o%fspl,1)
  in1 = size(spline_o%fspl,2)
  in2 = size(spline_o%fspl,3)

  ier = 57
  if(size(f,1).lt.in1) return
  if(size(f,2).lt.in2) return

  if(iexact) then
     ier = 58
     if(size(f,1).gt.in1) return
     if(size(f,2).gt.in2) return
  endif
  ier = 0
 
  !
  ! recompute min/max in case user changed the grid manually
  spline_o%x1max = maxval(spline_o%x1)
  spline_o%x2max = maxval(spline_o%x2)
  spline_o%x1min = minval(spline_o%x1)
  spline_o%x2min = minval(spline_o%x2)
 
  !
  ! set xpkg, useful for cloud and array interpolation
  imsg=0
  itol=0        ! range tolerance option
  ztol=5.e-7_ezspline_r8 ! range tolerance, if itol is set
  iper=0
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) iper=1
  call r8genxpkg(spline_o%n1,spline_o%x1(1),spline_o%x1pkg(1,1),&
       & iper,imsg,itol,ztol,spline_o%klookup1,ifail)
  if(ifail/=0) ier=27
  iper=0
  if(spline_o%ibctype2(1)==-1 .OR. spline_o%ibctype2(2)==-1) iper=1
  call r8genxpkg(spline_o%n2,spline_o%x2(1),spline_o%x2pkg(1,1),&
       & iper,imsg,itol,ztol,spline_o%klookup2,ifail)
  if(ifail/=0) ier=27
 
  spline_o%isReady = 0
 
  spline_o%fspl(1, 1:in1, 1:in2) = &
       &  f(1:in1, 1:in2)

  ! this fixes a VMS f90 compiler optimizer problem:
  if(ztol.eq.-1.2345d30) &
	write(6,*) 'spline_o%fspl(1,1,1) = ', spline_o%fspl(1,1,1)
 
  if (spline_o%isHybrid == 1) then
 
     call r8mkintrp2d( &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%hspline, spline_o%fspl(1,1,1), &
          &   in0,in1,in2, &
          &   spline_o%ibctype1(1), spline_o%bcval1min(1), &
          &   spline_o%ibctype1(2), spline_o%bcval1max(1), &
          &   spline_o%ibctype2(1), spline_o%bcval2min(1), &
          &   spline_o%ibctype2(2), spline_o%bcval2max(1), &
          &   ifail)

     spline_o%ilin1 = 0  ! Hybrid does not compute this
     spline_o%ilin2 = 0  ! Hybrid does not compute this

     if(ifail /= 0) then
        ier = 98
     else
        spline_o%isReady = 1
     endif

  else if (spline_o%isHermite == 0 .and. spline_o%isLinear == 0) then
 
     call r8mkbicub(             &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%fspl(1,1,1), spline_o%n1, &
          &   spline_o%ibctype1(1), spline_o%bcval1min(1), &
          &   spline_o%ibctype1(2), spline_o%bcval1max(1), &
          &   spline_o%ibctype2(1), spline_o%bcval2min(1), &
          &   spline_o%ibctype2(2), spline_o%bcval2max(1), &
          &   spline_o%ilin1, spline_o%ilin2, ifail)
 
     if(ifail /= 0) then
        ier = 98
     else
        spline_o%isReady = 1
     endif
 
  else if (spline_o%isLinear == 1) then

     spline_o%ilin1=0
     if(in1.gt.2) then
        if(spline_o%x1pkg(3,4).eq.0.0_ezspline_r8) spline_o%ilin1=1  ! evenly spaced grid
     else
        spline_o%ilin1=1
     end if

     spline_o%ilin2=0
     if(in2.gt.2) then
        if(spline_o%x2pkg(3,4).eq.0.0_ezspline_r8) spline_o%ilin2=1  ! evenly spaced grid
     else
        spline_o%ilin2=1
     end if

     spline_o%isReady = 1   ! no coefficient setup necessary

  else
     !
     ! Hermite polynomial coefficient setup based on Akima's method
     ! (df/dx..etc not required)
     !
     ipx = 0
     if (spline_o%ibctype1(1)==-1 .or. spline_o%ibctype1(2)==-1) then
        ipx=1
     else if (spline_o%ibctype1(1)<-1 .or. spline_o%ibctype1(1)>1 .or. &
          spline_o%ibctype1(2)<-1 .or. spline_o%ibctype1(2)>1 ) then
        ipx=99  ! an error...
     else if (spline_o%ibctype1(1)==1 .or. spline_o%ibctype1(2)==1) then
        ipx=2
        do jj=1,spline_o%n2
           if(spline_o%ibctype1(1)==1) then
              spline_o%fspl(2,1,jj)=spline_o%bcval1min(jj)
           else
              ! hand implemented default BC
              df1=(spline_o%fspl(1,2,jj)-spline_o%fspl(1,1,jj))/ &
                   (spline_o%x1(2)-spline_o%x1(1))
              df2=(spline_o%fspl(1,3,jj)-spline_o%fspl(1,2,jj))/ &
                   (spline_o%x1(3)-spline_o%x1(2))
              spline_o%fspl(2,1,jj)=(3*df1-df2)/2
           endif
           inum=spline_o%n1
           if(spline_o%ibctype1(2)==1) then
              spline_o%fspl(2,inum,jj)=spline_o%bcval1max(jj)
           else
              ! hand implemented default BC
              df1=(spline_o%fspl(1,inum,jj)-spline_o%fspl(1,inum-1,jj))/ &
                   (spline_o%x1(inum)-spline_o%x1(inum-1))
              df2=(spline_o%fspl(1,inum-1,jj)-spline_o%fspl(1,inum-2,jj))/ &
                   (spline_o%x1(inum-1)-spline_o%x1(inum-2))
              spline_o%fspl(2,inum,jj)=(3*df1-df2)/2
           endif
        enddo
     endif
     ipy = 0
     if (spline_o%ibctype2(1)==-1 .or. spline_o%ibctype2(2)==-1) then
        ipy=1
     else if (spline_o%ibctype2(1)<-1 .or. spline_o%ibctype2(1)>1 .or. &
          spline_o%ibctype2(2)<-1 .or. spline_o%ibctype2(2)>1 ) then
        ipy=99  ! an error...
     else if (spline_o%ibctype2(1)==1 .or. spline_o%ibctype2(2)==1) then
        ipy=2
        do ii=1,spline_o%n1
           if(spline_o%ibctype2(1)==1) then
              spline_o%fspl(3,ii,1)=spline_o%bcval2min(ii)
           else
              ! hand implemented default BC
              df1=(spline_o%fspl(1,ii,2)-spline_o%fspl(1,ii,1))/ &
                   (spline_o%x2(2)-spline_o%x2(1))
              df2=(spline_o%fspl(1,ii,3)-spline_o%fspl(1,ii,2))/ &
                   (spline_o%x2(3)-spline_o%x2(2))
              spline_o%fspl(3,ii,1)=(3*df1-df2)/2
           endif
           inum=spline_o%n2
           if(spline_o%ibctype2(2)==1) then
              spline_o%fspl(3,ii,inum)=spline_o%bcval2max(ii)
           else
              ! hand implemented default BC
              df1=(spline_o%fspl(1,ii,inum)-spline_o%fspl(1,ii,inum-1))/ &
                   (spline_o%x2(inum)-spline_o%x2(inum-1))
              df2=(spline_o%fspl(1,ii,inum-1)-spline_o%fspl(1,ii,inum-2))/ &
                   (spline_o%x2(inum-1)-spline_o%x2(inum-2))
              spline_o%fspl(3,ii,inum)=(3*df1-df2)/2
           endif
        enddo
     endif
     ifail = 0
     call r8akherm2p(spline_o%x1(1), spline_o%n1, &
          &         spline_o%x2(1), spline_o%n2, &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & spline_o%ilin1, spline_o%ilin2, &
          & ipx,ipy, ifail)
 
     if (ifail /=0 ) then
        ier = 91
     else
        spline_o%isReady = 1
     endif
 
  endif
 
 
end subroutine EZspline_setup2_r8
 
 
 
subroutine EZspline_setup3_r8(spline_o, f, ier, exact_dim)
  use ezspline_obj
  implicit none
  type(EZspline3_r8) spline_o
  real(ezspline_r8), dimension(:,:,:), intent(in) :: f
  ! ier:
  ! 0=ok
  ! 98=some error occurred in EZspline_setup
  integer, intent(out) :: ier

  logical, intent(in), OPTIONAL :: exact_dim  ! set T to require exact
  !  dimensioning match between f and spline_o%fspl; default is F

  logical :: iexact
  integer ifail
 
  integer :: ipx, ipy, ipz
  integer iper, imsg, itol, inum, ii, jj, in0, in1, in2, in3
  real(ezspline_r8) ztol, df1, df2
 
  !-------------------------
  iexact=.FALSE.
  if(present(exact_dim)) iexact = exact_dim
 
  if( .not.EZspline_allocated(spline_o) ) then
     ier = 98
     return
  endif
 
  in0 = size(spline_o%fspl,1)
  in1 = size(spline_o%fspl,2)
  in2 = size(spline_o%fspl,3)
  in3 = size(spline_o%fspl,4)

  ier = 57
  if(size(f,1).lt.in1) return
  if(size(f,2).lt.in2) return
  if(size(f,3).lt.in3) return

  if(iexact) then
     ier = 58
     if(size(f,1).gt.in1) return
     if(size(f,2).gt.in2) return
     if(size(f,3).gt.in3) return
  endif
  ier = 0
 
  !
  ! recompute min/max in case user changed the grid manually
  spline_o%x1max = maxval(spline_o%x1)
  spline_o%x2max = maxval(spline_o%x2)
  spline_o%x3max = maxval(spline_o%x3)
  spline_o%x1min = minval(spline_o%x1)
  spline_o%x2min = minval(spline_o%x2)
  spline_o%x3min = minval(spline_o%x3)
 
  !
  ! set xpkg, useful for cloud and array interpolation
  imsg=0
  itol=0        ! range tolerance option
  ztol=5.e-7_ezspline_r8 ! range tolerance, if itol is set
  iper=0
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) iper=1
  call r8genxpkg(spline_o%n1,spline_o%x1(1),spline_o%x1pkg(1,1),&
       & iper,imsg,itol,ztol,spline_o%klookup1,ifail)
  if(ifail/=0) ier=27
  iper=0
  if(spline_o%ibctype2(1)==-1 .OR. spline_o%ibctype2(2)==-1) iper=1
  call r8genxpkg(spline_o%n2,spline_o%x2(1),spline_o%x2pkg(1,1),&
       & iper,imsg,itol,ztol,spline_o%klookup2,ifail)
  if(ifail/=0) ier=27
  iper=0
  if(spline_o%ibctype3(1)==-1 .OR. spline_o%ibctype3(2)==-1) iper=1
  call r8genxpkg(spline_o%n3,spline_o%x3(1),spline_o%x3pkg(1,1),&
       & iper,imsg,itol,ztol,spline_o%klookup3,ifail)
  if(ifail/=0) ier=27
 
  spline_o%isReady = 0
 
  spline_o%fspl(1, 1:in1, 1:in2, 1:in3) = &
       &  f(1:in1, 1:in2, 1:in3)
 
  if (spline_o%isHybrid == 1) then
 
     call r8mkintrp3d( &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%x3(1), spline_o%n3, &
          &   spline_o%hspline, spline_o%fspl(1,1,1,1), &
          &   in0,in1,in2,in3, &
          &   spline_o%ibctype1(1), spline_o%bcval1min(1,1), &
          &   spline_o%ibctype1(2), spline_o%bcval1max(1,1), &
          &   spline_o%ibctype2(1), spline_o%bcval2min(1,1), &
          &   spline_o%ibctype2(2), spline_o%bcval2max(1,1), &
          &   spline_o%ibctype3(1), spline_o%bcval3min(1,1), &
          &   spline_o%ibctype3(2), spline_o%bcval3max(1,1), &
          &   ifail)

     spline_o%ilin1 = 0  ! Hybrid does not compute this
     spline_o%ilin2 = 0  ! Hybrid does not compute this
     spline_o%ilin3 = 0  ! Hybrid does not compute this

     if(ifail /= 0) then
        ier = 98
     else
        spline_o%isReady = 1
     endif

  else if (spline_o%isHermite == 0 .and. spline_o%isLinear == 0) then
 
     call r8mktricub(             &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%x3(1), spline_o%n3, &
          &   spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          &   spline_o%ibctype1(1), spline_o%bcval1min(1,1), &
          &   spline_o%ibctype1(2), spline_o%bcval1max(1,1), spline_o%n2, &
          &   spline_o%ibctype2(1), spline_o%bcval2min(1,1), &
          &   spline_o%ibctype2(2), spline_o%bcval2max(1,1), spline_o%n1, &
          &   spline_o%ibctype3(1), spline_o%bcval3min(1,1), &
          &   spline_o%ibctype3(2), spline_o%bcval3max(1,1), spline_o%n1, &
          &   spline_o%ilin1, spline_o%ilin2, spline_o%ilin3, ifail)
 
     if(ifail /= 0) then
        ier = 98
     else
        spline_o%isReady = 1
     endif
 
  else if (spline_o%isLinear == 1) then

     spline_o%ilin1=0
     if(in1.gt.2) then
        if(spline_o%x1pkg(3,4).eq.0.0_ezspline_r8) spline_o%ilin1=1  ! evenly spaced grid
     else
        spline_o%ilin1=1
     end if

     spline_o%ilin2=0
     if(in2.gt.2) then
        if(spline_o%x2pkg(3,4).eq.0.0_ezspline_r8) spline_o%ilin2=1  ! evenly spaced grid
     else
        spline_o%ilin2=1
     end if

     spline_o%ilin3=0
     if(in3.gt.2) then
        if(spline_o%x3pkg(3,4).eq.0.0_ezspline_r8) spline_o%ilin3=1  ! evenly spaced grid
     else
        spline_o%ilin3=1
     end if

     spline_o%isReady = 1   ! no coefficient setup necessary

  else
     !
     ! Hermite polynomial coefficient setup based on Akima's method
     ! (df/dx..etc not required)
     !
     ipx = 0
     if (spline_o%ibctype1(1)==-1 .or. spline_o%ibctype1(2)==-1) then
        ipx=1
     else if (spline_o%ibctype1(1)<-1 .or. spline_o%ibctype1(1)>1 .or. &
          spline_o%ibctype1(2)<-1 .or. spline_o%ibctype1(2)>1 ) then
        ipx=99  ! an error...
     else if (spline_o%ibctype1(1)==1 .or. spline_o%ibctype1(2)==1) then
        ipx=2
        do jj=1,spline_o%n3
           do ii=1,spline_o%n2
              if(spline_o%ibctype1(1)==1) then
                 spline_o%fspl(2,1,ii,jj)=spline_o%bcval1min(ii,jj)
              else
              ! hand implemented default BC
                 df1=(spline_o%fspl(1,2,ii,jj)-spline_o%fspl(1,1,ii,jj))/ &
                      (spline_o%x1(2)-spline_o%x1(1))
                 df2=(spline_o%fspl(1,3,ii,jj)-spline_o%fspl(1,2,ii,jj))/ &
                      (spline_o%x1(3)-spline_o%x1(2))
                 spline_o%fspl(2,1,ii,jj)=(3*df1-df2)/2
              endif
              inum=spline_o%n1
              if(spline_o%ibctype1(2)==1) then
                 spline_o%fspl(2,inum,ii,jj)=spline_o%bcval1max(ii,jj)
              else
                 ! hand implemented default BC
                 df1=(spline_o%fspl(1,inum,ii,jj)-spline_o%fspl(1,inum-1,ii,jj))/ &
                      (spline_o%x1(inum)-spline_o%x1(inum-1))
                 df2=(spline_o%fspl(1,inum-1,ii,jj)-spline_o%fspl(1,inum-2,ii,jj))/ &
                      (spline_o%x1(inum-1)-spline_o%x1(inum-2))
                 spline_o%fspl(2,inum,ii,jj)=(3*df1-df2)/2
              endif
           enddo
        enddo
     endif
     ipy = 0
     if (spline_o%ibctype2(1)==-1 .or. spline_o%ibctype2(2)==-1) then
        ipy=1
     else if (spline_o%ibctype2(1)<-1 .or. spline_o%ibctype2(1)>1 .or. &
          spline_o%ibctype2(2)<-1 .or. spline_o%ibctype2(2)>1 ) then
        ipy=99  ! an error...
     else if (spline_o%ibctype2(1)==1 .or. spline_o%ibctype2(2)==1) then
        ipy=2
        do jj=1,spline_o%n3
           do ii=1,spline_o%n1
              if(spline_o%ibctype2(1)==1) then
                 spline_o%fspl(3,ii,1,jj)=spline_o%bcval2min(ii,jj)
              else
              ! hand implemented default BC
                 df1=(spline_o%fspl(1,ii,2,jj)-spline_o%fspl(1,ii,1,jj))/ &
                      (spline_o%x2(2)-spline_o%x2(1))
                 df2=(spline_o%fspl(1,ii,3,jj)-spline_o%fspl(1,ii,2,jj))/ &
                      (spline_o%x2(3)-spline_o%x2(2))
                 spline_o%fspl(3,ii,1,jj)=(3*df1-df2)/2
              endif
              inum=spline_o%n2
              if(spline_o%ibctype2(2)==1) then
                 spline_o%fspl(3,ii,inum,jj)=spline_o%bcval2max(ii,jj)
              else
                 ! hand implemented default BC
                 df1=(spline_o%fspl(1,ii,inum,jj)-spline_o%fspl(1,ii,inum-1,jj))/ &
                      (spline_o%x2(inum)-spline_o%x2(inum-1))
                 df2=(spline_o%fspl(1,ii,inum-1,jj)-spline_o%fspl(1,ii,inum-2,jj))/ &
                      (spline_o%x2(inum-1)-spline_o%x2(inum-2))
                 spline_o%fspl(3,ii,inum,jj)=(3*df1-df2)/2
              endif
           enddo
        enddo
     endif
     ipz = 0
     if (spline_o%ibctype3(1)==-1 .or. spline_o%ibctype3(2)==-1) then
        ipz=1
     else if (spline_o%ibctype3(1)<-1 .or. spline_o%ibctype3(1)>1 .or. &
          spline_o%ibctype3(2)<-1 .or. spline_o%ibctype3(2)>1 ) then
        ipz=99  ! an error...
     else if (spline_o%ibctype3(1)==1 .or. spline_o%ibctype3(2)==1) then
        ipz=2
        do jj=1,spline_o%n2
           do ii=1,spline_o%n1
              if(spline_o%ibctype3(1)==1) then
                 spline_o%fspl(4,ii,jj,1)=spline_o%bcval3min(ii,jj)
              else
              ! hand implemented default BC
                 df1=(spline_o%fspl(1,ii,jj,2)-spline_o%fspl(1,ii,jj,1))/ &
                      (spline_o%x3(2)-spline_o%x3(1))
                 df2=(spline_o%fspl(1,ii,jj,3)-spline_o%fspl(1,ii,jj,2))/ &
                      (spline_o%x3(3)-spline_o%x3(2))
                 spline_o%fspl(4,ii,jj,1)=(3*df1-df2)/2
              endif
              inum=spline_o%n3
              if(spline_o%ibctype3(2)==1) then
                 spline_o%fspl(4,ii,jj,inum)=spline_o%bcval3max(ii,jj)
              else
                 ! hand implemented default BC
                 df1=(spline_o%fspl(1,ii,jj,inum)-spline_o%fspl(1,ii,jj,inum-1))/ &
                      (spline_o%x3(inum)-spline_o%x3(inum-1))
                 df2=(spline_o%fspl(1,ii,jj,inum-1)-spline_o%fspl(1,ii,jj,inum-2))/ &
                      (spline_o%x3(inum-1)-spline_o%x3(inum-2))
                 spline_o%fspl(4,ii,jj,inum)=(3*df1-df2)/2
              endif
           enddo
        enddo
     endif
     ifail = 0
     call r8akherm3p(spline_o%x1(1), spline_o%n1, &
          &         spline_o%x2(1), spline_o%n2, &
          &         spline_o%x3(1), spline_o%n3, &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & spline_o%ilin1, spline_o%ilin2, spline_o%ilin3, &
          & ipx,ipy,ipz, ifail)
 
     if (ifail /=0 ) then
        ier = 91
     else
        spline_o%isReady = 1
     endif
 
  endif
 
 
end subroutine EZspline_setup3_r8
!/////
! R4 !
!/////
subroutine EZspline_setup1_r4(spline_o, f, ier, exact_dim)
  use ezspline_obj
  implicit none
  type(EZspline1_r4) spline_o
  real(ezspline_r4), dimension(:), intent(in) :: f
  ! ier:
  ! 0=ok
  ! 98=some error occurred in EZspline_setup
  integer, intent(out) :: ier

  logical, intent(in), OPTIONAL :: exact_dim  ! set T to require exact
  !  dimensioning match between f and spline_o%fspl; default is F

  logical :: iexact
  integer ifail
 
  integer :: ipx
  integer iper, imsg, itol, inum, in1
  real(ezspline_r4) ztol, df1, df2
 
  !-------------------------
  iexact=.FALSE.
  if(present(exact_dim)) iexact = exact_dim
 
  if( .not.EZspline_allocated(spline_o) ) then
     ier = 98
     return
  endif
 
  in1 = size(spline_o%fspl,2)

  ier = 57
  if(size(f,1).lt.in1) return

  if(iexact) then
     ier = 58
     if(size(f,1).gt.in1) return
  endif
  ier = 0
 
  !
  ! recompute min/max in case user changed the grid manually
  spline_o%x1max = maxval(spline_o%x1)
  spline_o%x1min = minval(spline_o%x1)
 
  !
  ! set xpkg, useful for cloud and array interpolation
  imsg=0
  itol=0        ! range tolerance option
  ztol=5.e-7_ezspline_r4 ! range tolerance, if itol is set
  iper=0
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) iper=1
  call genxpkg(spline_o%n1, spline_o%x1(1), spline_o%x1pkg(1,1),&
       & iper,imsg,itol,ztol,spline_o%klookup1,ifail)
  if(ifail/=0) ier=27
 
  spline_o%isReady = 0
 
  spline_o%fspl(1, 1:in1) = &
       &  f(1:in1)
 
  if (spline_o%isHermite == 0 .and. spline_o%isLinear == 0) then
 
     call mkspline(             &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%fspl(1,1), &
          &   spline_o%ibctype1(1), spline_o%bcval1min, &
          &   spline_o%ibctype1(2), spline_o%bcval1max, &
          &   spline_o%ilin1, ifail)
 
     if(ifail /= 0) then
        ier = 98
     else
        spline_o%isReady = 1
     endif
 
  else if (spline_o%isLinear == 1) then

     spline_o%ilin1=0
     if(in1.gt.2) then
        if(spline_o%x1pkg(3,4).eq.0.0_ezspline_r4) spline_o%ilin1=1  ! evenly spaced grid
     else
        spline_o%ilin1=1
     end if

     spline_o%isReady = 1   ! no coefficient setup necessary

  else
     !
     ! Hermite polynomial coefficient setup based on Akima's method
     ! (df/dx..etc not required)
     !
     ipx = 0
     if (spline_o%ibctype1(1)==-1 .or. spline_o%ibctype1(2)==-1) then
        ipx=1
     else if (spline_o%ibctype1(1)<-1 .or. spline_o%ibctype1(1)>1 .or. &
          spline_o%ibctype1(2)<-1 .or. spline_o%ibctype1(2)>1 ) then
        ipx=99  ! an error...
     else if (spline_o%ibctype1(1)==1 .or. spline_o%ibctype1(2)==1) then
        ipx=2
        if(spline_o%ibctype1(1)==1) then
           spline_o%fspl(2,1)=spline_o%bcval1min
        else
           ! hand implemented default BC
           df1=(spline_o%fspl(1,2)-spline_o%fspl(1,1))/ &
                (spline_o%x1(2)-spline_o%x1(1))
           df2=(spline_o%fspl(1,3)-spline_o%fspl(1,2))/ &
                (spline_o%x1(3)-spline_o%x1(2))
           spline_o%fspl(2,1)=(3*df1-df2)/2
        endif
        inum=spline_o%n1
        if(spline_o%ibctype1(2)==1) then
           spline_o%fspl(2,inum)=spline_o%bcval1max
        else
           ! hand implemented default BC
           df1=(spline_o%fspl(1,inum)-spline_o%fspl(1,inum-1))/ &
                (spline_o%x1(inum)-spline_o%x1(inum-1))
           df2=(spline_o%fspl(1,inum-1)-spline_o%fspl(1,inum-2))/ &
                (spline_o%x1(inum-1)-spline_o%x1(inum-2))
           spline_o%fspl(2,inum)=(3*df1-df2)/2
        endif
     endif
     ifail = 0
     call akherm1p(spline_o%x1(1), spline_o%n1, &
          & spline_o%fspl(1,1), &
          & spline_o%ilin1, &
          & ipx, ifail)
 
     if (ifail /=0 ) then
        ier = 91
     else
        spline_o%isReady = 1
     endif
 
  endif
 
 
end subroutine EZspline_setup1_r4
 
 
subroutine EZspline_setup2_r4(spline_o, f, ier, exact_dim)
  use ezspline_obj
  implicit none
  type(EZspline2_r4) spline_o
  real(ezspline_r4), dimension(:,:), intent(in) :: f
  ! ier:
  ! 0=ok
  ! 98=some error occurred in EZspline_setup
  integer, intent(out) :: ier

  logical, intent(in), OPTIONAL :: exact_dim  ! set T to require exact
  !  dimensioning match between f and spline_o%fspl; default is F

  logical :: iexact
  integer ifail
 
  integer :: ipx, ipy
  integer iper, imsg, itol, inum, ii, jj, in0, in1, in2
  real(ezspline_r4) ztol, df1, df2
 
  !-------------------------
  iexact=.FALSE.
  if(present(exact_dim)) iexact = exact_dim
 
  if( .not.EZspline_allocated(spline_o) ) then
     ier = 98
     return
  endif
 
  in0 = size(spline_o%fspl,1)
  in1 = size(spline_o%fspl,2)
  in2 = size(spline_o%fspl,3)

  ier = 57
  if(size(f,1).lt.in1) return
  if(size(f,2).lt.in2) return

  if(iexact) then
     ier = 58
     if(size(f,1).gt.in1) return
     if(size(f,2).gt.in2) return
  endif
  ier = 0
 
  !
  ! recompute min/max in case user changed the grid manually
  spline_o%x1max = maxval(spline_o%x1)
  spline_o%x2max = maxval(spline_o%x2)
  spline_o%x1min = minval(spline_o%x1)
  spline_o%x2min = minval(spline_o%x2)
  !
  ! set xpkg, useful for cloud and array interpolation
  imsg=0
  itol=0        ! range tolerance option
  ztol=5.e-7_ezspline_r4 ! range tolerance, if itol is set
  iper=0
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) iper=1
  call genxpkg(spline_o%n1,spline_o%x1(1),spline_o%x1pkg(1,1),&
       & iper,imsg,itol,ztol,spline_o%klookup1,ifail)
  if(ifail/=0) ier=27
  iper=0
  if(spline_o%ibctype2(1)==-1 .OR. spline_o%ibctype2(2)==-1) iper=1
  call genxpkg(spline_o%n2,spline_o%x2(1),spline_o%x2pkg(1,1),&
       & iper,imsg,itol,ztol,spline_o%klookup2,ifail)
  if(ifail/=0) ier=27
 
  spline_o%isReady = 0
 
  spline_o%fspl(1, 1:in1, 1:in2) = &
       &  f(1:in1, 1:in2)
 
  if (spline_o%isHybrid == 1) then
 
     call mkintrp2d( &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%hspline, spline_o%fspl(1,1,1), &
          &   in0,in1,in2, &
          &   spline_o%ibctype1(1), spline_o%bcval1min(1), &
          &   spline_o%ibctype1(2), spline_o%bcval1max(1), &
          &   spline_o%ibctype2(1), spline_o%bcval2min(1), &
          &   spline_o%ibctype2(2), spline_o%bcval2max(1), &
          &   ifail)

     spline_o%ilin1 = 0  ! Hybrid does not compute this
     spline_o%ilin2 = 0  ! Hybrid does not compute this

     if(ifail /= 0) then
        ier = 98
     else
        spline_o%isReady = 1
     endif

  else if (spline_o%isHermite == 0 .and. spline_o%isLinear == 0) then
 
     call mkbicub(             &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%fspl(1,1,1), spline_o%n1, &
          &   spline_o%ibctype1(1), spline_o%bcval1min(1), &
          &   spline_o%ibctype1(2), spline_o%bcval1max(1), &
          &   spline_o%ibctype2(1), spline_o%bcval2min(1), &
          &   spline_o%ibctype2(2), spline_o%bcval2max(1), &
          &   spline_o%ilin1, spline_o%ilin2, ifail)
 
     if(ifail /= 0) then
        ier = 98
     else
        spline_o%isReady = 1
     endif
 
  else if (spline_o%isLinear == 1) then

     spline_o%ilin1=0
     if(in1.gt.2) then
        if(spline_o%x1pkg(3,4).eq.0.0_ezspline_r4) spline_o%ilin1=1  ! evenly spaced grid
     else
        spline_o%ilin1=1
     end if

     spline_o%ilin2=0
     if(in2.gt.2) then
        if(spline_o%x2pkg(3,4).eq.0.0_ezspline_r4) spline_o%ilin2=1  ! evenly spaced grid
     else
        spline_o%ilin2=1
     end if

     spline_o%isReady = 1   ! no coefficient setup necessary

  else
     !
     ! Hermite polynomial coefficient setup based on Akima's method
     ! (df/dx..etc not required)
     !
     ipx = 0
     if (spline_o%ibctype1(1)==-1 .or. spline_o%ibctype1(2)==-1) then
        ipx=1
     else if (spline_o%ibctype1(1)<-1 .or. spline_o%ibctype1(1)>1 .or. &
          spline_o%ibctype1(2)<-1 .or. spline_o%ibctype1(2)>1 ) then
        ipx=99  ! an error...
     else if (spline_o%ibctype1(1)==1 .or. spline_o%ibctype1(2)==1) then
        ipx=2
        do jj=1,spline_o%n2
           if(spline_o%ibctype1(1)==1) then
              spline_o%fspl(2,1,jj)=spline_o%bcval1min(jj)
           else
              ! hand implemented default BC
              df1=(spline_o%fspl(1,2,jj)-spline_o%fspl(1,1,jj))/ &
                   (spline_o%x1(2)-spline_o%x1(1))
              df2=(spline_o%fspl(1,3,jj)-spline_o%fspl(1,2,jj))/ &
                   (spline_o%x1(3)-spline_o%x1(2))
              spline_o%fspl(2,1,jj)=(3*df1-df2)/2
           endif
           inum=spline_o%n1
           if(spline_o%ibctype1(2)==1) then
              spline_o%fspl(2,inum,jj)=spline_o%bcval1max(jj)
           else
              ! hand implemented default BC
              df1=(spline_o%fspl(1,inum,jj)-spline_o%fspl(1,inum-1,jj))/ &
                   (spline_o%x1(inum)-spline_o%x1(inum-1))
              df2=(spline_o%fspl(1,inum-1,jj)-spline_o%fspl(1,inum-2,jj))/ &
                   (spline_o%x1(inum-1)-spline_o%x1(inum-2))
              spline_o%fspl(2,inum,jj)=(3*df1-df2)/2
           endif
        enddo
     endif
     ipy = 0
     if (spline_o%ibctype2(1)==-1 .or. spline_o%ibctype2(2)==-1) then
        ipy=1
     else if (spline_o%ibctype2(1)<-1 .or. spline_o%ibctype2(1)>1 .or. &
          spline_o%ibctype2(2)<-1 .or. spline_o%ibctype2(2)>1 ) then
        ipy=99  ! an error...
     else if (spline_o%ibctype2(1)==1 .or. spline_o%ibctype2(2)==1) then
        ipy=2
        do ii=1,spline_o%n1
           if(spline_o%ibctype2(1)==1) then
              spline_o%fspl(3,ii,1)=spline_o%bcval2min(ii)
           else
              ! hand implemented default BC
              df1=(spline_o%fspl(1,ii,2)-spline_o%fspl(1,ii,1))/ &
                   (spline_o%x2(2)-spline_o%x2(1))
              df2=(spline_o%fspl(1,ii,3)-spline_o%fspl(1,ii,2))/ &
                   (spline_o%x2(3)-spline_o%x2(2))
              spline_o%fspl(3,ii,1)=(3*df1-df2)/2
           endif
           inum=spline_o%n2
           if(spline_o%ibctype2(2)==1) then
              spline_o%fspl(3,ii,inum)=spline_o%bcval2max(ii)
           else
              ! hand implemented default BC
              df1=(spline_o%fspl(1,ii,inum)-spline_o%fspl(1,ii,inum-1))/ &
                   (spline_o%x2(inum)-spline_o%x2(inum-1))
              df2=(spline_o%fspl(1,ii,inum-1)-spline_o%fspl(1,ii,inum-2))/ &
                   (spline_o%x2(inum-1)-spline_o%x2(inum-2))
              spline_o%fspl(3,ii,inum)=(3*df1-df2)/2
           endif
        enddo
     endif
     ifail = 0
     call akherm2p(spline_o%x1(1), spline_o%n1, &
          &         spline_o%x2(1), spline_o%n2, &
          & spline_o%fspl(1,1,1), spline_o%n1, &
          & spline_o%ilin1, spline_o%ilin2, &
          & ipx,ipy, ifail)
 
     if (ifail /=0 ) then
        ier = 91
     else
        spline_o%isReady = 1
     endif
 
  endif
 
 
end subroutine EZspline_setup2_r4
 
 
 
subroutine EZspline_setup3_r4(spline_o, f, ier, exact_dim)
  use ezspline_obj
  implicit none
  type(EZspline3_r4) spline_o
  real(ezspline_r4), dimension(:,:,:), intent(in) :: f
  ! ier:
  ! 0=ok
  ! 98=some error occurred in EZspline_setup
  integer, intent(out) :: ier

  logical, intent(in), OPTIONAL :: exact_dim  ! set T to require exact
  !  dimensioning match between f and spline_o%fspl; default is F

  logical :: iexact
  integer ifail
 
  integer :: ipx, ipy, ipz
  integer iper, imsg, itol, inum, ii, jj, in0, in1, in2, in3
  real(ezspline_r4) ztol, df1, df2
 
  !-------------------------
  iexact=.FALSE.
  if(present(exact_dim)) iexact = exact_dim
 
  if( .not.EZspline_allocated(spline_o) ) then
     ier = 98
     return
  endif
 
  in0 = size(spline_o%fspl,1)
  in1 = size(spline_o%fspl,2)
  in2 = size(spline_o%fspl,3)
  in3 = size(spline_o%fspl,4)

  ier = 57
  if(size(f,1).lt.in1) return
  if(size(f,2).lt.in2) return
  if(size(f,3).lt.in3) return

  if(iexact) then
     ier = 58
     if(size(f,1).gt.in1) return
     if(size(f,2).gt.in2) return
     if(size(f,3).gt.in3) return
  endif
  ier = 0
 
  !
  ! recompute min/max in case user changed the grid manually
  spline_o%x1max = maxval(spline_o%x1)
  spline_o%x2max = maxval(spline_o%x2)
  spline_o%x3max = maxval(spline_o%x3)
  spline_o%x1min = minval(spline_o%x1)
  spline_o%x2min = minval(spline_o%x2)
  spline_o%x3min = minval(spline_o%x3)
 
  !
  ! set xpkg, useful for cloud and array interpolation
  imsg=0
  itol=0        ! range tolerance option
  ztol=5.e-7_ezspline_r4 ! range tolerance, if itol is set
  iper=0
  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) iper=1
  call genxpkg(spline_o%n1,spline_o%x1(1),spline_o%x1pkg(1,1),&
       & iper,imsg,itol,ztol,spline_o%klookup1,ifail)
  if(ifail/=0) ier=27
  iper=0
  if(spline_o%ibctype2(1)==-1 .OR. spline_o%ibctype2(2)==-1) iper=1
  call genxpkg(spline_o%n2,spline_o%x2(1),spline_o%x2pkg(1,1),&
       & iper,imsg,itol,ztol,spline_o%klookup2,ifail)
  if(ifail/=0) ier=27
  iper=0
  if(spline_o%ibctype3(1)==-1 .OR. spline_o%ibctype3(2)==-1) iper=1
  call genxpkg(spline_o%n3,spline_o%x3(1),spline_o%x3pkg(1,1),&
       & iper,imsg,itol,ztol,spline_o%klookup3,ifail)
  if(ifail/=0) ier=27
 
  spline_o%isReady = 0
 
  spline_o%fspl(1, 1:in1, 1:in2, 1:in3) = &
       &  f(1:in1, 1:in2, 1:in3)
 
  if (spline_o%isHybrid == 1) then
 
     call mkintrp3d( &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%x3(1), spline_o%n3, &
          &   spline_o%hspline, spline_o%fspl(1,1,1,1), &
          &   in0,in1,in2,in3, &
          &   spline_o%ibctype1(1), spline_o%bcval1min(1,1), &
          &   spline_o%ibctype1(2), spline_o%bcval1max(1,1), &
          &   spline_o%ibctype2(1), spline_o%bcval2min(1,1), &
          &   spline_o%ibctype2(2), spline_o%bcval2max(1,1), &
          &   spline_o%ibctype3(1), spline_o%bcval3min(1,1), &
          &   spline_o%ibctype3(2), spline_o%bcval3max(1,1), &
          &   ifail)

     spline_o%ilin1 = 0  ! Hybrid does not compute this
     spline_o%ilin2 = 0  ! Hybrid does not compute this
     spline_o%ilin3 = 0  ! Hybrid does not compute this

     if(ifail /= 0) then
        ier = 98
     else
        spline_o%isReady = 1
     endif

  else if (spline_o%isHermite == 0 .and. spline_o%isLinear == 0) then
 
     call mktricub(             &
          &   spline_o%x1(1), spline_o%n1, &
          &   spline_o%x2(1), spline_o%n2, &
          &   spline_o%x3(1), spline_o%n3, &
          &   spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          &   spline_o%ibctype1(1), spline_o%bcval1min(1,1), &
          &   spline_o%ibctype1(2), spline_o%bcval1max(1,1), spline_o%n2, &
          &   spline_o%ibctype2(1), spline_o%bcval2min(1,1), &
          &   spline_o%ibctype2(2), spline_o%bcval2max(1,1), spline_o%n1, &
          &   spline_o%ibctype3(1), spline_o%bcval3min(1,1), &
          &   spline_o%ibctype3(2), spline_o%bcval3max(1,1), spline_o%n1, &
          &   spline_o%ilin1, spline_o%ilin2, spline_o%ilin3, ifail)
 
     if(ifail /= 0) then
        ier = 98
     else
        spline_o%isReady = 1
     endif
 
  else if (spline_o%isLinear == 1) then

     spline_o%ilin1=0
     if(in1.gt.2) then
        if(spline_o%x1pkg(3,4).eq.0.0_ezspline_r4) spline_o%ilin1=1  ! evenly spaced grid
     else
        spline_o%ilin1=1
     end if

     spline_o%ilin2=0
     if(in2.gt.2) then
        if(spline_o%x2pkg(3,4).eq.0.0_ezspline_r4) spline_o%ilin2=1  ! evenly spaced grid
     else
        spline_o%ilin2=1
     end if

     spline_o%ilin3=0
     if(in3.gt.2) then
        if(spline_o%x3pkg(3,4).eq.0.0_ezspline_r4) spline_o%ilin3=1  ! evenly spaced grid
     else
        spline_o%ilin3=1
     end if

     spline_o%isReady = 1   ! no coefficient setup necessary

  else
     !
     ! Hermite polynomial coefficient setup based on Akima's method
     ! (df/dx..etc not required)
     !
     ipx = 0
     if (spline_o%ibctype1(1)==-1 .or. spline_o%ibctype1(2)==-1) then
        ipx=1
     else if (spline_o%ibctype1(1)<-1 .or. spline_o%ibctype1(1)>1 .or. &
          spline_o%ibctype1(2)<-1 .or. spline_o%ibctype1(2)>1 ) then
        ipx=99  ! an error...
     else if (spline_o%ibctype1(1)==1 .or. spline_o%ibctype1(2)==1) then
        ipx=2
        do jj=1,spline_o%n3
           do ii=1,spline_o%n2
              if(spline_o%ibctype1(1)==1) then
                 spline_o%fspl(2,1,ii,jj)=spline_o%bcval1min(ii,jj)
              else
              ! hand implemented default BC
                 df1=(spline_o%fspl(1,2,ii,jj)-spline_o%fspl(1,1,ii,jj))/ &
                      (spline_o%x1(2)-spline_o%x1(1))
                 df2=(spline_o%fspl(1,3,ii,jj)-spline_o%fspl(1,2,ii,jj))/ &
                      (spline_o%x1(3)-spline_o%x1(2))
                 spline_o%fspl(2,1,ii,jj)=(3*df1-df2)/2
              endif
              inum=spline_o%n1
              if(spline_o%ibctype1(2)==1) then
                 spline_o%fspl(2,inum,ii,jj)=spline_o%bcval1max(ii,jj)
              else
                 ! hand implemented default BC
                 df1=(spline_o%fspl(1,inum,ii,jj)-spline_o%fspl(1,inum-1,ii,jj))/ &
                      (spline_o%x1(inum)-spline_o%x1(inum-1))
                 df2=(spline_o%fspl(1,inum-1,ii,jj)-spline_o%fspl(1,inum-2,ii,jj))/ &
                      (spline_o%x1(inum-1)-spline_o%x1(inum-2))
                 spline_o%fspl(2,inum,ii,jj)=(3*df1-df2)/2
              endif
           enddo
        enddo
     endif
     ipy = 0
     if (spline_o%ibctype2(1)==-1 .or. spline_o%ibctype2(2)==-1) then
        ipy=1
     else if (spline_o%ibctype2(1)<-1 .or. spline_o%ibctype2(1)>1 .or. &
          spline_o%ibctype2(2)<-1 .or. spline_o%ibctype2(2)>1 ) then
        ipy=99  ! an error...
     else if (spline_o%ibctype2(1)==1 .or. spline_o%ibctype2(2)==1) then
        ipy=2
        do jj=1,spline_o%n3
           do ii=1,spline_o%n1
              if(spline_o%ibctype2(1)==1) then
                 spline_o%fspl(3,ii,1,jj)=spline_o%bcval2min(ii,jj)
              else
              ! hand implemented default BC
                 df1=(spline_o%fspl(1,ii,2,jj)-spline_o%fspl(1,ii,1,jj))/ &
                      (spline_o%x2(2)-spline_o%x2(1))
                 df2=(spline_o%fspl(1,ii,3,jj)-spline_o%fspl(1,ii,2,jj))/ &
                      (spline_o%x2(3)-spline_o%x2(2))
                 spline_o%fspl(3,ii,1,jj)=(3*df1-df2)/2
              endif
              inum=spline_o%n2
              if(spline_o%ibctype2(2)==1) then
                 spline_o%fspl(3,ii,inum,jj)=spline_o%bcval2max(ii,jj)
              else
                 ! hand implemented default BC
                 df1=(spline_o%fspl(1,ii,inum,jj)-spline_o%fspl(1,ii,inum-1,jj))/ &
                      (spline_o%x2(inum)-spline_o%x2(inum-1))
                 df2=(spline_o%fspl(1,ii,inum-1,jj)-spline_o%fspl(1,ii,inum-2,jj))/ &
                      (spline_o%x2(inum-1)-spline_o%x2(inum-2))
                 spline_o%fspl(3,ii,inum,jj)=(3*df1-df2)/2
              endif
           enddo
        enddo
     endif
     ipz = 0
     if (spline_o%ibctype3(1)==-1 .or. spline_o%ibctype3(2)==-1) then
        ipz=1
     else if (spline_o%ibctype3(1)<-1 .or. spline_o%ibctype3(1)>1 .or. &
          spline_o%ibctype3(2)<-1 .or. spline_o%ibctype3(2)>1 ) then
        ipz=99  ! an error...
     else if (spline_o%ibctype3(1)==1 .or. spline_o%ibctype3(2)==1) then
        ipz=2
        do jj=1,spline_o%n2
           do ii=1,spline_o%n1
              if(spline_o%ibctype3(1)==1) then
                 spline_o%fspl(4,ii,jj,1)=spline_o%bcval3min(ii,jj)
              else
              ! hand implemented default BC
                 df1=(spline_o%fspl(1,ii,jj,2)-spline_o%fspl(1,ii,jj,1))/ &
                      (spline_o%x3(2)-spline_o%x3(1))
                 df2=(spline_o%fspl(1,ii,jj,3)-spline_o%fspl(1,ii,jj,2))/ &
                      (spline_o%x3(3)-spline_o%x3(2))
                 spline_o%fspl(4,ii,jj,1)=(3*df1-df2)/2
              endif
              inum=spline_o%n3
              if(spline_o%ibctype3(2)==1) then
                 spline_o%fspl(4,ii,jj,inum)=spline_o%bcval3max(ii,jj)
              else
                 ! hand implemented default BC
                 df1=(spline_o%fspl(1,ii,jj,inum)-spline_o%fspl(1,ii,jj,inum-1))/ &
                      (spline_o%x3(inum)-spline_o%x3(inum-1))
                 df2=(spline_o%fspl(1,ii,jj,inum-1)-spline_o%fspl(1,ii,jj,inum-2))/ &
                      (spline_o%x3(inum-1)-spline_o%x3(inum-2))
                 spline_o%fspl(4,ii,jj,inum)=(3*df1-df2)/2
              endif
           enddo
        enddo
     endif
     ifail = 0
     call akherm3p(spline_o%x1(1), spline_o%n1, &
          &         spline_o%x2(1), spline_o%n2, &
          &         spline_o%x3(1), spline_o%n3, &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & spline_o%ilin1, spline_o%ilin2, spline_o%ilin3, &
          & ipx,ipy,ipz, ifail)
 
     if (ifail /=0 ) then
        ier = 91
     else
        spline_o%isReady = 1
     endif
 
  endif
 
 
end subroutine EZspline_setup3_r4
