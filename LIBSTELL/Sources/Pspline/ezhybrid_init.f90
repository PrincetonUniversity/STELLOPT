!/////
! R8 !
!/////
 
subroutine EZhybrid_init2_r8(spline_o, n1, n2, hspline, ier, &
     BCS1, BCS2)
  use ezspline_obj
  implicit none
  type(EZspline2_r8) spline_o
  integer, intent(in) :: n1, n2
  integer, intent(in) :: hspline(2)
  ! ier:
  ! 0=ok
  ! 1=allocation error
  ! 99=something strange happened in EZhybrid_init
  integer, intent(out) :: ier

  integer, intent(in), OPTIONAL :: BCS1(2), BCS2(2)

  integer i, iok, icoeff, idim1, idim2
  logical :: spline_present, hermite_present
 
  ier = 0
 
  if(EZspline_allocated(spline_o)) then
     ier = 100  ! allocated already
     return
  else
     call EZspline_preInit(spline_o)
  endif

  spline_o%n1 = n1
  spline_o%n2 = n2

  spline_o%klookup1 = 3    ! default lookup algorithm selection
  spline_o%klookup2 = 3    ! default lookup algorithm selection
 
  spline_o%isHermite = 0
  spline_o%isLinear = 0
  spline_o%isHybrid = 1
  spline_o%hspline = hspline

  spline_present=.FALSE.
  hermite_present=.FALSE.
  icoeff=1
  do i=1,2
     if(hspline(i).lt.-1) then
        ier=54
     else if(hspline(i).gt.2) then
        ier=54
     else if(hspline(i).eq.1) then
        hermite_present=.TRUE.
        icoeff=icoeff*2
     else if(hspline(i).eq.2) then
        spline_present=.TRUE.
        icoeff=icoeff*2
     endif
  enddo
  if((ier.eq.0).and.spline_present.and.hermite_present) ier=55
  if(ier.ne.0) return
 
  idim1 = n1
  if(hspline(1).eq.-1) idim1 = n1 - 1
  idim2 = n2
  if(hspline(2).eq.-1) idim2 = n2 - 1

  iok = 0
  allocate(spline_o%x1(n1), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x2(n2), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%fspl(icoeff,idim1,idim2), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval1min(idim2), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval1max(idim2), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval2min(idim1), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval2max(idim1), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x1pkg(n1, 4), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x2pkg(n2, 4), stat=iok)
  if(iok /= 0) ier = 1

  do i = 1, 2
 
     spline_o%bcval1min(1:idim2) = 0.0_ezspline_r8
     spline_o%bcval1max(1:idim2) = 0.0_ezspline_r8
     spline_o%ibctype1(i) = 0
     if(present(BCS1)) then
        select case(BCS1(i))
        case (-1)
           spline_o%ibctype1(i) = -1
        case (0)
           spline_o%ibctype1(i) = 0
        case (1)
           if(hspline(1).le.0) then
              spline_o%ibctype1(i) = 0
              ier = 56
           else
              spline_o%ibctype1(i) = 1
           endif
        case (2)
           if(hspline(1).le.1) then
              spline_o%ibctype1(i) = 0
              ier = 56
           else
              spline_o%ibctype1(i) = 2
           endif
        case default
           ier = 2
           spline_o%ibctype1(i) = 0
        end select
     endif
 
     spline_o%bcval2min(1:idim1) = 0.0_ezspline_r8
     spline_o%bcval2max(1:idim1) = 0.0_ezspline_r8
     spline_o%ibctype2(i) = 0
     if(present(BCS2)) then
        select case(BCS2(i))
        case (-1)
           spline_o%ibctype2(i) = -1
        case (0)
           spline_o%ibctype2(i) = 0
        case (1)
           if(hspline(2).le.0) then
              spline_o%ibctype2(i) = 0
              ier = 56
           else
              spline_o%ibctype2(i) = 1
           endif
        case (2)
           if(hspline(2).le.1) then
              spline_o%ibctype2(i) = 0
              ier = 56
           else
              spline_o%ibctype2(i) = 2
           endif
        case default
           ier = 2
           spline_o%ibctype2(i) = 0
        end select
     endif
 
  enddo

  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     spline_o%ibctype1(1)=-1
     spline_o%ibctype1(2)=-1
  endif
  if(spline_o%ibctype2(1)==-1 .OR. spline_o%ibctype2(2)==-1) then
     spline_o%ibctype2(1)=-1
     spline_o%ibctype2(2)=-1
  endif
 
  ! default grid
  spline_o%x1min = 0.0_ezspline_r8
  spline_o%x1max = 1.0_ezspline_r8
  if(spline_o%ibctype1(1).eq.-1) spline_o%x1max = ezspline_twopi_r8
 
  spline_o%x1 = spline_o%x1min + (spline_o%x1max - spline_o%x1min)* &
       & (/ (real(i-1,ezspline_r8)/real(spline_o%n1-1, ezspline_r8), i=1,spline_o%n1 ) /)
  spline_o%x2min = 0.0_ezspline_r8
  spline_o%x2max = 1.0_ezspline_r8
  if(spline_o%ibctype2(1).eq.-1) spline_o%x2max = ezspline_twopi_r8
  spline_o%x2 = spline_o%x2min + (spline_o%x2max - spline_o%x2min)* &
       & (/ (real(i-1,ezspline_r8)/real(spline_o%n2-1, ezspline_r8), i=1,spline_o%n2 ) /)
 
  spline_o%isReady = 0
 
  return
end subroutine EZhybrid_init2_r8
 
 
 
 
subroutine EZhybrid_init3_r8(spline_o, n1, n2, n3, hspline, ier, &
     BCS1, BCS2, BCS3)
  use ezspline_obj
  implicit none
  type(EZspline3_r8) spline_o
  integer, intent(in) :: n1, n2, n3
  integer, intent(in) :: hspline(3)
  ! ier:
  ! 0=ok
  ! 1=allocation error
  ! 99=something strange happened in EZhybrid_init
  integer, intent(out) :: ier

  integer, intent(in), OPTIONAL :: BCS1(2), BCS2(2), BCS3(2)

  integer i, iok, icoeff, idim1, idim2, idim3
  logical :: spline_present, hermite_present
 
  ier = 0
 
  if(EZspline_allocated(spline_o)) then
     ier = 100  ! allocated already
     return
  else
     call EZspline_preInit(spline_o)
  endif

  spline_o%n1 = n1
  spline_o%n2 = n2
  spline_o%n3 = n3
 
  spline_o%klookup1 = 3    ! default lookup algorithm selection
  spline_o%klookup2 = 3    ! default lookup algorithm selection
  spline_o%klookup3 = 3    ! default lookup algorithm selection
  
  spline_o%isHermite = 0
  spline_o%isLinear = 0
  spline_o%isHybrid = 1
  spline_o%hspline = hspline

  spline_present=.FALSE.
  hermite_present=.FALSE.
  icoeff=1
  do i=1,3
     if(hspline(i).lt.-1) then
        ier=54
     else if(hspline(i).gt.2) then
        ier=54
     else if(hspline(i).eq.1) then
        hermite_present=.TRUE.
        icoeff=icoeff*2
     else if(hspline(i).eq.2) then
        spline_present=.TRUE.
        icoeff=icoeff*2
     endif
  enddo
  if((ier.eq.0).and.spline_present.and.hermite_present) ier=55
  if(ier.ne.0) return
  
  idim1 = n1
  if(hspline(1).eq.-1) idim1 = n1 - 1
  idim2 = n2
  if(hspline(2).eq.-1) idim2 = n2 - 1
  idim3 = n3
  if(hspline(3).eq.-1) idim3 = n3 - 1

  iok = 0
  allocate(spline_o%x1(n1), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x2(n2), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x3(n3), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%fspl(icoeff,idim1,idim2,idim3), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval1min(idim2, idim3), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval1max(idim2, idim3), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval2min(idim1, idim3), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval2max(idim1, idim3), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval3min(idim1, idim2), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval3max(idim1, idim2), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x1pkg(n1, 4), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x2pkg(n2, 4), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x3pkg(n3, 4), stat=iok)
  if(iok /= 0) ier = 1

  do i = 1, 2
 
     spline_o%bcval1min(1:idim2, 1:idim3) = 0.0_ezspline_r8
     spline_o%bcval1max(1:idim2, 1:idim3) = 0.0_ezspline_r8
     spline_o%ibctype1(i) = 0
     if(present(BCS1)) then
        select case(BCS1(i))
        case (-1)
           spline_o%ibctype1(i) = -1
        case (0)
           spline_o%ibctype1(i) = 0
        case (1)
           if(hspline(1).le.0) then
              spline_o%ibctype1(i) = 0
              ier = 56
           else
              spline_o%ibctype1(i) = 1
           endif
        case (2)
           if(hspline(1).le.1) then
              spline_o%ibctype1(i) = 0
              ier = 56
           else
              spline_o%ibctype1(i) = 2
           endif
        case default
           ier = 2
           spline_o%ibctype1(i) = 0
        end select
     endif
 
     spline_o%bcval2min(1:idim1, 1:idim3) = 0.0_ezspline_r8
     spline_o%bcval2max(1:idim1, 1:idim3) = 0.0_ezspline_r8
     spline_o%ibctype2(i) = 0
     if(present(BCS2)) then
        select case(BCS2(i))
        case (-1)
           spline_o%ibctype2(i) = -1
        case (0)
           spline_o%ibctype2(i) = 0
        case (1)
           if(hspline(2).le.0) then
              spline_o%ibctype2(i) = 0
              ier = 56
           else
              spline_o%ibctype2(i) = 1
           endif
        case (2)
           if(hspline(2).le.1) then
              spline_o%ibctype2(i) = 0
              ier = 56
           else
              spline_o%ibctype2(i) = 2
           endif
        case default
           ier = 2
           spline_o%ibctype2(i) = 0
        end select
     endif
 
     spline_o%bcval3min(1:idim1, 1:idim2) = 0.0_ezspline_r8
     spline_o%bcval3max(1:idim1, 1:idim2) = 0.0_ezspline_r8
     spline_o%ibctype3(i) = 0
     if(present(BCS3)) then
        select case(BCS3(i))
        case (-1)
           spline_o%ibctype3(i) = -1
        case (0)
           spline_o%ibctype3(i) = 0
        case (1)
           if(hspline(3).le.0) then
              spline_o%ibctype3(i) = 0
              ier = 2
           else
              spline_o%ibctype3(i) = 1
           endif
        case (2)
           if(hspline(3).le.1) then
              spline_o%ibctype3(i) = 0
              ier = 2
           else
              spline_o%ibctype3(i) = 2
           endif
        case default
           ier = 2
           spline_o%ibctype3(i) = 0
        end select
     endif
 
  enddo

  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     spline_o%ibctype1(1)=-1
     spline_o%ibctype1(2)=-1
  endif
  if(spline_o%ibctype2(1)==-1 .OR. spline_o%ibctype2(2)==-1) then
     spline_o%ibctype2(1)=-1
     spline_o%ibctype2(2)=-1
  endif
  if(spline_o%ibctype3(1)==-1 .OR. spline_o%ibctype3(2)==-1) then
     spline_o%ibctype3(1)=-1
     spline_o%ibctype3(2)=-1
  endif
 
  ! default grid
  spline_o%x1min = 0.0_ezspline_r8
  spline_o%x1max = 1.0_ezspline_r8
  if(spline_o%ibctype1(1).eq.-1) spline_o%x1max = ezspline_twopi_r8
  spline_o%x1 = spline_o%x1min + (spline_o%x1max - spline_o%x1min)* &
       & (/ (real(i-1,ezspline_r8)/real(spline_o%n1-1, ezspline_r8), i=1,spline_o%n1 ) /)

  spline_o%x2min = 0.0_ezspline_r8
  spline_o%x2max = 1.0_ezspline_r8
  if(spline_o%ibctype2(1).eq.-1) spline_o%x2max = ezspline_twopi_r8
  spline_o%x2 = spline_o%x2min + (spline_o%x2max - spline_o%x2min)* &
       & (/ (real(i-1,ezspline_r8)/real(spline_o%n2-1, ezspline_r8), i=1,spline_o%n2 ) /)
 
  spline_o%x3min = 0.0_ezspline_r8
  spline_o%x3max = 1.0_ezspline_r8
  if(spline_o%ibctype3(1).eq.-1) spline_o%x3max = ezspline_twopi_r8
  spline_o%x3 = spline_o%x3min + (spline_o%x3max - spline_o%x3min)* &
       & (/ (real(i-1,ezspline_r8)/real(spline_o%n3-1, ezspline_r8), i=1,spline_o%n3 ) /)
 
  spline_o%isReady = 0
 
  return
end subroutine EZhybrid_init3_r8
 
!/////
! R4 !
!/////
 
subroutine EZhybrid_init2_r4(spline_o, n1, n2, hspline, ier, &
     BCS1, BCS2)
  use ezspline_obj
  implicit none
  type(EZspline2_r4) spline_o
  integer, intent(in) :: n1, n2
  integer, intent(in) :: hspline(2)
  ! ier:
  ! 0=ok
  ! 1=allocation error
  ! 99=something strange happened in EZhybrid_init
  integer, intent(out) :: ier

  integer, intent(in), OPTIONAL :: BCS1(2), BCS2(2)

  integer i, iok, icoeff, idim1, idim2
  logical :: spline_present, hermite_present
 
  ier = 0
 
  if(EZspline_allocated(spline_o)) then
     ier = 100  ! allocated already
     return
  else
     call EZspline_preInit(spline_o)
  endif
 
  spline_o%n1 = n1
  spline_o%n2 = n2
 
  spline_o%klookup1 = 3    ! default lookup algorithm selection
  spline_o%klookup2 = 3    ! default lookup algorithm selection
 
  spline_o%isHermite = 0
  spline_o%isLinear = 0
  spline_o%isHybrid = 1
  spline_o%hspline = hspline

  spline_present=.FALSE.
  hermite_present=.FALSE.
  icoeff=1
  do i=1,2
     if(hspline(i).lt.-1) then
        ier=54
     else if(hspline(i).gt.2) then
        ier=54
     else if(hspline(i).eq.1) then
        hermite_present=.TRUE.
        icoeff=icoeff*2
     else if(hspline(i).eq.2) then
        spline_present=.TRUE.
        icoeff=icoeff*2
     endif
  enddo
  if((ier.eq.0).and.spline_present.and.hermite_present) ier=55
  if(ier.ne.0) return
 
  idim1 = n1
  if(hspline(1).eq.-1) idim1 = n1 - 1
  idim2 = n2
  if(hspline(2).eq.-1) idim2 = n2 - 1
 
  iok = 0
  allocate(spline_o%x1(n1), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x2(n2), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%fspl(icoeff,idim1,idim2), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval1min(idim2), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval1max(idim2), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval2min(idim1), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval2max(idim1), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x1pkg(n1, 4), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x2pkg(n2, 4), stat=iok)
  if(iok /= 0) ier = 1

  do i = 1, 2
 
     spline_o%bcval1min(1:idim2) = 0.0
     spline_o%bcval1max(1:idim2) = 0.0
     spline_o%ibctype1(i) = 0
     if(present(BCS1)) then
        select case(BCS1(i))
        case (-1)
           spline_o%ibctype1(i) = -1
        case (0)
           spline_o%ibctype1(i) = 0
        case (1)
           if(hspline(1).le.0) then
              spline_o%ibctype1(i) = 0
              ier = 56
           else
              spline_o%ibctype1(i) = 1
           endif
        case (2)
           if(hspline(1).le.1) then
              spline_o%ibctype1(i) = 0
              ier = 56
           else
              spline_o%ibctype1(i) = 2
           endif
        case default
           ier = 2
           spline_o%ibctype1(i) = 0
        end select
     endif
 
     spline_o%bcval2min(1:idim1) = 0.0
     spline_o%bcval2max(1:idim1) = 0.0
     spline_o%ibctype2(i) = 0
     if(present(BCS2)) then
        select case(BCS2(i))
        case (-1)
           spline_o%ibctype2(i) = -1
        case (0)
           spline_o%ibctype2(i) = 0
        case (1)
           if(hspline(2).le.0) then
              spline_o%ibctype2(i) = 0
              ier = 56
           else
              spline_o%ibctype2(i) = 1
           endif
        case (2)
           if(hspline(2).le.1) then
              spline_o%ibctype2(i) = 0
              ier = 56
           else
              spline_o%ibctype2(i) = 2
           endif
        case default
           ier = 2
           spline_o%ibctype2(i) = 0
        end select
     endif
 
  enddo

  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     spline_o%ibctype1(1)=-1
     spline_o%ibctype1(2)=-1
  endif
  if(spline_o%ibctype2(1)==-1 .OR. spline_o%ibctype2(2)==-1) then
     spline_o%ibctype2(1)=-1
     spline_o%ibctype2(2)=-1
  endif
 
  ! default grid
  spline_o%x1min = 0.0_ezspline_r4
  spline_o%x1max = 1.0_ezspline_r4
  if(spline_o%ibctype1(1).eq.-1) spline_o%x1max = ezspline_twopi_r4
  spline_o%x1 = spline_o%x1min + (spline_o%x1max - spline_o%x1min)* &
       & (/ (real(i-1,ezspline_r4)/real(spline_o%n1-1, ezspline_r4), i=1,spline_o%n1 ) /)

  spline_o%x2min = 0.0_ezspline_r4
  spline_o%x2max = 1.0_ezspline_r4
  if(spline_o%ibctype2(1).eq.-1) spline_o%x2max = ezspline_twopi_r4
  spline_o%x2 = spline_o%x2min + (spline_o%x2max - spline_o%x2min)* &
       & (/ (real(i-1,ezspline_r4)/real(spline_o%n2-1, ezspline_r4), i=1,spline_o%n2 ) /)
 
 
  spline_o%isReady = 0
 
  return
end subroutine EZhybrid_init2_r4
 
 
 
 
subroutine EZhybrid_init3_r4(spline_o, n1, n2, n3, hspline, ier, &
     BCS1, BCS2, BCS3)
  use ezspline_obj
  implicit none
  type(EZspline3_r4) spline_o
  integer, intent(in) :: n1, n2, n3
  integer, intent(in) :: hspline(3)
  ! ier:
  ! 0=ok
  ! 1=allocation error
  ! 99=something strange happened in EZhybrid_init
  integer, intent(out) :: ier

  integer, intent(in), OPTIONAL :: BCS1(2), BCS2(2), BCS3(2)

  integer i, iok, icoeff, idim1, idim2, idim3
  logical :: spline_present, hermite_present
 
  ier = 0
 
  if(EZspline_allocated(spline_o)) then
     ier = 100  ! allocated already
     return
  else
     call EZspline_preInit(spline_o)
  endif
 
  spline_o%n1 = n1
  spline_o%n2 = n2
  spline_o%n3 = n3
 
  spline_o%klookup1 = 3    ! default lookup algorithm selection
  spline_o%klookup2 = 3    ! default lookup algorithm selection
  spline_o%klookup3 = 3    ! default lookup algorithm selection
 
  spline_o%isHermite = 0
  spline_o%isLinear = 0
  spline_o%isHybrid = 1
  spline_o%hspline = hspline

  spline_present=.FALSE.
  hermite_present=.FALSE.
  icoeff=1
  do i=1,3
     if(hspline(i).lt.-1) then
        ier=54
     else if(hspline(i).gt.2) then
        ier=54
     else if(hspline(i).eq.1) then
        hermite_present=.TRUE.
        icoeff=icoeff*2
     else if(hspline(i).eq.2) then
        spline_present=.TRUE.
        icoeff=icoeff*2
     endif
  enddo
  if((ier.eq.0).and.spline_present.and.hermite_present) ier=55
  if(ier.ne.0) return
 
  idim1 = n1
  if(hspline(1).eq.-1) idim1 = n1 - 1
  idim2 = n2
  if(hspline(2).eq.-1) idim2 = n2 - 1
  idim3 = n3
  if(hspline(3).eq.-1) idim3 = n3 - 1
 
  iok = 0
  allocate(spline_o%x1(n1), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x2(n2), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x3(n3), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%fspl(icoeff,idim1,idim2,idim3), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval1min(idim2, idim3), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval1max(idim2, idim3), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval2min(idim1, idim3), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval2max(idim1, idim3), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval3min(idim1, idim2), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%bcval3max(idim1, idim2), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x1pkg(n1, 4), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x2pkg(n2, 4), stat=iok)
  if(iok /= 0) ier = 1
  allocate(spline_o%x3pkg(n3, 4), stat=iok)
  if(iok /= 0) ier = 1

  do i = 1, 2
 
     spline_o%bcval1min(1:idim2, 1:idim3) = 0.0
     spline_o%bcval1max(1:idim2, 1:idim3) = 0.0
     spline_o%ibctype1(i) = 0
     if(present(BCS1)) then
        select case(BCS1(i))
        case (-1)
           spline_o%ibctype1(i) = -1
        case (0)
           spline_o%ibctype1(i) = 0
        case (1)
           if(hspline(1).le.0) then
              spline_o%ibctype1(i) = 0
              ier = 56
           else
              spline_o%ibctype1(i) = 1
           endif
        case (2)
           if(hspline(1).le.1) then
              spline_o%ibctype1(i) = 0
              ier = 56
           else
              spline_o%ibctype1(i) = 2
           endif
        case default
           ier = 2
           spline_o%ibctype1(i) = 0
        end select
     endif
 
     spline_o%bcval2min(1:idim1, 1:idim3) = 0.0
     spline_o%bcval2max(1:idim1, 1:idim3) = 0.0
     spline_o%ibctype2(i) = 0
     if(present(BCS2)) then
        select case(BCS2(i))
        case (-1)
           spline_o%ibctype2(i) = -1
        case (0)
           spline_o%ibctype2(i) = 0
        case (1)
           if(hspline(2).le.0) then
              spline_o%ibctype2(i) = 0
              ier = 56
           else
              spline_o%ibctype2(i) = 1
           endif
        case (2)
           if(hspline(2).le.1) then
              spline_o%ibctype2(i) = 0
              ier = 56
           else
              spline_o%ibctype2(i) = 2
           endif
        case default
           ier = 2
           spline_o%ibctype2(i) = 0
        end select
     endif
 
     spline_o%bcval3min(1:idim1, 1:idim2) = 0.0
     spline_o%bcval3max(1:idim1, 1:idim2) = 0.0
     spline_o%ibctype3(i) = 0
     if(present(BCS3)) then
        select case(BCS3(i))
        case (-1)
           spline_o%ibctype3(i) = -1
        case (0)
           spline_o%ibctype3(i) = 0
        case (1)
           if(hspline(3).le.0) then
              spline_o%ibctype3(i) = 0
              ier = 2
           else
              spline_o%ibctype3(i) = 1
           endif
        case (2)
           if(hspline(3).le.1) then
              spline_o%ibctype3(i) = 0
              ier = 2
           else
              spline_o%ibctype3(i) = 2
           endif
        case default
           ier = 2
           spline_o%ibctype3(i) = 0
        end select
     endif
 
  enddo

  if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
     spline_o%ibctype1(1)=-1
     spline_o%ibctype1(2)=-1
  endif
  if(spline_o%ibctype2(1)==-1 .OR. spline_o%ibctype2(2)==-1) then
     spline_o%ibctype2(1)=-1
     spline_o%ibctype2(2)=-1
  endif
  if(spline_o%ibctype3(1)==-1 .OR. spline_o%ibctype3(2)==-1) then
     spline_o%ibctype3(1)=-1
     spline_o%ibctype3(2)=-1
  endif
 
  ! default grid
  spline_o%x1min = 0.0_ezspline_r4
  spline_o%x1max = 1.0_ezspline_r4
  if(spline_o%ibctype1(1).eq.-1) spline_o%x1max = ezspline_twopi_r4
  spline_o%x1 = spline_o%x1min + (spline_o%x1max - spline_o%x1min)* &
       & (/ (real(i-1,ezspline_r4)/real(spline_o%n1-1, ezspline_r4), i=1,spline_o%n1 ) /)

  spline_o%x2min = 0.0_ezspline_r4
  spline_o%x2max = 1.0_ezspline_r4
  if(spline_o%ibctype2(1).eq.-1) spline_o%x2max = ezspline_twopi_r4
  spline_o%x2 = spline_o%x2min + (spline_o%x2max - spline_o%x2min)* &
       & (/ (real(i-1,ezspline_r4)/real(spline_o%n2-1, ezspline_r4), i=1,spline_o%n2 ) /)
 
  spline_o%x3min = 0.0_ezspline_r4
  spline_o%x3max = 1.0_ezspline_r4
  if(spline_o%ibctype3(1).eq.-1) spline_o%x3max = ezspline_twopi_r4
  spline_o%x3 = spline_o%x3min + (spline_o%x3max - spline_o%x3min)* &
       & (/ (real(i-1,ezspline_r4)/real(spline_o%n3-1, ezspline_r4), i=1,spline_o%n3 ) /)
 
  spline_o%isReady = 0
 
  return
end subroutine EZhybrid_init3_r4
 
