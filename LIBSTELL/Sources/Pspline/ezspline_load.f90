!/////
! R8 !
!/////
!
! 1-D
!
 
!  DMC modifications -- spl_name optional argument.
!  a single file can contain more than one spline (or not).  If it does
!  contain more than one spline, then each one must have a distinct name,
!  which must be specified in the call to tell which one to read this time.

!  fullsv -- ezspline_save now has a "full_save" option which means, save
!  the spline including coefficients.  ezspline_load is modified to look
!  for and read fullsv data if it is present; if not, the spline coefficients
!  are recalculated.
 
   subroutine EZspline_load1_r8(spline_o, filename, ier, spl_name)
     use ezspline_obj
     use ezcdf
     implicit none
     type(EZspline1_r8) :: spline_o
     character*(*) :: filename
     ! ier:
     ! 21=could not load spline object from file filename
     ! 22=loaded spline object from file filename but failed at coefficient set-up
     integer, intent(out) :: ier

     character*(*), intent(in), OPTIONAL :: spl_name  ! specify name of spline
     !  to be read-- in case file contains more than one.

     integer ncid, ifail
     integer dimlens(3)
     character(2), parameter :: real8='R8'
     character(3), parameter :: int='INT'
     integer n1, BCS1(2)
     real(ezspline_r8), dimension(:), allocatable :: f
 
     character*4 :: data_type

     character*32 zpre
     logical :: fullsv

     if(present(spl_name)) then
        zpre = spl_name
     else
        zpre = ' '
     endif

     ier = 0
     call cdfOpn(ncid, filename, 'r') ! no error flag??
     if(ncid==0) then
        ier=43
        return
     endif
 
     ! check if n1 is present; check if spline object has correct rank
     !                         n2 should NOT be present.

     call cdfInqVar(ncid, trim(zpre)//'n1', dimlens, data_type, ifail)
     if(ifail.ne.0) then
        ier = 21
        return
     endif
     call cdfInqVar(ncid, trim(zpre)//'n2', dimlens, data_type, ifail)
     if(ifail.eq.0) then
        ier = 20
        return
     endif

     ! check if fullsv was in effect -- if so, x1pkg will be found.

     fullsv=.FALSE.
     call cdfInqVar(ncid, trim(zpre)//'x1pkg', dimlens, data_type, ifail)
     if(ifail.eq.0) then
        fullsv=.TRUE.
     endif

     if( spline_o%nguard /= 123456789 ) call EZspline_preInit(spline_o)

     if( ezspline_allocated(spline_o) ) call EZspline_free1_r8(spline_o, ier)

     call cdfGetVar(ncid, trim(zpre)//'n1', n1, ifail)
     call cdfGetVar(ncid, trim(zpre)//'isLinear', spline_o%isLinear, ifail)
     if(ifail.ne.0) spline_o%isLinear=0

     if(spline_o%isLinear.eq.1) then
        call EZlinear_init1_r8(spline_o, n1, ier)
     else
        BCS1 = (/0, 0/)
        call EZspline_init1_r8(spline_o, n1, BCS1, ier)
     endif
     if(ier.ne.0) return

     call cdfGetVar(ncid, trim(zpre)//'klookup1', spline_o%klookup1, ifail)
     if(ifail.ne.0) spline_o%klookup1=-3  ! the old default
 
     call cdfGetVar(ncid, trim(zpre)//'isHermite', spline_o%isHermite, ifail)
     call cdfGetVar(ncid, trim(zpre)//'isReady', spline_o%isReady, ifail)
     call cdfGetVar(ncid, trim(zpre)//'ibctype1', spline_o%ibctype1, ifail)
     call cdfGetVar(ncid, trim(zpre)//'x1', spline_o%x1, ifail)
     call cdfGetVar(ncid, trim(zpre)//'bcval1min', spline_o%bcval1min, ifail)
     call cdfGetVar(ncid, trim(zpre)//'bcval1max', spline_o%bcval1max, ifail)
     if(ifail/=0) then
        ier=44
        return
     endif
 
     if(.not.fullsv) then
        allocate(f(n1))
        if(spline_o%isReady==1) then
           call cdfGetVar(ncid, trim(zpre)//'f', f, ifail)
        else
           ier = 92
        endif
     else
        call cdfGetVar(ncid, trim(zpre)//'x1min', spline_o%x1min, ifail)
        call cdfGetVar(ncid, trim(zpre)//'x1max', spline_o%x1max, ifail)
        call cdfGetVar(ncid, trim(zpre)//'ilin1', spline_o%ilin1, ifail)
        call cdfGetVar(ncid, trim(zpre)//'x1pkg', spline_o%x1pkg, ifail)
        call cdfGetVar(ncid, trim(zpre)//'fspl', spline_o%fspl, ifail)
     endif

     if(ifail /=0) then
        ier = 21
        return
     endif
 
     call cdfCls(ncid)
  
     if(.not.fullsv) then
        call EZspline_setup1_r8x(spline_o, f, ier)
        deallocate(f)
        if(ifail /=0) ier = 22
     endif
 
   end subroutine EZspline_load1_r8
 
 
 
 
 
!!
!! 2-D
!!
 
 
   subroutine EZspline_load2_r8(spline_o, filename, ier, spl_name)
     use ezspline_obj
     use ezcdf
     implicit none
     type(EZspline2_r8) :: spline_o
     character*(*) :: filename
     ! ier:
     ! 20=attempt to load spline object with incorrect rank
     ! 21=could not load spline object from file filename
     ! 22=loaded spline object from file filename but failed at coefficient set-up
     integer, intent(out) :: ier

     character*(*), intent(in), OPTIONAL :: spl_name  ! specify name of spline
     !  to be read-- in case file contains more than one.

     integer ncid, ifail, in0, in1, in2
     integer dimlens(3)
     character(2), parameter :: real8='R8'
     character(3), parameter :: int='INT'
     integer n1, n2, BCS1(2), BCS2(2), hspline(2)
     real(ezspline_r8), dimension(:,:), allocatable :: f
 
     character*4 :: data_type
 
     character*32 zpre
     logical :: fullsv

     if(present(spl_name)) then
        zpre = spl_name
     else
        zpre = ' '
     endif

     ier = 0
     call cdfOpn(ncid, filename, 'r') ! no error flag??
     if(ncid==0) then
        ier=43
        return
     endif
 
     ! check if n1 is present; check if spline object has correct rank
     !                         n2 should be present but not n3...

     call cdfInqVar(ncid, trim(zpre)//'n1', dimlens, data_type, ifail)
     if(ifail.ne.0) then
        ier = 21
        return
     endif
     call cdfInqVar(ncid, trim(zpre)//'n2', dimlens, data_type, ifail)
     if(ifail.ne.0) then
        ier = 20
        return
     endif
     call cdfInqVar(ncid, trim(zpre)//'n3', dimlens, data_type, ifail)
     if(ifail.eq.0) then
        ier = 20
        return
     endif
 
     ! check if fullsv was in effect -- if so, x1pkg will be found.

     fullsv=.FALSE.
     call cdfInqVar(ncid, trim(zpre)//'x1pkg', dimlens, data_type, ifail)
     if(ifail.eq.0) then
        fullsv=.TRUE.
     endif

     if(spline_o%nguard /= 123456789 ) call EZspline_preInit(spline_o)

     if( ezspline_allocated(spline_o) ) call EZspline_free2_r8(spline_o, ier)

     call cdfGetVar(ncid, trim(zpre)//'n1', n1, ifail)
     call cdfGetVar(ncid, trim(zpre)//'n2', n2, ifail)
     call cdfGetVar(ncid, trim(zpre)//'isLinear', spline_o%isLinear, ifail)
     if(ifail.ne.0) spline_o%isLinear=0
     call cdfGetVar(ncid, trim(zpre)//'isHybrid', spline_o%isHybrid, ifail)
     if(ifail.ne.0) spline_o%isHybrid=0
     spline_o%hspline = 0

     if(spline_o%isLinear.eq.1) then
        call EZlinear_init2_r8(spline_o, n1, n2, ier)
     else if(spline_o%isHybrid.eq.1) then
        call cdfGetVar(ncid, trim(zpre)//'hspline', hspline, ifail)
        call EZhybrid_init2_r8x(spline_o, n1, n2, hspline, ier)
     else
        BCS1 = (/0, 0/); BCS2 = (/0, 0/)
        call EZspline_init2_r8(spline_o, n1, n2, BCS1, BCS2, ier)
     endif
     if(ier.ne.0) return

     call cdfGetVar(ncid, trim(zpre)//'klookup1', spline_o%klookup1, ifail)
     if(ifail.ne.0) spline_o%klookup1=-3  ! the old default
     call cdfGetVar(ncid, trim(zpre)//'klookup2', spline_o%klookup2, ifail)
     if(ifail.ne.0) spline_o%klookup2=-3  ! the old default
 
     call cdfGetVar(ncid, trim(zpre)//'isHermite', spline_o%isHermite, ifail)
     call cdfGetVar(ncid, trim(zpre)//'isReady', spline_o%isReady, ifail)
     call cdfGetVar(ncid, trim(zpre)//'ibctype1', spline_o%ibctype1, ifail)
     call cdfGetVar(ncid, trim(zpre)//'ibctype2', spline_o%ibctype2, ifail)
     call cdfGetVar(ncid, trim(zpre)//'x1', spline_o%x1, ifail)
     call cdfGetVar(ncid, trim(zpre)//'x2', spline_o%x2, ifail)
     call cdfGetVar(ncid, trim(zpre)//'bcval1min', spline_o%bcval1min, ifail)
     call cdfGetVar(ncid, trim(zpre)//'bcval1max', spline_o%bcval1max, ifail)
     call cdfGetVar(ncid, trim(zpre)//'bcval2min', spline_o%bcval2min, ifail)
     call cdfGetVar(ncid, trim(zpre)//'bcval2max', spline_o%bcval2max, ifail)
     if(ifail/=0) then
        ier=44
        return
     endif

     in0=size(spline_o%fspl,1)
     in1=size(spline_o%fspl,2)
     in2=size(spline_o%fspl,3)
 
     if(.not.fullsv) then
        allocate(f(in1,in2))
        if(spline_o%isReady==1) then
           call cdfGetVar(ncid, trim(zpre)//'f', f, ifail)
        else
           ier = 92
        endif
     else
        call cdfGetVar(ncid, trim(zpre)//'x1min', spline_o%x1min, ifail)
        call cdfGetVar(ncid, trim(zpre)//'x1max', spline_o%x1max, ifail)
        call cdfGetVar(ncid, trim(zpre)//'ilin1', spline_o%ilin1, ifail)
        call cdfGetVar(ncid, trim(zpre)//'x2min', spline_o%x2min, ifail)
        call cdfGetVar(ncid, trim(zpre)//'x2max', spline_o%x2max, ifail)
        call cdfGetVar(ncid, trim(zpre)//'ilin2', spline_o%ilin2, ifail)
        call cdfGetVar(ncid, trim(zpre)//'x1pkg', spline_o%x1pkg, ifail)
        call cdfGetVar(ncid, trim(zpre)//'x2pkg', spline_o%x2pkg, ifail)
        call cdfGetVar(ncid, trim(zpre)//'fspl', spline_o%fspl, ifail)
     endif

     if(ifail /=0) then
        ier = 21
        return
     endif
 
     call cdfCls(ncid)
 
     if(.not.fullsv) then
        call EZspline_setup2_r8x(spline_o, f, ier)
        deallocate(f)
        if(ifail /=0) ier = 22
     endif
 
   end subroutine EZspline_load2_r8
 
 
!!!
!!! 3-D
!!!
 
 
   subroutine EZspline_load3_r8(spline_o, filename, ier, spl_name)
     use ezspline_obj
     use ezcdf
     implicit none
     type(EZspline3_r8) :: spline_o
     character*(*) :: filename
     ! ier:
     ! 21=could not load spline object from file filename
     ! 22=loaded spline object from file filename but failed at coefficient set-up
     integer, intent(out) :: ier

     character*(*), intent(in), OPTIONAL :: spl_name  ! specify name of spline
     !  to be read-- in case file contains more than one.

     integer ncid, ifail, in0, in1, in2, in3
     integer dimlens(3)
     character(2), parameter :: real8='R8'
     character(3), parameter :: int='INT'
     integer n1, n2, n3, BCS1(2), BCS2(2), BCS3(2), hspline(3)
     real(ezspline_r8), dimension(:,:,:), allocatable :: f
  
     character*4 :: data_type

     character*32 zpre
     logical :: fullsv

     if(present(spl_name)) then
        zpre = spl_name
     else
        zpre = ' '
     endif

     ier = 0
     call cdfOpn(ncid, filename, 'r') ! no error flag??
     if(ncid==0) then
        ier=43
        return
     endif
 
     ! check if n1 is present; check if spline object has correct rank
     !                         n2 and n3 should also be present.

     call cdfInqVar(ncid, trim(zpre)//'n1', dimlens, data_type, ifail)
     if(ifail.ne.0) then
        ier = 21
        return
     endif
     call cdfInqVar(ncid, trim(zpre)//'n2', dimlens, data_type, ifail)
     if(ifail.ne.0) then
        ier = 20
        return
     endif
     call cdfInqVar(ncid, trim(zpre)//'n3', dimlens, data_type, ifail)
     if(ifail.ne.0) then
        ier = 20
        return
     endif

     ! check if fullsv was in effect -- if so, x1pkg will be found.

     fullsv=.FALSE.
     call cdfInqVar(ncid, trim(zpre)//'x1pkg', dimlens, data_type, ifail)
     if(ifail.eq.0) then
        fullsv=.TRUE.
     endif
 
     if(spline_o%nguard /= 123456789 ) call EZspline_preInit(spline_o)

     if( ezspline_allocated(spline_o) ) call EZspline_free3_r8(spline_o, ier)

     call cdfGetVar(ncid, trim(zpre)//'n1', n1, ifail)
     call cdfGetVar(ncid, trim(zpre)//'n2', n2, ifail)
     call cdfGetVar(ncid, trim(zpre)//'n3', n3, ifail)
     call cdfGetVar(ncid, trim(zpre)//'isLinear', spline_o%isLinear, ifail)
     if(ifail.ne.0) spline_o%isLinear=0
     call cdfGetVar(ncid, trim(zpre)//'isHybrid', spline_o%isHybrid, ifail)
     if(ifail.ne.0) spline_o%isHybrid=0
     spline_o%hspline = 0

     if(spline_o%isLinear.eq.1) then
        call EZlinear_init3_r8(spline_o, n1, n2, n3, ier)
     else if(spline_o%isHybrid.eq.1) then
        call cdfGetVar(ncid, trim(zpre)//'hspline', hspline, ifail)
        call EZhybrid_init3_r8x(spline_o, n1, n2, n3, hspline, ier)
     else
        BCS1 = (/0, 0/); BCS2 = (/0, 0/); BCS3 = (/0, 0/)
        call EZspline_init3_r8(spline_o, n1, n2, n3, BCS1, BCS2, BCS3, ier)
     endif
     if(ier.ne.0) return

     call cdfGetVar(ncid, trim(zpre)//'klookup1', spline_o%klookup1, ifail)
     if(ifail.ne.0) spline_o%klookup1=-3  ! the old default
     call cdfGetVar(ncid, trim(zpre)//'klookup2', spline_o%klookup2, ifail)
     if(ifail.ne.0) spline_o%klookup2=-3  ! the old default
     call cdfGetVar(ncid, trim(zpre)//'klookup3', spline_o%klookup3, ifail)
     if(ifail.ne.0) spline_o%klookup3=-3  ! the old default
 
     call cdfGetVar(ncid, trim(zpre)//'isHermite', spline_o%isHermite, ifail)
     call cdfGetVar(ncid, trim(zpre)//'isReady', spline_o%isReady, ifail)
     call cdfGetVar(ncid, trim(zpre)//'ibctype1', spline_o%ibctype1, ifail)
     call cdfGetVar(ncid, trim(zpre)//'ibctype2', spline_o%ibctype2, ifail)
     call cdfGetVar(ncid, trim(zpre)//'ibctype3', spline_o%ibctype3, ifail)
     call cdfGetVar(ncid, trim(zpre)//'x1', spline_o%x1, ifail)
     call cdfGetVar(ncid, trim(zpre)//'x2', spline_o%x2, ifail)
     call cdfGetVar(ncid, trim(zpre)//'x3', spline_o%x3, ifail)
     call cdfGetVar(ncid, trim(zpre)//'bcval1min', spline_o%bcval1min, ifail)
     call cdfGetVar(ncid, trim(zpre)//'bcval1max', spline_o%bcval1max, ifail)
     call cdfGetVar(ncid, trim(zpre)//'bcval2min', spline_o%bcval2min, ifail)
     call cdfGetVar(ncid, trim(zpre)//'bcval2max', spline_o%bcval2max, ifail)
     call cdfGetVar(ncid, trim(zpre)//'bcval3min', spline_o%bcval3min, ifail)
     call cdfGetVar(ncid, trim(zpre)//'bcval3max', spline_o%bcval3max, ifail)
     if(ifail/=0) then
        ier=44
        return
     endif
 
     in0=size(spline_o%fspl,1)
     in1=size(spline_o%fspl,2)
     in2=size(spline_o%fspl,3)
     in3=size(spline_o%fspl,4)

     if(.not.fullsv) then
        allocate(f(in1,in2,in3))
        if(spline_o%isReady==1) then
           call cdfGetVar(ncid, trim(zpre)//'f', f, ifail)
        else
           ier = 92
        endif
     else
        call cdfGetVar(ncid, trim(zpre)//'x1min', spline_o%x1min, ifail)
        call cdfGetVar(ncid, trim(zpre)//'x1max', spline_o%x1max, ifail)
        call cdfGetVar(ncid, trim(zpre)//'ilin1', spline_o%ilin1, ifail)
        call cdfGetVar(ncid, trim(zpre)//'x2min', spline_o%x2min, ifail)
        call cdfGetVar(ncid, trim(zpre)//'x2max', spline_o%x2max, ifail)
        call cdfGetVar(ncid, trim(zpre)//'ilin2', spline_o%ilin2, ifail)
        call cdfGetVar(ncid, trim(zpre)//'x3min', spline_o%x3min, ifail)
        call cdfGetVar(ncid, trim(zpre)//'x3max', spline_o%x3max, ifail)
        call cdfGetVar(ncid, trim(zpre)//'ilin3', spline_o%ilin3, ifail)
        call cdfGetVar(ncid, trim(zpre)//'x1pkg', spline_o%x1pkg, ifail)
        call cdfGetVar(ncid, trim(zpre)//'x2pkg', spline_o%x2pkg, ifail)
        call cdfGetVar(ncid, trim(zpre)//'x3pkg', spline_o%x3pkg, ifail)
        call ezspline_cdfget3(ncid, trim(zpre)//'fspl', spline_o%fspl, &
             in0*in1, in2, in3, ifail)
     endif

     if(ifail /=0) then
        ier = 21
        return
     endif
 
     call cdfCls(ncid)
 
     if(.not.fullsv) then
        call EZspline_setup3_r8x(spline_o, f, ier)
        deallocate(f)
        if(ifail /=0) ier = 22
     endif
 
   end subroutine EZspline_load3_r8
!/////
! R4 !
!/////
!
! 1-D
!
 
 
   subroutine EZspline_load1_r4(spline_o, filename, ier, spl_name)
     use ezspline_obj
     use ezcdf
     implicit none
     type(EZspline1_r4) :: spline_o
     character*(*) :: filename
     ! ier:
     ! 21=could not load spline object from file filename
     ! 22=loaded spline object from file filename but failed at coefficient set-up
     integer, intent(out) :: ier

     character*(*), intent(in), OPTIONAL :: spl_name  ! specify name of spline
     !  to be read-- in case file contains more than one.

     integer ncid, ifail
     integer dimlens(3)
     character(2), parameter :: real8='R8'
     character(3), parameter :: int='INT'
     integer n1, BCS1(2)
     real(ezspline_r4), dimension(:), allocatable :: f
 
     character*4 :: data_type
 
     character*32 zpre
     logical :: fullsv

     if(present(spl_name)) then
        zpre = spl_name
     else
        zpre = ' '
     endif

     ier = 0
     call cdfOpn(ncid, filename, 'r') ! no error flag??
     if(ncid==0) then
        ier=43
        return
     endif
 
     ! check if n1 is present; check if spline object has correct rank
     !                         n2 should NOT be present.

     call cdfInqVar(ncid, trim(zpre)//'n1', dimlens, data_type, ifail)
     if(ifail.ne.0) then
        ier = 21
        return
     endif
     call cdfInqVar(ncid, trim(zpre)//'n2', dimlens, data_type, ifail)
     if(ifail.eq.0) then
        ier = 20
        return
     endif

     ! check if fullsv was in effect -- if so, x1pkg will be found.

     fullsv=.FALSE.
     call cdfInqVar(ncid, trim(zpre)//'x1pkg', dimlens, data_type, ifail)
     if(ifail.eq.0) then
        fullsv=.TRUE.
     endif
  
     if(spline_o%nguard /= 123456789 ) call EZspline_preInit(spline_o)

     if( ezspline_allocated(spline_o) ) call EZspline_free1_r4(spline_o, ier)

     call cdfGetVar(ncid, trim(zpre)//'n1', n1, ifail)
     call cdfGetVar(ncid, trim(zpre)//'isLinear', spline_o%isLinear, ifail)
     if(ifail.ne.0) spline_o%isLinear=0

     if(spline_o%isLinear.eq.1) then
        call EZlinear_init1_r4(spline_o, n1, ier)
     else
        BCS1 = (/0, 0/)
        call EZspline_init1_r4(spline_o, n1, BCS1, ier)
     endif
     if(ier.ne.0) return

     call cdfGetVar(ncid, trim(zpre)//'klookup1', spline_o%klookup1, ifail)
     if(ifail.ne.0) spline_o%klookup1=-3  ! the old default
 
     call cdfGetVar(ncid, trim(zpre)//'isHermite', spline_o%isHermite, ifail)
     call cdfGetVar(ncid, trim(zpre)//'isReady', spline_o%isReady, ifail)
     call cdfGetVar(ncid, trim(zpre)//'ibctype1', spline_o%ibctype1, ifail)
     call cdfGetVar(ncid, trim(zpre)//'x1', spline_o%x1, ifail)
     call cdfGetVar(ncid, trim(zpre)//'bcval1min', spline_o%bcval1min, ifail)
     call cdfGetVar(ncid, trim(zpre)//'bcval1max', spline_o%bcval1max, ifail)
     if(ifail/=0) then
        ier=44
        return
     endif
 
     if(.not.fullsv) then
        allocate(f(n1))
        if(spline_o%isReady==1) then
           call cdfGetVar(ncid, trim(zpre)//'f', f, ifail)
        else
           ier = 92
        endif
     else
        call cdfGetVar(ncid, trim(zpre)//'x1min', spline_o%x1min, ifail)
        call cdfGetVar(ncid, trim(zpre)//'x1max', spline_o%x1max, ifail)
        call cdfGetVar(ncid, trim(zpre)//'ilin1', spline_o%ilin1, ifail)
        call cdfGetVar(ncid, trim(zpre)//'x1pkg', spline_o%x1pkg, ifail)
        call cdfGetVar(ncid, trim(zpre)//'fspl', spline_o%fspl, ifail)
     endif

     if(ifail /=0) then
        ier = 21
        return
     endif
 
     call cdfCls(ncid)
 
     if(.not.fullsv) then
        call EZspline_setup1_r4x(spline_o, f, ier)
        deallocate(f)
        if(ifail /=0) ier = 22
     endif
 
   end subroutine EZspline_load1_r4
 
 
 
 
 
!!
!! 2-D
!!
 
 
   subroutine EZspline_load2_r4(spline_o, filename, ier, spl_name)
     use ezspline_obj
     use ezcdf
     implicit none
     type(EZspline2_r4) :: spline_o
     character*(*) :: filename
     ! ier:
     ! 20=attempt to load spline object with incorrect rank
     ! 21=could not load spline object from file filename
     ! 22=loaded spline object from file filename but failed at coefficient set-up
     integer, intent(out) :: ier

     character*(*), intent(in), OPTIONAL :: spl_name  ! specify name of spline
     !  to be read-- in case file contains more than one.

     integer ncid, ifail, in0, in1, in2
     integer dimlens(3)
     character(2), parameter :: real8='R8'
     character(3), parameter :: int='INT'
     integer n1, n2, BCS1(2), BCS2(2), hspline(2)
     real(ezspline_r4), dimension(:,:), allocatable :: f
 
     character*4 :: data_type
 
     character*32 zpre
     logical :: fullsv

     if(present(spl_name)) then
        zpre = spl_name
     else
        zpre = ' '
     endif

     ier = 0
     call cdfOpn(ncid, filename, 'r') ! no error flag??
     if(ncid==0) then
        ier=43
        return
     endif
 
     ! check if n1 is present; check if spline object has correct rank
     !                         n2 should be present but not n3...

     call cdfInqVar(ncid, trim(zpre)//'n1', dimlens, data_type, ifail)
     if(ifail.ne.0) then
        ier = 21
        return
     endif
     call cdfInqVar(ncid, trim(zpre)//'n2', dimlens, data_type, ifail)
     if(ifail.ne.0) then
        ier = 20
        return
     endif
     call cdfInqVar(ncid, trim(zpre)//'n3', dimlens, data_type, ifail)
     if(ifail.eq.0) then
        ier = 20
        return
     endif

     ! check if fullsv was in effect -- if so, x1pkg will be found.

     fullsv=.FALSE.
     call cdfInqVar(ncid, trim(zpre)//'x1pkg', dimlens, data_type, ifail)
     if(ifail.eq.0) then
        fullsv=.TRUE.
     endif
  
     if(spline_o%nguard /= 123456789 ) call EZspline_preInit(spline_o)

     if( ezspline_allocated(spline_o) ) call EZspline_free2_r4(spline_o, ier)

     call cdfGetVar(ncid, trim(zpre)//'n1', n1, ifail)
     call cdfGetVar(ncid, trim(zpre)//'n2', n2, ifail)
     call cdfGetVar(ncid, trim(zpre)//'isLinear', spline_o%isLinear, ifail)
     if(ifail.ne.0) spline_o%isLinear=0
     call cdfGetVar(ncid, trim(zpre)//'isHybrid', spline_o%isHybrid, ifail)
     if(ifail.ne.0) spline_o%isHybrid=0
     spline_o%hspline = 0

     if(spline_o%isLinear.eq.1) then
        call EZlinear_init2_r4(spline_o, n1, n2, ier)
     else if(spline_o%isHybrid.eq.1) then
        call cdfGetVar(ncid, trim(zpre)//'hspline', hspline, ifail)
        call EZhybrid_init2_r4x(spline_o, n1, n2, hspline, ier)
     else
        BCS1 = (/0, 0/); BCS2 = (/0, 0/)
        call EZspline_init2_r4(spline_o, n1, n2, BCS1, BCS2, ier)
     endif
     if(ier.ne.0) return

     call cdfGetVar(ncid, trim(zpre)//'klookup1', spline_o%klookup1, ifail)
     if(ifail.ne.0) spline_o%klookup1=-3  ! the old default
     call cdfGetVar(ncid, trim(zpre)//'klookup2', spline_o%klookup2, ifail)
     if(ifail.ne.0) spline_o%klookup2=-3  ! the old default
 
     call cdfGetVar(ncid, trim(zpre)//'isHermite', spline_o%isHermite, ifail)
     call cdfGetVar(ncid, trim(zpre)//'isReady', spline_o%isReady, ifail)
     call cdfGetVar(ncid, trim(zpre)//'ibctype1', spline_o%ibctype1, ifail)
     call cdfGetVar(ncid, trim(zpre)//'ibctype2', spline_o%ibctype2, ifail)
     call cdfGetVar(ncid, trim(zpre)//'x1', spline_o%x1, ifail)
     call cdfGetVar(ncid, trim(zpre)//'x2', spline_o%x2, ifail)
     call cdfGetVar(ncid, trim(zpre)//'bcval1min', spline_o%bcval1min, ifail)
     call cdfGetVar(ncid, trim(zpre)//'bcval1max', spline_o%bcval1max, ifail)
     call cdfGetVar(ncid, trim(zpre)//'bcval2min', spline_o%bcval2min, ifail)
     call cdfGetVar(ncid, trim(zpre)//'bcval2max', spline_o%bcval2max, ifail)
     if(ifail/=0) then
        ier=44
        return
     endif
 
     in0=size(spline_o%fspl,1)
     in1=size(spline_o%fspl,2)
     in2=size(spline_o%fspl,3)
 
     if(.not.fullsv) then
        allocate(f(in1,in2))
        if(spline_o%isReady==1) then
           call cdfGetVar(ncid, trim(zpre)//'f', f, ifail)
        else
           ier = 92
        endif
     else
        call cdfGetVar(ncid, trim(zpre)//'x1min', spline_o%x1min, ifail)
        call cdfGetVar(ncid, trim(zpre)//'x1max', spline_o%x1max, ifail)
        call cdfGetVar(ncid, trim(zpre)//'ilin1', spline_o%ilin1, ifail)
        call cdfGetVar(ncid, trim(zpre)//'x2min', spline_o%x2min, ifail)
        call cdfGetVar(ncid, trim(zpre)//'x2max', spline_o%x2max, ifail)
        call cdfGetVar(ncid, trim(zpre)//'ilin2', spline_o%ilin2, ifail)
        call cdfGetVar(ncid, trim(zpre)//'x1pkg', spline_o%x1pkg, ifail)
        call cdfGetVar(ncid, trim(zpre)//'x2pkg', spline_o%x2pkg, ifail)
        call cdfGetVar(ncid, trim(zpre)//'fspl', spline_o%fspl, ifail)
     endif

     if(ifail /=0) then
        ier = 21
        return
     endif
 
     call cdfCls(ncid)
 
     if(.not.fullsv) then
        call EZspline_setup2_r4x(spline_o, f, ier)
        deallocate(f)
        if(ifail /=0) ier = 22
     endif
 
   end subroutine EZspline_load2_r4
 
 
!!!
!!! 3-D
!!!
 
 
   subroutine EZspline_load3_r4(spline_o, filename, ier, spl_name)
     use ezspline_obj
     use ezcdf
     implicit none
     type(EZspline3_r4) :: spline_o
     character*(*) :: filename
     ! ier:
     ! 21=could not load spline object from file filename
     ! 22=loaded spline object from file filename but failed at coefficient set-up
     integer, intent(out) :: ier

     character*(*), intent(in), OPTIONAL :: spl_name  ! specify name of spline
     !  to be read-- in case file contains more than one.

     integer ncid, ifail, in0, in1, in2, in3
     integer dimlens(3)
     character(2), parameter :: real8='R8'
     character(3), parameter :: int='INT'
     integer n1, n2, n3, BCS1(2), BCS2(2), BCS3(2), hspline(3)
     real(ezspline_r4), dimension(:,:,:), allocatable :: f
 
     character*4 :: data_type

     character*32 zpre
     logical :: fullsv

     if(present(spl_name)) then
        zpre = spl_name
     else
        zpre = ' '
     endif

     ier = 0
     call cdfOpn(ncid, filename, 'r') ! no error flag??
     if(ncid==0) then
        ier=43
        return
     endif
 
     ! check if n1 is present; check if spline object has correct rank
     !                         n2 and n3 should also be present.

     call cdfInqVar(ncid, trim(zpre)//'n1', dimlens, data_type, ifail)
     if(ifail.ne.0) then
        ier = 21
        return
     endif
     call cdfInqVar(ncid, trim(zpre)//'n2', dimlens, data_type, ifail)
     if(ifail.ne.0) then
        ier = 20
        return
     endif
     call cdfInqVar(ncid, trim(zpre)//'n3', dimlens, data_type, ifail)
     if(ifail.ne.0) then
        ier = 20
        return
     endif

     ! check if fullsv was in effect -- if so, x1pkg will be found.

     fullsv=.FALSE.
     call cdfInqVar(ncid, trim(zpre)//'x1pkg', dimlens, data_type, ifail)
     if(ifail.eq.0) then
        fullsv=.TRUE.
     endif
  
     if(spline_o%nguard /= 123456789 ) call EZspline_preInit(spline_o)

     if( ezspline_allocated(spline_o) ) call EZspline_free3_r4(spline_o, ier)

     call cdfGetVar(ncid, trim(zpre)//'n1', n1, ifail)
     call cdfGetVar(ncid, trim(zpre)//'n2', n2, ifail)
     call cdfGetVar(ncid, trim(zpre)//'n3', n3, ifail)
     call cdfGetVar(ncid, trim(zpre)//'isLinear', spline_o%isLinear, ifail)
     if(ifail.ne.0) spline_o%isLinear=0
     call cdfGetVar(ncid, trim(zpre)//'isHybrid', spline_o%isHybrid, ifail)
     if(ifail.ne.0) spline_o%isHybrid=0
     spline_o%hspline = 0

     if(spline_o%isLinear.eq.1) then
        call EZlinear_init3_r4(spline_o, n1, n2, n3, ier)
     else if(spline_o%isHybrid.eq.1) then
        call cdfGetVar(ncid, trim(zpre)//'hspline', hspline, ifail)
        call EZhybrid_init3_r4x(spline_o, n1, n2, n3, hspline, ier)
     else
        BCS1 = (/0, 0/); BCS2 = (/0, 0/); BCS3 = (/0, 0/)
        call EZspline_init3_r4(spline_o, n1, n2, n3, BCS1, BCS2, BCS3, ier)
     endif
     if(ier.ne.0) return

     call cdfGetVar(ncid, trim(zpre)//'klookup1', spline_o%klookup1, ifail)
     if(ifail.ne.0) spline_o%klookup1=-3  ! the old default
     call cdfGetVar(ncid, trim(zpre)//'klookup2', spline_o%klookup2, ifail)
     if(ifail.ne.0) spline_o%klookup2=-3  ! the old default
     call cdfGetVar(ncid, trim(zpre)//'klookup3', spline_o%klookup3, ifail)
     if(ifail.ne.0) spline_o%klookup3=-3  ! the old default
 
     call cdfGetVar(ncid, trim(zpre)//'isHermite', spline_o%isHermite, ifail)
     call cdfGetVar(ncid, trim(zpre)//'isReady', spline_o%isReady, ifail)
     call cdfGetVar(ncid, trim(zpre)//'ibctype1', spline_o%ibctype1, ifail)
     call cdfGetVar(ncid, trim(zpre)//'ibctype2', spline_o%ibctype2, ifail)
     call cdfGetVar(ncid, trim(zpre)//'ibctype3', spline_o%ibctype3, ifail)
     call cdfGetVar(ncid, trim(zpre)//'x1', spline_o%x1, ifail)
     call cdfGetVar(ncid, trim(zpre)//'x2', spline_o%x2, ifail)
     call cdfGetVar(ncid, trim(zpre)//'x3', spline_o%x3, ifail)
     call cdfGetVar(ncid, trim(zpre)//'bcval1min', spline_o%bcval1min, ifail)
     call cdfGetVar(ncid, trim(zpre)//'bcval1max', spline_o%bcval1max, ifail)
     call cdfGetVar(ncid, trim(zpre)//'bcval2min', spline_o%bcval2min, ifail)
     call cdfGetVar(ncid, trim(zpre)//'bcval2max', spline_o%bcval2max, ifail)
     call cdfGetVar(ncid, trim(zpre)//'bcval3min', spline_o%bcval3min, ifail)
     call cdfGetVar(ncid, trim(zpre)//'bcval3max', spline_o%bcval3max, ifail)
     if(ifail/=0) then
        ier=44
        return
     endif
 
     in0=size(spline_o%fspl,1)
     in1=size(spline_o%fspl,2)
     in2=size(spline_o%fspl,3)
     in3=size(spline_o%fspl,4)
 
     if(.not.fullsv) then
        allocate(f(in1,in2,in3))
        if(spline_o%isReady==1) then
           call cdfGetVar(ncid, trim(zpre)//'f', f, ifail)
        else
           ier = 92
        endif
     else
        call cdfGetVar(ncid, trim(zpre)//'x1min', spline_o%x1min, ifail)
        call cdfGetVar(ncid, trim(zpre)//'x1max', spline_o%x1max, ifail)
        call cdfGetVar(ncid, trim(zpre)//'ilin1', spline_o%ilin1, ifail)
        call cdfGetVar(ncid, trim(zpre)//'x2min', spline_o%x2min, ifail)
        call cdfGetVar(ncid, trim(zpre)//'x2max', spline_o%x2max, ifail)
        call cdfGetVar(ncid, trim(zpre)//'ilin2', spline_o%ilin2, ifail)
        call cdfGetVar(ncid, trim(zpre)//'x3min', spline_o%x3min, ifail)
        call cdfGetVar(ncid, trim(zpre)//'x3max', spline_o%x3max, ifail)
        call cdfGetVar(ncid, trim(zpre)//'ilin3', spline_o%ilin3, ifail)
        call cdfGetVar(ncid, trim(zpre)//'x1pkg', spline_o%x1pkg, ifail)
        call cdfGetVar(ncid, trim(zpre)//'x2pkg', spline_o%x2pkg, ifail)
        call cdfGetVar(ncid, trim(zpre)//'x3pkg', spline_o%x3pkg, ifail)
        call ezspline_cdfget3(ncid, trim(zpre)//'fspl', spline_o%fspl, &
             in0*in1, in2, in3, ifail)
     endif

     if(ifail /=0) then
        ier = 21
        return
     endif
 
     call cdfCls(ncid)
 
     if(.not.fullsv) then
        call EZspline_setup3_r4x(spline_o, f, ier)
        deallocate(f)
        if(ifail /=0) ier = 22
     endif
 
   end subroutine EZspline_load3_r4
