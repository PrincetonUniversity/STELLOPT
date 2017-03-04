!/////
! R8 !
!/////
!
!  MOD DMC Apr 2007 -- support Hybrid spline objects.
!   these objects can have "zonal step function" dimensions 1 less than
!   the grid dimension; also, the number of spline coefficients varies
!   depending on how many dimensions are using cubic interpolation.  So,
!   local variables in0, in1, in2, in3 are defined to distinguish these
!   dimensions from the grid dimensions spline_o%n1, etc.
!
! 1-D
!

!DEC$ IF DEFINED (NETCDF)
   subroutine EZspline_save1_r8(spline_o, filename, ier, &
        spl_name, fullsave)
     use EZspline_obj
     use EZcdf
     implicit none
     type(EZspline1_r8) :: spline_o
     character*(*) :: filename
     ! ier:
     ! 17=could not save spline object in file filename
     integer, intent(out) :: ier

     !  if SPL_NAME is set, APPEND to file instead of creating new; prepend
     !  name to all NetCDF data items.  This allows one file to contain
     !  multiple items. (new dmc Mar 2006)
     character*(*), intent(in), OPTIONAL :: SPL_NAME  ! optional spline "name"

     !  if FULLSAVE is set .TRUE., save derived coefficients along with data
     !  this saves recalculating the coefficients when file is read (Mar 2006)
     logical, intent(in), OPTIONAL :: FULLSAVE

     integer ncid, ifail, in0, in1
     integer dimlens(3)
     character(2), parameter :: real8='R8'
     character(3), parameter :: int='INT'

     character*32 zpre
     logical :: fullsv,imodify
     
     ier = 0
     if(spline_o%isReady /= 1) then
        ier = 94
        return
     endif

     in0 = size(spline_o%fspl,1)
     in1 = size(spline_o%fspl,2)

     if(present(spl_name)) then
        call ezspline_spl_name_chk(spl_name,ier)
        if(ier.ne.0) return
        zpre=spl_name
     else
        zpre=' '
     endif

     if(present(fullsave)) then
        fullsv=fullsave
     else
        fullsv=.FALSE.
     endif

     if(zpre.eq.' ') then
        call cdfOpn(ncid, filename, 'w') ! no error flag??
        imodify=.FALSE.
     else
        call cdfOpn(ncid, filename, 'v') ! Append if exists, else create new
        imodify=.TRUE.
     endif
     if(ncid==0) then
        ier = 39
        return
     endif

     dimlens = (/1, 1, 1/) ! scalar
     call ezspline_defVar(ncid, trim(zpre)//'n1', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'klookup1', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'isHermite', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'isLinear', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'isReady', dimlens, int, imodify, ier)
     dimlens = (/2, 1, 1/) 
     call ezspline_defVar(ncid, trim(zpre)//'ibctype1', dimlens, int, imodify, ier)
     dimlens = (/spline_o%n1, 1, 1/)
     call ezspline_defVar(ncid, trim(zpre)//'x1', dimlens, real8, imodify, ier)
     dimlens = (/1, 1, 1/)
     call ezspline_defVar(ncid, trim(zpre)//'bcval1min', dimlens, real8, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'bcval1max', dimlens, real8, imodify, ier)
     !
     ! save only the function values at the nodes
     !
     dimlens = (/spline_o%n1, 1, 1/)
     if(.not.fullsv) then
        call ezspline_defVar(ncid, trim(zpre)//'f', dimlens, real8, imodify, ier)
     else
        dimlens = (/1, 1, 1/) ! scalar
        call ezspline_defVar(ncid, trim(zpre)//'x1min', dimlens, real8, imodify, ier)
        call ezspline_defVar(ncid, trim(zpre)//'x1max', dimlens, real8, imodify, ier)
        call ezspline_defVar(ncid, trim(zpre)//'ilin1', dimlens, int, imodify, ier)
        dimlens = (/spline_o%n1, 4, 1/)
        call ezspline_defVar(ncid, trim(zpre)//'x1pkg', dimlens, real8, imodify, ier)
        dimlens = (/2, spline_o%n1, 1/)
        call ezspline_defVar(ncid, trim(zpre)//'fspl', dimlens, real8, imodify, ier)
     endif
     if(ier.ne.0) return

     call cdfPutVar(ncid, trim(zpre)//'n1', spline_o%n1, ifail)
     call cdfPutVar(ncid, trim(zpre)//'klookup1', spline_o%klookup1, ifail)
     call cdfPutVar(ncid, trim(zpre)//'isHermite', spline_o%isHermite, ifail)
     call cdfPutVar(ncid, trim(zpre)//'isLinear', spline_o%isLinear, ifail)
     call cdfPutVar(ncid, trim(zpre)//'isReady', spline_o%isReady, ifail)
     call cdfPutVar(ncid, trim(zpre)//'ibctype1', spline_o%ibctype1, ifail)
     call cdfPutVar(ncid, trim(zpre)//'x1', spline_o%x1, ifail)
     call cdfPutVar(ncid, trim(zpre)//'bcval1min', spline_o%bcval1min, ifail)
     call cdfPutVar(ncid, trim(zpre)//'bcval1max', spline_o%bcval1max, ifail)
     if(ifail/=0) then
        ier=40
        return
     endif

     if(spline_o%isReady == 1) then
        if(.not.fullsv) then
           call cdfPutVar(ncid, trim(zpre)//'f', spline_o%fspl(1,:), ifail)
        else
           call cdfPutVar(ncid, trim(zpre)//'x1min', spline_o%x1min, ifail)
           call cdfPutVar(ncid, trim(zpre)//'x1max', spline_o%x1max, ifail)
           call cdfPutVar(ncid, trim(zpre)//'ilin1', spline_o%ilin1, ifail)
           call cdfPutVar(ncid, trim(zpre)//'x1pkg', spline_o%x1pkg, ifail)
           call cdfPutVar(ncid, trim(zpre)//'fspl', spline_o%fspl, ifail)
        endif
     else
        ier = 93
     endif

     call cdfCls(ncid)

     if(ifail /=0) ier = 17

     return
   end subroutine EZspline_save1_r8

!! 
!! 2-D
!!


   subroutine EZspline_save2_r8(spline_o, filename, ier, &
        spl_name, fullsave)
     use EZspline_obj
     use EZcdf
     implicit none
     type(EZspline2_r8) :: spline_o
     character*(*) :: filename
     ! ier:
     ! 17=could not save spline object in file filename
     integer, intent(out) :: ier

     !  if SPL_NAME is set, APPEND to file instead of creating new; prepend
     !  name to all NetCDF data items.  This allows one file to contain
     !  multiple items. (new dmc Mar 2006)
     character*(*), intent(in), OPTIONAL :: SPL_NAME  ! optional spline "name"

     !  if FULLSAVE is set .TRUE., save derived coefficients along with data
     !  this saves recalculating the coefficients when file is read (Mar 2006)
     logical, intent(in), OPTIONAL :: FULLSAVE

     integer ncid, ifail, in0, in1, in2
     integer dimlens(3)
     character(2), parameter :: real8='R8'
     character(3), parameter :: int='INT'

     character*32 zpre
     logical :: fullsv,imodify
     
     ier = 0
     if(spline_o%isReady /= 1) then
        ier = 94
        return
     endif

     in0 = size(spline_o%fspl,1)
     in1 = size(spline_o%fspl,2)
     in2 = size(spline_o%fspl,3)

     if(present(spl_name)) then
        call ezspline_spl_name_chk(spl_name,ier)
        if(ier.ne.0) return
        zpre=spl_name
     else
        zpre=' '
     endif

     if(present(fullsave)) then
        fullsv=fullsave
     else
        fullsv=.FALSE.
     endif

     if(zpre.eq.' ') then
        call cdfOpn(ncid, filename, 'w') ! no error flag??
        imodify=.FALSE.
     else
        call cdfOpn(ncid, filename, 'v') ! Append if exists, else create new
        imodify=.TRUE.
     endif
     if(ncid==0) then
        ier = 39
        return
     endif

     dimlens = (/1, 1, 1/) ! scalar
     call ezspline_defVar(ncid, trim(zpre)//'n1', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'n2', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'klookup1', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'klookup2', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'isHermite', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'isLinear', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'isHybrid', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'isReady', dimlens, int, imodify, ier)
     dimlens = (/2, 1, 1/) 
     call ezspline_defvar(ncid, trim(zpre)//'hspline', dimlens, int, imodify, ier)
     dimlens = (/2, 1, 1/) 
     call ezspline_defVar(ncid, trim(zpre)//'ibctype1', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'ibctype2', dimlens, int, imodify, ier)
     dimlens = (/spline_o%n1, 1, 1/)
     call ezspline_defVar(ncid, trim(zpre)//'x1', dimlens, real8, imodify, ier)
     dimlens = (/spline_o%n2, 1, 1/)
     call ezspline_defVar(ncid, trim(zpre)//'x2', dimlens, real8, imodify, ier)
     dimlens = (/in2, 1, 1/)
     call ezspline_defVar(ncid, trim(zpre)//'bcval1min', dimlens, real8, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'bcval1max', dimlens, real8, imodify, ier)
     dimlens = (/in1, 1, 1/)
     call ezspline_defVar(ncid, trim(zpre)//'bcval2min', dimlens, real8, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'bcval2max', dimlens, real8, imodify, ier)
     !
     ! save only the function values at the nodes
     !
     dimlens = (/in1, in2, 1/)
     if(.not.fullsv) then
        call ezspline_defVar(ncid, trim(zpre)//'f', dimlens, real8, imodify, ier)
     else
        dimlens = (/1, 1, 1/) ! scalar
        call ezspline_defVar(ncid, trim(zpre)//'x1min', dimlens, real8, imodify, ier)
        call ezspline_defVar(ncid, trim(zpre)//'x1max', dimlens, real8, imodify, ier)
        call ezspline_defVar(ncid, trim(zpre)//'ilin1', dimlens, int, imodify, ier)
        call ezspline_defVar(ncid, trim(zpre)//'x2min', dimlens, real8, imodify, ier)
        call ezspline_defVar(ncid, trim(zpre)//'x2max', dimlens, real8, imodify, ier)
        call ezspline_defVar(ncid, trim(zpre)//'ilin2', dimlens, int, imodify, ier)
        dimlens = (/spline_o%n1, 4, 1/)
        call ezspline_defVar(ncid, trim(zpre)//'x1pkg', dimlens, real8, imodify, ier)
        dimlens = (/spline_o%n2, 4, 1/)
        call ezspline_defVar(ncid, trim(zpre)//'x2pkg', dimlens, real8, imodify, ier)
        dimlens = (/in0, in1, in2/)
        call ezspline_defVar(ncid, trim(zpre)//'fspl', dimlens, real8, imodify, ier)
     endif
     if(ier.ne.0) return

     call cdfPutVar(ncid, trim(zpre)//'n1', spline_o%n1, ifail)
     call cdfPutVar(ncid, trim(zpre)//'n2', spline_o%n2, ifail)
     call cdfPutVar(ncid, trim(zpre)//'klookup1', spline_o%klookup1, ifail)
     call cdfPutVar(ncid, trim(zpre)//'klookup2', spline_o%klookup2, ifail)
     call cdfPutVar(ncid, trim(zpre)//'isHermite', spline_o%isHermite, ifail)
     call cdfPutVar(ncid, trim(zpre)//'isLinear', spline_o%isLinear, ifail)
     call cdfPutVar(ncid, trim(zpre)//'isHybrid', spline_o%isHybrid, ifail)
     call cdfPutVar(ncid, trim(zpre)//'isReady', spline_o%isReady, ifail)
     call cdfPutVar(ncid, trim(zpre)//'hspline', spline_o%hspline, ifail)
     call cdfPutVar(ncid, trim(zpre)//'ibctype1', spline_o%ibctype1, ifail)
     call cdfPutVar(ncid, trim(zpre)//'ibctype2', spline_o%ibctype2, ifail)
     call cdfPutVar(ncid, trim(zpre)//'x1', spline_o%x1, ifail)
     call cdfPutVar(ncid, trim(zpre)//'x2', spline_o%x2, ifail)
     call cdfPutVar(ncid, trim(zpre)//'bcval1min', spline_o%bcval1min, ifail)
     call cdfPutVar(ncid, trim(zpre)//'bcval1max', spline_o%bcval1max, ifail)
     call cdfPutVar(ncid, trim(zpre)//'bcval2min', spline_o%bcval2min, ifail)
     call cdfPutVar(ncid, trim(zpre)//'bcval2max', spline_o%bcval2max, ifail)
     if(ifail/=0) then
        ier=40
        return
     endif

     if(spline_o%isReady == 1) then
        if(.not.fullsv) then
           call cdfPutVar(ncid, trim(zpre)//'f', spline_o%fspl(1,:,:), ifail)
        else
           call cdfPutVar(ncid, trim(zpre)//'x1min', spline_o%x1min, ifail)
           call cdfPutVar(ncid, trim(zpre)//'x1max', spline_o%x1max, ifail)
           call cdfPutVar(ncid, trim(zpre)//'ilin1', spline_o%ilin1, ifail)
           call cdfPutVar(ncid, trim(zpre)//'x2min', spline_o%x2min, ifail)
           call cdfPutVar(ncid, trim(zpre)//'x2max', spline_o%x2max, ifail)
           call cdfPutVar(ncid, trim(zpre)//'ilin2', spline_o%ilin2, ifail)
           call cdfPutVar(ncid, trim(zpre)//'x1pkg', spline_o%x1pkg, ifail)
           call cdfPutVar(ncid, trim(zpre)//'x2pkg', spline_o%x2pkg, ifail)
           call cdfPutVar(ncid, trim(zpre)//'fspl', spline_o%fspl, ifail)
        endif
     else
        ier = 93
     endif

     call cdfCls(ncid)

     if(ifail /=0) ier = 17

     return
   end subroutine EZspline_save2_r8





!!! 
!!! 3-D
!!!



   subroutine EZspline_save3_r8(spline_o, filename, ier, &
        spl_name, fullsave)
     use EZspline_obj
     use EZcdf
     implicit none
     type(EZspline3_r8) :: spline_o
     character*(*) :: filename
     ! ier:
     ! 17=could not save spline object in file filename
     integer, intent(out) :: ier

     !  if SPL_NAME is set, APPEND to file instead of creating new; prepend
     !  name to all NetCDF data items.  This allows one file to contain
     !  multiple items. (new dmc Mar 2006)
     character*(*), intent(in), OPTIONAL :: SPL_NAME  ! optional spline "name"

     !  if FULLSAVE is set .TRUE., save derived coefficients along with data
     !  this saves recalculating the coefficients when file is read (Mar 2006)
     logical, intent(in), OPTIONAL :: FULLSAVE

     integer ncid, ifail, in0, in1, in2, in3
     integer dimlens(3)
     character(2), parameter :: real8='R8'
     character(3), parameter :: int='INT'

     character*32 zpre
     logical :: fullsv,imodify
     
     ier = 0
     if(spline_o%isReady /= 1) then
        ier = 94
        return
     endif

     in0 = size(spline_o%fspl,1)
     in1 = size(spline_o%fspl,2)
     in2 = size(spline_o%fspl,3)
     in3 = size(spline_o%fspl,4)

     if(present(spl_name)) then
        call ezspline_spl_name_chk(spl_name,ier)
        if(ier.ne.0) return
        zpre=spl_name
     else
        zpre=' '
     endif

     if(present(fullsave)) then
        fullsv=fullsave
     else
        fullsv=.FALSE.
     endif

     if(zpre.eq.' ') then
        call cdfOpn(ncid, filename, 'w') ! no error flag??
        imodify=.FALSE.
     else
        call cdfOpn(ncid, filename, 'v') ! Append if exists, else create new
        imodify=.TRUE.
     endif
     if(ncid==0) then
        ier = 39
        return
     endif

     dimlens = (/1, 1, 1/) ! scalar
     call ezspline_defVar(ncid, trim(zpre)//'n1', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'n2', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'n3', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'klookup1', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'klookup2', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'klookup3', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'isHermite', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'isLinear', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'isHybrid', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'isReady', dimlens, int, imodify, ier)
     dimlens = (/3, 1, 1/) 
     call ezspline_defvar(ncid, trim(zpre)//'hspline', dimlens, int, imodify, ier)
     dimlens = (/2, 1, 1/) 
     call ezspline_defVar(ncid, trim(zpre)//'ibctype1', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'ibctype2', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'ibctype3', dimlens, int, imodify, ier)
     dimlens = (/spline_o%n1, 1, 1/)
     call ezspline_defVar(ncid, trim(zpre)//'x1', dimlens, real8, imodify, ier)
     dimlens = (/spline_o%n2, 1, 1/)
     call ezspline_defVar(ncid, trim(zpre)//'x2', dimlens, real8, imodify, ier)
     dimlens = (/spline_o%n3, 1, 1/)
     call ezspline_defVar(ncid, trim(zpre)//'x3', dimlens, real8, imodify, ier)
     dimlens = (/in2, in3, 1/)
     call ezspline_defVar(ncid, trim(zpre)//'bcval1min', dimlens, real8, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'bcval1max', dimlens, real8, imodify, ier)
     dimlens = (/in1, in3, 1/)
     call ezspline_defVar(ncid, trim(zpre)//'bcval2min', dimlens, real8, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'bcval2max', dimlens, real8, imodify, ier)
     dimlens = (/in1, in2, 1/)
     call ezspline_defVar(ncid, trim(zpre)//'bcval3min', dimlens, real8, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'bcval3max', dimlens, real8, imodify, ier)
     !
     ! save only the function values at the nodes
     !
     dimlens = (/in1, in2, in3/)
     if(.not.fullsv) then
        call ezspline_defVar(ncid, trim(zpre)//'f', dimlens, real8, imodify, ier)
     else
        dimlens = (/1, 1, 1/) ! scalar
        call ezspline_defVar(ncid, trim(zpre)//'x1min', dimlens, real8, imodify, ier)
        call ezspline_defVar(ncid, trim(zpre)//'x1max', dimlens, real8, imodify, ier)
        call ezspline_defVar(ncid, trim(zpre)//'ilin1', dimlens, int, imodify, ier)
        call ezspline_defVar(ncid, trim(zpre)//'x2min', dimlens, real8, imodify, ier)
        call ezspline_defVar(ncid, trim(zpre)//'x2max', dimlens, real8, imodify, ier)
        call ezspline_defVar(ncid, trim(zpre)//'ilin2', dimlens, int, imodify, ier)
        call ezspline_defVar(ncid, trim(zpre)//'x3min', dimlens, real8, imodify, ier)
        call ezspline_defVar(ncid, trim(zpre)//'x3max', dimlens, real8, imodify, ier)
        call ezspline_defVar(ncid, trim(zpre)//'ilin3', dimlens, int, imodify, ier)
        dimlens = (/spline_o%n1, 4, 1/)
        call ezspline_defVar(ncid, trim(zpre)//'x1pkg', dimlens, real8, imodify, ier)
        dimlens = (/spline_o%n2, 4, 1/)
        call ezspline_defVar(ncid, trim(zpre)//'x2pkg', dimlens, real8, imodify, ier)
        dimlens = (/spline_o%n3, 4, 1/)
        call ezspline_defVar(ncid, trim(zpre)//'x3pkg', dimlens, real8, imodify, ier)
        dimlens = (/in0*in1, in2, in3/)
        call ezspline_defVar(ncid, trim(zpre)//'fspl', dimlens, real8, imodify, ier)
     endif
     if(ier.ne.0) return

     call cdfPutVar(ncid, trim(zpre)//'n1', spline_o%n1, ifail)
     call cdfPutVar(ncid, trim(zpre)//'n2', spline_o%n2, ifail)
     call cdfPutVar(ncid, trim(zpre)//'n3', spline_o%n3, ifail)
     call cdfPutVar(ncid, trim(zpre)//'klookup1', spline_o%klookup1, ifail)
     call cdfPutVar(ncid, trim(zpre)//'klookup2', spline_o%klookup2, ifail)
     call cdfPutVar(ncid, trim(zpre)//'klookup3', spline_o%klookup3, ifail)
     call cdfPutVar(ncid, trim(zpre)//'isHermite', spline_o%isHermite, ifail)
     call cdfPutVar(ncid, trim(zpre)//'isLinear', spline_o%isLinear, ifail)
     call cdfPutVar(ncid, trim(zpre)//'isHybrid', spline_o%isHybrid, ifail)
     call cdfPutVar(ncid, trim(zpre)//'isReady', spline_o%isReady, ifail)
     call cdfPutVar(ncid, trim(zpre)//'hspline', spline_o%hspline, ifail)
     call cdfPutVar(ncid, trim(zpre)//'ibctype1', spline_o%ibctype1, ifail)
     call cdfPutVar(ncid, trim(zpre)//'ibctype2', spline_o%ibctype2, ifail)
     call cdfPutVar(ncid, trim(zpre)//'ibctype3', spline_o%ibctype3, ifail)
     call cdfPutVar(ncid, trim(zpre)//'x1', spline_o%x1, ifail)
     call cdfPutVar(ncid, trim(zpre)//'x2', spline_o%x2, ifail)
     call cdfPutVar(ncid, trim(zpre)//'x3', spline_o%x3, ifail)
     call cdfPutVar(ncid, trim(zpre)//'bcval1min', spline_o%bcval1min, ifail)
     call cdfPutVar(ncid, trim(zpre)//'bcval1max', spline_o%bcval1max, ifail)
     call cdfPutVar(ncid, trim(zpre)//'bcval2min', spline_o%bcval2min, ifail)
     call cdfPutVar(ncid, trim(zpre)//'bcval2max', spline_o%bcval2max, ifail)
     call cdfPutVar(ncid, trim(zpre)//'bcval3min', spline_o%bcval3min, ifail)
     call cdfPutVar(ncid, trim(zpre)//'bcval3max', spline_o%bcval3max, ifail)
     if(ifail/=0) then
        ier=40
        return
     endif

     if(spline_o%isReady == 1) then
        if(.not.fullsv) then
           call cdfPutVar(ncid, trim(zpre)//'f', spline_o%fspl(1,:,:,:), ifail)
        else
           call cdfPutVar(ncid, trim(zpre)//'x1min', spline_o%x1min, ifail)
           call cdfPutVar(ncid, trim(zpre)//'x1max', spline_o%x1max, ifail)
           call cdfPutVar(ncid, trim(zpre)//'ilin1', spline_o%ilin1, ifail)
           call cdfPutVar(ncid, trim(zpre)//'x2min', spline_o%x2min, ifail)
           call cdfPutVar(ncid, trim(zpre)//'x2max', spline_o%x2max, ifail)
           call cdfPutVar(ncid, trim(zpre)//'ilin2', spline_o%ilin2, ifail)
           call cdfPutVar(ncid, trim(zpre)//'x3min', spline_o%x3min, ifail)
           call cdfPutVar(ncid, trim(zpre)//'x3max', spline_o%x3max, ifail)
           call cdfPutVar(ncid, trim(zpre)//'ilin3', spline_o%ilin3, ifail)
           call cdfPutVar(ncid, trim(zpre)//'x1pkg', spline_o%x1pkg, ifail)
           call cdfPutVar(ncid, trim(zpre)//'x2pkg', spline_o%x2pkg, ifail)
           call cdfPutVar(ncid, trim(zpre)//'x3pkg', spline_o%x3pkg, ifail)
           call ezspline_cdfput3(ncid, trim(zpre)//'fspl', spline_o%fspl, &
                in0*in1, in2, in3, ifail)
        endif
     else
        ier = 93
     endif

     call cdfCls(ncid)

     if(ifail /=0) ier = 17

     return
   end subroutine EZspline_save3_r8
!/////
! R4 !
!/////
!
! 1-D
!


   subroutine EZspline_save1_r4(spline_o, filename, ier, &
        spl_name, fullsave)
     use EZspline_obj
     use EZcdf
     implicit none
     type(EZspline1_r4) :: spline_o
     character*(*) :: filename
     ! ier:
     ! 17=could not save spline object in file filename
     integer, intent(out) :: ier

     !  if SPL_NAME is set, APPEND to file instead of creating new; prepend
     !  name to all NetCDF data items.  This allows one file to contain
     !  multiple items. (new dmc Mar 2006)
     character*(*), intent(in), OPTIONAL :: SPL_NAME  ! optional spline "name"

     !  if FULLSAVE is set .TRUE., save derived coefficients along with data
     !  this saves recalculating the coefficients when file is read (Mar 2006)
     logical, intent(in), OPTIONAL :: FULLSAVE

     integer ncid, ifail, in0, in1
     integer dimlens(3)
     character(2), parameter :: real4='R4'
     character(3), parameter :: int='INT'

     character*32 zpre
     logical :: fullsv,imodify
     
     ier = 0
     if(spline_o%isReady /= 1) then
        ier = 94
        return
     endif

     in0 = size(spline_o%fspl,1)
     in1 = size(spline_o%fspl,2)

     if(present(spl_name)) then
        call ezspline_spl_name_chk(spl_name,ier)
        if(ier.ne.0) return
        zpre=spl_name
     else
        zpre=' '
     endif

     if(present(fullsave)) then
        fullsv=fullsave
     else
        fullsv=.FALSE.
     endif

     if(zpre.eq.' ') then
        call cdfOpn(ncid, filename, 'w') ! no error flag??
        imodify=.FALSE.
     else
        call cdfOpn(ncid, filename, 'v') ! Append if exists, else create new
        imodify=.TRUE.
     endif
     if(ncid==0) then
        ier = 39
        return
     endif

     dimlens = (/1, 1, 1/) ! scalar
     call ezspline_defVar(ncid, trim(zpre)//'n1', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'klookup1', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'isHermite', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'isLinear', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'isReady', dimlens, int, imodify, ier)
     dimlens = (/2, 1, 1/) 
     call ezspline_defVar(ncid, trim(zpre)//'ibctype1', dimlens, int, imodify, ier)
     dimlens = (/spline_o%n1, 1, 1/)
     call ezspline_defVar(ncid, trim(zpre)//'x1', dimlens, real4, imodify, ier)
     dimlens = (/1, 1, 1/)
     call ezspline_defVar(ncid, trim(zpre)//'bcval1min', dimlens, real4, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'bcval1max', dimlens, real4, imodify, ier)
     !
     ! save only the function values at the nodes
     !
     dimlens = (/spline_o%n1, 1, 1/)
     if(.not.fullsv) then
        call ezspline_defVar(ncid, trim(zpre)//'f', dimlens, real4, imodify, ier)
     else
        dimlens = (/1, 1, 1/) ! scalar
        call ezspline_defVar(ncid, trim(zpre)//'x1min', dimlens, real4, imodify, ier)
        call ezspline_defVar(ncid, trim(zpre)//'x1max', dimlens, real4, imodify, ier)
        call ezspline_defVar(ncid, trim(zpre)//'ilin1', dimlens, int, imodify, ier)
        dimlens = (/spline_o%n1, 4, 1/)
        call ezspline_defVar(ncid, trim(zpre)//'x1pkg', dimlens, real4, imodify, ier)
        dimlens = (/2, spline_o%n1, 1/)
        call ezspline_defVar(ncid, trim(zpre)//'fspl', dimlens, real4, imodify, ier)
     endif
     if(ier.ne.0) return

     call cdfPutVar(ncid, trim(zpre)//'n1', spline_o%n1, ifail)
     call cdfPutVar(ncid, trim(zpre)//'klookup1', spline_o%klookup1, ifail)
     call cdfPutVar(ncid, trim(zpre)//'isHermite', spline_o%isHermite, ifail)
     call cdfPutVar(ncid, trim(zpre)//'isLinear', spline_o%isLinear, ifail)
     call cdfPutVar(ncid, trim(zpre)//'isReady', spline_o%isReady, ifail)
     call cdfPutVar(ncid, trim(zpre)//'ibctype1', spline_o%ibctype1, ifail)
     call cdfPutVar(ncid, trim(zpre)//'x1', spline_o%x1, ifail)
     call cdfPutVar(ncid, trim(zpre)//'bcval1min', spline_o%bcval1min, ifail)
     call cdfPutVar(ncid, trim(zpre)//'bcval1max', spline_o%bcval1max, ifail)
     if(ifail/=0) then
        ier=40
        return
     endif

     if(spline_o%isReady == 1) then
        if(.not.fullsv) then
           call cdfPutVar(ncid, trim(zpre)//'f', spline_o%fspl(1,:), ifail)
        else
           call cdfPutVar(ncid, trim(zpre)//'x1min', spline_o%x1min, ifail)
           call cdfPutVar(ncid, trim(zpre)//'x1max', spline_o%x1max, ifail)
           call cdfPutVar(ncid, trim(zpre)//'ilin1', spline_o%ilin1, ifail)
           call cdfPutVar(ncid, trim(zpre)//'x1pkg', spline_o%x1pkg, ifail)
           call cdfPutVar(ncid, trim(zpre)//'fspl', spline_o%fspl, ifail)
        endif
     else
        ier = 93
     endif

     call cdfCls(ncid)

     if(ifail /=0) ier = 17

     return
   end subroutine EZspline_save1_r4

!! 
!! 2-D
!!


   subroutine EZspline_save2_r4(spline_o, filename, ier, &
        spl_name, fullsave)
     use EZspline_obj
     use EZcdf
     implicit none
     type(EZspline2_r4) :: spline_o
     character*(*) :: filename
     ! ier:
     ! 17=could not save spline object in file filename
     integer, intent(out) :: ier

     !  if SPL_NAME is set, APPEND to file instead of creating new; prepend
     !  name to all NetCDF data items.  This allows one file to contain
     !  multiple items. (new dmc Mar 2006)
     character*(*), intent(in), OPTIONAL :: SPL_NAME  ! optional spline "name"

     !  if FULLSAVE is set .TRUE., save derived coefficients along with data
     !  this saves recalculating the coefficients when file is read (Mar 2006)
     logical, intent(in), OPTIONAL :: FULLSAVE

     integer ncid, ifail, in0, in1, in2
     integer dimlens(3)
     character(2), parameter :: real4='R4'
     character(3), parameter :: int='INT'

     character*32 zpre
     logical :: fullsv,imodify
     
     ier = 0
     if(spline_o%isReady /= 1) then
        ier = 94
        return
     endif

     in0 = size(spline_o%fspl,1)
     in1 = size(spline_o%fspl,2)
     in2 = size(spline_o%fspl,3)

     if(present(spl_name)) then
        call ezspline_spl_name_chk(spl_name,ier)
        if(ier.ne.0) return
        zpre=spl_name
     else
        zpre=' '
     endif

     if(present(fullsave)) then
        fullsv=fullsave
     else
        fullsv=.FALSE.
     endif

     if(zpre.eq.' ') then
        call cdfOpn(ncid, filename, 'w') ! no error flag??
        imodify=.FALSE.
     else
        call cdfOpn(ncid, filename, 'v') ! Append if exists, else create new
        imodify=.TRUE.
     endif
     if(ncid==0) then
        ier = 39
        return
     endif

     dimlens = (/1, 1, 1/) ! scalar
     call ezspline_defVar(ncid, trim(zpre)//'n1', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'n2', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'klookup1', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'klookup2', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'isHermite', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'isLinear', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'isHybrid', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'isReady', dimlens, int, imodify, ier)
     dimlens = (/2, 1, 1/) 
     call ezspline_defvar(ncid, trim(zpre)//'hspline', dimlens, int, imodify, ier)
     dimlens = (/2, 1, 1/) 
     call ezspline_defVar(ncid, trim(zpre)//'ibctype1', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'ibctype2', dimlens, int, imodify, ier)
     dimlens = (/spline_o%n1, 1, 1/)
     call ezspline_defVar(ncid, trim(zpre)//'x1', dimlens, real4, imodify, ier)
     dimlens = (/spline_o%n2, 1, 1/)
     call ezspline_defVar(ncid, trim(zpre)//'x2', dimlens, real4, imodify, ier)
     dimlens = (/in2, 1, 1/)
     call ezspline_defVar(ncid, trim(zpre)//'bcval1min', dimlens, real4, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'bcval1max', dimlens, real4, imodify, ier)
     dimlens = (/in1, 1, 1/)
     call ezspline_defVar(ncid, trim(zpre)//'bcval2min', dimlens, real4, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'bcval2max', dimlens, real4, imodify, ier)
     !
     ! save only the function values at the nodes
     !
     dimlens = (/in1, in2, 1/)
     if(.not.fullsv) then
        call ezspline_defVar(ncid, trim(zpre)//'f', dimlens, real4, imodify, ier)
     else
        dimlens = (/1, 1, 1/) ! scalar
        call ezspline_defVar(ncid, trim(zpre)//'x1min', dimlens, real4, imodify, ier)
        call ezspline_defVar(ncid, trim(zpre)//'x1max', dimlens, real4, imodify, ier)
        call ezspline_defVar(ncid, trim(zpre)//'ilin1', dimlens, int, imodify, ier)
        call ezspline_defVar(ncid, trim(zpre)//'x2min', dimlens, real4, imodify, ier)
        call ezspline_defVar(ncid, trim(zpre)//'x2max', dimlens, real4, imodify, ier)
        call ezspline_defVar(ncid, trim(zpre)//'ilin2', dimlens, int, imodify, ier)
        dimlens = (/spline_o%n1, 4, 1/)
        call ezspline_defVar(ncid, trim(zpre)//'x1pkg', dimlens, real4, imodify, ier)
        dimlens = (/spline_o%n2, 4, 1/)
        call ezspline_defVar(ncid, trim(zpre)//'x2pkg', dimlens, real4, imodify, ier)
        dimlens = (/in0, in1, in2/)
        call ezspline_defVar(ncid, trim(zpre)//'fspl', dimlens, real4, imodify, ier)
     endif
     if(ier.ne.0) return

     call cdfPutVar(ncid, trim(zpre)//'n1', spline_o%n1, ifail)
     call cdfPutVar(ncid, trim(zpre)//'n2', spline_o%n2, ifail)
     call cdfPutVar(ncid, trim(zpre)//'klookup1', spline_o%klookup1, ifail)
     call cdfPutVar(ncid, trim(zpre)//'klookup2', spline_o%klookup2, ifail)
     call cdfPutVar(ncid, trim(zpre)//'isHermite', spline_o%isHermite, ifail)
     call cdfPutVar(ncid, trim(zpre)//'isLinear', spline_o%isLinear, ifail)
     call cdfPutVar(ncid, trim(zpre)//'isHybrid', spline_o%isHybrid, ifail)
     call cdfPutVar(ncid, trim(zpre)//'isReady', spline_o%isReady, ifail)
     call cdfPutVar(ncid, trim(zpre)//'hspline', spline_o%hspline, ifail)
     call cdfPutVar(ncid, trim(zpre)//'ibctype1', spline_o%ibctype1, ifail)
     call cdfPutVar(ncid, trim(zpre)//'ibctype2', spline_o%ibctype2, ifail)
     call cdfPutVar(ncid, trim(zpre)//'x1', spline_o%x1, ifail)
     call cdfPutVar(ncid, trim(zpre)//'x2', spline_o%x2, ifail)
     call cdfPutVar(ncid, trim(zpre)//'bcval1min', spline_o%bcval1min, ifail)
     call cdfPutVar(ncid, trim(zpre)//'bcval1max', spline_o%bcval1max, ifail)
     call cdfPutVar(ncid, trim(zpre)//'bcval2min', spline_o%bcval2min, ifail)
     call cdfPutVar(ncid, trim(zpre)//'bcval2max', spline_o%bcval2max, ifail)
     if(ifail/=0) then
        ier=40
        return
     endif

     if(spline_o%isReady == 1) then
        if(.not.fullsv) then
           call cdfPutVar(ncid, trim(zpre)//'f', spline_o%fspl(1,:,:), ifail)
        else
           call cdfPutVar(ncid, trim(zpre)//'x1min', spline_o%x1min, ifail)
           call cdfPutVar(ncid, trim(zpre)//'x1max', spline_o%x1max, ifail)
           call cdfPutVar(ncid, trim(zpre)//'ilin1', spline_o%ilin1, ifail)
           call cdfPutVar(ncid, trim(zpre)//'x2min', spline_o%x2min, ifail)
           call cdfPutVar(ncid, trim(zpre)//'x2max', spline_o%x2max, ifail)
           call cdfPutVar(ncid, trim(zpre)//'ilin2', spline_o%ilin2, ifail)
           call cdfPutVar(ncid, trim(zpre)//'x1pkg', spline_o%x1pkg, ifail)
           call cdfPutVar(ncid, trim(zpre)//'x2pkg', spline_o%x2pkg, ifail)
           call cdfPutVar(ncid, trim(zpre)//'fspl', spline_o%fspl, ifail)
        endif
     else
        ier = 93
     endif

     call cdfCls(ncid)

     if(ifail /=0) ier = 17

     return
   end subroutine EZspline_save2_r4





!!! 
!!! 3-D
!!!



   subroutine EZspline_save3_r4(spline_o, filename, ier, &
        spl_name, fullsave)
     use EZspline_obj
     use EZcdf
     implicit none
     type(EZspline3_r4) :: spline_o
     character*(*) :: filename
     ! ier:
     ! 17=could not save spline object in file filename
     integer, intent(out) :: ier

     !  if SPL_NAME is set, APPEND to file instead of creating new; prepend
     !  name to all NetCDF data items.  This allows one file to contain
     !  multiple items. (new dmc Mar 2006)
     character*(*), intent(in), OPTIONAL :: SPL_NAME  ! optional spline "name"

     !  if FULLSAVE is set .TRUE., save derived coefficients along with data
     !  this saves recalculating the coefficients when file is read (Mar 2006)
     logical, intent(in), OPTIONAL :: FULLSAVE

     integer ncid, ifail, in0, in1, in2, in3
     integer dimlens(3)
     character(2), parameter :: real4='R4'
     character(3), parameter :: int='INT'

     character*32 zpre
     logical :: fullsv,imodify
     
     ier = 0
     if(spline_o%isReady /= 1) then
        ier = 94
        return
     endif

     in0 = size(spline_o%fspl,1)
     in1 = size(spline_o%fspl,2)
     in2 = size(spline_o%fspl,3)
     in3 = size(spline_o%fspl,4)

     if(present(spl_name)) then
        call ezspline_spl_name_chk(spl_name,ier)
        if(ier.ne.0) return
        zpre=spl_name
     else
        zpre=' '
     endif

     if(present(fullsave)) then
        fullsv=fullsave
     else
        fullsv=.FALSE.
     endif

     if(zpre.eq.' ') then
        call cdfOpn(ncid, filename, 'w') ! no error flag??
        imodify=.FALSE.
     else
        call cdfOpn(ncid, filename, 'v') ! Append if exists, else create new
        imodify=.TRUE.
     endif
     if(ncid==0) then
        ier = 39
        return
     endif

     dimlens = (/1, 1, 1/) ! scalar
     call ezspline_defVar(ncid, trim(zpre)//'n1', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'n2', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'n3', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'klookup1', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'klookup2', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'klookup3', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'isHermite', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'isLinear', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'isHybrid', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'isReady', dimlens, int, imodify, ier)
     dimlens = (/3, 1, 1/) 
     call ezspline_defvar(ncid, trim(zpre)//'hspline', dimlens, int, imodify, ier)
     dimlens = (/2, 1, 1/) 
     call ezspline_defVar(ncid, trim(zpre)//'ibctype1', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'ibctype2', dimlens, int, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'ibctype3', dimlens, int, imodify, ier)
     dimlens = (/spline_o%n1, 1, 1/)
     call ezspline_defVar(ncid, trim(zpre)//'x1', dimlens, real4, imodify, ier)
     dimlens = (/spline_o%n2, 1, 1/)
     call ezspline_defVar(ncid, trim(zpre)//'x2', dimlens, real4, imodify, ier)
     dimlens = (/spline_o%n3, 1, 1/)
     call ezspline_defVar(ncid, trim(zpre)//'x3', dimlens, real4, imodify, ier)
     dimlens = (/in2, in3, 1/)
     call ezspline_defVar(ncid, trim(zpre)//'bcval1min', dimlens, real4, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'bcval1max', dimlens, real4, imodify, ier)
     dimlens = (/in1, in3, 1/)
     call ezspline_defVar(ncid, trim(zpre)//'bcval2min', dimlens, real4, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'bcval2max', dimlens, real4, imodify, ier)
     dimlens = (/in1, in2, 1/)
     call ezspline_defVar(ncid, trim(zpre)//'bcval3min', dimlens, real4, imodify, ier)
     call ezspline_defVar(ncid, trim(zpre)//'bcval3max', dimlens, real4, imodify, ier)
     !
     ! save only the function values at the nodes
     !
     dimlens = (/in1, in2, in3/)
     if(.not.fullsv) then
        call ezspline_defVar(ncid, trim(zpre)//'f', dimlens, real4, imodify, ier)
     else
        dimlens = (/1, 1, 1/) ! scalar
        call ezspline_defVar(ncid, trim(zpre)//'x1min', dimlens, real4, imodify, ier)
        call ezspline_defVar(ncid, trim(zpre)//'x1max', dimlens, real4, imodify, ier)
        call ezspline_defVar(ncid, trim(zpre)//'ilin1', dimlens, int, imodify, ier)
        call ezspline_defVar(ncid, trim(zpre)//'x2min', dimlens, real4, imodify, ier)
        call ezspline_defVar(ncid, trim(zpre)//'x2max', dimlens, real4, imodify, ier)
        call ezspline_defVar(ncid, trim(zpre)//'ilin2', dimlens, int, imodify, ier)
        call ezspline_defVar(ncid, trim(zpre)//'x3min', dimlens, real4, imodify, ier)
        call ezspline_defVar(ncid, trim(zpre)//'x3max', dimlens, real4, imodify, ier)
        call ezspline_defVar(ncid, trim(zpre)//'ilin3', dimlens, int, imodify, ier)
        dimlens = (/spline_o%n1, 4, 1/)
        call ezspline_defVar(ncid, trim(zpre)//'x1pkg', dimlens, real4, imodify, ier)
        dimlens = (/spline_o%n2, 4, 1/)
        call ezspline_defVar(ncid, trim(zpre)//'x2pkg', dimlens, real4, imodify, ier)
        dimlens = (/spline_o%n3, 4, 1/)
        call ezspline_defVar(ncid, trim(zpre)//'x3pkg', dimlens, real4, imodify, ier)
        dimlens = (/in0*in1, in2, in3/)
        call ezspline_defVar(ncid, trim(zpre)//'fspl', dimlens, real4, imodify, ier)
     endif
     if(ier.ne.0) return

     call cdfPutVar(ncid, trim(zpre)//'n1', spline_o%n1, ifail)
     call cdfPutVar(ncid, trim(zpre)//'n2', spline_o%n2, ifail)
     call cdfPutVar(ncid, trim(zpre)//'n3', spline_o%n3, ifail)
     call cdfPutVar(ncid, trim(zpre)//'klookup1', spline_o%klookup1, ifail)
     call cdfPutVar(ncid, trim(zpre)//'klookup2', spline_o%klookup2, ifail)
     call cdfPutVar(ncid, trim(zpre)//'klookup3', spline_o%klookup3, ifail)
     call cdfPutVar(ncid, trim(zpre)//'isHermite', spline_o%isHermite, ifail)
     call cdfPutVar(ncid, trim(zpre)//'isLinear', spline_o%isLinear, ifail)
     call cdfPutVar(ncid, trim(zpre)//'isHybrid', spline_o%isHybrid, ifail)
     call cdfPutVar(ncid, trim(zpre)//'isReady', spline_o%isReady, ifail)
     call cdfPutVar(ncid, trim(zpre)//'hspline', spline_o%hspline, ifail)
     call cdfPutVar(ncid, trim(zpre)//'ibctype1', spline_o%ibctype1, ifail)
     call cdfPutVar(ncid, trim(zpre)//'ibctype2', spline_o%ibctype2, ifail)
     call cdfPutVar(ncid, trim(zpre)//'ibctype3', spline_o%ibctype3, ifail)
     call cdfPutVar(ncid, trim(zpre)//'x1', spline_o%x1, ifail)
     call cdfPutVar(ncid, trim(zpre)//'x2', spline_o%x2, ifail)
     call cdfPutVar(ncid, trim(zpre)//'x3', spline_o%x3, ifail)
     call cdfPutVar(ncid, trim(zpre)//'bcval1min', spline_o%bcval1min, ifail)
     call cdfPutVar(ncid, trim(zpre)//'bcval1max', spline_o%bcval1max, ifail)
     call cdfPutVar(ncid, trim(zpre)//'bcval2min', spline_o%bcval2min, ifail)
     call cdfPutVar(ncid, trim(zpre)//'bcval2max', spline_o%bcval2max, ifail)
     call cdfPutVar(ncid, trim(zpre)//'bcval3min', spline_o%bcval3min, ifail)
     call cdfPutVar(ncid, trim(zpre)//'bcval3max', spline_o%bcval3max, ifail)
     if(ifail/=0) then
        ier=40
        return
     endif

     if(spline_o%isReady == 1) then
        if(.not.fullsv) then
           call cdfPutVar(ncid, trim(zpre)//'f', spline_o%fspl(1,:,:,:), ifail)
        else
           call cdfPutVar(ncid, trim(zpre)//'x1min', spline_o%x1min, ifail)
           call cdfPutVar(ncid, trim(zpre)//'x1max', spline_o%x1max, ifail)
           call cdfPutVar(ncid, trim(zpre)//'ilin1', spline_o%ilin1, ifail)
           call cdfPutVar(ncid, trim(zpre)//'x2min', spline_o%x2min, ifail)
           call cdfPutVar(ncid, trim(zpre)//'x2max', spline_o%x2max, ifail)
           call cdfPutVar(ncid, trim(zpre)//'ilin2', spline_o%ilin2, ifail)
           call cdfPutVar(ncid, trim(zpre)//'x3min', spline_o%x3min, ifail)
           call cdfPutVar(ncid, trim(zpre)//'x3max', spline_o%x3max, ifail)
           call cdfPutVar(ncid, trim(zpre)//'ilin3', spline_o%ilin3, ifail)
           call cdfPutVar(ncid, trim(zpre)//'x1pkg', spline_o%x1pkg, ifail)
           call cdfPutVar(ncid, trim(zpre)//'x2pkg', spline_o%x2pkg, ifail)
           call cdfPutVar(ncid, trim(zpre)//'x3pkg', spline_o%x3pkg, ifail)
           call ezspline_cdfput3(ncid, trim(zpre)//'fspl', spline_o%fspl, &
                in0*in1, in2, in3, ifail)
        endif
     else
        ier = 93
     endif

     call cdfCls(ncid)

     if(ifail /=0) ier = 17

     return
   end subroutine EZspline_save3_r4

   subroutine ezspline_spl_name_chk(spl_name,ier)
     character*(*), intent(in) :: spl_name
     integer, intent(out) :: ier

     !----------------
     integer :: ic,ilen,indx
     character*63 :: zlegal = &
          '0123456789_abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
     !----------------
     ier=0

     ilen = len(trim(spl_name))

     if(ilen.le.0) then
        ier=50
     else if(ilen.gt.20) then
        ier=51
     else
        indx=index(zlegal,spl_name(1:1))
        if(indx.le.10) then
           ier=52
        else
           do ic=2,ilen
              indx=index(zlegal,spl_name(ic:ic))
              if(indx.le.0) ier=52
           enddo
        endif
     endif

   end subroutine ezspline_spl_name_chk

   subroutine ezspline_defVar(ncid, name, dimlens, ztype, imodify, ier)

     use EZcdf
     implicit NONE

     ! inquire about variable with logic specific to ezspline_save

     integer, intent(in) :: ncid          ! NetCDF handle
     character*(*), intent(in) :: name    ! item name to define
     integer, intent(in) :: dimlens(3)    ! dimensioning information
     character*(*), intent(in) :: ztype   ! data type e.g. "R8"
     logical, intent(in) :: imodify       ! T if item may already exist

     integer, intent(inout) :: ier        ! status code, set iff a real error
     ! is detected...

     !-----------------------------
     character*10 :: ztype_loc
     integer :: dimlens_loc(3),id1,id2,ii
     integer :: ifail,ivarid
     !-----------------------------

     if(.not.imodify) then

        ! simply define the quantity...

        call cdfDefVar(ncid, name, dimlens, ztype, ifail)
        if(ifail.ne.0) ier=42

     else

        ifail = nf_inq_varid(ncid, name, ivarid)
        if(ifail.ne.0) then

           ! assume item does not exist, so, just define it...

           call cdfDefVar(ncid, name, dimlens, ztype, ifail)
           if(ifail.ne.0) ier=42

        else

           ! item DOES exist; require dimension & type match...

           call cdfInqVar(ncid, name, dimlens_loc, ztype_loc, ifail)
           if(ztype.ne.ztype_loc) then
              ier=53

           else
              do ii=1,3
                 id1=max(1,dimlens(ii))
                 id2=max(1,dimlens_loc(ii))
                 if(id1.ne.id2) ier=53
              enddo
           endif

           !  if ier= 0, no cdfDefVar is needed: object already defined
           !  and attributes are correct.

        endif
     endif

   end subroutine ezspline_defVar

!DEC$ ENDIF
