module ezspline_obj
  use ezspline_type
  interface EZspline_preInit
     !
     ! usage:  call EZspline_preinit(spline_o)
     !   ...where spline_o is a 1d, 2d, or 3d spline object of R4 or R8
     !   precision.

     module procedure &
          EZspline_preInit1_r8, &
          EZspline_preInit2_r8, &
          EZspline_preInit3_r8, &
          EZspline_preInit1_r4, &
          EZspline_preInit2_r4, &
          EZspline_preInit3_r4

  end interface

  interface EZspline_allocated
     ! logical function returns TRUE if allocated, FALSE otherwise
     !
     ! usage:  
     !   logical :: answer
     !   answer = EZspline_allocated(spline_o)
     !   ...where spline_o is a 1d, 2d, or 3d spline object of R4 or R8
     !   precision.

     module procedure &
          EZspline_allocated1_r8, &
          EZspline_allocated2_r8, &
          EZspline_allocated3_r8, &
          EZspline_allocated1_r4, &
          EZspline_allocated2_r4, &
          EZspline_allocated3_r4
  end interface

  contains

    subroutine EZspline_preInit1_r8(spline_o)
      use ezspline_type
      type(EZspline1_r8) :: spline_o
      spline_o%nguard=123456789
      spline_o%isReady=0
      spline_o%ibctype1=0 
    end subroutine EZspline_preInit1_r8

    subroutine EZspline_preInit2_r8(spline_o)
      use ezspline_type
      type(EZspline2_r8) :: spline_o
      spline_o%nguard=123456789
      spline_o%isReady=0
      spline_o%ibctype1=0 ; spline_o%ibctype2=0
    end subroutine EZspline_preInit2_r8

    subroutine EZspline_preInit3_r8(spline_o)
      use ezspline_type
      type(EZspline3_r8) spline_o
      spline_o%nguard=123456789
      spline_o%isReady=0
      spline_o%ibctype1=0 ; spline_o%ibctype2=0 ; spline_o%ibctype3=0
    end subroutine EZspline_preInit3_r8

    subroutine EZspline_preInit1_r4(spline_o)
      use ezspline_type
      type(EZspline1_r4) spline_o
      spline_o%nguard=123456789
      spline_o%isReady=0
      spline_o%ibctype1=0 
    end subroutine EZspline_preInit1_r4

    subroutine EZspline_preInit2_r4(spline_o)
      use ezspline_type
      type(EZspline2_r4) spline_o
      spline_o%nguard=123456789
      spline_o%isReady=0
      spline_o%ibctype1=0 ; spline_o%ibctype2=0
    end subroutine EZspline_preInit2_r4

    subroutine EZspline_preInit3_r4(spline_o)
      use ezspline_type
      type(EZspline3_r4) spline_o
      spline_o%nguard=123456789
      spline_o%isReady=0
      spline_o%ibctype1=0 ; spline_o%ibctype2=0 ; spline_o%ibctype3=0
    end subroutine EZspline_preInit3_r4

    logical function EZspline_allocated1_r8(spline_o)
      use ezspline_type
      type(EZspline1_r8) spline_o
      EZspline_allocated1_r8 = allocated(spline_o%fspl) &
           & .and. allocated(spline_o%x1) .and. allocated(spline_o%x1pkg) &
           & .and. (spline_o%nguard == 123456789) ! check that ezspline_init has been called
    end function EZspline_allocated1_r8

    logical function EZspline_allocated2_r8(spline_o)
      use ezspline_type
      type(EZspline2_r8) spline_o
      EZspline_allocated2_r8 = allocated(spline_o%fspl) &
           & .and. allocated(spline_o%x1) .and. allocated(spline_o%x1pkg) &
           & .and. allocated(spline_o%x2) .and. allocated(spline_o%x2pkg) &
           & .and. (spline_o%nguard == 123456789) ! check that ezspline_init has been called
    end function EZspline_allocated2_r8

    logical function EZspline_allocated3_r8(spline_o)
      use ezspline_type
      type(EZspline3_r8) spline_o
      EZspline_allocated3_r8 = allocated(spline_o%fspl) &
           & .and. allocated(spline_o%x1) .and. allocated(spline_o%x1pkg) &
           & .and. allocated(spline_o%x2) .and. allocated(spline_o%x2pkg) &
           & .and. allocated(spline_o%x3) .and. allocated(spline_o%x3pkg) &
           & .and. (spline_o%nguard == 123456789) ! check that ezspline_init has been called
    end function EZspline_allocated3_r8

    logical function EZspline_allocated1_r4(spline_o)
      use ezspline_type
      type(EZspline1_r4) spline_o
      EZspline_allocated1_r4 = allocated(spline_o%fspl) &
           & .and. allocated(spline_o%x1) .and. allocated(spline_o%x1pkg) &
           & .and. (spline_o%nguard == 123456789) ! check that ezspline_init has been called
    end function EZspline_allocated1_r4

    logical function EZspline_allocated2_r4(spline_o)
      use ezspline_type
      type(EZspline2_r4) spline_o
      EZspline_allocated2_r4 = allocated(spline_o%fspl) &
           & .and. allocated(spline_o%x1) .and. allocated(spline_o%x1pkg) &
           & .and. allocated(spline_o%x2) .and. allocated(spline_o%x2pkg) &
           & .and. (spline_o%nguard == 123456789) ! check that ezspline_init has been called
    end function EZspline_allocated2_r4

    logical function EZspline_allocated3_r4(spline_o)
      use ezspline_type
      type(EZspline3_r4) spline_o
      EZspline_allocated3_r4 = allocated(spline_o%fspl) &
           & .and. allocated(spline_o%x1) .and. allocated(spline_o%x1pkg) &
           & .and. allocated(spline_o%x2) .and. allocated(spline_o%x2pkg) &
           & .and. allocated(spline_o%x3) .and. allocated(spline_o%x3pkg) &
           & .and. (spline_o%nguard == 123456789) ! check that ezspline_init has been called
    end function EZspline_allocated3_r4

    subroutine ezmake_ict1(i,ict)
      !  (private utility for ezspline derivative2 subroutines)
      !  make ict(1:6) array
      !  for higher derivatives; d[i1+i2]f/dx[i1]dy[i2]
      !  expecting i in range [0:3] (NOT CHECKED)

      implicit NONE
      integer, intent(in) :: i
      integer, intent(out) :: ict(3)

      if(i.eq.0) then
         ict = (/1, 0, 0 /) ! seek f @ (p1)
      else if(i.eq.1) then
         ict = (/0, 1, 0 /) ! df/dx
      else if(i.eq.2) then
         ict = (/0, 0, 1 /) ! d2f/dx2
      else
         ict = (/3, 0, 0 /) ! d3f/dx3
      endif

    end subroutine ezmake_ict1

    subroutine ezmake_ict2(i1,i2,ict)
      !  (private utility for ezspline derivative2 subroutines)
      !  make ict(1:6) array
      !  for higher derivatives; d[i1+i2]f/dx[i1]dy[i2]
      !  expecting i1 & i2 in range [0:3] (NOT CHECKED)

      implicit NONE
      integer, intent(in) :: i1,i2
      integer, intent(out) :: ict(6)

      integer :: imark,isum,iii

      !  this generates the control argument needed by evbicub & similar
      !  routines...
      !----------------------

      isum = i1+i2
      ict(1)=isum

      imark=0

      if(isum.eq.0) then
         ict = (/1, 0, 0, 0, 0, 0 /) ! seek f @ (p1, p2, p3)
      else if(isum.eq.1) then
         if(i1.eq.1) then
            ict = (/0, 1, 0, 0, 0, 0 /) ! df/dx
         else
            ict = (/0, 0, 1, 0, 0, 0 /) ! df/dy
         endif
      else if(isum.eq.2) then
         if(i1.eq.2) then
            ict = (/0, 0, 0, 1, 0, 0 /) ! d2f/dx2
         else if(i2.eq.2) then
            ict = (/0, 0, 0, 0, 1, 0 /) ! d2f/dy2
         else
            ict = (/0, 0, 0, 0, 0, 1 /) ! d2f/dxdy
         endif
      else if(isum.eq.3) then
         if(i1.eq.3) then
            imark=2  ! fxxx
         else if(i1.eq.2) then
            imark=3  ! fxxy
         else if(i1.eq.1) then
            imark=4  ! fxyy
         else
            imark=5  ! fyyy
         endif
      else if(isum.eq.4) then
         if(i1.eq.3) then
            imark=2  ! fxxxy
         else if(i2.eq.3) then
            imark=4  ! fxyyy
         else
            imark=3  ! fxxyy
         endif
      else if(isum.eq.5) then
         if(i1.eq.3) then
            imark=2  ! fxxxyy
         else if(i2.eq.3) then
            imark=3  ! fxxyyy
         endif
      endif

      !  isum=6 --> fxxxyyy

      if(isum.gt.2) then
         do iii=2,6
            if(iii.eq.imark) then
               ict(iii)=1
            else
               ict(iii)=0
            endif
         enddo
      endif

    end subroutine ezmake_ict2

    subroutine ezmake_ict3(i1,i2,i3,ict)
      !  (private utility for ezspline derivative3 subroutines)
      !  make ict(1:10) array
      !  for higher derivatives; d[i1+i2+i3]f/dx[i1]dy[i2]dz[i3]
      !  i1 & i2 & i3 in range [0:3] (NOT CHECKED)

      implicit NONE
      integer, intent(in) :: i1,i2,i3
      integer, intent(out) :: ict(10)

      integer :: imark,isum,iii

      !  this generates the control argument needed by evtricub & similar
      !  routines...
      !----------------------

      isum = i1+i2+i3
      if(max(i1,i2,i3).eq.3) then
         isum=-isum
      endif
      ict(1)=isum

      imark=0

      if(isum.eq.0) then
         ict = (/1, 0, 0, 0, 0, 0, 0, 0, 0, 0 /) ! seek f @ (p1, p2, p3)

      else if(isum.eq.1) then
!  1st derivatives
         if(i1.eq.1) then
            ict = (/0, 1, 0, 0, 0, 0, 0, 0, 0, 0 /)  ! df/dx
         else if(i2.eq.1) then
            ict = (/0, 0, 1, 0, 0, 0, 0, 0, 0, 0 /)  ! df/dy
         else
            ict = (/0, 0, 0, 1, 0, 0, 0, 0, 0, 0 /)  ! df/dz
         endif

      else if(isum.eq.2) then
!  2nd derivatives -- legacy ordering; x-precedence ordering for all
!  higher derivatives...

         if(i1.eq.2) then
            ict = (/0, 0, 0, 0, 1, 0, 0, 0, 0, 0 /)  ! d2f/dx2
         else if(i2.eq.2) then
            ict = (/0, 0, 0, 0, 0, 1, 0, 0, 0, 0 /)  ! d2f/dy2
         else if(i3.eq.2) then
            ict = (/0, 0, 0, 0, 0, 0, 1, 0, 0, 0 /)  ! d2f/dz2
         else if(i3.eq.0) then
            ict = (/0, 0, 0, 0, 0, 0, 0, 1, 0, 0 /)  ! d2f/dxdy
         else if(i2.eq.0) then
            ict = (/0, 0, 0, 0, 0, 0, 0, 0, 1, 0 /)  ! d2f/dxdz
         else
            ict = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 1 /)  ! d2f/dydz
         endif
        
      else if(isum.eq.3) then
!  3rd derivative, continuous: max(i1,i2,i3)<3
         if(i1.eq.2) then
            if(i2.eq.1) then
               imark=2     ! fxxy
            else
               imark=3     ! fxxz
            endif
         else if(i1.eq.1) then
            if(i2.eq.2) then
               imark=4     ! fxyy
            else if(i2.eq.1) then
               imark=5     ! fxyz
            else
               imark=6     ! fxzz
            endif
         else
            if(i2.eq.2) then
               imark=7     ! fyyz
            else
               imark=8     ! fyzz
            endif
         endif

      else if(isum.eq.-3) then
!  3rd derivative
         if(i1.eq.3) then
            imark=2        ! fxxx
         else if(i2.eq.3) then
            imark=3        ! fyyy
         else if(i3.eq.3) then
            imark=4        ! fzzz
         endif

      else if(isum.eq.4) then
!  4th derivative, continuous: max(i1,i2,i3)<3
         if(i1.eq.2) then
            if(i2.eq.2) then
               imark=2     ! fxxyy
            else if(i2.eq.1) then
               imark=3     ! fxxyz
            else
               imark=4     ! fxxzz
            endif
         else if(i1.eq.1) then
            if(i2.eq.2) then
               imark=5     ! fxyyz
            else
               imark=6     ! fxyzz
            endif
         else
            imark=7        ! fyyzz
         endif

      else if(isum.eq.-4) then
!  4th derivative
         if(i1.eq.3) then
            if(i2.eq.1) then
               imark=2     ! fxxxy
            else
               imark=3     ! fxxxz
            endif
         else if(i1.eq.1) then
            if(i2.eq.3) then
               imark=4     ! fxyyy
            else
               imark=5     ! fxzzz
            endif
         else
            if(i2.eq.3) then
               imark=6     ! fyyyz
            else
               imark=7     ! fyzzz
            endif
         endif

      else if(isum.eq.5) then
!  5th derivative, continuous: max(i1,i2,i3)<3
         if(i3.eq.1) then
            imark=2     ! fxxyyz
         else if(i2.eq.1) then
            imark=3     ! fxxyzz
         else
            imark=4     ! fxyyzz
         endif

      else if(isum.eq.-5) then
!  5th derivative
         if(i1.eq.3) then
            if(i2.eq.2) then
               imark=2  ! fxxxyy
            else if(i2.eq.1) then
               imark=3  ! fxxxyz
            else
               imark=4  ! fxxxzz
            endif
         else if(i1.eq.2) then
            if(i2.eq.3) then
               imark=5  ! fxxyyy
            else
               imark=6  ! fxxzzz
            endif
         else if(i1.eq.1) then
            if(i2.eq.3) then
               imark=7  ! fxyyyz
            else
               imark=8  ! fxyzzz
            endif
         else
            if(i2.eq.3) then
               imark=9  ! fyyyzz
            else
               imark=10 ! fyyzzz
            endif
         endif

!  isum=6 --> fxxyyzz  (i1=i2=i3=2)
      else if(isum.eq.-6) then
!  6th derivative
         if(i1.eq.3) then
            if(i2.eq.3) then
               imark=2  ! fxxxyyy
            else if(i2.eq.2) then
               imark=3  ! fxxxyyz
            else if(i2.eq.1) then
               imark=4  ! fxxxyzz
            else
               imark=5  ! fxxxzzz
            endif
         else if(i1.eq.2) then
            if(i2.eq.3) then
               imark=6  ! fxxyyyz
            else if(i2.eq.1) then
               imark=7  ! fxxyzzz
            endif
         else if(i1.eq.1) then
            if(i2.eq.3) then
               imark=8  ! fxyyyzz
            else
               imark=9  ! fxyyzzz
            endif
         else
            imark=10    ! fyyyzzz
         endif

!  isum=7 not possible
      else if(isum.eq.-7) then
!  7th derivative
         if(i1.eq.3) then
            if(i2.eq.3) then
               imark=2  ! fxxxyyyz
            else if(i2.eq.2) then
               imark=3  ! fxxxyyzz
            else
               imark=4  ! fxxxyzzz
            endif
         else if(i1.eq.2) then
            if(i2.eq.3) then
               imark=5  ! fxxyyyzz
            else
               imark=6  ! fxxyyzzz
            endif
         else
            imark=7     ! fxyyyzzz
         endif

!  isum=8 not possible
      else if(isum.eq.-8) then
!  8th derivative
         if(i3.eq.2) then 
            imark=2  ! fxxxyyyzz
         else if(i2.eq.2) then
            imark=3  ! fxxxyyzzz
         else
            imark=4  ! fxxyyyzzz
         endif

!  isum=9 not possible
!  isum=-9 --> fxxxyyyzzz

      endif

      if(abs(isum).gt.2) then
         do iii=2,10
            if(iii.eq.imark) then
               ict(iii)=1
            else
               ict(iii)=0
            endif
         enddo
      endif

    end subroutine ezmake_ict3

end module EZspline_obj
