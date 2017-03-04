!/////
! R8 !
!/////
   subroutine EZspline_isGridRegular1_r8(spline_o, ier)
     use EZspline_obj
     implicit none
     type(EZspline1_r8) :: spline_o
     ! ier:
     ! 14=x1 grid is not strictly increasing
     integer, intent(out) :: ier
     integer i
     
     ier = 0
     do i=1, spline_o%n1-1
        if(spline_o%x1(i+1) <= spline_o%x1(i)) then
           ier = 14
           return
        endif
     enddo

   end subroutine EZspline_isGridRegular1_r8

   subroutine EZspline_isGridRegular2_r8(spline_o, ier)
     use EZspline_obj
     implicit none
     type(EZspline2_r8) :: spline_o
     ! ier:
     ! 14=x1 grid is not strictly increasing
     ! 15=x2 grid is not strictly increasing
     integer, intent(out) :: ier
     integer i
     
     ier = 0
     do i=1, spline_o%n1-1
        if(spline_o%x1(i+1) <= spline_o%x1(i)) then
           ier = 14
           return
        endif
     enddo
     do i=1, spline_o%n2-1
        if(spline_o%x2(i+1) <= spline_o%x2(i)) then
           ier = 15
           return
        endif
     enddo

   end subroutine EZspline_isGridRegular2_r8


   subroutine EZspline_isGridRegular3_r8(spline_o, ier)
     use EZspline_obj
     implicit none
     type(EZspline3_r8) :: spline_o
     ! ier:
     ! 14=x1 grid is not strictly increasing
     ! 15=x2 grid is not strictly increasing
     ! 16=x3 grid is not strictly increasing
     integer, intent(out) :: ier
     integer i
     
     ier = 0
     do i=1, spline_o%n1-1
        if(spline_o%x1(i+1) <= spline_o%x1(i)) then
           ier = 14
           return
        endif
     enddo
     do i=1, spline_o%n2-1
        if(spline_o%x2(i+1) <= spline_o%x2(i)) then
           ier = 15
           return
        endif
     enddo
     do i=1, spline_o%n3-1
        if(spline_o%x3(i+1) <= spline_o%x3(i)) then
           ier = 16
           return
        endif
     enddo

   end subroutine EZspline_isGridRegular3_r8
!/////
! R4 !
!/////
   subroutine EZspline_isGridRegular1_r4(spline_o, ier)
     use EZspline_obj
     implicit none
     type(EZspline1_r4) :: spline_o
     ! ier:
     ! 14=x1 grid is not strictly increasing
     integer, intent(out) :: ier
     integer i
     
     ier = 0
     do i=1, spline_o%n1-1
        if(spline_o%x1(i+1) <= spline_o%x1(i)) then
           ier = 14
           return
        endif
     enddo

   end subroutine EZspline_isGridRegular1_r4

   subroutine EZspline_isGridRegular2_r4(spline_o, ier)
     use EZspline_obj
     implicit none
     type(EZspline2_r4) :: spline_o
     ! ier:
     ! 14=x1 grid is not strictly increasing
     ! 15=x2 grid is not strictly increasing
     integer, intent(out) :: ier
     integer i
     
     ier = 0
     do i=1, spline_o%n1-1
        if(spline_o%x1(i+1) <= spline_o%x1(i)) then
           ier = 14
           return
        endif
     enddo
     do i=1, spline_o%n2-1
        if(spline_o%x2(i+1) <= spline_o%x2(i)) then
           ier = 15
           return
        endif
     enddo

   end subroutine EZspline_isGridRegular2_r4


   subroutine EZspline_isGridRegular3_r4(spline_o, ier)
     use EZspline_obj
     implicit none
     type(EZspline3_r4) :: spline_o
     ! ier:
     ! 14=x1 grid is not strictly increasing
     ! 15=x2 grid is not strictly increasing
     ! 16=x3 grid is not strictly increasing
     integer, intent(out) :: ier
     integer i
     
     ier = 0
     do i=1, spline_o%n1-1
        if(spline_o%x1(i+1) <= spline_o%x1(i)) then
           ier = 14
           return
        endif
     enddo
     do i=1, spline_o%n2-1
        if(spline_o%x2(i+1) <= spline_o%x2(i)) then
           ier = 15
           return
        endif
     enddo
     do i=1, spline_o%n3-1
        if(spline_o%x3(i+1) <= spline_o%x3(i)) then
           ier = 16
           return
        endif
     enddo

   end subroutine EZspline_isGridRegular3_r4
