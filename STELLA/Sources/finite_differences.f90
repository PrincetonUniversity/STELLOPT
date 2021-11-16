module finite_differences

  implicit none

  public :: first_order_upwind
  public :: third_order_upwind
  public :: fifth_order_upwind
  public :: third_order_upwind_zed
  public :: first_order_upwind_zed
  public :: second_order_centered
  public :: fourth_order_centered
  public :: second_order_centered_zed
  public :: four_point_triangle
  public :: fd3pt, fd5pt
  public :: d2_3pt
  public :: fd_variable_upwinding_vpa
  public :: fd_variable_upwinding_zed
  public :: fd_cell_centres_zed, cell_centres_zed

  interface fd3pt
     module procedure fd3pt_real
     module procedure fd3pt_real_array
     module procedure fd3pt_complex_array
  end interface

  interface fd5pt
     module procedure fd5pt_real
     module procedure fd5pt_array
  end interface

  interface first_order_upwind
     module procedure first_order_upwind_real
     module procedure first_order_upwind_complex
  end interface

  interface third_order_upwind
     module procedure third_order_upwind_complex
     module procedure third_order_upwind_real
  end interface

  interface fifth_order_upwind
     module procedure fifth_order_upwind_complex
     module procedure fifth_order_upwind_real
  end interface

  interface tridag
     module procedure tridag_real
     module procedure tridag_complex
  end interface

  interface second_order_centered
     module procedure second_order_centered_real
     module procedure second_order_centered_complex
  end interface

  interface four_point_triangle
     module procedure  four_point_triangle_real
     module procedure  four_point_triangle_complex
  end interface

  interface fourth_order_centered
     module procedure fourth_order_centered_real
     module procedure fourth_order_centered_complex
  end interface

  interface second_order_centered_zed
     module procedure second_order_centered_zed_real
     module procedure second_order_centered_zed_complex
  end interface

  interface d2_3pt
     module procedure d2_3pt_real
     module procedure d2_3pt_complex
  end interface

contains

  subroutine first_order_upwind_real (llim, f, del, sgn, df)
    
    implicit none
    
    integer, intent (in) :: llim
    real, dimension (llim:), intent (in) :: f
    real, intent (in) :: del
    integer, intent (in) :: sgn
    real, dimension (llim:), intent (out) :: df
    
    integer :: i, n, istart, iend

    n = size(f)

    if (sgn == -1) then
       istart = llim
       iend = llim+n-1
    else
       istart = llim+n-1
       iend = llim
    end if

    ! zero BC, 1st order accurate upwind
    df(istart) = -f(istart)*sgn/del
    do i = istart-sgn, iend, -sgn
       df(i) = sgn*(f(i+sgn)-f(i))/del
    end do

  end subroutine first_order_upwind_real
  
  subroutine first_order_upwind_complex (llim, f, del, sgn, df)
    
    implicit none
    
    integer, intent (in) :: llim
    complex, dimension (llim:), intent (in) :: f
    real, intent (in) :: del
    integer, intent (in) :: sgn
    complex, dimension (llim:), intent (out) :: df
    
    integer :: i, n, istart, iend

    n = size(f)

    if (sgn == -1) then
       istart = llim
       iend = llim+n-1
    else
       istart = llim+n-1
       iend = llim
    end if

    ! zero BC, 1st order accurate upwind
    df(istart) = -f(istart)*sgn/del
    do i = istart-sgn, iend, -sgn
       df(i) = sgn*(f(i+sgn)-f(i))/del
    end do

  end subroutine first_order_upwind_complex
  
  subroutine third_order_upwind_complex (llim, f, del, sgn, df)
    
    implicit none
    
    integer, intent (in) :: llim
    complex, dimension (llim:), intent (in) :: f
    real, intent (in) :: del
    integer, intent (in) :: sgn
    complex, dimension (llim:), intent (out) :: df
    
    integer :: i, n, istart, iend

    n = size(f)
    if (sgn == -1) then
       istart = llim
       iend = llim+n-1
    else
       istart = llim+n-1
       iend = llim
    end if

    i = istart-sgn
    ! zero BC, 1st order accurate upwind
    df(istart) = -f(istart)*sgn/del

    ! zero BC, 3rd order accurate upwind
    df(i) = -sgn*(2.*f(i-sgn)+3.*f(i)-6.*f(i+sgn))/(6.*del)

    ! 1st order accurate upwind
    df(iend) = sgn*(f(iend+sgn)-f(iend))/del

    ! 3rd order accurate upwind
    do i = istart-2*sgn, iend+sgn, -sgn
       df(i) = -sgn*(2.*f(i-sgn)+3*f(i)-6.*f(i+sgn)+f(i+2*sgn))/(6.*del)
    end do

  end subroutine third_order_upwind_complex

  subroutine third_order_upwind_real (llim, f, del, sgn, df)
    
    implicit none
    
    integer, intent (in) :: llim
    real, dimension (llim:), intent (in) :: f
    real, intent (in) :: del
    integer, intent (in) :: sgn
    real, dimension (llim:), intent (out) :: df
    
    integer :: i, n, istart, iend

    n = size(f)
    if (sgn == -1) then
       istart = llim
       iend = llim+n-1
    else
       istart = llim+n-1
       iend = llim
    end if

    i = istart-sgn
    ! zero BC, 1st order accurate upwind
    df(istart) = -f(istart)*sgn/del
    ! zero BC, 3rd order accurate upwind
    df(i) = -sgn*(2.*f(i-sgn)+3.*f(i)-6.*f(i+sgn))/(6.*del)

    ! 1st order accurate upwind
    df(iend) = sgn*(f(iend+sgn)-f(iend))/del

    ! 3rd order accurate upwind
    do i = istart-2*sgn, iend+sgn, -sgn
       df(i) = -sgn*(2.*f(i-sgn)+3*f(i)-6.*f(i+sgn)+f(i+2*sgn))/(6.*del)
    end do

  end subroutine third_order_upwind_real

  subroutine fifth_order_upwind_complex (llim, f, del, sgn, df)
    
    implicit none
    
    integer, intent (in) :: llim
    complex, dimension (llim:), intent (in) :: f
    real, intent (in) :: del
    integer, intent (in) :: sgn
    complex, dimension (llim:), intent (out) :: df
    
    integer :: i, n, istart, iend

    n = size(f)
    if (sgn == -1) then
       istart = llim
       iend = llim+n-1
    else
       istart = llim+n-1
       iend = llim
    end if

    ! zero BC, 1st order accurate upwind
    df(istart) = -f(istart)*sgn/del
    ! zero BC, 3rd order accurate upwind
    i = istart-sgn
    df(i) = -sgn*(2.*f(i-sgn)+3.*f(i)-6.*f(i+sgn))/(6.*del)
    ! zero BC, 5th order accurate upwind
    i = istart-2*sgn
    df(i) = -sgn*(-3.*f(i-2*sgn)+30.*f(i-sgn)+20.*f(i)-60.*f(i+sgn)+15.*f(i+2*sgn))/(60.*del)

    ! 1st order accurate upwind
    df(iend) = sgn*(f(iend+sgn)-f(iend))/del
    ! 3rd order accurate upwind
    df(iend+sgn) = -sgn*(2.*f(iend)+3*f(iend+sgn)-6.*f(iend+2*sgn)+f(iend+3*sgn))/(6.*del)

    ! 5th order accurate upwind
    do i = istart-3*sgn, iend+2*sgn, -sgn
       df(i) = -sgn*(-3.*f(i-2*sgn)+30.*f(i-sgn)+20.*f(i)-60.*f(i+sgn)+15.*f(i+2*sgn)-2.*f(i+3*sgn)) &
                /(60.*del)
    end do

  end subroutine fifth_order_upwind_complex

  subroutine fifth_order_upwind_real (llim, f, del, sgn, df)
    
    implicit none
    
    integer, intent (in) :: llim
    real, dimension (llim:), intent (in) :: f
    real, intent (in) :: del
    integer, intent (in) :: sgn
    real, dimension (llim:), intent (out) :: df
    
    integer :: i, n, istart, iend

    n = size(f)
    if (sgn == -1) then
       istart = llim
       iend = llim+n-1
    else
       istart = llim+n-1
       iend = llim
    end if

    ! zero BC, 1st order accurate upwind
    df(istart) = -f(istart)*sgn/del
    ! zero BC, 3rd order accurate upwind
    i = istart-sgn
    df(i) = -sgn*(2.*f(i-sgn)+3.*f(i)-6.*f(i+sgn))/(6.*del)
    ! zero BC, 5th order accurate upwind
    i = istart-2*sgn
    df(i) = -sgn*(-3.*f(i-2*sgn)+30.*f(i-sgn)+20.*f(i)-60.*f(i+sgn)+15.*f(i+2*sgn))/(60.*del)

    ! 1st order accurate upwind
    df(iend) = -sgn*(f(iend)-f(iend+sgn))/del
    ! 3rd order accurate upwind
    df(iend+sgn) = -sgn*(2.*f(iend)+3*f(iend+sgn)-6.*f(iend+2*sgn)+f(iend+3*sgn))/(6.*del)

    ! 5th order accurate upwind
    do i = istart-3*sgn, iend+2*sgn, -sgn
       df(i) = -sgn*(-3.*f(i-2*sgn)+30.*f(i-sgn)+20.*f(i)-60.*f(i+sgn)+15.*f(i+2*sgn)-2.*f(i+3*sgn)) &
                /(60.*del)
    end do

  end subroutine fifth_order_upwind_real

  subroutine third_order_upwind_zed (llim, iseg, nseg, f, del, sgn, fl, fr, periodic, df)

    implicit none
    
    integer, intent (in) :: llim, iseg, nseg
    complex, dimension (llim:), intent (in) :: f
    real, intent (in) :: del
    integer, intent (in) :: sgn
    logical, intent (in) :: periodic
    complex, dimension (:), intent (in) :: fl, fr
    complex, dimension (llim:), intent (out) :: df

    integer :: i, istart, iend, ulim

    ulim = size(f)+llim-1    
    ! if sgn > 0, then stream speed is negative
    ! so sweep from more positive to more negative zed
    if (sgn > 0) then
       if (iseg == nseg .and..not.periodic) then
          i = ulim
          df(i) = -f(i)/del
          i = ulim-1
          df(i) = -(2.*f(i-1)+3.*f(i)-6.*f(i+1))/(6.*del)
       else
          i = ulim
          df(i) = -(2.*f(i-1)+3.*f(i)-6.*fr(1)+fr(2))/(6.*del)
          i = ulim-1
          df(i) = -(2.*f(i-1)+3.*f(i)-6.*f(i+1)+fr(1))/(6.*del)
       end if
       if (iseg == 1.and..not.periodic) then
          i = llim
          df(i) = (f(i+1)-f(i))/del
       else
          i = llim
          df(i) = -(2.*fl(2)+3*f(i)-6.*f(i+1)+f(i+2))/(6.*del)
       end if
       istart = ulim
       iend = llim
    else
       if (iseg == 1.and..not.periodic) then
          i = llim
          df(i) = f(i)/del
          i = llim+1
          df(i) = (2.*f(i+1)+3.*f(i)-6.*f(i-1))/(6.*del)
       else
          i = llim
          df(i) = (2.*f(i+1)+3*f(i)-6.*fl(2)+fl(1))/(6.*del)
          i = llim+1
          df(i) = (2.*f(i+1)+3*f(i)-6.*f(i-1)+fl(2))/(6.*del)
       end if
       if (iseg == nseg.and..not.periodic) then
          i = ulim
          df(i) = (f(i)-f(i-1))/del
       else
          i = ulim
          df(i) = (2.*fr(1)+3*f(i)-6.*f(i-1)+f(i-2))/(6.*del)
       end if
       istart = llim
       iend = ulim
    end if

    ! 3rd order accurate upwind
    do i = istart-2*sgn, iend+sgn, -sgn
       df(i) = -sgn*(2.*f(i-sgn)+3*f(i)-6.*f(i+sgn)+f(i+2*sgn))/(6.*del)
    end do

  end subroutine third_order_upwind_zed

  subroutine first_order_upwind_zed (llim, iseg, nseg, f, del, sgn, fl, fr, periodic, df)

    implicit none
    
    integer, intent (in) :: llim, iseg, nseg
    complex, dimension (llim:), intent (in) :: f
    real, intent (in) :: del
    integer, intent (in) :: sgn
    logical, intent (in) :: periodic
    complex, dimension (:), intent (in) :: fl, fr
    complex, dimension (llim:), intent (out) :: df

    integer :: i, istart, iend, ulim

    ulim = size(f)+llim-1    
    ! if sgn > 0, then stream speed is negative
    ! so sweep from more positive to more negative zed
    if (sgn > 0) then
       if (iseg == nseg.and..not.periodic) then
          i = ulim
          df(i) = -f(i)/del
          i = ulim-1
          df(i) = (f(i+1)-f(i))/del
       else
          i = ulim
          df(i) = (fr(1)-f(i))/del
          i = ulim-1
          df(i) = (f(i+1)-f(i))/del
       end if
       i = llim
       df(i) = (f(i+1)-f(i))/del
       istart = ulim
       iend = llim
    else
       if (iseg == 1.and..not.periodic) then
          i = llim
          df(i) = f(i)/del
          i = llim+1
          df(i) = (f(i)-f(i-1))/del
       else
          i = llim
          df(i) = (f(i)-fl(2))/del
          i = llim+1
          df(i) = (f(i)-f(i-1))/del
       end if
       i = ulim
       df(i) = (f(i)-f(i-1))/del
       istart = llim
       iend = ulim
    end if

    ! 3rd order accurate upwind
    do i = istart-2*sgn, iend+sgn, -sgn
       df(i) = sgn*(f(i+sgn)-f(i))/del
    end do

  end subroutine first_order_upwind_zed

  subroutine second_order_centered_real (llim, f, del, df)
    
    implicit none
    
    integer, intent (in) :: llim
    real, dimension (llim:), intent (in) :: f
    real, intent (in) :: del
    real, dimension (llim:), intent (out) :: df
    
    integer :: i, n, istart, iend

    n = size(f)
    istart = llim
    iend = llim+n-1

    ! zero BC
    df(istart) = f(istart+1) / (2.*del)
    df(iend)   =-f(iend-1)   / (2.*del) 


    ! 2nd order accurate centered
    do i = istart+1, iend-1
       df(i) = (f(i+1) - f(i-1)) / (2.*del)
    end do

  end subroutine second_order_centered_real

  subroutine second_order_centered_complex (llim, f, del, df)
    
    implicit none
    
    integer, intent (in) :: llim
    complex, dimension (llim:), intent (in) :: f
    real, intent (in) :: del
    complex, dimension (llim:), intent (out) :: df
    
    integer :: i, n, istart, iend

    n = size(f)
    istart = llim
    iend = llim+n-1

    ! zero BC
    df(istart) = f(istart+1) / (2.*del)
    df(iend)   =-f(iend-1)   / (2.*del) 


    ! 2nd order accurate centered
    do i = istart+1, iend-1
       df(i) = (f(i+1) - f(i-1)) / (2.*del)
    end do

  end subroutine second_order_centered_complex

  subroutine four_point_triangle_real (llim, f, del, df)
    
    implicit none
    
    integer, intent (in) :: llim
    real, dimension (llim:), intent (in) :: f
    real, intent (in) :: del
    real, dimension (llim:), intent (out) :: df
    
    integer :: i, n, istart, iend

    n = size(f)
    istart = llim
    iend = llim+n-1

    ! zero BC
    i=istart
    df(i) = f(i+1) / (2.0*del)
    i=istart+1
    df(i) = (f(i+1)-f(i-1)) / (2.0*del)
    i=istart+2
    df(i) = (-2.*f(i+3)+9.*f(i+1)-9.*f(i-1)) / (18.0*del)

    i=iend
    df(i) = -f(i-1) / (2.0*del)
    i=iend-1
    df(i) = (f(i+1)-f(i-1)) / (2.0*del)
    i=iend-2
    df(i) = (9.*f(i+1)-9.*f(i-1)+2.*f(i-3)) / (18.0*del)


    ! 2nd order accurate centered
    do i = istart+3, iend-3
      df(i) = (-2.*f(i+3)+9.*f(i+1)-9.*f(i-1)+2.*f(i-3)) / (18.0*del)
    end do

  end subroutine four_point_triangle_real

  subroutine four_point_triangle_complex (llim, f, del, df)
    
    implicit none
    
    integer, intent (in) :: llim
    complex, dimension (llim:), intent (in) :: f
    real, intent (in) :: del
    complex, dimension (llim:), intent (out) :: df
    
    integer :: i, n, istart, iend

    n = size(f)
    istart = llim
    iend = llim+n-1

    ! zero BC
    i=istart
    df(i) = f(i+1) / (2.0*del)
    i=istart+1
    df(i) = (f(i+1)-f(i-1)) / (2.0*del)
    i=istart+2
    df(i) = (-2.*f(i+3)+9.*f(i+1)-9.*f(i-1)) / (18.0*del)

    i=iend
    df(i) = -f(i-1) / (2.0*del)
    i=iend-1
    df(i) = (f(i+1)-f(i-1)) / (2.0*del)
    i=iend-2
    df(i) = (9.*f(i+1)-9.*f(i-1)+2.*f(i-3)) / (18.0*del)


    ! 2nd order accurate centered
    do i = istart+3, iend-3
      df(i) = (-2.*f(i+3)+9.*f(i+1)-9.*f(i-1)+2.*f(i-3)) / (18.0*del)
    end do

  end subroutine four_point_triangle_complex

  subroutine fourth_order_centered_real (llim, f, del, df)
    
    implicit none
    
    integer, intent (in) :: llim
    real, dimension (llim:), intent (in) :: f
    real, intent (in) :: del
    real, dimension (llim:), intent (out) :: df
    
    integer :: i, n, istart, iend

    n = size(f)
    istart = llim
    iend = llim+n-1

    ! zero BC
    ! 2nd order accurate centered
    df(istart) = f(istart + 1) / (2.*del)
    df(iend)   =-f(iend   - 1) / (2.*del) 

    ! 4th order accurate centered
    df(istart+1) = (f(istart+3) - 8.*f(istart+2) + 8.*f(istart)) / (12.*del)
    df(iend-1)   = (-8*f(iend)  + 8.*f(iend-2) - f(iend-3))      / (12.*del) 


    ! 4th order accurate centered
    do i = istart+2, iend-2
       df(i) = (f(i+2)-8.*f(i+1) + 8.*f(i-1)-f(i-2)) / (12.*del)
    end do

  end subroutine fourth_order_centered_real

  subroutine fourth_order_centered_complex (llim, f, del, df)
    
    implicit none
    
    integer, intent (in) :: llim
    complex, dimension (llim:), intent (in) :: f
    real, intent (in) :: del
    complex, dimension (llim:), intent (out) :: df
    
    integer :: i, n, istart, iend

    n = size(f)
    istart = llim
    iend = llim+n-1

    ! zero BC
    ! 2nd order accurate centered
    df(istart) = f(istart + 1) / (2.*del)
    df(iend)   =-f(iend   - 1) / (2.*del) 

    ! 4th order accurate centered
    df(istart+1) = (f(istart+3) - 8.*f(istart+2) + 8.*f(istart))/ (12.*del)
    df(iend-1)   = (-8*f(iend) + 8.*f(iend-2) - f(iend-3))      / (12.*del) 


    ! 4th order accurate centered
    do i = istart+2, iend-2
       df(i) = (f(i+2)-8.*f(i+1) + 8.*f(i-1)-f(i-2)) / (12.*del)
    end do

  end subroutine fourth_order_centered_complex

  subroutine second_order_centered_zed_real (llim, iseg, nseg, f, del, sgn, fl, fr, periodic, df)

    implicit none
    
    integer, intent (in) :: llim, iseg, nseg
    real, dimension (llim:), intent (in) :: f
    integer, intent (in) :: sgn
    real, intent (in) :: del
    real, dimension (:), intent (in) :: fl, fr
    logical, intent (in) :: periodic
    real, dimension (llim:), intent (out) :: df

    integer :: i, ulim

    ulim = size(f)+llim-1    

    i = llim
    if (iseg == 1 .and. sgn>0 .and..not.periodic) then
       ! sgn > 0 corresponds to negative advection speed
       ! upwind at boundary requires taking information from right
       df(i) = (f(i+1)-f(i))/del
    else
       df(i) = 0.5*(f(i+1)-fl(2))/del
    end if

    i = ulim
    if (iseg == nseg .and. sgn<0 .and..not.periodic) then
       ! sgn < 0 corresponds to positive advection speed
       ! upwind at boundary requires taking information from left
       df(i) = (f(i)-f(i-1))/del
    else
       df(i) = 0.5*(fr(1)-f(i-1))/del
    end if
    
    do i = llim+1, ulim-1
       df(i) = 0.5*(f(i+1)-f(i-1))/del
    end do

  end subroutine second_order_centered_zed_real

  subroutine second_order_centered_zed_complex (llim, iseg, nseg, f, del, sgn, fl, fr, periodic, df)

    implicit none
    
    integer, intent (in) :: llim, iseg, nseg
    complex, dimension (llim:), intent (in) :: f
    integer, intent (in) :: sgn
    real, intent (in) :: del
    complex, dimension (:), intent (in) :: fl, fr
    logical, intent (in) :: periodic
    complex, dimension (llim:), intent (out) :: df

    integer :: i, ulim

    ulim = size(f)+llim-1    

    i = llim
    if (iseg == 1 .and. sgn>0 .and..not.periodic) then
       ! sgn > 0 corresponds to negative advection speed
       ! upwind at boundary requires taking information from right
       df(i) = (f(i+1)-f(i))/del
    else
       df(i) = 0.5*(f(i+1)-fl(2))/del
    end if

    i = ulim
    if (iseg == nseg .and. sgn<0 .and..not.periodic) then
       ! sgn < 0 corresponds to positive advection speed
       ! upwind at boundary requires taking information from left
       df(i) = (f(i)-f(i-1))/del
    else
       df(i) = 0.5*(fr(1)-f(i-1))/del
    end if
    
    do i = llim+1, ulim-1
       df(i) = 0.5*(f(i+1)-f(i-1))/del
    end do

  end subroutine second_order_centered_zed_complex

  subroutine second_order_centered_vpa (llim, f, del, df)

    implicit none
    
    integer, intent (in) :: llim
    complex, dimension (llim:), intent (in) :: f
    real, intent (in) :: del
    complex, dimension (llim:), intent (out) :: df

    integer :: i, ulim

    ulim = size(f)+llim-1    

    i = llim
    df(i) = 0.5*f(i+1)/del

    i = ulim
    df(i) = -0.5*f(i-1)/del

    do i = llim+1, ulim-1
       df(i) = 0.5*(f(i+1)-f(i-1))/del
    end do

  end subroutine second_order_centered_vpa

  subroutine fd_cell_centres_zed (llim, f, del, sgn, fl, fr, df)

    implicit none

    integer, intent (in) :: llim
    complex, dimension (llim:), intent (in) :: f
    real, intent (in) :: del
    integer, intent (in) :: sgn
    complex, intent (in) :: fl, fr
    complex, dimension (llim:), intent (out) :: df

    integer :: i, ulim

    ulim = size(f)+llim-1
       
    if (sgn > 0) then
       ! if sgn > 0, then stream speed is negative
       ! so sweep from more positive to more negative zed
       i = ulim
       df(i) = (fr-f(i))/del
       do i = ulim-1, llim, -1
          df(i) = (f(i+1)-f(i))/del
       end do
    else
       ! if sgn < 0, then stream speed is positive
       ! so sweep from more negative to more positive zed
       i = llim
       df(i) = (f(i)-fl)/del
       do i = llim+1, ulim
          df(i) = (f(i)-f(i-1))/del
       end do
    end if

  end subroutine fd_cell_centres_zed

  ! cell_centres_zed takes f at z grid locations
  ! and returns f at cell centres
  ! (with possible offset due to upwinding)
  subroutine cell_centres_zed (llim, f, upwnd, sgn, fl, fr, fc)

    implicit none

    integer, intent (in) :: llim
    complex, dimension (llim:), intent (in) :: f
    real, intent (in) :: upwnd
    integer, intent (in) :: sgn
    complex, intent (in) :: fl, fr
    complex, dimension (llim:), intent (out) :: fc

    integer :: i, ulim

    ulim = size(f)+llim-1

    if (sgn > 0) then
       ! if sgn > 0, then stream speed is negative
       ! so sweep from more positive to more negative zed
       i = ulim
       fc(i) = 0.5*((1.-upwnd)*fr + (1.+upwnd)*f(i))
       do i = ulim-1, llim, -1
          fc(i) = 0.5*((1.-upwnd)*f(i+1) + (1.+upwnd)*f(i))
       end do
    else
       ! if sgn < 0, then stream speed is positive
       ! so sweep from more negative to more positive zed
       i = llim
       fc(i) = 0.5*((1.+upwnd)*f(i)+(1.-upwnd)*fl)
       do i = llim+1, ulim
          fc(i) = 0.5*((1.+upwnd)*f(i)+(1.-upwnd)*f(i-1))
       end do
    end if

  end subroutine cell_centres_zed

  subroutine fd_variable_upwinding_zed (llim, iseg, nseg, f, del, sgn, upwnd, fl, fr, periodic, df)

    implicit none

    integer, intent (in) :: llim, iseg, nseg
    complex, dimension (llim:), intent (in) :: f
    real, intent (in) :: del, upwnd
    integer, intent (in) :: sgn
    complex, dimension (:), intent (in) :: fl, fr
    logical, intent (in) :: periodic
    complex, dimension (llim:), intent (out) :: df

    integer :: i, istart, iend, ulim

    ! if upwnd is zero or if vpa=0, then use centered differences
    if (abs(upwnd) < epsilon(0.) .or. sgn == 0) then
       call second_order_centered_zed (llim, iseg, nseg, f, del, sgn, fl, fr, periodic, df)
    else
       ulim = size(f)+llim-1
       
       ! if sgn > 0, then stream speed is negative
       ! so sweep from more positive to more negative zed
       if (sgn > 0) then
          if (iseg == nseg.and..not.periodic) then
             i = ulim
             df(i) = (0.5*(upwnd-1.)*f(i-1)-upwnd*f(i))/del
          else
             i = ulim
             df(i) = (0.5*(upwnd-1.)*f(i-1)-upwnd*f(i)+0.5*(1.+upwnd)*fr(1))/del
          end if
          if (iseg == 1.and..not.periodic) then
             i = llim
             ! at left boundary, must upwind fully as no info for f(i-1)
             df(i) = (f(i+1)-f(i))/del
          else
             i = llim
             df(i) = (0.5*(1.+upwnd)*f(i+1)-upwnd*f(i)+0.5*(upwnd-1.)*fl(2))/del
          end if
          istart = ulim
          iend = llim
       else
          if (iseg == 1.and..not.periodic) then
             i = llim
             df(i) = (0.5*(1.-upwnd)*f(i+1)+upwnd*f(i))/del
          else
             i = llim
             df(i) = (0.5*(1.-upwnd)*f(i+1)+upwnd*f(i)-0.5*(1.+upwnd)*fl(2))/del
          end if
          if (iseg == nseg.and..not.periodic) then
             i = ulim
             ! if at rightmost zed, have no info for f(i+1) so must fully upwind
             df(i) = (f(i)-f(i-1))/del
          else
             i = ulim
             df(i) = (0.5*(1.-upwnd)*fr(1)+upwnd*f(i)-0.5*(1.+upwnd)*f(i-1))/del
          end if
          istart = llim
          iend = ulim
       end if

       ! mixed 2nd order centered and 1st order upwind scheme
       do i = istart-sgn, iend+sgn, -sgn
          df(i) = sgn*(0.5*(1.+upwnd)*f(i+sgn) - upwnd*f(i) + 0.5*(upwnd-1.)*f(i-sgn))/del
       end do

    end if

  end subroutine fd_variable_upwinding_zed

  subroutine fd_variable_upwinding_vpa (llim, f, del, sgn, upwnd, df)

    implicit none

    integer, intent (in) :: llim
    complex, dimension (llim:), intent (in) :: f
    real, intent (in) :: del, upwnd
    integer, intent (in) :: sgn
    complex, dimension (llim:), intent (out) :: df

    integer :: i, istart, iend, ulim

    ! if upwnd is zero or if z=0, then use centered differences
    if (abs(upwnd) < epsilon(0.) .or. sgn == 0) then
       call second_order_centered_vpa (llim, f, del, df)
    else
       ulim = size(f)+llim-1

       ! if sgn > 0, then stream speed is negative
       ! so sweep from more positive to more negative zed
       if (sgn > 0) then
          istart = ulim
          iend = llim
       else
          istart = llim
          iend = ulim
       end if

       ! zero_bc assumes that g -> zero beyond grid
       ! boundaries in vpa
       df(istart) = sgn*(0.5*(upwnd-1.0)*f(istart-sgn)-upwnd*f(istart))/del

       ! as do not have info beyond grid boundary at end of sweep
       ! use pure upwinding
       df(iend) = sgn*(f(iend+sgn)-f(iend))/del

       ! mixed centered and 1st order upwind scheme
       do i = istart-sgn, iend+sgn, -sgn
          df(i) = sgn*(0.5*(1.+upwnd)*f(i+sgn) - upwnd*f(i) + 0.5*(upwnd-1.)*f(i-sgn))/del
       end do

    end if

  end subroutine fd_variable_upwinding_vpa

  ! only good for equally-spaced grid-pts
  subroutine fd3pt_real (prof, profgrad, dr)
    
    implicit none
    
    real, dimension (:), intent (in) :: prof
    real, dimension (:), intent (out) :: profgrad
    real, intent (in) :: dr
    
    integer :: ix, npts
    real, dimension (:), allocatable :: aa, bb, cc
    
    npts = size(prof)
    allocate (aa(npts), bb(npts), cc(npts))
    
    aa = 1.0 ; bb = 4.0 ; cc = 1.0
    aa(1) = 0.0 ; bb(1) = 0.5 ; cc(1) = 0.5
    aa(npts) = 0.5 ; bb(npts) = 0.5 ; cc(npts) = 0.0
    
    do ix = 2, npts-1
       profgrad(ix) = 3.0 * (prof(ix+1) - prof(ix-1)) / dr
    end do
    profgrad(1) = (prof(2)-prof(1))/dr
    profgrad(npts) = (prof(npts)-prof(npts-1))/dr
    
    call tridag (aa, bb, cc, profgrad)
    
    deallocate (aa, bb, cc)
    
  end subroutine fd3pt_real
  
  subroutine fd3pt_real_array (prof, profgrad, dr)
    
    implicit none
    
    real, dimension (:), intent (in) :: prof, dr
    real, dimension (:), intent (out) :: profgrad
    
    integer :: ix, npts
    real :: a, b, c
    
    npts = size(prof)
    
    do ix = 2, npts-1
       profgrad(ix) = ((prof(ix)-prof(ix-1))*dr(ix)/dr(ix-1) &
            + (prof(ix+1)-prof(ix))*dr(ix-1)/dr(ix)) / (dr(ix-1)+dr(ix))
    end do
    ix = 1
    a = -(2.*dr(1) + dr(2))/(dr(1)*(dr(1)+dr(2)))
    b = (dr(1)+dr(2))/(dr(1)*dr(2))
    c = -dr(1)/(dr(2)*(dr(1)+dr(2)))
    profgrad(1) = a*prof(1)+b*prof(2)+c*prof(3)
    ix = npts
    a = dr(npts-1)/(dr(npts-2)*(dr(npts-2)+dr(npts-1)))
    b = -(dr(npts-1)+dr(npts-2))/(dr(npts-2)*dr(npts-1))
    c = (2.*dr(npts-1)+dr(npts-2))/(dr(npts-1)*(dr(npts-1)+dr(npts-2)))
    profgrad(npts) = a*prof(npts-2) + b*prof(npts-1) + c*prof(npts)
    
  end subroutine fd3pt_real_array

  subroutine fd3pt_complex_array (prof, profgrad, dr)
    
    implicit none
    
    complex, dimension (:), intent (in) :: prof
    real, dimension (:), intent (in) :: dr
    complex, dimension (:), intent (out) :: profgrad
    
    integer :: ix, npts
    real :: a, b, c
    
    npts = size(prof)
    
    do ix = 2, npts-1
       profgrad(ix) = ((prof(ix)-prof(ix-1))*dr(ix)/dr(ix-1) &
            + (prof(ix+1)-prof(ix))*dr(ix-1)/dr(ix)) / (dr(ix-1)+dr(ix))
    end do
    ix = 1
    a = -(2.*dr(1) + dr(2))/(dr(1)*(dr(1)+dr(2)))
    b = (dr(1)+dr(2))/(dr(1)*dr(2))
    c = -dr(1)/(dr(2)*(dr(1)+dr(2)))
    profgrad(1) = a*prof(1)+b*prof(2)+c*prof(3)
    ix = npts
    a = dr(npts-1)/(dr(npts-2)*(dr(npts-2)+dr(npts-1)))
    b = -(dr(npts-1)+dr(npts-2))/(dr(npts-2)*dr(npts-1))
    c = (2.*dr(npts-1)+dr(npts-2))/(dr(npts-1)*(dr(npts-1)+dr(npts-2)))
    profgrad(npts) = a*prof(npts-2) + b*prof(npts-1) + c*prof(npts)
    
  end subroutine fd3pt_complex_array

  ! boundary points are 2nd-order accurate (2-pt compact difference)
  ! next to boundary points are 4th-order accurate (2-pt centered compact difference)
  ! interior points are 6th-order accurate (4-pt centered compact difference)
  subroutine fd5pt_real (prof, profgrad, dr)

    implicit none

    real, dimension (:), intent (in) :: prof
    real, dimension (:), intent (out) :: profgrad
    real, intent (in) :: dr

    integer :: ix, npts
    real, dimension (:), allocatable :: aa, bb, cc

    npts = size(prof)
    allocate (aa(npts), bb(npts), cc(npts))

    aa = 1.0 ; bb = 3.0 ; cc = 1.0
    aa(1) = 0.0 ; bb(1) = 0.5 ; cc(1) = 0.5
    aa(2) = 1.0 ; bb(2) = 4.0 ; cc(2) = 1.0
    aa(npts-1) = 1.0 ; bb(npts-1) = 4.0 ; cc(npts-1) = 1.0
    aa(npts) = 0.5 ; bb(npts) = 0.5 ; cc(npts) = 0.0

    do ix = 3, npts-2
       profgrad(ix) = (7.*(prof(ix+1) - prof(ix-1)) + 0.25*(prof(ix+2)-prof(ix-2))) / (3.*dr)
    end do
    profgrad(1) = (prof(2)-prof(1))/dr
    profgrad(2) = 3.0*(prof(3) - prof(1))/dr
    profgrad(npts-1) = 3.0*(prof(npts) - prof(npts-2))/dr
    profgrad(npts) = (prof(npts)-prof(npts-1))/dr

    call tridag (aa, bb, cc, profgrad)

    deallocate (aa, bb, cc)

  end subroutine fd5pt_real

  ! boundary points are 2nd-order accurate (2-pt compact difference)
  ! next to boundary points are 4th-order accurate (2-pt centered compact difference)
  ! interior points are 6th-order accurate (4-pt centered compact difference)
  subroutine fd5pt_array (prof, profgrad, dr)

    implicit none

    real, dimension (:), intent (in) :: prof, dr
    real, dimension (:), intent (out) :: profgrad

    integer :: ix, npts
    real, dimension (:), allocatable :: aa, bb, cc

    npts = size(prof)
    allocate (aa(npts), bb(npts), cc(npts))

    aa = 1.0 ; bb = 3.0 ; cc = 1.0
    aa(1) = 0.0 ; bb(1) = 0.5 ; cc(1) = 0.5
    aa(2) = 1.0 ; bb(2) = 4.0 ; cc(2) = 1.0
    aa(npts-1) = 1.0 ; bb(npts-1) = 4.0 ; cc(npts-1) = 1.0
    aa(npts) = 0.5 ; bb(npts) = 0.5 ; cc(npts) = 0.0

    do ix = 3, npts-2
       profgrad(ix) = (7.*(prof(ix+1) - prof(ix-1)) + 0.25*(prof(ix+2)-prof(ix-2))) / (3.*dr(ix))
    end do
    profgrad(1) = (prof(2)-prof(1))/dr(1)
    profgrad(2) = 3.0*(prof(3) - prof(1))/dr(2)
    profgrad(npts-1) = 3.0*(prof(npts) - prof(npts-2))/dr(npts-1)
    profgrad(npts) = (prof(npts)-prof(npts-1))/dr(npts)

    call tridag (aa, bb, cc, profgrad)

    deallocate (aa, bb, cc)

  end subroutine fd5pt_array


  ! second derivative using centered differences
  ! second order accurate
  subroutine d2_3pt_real (f, d2f, dr)

    implicit none

    real, dimension (:), intent (in) :: f
    real, dimension (:), intent (in) :: dr
    real, dimension (:), intent (out) :: d2f

    real :: a, b, c, d
    integer :: i, n

    n = size(f)

    do i = 2, n-1
       a = 2./(dr(i-1)*(dr(i)+dr(i-1)))
       b = -2./(dr(i-1)*dr(i))
       c = 2./(dr(i)*(dr(i)+dr(i-1)))
       d2f(i) = a*f(i-1)+b*f(i)+c*f(i+1)
    end do
    i = 1
    a = (6.*dr(1)+4.*dr(2)+2.*dr(3)) &
         / (dr(1)*(dr(1)+dr(2))*(dr(1)+dr(2)+dr(3)))
    b = -(4.*(dr(1)+dr(2))+2.*dr(3)) &
         / (dr(1)*dr(2)*(dr(2)+dr(3)))
    c = (4.*dr(1)+2.*(dr(2)+dr(3))) &
         / ((dr(2)+dr(1))*dr(2)*dr(3))
    d = -(4.*dr(1)+2.*dr(2)) &
         / ((dr(3)+dr(2)+dr(1))*(dr(3)+dr(2))*dr(3))
    d2f(i) = a*f(1)+b*f(2)+c*f(3)+d*f(4)
    i = n
    a = -(4.*dr(n-1)+2.*dr(n-2)) &
         / (dr(n-3)*(dr(n-3)+dr(n-2))*(dr(n-3)+dr(n-2)+dr(n-1)))
    b = (4.*dr(n-1)+2.*(dr(n-2)+dr(n-3))) &
         / (dr(n-3)*dr(n-2)*(dr(n-2)+dr(n-1)))
    c = -(4.*(dr(n-1)+dr(n-2))+2.*dr(n-3)) &
         / (dr(n-1)*dr(n-2)*(dr(n-2)+dr(n-3)))
    d = (6.*dr(n-1)+4.*dr(n-2)+2.*dr(n-3)) &
         / (dr(n-1)*(dr(n-1)+dr(n-2))*(dr(n-1)+dr(n-2)+dr(n-3)))
    d2f(i) = a*f(n-3)+b*f(n-2)+c*f(n-1)+d*f(n)

!     ! FLAG -- this is a hack
!     ! do not anticipate needing 2nd derivatives
!     ! at first and last grid points
!     d2f(1) = d2f(2)
!     d2f(n) = d2f(n-1)

  end subroutine d2_3pt_real

  subroutine d2_3pt_complex (f, d2f, dr)

    implicit none

    complex, dimension (:), intent (in) :: f
    real, dimension (:), intent (in) :: dr
    complex, dimension (:), intent (out) :: d2f

    real :: a, b, c, d
    integer :: i, n

    n = size(f)

    do i = 2, n-1
       a = 2./(dr(i-1)*(dr(i)+dr(i-1)))
       b = -2./(dr(i-1)*dr(i))
       c = 2./(dr(i)*(dr(i)+dr(i-1)))
       d2f(i) = a*f(i-1)+b*f(i)+c*f(i+1)
    end do
    i = 1
    a = (6.*dr(1)+4.*dr(2)+2.*dr(3)) &
         / (dr(1)*(dr(1)+dr(2))*(dr(1)+dr(2)+dr(3)))
    b = -(4.*(dr(1)+dr(2))+2.*dr(3)) &
         / (dr(1)*dr(2)*(dr(2)+dr(3)))
    c = (4.*dr(1)+2.*(dr(2)+dr(3))) &
         / ((dr(2)+dr(1))*dr(2)*dr(3))
    d = -(4.*dr(1)+2.*dr(2)) &
         / ((dr(3)+dr(2)+dr(1))*(dr(3)+dr(2))*dr(3))
    d2f(i) = a*f(1)+b*f(2)+c*f(3)+d*f(4)
    i = n
    a = -(4.*dr(n-1)+2.*dr(n-2)) &
         / (dr(n-3)*(dr(n-3)+dr(n-2))*(dr(n-3)+dr(n-2)+dr(n-1)))
    b = (4.*dr(n-1)+2.*(dr(n-2)+dr(n-3))) &
         / (dr(n-3)*dr(n-2)*(dr(n-2)+dr(n-1)))
    c = -(4.*(dr(n-1)+dr(n-2))+2.*dr(n-3)) &
         / (dr(n-1)*dr(n-2)*(dr(n-2)+dr(n-3)))
    d = (6.*dr(n-1)+4.*dr(n-2)+2.*dr(n-3)) &
         / (dr(n-1)*(dr(n-1)+dr(n-2))*(dr(n-1)+dr(n-2)+dr(n-3)))
    d2f(i) = a*f(n-3)+b*f(n-2)+c*f(n-1)+d*f(n)

  end subroutine d2_3pt_complex

  subroutine tridag_real (aa, bb, cc, sol)
    
    implicit none
    
    real, dimension (:), intent (in) :: aa, bb, cc
    real, dimension (:), intent (in out) :: sol
    
    integer :: ix, npts
    real :: bet
    
    real, dimension (:), allocatable :: gam
    
    npts = size(aa)
    allocate (gam(npts))
    
    bet = bb(1)
    sol(1) = sol(1)/bet
    
    do ix = 2, npts
       gam(ix) = cc(ix-1)/bet
       bet = bb(ix) - aa(ix)*gam(ix)
       if (bet == 0.0) write (*,*) 'tridiagonal solve failed'
       sol(ix) = (sol(ix)-aa(ix)*sol(ix-1))/bet
    end do

    do ix = npts-1, 1, -1
       sol(ix) = sol(ix) - gam(ix+1)*sol(ix+1)
    end do

    deallocate (gam)

  end subroutine tridag_real

  subroutine tridag_complex (llim, aa, bb, cc, sol)
    
    implicit none
    
    integer, intent (in) :: llim
    real, dimension (llim:), intent (in) :: aa, bb, cc
    complex, dimension (llim:), intent (in out) :: sol
    
    integer :: ix, npts
    real :: bet
    
    real, dimension (:), allocatable :: gam
    
    npts = size(bb)
    allocate (gam(llim:llim+npts-1))
    
    bet = bb(llim)
    sol(llim) = sol(llim)/bet
    
    do ix = llim+1, llim+npts-1
       gam(ix) = cc(ix-1)/bet
       bet = bb(ix) - aa(ix)*gam(ix)
       if (bet == 0.0) write (*,*) 'tridiagonal solve failed'
       sol(ix) = (sol(ix)-aa(ix)*sol(ix-1))/bet
    end do

    do ix = llim+npts-2, llim, -1
       sol(ix) = sol(ix) - gam(ix+1)*sol(ix+1)
    end do

    deallocate (gam)

  end subroutine tridag_complex
  
end module finite_differences
