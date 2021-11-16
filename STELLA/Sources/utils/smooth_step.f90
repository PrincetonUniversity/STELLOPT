module smooth_step

  implicit none

  public :: smoothstep

contains

  pure function smoothstep (x,N,minV,maxV)

    implicit none

    real :: smoothstep
    real, intent (in) :: x
    integer, intent (in) :: N 
    real, optional, intent (in) :: minV
    real, optional, intent (in) :: maxV
    real :: minVl 
    real :: maxVl
    real :: dV

    minVl = 0
    maxVl = 1

    if(present(minV)) minVl = minV
    if(present(maxV)) maxVl = maxV

    dV = maxVl - minVl

    if(x .le. 0) then
      smoothstep = minVl
      return
    endif
    if(x .ge. 1) then
      smoothstep = maxVl
      return
    endif

    select case (N)
      case(0)
        smoothstep = minVl + dV*smoothstep0(x)
      case(1)
        smoothstep = minVl + dV*smoothstep1(x)
      case(2)
        smoothstep = minVl + dV*smoothstep2(x)
      case default
        smoothstep = minVl + dV*smoothstepN(x,N)
    end select  
  end function smoothstep

  pure function smoothstep0 (x)
    implicit none
    real :: smoothstep0
    real, intent (in) :: x
    smoothstep0 = x
  end function smoothstep0

  pure function smoothstep1 (x)
    implicit none
    real :: smoothstep1
    real, intent (in) :: x
    smoothstep1 = x*x*(3-2*x)
  end function smoothstep1

  pure function smoothstep2 (x)
    implicit none
    real :: smoothstep2
    real, intent (in) :: x
    smoothstep2 = x*x*x*(x*(x*6-15)+10)
  end function smoothstep2

  pure function smoothstepN (x,N)
    implicit none
    real  :: smoothstepN
    real :: sumV
    real :: xp
    real, intent (in) :: x
    integer, intent (in) :: N
    integer :: i

    if(N < 0) then
      smoothstepN = 0.5*sin(3.14159265358979*(x-0.5))+0.5
      return
    endif

    xp = x
    if(x > 0.5) xp = 1.0-x
    sumV = 0
    do i = N, 0, -1
      sumV = xp*sumV + &
             pascalTriangle(-N-1,i) * &
             pascalTriangle(2*N + 1, N-i)
    end do
    sumV = xp**(N+1)*sumV
    smoothstepN = sumV
    if(x > 0.5) smoothstepN = 1-sumV
  end function smoothstepN

  pure function pascalTriangle (a,b)
    implicit none
    integer :: pascalTriangle
    integer, intent (in) :: a,b
    integer i
    pascalTriangle = 1
    do i = 0, b-1
      pascalTriangle = pascalTriangle*(a-i)/(i+1)
    end do
  end function pascalTriangle

end module smooth_step
