module sort

  implicit none

  public :: sort_array_ascending
  public :: unsort_array_ascending

contains

  subroutine sort_array_ascending (array, sort_index)

    implicit none

    real, dimension (:), intent (in out) :: array
    integer, dimension (:), intent (out) :: sort_index

    integer :: i, j, n
    real, dimension (:), allocatable :: tmp

    n = size(array)

    sort_index = 1
    do i = 1, n-1
       do j = i+1, n
          if (array(i) > array(j)) then
             sort_index(i) = sort_index(i) + 1
          else
             sort_index(j) = sort_index(j) + 1
          end if
       end do
    end do

    allocate (tmp(n))
    tmp = array
    do i = 1, n
       array(sort_index(i)) = tmp(i)
    end do
    deallocate(tmp)

  end subroutine sort_array_ascending

  subroutine unsort_array_ascending (array, sort_index)

    implicit none

    real, dimension (:), intent (in out) :: array
    integer, dimension (:), intent (in) :: sort_index

    integer :: i, n
    real, dimension (:), allocatable :: tmp

    n = size(array)
    
    allocate (tmp(n))
    tmp = array
    do i = 1, n
       array(i) = tmp(sort_index(i))
    end do
    deallocate (tmp)

  end subroutine unsort_array_ascending

end module sort
