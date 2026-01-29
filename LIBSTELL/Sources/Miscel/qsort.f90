!------------------------------------------------------------------------------
!     Module:        qsort
!     Author:        Lucas van Ham (lucas.van.ham@ipp.mpg.de)
!     Date:          January 2026
!     Description:   Quicksorts an array of indices based on an array of values
!                    (i.e. keeps track of permutations)
!                    For simple quicksort (sort values based on values), use
!                    DLASRT (LAPACK)
!------------------------------------------------------------------------------
    MODULE qsort
!------------------------------------------------------------------------------
!     Functions
!        partition: partitions input array around pivot point
!        quicksort: the actual quicksorting function
!------------------------------------------------------------------------------
    CONTAINS

    FUNCTION partition(array, indices, dex_left, dex_right) result(dex_pivot)
    !-----------------------------------------------------------------------
    ! param[in]:    array:     array to sort
    ! param[inout]: indices:   indices of incoming array
    ! param[in]:    dex_left:  start index of section
    ! param[in]:    dex_right: end index of section
    ! param[out]:   dex_pivot: index of pivot element
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(in) :: array(:)
    INTEGER, INTENT(inout) :: indices(:)
    INTEGER, INTENT(in)  :: dex_left, dex_right
    INTEGER :: dex_pivot
    INTEGER :: i, temp
    DOUBLE PRECISION :: pivot
    
    pivot = array(indices(dex_right))
    dex_pivot = dex_left - 1

    DO i = dex_left, dex_right-1
      IF (array(indices(i)).LE.pivot) THEN
        dex_pivot = dex_pivot + 1
        temp = indices(i)
        indices(i) = indices(dex_pivot)
        indices(dex_pivot) = temp
      END IF
    END DO

    temp = indices(dex_pivot+1)
    indices(dex_pivot+1) = indices(dex_right)
    indices(dex_right) = temp

    dex_pivot = dex_pivot+1    

    END FUNCTION partition

    RECURSIVE SUBROUTINE quicksort(array, indices, dex_left, dex_right)
    !-----------------------------------------------------------------------
    ! param[in]:    array:     array to sort
    ! param[inout]: indices:   indices of incoming array
    ! param[in]:    dex_left:  start index of section
    ! param[in]:    dex_right: end index of section
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    
    DOUBLE PRECISION, INTENT(in) :: array(:)
    INTEGER, INTENT(inout) :: indices(:)
    INTEGER, INTENT(in)  :: dex_left, dex_right
    INTEGER :: dex_pivot

    IF (dex_left.GE.dex_right) RETURN
    dex_pivot = partition(array, indices, dex_left, dex_right)

    quicksort(array, indices, dex_left, dex_pivot-1)
    quicksort(array, indices, dex_pivot+1,dex_right)

    END SUBROUTINE quicksort

    END MODULE qsort 