subroutine sort(n, a, idx)
    implicit none
    integer :: n, i, j, d
    integer :: idx(n)
    real :: a(n), x
 
    do i = 1, n
       idx(i) = i
    end do

    do i = 2, n
        x = a(i)
        d = idx(i)
        j = i - 1
        do while (j >= 1)
            if (a(j) <= x) exit
            idx(j+1) = idx(j)
            a(j + 1) = a(j)
            j = j - 1
        end do
        idx(j+1) = d
        a(j + 1) = x
    end do
    return
end subroutine