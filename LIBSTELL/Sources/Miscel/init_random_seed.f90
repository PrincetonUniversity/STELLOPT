SUBROUTINE init_random_seed(nseed)
  INTEGER, OPTIONAL, INTENT(IN) :: nseed
  INTEGER :: i, n, clock
  INTEGER, ALLOCATABLE :: seed(:)
    CALL random_seed(size=n)
    ALLOCATE(seed(n))
    CALL system_clock(count=clock)
IF (PRESENT(nseed) .and. nseed .ne. -1) THEN
    seed = nseed
ELSE
    seed = clock + 37 * (/ (i-1,i=1,n)/)
END IF
CALL random_seed(put=seed)
DEALLOCATE(seed)
END SUBROUTINE init_random_seed
