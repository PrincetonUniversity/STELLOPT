      SUBROUTINE init_random_seed
      INTEGER :: i, n, clock
      INTEGER, ALLOCATABLE :: seed(:)
     
      CALL random_seed(size=n)
      ALLOCATE(seed(n))
      
      CALL system_clock(count=clock)
      seed = clock + 37 * (/ (i-1,i=1,n)/)
      CALL random_seed(put=seed)
      DEALLOCATE(seed)
      END SUBROUTINE init_random_seed
