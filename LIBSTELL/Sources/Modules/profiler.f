!*******************************************************************************
!>  @file profiler.f
!>  @brief Contains module @ref profiler.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines functions for measuring an tabulating performance of function and
!>  subroutine calls. These routines are only active when the PROFILE_ON macro
!>  is defined.
!*******************************************************************************

      MODULE profiler
      USE stel_kinds

      IMPLICIT NONE

!*******************************************************************************
!  profiler module parameters
!*******************************************************************************
!>  Max string length
      INTEGER, PARAMETER :: profiler_string_size = 68
!>  Max number of buckets
      INTEGER, PARAMETER :: profiler_bucket_size = 1000

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) bucket
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Full table of profiled functions.
!-------------------------------------------------------------------------------
      TYPE profiler_bucket
!>  OpenMP lock variable.
!$       INTEGER (kind=8)                     :: lock
!>  Number of calls to a function.
         INTEGER                              :: number_of_calls = 0
!>  Total time spent in that function.
         REAL (rprec)                         :: total_time = 0.0
!>  Total time spent in that function.
         REAL (rprec)                         :: average_time = 0.0
!>  Symbol name of that function.
         CHARACTER (len=profiler_string_size) :: symbol_name = ''
      END TYPE

!*******************************************************************************
!  profiler module variables
!*******************************************************************************
#if PROFILE_ON
!>  Array of buckets to hold the values.
      TYPE (profiler_bucket), DIMENSION(profiler_bucket_size), SAVE            &
     &   :: buckets
#endif

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref profiler_table object.
!>
!>  Allocates memory for a profiler table.
!>
!>  @returns A pointer to a constructed @ref profiler_table object.
!-------------------------------------------------------------------------------
      SUBROUTINE profiler_construct()

      IMPLICIT NONE

!  local variables
!$    INTEGER :: i

!  Start of executable code
#if PROFILE_ON
!$    DO i = 1, profiler_bucket_size
!$       CALL OMP_INIT_LOCK(buckets(i)%lock)
!$    END DO
#endif
      END SUBROUTINE

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref profiler_table object.
!>
!>  Deallocates memory for a @ref profiler_table object.
!>
!>  @param[inout] this A @ref profiler_table instance.
!-------------------------------------------------------------------------------
      SUBROUTINE profiler_destruct()

      IMPLICIT NONE

!  local variables
!$    INTEGER :: i

!  Start of executable code
#if PROFILE_ON
!$    DO i = 1, profiler_bucket_size
!$       CALL OMP_DESTROY_LOCK(buckets(i)%lock)
!$    END DO
#endif
      END SUBROUTINE

!*******************************************************************************
!  SETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Gets the end time of profiled function.
!>
!>  Sets the end time of profiled function and adds the result to the table. The
!>  symbol name is hashed to determine the index of the array it lives in. If
!>  the index is full, the next instance is used. When profiling is turned off,
!>  this does nothing. To avoid race conditons in multi threaded code, this must
!>  only be run by one thread or process at a time.
!>
!>  @param[in]    symbol_name Name of the symbol to profile.
!>  @param[in]    start_time  Starting time of the symbol call.
!-------------------------------------------------------------------------------
      SUBROUTINE profiler_set_stop_time(symbol_name, start_time)

!  Declare Arguments
      CHARACTER (len=*), INTENT(in)        :: symbol_name
      REAL (rprec), INTENT(in)             :: start_time

!  local variables
      INTEGER                              :: i, index
      REAL (rprec)                         :: time

!  Start of executable code.
#if PROFILE_ON
      CALL second0(time)
      time = time - start_time
      index = profiler_hash_function(symbol_name)

      DO i = index, profiler_bucket_size
!$       CALL OMP_SET_LOCK(buckets(i)%lock)
         IF (buckets(i)%symbol_name .eq. '') THEN
            buckets(i)%symbol_name = symbol_name
            buckets(i)%number_of_calls = 1
            buckets(i)%total_time = time
            buckets(i)%average_time = time
!$          CALL OMP_UNSET_LOCK(buckets(i)%lock)
            RETURN
         ELSE IF (buckets(i)%symbol_name .eq. symbol_name) THEN
            buckets(i)%number_of_calls = buckets(i)%number_of_calls + 1
            buckets(i)%total_time = buckets(i)%total_time + time
            buckets(i)%average_time = buckets(i)%total_time                    &
     &                              / buckets(i)%number_of_calls
!$          CALL OMP_UNSET_LOCK(buckets(i)%lock)
            RETURN
         END IF
!$       CALL OMP_UNSET_LOCK(buckets(i)%lock)
      END DO

      DO i = 1, index - 1
!$       CALL OMP_SET_LOCK(buckets(i)%lock)
         IF (buckets(i)%symbol_name .eq. '') THEN
            buckets(i)%symbol_name = symbol_name
            buckets(i)%number_of_calls = 1
            buckets(i)%total_time = time
            buckets(i)%average_time = time
!$          CALL OMP_UNSET_LOCK(buckets(i)%lock)
            RETURN
         ELSE IF (buckets(i)%symbol_name .eq. symbol_name) THEN
            buckets(i)%number_of_calls = buckets(i)%number_of_calls + 1
            buckets(i)%total_time = buckets(i)%total_time + time
            buckets(i)%average_time = buckets(i)%total_time                    &
     &                              / buckets(i)%number_of_calls
!$          CALL OMP_UNSET_LOCK(buckets(i)%lock)
            RETURN
         END IF
!$       CALL OMP_UNSET_LOCK(buckets(i)%lock)
      END DO

      WRITE (*,*) 'Profile table full: ' // TRIM(symbol_name) //               &
     &            ' could not be profiled.'
#endif

       END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Gets the start time of profiled function.
!>
!>  Gets a starting time of a function. When profiling is turned off, this
!>  returns 0 instead.
!>
!>  @returns The start time of a function.
!-------------------------------------------------------------------------------
      FUNCTION profiler_get_start_time()

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: profiler_get_start_time

!  Start of executable code
#if PROFILE_ON
      CALL second0(profiler_get_start_time)
#else
      profiler_get_start_time = 0.0
#endif

      END FUNCTION

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Sorts the profile table based on the average call time.
!>
!>  A recursive merge sort algorithm.
!>
!>  @param[in] low_index  Lower index to sort.
!>  @param[in] high_index Higher index to sort.
!-------------------------------------------------------------------------------
      RECURSIVE SUBROUTINE profiler_sort(low_index, high_index)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER, INTENT(in)                               :: low_index
      INTEGER, INTENT(in)                               :: high_index

!  local variables
      INTEGER                                           :: mid_index
      INTEGER                                           :: i, i1, i2
      TYPE (profiler_bucket)                            :: swap_bucket
      TYPE (profiler_bucket), DIMENSION(:), ALLOCATABLE :: temp_buckets

!  Start of executable code
#if PROFILE_ON
      IF (high_index - low_index .lt. 2) THEN
!  Since the table has been reduced to two elements, sort the sub table.
         IF (buckets(high_index)%average_time .gt.                             &
     &       buckets(low_index)%average_time) THEN
!  Swap the values
            swap_bucket = buckets(low_index)
            buckets(low_index) = buckets(high_index)
            buckets(high_index) = swap_bucket
         END IF
         RETURN
      END IF

!  Split the table in the middle and sort the sub tables.
      mid_index = (high_index + low_index)/2
      CALL profiler_sort(low_index, mid_index - 1)
      CALL profiler_sort(mid_index, high_index)

!  Merge the two tables.
      i1 = low_index
      i2 = mid_index

      ALLOCATE(temp_buckets(low_index:high_index))

      DO i = low_index, high_index
         IF (buckets(i1)%average_time .gt.                                     &
     &       buckets(i2)%average_time) THEN
!  Take the value of the lower table and increment the lower index
            temp_buckets(i) = buckets(i1)
            i1 = i1 + 1
            IF (i1 .gt. mid_index - 1) THEN
!  Merged the last value in the lower table the remining values are all from the
!  upper.
               temp_buckets(i + 1:high_index) = buckets(i2:high_index)
               EXIT
            END IF
         ELSE
!  Take the value of the upper table and increment the upper index
            temp_buckets(i) = buckets(i2)
            i2 = i2 + 1
            IF (i2 .gt. high_index) THEN
!  Merged the last value in the upper tabel the remining values are all from the
!  lower.
               temp_buckets(i + 1:high_index) =                                &
     &            buckets(i1:mid_index - 1)
               EXIT
            END IF
         END IF
      END DO

      buckets(low_index:high_index) = temp_buckets(low_index:high_index)

      DEALLOCATE(temp_buckets)
#endif
      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Computes a hash for the symbol name.
!>
!>  The hash is computed using the djb2 algorithm.
!>
!>      hash_0 = 5381
!>      hash_i+1 = hash_i * 33 ^ c_i
!>
!>  The magic numbers of 5381 and 33 help the hash function avoid collision and
!>  good distribution.
!>
!>  @param[in] symbol_name Name of the symbol to hash.
!>  @returns The hash value of the symbol.
!-------------------------------------------------------------------------------
      FUNCTION profiler_hash_function(symbol_name)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                       :: profiler_hash_function
      CHARACTER (len=*), INTENT(in) :: symbol_name

!  local variables
      INTEGER (kind=8)              :: hash
      INTEGER                       :: i

!  Start of executable code
      hash = 5381
      DO i = 1, LEN_TRIM(symbol_name)
         hash = (LSHIFT(hash, 5) + hash) + ICHAR(symbol_name(i:i))
      END DO

      profiler_hash_function = MOD(ABS(hash), profiler_bucket_size) + 1

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write out the profiled data to an output file.
!>
!>  This loops through the buckets and prints out the profiled data. Any empty
!>  buckets get skiped.
!>
!>  @param[in] iou  Input/output unit of the output file.
!-------------------------------------------------------------------------------
      SUBROUTINE profiler_write(iou)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER, INTENT(in)               :: iou

!  local variables
      INTEGER                           :: i

!  Start of executable code
#if PROFILE_ON
      WRITE (iou,*)
      WRITE (iou,*) ' *** Code Profile Information'

      CALL profiler_sort(1, profiler_bucket_size)

      WRITE (iou,1000)

      DO i = 1, profiler_bucket_size
         IF (buckets(i)%symbol_name .ne. '') THEN
            WRITE (iou, 1001) buckets(i)%symbol_name,                          &
     &                        buckets(i)%average_time,                         &
     &                        buckets(i)%total_time,                           &
     &                        buckets(i)%number_of_calls
         END IF
      END DO

1000  FORMAT('Symbol Name',59x,'Average Time',2x,'TotalTime',5x,               &
     &       'Num Calls')
1001  FORMAT(a68,2x,es12.5,2x,es12.5,2x,i9)
#endif

      END SUBROUTINE

      END MODULE
