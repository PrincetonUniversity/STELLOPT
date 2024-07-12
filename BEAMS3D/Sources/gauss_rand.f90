      SUBROUTINE gauss_rand(n,x)
!  Subroutine to fill an array with gaussian random variables of
!  unit variance. Pass size of array, so that don't have to 
!  fuss with interfaces

      USE stel_kinds, ONLY: rprec

      INTEGER, INTENT(in) :: n
      REAL(rprec), DIMENSION(n), INTENT(inout) :: x

      INTEGER :: np1o2, i, j
      REAL(rprec), DIMENSION(2) :: v
      REAL(rprec) :: rsq, fac

!  start of executable code
      np1o2 = (n + 1) / 2
      j = 0
      DO i = 1,np1o2
         j = j + 1
100      CALL RANDOM_NUMBER(v)
         v = 2. * v - 1.
         rsq = v(1) * v(1) + v(2) * v(2)
         IF ((rsq .ge. 1.) .or. (rsq .eq. 0.)) GO TO 100
         fac = SQRT(-2. * LOG(rsq) / rsq)
         x(j) = v(1) * fac
         j = j + 1
         IF (j .le. n) THEN
            x(j) = v(2) * fac
         ENDIF
      END DO
      RETURN
      END SUBROUTINE gauss_rand

          subroutine init_random_seed()
            use iso_fortran_env, only: int64
            implicit none
            integer, allocatable :: seed(:)
            integer :: i, n, un, istat, dt(8), pid
            integer(int64) :: t
            INTEGER*4 :: getpid ! Need for ifort
          
            call random_seed(size = n)
            allocate(seed(n))
            ! First try if the OS provides a random number generator
            open(newunit=un, file="/dev/urandom", access="stream", &
                 form="unformatted", action="read", status="old", iostat=istat)
            if (istat == 0) then
               read(un) seed
               close(un)
            else
               ! Fallback to XOR:ing the current time and pid. The PID is
               ! useful in case one launches multiple instances of the same
               ! program in parallel.
               call system_clock(t)
               if (t == 0) then
                  call date_and_time(values=dt)
                  t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
                       + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
                       + dt(3) * 24_int64 * 60 * 60 * 1000 &
                       + dt(5) * 60 * 60 * 1000 &
                       + dt(6) * 60 * 1000 + dt(7) * 1000 &
                       + dt(8)
               end if
               pid = getpid()
               t = ieor(t, int(pid, kind(t)))
               do i = 1, n
                  seed(i) = lcg(t)
               end do
            end if
            call random_seed(put=seed)
            DEALLOCATE(seed)
          contains
            ! This simple PRNG might not be good enough for real work, but is
            ! sufficient for seeding a better PRNG.
            function lcg(s)
              integer :: lcg
              integer(int64) :: s
              if (s == 0) then
                 s = 104729
              else
                 s = mod(s, 4294967296_int64)
              end if
              s = mod(s * 279470273_int64, 4294967291_int64)
              lcg = int(mod(s, int(huge(0), int64)), kind(0))
            end function lcg
          end subroutine init_random_seed
