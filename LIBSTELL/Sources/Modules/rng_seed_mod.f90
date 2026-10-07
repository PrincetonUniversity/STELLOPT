!-----------------------------------------------------------------------
!     Module:        rng_seed_mod
!     Authors:       D. Kulla (david.kulla@ipp.mpg.de)
!     Date:          08/06/2026
!     Description:   This module seeds the intrinsic RANDOM_NUMBER
!                    generator.  It replaces the former stand alone
!                    init_random_seed subroutine and adds an optional
!                    user supplied seed so that runs which draw random
!                    numbers (collision operators, beam deposition,
!                    field line diffusion) can be reproduced exactly.
!
!                    Every MPI rank must draw an independent sequence,
!                    otherwise the collisional noise is correlated
!                    across ranks.  The user seed is therefore mixed
!                    with a stream identifier (the MPI rank) before it
!                    is handed to RANDOM_SEED.  Results are reproducible
!                    for a given number of ranks, not across different
!                    rank counts.
!-----------------------------------------------------------------------
      MODULE rng_seed_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE iso_fortran_env, ONLY: int64

!-----------------------------------------------------------------------
!     Module Variables
!          RNG_SEED_RANDOM   Value of a seed which requests a
!                            non-reproducible (entropy based) stream.
!-----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, PARAMETER :: RNG_SEED_RANDOM = -1

      PRIVATE :: lcg

!-----------------------------------------------------------------------
!     Subroutines
!         init_rng_seed:   Seeds RANDOM_NUMBER
!-----------------------------------------------------------------------
      CONTAINS

!-----------------------------------------------------------------------
!     Subroutine:    init_rng_seed
!     Description:   Seeds the intrinsic random number generator.  If
!                    user_seed is negative (RNG_SEED_RANDOM) the seed is
!                    taken from /dev/urandom, falling back to the system
!                    clock XORed with the PID.  Otherwise a
!                    deterministic seed is built from user_seed and
!                    stream_id.
!     Arguments:
!         user_seed  Seed from the input namelist, <0 for random
!         stream_id  Stream identifier, pass the MPI rank (default 0)
!-----------------------------------------------------------------------
      SUBROUTINE init_rng_seed(user_seed, stream_id)
      IMPLICIT NONE
      INTEGER, INTENT(in)           :: user_seed
      INTEGER, INTENT(in), OPTIONAL :: stream_id
      INTEGER :: i, n, un, istat, mystream, dt(8), pid
      INTEGER(int64) :: t
      INTEGER, ALLOCATABLE :: seed(:)
      INTEGER*4 :: getpid ! Need for ifort

      mystream = 0
      IF (PRESENT(stream_id)) mystream = stream_id

      CALL RANDOM_SEED(size = n)
      ALLOCATE(seed(n))

      IF (user_seed >= 0) THEN
         ! Deterministic.  Mix the user seed with the stream id so that
         ! neighbouring ranks do not start from neighbouring states.
         t = INT(user_seed, int64) + INT(mystream, int64) * 2654435761_int64
         DO i = 1, n
            seed(i) = lcg(t)
         END DO
      ELSE
         ! First try if the OS provides a random number generator
         OPEN(newunit=un, file="/dev/urandom", access="stream", &
              form="unformatted", action="read", status="old", iostat=istat)
         IF (istat == 0) THEN
            READ(un) seed
            CLOSE(un)
         ELSE
            ! Fallback to XOR:ing the current time and pid. The PID is
            ! useful in case one launches multiple instances of the same
            ! program in parallel.
            CALL SYSTEM_CLOCK(t)
            IF (t == 0) THEN
               CALL DATE_AND_TIME(values=dt)
               t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
                    + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
                    + dt(3) * 24_int64 * 60 * 60 * 1000 &
                    + dt(5) * 60 * 60 * 1000 &
                    + dt(6) * 60 * 1000 + dt(7) * 1000 &
                    + dt(8)
            END IF
            pid = getpid()
            t = IEOR(t, INT(pid, KIND(t)))
            t = IEOR(t, INT(mystream, KIND(t)) * 2654435761_int64)
            DO i = 1, n
               seed(i) = lcg(t)
            END DO
         END IF
      END IF

      CALL RANDOM_SEED(put=seed)
      DEALLOCATE(seed)

      END SUBROUTINE init_rng_seed

!-----------------------------------------------------------------------
!     Function:      lcg
!     Description:   Simple linear congruential generator.  Not good
!                    enough for real work, but sufficient for seeding a
!                    better PRNG.  All intermediates stay below 2**61 so
!                    no integer overflow can occur.
!-----------------------------------------------------------------------
      INTEGER FUNCTION lcg(s)
      IMPLICIT NONE
      INTEGER(int64), INTENT(inout) :: s
      IF (s == 0) THEN
         s = 104729
      ELSE
         s = MOD(s, 4294967296_int64)
      END IF
      s = MOD(s * 279470273_int64, 4294967291_int64)
      lcg = INT(MOD(s, INT(HUGE(0), int64)), KIND(0))
      END FUNCTION lcg

      END MODULE rng_seed_mod
