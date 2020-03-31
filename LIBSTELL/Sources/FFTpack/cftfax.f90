      SUBROUTINE cftfax_g(n, ifax, trigs)
      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER n
      INTEGER, DIMENSION(13) :: ifax
      REAL(rprec), DIMENSION(*) :: trigs
#if !defined(CRAY) || defined(LONESTAR) || defined(MCURIE)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: k
!-----------------------------------------------
!
! on input     n
!               the length of each complex transform to be performed
!
!               n must be greater than 1 and contain no prime
!               factors greater than 5.
!
! on output    ifax
!               ifax(1)
!                 the number of factors chosen or -99 in CASE of error
!               ifax(2) thru ifax( ifax(1)+1 )
!                 the factors of n in the followin order:  appearing
!                 first are as many factors of 4 as can be obtained.
!                 subsequent factors are primes, and appear in
!                 ascending order, except for multiple factors.
!
!              trigs
!               2n SIN and COS values for USE by the transform routine
!
      CALL fact_g (n, ifax)
      k = ifax(1)
      IF (k<1 .or. ifax(k+1)>5) ifax(1) = -99
      IF (ifax(1) .le. 0) STOP 'IFAX(1) <= 0 in CFTFAX'
      CALL cftrig_g (n, trigs)
#else
      CALL cftfax (n, ifax, trigs)
#endif
      END SUBROUTINE cftfax_g
