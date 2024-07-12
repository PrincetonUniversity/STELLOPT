      SUBROUTINE cftfax_g(n, ifax, trigs)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER n
      INTEGER, DIMENSION(13) :: ifax
      REAL(rprec), DIMENSION(*) :: trigs
!DEC$ IF .NOT.DEFINED(CRAY) .OR. DEFINED(LONESTAR) .OR. DEFINED(MCURIE)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: k
C-----------------------------------------------
c
c on input     n
c               the length of each complex transform to be performed
c
c               n must be greater than 1 and contain no prime
c               factors greater than 5.
c
c on output    ifax
c               ifax(1)
c                 the number of factors chosen or -99 in CASE of error
c               ifax(2) thru ifax( ifax(1)+1 )
c                 the factors of n in the followin order:  appearing
c                 first are as many factors of 4 as can be obtained.
c                 subsequent factors are primes, and appear in
c                 ascending order, except for multiple factors.
c
c              trigs
c               2n SIN and COS values for USE by the transform routine
c
      CALL fact_g (n, ifax)
      k = ifax(1)
      IF (k<1 .or. ifax(k+1)>5) ifax(1) = -99
      IF (ifax(1) .le. 0) STOP 'IFAX(1) <= 0 in CFTFAX'
      CALL cftrig_g (n, trigs)
!DEC$ ELSE
      CALL cftfax (n, ifax, trigs)
!DEC$ ENDIF
      END SUBROUTINE cftfax_g
