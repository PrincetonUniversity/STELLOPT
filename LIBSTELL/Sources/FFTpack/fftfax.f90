      SUBROUTINE fftfax_g(n, ifax, trigs)
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
      INTEGER :: mode=3, i
!-----------------------------------------------
!
! mode 3 is used for REAL/half-complex transforms.  it is possible
! to do complex/complex transforms with other values of mode, but
! documentation of the details were not available when this routine
! was written.
!
      CALL fax (ifax, n, mode)
      i = ifax(1)
      IF (ifax(i+1)>5 .or. n<=4) ifax(1) = -99
      IF (ifax(1) <= 0) STOP 'IFAX(1) <= 0 in fftfax'
      CALL fftrig_g (trigs, n, mode)
#else
      CALL fftfax (n, ifax, trigs)
#endif
      END SUBROUTINE fftfax_g
