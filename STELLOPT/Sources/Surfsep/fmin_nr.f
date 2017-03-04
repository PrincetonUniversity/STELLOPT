      FUNCTION fmin_nr (x)
!  (c) copr. 1986-92 numerical recipes software

      USE stel_kinds
      USE newtv
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(*) :: x
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec) :: fmin_nr, r1
!-----------------------------------------------
!
!     USEs funcv
!
      CALL funcv (nn, x, fvec)
      r1 = DOT_PRODUCT(fvec(:nn),fvec(:nn))
      fmin_nr = r1/2

      END FUNCTION fmin_nr
