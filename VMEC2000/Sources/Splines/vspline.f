      MODULE vspline
      USE vsvd0
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, DIMENSION(1) :: jspmin
      INTEGER :: iknots
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: hthom, ythom,
     1   y2thom, pknots, hstark, y2stark, ystark, sknots
C-----------------------------------------------
      END MODULE vspline
