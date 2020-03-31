      SUBROUTINE sort_data (x, index_array, n)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
c-----------------------------------------------
      INTEGER n
      INTEGER, DIMENSION(n) :: index_array
      REAL(rprec), DIMENSION(n) :: x
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, j
      INTEGER, DIMENSION(1) :: isamax
      REAL(rprec), DIMENSION(n) :: dumx
C-----------------------------------------------
!
!       RETURNS INDEX(I) ARRAY, SO THAT X(INDEX(I)) IS SORTED, I=1,N
!       RETURNS Xin(INDEX(I)) = Xout(I)
!

      DO i = n, 1, -1
         isamax = MAXLOC(ABS(x))
         j = isamax(1)
         dumx(i) = x(j)
         x(j) = 0.
         index_array(i) = j
      END DO

      x = dumx

      END SUBROUTINE sort_data
