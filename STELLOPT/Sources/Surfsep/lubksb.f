      SUBROUTINE lubksb(a, n, np, indx, b)
!  (c) copr. 1986-92 numerical recipes software

      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER n, np
      INTEGER, DIMENSION(n) :: indx
      REAL(rprec), DIMENSION(np,np) :: a
      REAL(rprec), DIMENSION(n) :: b
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, ii, ll
      REAL(rprec) :: sumn
!-----------------------------------------------
!
!
      ii = 0
      DO i = 1, n
         ll = indx(i)
         sumn = b(ll)
         b(ll) = b(i)
         IF (ii .ne. 0) THEN
            sumn = sumn - SUM(a(i,ii:i-1)*b(ii:i-1))
         ELSE IF (sumn .ne. 0._dp) THEN
            ii = i
         ENDIF
         b(i) = sumn
      END DO
      DO i = n, 1, -1
         sumn = b(i)
         sumn = sumn - SUM(a(i,i+1:n)*b(i+1:n))
         b(i) = sumn/a(i,i)
      END DO

      END SUBROUTINE lubksb
