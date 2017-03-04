      SUBROUTINE ludcmp(a, n, np, indx, d)
!  (c) copr. 1986-92 numerical recipes software

      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER n, np
      REAL(rprec) d
      INTEGER, DIMENSION(n) :: indx
      REAL(rprec), DIMENSION(np,np) :: a
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: nmax = 500
      REAL(rprec), PARAMETER :: tiny_no = 1.0e-20_dp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, imax, j, k
      REAL(rprec) :: aamax, dum, sumn
      REAL(rprec), DIMENSION(nmax) :: vv
!-----------------------------------------------
!
      d = 1
      vv = 0
      DO i = 1, n
         aamax = MAXVAL(ABS(a(i,:n)))
         IF (aamax == 0._dp) THEN
            WRITE (*, '(2a)') 'PAUSE ', 'singular matrix in ludcmp'
         ELSE
            vv(i) = 1._dp/aamax
         ENDIF
      END DO
      DO j = 1, n
         DO i = 1, j - 1
            sumn = a(i,j)
            sumn = sumn - SUM(a(i,:i-1)*a(:i-1,j))
            a(i,j) = sumn
         END DO
         aamax = 0
         DO i = j, n
            sumn = a(i,j)
            sumn = sumn - SUM(a(i,:j-1)*a(:j-1,j))
            a(i,j) = sumn
            dum = vv(i)*ABS(sumn)
            IF (dum >= aamax) THEN
               imax = i
               aamax = dum
            ENDIF
         END DO
         IF (j /= imax) THEN
            DO k = 1, n
               dum = a(imax,k)
               a(imax,k) = a(j,k)
               a(j,k) = dum
            END DO
            d = -d
            vv(imax) = vv(j)
         ENDIF
         indx(j) = imax
         IF (a(j,j) == 0._dp) a(j,j) = tiny_no
         IF (j /= n) THEN
            dum = 1._dp/a(j,j)
            a(j+1:n,j) = a(j+1:n,j)*dum
         ENDIF
      END DO

      END SUBROUTINE ludcmp
