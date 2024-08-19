      SUBROUTINE rz_surf(theta, phi, rw, zw, kmodes, rmn, zmn, m_num,
     1   n_num, nfp)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE stel_constants
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER kmodes, nfp
      INTEGER, DIMENSION(kmodes) :: m_num, n_num
      REAL(rprec) :: theta, phi, rw, zw
      REAL(rprec), DIMENSION(kmodes) :: rmn, zmn
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i
!-----------------------------------------------

      rw = zero
      zw = zero

! 2/16/98 WHM The NESCOIL sign convention is used here

      DO i=1,kmodes
         rw = rw + rmn(i)*COS(m_num(i)*theta+nfp*n_num(i)*phi)
         zw = zw + zmn(i)*SIN(m_num(i)*theta+nfp*n_num(i)*phi)
      END DO

      END SUBROUTINE rz_surf
