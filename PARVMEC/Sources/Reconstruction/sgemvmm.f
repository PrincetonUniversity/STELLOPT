      SUBROUTINE sgemvmm(amat_i, amat_p, amatsq, b, bsq, wten,
     1   mdata, niota, npres, nots)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER mdata, niota, npres, nots
      REAL(rprec), DIMENSION(niota,mdata) :: amat_i
      REAL(rprec), DIMENSION(npres,mdata) :: amat_p
      REAL(rprec), DIMENSION(nots,nots) :: amatsq
      REAL(rprec), DIMENSION(mdata) :: b
      REAL(rprec), DIMENSION(nots) :: bsq, wten
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0.0_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, j, ioff, joff
C-----------------------------------------------

      amatsq = zero
      bsq    = zero

!
!       INITIALIZE IOTA, PRESSURE DIAGONAL ELEMENTS (ALREADY 'SQUARED')
!
!
!       COMPUTE LOWER DIAGONAL HALF OF SQUARE OF MATRIX
!       A(trans)*A and A(trans)*B
!       BY ADDING CONTRIBUTIONS FROM EXTERNAL MAGNETICS SIGNALS
!

!
!       FIRST UPPER LEFT NIOTA X NIOTA BLOCK
!
      DO i = 1, niota
         bsq(i) = bsq(i) + SUM(b*amat_i(i,:))
         DO j = 1, i
            amatsq(i,j) = amatsq(i,j) + SUM(amat_i(i,:)*amat_i(j,:))
         END DO
      END DO

!
!       LOWER NPRES X NIOTA BLOCK, NPRES X NPRES BLOCK
!
      DO i = 1, npres
         ioff = i + niota
         bsq(ioff) = bsq(ioff) + SUM(b*amat_p(i,:))
         DO j = 1, niota
            amatsq(ioff,j) = amatsq(ioff,j) +
     1                       SUM(amat_p(i,:)*amat_i(j,:))
         END DO
         DO j = 1, i
            joff = j + niota
            amatsq(ioff,joff) = amatsq(ioff,joff) +
     1                          SUM(amat_p(i,:)*amat_p(j,:))
         END DO
      END DO

      DO i = 1, nots
         wten(i) = amatsq(i,i)
         amatsq(1:i-1,i) = amatsq(i,1:i-1)
      END DO


      END SUBROUTINE sgemvmm
