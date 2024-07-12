

!========================================================================
!  The following are supplemental subroutines from various libraries
!========================================================================

      SUBROUTINE SIMPUN(XX, FX, NX, AX)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      use stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER NX
      REAL(rprec), DIMENSION(NX) :: XX, FX, AX
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero=0, two=2, six=6
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: II, IX, IC
      REAL(rprec) :: D1, D2, D3, A2, A3
C-----------------------------------------------
c---
c  II>0 then ax(i) = Integral(f(x)dx) from x=xx(1) to x=xx(i). Here f(x)=fx(i)
c  II<0 then ax(i) = Integral(f(x)dx) from x=xx(nx) to x=xx(i)
c  PROGRAM AUTHOR      J. BARISH,
c  COMPUTING TECHNOLOGY CENTER, UNION CARBIDE CORP., NUCLEAR DIV.,
c  OAK RIDGE, TENN.
c---
c
      II = 1
      IF (II < 0) GO TO 30
      AX(1) = ZERO
      DO IX = 2, NX, 2
         D1 = XX(IX) - XX(IX-1)
         AX(IX) = AX(IX-1) + D1/TWO*(FX(IX)+FX(IX-1))
         IF (NX .eq. IX) EXIT
         D2 = XX(IX+1) - XX(IX-1)
         D3 = D2/D1
         A2 = D3/SIX*D2**2/(XX(IX+1)-XX(IX))
         A3 = D2/TWO - A2/D3
         AX(IX+1)=AX(IX-1)+(D2-A2-A3)*FX(IX-1)+A2*FX(IX)+A3*FX(IX+1)
      END DO
   20 CONTINUE
      RETURN
   30 CONTINUE
      AX(NX) = ZERO
      DO IX = 2, NX, 2
         IC = NX + 1 - IX
         D1 = XX(IC+1) - XX(IC)
         AX(IC) = AX(IC+1) + D1/TWO*(FX(IC+1)+FX(IC))
         IF (NX .eq. IX) GO TO 20
         D2 = XX(IC+1) - XX(IC-1)
         D3 = D2/(XX(IC)-XX(IC-1))
         A2 = D3/SIX*D2**2/D1
         A3 = D2/TWO - A2/D3
         AX(IC-1)=AX(IC+1)+(D2-A2-A3)*FX(IC-1)+A2*FX(IC)+A3*FX(IC+1)
      END DO

      END SUBROUTINE SIMPUN
