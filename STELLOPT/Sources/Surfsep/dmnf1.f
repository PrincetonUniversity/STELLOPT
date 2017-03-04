      SUBROUTINE dmnf1(u, v, SUMr, SUMz)
      USE stel_kinds
      USE cmnf1
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec) u, v, SUMr, SUMz
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec) :: pi2
!-----------------------------------------------
!
!     constants
      pi2 = 4*ASIN(1._dp)
!
      SUMr = SUM(rmn1(:nmn1)*COS(pi2*(xm1(:nmn1)*u+xn1(:nmn1)*v)))
      SUMz = SUM(zmn1(:nmn1)*SIN(pi2*(xm1(:nmn1)*u+xn1(:nmn1)*v)))
!
      END SUBROUTINE dmnf1
