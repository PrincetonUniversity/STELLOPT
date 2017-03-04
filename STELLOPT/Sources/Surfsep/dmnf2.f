      SUBROUTINE dmnf2(u, v, SUMr, SUMz)
      USE stel_kinds
      USE cmnf2
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
!
!     constants
      pi2 = 4*ASIN(1._dp)
!
      SUMr = SUM(rmn2(:nmn2)*COS(pi2*(xm2(:nmn2)*u+xn2(:nmn2)*v)))
      SUMz = SUM(zmn2(:nmn2)*SIN(pi2*(xm2(:nmn2)*u+xn2(:nmn2)*v)))

      END SUBROUTINE dmnf2
