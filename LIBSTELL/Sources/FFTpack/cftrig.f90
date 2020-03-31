      SUBROUTINE cftrig_g (n, trigs)
      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER n
      REAL(rprec), DIMENSION(*) :: trigs
#if !defined(CRAY) || defined(LONESTAR) || defined(MCURIE)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: one = 1, two = 2, p5 = .5_dp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: l, i
      REAL(rprec) :: pi, del, angle
!-----------------------------------------------
      pi = two*ASIN(one)
      del = (pi + pi)/n
      l = n + n
      DO 10 i=1,l,2
        angle=(p5*del)*(i-1)
        trigs(i)=COS(angle)
        trigs(i+1)=SIN(angle)
   10 CONTINUE
#else
      CALL cftrig (n, trigs)
#endif
      END SUBROUTINE cftrig_g
