      SUBROUTINE cftrig_g (n, trigs)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER n
      REAL(rprec), DIMENSION(*) :: trigs
!DEC$ IF .NOT.DEFINED (CRAY) .OR. DEFINED(LONESTAR) .OR. DEFINED(MCURIE)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: one = 1, two = 2, p5 = .5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: l, i
      REAL(rprec) :: pi, del, angle
C-----------------------------------------------
      pi = two*ASIN(one)
      del = (pi + pi)/n
      l = n + n
      DO 10 i=1,l,2
        angle=(p5*del)*(i-1)
        trigs(i)=COS(angle)
        trigs(i+1)=SIN(angle)
   10 CONTINUE
!DEC$ ELSE
      CALL cftrig (n, trigs)
!DEC$ ENDIF
      END SUBROUTINE cftrig_g
