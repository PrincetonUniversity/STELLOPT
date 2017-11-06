      SUBROUTINE pofs(re, ro, ns, rthom, itse)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER ns, itse
      REAL(rprec), DIMENSION(ns) :: re, ro
      REAL(rprec), DIMENSION(*) :: rthom
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: j, js, i
      REAL(rprec), DIMENSION(ns) :: s
      REAL(rprec) :: sqrjs, rt2, ds
C-----------------------------------------------
!
!       Interpolate Rmid = Re + Ro to get R(s)
!
!       THIS ROUTINE INTERPOLATES THE RTHOM "S-SPACE" ARRAY
!       ONTO THE INSTANTANEOUS RMID ARRAY
!       ON INPUT, RTHOM IS THE S-VALUE ARRAY FOR THOMPSON DATA
!       ON OUTPUT,RTHOM IS THE CORRESPONDING (THETA=0) RMID ARRAY
!
      DO j = 1, ns
         s(j) = REAL(j - 1,rprec)/(ns - 1)
      END DO

      js = 1
      DO i = 1, itse
         rt2 = rthom(i)
  100    CONTINUE
         IF (rt2.ge.s(js) .and. rt2.le.s(js+1)) THEN
            ds = (rt2 - s(js))/(s(js+1)-s(js))
            sqrjs = SQRT(rt2)
            rthom(i) = re(js) + (re(js+1)-re(js))*ds + sqrjs*
     1         (ro(js)+(ro(js+1)-ro(js))*ds)
         ELSE
            js = js + 1
            IF (js < ns) GOTO 100
         ENDIF
      END DO

      END SUBROUTINE pofs
