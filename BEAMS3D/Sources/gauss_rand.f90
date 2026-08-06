      SUBROUTINE gauss_rand(n,x)
!  Subroutine to fill an array with gaussian random variables of
!  unit variance. Pass size of array, so that don't have to 
!  fuss with interfaces

      USE stel_kinds, ONLY: rprec

      INTEGER, INTENT(in) :: n
      REAL(rprec), DIMENSION(n), INTENT(inout) :: x

      INTEGER :: np1o2, i, j
      REAL(rprec), DIMENSION(2) :: v
      REAL(rprec) :: rsq, fac

!  start of executable code
      np1o2 = (n + 1) / 2
      j = 0
      DO i = 1,np1o2
         j = j + 1
100      CALL RANDOM_NUMBER(v)
         v = 2. * v - 1.
         rsq = v(1) * v(1) + v(2) * v(2)
         IF ((rsq .ge. 1.) .or. (rsq .eq. 0.)) GO TO 100
         fac = SQRT(-2. * LOG(rsq) / rsq)
         x(j) = v(1) * fac
         j = j + 1
         IF (j .le. n) THEN
            x(j) = v(2) * fac
         ENDIF
      END DO
      RETURN
      END SUBROUTINE gauss_rand
