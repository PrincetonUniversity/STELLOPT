      SUBROUTINE power_to_legendre(n, a, b, ac, tc)
      USE stel_kinds
      IMPLICIT NONE
      INTEGER, INTENT(IN):: n
      REAL(rprec), DIMENSION(0:n), INTENT(IN):: ac
      REAL(rprec), DIMENSION(0:n), INTENT(OUT):: tc
      REAL(rprec), DIMENSION(0:n,0:n), INTENT(IN):: a, b
      INTEGER:: i, j, k
!------------------------------------------------------------------
!      Given the following notation:
!
!           TC == (tc(1), ...tc(n))==>  vector of coefficients
!                 for Legendre series in [-1,1]
!           AC == (ac(1), ...ac(n))==>  vector of coefficients
!                 for power series in [0,1]
!      THEN:
!                        TC = AC* B * A
!------------------------------------------------------------------
      DO i = 0, n
         tc(i) = 0
         DO j= 0, n
           DO k = 0, n
             tc(i) = tc(i) + ac(j) * b(j,k) * a(k,i)
           END DO
         END DO
      END DO

      END SUBROUTINE power_to_legendre
