      SUBROUTINE legendre_to_power(n, a_inv, b_inv, tc, ac)
      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN):: n
      REAL(rprec), DIMENSION(0:n), INTENT(IN):: tc
      REAL(rprec), DIMENSION(0:n), INTENT(OUT):: ac
      REAL(rprec), DIMENSION(0:n,0:n), INTENT(IN):: a_inv, b_inv
!---------------------------------------------------------------------
      INTEGER:: i, j, k
!---------------------------------------------------------------------
!     Given the following notation:
!
!           AC == (ac(1), ...ac(n))==>  vector of coefficients for
!                 power series in [0,1]
!           TC == (tc(1), ...tc(n))==>  vector of coefficients for
!                 Legendre series in [-1,1]
!     THEN:
!                        AC = TC* A_INV * B_INV
!----------------------------------------------------------------------
      DO i = 0, n
         ac(i) = 0
         DO j= 0, n
           DO k = 0, n
              ac(i) = ac(i) + tc(j) * a_inv(j,k) * b_inv(k,i)
           ENDDO
         ENDDO
      ENDDO

      END SUBROUTINE legendre_to_power
