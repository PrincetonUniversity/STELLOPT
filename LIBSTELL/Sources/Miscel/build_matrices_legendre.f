      SUBROUTINE build_matrices_legendre(n, a, b, a_inv, b_inv)
      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: n
      REAL(rprec), INTENT(OUT), DIMENSION(0:n,0:n) :: a, b,
     1  a_inv, b_inv
!-----------------------------------------------
      INTEGER :: idx, i, j
      REAL(rprec) :: factorial
      REAL(rprec), DIMENSION(13), PARAMETER :: norm_legend =
     1    ( / 1, 1, 3, 5, 35, 63, 231, 429, 6435,
     2        12155, 46189, 88179, 676039 / )
      REAL(rprec), DIMENSION(13), PARAMETER ::  norm_inv_legend =
     1    ( / 1, 1, 2, 2, 8, 8, 16, 16, 128,
     2        128, 256, 256, 1024 / )
      REAL(rprec), DIMENSION(13*13), PARAMETER :: b_legend =
     1  ( / 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     2  0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     3  1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     4  0, 3, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     5  7, 0, 20, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0,
     6  0, 27, 0, 28, 0, 8, 0, 0, 0, 0, 0, 0, 0,
     7  33, 0, 110, 0, 72, 0, 16, 0, 0, 0, 0, 0, 0,
     8  0, 143, 0, 182, 0, 88, 0, 16, 0, 0, 0, 0, 0,
     9  715, 0, 2600, 0, 2160, 0, 832, 0, 128, 0, 0, 0, 0,
     A  0, 3315, 0, 4760, 0, 2992, 0, 960, 0, 128, 0, 0, 0,
     B  4199, 0, 16150, 0, 15504, 0, 7904, 0, 2176, 0,
     C  256, 0, 0,
     D  0, 20349, 0, 31654, 0, 23408, 0, 10080, 0, 2432,
     E  0, 256, 0,
     F  52003, 0, 208012, 0, 220248, 0, 133952, 0, 50048,
     G  0, 10752, 0, 1024 / )
      REAL(rprec), DIMENSION(13*13), PARAMETER :: b_legend_inv =
     1  ( / 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     2  0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     3  -1, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     4  0, -3, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     5  3, 0, -30, 0, 35, 0, 0, 0, 0, 0, 0, 0, 0,
     6  0, 15, 0, -70, 0, 63, 0, 0, 0, 0, 0, 0, 0,
     7  -5, 0, 105, 0, -315, 0, 231, 0, 0, 0, 0, 0, 0,
     8  0, -35, 0, 315, 0, -693, 0, 429, 0, 0, 0, 0, 0,
     9  35, 0, -1260, 0, 6930, 0, -12012, 0, 6435, 0, 0,
     A  0, 0,
     B  0, 315, 0, -4620, 0, 18018, 0, -25740, 0, 12155,
     C  0, 0, 0,
     D  -63, 0, 3465, 0, -30030, 0, 90090, 0, -109395, 0,
     E  46189, 0, 0,
     F  0, -693, 0, 15015, 0, -90090, 0, 218790, 0, -230945,
     G  0, 88179, 0,
     H  231, 0, -18018, 0, 225225, 0, -1021020, 0, 2078505,
     I  0, -1939938, 0, 676039 / )

!--------------------------------------------------------------------------------------------------
!      P_n(y):: Legendre polynomial in [-1,1]. Matrix coefficients from
!      "Handbook of Mathematical Functions", Abramowitz and Stegun, 1973, p.795.
!      If higher degree polyn. needed, then data set must be extended.
!
!      DESCRIPTION OF MATRICES:
!
!      A, A_inv: convert legendre Pol. to and from y-powers in [-1,1]
!
!              norm_legend(m)*y**m = sum_n=0^m (b_legend{mn} L_n)
!              norm_legend_inv(m)*L_m = sum_n=0^N (b_legend_inv{mn} y**n)
!
!           ==>   A_{mn} = b_legend_{mn}/norm_legend(m)                IF  m<=n  or  0. otherwise
!           ==>   A_inv_{mn} = b_legend_inv_{mn}/norm_inv_legend(m)    IF  m<=n  or  0. otherwise
!
!      B, B_inv: convert y-powers in [-1,1] to and from x-powers in [0,1]
!
!             x = (y + 1)/2;   x**m =  sum_n=0^m  2**(-m) m!/(n!(m-n)!)  y**n
!             y = 2*x - 1;     y**m =  sum_n=0^m  2**n (-1)**(m-n) m!/(n!(m-n)!)  x**n
!
!           ==>   B_{mn} =  2**(-m) m!/(n!(m-n)!)                     IF  m<=n  or  0. otherwise
!           ==>   B_inv_{mn} =  2**n (-1)**(m-n) m!/(n!(m-n)!)        IF  m<=n  or  0. otherwise
!
!
!      NOTICE that:   A_inv = (A)**(-1);  B_inv = (B)**(-1)
!--------------------------------------------------------------------------------------------------

      IF (n > 12) STOP 'N(legendre) CANNOT be larger than 12!!!'

      a = 0
      a_inv = 0
      b = 0
      b_inv = 0

       DO i = 0, n
         idx = 13*i
         DO j = 0, i
           a(i,j) = b_legend(idx+j+1)/norm_legend(i+1)
           a_inv(i,j) = b_legend_inv(idx+j+1)/norm_inv_legend(i+1)
           b(i,j) = 2._dp**(-i) * factorial(i)/
     1       (factorial(j)*factorial(i-j))
           b_inv(i,j) = 2._dp**j*(-1._dp)**(i-j)* factorial(i)/
     1       (factorial(j)*factorial(i-j))
         END DO
       END DO

      END SUBROUTINE build_matrices_legendre
