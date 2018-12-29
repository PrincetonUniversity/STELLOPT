      MODULE legendre_profile
! A collection of functions related to Legendre polynomials profiles
! Based loosely on the originally code that appeared in the version 1
! of the stellopt. It was ported into STELLOPTV2, and since then, VMEC,
! LIBSTELL, and the STELLOPT codes have been reorgized significantly.
! Much of the original code has been migrated to this file, but not
! all of it has been tested for backwards-compatibility with 
! the original STELLOPT code.
! Development: legedre_ip: legendre polynomials for current density
!                          integrated Leg.polys for total current out
! To do:  legendre_I, legendre_p, legedre_p0, legendre_iota
! 
      USE stel_kinds        
      IMPLICIT NONE

      INTEGER :: n_leg     
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: a_leg, b_leg,
     1   a_leg_inv, b_leg_inv, c_leg_int 
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: tc, ti, tm  

      CONTAINS
     
      SUBROUTINE build_matrices_legendre(n, a, b, a_inv, b_inv)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: n
      REAL(rprec), INTENT(OUT), DIMENSION(0:n,0:n) :: a, b
      REAL(rprec), INTENT(OUT), DIMENSION(0:n,0:n) :: a_inv, b_inv
!-----------------------------------------------
      INTEGER :: idx, idx_2, i, j
      REAL(rprec) :: factorial
      REAL(rprec), DIMENSION(13), PARAMETER :: norm_legend = 
     1     ( / 1, 1, 3, 5, 35, 63, 231, 429, 6435, 
     2         12155, 46189, 88179, 676039 / )
      REAL(rprec), DIMENSION(13), PARAMETER ::  norm_inv_legend = 
     1     ( / 1, 1, 2, 2, 8, 8, 16, 16, 128, 
     2         128, 256, 256, 1024 / )
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
!         
!--------------------------------------------------------------------------------------------------

      IF (n > 12) STOP 'N(legendre) CANNOT be larger than 12!!!'

      DO i = 0, n
        idx = 13*i
        DO j = 0, i
          a(i,j) = b_legend(idx+j+1)/norm_legend(i+1)
          a_inv(i,j) = b_legend_inv(idx+j+1)/norm_inv_legend(i+1)
          b(i,j) = 2._dp**(-i) * factorial(i)/
     1      (factorial(j)*factorial(i-j))
          b_inv(i,j) = 2._dp**j*(-1._dp)**(i-j)* factorial(i)/
     1      (factorial(j)*factorial(i-j))
        END DO
      END DO

      END SUBROUTINE build_matrices_legendre


      SUBROUTINE build_matrices_legendre_int(n, c)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: n
      REAL(rprec), INTENT(OUT), DIMENSION(0:n,0:(n+1)) :: c
!-----------------------------------------------
      INTEGER :: idx, i, j
!--------------------------------------------------------------------------------------------------
!      C_legend_int: Pre-evaluated integrals of the Legendre polynomials
!      from 0<=s<=x. Coefficients calculated with Mathematica
!           {Integrate[LegendreP[0, 2*x - 1], x],
!            Integrate[LegendreP[1, 2*x - 1], x],
!            Integrate[LegendreP[2, 2*x - 1], x],
!            Integrate[LegendreP[3, 2*x - 1], x],
!            Integrate[LegendreP[4, 2*x - 1], x],
!            Integrate[LegendreP[5, 2*x - 1], x],
!            Integrate[LegendreP[6, 2*x - 1], x],
!            Integrate[LegendreP[7, 2*x - 1], x],
!            Integrate[LegendreP[8, 2*x - 1], x],
!            Integrate[LegendreP[9, 2*x - 1], x],
!            Integrate[LegendreP[10, 2*x - 1], x],
!            Integrate[LegendreP[11, 2*x - 1], x],
!            Integrate[LegendreP[12, 2*x - 1], x]}
!           CoefficientList[%, x]
!         
!--------------------------------------------------------------------------------------------------
      REAL(rprec), DIMENSION(13*14), PARAMETER :: c_legend_int = 
     1  ( / 0,  1 ,  0,     0,      0,     0,    0, 0, 0, 0, 0, 0, 0, 0,
     2      0, -1 ,  1,     0,      0,     0,    0, 0, 0, 0, 0, 0, 0, 0,
     3      0,  1 , -3,     2,      0,     0,    0, 0, 0, 0, 0, 0, 0, 0,
     4      0, -1 ,  6,   -10,      5,     0,    0, 0, 0, 0, 0, 0, 0, 0,
     5      0,  1, -10,    30,    -35,    14,    0, 0, 0, 0, 0, 0, 0, 0,
     6      0, -1,  15,   -70,    140,  -126,   42, 0, 0, 0, 0, 0, 0, 0,
     7      0,  1, -21,   140,   -420,   630, -462, 132, 
     8              0, 0, 0, 0, 0, 0,
     9      0, -1,  28,  -252,   1050, -2310,    2772,   -1716,
     A            429,        0, 0, 0, 0, 0,
     B      0,  1, -36,   420,  -2310,  6930,  -12012,   12012,
     C          -6435,     1430,        0, 0, 0, 0,
     D      0, -1,  45,  -660,   4620, -18018,   42042,  -60060,
     E          51480,   -24310,     4862,       0, 0, 0,
     F      0,  1, -55,   990,  -8580,  42042, -126126,  240240,
     G        -291720,   218790,   -92378,   16796,        0, 0,
     H      0, -1,  66, -1430,  15015, -90090,  336336, -816816,
     I        1312740, -1385670,   923780, -352716,    58786,      0,
     J      0,  1, -78,  2002, -25025, 180180, -816816, 2450448,
     K       -4988412,  6928350, -6466460, 3879876, -1352078, 208012 / )


      IF (n > 12) STOP 'N(legendre) CANNOT be larger than 12!!!'

      ! Initialize c to zero
      c = 0
      DO i = 0, n
        idx = 14*i
        DO j = 0, (i+1)
          c(i,j) = c_legend_int(idx+j+1)
        END DO
      END DO

      END SUBROUTINE build_matrices_legendre_int


      SUBROUTINE legendre_poly_int(n_leg, c_leg_int, xx, yy, aa, 
     1                              i_start, i_end)

      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(IN)  ::  n_leg, i_start, i_end
      REAL(rprec), INTENT(IN)  ::  xx
      REAL(rprec), INTENT(OUT) :: yy
      REAL(rprec), DIMENSION(0:n_leg, 0:n_leg+1) :: c_leg_int
      REAL(rprec), DIMENSION(0:n_leg), INTENT(IN) :: aa
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER  ::  i, j, k
      REAL(rprec) :: sum_1, sum_2
C-----------------------------------------------
!       Legendre polynomial integration routine.
!       n_gen:  number of legendre polynomials
!       c_leg_int: matrix of integrated Legendre poly, in terms of
!       powers`
!       aa: array of coefficients
!       computed from call to SPLINE
!       xx : value at which yy = F(X) is to be computed from Legendre
!       polynomials
!       i_start: The index of the first element that is to be used for
!       the Legendre polynomial series summation.
!       i_end: The index of the last element to be used in the Legendre
!       summation.
!---------------------------------------------------------------------
!     Given the following notation:
!
!           aa: (aa(0), ...aa(n_gen-1))==>  vector of coefficients for
!                 integrated legendre_ip series in [0,1]
!           c_leg_int: c_mn: Each row: The coefficients for a 
!                 polyinomial expansion of the integration of the nth
!                 Legendre polynomial[-1,1], evaluatad at [2*xx-1]
!     THEN:
!           pcurr = sum_i (aa(i) * sum_j (c_leg_int(i,j)*xx**j)
!     Coefficients for c_leg_int calculated via Mathematica.
!----------------------------------------------------------------------
!!  Legendre polynomials for I-prime
!!  I(s) = Sum(i,0,n)[ac(i)*Integra(LegendreP(i, 2*s-1), ds, 0<s<1)~s
!!  a. Curtor sets the overall total integrated current.
!!  b. The integral of 0th Legendre polynomial evaluated at (s=1) is 1
!!  c. The integral of all higher-order (>0) Legendre polynomials,
!!     evaluated at (s=1) is 0.
!!  So, the 0th component should be set to 1 and the normalization (curtor)
!!      sets the overall current, or, set the 0th component to curtor.

      yy = 0

      DO i = i_start, i_end
         sum_2 = 0
         DO j = i_end+1, i_start, -1
            sum_2 = sum_2 + c_leg_int(i,j) * xx**j
         END DO
         yy = yy + aa(i) * sum_2
      END DO

      END SUBROUTINE legendre_poly_int


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


      END MODULE legendre_profile                            

