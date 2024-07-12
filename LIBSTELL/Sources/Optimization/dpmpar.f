      FUNCTION dpmpar (i)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: i
      REAL(rprec) :: dpmpar
C-----------------------------------------------
c
c     FUNCTION dpmpar
c
c     this FUNCTION provides single (double) precision machine parameters
c     when the appropriate set of data statements is activated (by
c     removing the c from column 1) and ALL other data statements are
c     rendered inactive. most of the PARAMETER values were obtained
c     from the corresponding bell laboratories port library FUNCTION.
c
c     the FUNCTION statement is
c
c       FUNCTION dpmpar(i)
c
c     WHERE
c
c       i is an integer input variable set to 1, 2, or 3 which
c         selects the desired machine parameter. if the machine has
c         t base b digits and its smallest and largest exponents are
c         emin and emax, respectively, then these parameters are
c
c         dpmpar(1) = b**(1 - t), the machine precision,
c
c         dpmpar(2) = b**(emin - 1), the smallest magnitude,
c
c         dpmpar(3) = b**emax*(1 - b**(-t)), the largest magnitude.
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     SINGLE PRECISION
c     data rmach(1)/1.490116100E-8/
c     data rmach(2)/1.4693679000E-30/
c     data rmach(3)/1.701411800E+30/
c     DOUBLE PRECISION - IEEE
c     data dmach(1) /2.22044604926d-16/
c     data dmach(2) /2.22507385852d-308/
c     data dmach(3) /1.79769313485d+308/
c     modified for f90 (sph, august 1997)
c
      SELECT CASE(i)
      CASE(:1)
        dpmpar = EPSILON(dpmpar)      !2.22044604926e-16_dp
      CASE(2)
        dpmpar = TINY(dpmpar)         !2.22507385852e-308_dp
      CASE(3:)
        dpmpar = HUGE(dpmpar)         !1.79769313485e+308_dp
      END SELECT

      END FUNCTION dpmpar
