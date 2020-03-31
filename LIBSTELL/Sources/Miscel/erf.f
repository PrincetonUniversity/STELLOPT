      FUNCTION ERF(x)
      USE stel_kinds, ONLY: rprec
!-------------------------------------------------------------------------------
!ERF evaluates the error function erf(x)
!References:
!  M.Abramowitz, L.A.Stegun, Handbook of Math. Functions, p. 299
!  W.A.Houlberg 7/2003
!Comments:
!  The error function can consume a lot of time deep in the multiple species
!    loop so a very efficient calculation is called for
!  Time consumption is much more critical than accuracy as suggested by T.Amano
!  A three term expansion from Abramowitz is not sufficiently accurate because
!    it generates viscosities with singularities in the vicinity of the BP-PS
!    transition for ions
!-------------------------------------------------------------------------------

!Declaration of input variables
      REAL(rprec), INTENT(IN) :: x                   !argument of error function [-]

!Declaration of output variables
      REAL(rprec) ::  ERF                            !value of error function [-]

!Declaration of local variables
      REAL(rprec) :: t

      REAL(rprec), PARAMETER :: one=1, p=0.3275911_rprec,
     1                          a1=0.254829592_rprec,
     2                          a2=-0.284496736_rprec,     
     3                          a3=1.421413741_rprec,
     4                          a4=-1.453152027_rprec,
     5                          a5=1.061405429_rprec
!-------------------------------------------------------------------------------
!Apply fit
!-------------------------------------------------------------------------------
      t=one/(one + p*x)
      ERF = one - ( a1       
     1          + ( a2
     2          + ( a3
     3          + ( a4 
     4          +   a5*t)*t)*t)*t)*t*EXP(-x**2)

      END FUNCTION ERF
