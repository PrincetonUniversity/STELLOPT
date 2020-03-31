      FUNCTION ga_evaluate(fcn, nopt, fvec, n, x, iflag, nfev)
      USE stel_kinds
      INTEGER :: n, nfev
      REAL(rprec) :: ga_evaluate
      REAL(rprec), DIMENSION(n) :: x
      INTEGER :: nopt
      REAL(rprec), DIMENSION(nopt) :: fvec
      EXTERNAL  fcn

      INTEGER :: iflag

      CALL fcn(nopt, n, x, fvec, iflag, nfev)
      ga_evaluate = -SUM(fvec(:nopt)**2)

      END FUNCTION ga_evaluate
