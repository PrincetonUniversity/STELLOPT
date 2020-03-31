
MODULE neo_rhsbo
  INTERFACE rhs_bo
     SUBROUTINE rhs_bo1(phi,y,dery)
       USE neo_precision
       REAL(dp),                             INTENT(in)  ::  phi
       REAL(dp), DIMENSION(:),               INTENT(in)  ::  y
       REAL(dp), DIMENSION(SIZE(y)), TARGET, INTENT(out) ::  dery
     END SUBROUTINE rhs_bo1
  END INTERFACE
END MODULE neo_rhsbo
