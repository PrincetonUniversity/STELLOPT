
MODULE neo_rhscur
  INTERFACE rhs_cur
     SUBROUTINE rhs_cur1(phi,y,dery)
       USE neo_precision
       REAL(dp),                             INTENT(in)  ::  phi
       REAL(dp), DIMENSION(:), TARGET,       INTENT(in)  ::  y
       REAL(dp), DIMENSION(SIZE(y)), TARGET, INTENT(out) ::  dery
     END SUBROUTINE rhs_cur1
  END INTERFACE
END MODULE neo_rhscur
