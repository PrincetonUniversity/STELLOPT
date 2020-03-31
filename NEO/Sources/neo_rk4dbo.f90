
MODULE neo_rk4dbo
  INTERFACE rk4d_bo
     SUBROUTINE rk4d_bo1(x,y,h)
       USE neo_precision
       REAL(dp),               INTENT(inout) ::  x
       REAL(dp),               INTENT(in)    ::  h
       REAL(dp), DIMENSION(:), INTENT(inout) ::  y
     END SUBROUTINE rk4d_bo1
  END INTERFACE
END MODULE neo_rk4dbo
