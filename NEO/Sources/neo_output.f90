
MODULE neo_output
  USE neo_precision
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE ::  epspar
  REAL(kind=dp)                            ::  epstot,ctrone,ctrtot
  REAL(kind=dp)                            ::  bareph,barept,drdpsi
  REAL(kind=dp)                            ::  yps
  INTEGER                                  ::  nintfp
  INTEGER                                  ::  ierr
  REAL(kind=dp)                            ::  lambda_b
  REAL(kind=dp)                            ::  lambda_b1, lambda_b2
  REAL(kind=dp)                            ::  lambda_ps1, lambda_ps2
  REAL(kind=dp)                            ::  avnabpsi
END MODULE neo_output
