
MODULE partpa_cur
  USE neo_precision
! Exchange between flint_cur and rhs_cur
  USE sizey_cur
  REAL(kind=dp)                            :: bmod0
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: y_part, yfac, sqyfac
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: k_fac1, k_fac2
END MODULE partpa_cur
