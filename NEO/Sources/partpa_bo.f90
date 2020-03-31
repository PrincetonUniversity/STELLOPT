
MODULE partpa_bo
  USE neo_precision
! Exchange between flint_bo and rhs_bo
  USE sizey_bo
  INTEGER                                  :: ipmax
  INTEGER,       DIMENSION(:), ALLOCATABLE :: isw,ipa,icount
  REAL(kind=dp)                            :: pard0,bmod0
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: eta
END MODULE partpa_bo
