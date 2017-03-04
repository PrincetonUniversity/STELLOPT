!=====================================================
SUBROUTINE spfper(np1,amx1,amx2,amx3)

! Helper routine for splfi

  USE neo_precision

  IMPLICIT NONE

  INTEGER,                       INTENT(in)  :: np1
  REAL(kind=dp), DIMENSION(np1), INTENT(out) :: amx1, amx2, amx3
  REAL(kind=dp)                              :: beta, ss
  INTEGER                                    :: n, n1, i, i1

  n = np1-1

  n1 = n-1
  amx1(1) = 2
  amx2(1) = 0.5_dp
  amx3(1) = 0.5_dp
  amx1(2) = sqrt(15._dp)/2
  amx2(2) = ONE/amx1(2)
  amx3(2) = -.25_dp/amx1(2)
  beta = 3.75_dp
  DO i = 3,n1
     i1 = i-1
     beta = 4-ONE/beta
     amx1(i) = sqrt(beta)
     amx2(i) = ONE/amx1(i)
     amx3(i) = -amx3(i1)/amx1(i)/amx1(i1)
  END DO
  amx3(n1) = amx3(n1)+ONE/amx1(n1)
  amx2(n1) = amx3(n1)
  ss = 0
  DO i = 1,n1
     ss = ss+amx3(i)*amx3(i)
  END DO
  amx1(n) = sqrt(4-ss)

  RETURN
END SUBROUTINE spfper
