      SUBROUTINE printmatrix(amatp, amati, data_in, 
     1                       i, i1, ip, is, type_in)
      USE vparams, ONLY: rprec, nthreed
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: i, i1
      REAL(rprec)  :: data_in
      CHARACTER(LEN=10) :: type_in
      INTEGER, DIMENSION(*) :: ip, is
      REAL(rprec), DIMENSION(*) :: amatp, amati
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: k
C-----------------------------------------------

      WRITE (nthreed, 10) i, type_in, i1, data_in, 
     1    (amatp(ip(k)),k=1,5), (amati(is(k)),k=1,5)

   10 FORMAT(1x,i3,a10,i2,')',1p,11e10.2)

      END SUBROUTINE printmatrix
