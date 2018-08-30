      MODULE vac_persistent
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, DIMENSION(:), ALLOCATABLE :: imirr
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: sinper, cosper,
     1   sinuv, cosuv, tanu, tanv, xmpot, xnpot, csign
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: sinu, cosu,
     1 sinv, cosv, sinui, cosui, sinu1, cosu1, sinv1, cosv1
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: cmns
      END MODULE vac_persistent
