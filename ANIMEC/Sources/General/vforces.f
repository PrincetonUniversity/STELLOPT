      MODULE vforces
      USE stel_kinds, ONLY: rprec
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), DIMENSION(:), ALLOCATABLE, TARGET ::
     1    armn, azmn, brmn, bzmn, crmn, czmn, blmn, clmn
      REAL(rprec), POINTER, DIMENSION(:) ::
     1    armn_e, armn_o, azmn_e, azmn_o,
     2    brmn_e, brmn_o, bzmn_e, bzmn_o,
     3    crmn_e, crmn_o, czmn_e, czmn_o, blmn_e,
     4    blmn_o, clmn_e, clmn_o
c-----------------------------------------------
      END MODULE vforces
