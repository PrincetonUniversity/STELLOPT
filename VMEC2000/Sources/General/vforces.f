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
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::
     1    parmn, pazmn, pbrmn, pbzmn, pcrmn, pczmn, pblmn, pclmn
      REAL(rprec), POINTER, DIMENSION(:,:) ::
     1    parmn_e, parmn_o, pazmn_e, pazmn_o,
     2    pbrmn_e, pbrmn_o, pbzmn_e, pbzmn_o,
     3    pcrmn_e, pcrmn_o, pczmn_e, pczmn_o, 
     4    pblmn_e, pblmn_o, pclmn_e, pclmn_o
c-----------------------------------------------
      END MODULE vforces
