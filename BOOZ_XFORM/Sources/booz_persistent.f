      MODULE booz_persistent
      USE stel_kinds
C     ONLY ARRAYS SMALL ENOUGH TO SAVE BETWEEN INTERNAL CALLS SHOULD BE STORED HERE
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: nu2_b, nu3_b
      REAL(rprec), DIMENSION(:), ALLOCATABLE ::
     1   sfull, scl, thgrd, ztgrd, xm, xn, xm_nyq, xn_nyq, xnb, xmb
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::
     1   cosm_b, cosn_b, sinm_b, sinn_b, rmnc, zmns, lmns,
     2   rmns, zmnc, lmnc
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::
     1   cosm_nyq, cosn_nyq, sinm_nyq, sinn_nyq
C-----------------------------------------------
      END MODULE booz_persistent
