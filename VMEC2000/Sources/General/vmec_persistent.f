      MODULE vmec_persistent
      USE stel_kinds, ONLY: rprec
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, DIMENSION(:), ALLOCATABLE :: ixm, jmin3
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: cosmu, sinmu,
     1   cosmum, sinmum, cosmumi, sinmumi,
     2   cosnv, sinnv, cosnvn, sinnvn
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: cosmui, sinmui, 
     1             cosmui3, cosmumi3
      REAL(rprec), DIMENSION(:), ALLOCATABLE, TARGET ::
     1   xm, xn, xm_nyq, xn_nyq
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: cos01, sin01
c-----------------------------------------------
      END MODULE vmec_persistent
