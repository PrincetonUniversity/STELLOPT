      MODULE booz_params
      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: unit_booz = 20
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: jsize
      INTEGER, ALLOCATABLE, DIMENSION(:) :: jlist
      INTEGER, DIMENSION(0:1) :: ntorsum
      INTEGER :: nfp, ns, mpol, mpol1, ntor, mnmax
      INTEGER :: mpol_nyq, ntor_nyq, mnmax_nyq
      INTEGER :: nu_boz, nv_boz, nunv, mboz, nboz, mnboz
      REAL(rprec), DIMENSION(:), ALLOCATABLE ::
     1  hiota, phip, pres, beta_vol, phi, buco, bvco, gpsi, ipsi
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: pmns, pmnc
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: bmncb,
     1  rmncb, zmnsb, pmnsb, gmncb, bmod_b
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: bmnsb,
     1  rmnsb, zmncb, pmncb, gmnsb
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: bsubumnc, 
     1  bsubvmnc, bsubumns, bsubvmns, bmodmnc, bmodmns 
      REAL(rprec) :: ohs
      LOGICAL :: lscreen, lasym_b
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: lsurf_boz
      END MODULE booz_params
