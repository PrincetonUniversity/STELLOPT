      MODULE vacmod
      USE vacmod0
      USE vac_persistent
      USE vmec_input, ONLY: lasym
      USE vmec_params, ONLY: signgs
      USE vparams, ONLY: zero, one, c2p0, cp5
      USE mgrid_mod, ONLY: nr0b, np0b, nz0b, 
     1   rminb, zminb, rmaxb, zmaxb, delrb, delzb
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: p5 = cp5, two = c2p0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: nfper, nvper
      REAL(rprec), DIMENSION(:), ALLOCATABLE, TARGET :: potvac
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: bvecsav, amatsav,
     1   bexni, brv, bphiv, bzv, bsqvac, bsqvac0, r1b, rub, rvb, z1b,
     2   zub, zvb, bexu, bexv, bexn, auu, auv, avv, snr, snv, snz, drv,
     3   guu_b, guv_b, gvv_b, rzb2, rcosuv, rsinuv,
     5   bredge, bpedge, bzedge
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: raxis_nestor, 
     1                                          zaxis_nestor
      REAL(rprec) :: bsubvvac, pi2,
     2   pi3, pi4, alp, alu, alv, alvp, onp, onp2
      END MODULE vacmod
