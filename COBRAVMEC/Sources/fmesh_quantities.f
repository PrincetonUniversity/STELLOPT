      MODULE fmesh_quantities
      USE stel_kinds
      USE general_dimensions
      USE readin_data, ONLY: xn_v, xm_v, xm_vnyq, xn_vnyq
      IMPLICIT NONE
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: iotaf, phipf,
     1   presf, mercierf, iotapf, prespf, radios
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: rmncf, zmnsf,
     1   lmnsf, bmncf, bsupvmncf, bsupumncf, rmncpf, zmnspf,
     2   lmnspf, bmncpf
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: rmnsf, zmncf,           ! 110909 RS = for ASYMMETRIC input
     1   lmncf, bmnsf, bsupvmnsf, bsupumnsf, rmnspf, zmncpf,
     2   lmncpf, bmnspf
      END MODULE fmesh_quantities
