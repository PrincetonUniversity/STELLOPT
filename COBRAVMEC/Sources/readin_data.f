      MODULE readin_data
      USE stel_kinds
      USE general_dimensions
      INTEGER,DIMENSION(:),ALLOCATABLE :: list
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: hiota, hphip,
     1   hpres
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: xn_v, xm_v
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: xn_vnyq, xm_vnyq
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: lmnsh, bmnch,
     1   bsupvmnch, bsupumnch
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: lmnch, bmnsh,               ! RS110909 - For ASYMMETRIC input
     1   bsupvmnsh, bsupumnsh
      LOGICAL :: lscreen
      END MODULE readin_data
