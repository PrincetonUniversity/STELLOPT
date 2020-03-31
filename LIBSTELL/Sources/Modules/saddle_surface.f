      MODULE saddle_surface
      USE stel_constants

      INTEGER, PARAMETER :: nsurf = 500
      INTEGER :: numsurf_sad, nopt_wsurf
      INTEGER, DIMENSION(nsurf) :: m_sad, n_sad
      INTEGER :: mpol_opt, ntor_opt, irho_bdy, irm0_bdy, izm0_bdy
      REAL(rprec), DIMENSION(nsurf) :: rmn_sad, zmn_sad
      INTEGER, DIMENSION(:), ALLOCATABLE :: nbrho_opt, mbrho_opt
      INTEGER, ALLOCATABLE :: nrz0_opt(:)
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: rbc,zbs,rhobc,delta_mn

      END MODULE saddle_surface
