
MODULE neo_spline
! Spline arrays
  USE neo_precision
  INTEGER, PARAMETER       ::   mt = 1
  INTEGER, PARAMETER       ::   mp = 1
  INTEGER                  ::   theta_ind, phi_ind
  INTEGER                  ::   ierr

  REAL(kind=dp) ::  theta_d, phi_d

! Spline array for modb
  REAL(kind=dp),    DIMENSION(:,:,:,:), ALLOCATABLE :: b_spl
! Spline array for geodesic curviture
  REAL(kind=dp),    DIMENSION(:,:,:,:), ALLOCATABLE :: k_spl
! Spline array for sqrg11
  REAL(kind=dp),    DIMENSION(:,:,:,:), ALLOCATABLE :: g_spl
! Spline array for parallel derivative
  REAL(kind=dp),    DIMENSION(:,:,:,:), ALLOCATABLE :: p_spl
! Spline array for quasi-toroidal phi component of b
  REAL(kind=dp),    DIMENSION(:,:,:,:), ALLOCATABLE :: q_spl
END MODULE neo_spline
