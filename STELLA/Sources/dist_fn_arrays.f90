!> A container for the arrays that are used to store the distribution function among other things.
!!  These need to be accessible at a lower dependency level than the dist_fn module itself.
!! These arrays are allocated in the function dist_fn::allocate_arrays. 

module dist_fn_arrays

  public :: gnew, gold
  public :: g0, g1, g2, g3
  public :: g_krook, g_proj
  public :: gvmu
  public :: kperp2, dkperp2dr
  public :: wstar, wstarp
  public :: wdriftx_g, wdrifty_g
  public :: wdriftx_phi, wdrifty_phi
  public :: wdriftpx_g, wdriftpy_g
  public :: wdriftpx_phi, wdriftpy_phi

  ! dist fn
  complex, dimension (:,:,:,:,:), allocatable :: gnew, gold
  ! (naky, nakx, -nzgrid:nzgrid, ntubes, -vmu-layout-)

  complex, dimension (:,:,:,:,:), allocatable :: g0, g1, g2, g3
  ! (naky, nakx, -nzgrid:nzgrid, ntubes, -vmu-layout-)

  ! needed to implement time-delayed source when using Krook operator
  complex, dimension (:,:,:,:), allocatable :: g_krook
  ! (nakx, -nzgrid:nzgrid, ntubes, -vmu-layout-)

  ! needed to implement time-delayed source when using projection method
  complex, dimension (:,:,:,:), allocatable :: g_proj
  ! (nakx, -nzgrid:nzgrid, ntubes, -vmu-layout-)

  complex, dimension (:,:,:), allocatable :: gvmu
  ! (nvpa, nmu, -kxkyz-layout-)

  real, dimension (:,:,:), allocatable :: wstar, wstarp
  ! (nalpha, -nzgrid:nzgrid, -vmu-layout-)

  real, dimension (:,:,:), allocatable :: wdriftx_g, wdrifty_g
  real, dimension (:,:,:), allocatable :: wdriftx_phi, wdrifty_phi

  real, dimension (:,:,:), allocatable :: wdriftpx_g, wdriftpy_g
  real, dimension (:,:,:), allocatable :: wdriftpx_phi, wdriftpy_phi
  ! (nalpha, -nzgrid:nzgrid, -vmu-layout-)

  real, dimension (:,:,:,:), allocatable :: kperp2, dkperp2dr
  ! (naky, nakx, nalpha, -nzgrid:nzgrid)
  ! note: dkperp2dr is divided by kperp2

end module dist_fn_arrays
