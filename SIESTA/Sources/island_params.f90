      MODULE island_params
      USE stel_kinds
!
!     CONTAINS VARIABLE DECLARATIONS FOR VMECPP (ISLAND) CODE
!     THAT WILL BE SHARED THROUGHOUT THE PROJECT
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
	  INTEGER :: ns_i, nu_i, nv_i, nuv_i                                 ! size of s, u, v dimensions of (island) metric tensors
      INTEGER :: mpol_i, ntor_i
!      INTEGER :: mpol32_i, ntor32_i                                      ! extended grid dims
      INTEGER :: nfp_i, mnmax_i
      INTEGER :: nsh                                                     ! Number of points in half mesh  
      REAL(rprec) :: hs_i, ohs_i, dnorm_i, gnorm_i, wb_i, wp_i
      REAL(rprec), PARAMETER :: gamma = 5._dp/3._dp                      ! Adiabatic constant
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: cosmu, cosmum,        &
     &    cosmui, sinmu, sinmum, sinmui, cosnv, cosnvn, sinnv, sinnvn
      REAL(rprec), ALLOCATABLE, DIMENSION(:)   :: phipf_i, chipf_i, presf_i, vp_f
!-----------------------------------------------

      END MODULE island_params
