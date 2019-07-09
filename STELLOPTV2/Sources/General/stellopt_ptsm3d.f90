!----------------------------------------------------------------------
!     Subroutine:     stellopt_ptsm3d
!     Authors:        B.J. Faber (bfaber@wisc.edu)
!     Date:           7 Feb 2018
!     Description:    T
!
!----------------------------------------------------------------------
      SUBROUTINE stellopt_ptsm3d(lscreen,iflag)

!----------------------------------------------------------------------
!       Libraries
!----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_input_mod
      USE stellopt_vars
      USE stellopt_targets
      USE equil_vals
      USE equil_utils
      USE booz_persistent
      USE read_boozer_mod
      USE EZspline_obj
      USE EZspline

      !VMEC2PEST
      USE interfaces
      
      !PTSM3D Files
      !USE PTSM3D_setup
      !USE PTSM3D_geom
      !USE PTSM3D_itg
      !USE PTSM3D_triplets
      !USE PTSM3D_targets

       
      IMPLICIT NONE
      LOGICAL, INTENT(in)     :: lscreen
      INTEGER, INTENT(inout)  :: iflag

!----------------------------------------------------------------------
!     Local variables
!
!----------------------------------------------------------------------
!      INTEGER :: maxPnt, nalpha0_, i, iunit, ier, ncnt, periods 
!      INTEGER :: j, k, global_npol, nzgrid, nalpha, vmec_option
!      REAL(rprec) :: dpdx, maxTheta, pval, pprime
!      REAL(rprec) :: zeta_center, s_used
!      REAL(rprec) :: periodszeta
!      REAL(rprec), DIMENSION(:), ALLOCATABLE :: th
!      character(len=128) :: temp_str, gist_filename, num_str, vmec2sfl_geom_file
!      LOGICAL :: uflag, res!, verbose
! 
!      REAL(rprec), PARAMETER :: zero   = 0.0_rprec
!      REAL(rprec), PARAMETER :: one    = 1.0_rprec
!      REAL(rprec), PARAMETER :: two    = 2.0_rprec

!----------------------------------------------------------------------
!     Replicate the GIST geometry calculations done in stellopt_txport
!     Do this for each (kx,ky) pair focusing only on a range defined
!     by theta_k and local_npol
!----------------------------------------------------------------------
      real(rprec), dimension(:), allocatable :: surfaces, data_arr
      integer :: nx2, nx3, i, j
      character(len=128) :: x3_coord, norm_type, grid_type
      real(rprec) :: nfpi, x3_center
      character(len=16), dimension(7) :: ptsm3d_geom_strings
      IF (lscreen) WRITE(6,'(a)') &
      &  ' -------------------------  BEGIN PTSM3D CALCULATION &
      & ------------------------ '
      ptsm3d_geom_strings(1) = 'bmag'
      ptsm3d_geom_strings(2) = 'jac'
      ptsm3d_geom_strings(3) = 'g11'
      ptsm3d_geom_strings(4) = 'g12'
      ptsm3d_geom_strings(5) = 'g22'
      ptsm3d_geom_strings(6) = 'curv_drift_x1'
      ptsm3d_geom_strings(7) = 'curv_drift_x2'

      ! Move this to the PTSM3D namelist
      if(allocated(surfaces)) deallocate(surfaces)
      allocate(surfaces(1))
      surfaces(1) = 0.5
      nx2 = 1
      nx3 = 64
      allocate(data_arr(nx3+1))
      x3_center = 0.0
      nfpi = 5.0
      norm_type = "minor_r"
      grid_type = "gene"
      x3_coord = "theta"
      call vmec2pest_stellopt_interface(surfaces,nx2,nx3,x3_center,&
        &trim(x3_coord),nfpi,trim(norm_type),trim(grid_type))

      do j = 1,7
        call get_pest_data_interface(0,0,trim(ptsm3d_geom_strings(j)),1,nx3+1,data_arr)

      end do

      deallocate(surfaces,data_arr)
      
      END SUBROUTINE stellopt_ptsm3d
