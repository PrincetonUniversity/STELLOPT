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
      USE read_wout_mod

      !VMECTOOLS
      USE interfaces
      
      !PTSM3D Files
      USE ptsm3d_setup
      USE ptsm3d_geom
      USE ptsm3d_itg
      USE ptsm3d_triplets
      USE ptsm3d_targets

       
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
      !real(rprec), parameter :: pi = 3.1415926535897932846264338327950d+0 
      real(rprec), dimension(:), allocatable :: surfaces, data_arr
      integer :: nx2, nx3, i, j
      character(len=128) :: x3_coord, norm_type, grid_type
      real(rprec) :: nfpi, x3_center, max_z
      real(rprec), dimension(1) :: s0_data, vmec_data
      character(len=16), dimension(8) :: ptsm3d_geom_strings
      IF (lscreen) WRITE(6,'(a)') &
      &  ' -------------------------  BEGIN PTSM3D CALCULATION &
      & ------------------------ '
      ptsm3d_geom_strings(4) = 'bmag'
      ptsm3d_geom_strings(5) = 'jac'
      ptsm3d_geom_strings(1) = 'g11'
      ptsm3d_geom_strings(2) = 'g12'
      ptsm3d_geom_strings(3) = 'g22'
      ptsm3d_geom_strings(6) = 'curv_drift_x2'
      ptsm3d_geom_strings(7) = 'curv_drift_x1'
      ptsm3d_geom_strings(8) = 'x3'

      call vmec2pest_stellopt_interface(&
        & (/s0/),1,1,0d+0,"zeta",1d+0,"minor_r","gene")

      call get_pest_data_interface(0,0,"shat",0,1,s0_data)
      shat = s0_data(1)
      ! Move this to the PTSM3D namelist
      if(allocated(surfaces)) deallocate(surfaces)
      allocate(surfaces(1))
      surfaces(1) = s0
      nx2 = 1
      max_z = 1.0/(abs(shat)*dky)
      call get_pest_data_interface(0,0,"nfp",-1,1,vmec_data)
      nfpi = (real(ceiling(max_z/pi)) + local_npol)*vmec_data(1)

      nx3 = (ceiling(max_z/pi)+local_npol)*points_per_turn
      nz = nx3+1
      li1=0
      li2=nx3
      allocate(data_arr(nx3+1))
      x3_center = 0.0
      norm_type = "minor_r"
      grid_type = "gene"
      x3_coord = "theta"

      allocate(geom(8,nx3+1))
      call ptsm3d_initialize_geom
      call vmec2pest_stellopt_interface(surfaces,nx2,nx3,x3_center,&
        &trim(x3_coord),nfpi,trim(norm_type),trim(grid_type))

      do j = 1,8
        call get_pest_data_interface(0,0,trim(ptsm3d_geom_strings(j)),1,nx3+1,data_arr)
        geom(j,:) = data_arr
      end do

      call ptsm3d_set_norms
      call ptsm3d_initialize_itg_solve
      call ptsm3d_itg_solve
      !print *, omdw(:,:,0)
      call ptsm3d_initialize_triplets
      call ptsm3d_compute_triplets
      call ptsm3d_compute_targets
      if (opt_target .eq. 'zf') then
        ptsm3d_target = target_12f
      elseif (opt_target .eq. 'nzf') then
        ptsm3d_target = target_qst
      endif
      if (lscreen) write(6,"((A),2x,F12.7)") "PTSM3D target: ",ptsm3d_target

      call ptsm3d_finalize_triplets
      call ptsm3d_finalize_itg_solve
      call ptsm3d_finalize_geom

      !ptsm3d_finalize_geom deallocates geom
      deallocate(surfaces,data_arr)
      
      END SUBROUTINE stellopt_ptsm3d
