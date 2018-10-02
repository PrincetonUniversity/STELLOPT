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

      !VMEC2SFL
      USE vmec2sfl_vars_mod, vshat => shat
      USE vmec2sfl_mod
      
      !PTSM3D Files
      USE PTSM3D_setup
      USE PTSM3D_geom
      USE PTSM3D_itg
      USE PTSM3D_triplets
      USE PTSM3D_targets

       
      IMPLICIT NONE
      LOGICAL, INTENT(in)     :: lscreen
      INTEGER, INTENT(inout)  :: iflag

!----------------------------------------------------------------------
!     Local variables
!
!----------------------------------------------------------------------
      INTEGER :: maxPnt, nalpha0_, i, iunit, ier, ncnt, periods 
      INTEGER :: j, k, global_npol, nzgrid, nalpha, vmec_option
      REAL(rprec) :: dpdx, maxTheta, pval, pprime
      REAL(rprec) :: zeta_center, s_used
      REAL(rprec) :: periodszeta
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: th
      character(len=128) :: temp_str, gist_filename, num_str, vmec2sfl_geom_file
      LOGICAL :: uflag, res!, verbose
 
      REAL(rprec), PARAMETER :: zero   = 0.0_rprec
      REAL(rprec), PARAMETER :: one    = 1.0_rprec
      REAL(rprec), PARAMETER :: two    = 2.0_rprec

!----------------------------------------------------------------------
!     Replicate the GIST geometry calculations done in stellopt_txport
!     Do this for each (kx,ky) pair focusing only on a range defined
!     by theta_k and local_npol
!----------------------------------------------------------------------
      IF (lscreen) WRITE(6,'(a)') &
      &  ' -------------------------  BEGIN PTSM3D CALCULATION &
      & ------------------------ '

      nz0 = points_per_turn
      
      vmec_option = 0
      verbose = .false.

      call get_surface_quantities(s0,vmec_option)
      q0 = safety_factor_q
      shat = vshat
      maxTheta = abs(1.0/shat/dky) + 2*local_npol*pi
      global_npol = nint(maxTheta/pi)
      !make sure global_npol is even
      if (modulo(global_npol,2) == 1) global_npol = global_npol + 1
      maxTheta = global_npol*pi
      maxPnt = global_npol*nz0

      periods = float(global_npol)*nfp
      nzgrid = maxPnt/2
      periodszeta = periods*q0
      zeta_center = 0.0
      nalpha=1

      ALLOCATE(alpha(nalpha))
      ALLOCATE(zeta(-nzgrid:nzgrid))
      ALLOCATE(bmag(nalpha, -nzgrid:nzgrid))
      ALLOCATE(gradpar(nalpha, -nzgrid:nzgrid))
      ALLOCATE(gds2(nalpha, -nzgrid:nzgrid))
      ALLOCATE(gds21(nalpha, -nzgrid:nzgrid))
      ALLOCATE(gds22(nalpha, -nzgrid:nzgrid))
      ALLOCATE(gbdrift(nalpha, -nzgrid:nzgrid))
      ALLOCATE(gbdrift0(nalpha, -nzgrid:nzgrid))
      ALLOCATE(cvdrift(nalpha, -nzgrid:nzgrid))
      ALLOCATE(cvdrift0(nalpha, -nzgrid:nzgrid))
      ALLOCATE(jac_gist_inv(nalpha, -nzgrid:nzgrid))
      ALLOCATE(d_B_d_par(nalpha, -nzgrid:nzgrid))
      ALLOCATE(th(0:2*nzgrid))
      
      call vmec2sfl(nalpha,nzgrid,zeta_center,periodszeta)

      nz = 2*nzgrid+1
      li1 = 0
      li2 = 2*nzgrid
      CALL PTSM3D_initialize_geom
      if (allocated(geom)) deallocate(geom)
      allocate(geom(9,li1:li2))

      i = 1 !eventually we may allow for multiple s values
      geom(1,:) = gds22(i,:)/shat**2   ! g11
      geom(2,:) = gds21(i,:)/shat      ! g21
      geom(3,:) = gds2(i,:)            ! g22
      geom(4,:) = Bmag(i,:)            ! modB
      geom(5,:) = 1.0/abs(jac_gist_inv(i,:))  ! Jacobian
      geom(6,:) = Bmag(i,:)/2 * cvdrift(i,:)       ! L2
      geom(7,:) = -Bmag(i,:)/2/shat * cvdrift0(i,:) ! L1
      geom(8,:) = d_B_d_par(i,:)            ! dBdpar

      !zeta_center is 0 always, so assume theta_center is 0 too
      !Then conversion between theta and zeta is simply
      !theta = zeta/q
      th = zeta/q0
      geom(9,:) = th 

      temp_str = TRIM('gist_genet_'//TRIM(proc_string))
      ncnt = 0
      uflag = .false.
      DO WHILE (uflag .eqv. .false.)
        WRITE(num_str,"(I5.5)") ncnt 
        gist_filename = TRIM(temp_str)//&
        & '.' //TRIM(ADJUSTL(num_str))  
        INQUIRE(FILE=TRIM(gist_filename),EXIST=res)
        IF (res) THEN
          ncnt = ncnt + 1
        ELSE
          uflag = .true.
        END IF
      END DO 
      if (write_gist) then
        iunit = 50000+iflag+100*ncnt
        CALL safe_open(iunit,iflag,TRIM(gist_filename),&
          & 'unknown','formatted')
        WRITE(iunit,'(A)') '&PARAMETERS'
        WRITE(iunit,"(A,F12.7)") "s0 = ", normalized_toroidal_flux_used
        WRITE(iunit,"(A,F12.7)") "minor_a = ", L_reference
        WRITE(iunit,"(A,F12.7)") "Bref = ", B_reference
        WRITE(iunit,"(A,F12.7)") "q0 = ",ABS(q0)
        WRITE(iunit,"(A,F12.7)") "shat = ",shat 
        WRITE(iunit,"(A,I5)") "gridpoints = ",nzgrid*2+1
        WRITE(iunit,"(A,I7)") "n_pol = ",periods/nfp
        WRITE(iunit,"(A)") "/"
      end if

      DO j = li1,li2 
         !WRITE(iunit, "(9ES22.12E3)") gds2(i,j), gds21(i,j), &
         !    gds22(i,j), bmag(i,j), jac_gist_inv(i,j), &
         !    cvdrift(i,j), cvdrift0(i,j), gradpar(i,j), zeta(j)
         WRITE(iunit, "(9ES22.12E3)") geom(1,j), geom(2,j), geom(3,j), geom(4,j), &
                     geom(5,j), geom(6,j), geom(7,j), geom(8,j), geom(9,j)
      END DO
    
      IF (write_gist) CLOSE(iunit)

      CALL PTSM3D_set_norms 

      CALL PTSM3D_initialize_itg_solve
      ! Solve ITG dispersion relation for each (kx,ky)
      !CALL PTSM3D_itg_solve(j,k)
      CALL PTSM3D_itg_solve

      ! Call the rest of the PTSM3D functions
      CALL PTSM3D_initialize_triplets

      CALL PTSM3D_compute_triplets

      CALL PTSM3D_compute_targets
      
      IF (opt_target == 'zf') THEN
        ptsm3d_target = target_12f
        !IF (lscreen) WRITE(6,"(2A,F12.7)"),&
        WRITE(6,"(2A,F12.7)"),&
          & TRIM(gist_filename),&
          !& TRIM(TRIM(proc_string)//"."//TRIM(num_str)),&
          &", TARGET_12F  : ",target_12f 
      END IF
      IF (opt_target == 'nzf') THEN
        ptsm3d_target = target_qst
        !IF (lscreen) WRITE(6,"(A,F12.7)"),"TARGET_QST  : ",target_qst 
        WRITE(6,"(2A,F12.7)"),&
          & TRIM(gist_filename),&
         ! & TRIM(TRIM(proc_string)//"."//TRIM(num_str)),&
          & " TARGET_QST  : ",target_qst
      END IF
      IF (opt_target == 'combo') THEN
        ptsm3d_target = target_12f+target_qst 
        !IF (lscreen) WRITE(6,"(A,F12.7)"),"PTSM3D_TARGET   : ",ptsm3d_target
        WRITE(6,"(2A,F12.7)"),&
          !& TRIM(TRIM(proc_string)//"."//TRIM(num_str)),&
          & TRIM(gist_filename),&
          & "PTSM3D_TARGET   : ",ptsm3d_target
      END IF

      IF (lscreen) WRITE(6,'(a)') &
      &  ' -------------------------  END PTSM3D CALCULATION &
      & --------------------------'
      
      DEALLOCATE(alpha)
      DEALLOCATE(zeta)
      DEALLOCATE(bmag)
      DEALLOCATE(gradpar)
      DEALLOCATE(gds2)
      DEALLOCATE(gds21)
      DEALLOCATE(gds22)
      DEALLOCATE(gbdrift)
      DEALLOCATE(gbdrift0)
      DEALLOCATE(cvdrift)
      DEALLOCATE(cvdrift0)
      DEALLOCATE(jac_gist_inv)
      DEALLOCATE(d_B_d_par)
      DEALLOCATE(th)

      END SUBROUTINE stellopt_ptsm3d
