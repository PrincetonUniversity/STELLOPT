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
      USE vmec2gs2_mod
      
      !PTSM3D Files
      USE PTSM3D_par
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
      REAL(rprec) :: a, s, Ba, Fa, iot, iotp, q, qprim
      REAL(rprec) :: dpdx, maxTheta, pval, pprime
      REAL(rprec) :: zeta_center, s_used
      REAL(rprec) :: L_ref, B_ref, periodszeta
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: g11,g12,g22,Bhat,abs_jac,L1,L2, &
                                                dBdt,zeta,alpha,th
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: bmag, gradpar, gds2, gds21, &
         gds22, gbdrift, gbdrift0, cvdrift, cvdrift0, jac_gist_inv, d_B_d_par
      character(len=128) :: temp_str, gist_filename, num_str
      LOGICAL :: uflag, res, verbose
 
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

      s = s0 
      Fa = phiedge
      a = Aminor
      Ba = ABS(Fa/(pi*a**2))
      CALL EZspline_interp(iota_spl,s,iot,ier)
      CALL EZspline_derivative(iota_spl,1,s,iotp,ier)
      !iot = iota_b(ik)
      !iotp = 0.5*(iota_b(ik+1)-iota_b(ik-1))/drho
      q = one/iot
      qprim = -iotp*q**2
      shat = two*s/q*qprim
      CALL get_equil_p(s,pval,ier,pprime)
      dpdx = -4.0_rprec*sqrt(s)/Ba**2 * pprime*mu0
      !maxPnt = nz0*local_npol
      
      q0 = q
      qprime = q0*shat/(2.0*s0) 
      Bref = Ba
      minor_a = a
      dtheta = pi2/real(nz0,rprec)
      maxTheta = nint(abs(1.0/qprim*1.0/dky))+local_npol*pi
      maxPnt = nint(2.0/dtheta*maxTheta)
      global_npol = nint(maxTheta/pi)
      
      
      periods = float(global_npol)
      nzgrid = maxPnt
      zeta_center = 0.0
      vmec_option = 0
      verbose = .false.
      nalpha = 1
      
      !debugging
      nzgrid = 6400
      periods = 400
      periodszeta = periods*q

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
      ALLOCATE(g11(0:2*nzgrid))
      ALLOCATE(g12(0:2*nzgrid))
      ALLOCATE(g22(0:2*nzgrid))
      ALLOCATE(Bhat(0:2*nzgrid))
      ALLOCATE(abs_jac(0:2*nzgrid))
      ALLOCATE(L1(0:2*nzgrid))
      ALLOCATE(L2(0:2*nzgrid))
      ALLOCATE(dBdt(0:2*nzgrid))
      ALLOCATE(th(0:2*nzgrid))

      CALL PTSM3D_initialize_geom(2*nzgrid+1)

      CALL PTSM3D_set_norms 

      CALL PTSM3D_initialize_itg_solve



      !Call Matt's program to calculate the quantities on a gs2 grid
      CALL vmec2gs2('wout_'//TRIM(PROC_STRING)//'.nc', nalpha, nzgrid,&
             zeta_center, periodszeta, s, vmec_option, verbose, s_used, q, shat, &
             L_ref, B_ref, alpha, zeta, bmag, gradpar, gds2, gds21, gds22, & 
             gbdrift, gbdrift0, cvdrift, cvdrift0, jac_gist_inv, d_B_d_par)


      ! Convert to GIST-GENE coordinates

      i = 1 !eventually we may allow for multiple s values
      g11 = gds22(i,:)/shat**2
      g12 = gds21(i,:)/shat
      g22 = gds2(i,:)
      Bhat = Bmag(i,:)
      abs_jac = 1.0/abs(jac_gist_inv(i,:))
      L2 = Bhat/2 * cvdrift(i,:)
      L1 = -Bhat/2/shat * cvdrift0(i,:)
      dBdt = d_B_d_par(i,:)

      !Set the parameters in PTSM3D
      gxx = g11
      gxy = g12
      gyy = g22
      modB = Bhat
      jac = abs_jac
      dBdy = L2
      dBdx = L1
      


      !zeta_center is 0 always, so assume theta_center is 0 too
      !Then conversion between theta and zeta is simply
      !theta = zeta/q
      th = zeta/q



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
        WRITE(iunit,"(A,F12.7)") "s0 = ",s_used
        WRITE(iunit,"(A,F12.7)") "minor_a = ", L_ref
        WRITE(iunit,"(A,F12.7)") "Bref = ", B_ref
        WRITE(iunit,"(A,F12.7)") "q0 = ",ABS(q)
        WRITE(iunit,"(A,F12.7)") "shat = ",shat 
        WRITE(iunit,"(A,I5)") "gridpoints = ",nzgrid*2+1
        WRITE(iunit,"(A,I7)") "n_pol = ",periods/nfp
        WRITE(iunit,"(A)") "/"
      end if

      DO j = 0,2*nzgrid
         !WRITE(iunit, "(9ES22.12E3)") gds2(i,j), gds21(i,j), &
         !    gds22(i,j), bmag(i,j), jac_gist_inv(i,j), &
         !    cvdrift(i,j), cvdrift0(i,j), gradpar(i,j), zeta(j)
         WRITE(iunit, "(9ES22.12E3)") g11(j), g12(j), g22(j), Bhat(j), &
                     abs_jac(j), L2(j), L1(j), dBdt(j), th(j)
      END DO
    
      IF (write_gist) CLOSE(iunit)

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
          & TRIM(TRIM(proc_string)//"."//TRIM(num_str)),&
          &", TARGET_12F  : ",target_12f 
      END IF
      IF (opt_target == 'nzf') THEN
        ptsm3d_target = target_qst
        !IF (lscreen) WRITE(6,"(A,F12.7)"),"TARGET_QST  : ",target_qst 
        WRITE(6,"(2A,F12.7)"),&
          & TRIM(TRIM(proc_string)//"."//TRIM(num_str)),&
          & "TARGET_QST  : ",target_qst
      END IF
      IF (opt_target == 'combo') THEN
        ptsm3d_target = target_12f+target_qst 
        !IF (lscreen) WRITE(6,"(A,F12.7)"),"PTSM3D_TARGET   : ",ptsm3d_target
        WRITE(6,"(2A,F12.7)"),&
          & TRIM(TRIM(proc_string)//"."//TRIM(num_str)),&
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
      DEALLOCATE(g11)
      DEALLOCATE(g12)
      DEALLOCATE(g22)
      DEALLOCATE(Bhat)
      DEALLOCATE(abs_jac)
      DEALLOCATE(L1)
      DEALLOCATE(L2)
      DEALLOCATE(dBdt)
      DEALLOCATE(th)


      END SUBROUTINE stellopt_ptsm3d
