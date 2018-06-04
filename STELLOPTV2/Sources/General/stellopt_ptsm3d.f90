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
      INTEGER :: maxPnt, nalpha0_, ialpha, i, iunit, ik, ier, ncnt 
      INTEGER :: j, k, global_npol, nzgrid, nalpha, vmec_option
      REAL(rprec) :: a, s, Ba, Fa, iot, iotp, qprim, pval, pprime
      REAL(rprec) :: dalpha, alpha0_start_, phi0, th, jac1, c
      REAL(rprec) :: gaa, gat, gst, gsa, gss, thetastar
      REAL(rprec) :: dBds, dBda, q, dpdx, sqrtg, absb, u, v, R, Z
      REAL(rprec) :: temp1, temp2, temp3, abserr, alpha0_end
      REAL(rprec) :: alpha0_start, maxTheta, zeta_center, s_used
      REAL(rprec) :: L_ref, B_ref, periods
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: g11,g12,g22,Bhat,abs_jac,L1,L2, &
                                                dBdt,zeta,alpha
      REAL(rprec), DIMENSION(3) :: sflCrd0,sflCrd, sflCrd_sav, gradS,gradThetaStar,&
                                   gradPhi,mag,gradAlpha, wrk, gradB,R_grad,Z_grad,&
                                   esubs, esubu, esubv, es, eu, ev, gradlam, ea, et
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
      ALLOCATE(g11(-nzgrid:nzgrid))
      ALLOCATE(g12(-nzgrid:nzgrid))
      ALLOCATE(g22(-nzgrid:nzgrid))
      ALLOCATE(Bhat(-nzgrid:nzgrid))
      ALLOCATE(abs_jac(-nzgrid:nzgrid))
      ALLOCATE(L1(-nzgrid:nzgrid))
      ALLOCATE(L2(-nzgrid:nzgrid))
      ALLOCATE(dBdt(-nzgrid:nzgrid))




      CALL vmec2gs2('wout_'//TRIM(PROC_STRING)//'.nc', nalpha, nzgrid,&
             zeta_center, periods, s, vmec_option, verbose, s_used, q, shat, &
             L_ref, B_ref, alpha, zeta, bmag, gradpar, gds2, gds21, gds22, & 
             gbdrift, gbdrift0, cvdrift, cvdrift0, jac_gist_inv, d_B_d_par)

      CALL PTSM3D_initialize_geom(nzgrid)

      CALL PTSM3D_set_norms 

      CALL PTSM3D_initialize_itg_solve

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
        WRITE(iunit,"(A,I5)") "gridpoints = ",nzgrid
        WRITE(iunit,"(A,F12.7)") "n_pol = ",periods
        WRITE(iunit,"(A)") "/"
      end if

      CALL PTSM3D_initialize_geom(maxPnt)

      CALL PTSM3D_set_norms 

      CALL PTSM3D_initialize_itg_solve

      
!      sflCrd0(1) = s
!      sflCrd0(2) = 0.0 
!      phi0 = 0.0
!      sflCrd0(3) = phi0
!      DO i = 1, maxPnt ! Loop over field line
!        th = -maxTheta + (i-1)*dtheta
!        sflCrd(1) = sflCrd0(1)
!        sflCrd(2) = th
!        sflCrd(3) = sflCrd0(3) + q*(th-sflCrd0(2))
!        ! Get Metric Elements
!        sflCrd_sav = sflCrd
!        CALL pest2vmec(sflCrd) ! Returns on 2pi grid
!        u = sflCrd(2)
!        v = sflCrd(3)
!        IF (u < 0) THEN
!           u = -MOD(ABS(u),pi2)
!           u = u + pi2
!        END IF
!        IF (v < 0) THEN
!           v = -MOD(ABS(v),pi2)
!           v = v + pi2
!        END IF
!        IF (u > pi2) u = MOD(u,pi2)
!        IF (v > pi2) v = MOD(v,pi2)
!        CALL get_equil_RZ(sflCrd(1),u,v,R,Z,ier,R_GRAD=R_grad,Z_GRAD=Z_grad)
!        ! Get equil_RZ returns dR/drho and dZ/drho
!        !   dR/ds=(0.5/rho)*dR/drho=(0.5/sqrt(s))*dR/drho
!        R_grad(3) = 0.5*R_grad(3)/SQRT(sflCrd(1))
!        Z_grad(3) = 0.5*Z_grad(3)/SQRT(sflCrd(1))
!        ! e_s
!        esubs(1) = R_grad(3)
!        esubs(2) = zero
!        esubs(3) = Z_grad(3)
!        ! e_u
!        esubu(1) = R_grad(1)
!        esubu(2) = zero
!        esubu(3) = Z_grad(1)
!        ! e_v
!        esubv(1) = R_grad(2)
!        esubv(2) = one
!        esubv(3) = Z_grad(2) 
!        esubv(1) = esubv(1)*nfp
!        esubv(3) = esubv(3)*nfp
!        ! Cylindrical Coordianates
!        !CALL EZspline_interp(R_spl,u,v,sflCrd(1),R,iflag)
!        esubs(2) = esubs(2)*R
!        esubu(2) = esubu(2)*R
!        esubv(2) = esubv(2)*R
!        ! sqrt(g) = R*(Ru*Zs-Rs*Zu)
!        sqrtg = R*(R_grad(1)*Z_grad(3)-R_grad(3)*Z_grad(1))
!        ! e^i = (e_j x e_k)/sqrt(g)
!        !CALL EZspline_interp(G_spl,u,v,sflCrd(1),sqrtg,iflag)
!        CALL cross_product(esubu,esubv,es)
!        CALL cross_product(esubv,esubs,eu)
!        CALL cross_product(esubs,esubu,ev)
!        es = es/sqrtg
!        eu = eu/sqrtg
!        ev = ev/sqrtg
!        ! Get Field (before we adjust eu so gradB is correct)
!        CALL get_equil_Bflx(sflCrd(1),u,v,temp1,temp2,temp3,ier,absb,gradB)
!        gradB = gradB(3)*es + gradB(1)*eu + gradB(2)*ev*nfp
!        ! Now Adjust e^u for lambda
!        !IF (pest) THEN
!           CALL get_equil_L(sflCrd(1),u,v,temp1,ier,gradlam)
!           eu = eu + gradlam(3)*es + gradlam(1)*eu + gradlam(2)*ev*nfp
!        !ELSE
!        !END IF
!        ! Now do some calculations
!        sflCrd = sflCrd_sav
!        thetastar = sflCrd(2)
!        gradS = es
!        gradThetaStar = eu
!        gradPhi = ev
!        gradAlpha = qprim * thetaStar*gradS+q*gradThetaStar-gradPhi
!        alpha     = q*thetaStar - sflCrd(3)
!        ! Metrice and Jacobian
!        gss = DOT_PRODUCT(gradS,gradS)
!        gsa = DOT_PRODUCT(gradS,gradAlpha)
!        gst = DOT_PRODUCT(gradS,gradThetaStar)
!        gaa = DOT_PRODUCT(gradAlpha,gradAlpha)
!        gat = DOT_PRODUCT(gradAlpha,gradThetaStar)
!        CALL cross_product(gradS,gradAlpha,wrk)
!        jac1 = one/DOT_PRODUCT(wrk,gradThetaStar)
!        
!        Bhat = absb/Ba
!        modB(i-1) = Bhat
!
!        g11 = gss*a**2/(4*s)
!        !g12(ialpha,i) = gsa*a**2*iot/2*sloc_fac
!        g12 = gsa*a**2*iot/2
!        g22 = (Bhat**2+g12**2)/g11
!        abs_jac = ABS(jac1*2*q/a**3)
!        gxx(i-1) = g11
!        gxy(i-1) = g12
!        gyy(i-1) = g22
!        jac(i-1) = minor_a/(Bref*q0)*abs_jac 
!                 
!        CALL cross_product(gradAlpha,gradThetaStar, es) 
!        CALL cross_product(gradThetaStar,gradS,     ea) 
!        CALL cross_product(gradS,gradAlpha,         et)  
!        es(:) = es(:)*jac1
!        ea(:) = ea(:)*jac1
!        et(:) = et(:)*jac1
!        gradB = gradB/Ba
!        dBds = DOT_PRODUCT(gradB,es)
!        dBda = DOT_PRODUCT(gradB,ea)
!        dBdt = DOT_PRODUCT(gradB,et)
!        
!        c = iot*iot*a**4
!        L1 = q/sqrt(s)*(dBda + c*(gss*gat-gsa*gst)*dBdt/(4*Bhat**2))
!        L2 = two*sqrt(s)*(dBds + c*(gaa*gst-gsa*gat)*dBdt/(4*Bhat**2))
!        dBdx(i-1) = L1
!        dBdy(i-1) = L2
!        IF (write_gist) THEN
!          WRITE(iunit,"(9ES20.10)") g11,g12,g22,Bhat,&
!          & abs_jac, L2, L1, th, 0.0; CALL FLUSH(iunit)
!        END IF
!        
!      ENDDO ! End loop over field line
      
      i = 1 !eventually we may allow for multiple s values
      g11 = gds22(i,:)/shat**2
      g12 = gds21(i,:)/shat
      g22 = gds2(i,:)
      Bhat = Bmag(i,:)
      abs_jac = 1.0/abs(jac_gist_inv(i,:))
      L2 = Bhat/2 * cvdrift(i,:)
      L1 = Bhat/2/shat**2 * cvdrift0(i,:)
      dBdt = d_B_d_par(i,:)
     
      DO j = -nzgrid,nzgrid
         !WRITE(iunit, "(9ES22.12E3)") gds2(i,j), gds21(i,j), &
         !    gds22(i,j), bmag(i,j), jac_gist_inv(i,j), &
         !    cvdrift(i,j), cvdrift0(i,j), gradpar(i,j), zeta(j)
         WRITE(iunit, "(9ES22.12E3)") g11(j), g12(j), g22(j), Bhat(j), &
                     abs_jac(j), L2(j), L1(j), dBdt(j), zeta(j)
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

      CALL PTSM3D_finalize_triplets

      CALL PTSM3D_finalize_itg_solve

      CALL PTSM3D_finalize_geom

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


      END SUBROUTINE stellopt_ptsm3d
