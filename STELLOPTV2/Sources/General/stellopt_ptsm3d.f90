!----------------------------------------------------------------------
!     Subroutine:     stellopt_ptsm3d
!     Authors:        B.J. Faber (bfaber@wisc.edu)
!     Date:           7 Feb 2018
!     Description:    T
!
!----------------------------------------------------------------------
      SUBROUTINE stellopt_ptsm3d

!----------------------------------------------------------------------
!       Libraries
!----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellop_input_mod
      USE stellopt_vars
      USE equil_vals
      USE equil_utils
      
      !PTSM3D Files
      USE PTSM3D_par
      USE PTSM3D_geom
      USE PTSM3D_itg
      USE PTSM3D_triplets
      USE PTSM3D_target

       
      IMPLICIT NONE
      LOGICAL, INTENT(in)     :: lscreen
      LOGICAL, INTENT(inout)  :: iflag

!----------------------------------------------------------------------
!     Local variables
!
!----------------------------------------------------------------------
      INTEGER :: maxPnt, nalpha0_, ialpha, i, iunit, ik, ier 
      REAL(rprec) :: a, s, Ba, Fa, iot,iotp,qprim,&
                     pval, pprime, dalpha, alpha0_start_, phi0, &
                     th, jac1, c, &
                     gaa, gat, gst, gsa, gss, alpha, thetastar, &
                     dBds, dBda, q, dpdx, sqrtg, modb, &
                     u, v, R, Z, &
                     temp1,temp2,temp3, abserr alpha0_end, alpha0_start
      REAL(rprec), DIMENSION(3) :: sflCrd0,sflCrd, sflCrd_sav, gradS,gradThetaStar,&
                                   gradPhi,mag,gradAlpha, wrk, gradB,R_grad,Z_grad,&
                                   esubs, esubu, esubv, es, eu, ev, gradlam, ea, et
 
      REAL(rprec), PARAMETER :: zero   = 0.0_rprec
      REAL(rprec), PARAMETER :: one    = 1.0_rprec
      REAL(rprec), PARAMETER :: two    = 2.0_rprec

!----------------------------------------------------------------------
!     Replicate the GIST geometry calculations done in stellopt_txport
!     Do this for each (kx,ky) pair focusing only on a range defined
!     by theta_k and local_npol
!----------------------------------------------------------------------

      CALL PTSM3D_read_parameters

      Fa = phiedge
      a = Aminor
      Ba = ABS(Fa/(pi*a**2))
      CALL EZspline_interp(iota_spl,s,iot,iflag)
      CALL EZspline_derivative(iota_spl,1,s,iotp,iflag)
      !iot = iota_b(ik)
      !iotp = 0.5*(iota_b(ik+1)-iota_b(ik-1))/drho
      q = one/iot
      qprim = -iotp*q**2
      shat = two*s/q*qprim
      CALL get_equil_p(s,pval,iflag,pprime)
      dpdx = -4.0_rprec*sqrt(s)/Ba**2 * pprime*mu0
      maxPnt = nz0*local_npol

      q0 = q
      qprime = q0*shat/(2.0*s0) 
      Bref = Ba
      minor_a = a
      s0 = s
      CALL PTSM3D_intialize_geom
      CALL PTSM3D_set_norms 
      CALL PTSM3D_initialize_itg_solve

      DO k=lk1+1,lk2
        DO j=lj1,lj2
          sflCrd0(1) = s
          sflCrd0(2) = 1.0/(2.0*s0*qprime)*kx(j)/ky(k)
          !dtheta = pi2*local_npol/maxPnt
          phi0 = 0.0
          sflCrd0(3) = phi0
          DO i = 1, maxPnt ! Loop over field line
            th = -pi*local_npol + (i-1)*dtheta
            sflCrd(1) = sflCrd0(1)
            sflCrd(2) = th
            sflCrd(3) = sflCrd0(3) + q*(th-sflCrd0(2))
            ! Get Metric Elements
            sflCrd_sav = sflCrd
            CALL pest2vmec(sflCrd) ! Returns on 2pi grid
            u = sflCrd(2)
            v = sflCrd(3)
            IF (u < 0) THEN
               u = -MOD(ABS(u),pi2)
               u = u + pi2
            END IF
            IF (v < 0) THEN
               v = -MOD(ABS(v),pi2)
               v = v + pi2
            END IF
            IF (u > pi2) u = MOD(u,pi2)
            IF (v > pi2) v = MOD(v,pi2)
            CALL get_equil_RZ(sflCrd(1),u,v,R,Z,iflag,R_GRAD=R_grad,Z_GRAD=Z_grad)
            ! Get equil_RZ returns dR/drho and dZ/drho
            !   dR/ds=(0.5/rho)*dR/drho=(0.5/sqrt(s))*dR/drho
            R_grad(3) = 0.5*R_grad(3)/SQRT(sflCrd(1))
            Z_grad(3) = 0.5*Z_grad(3)/SQRT(sflCrd(1))
            ! e_s
            esubs(1) = R_grad(3)
            esubs(2) = zero
            esubs(3) = Z_grad(3)
            ! e_u
            esubu(1) = R_grad(1)
            esubu(2) = zero
            esubu(3) = Z_grad(1)
            ! e_v
            esubv(1) = R_grad(2)
            esubv(2) = one
            esubv(3) = Z_grad(2) 
            esubv(1) = esubv(1)*nfp
            esubv(3) = esubv(3)*nfp
            ! Cylindrical Coordianates
            !CALL EZspline_interp(R_spl,u,v,sflCrd(1),R,iflag)
            esubs(2) = esubs(2)*R
            esubu(2) = esubu(2)*R
            esubv(2) = esubv(2)*R
            ! sqrt(g) = R*(Ru*Zs-Rs*Zu)
            sqrtg = R*(R_grad(1)*Z_grad(3)-R_grad(3)*Z_grad(1))
            ! e^i = (e_j x e_k)/sqrt(g)
            !CALL EZspline_interp(G_spl,u,v,sflCrd(1),sqrtg,iflag)
            CALL cross_product(esubu,esubv,es)
            CALL cross_product(esubv,esubs,eu)
            CALL cross_product(esubs,esubu,ev)
            es = es/sqrtg
            eu = eu/sqrtg
            ev = ev/sqrtg
            ! Get Field (before we adjust eu so gradB is correct)
            CALL get_equil_Bflx(sflCrd(1),u,v,temp1,temp2,temp3,iflag,modb,gradB)
            gradB = gradB(3)*es + gradB(1)*eu + gradB(2)*ev*nfp
            ! Now Adjust e^u for lambda
            !IF (pest) THEN
               CALL get_equil_L(sflCrd(1),u,v,temp1,iflag,gradlam)
               eu = eu + gradlam(3)*es + gradlam(1)*eu + gradlam(2)*ev*nfp
            !ELSE
            !END IF
            ! Now do some calculations
            sflCrd = sflCrd_sav
            thetastar = sflCrd(2)
            gradS = es
            gradThetaStar = eu
            gradPhi = ev
            gradAlpha = qprim * thetaStar*gradS+q*gradThetaStar-gradPhi
            alpha     = q*thetaStar - sflCrd(3)
            ! Metrice and Jacobian
            gss = DOT_PRODUCT(gradS,gradS)
            gsa = DOT_PRODUCT(gradS,gradAlpha)
            gst = DOT_PRODUCT(gradS,gradThetaStar)
            gaa = DOT_PRODUCT(gradAlpha,gradAlpha)
            gat = DOT_PRODUCT(gradAlpha,gradThetaStar)
            CALL cross_product(gradS,gradAlpha,wrk)
            jac1 = one/DOT_PRODUCT(wrk,gradThetaStar)
            
            Bhat(i) = modb/Ba
            modB(i-1) = Bhat
 
            g11(i) = gss*a**2/(4*s)
            !g12(ialpha,i) = gsa*a**2*iot/2*sloc_fac
            g12(i) = gsa*a**2*iot/2
            g22(i) = (Bhat(i)**2+g12(i)**2)/g11(i)
            abs_jac(i) = ABS(jac1*2*q/a**3)
            gxx(i-1) = g11(i)
            gxy(i-1) = g12(i)
            gyy(i-1) = g22(i)
            jac(i-1) = minor_a/(Bref*q0)*abs_jac(i) 
                     
            CALL cross_product(gradAlpha,gradThetaStar, es) 
            CALL cross_product(gradThetaStar,gradS,     ea) 
            CALL cross_product(gradS,gradAlpha,         et)  
            es(:) = es(:)*jac1
            ea(:) = ea(:)*jac1
            et(:) = et(:)*jac1
            gradB = gradB/Ba
            dBds = DOT_PRODUCT(gradB,es)
            dBda = DOT_PRODUCT(gradB,ea)
            dBdt(i) = DOT_PRODUCT(gradB,et)
            
            c = iot*iot*a**4
            L1(i) = q/sqrt(s)*(dBda + c*(gss*gat-gsa*gst)*dBdt(i)/(4*Bhat(i)**2))
            L2(i) = two*sqrt(s)*(dBds + c*(gaa*gst-gsa*gat)*dBdt(i)/(4*Bhat(i)**2))
            dBdx(i-1) = L1(i)
            dBdy(i-1) = L2(i)
 
          ENDDO ! End loop over field line
 
          ! Solve ITG dispersion relation for each (kx,ky)
          CALL PTSM3D_itg_solve(j,k)
        ENDDO
      ENDDO ! End loop over (kx,ky)

      ! Call the rest of the PTSM3D functions
      CALL PTSM3D_initialize_triplets

      CALL PTSM3D_compute_triplets

      CALL PTSM3D_compute_target
      
      IF (opt_target == 'zf') target_ptsm3d = 1.0/target_12f
      IF (opt_target == 'nzf') target_ptsm3d = 1.0/target_qst
      IF (opt_target == 'combo') target_ptsm3d = &
        & 1.0/(target_12f+target_qst) 

      CALL PTSM3D_finalize_triplets

      CALL PTSM3D_finalize_itg_solve

      CALL PTSM3D_finalize_geom

      END SUBROUTINE stellopt_ptsm3d
