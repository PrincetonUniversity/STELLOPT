!-----------------------------------------------------------------------
!     Subroutine:    stellopt_load_equil
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          06/06/2012
!     Description:   This subroutine handles loading the various
!                    equilibrium parameters produced by the various
!                    equilibrium codes.  Essentailly we're defining
!                    an equilibrium object we can use to evaluate
!                    the properties of
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_load_equil(lscreen,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_input_mod
      USE stellopt_vars
      USE equil_vals
      USE equil_utils
      ! VMEC DATA
      USE read_wout_mod, ONLY: aspect_vmec => aspect, beta_vmec => betatot,&
                               betap_vmec => betapol, betat_vmec => betator,&
                               curtor_vmec => Itor, phi_vmec => phi, &
                               rbtor_vmec => RBtor0, Rmajor_vmec => Rmajor, &
                               Aminor_vmec => Aminor, &
                               ns_vmec => ns, volume_vmec => Volume, &
                               wp_vmec => wp, pres_vmec => presf, &
                               vp_vmec => vp, presh_vmec => pres, &
                               jdotb_vmec => jdotb, &
                               iota_vmec => iotaf, rmnc_vmec => rmnc, &
                               rmns_vmec => rmns, zmnc_vmec => zmnc, &
                               gmns_vmec => gmns, gmnc_vmec => gmnc, &
                               lmns_vmec => lmns, lmnc_vmec => lmnc, &
                               bmns_vmec => bmns, bmnc_vmec => bmnc, &
                               zmns_vmec => zmns, mnmax_vmec => mnmax,&
                               bsupumnc_vmec => bsupumnc, bsupvmnc_vmec => bsupvmnc, &
                               bsupumns_vmec => bsupumns, bsupvmns_vmec => bsupvmns, &
                               xm_vmec => xm, xn_vmec => xn, &
                               xm_nyq_vmec => xm_nyq, xn_nyq_vmec => xn_nyq, &
                               lasym_vmec => lasym, mpol_vmec => mpol,&
                               ntor_vmec => ntor, nfp_vmec => nfp, &
                               extcur_vmec => extcur, &
                               read_wout_file, read_wout_deallocate, &
                               jcurv_vmec => jcurv, rmax_vmec => rmax_surf, &
                               rmin_vmec => rmin_surf, zmax_vmec => zmax_surf, &
                               omega_vmec => omega, machsq_vmec => machsq, &
                               tpotb_vmec => tpotb, b0_vmec => b0
      USE vmec_input,    ONLY: rbc_vmec => rbc, rbs_vmec => rbs, &
                               zbc_vmec => zbc, zbs_vmec => zbs, &
                               raxis_cc_vmec => raxis_cc, raxis_cs_vmec => raxis_cs, &
                               zaxis_cc_vmec => zaxis_cc, zaxis_cs_vmec => zaxis_cs, &
                               lfreeb_vmec => lfreeb
                               
      USE mgrid_mod,     ONLY: nextcur, rminb, rmaxb, zmaxb
      !USE vparams, ONLY: mu0
      ! SPEC
      ! PIES
      ! SIESTA
      ! LIBSTELL
      USE safe_open_mod, ONLY: safe_open
      USE EZspline
      
!-----------------------------------------------------------------------
!     Subroutine Parameters
!        iflag         Error flag
!----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL, INTENT(in)    :: lscreen
      INTEGER, INTENT(inout) :: iflag
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        iunit       File unit number
!----------------------------------------------------------------------
      INTEGER ::  ier, iunit,nvar_in
      INTEGER ::  nu, nv, u, v, mn, dex, mnmax_temp
      INTEGER :: im,in
      INTEGER, ALLOCATABLE :: xm_temp(:), xn_temp(:)
      REAL(rprec) :: temp, s_temp, u_temp, phi_temp
      REAL(rprec), ALLOCATABLE :: xu(:), xv(:)
      REAL(rprec), ALLOCATABLE :: rmnc_temp(:,:), zmns_temp(:,:), lmns_temp(:,:)
      REAL(rprec), ALLOCATABLE :: rmns_temp(:,:), zmnc_temp(:,:), lmnc_temp(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: Vol(:)
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      ier = 0
      SELECT CASE (TRIM(equil_type))
         CASE('vmec2000','animec','flow','satire','parvmec','paravmec','vboot')
            ! Read the VMEC output
            CALL read_wout_deallocate
            CALL read_wout_file(TRIM(proc_string),ier)
            IF (ier .ne. 0) RETURN
            ! Check for grid size
            IF (lfreeb_vmec .and. ((rmax_vmec >= rmaxb) .or. (rmin_vmec <= rminb) .or. (zmax_vmec >= zmaxb))) THEN
               IF (lscreen) WRITE(6,'(A)')   '!!!!!  VMEC Solution exceeds Vacuum Grid Size  !!!!!'
               iflag = -1
               RETURN
            END IF
            ! Logical Values
            lasym = lasym_vmec
            ! Scalar Values
            aspect  = aspect_vmec
            beta    = beta_vmec
            betap   = betap_vmec
            betat   = betat_vmec
            curtor  = curtor_vmec
            phiedge = phi_vmec(ns_vmec)
            volume  = volume_vmec
            rmajor  = Rmajor_vmec
            aminor  = Aminor_vmec
            !wp      = wp_vmec/mu0
            drho    = 1./ns_vmec
            nrad    = ns_vmec
            wp      = 1.5_rprec*pi2*pi2*SUM(vp_vmec(2:nrad)*presh_vmec(2:nrad))/(nrad-1) ! Old STELLOPT Way
            rbtor   = rbtor_vmec
            Baxis   = 0.5*SUM(3*bmnc_vmec(:,2)-bmnc_vmec(:,3),1) ! Asymmetric part zero at phi=0
            ! May need to create some radial arrays
            IF (ALLOCATED(rho)) DEALLOCATE(rho)
            ALLOCATE(rho(ns_vmec))
            FORALL(u=1:ns_vmec) rho(u) = REAL(u-1)/REAL(ns_vmec-1)
            nfp     = nfp_vmec
            ! Get the external currents
            IF (ALLOCATED(extcur)) DEALLOCATE(extcur)
            ALLOCATE(extcur(nextcur))
            extcur(1:nextcur) = extcur_vmec(1:nextcur)
            ! Radial Splines Arrays
            IF (EZspline_allocated(pres_spl)) CALL EZspline_free(pres_spl,iflag)
            IF (EZspline_allocated(iota_spl)) CALL EZspline_free(iota_spl,iflag)
            IF (EZspline_allocated(ip_spl)) CALL EZspline_free(ip_spl,iflag)
            IF (EZspline_allocated(jdotb_spl)) CALL EZspline_free(jdotb_spl,iflag)
            IF (EZspline_allocated(jcurv_spl)) CALL EZspline_free(jcurv_spl,iflag)
            IF (EZspline_allocated(V_spl)) CALL EZspline_free(V_spl,iflag)
            CALL EZspline_init(pres_spl,ns_vmec,bcs0,iflag)
            CALL EZspline_init(iota_spl,ns_vmec,bcs0,iflag)
            CALL EZspline_init(jdotb_spl,ns_vmec,bcs0,iflag)
            CALL EZspline_init(jcurv_spl,ns_vmec,bcs0,iflag)
            CALL EZspline_init(V_spl,ns_vmec,bcs0,iflag)
            !CALL EZspline_init(ip_spl,ns_vmec,bcs0,iflag)
            pres_spl%x1 = rho
            iota_spl%x1 = rho
            jdotb_spl%x1 = rho
            jcurv_spl%x1 = rho
            V_spl%x1 = rho
            !ip_spl%x1 = phi_vmec
            pres_spl%isHermite = 1
            iota_spl%isHermite = 1
            jdotb_spl%isHermite = 1
            jcurv_spl%isHermite = 1
            V_spl%isHermite = 1
            !ip_spl%isHermite = 1
            CALL EZspline_setup(pres_spl,pres_vmec,ier)
            CALL EZspline_setup(iota_spl,iota_vmec,ier)
            CALL EZspline_setup(jdotb_spl,jdotb_vmec,ier)
            CALL EZspline_setup(jcurv_spl,jcurv_vmec,ier)
            ALLOCATE(Vol(ns_vmec))
            FORALL(u=1:ns_vmec) Vol(u) = SUM(vp_vmec(1:u))
            CALL EZspline_setup(V_spl,Vol,ier)
            DEALLOCATE(Vol)
            ! Handle the FLOW arrays
            IF (TRIM(equil_type)=='flow' .or. TRIM(equil_type)=='satire') THEN
               mach0 = SQRT(machsq_vmec)
               IF (EZspline_allocated(omega_spl)) CALL EZspline_free(omega_spl,iflag)
               CALL EZspline_init(omega_spl,ns_vmec,bcs0,iflag)
               omega_spl%x1 = rho
               omega_spl%isHermite = 1
               CALL EZspline_setup(omega_spl,omega_vmec,ier)
            END IF
            !CALL EZspline_setup(ip_spl,pres_vmec,ier)
            ! Get the realspace R and Z and metric elements
            nu = 8 * mpol_vmec + 1
            nu = 2 ** CEILING(log(REAL(nu))/log(2.0_rprec))
            nv = 4 * ntor_vmec + 5                                      ! Use at least 5 toroidal points
            nv = 2 ** CEILING(log(REAL(nv))/log(2.0_rprec)) + 1  ! Odd so we get nfp/2 plane
            ! Handle Nyquist issues
            IF (SIZE(xm_nyq_vmec) > SIZE(xm_vmec)) THEN
               mnmax_temp = SIZE(xm_nyq_vmec)
            ELSE
               mnmax_temp = mnmax_vmec
            END IF
            ALLOCATE(xm_temp(mnmax_temp),xn_temp(mnmax_temp))
            ALLOCATE(rmnc_temp(mnmax_temp,ns_vmec), zmns_temp(mnmax_temp,ns_vmec), lmns_temp(mnmax_temp,ns_vmec))
            rmnc_temp=0; zmns_temp=0; lmns_temp = 0;
            IF (lasym_vmec) THEN
                ALLOCATE(rmns_temp(mnmax_temp,ns_vmec), zmnc_temp(mnmax_temp,ns_vmec), lmnc_temp(mnmax_temp,ns_vmec))
                rmns_temp=0; zmnc_temp=0; lmnc_temp = 0;
            END IF
            IF (SIZE(xm_nyq_vmec) > SIZE(xm_vmec)) THEN
               xm_temp(1:mnmax_temp) = xm_nyq_vmec(1:mnmax_temp)
               xn_temp(1:mnmax_temp) = xn_nyq_vmec(1:mnmax_temp)
               DO u = 1,mnmax_temp
                  DO v = 1, mnmax_vmec
                     IF ((xm_vmec(v) .eq. xm_nyq_vmec(u)) .and. (xn_vmec(v) .eq. xn_nyq_vmec(u))) THEN
                        rmnc_temp(u,1:ns_vmec) = rmnc_vmec(v,1:ns_vmec)
                        zmns_temp(u,1:ns_vmec) = zmns_vmec(v,1:ns_vmec)
                        lmns_temp(u,1:ns_vmec) = lmns_vmec(v,1:ns_vmec)
                        IF (lasym_vmec) THEN
                           rmns_temp(u,1:ns_vmec) = rmns_vmec(v,1:ns_vmec)
                           zmnc_temp(u,1:ns_vmec) = zmnc_vmec(v,1:ns_vmec)
                           lmnc_temp(u,1:ns_vmec) = lmnc_vmec(v,1:ns_vmec)
                        END IF
                     END IF
                  END DO
               END DO
            ELSE
               xm_temp(1:mnmax_temp) = xm_vmec(1:mnmax_temp)
               xn_temp(1:mnmax_temp) = xn_vmec(1:mnmax_temp)
               rmnc_temp(1:mnmax_temp,1:ns_vmec) = rmnc_vmec(1:mnmax_temp,1:ns_vmec)
               zmns_temp(1:mnmax_temp,1:ns_vmec) = zmns_vmec(1:mnmax_temp,1:ns_vmec)
               lmns_temp(1:mnmax_temp,1:ns_vmec) = lmns_vmec(1:mnmax_temp,1:ns_vmec)
               IF (lasym_vmec) THEN
                  rmns_temp(1:mnmax_temp,1:ns_vmec) = rmns_vmec(1:mnmax_temp,1:ns_vmec)
                  zmnc_temp(1:mnmax_temp,1:ns_vmec) = zmnc_vmec(1:mnmax_temp,1:ns_vmec)
                  lmnc_temp(1:mnmax_temp,1:ns_vmec) = lmnc_vmec(1:mnmax_temp,1:ns_vmec)
               END IF
            END IF
            ! Half to full grid
            bsupumnc_vmec(:,1) = (3*bsupumnc_vmec(:,2) - bsupumnc_vmec(:,3))*0.5D+00
            bsupvmnc_vmec(:,1) = (3*bsupvmnc_vmec(:,2) - bsupvmnc_vmec(:,3))*0.5D+00
            gmnc_vmec(:,1) = (3*gmnc_vmec(:,2) - gmnc_vmec(:,3))*0.5D+00
            lmns_temp(:,1) = (3*lmns_temp(:,2) - lmns_temp(:,3))*0.5D+00
            FORALL(mn = 1:mnmax_temp) bsupumnc_vmec(mn,2:ns_vmec-1) = 0.5*(bsupumnc_vmec(mn,2:ns_vmec-1) + bsupumnc_vmec(mn,3:ns_vmec))
            FORALL(mn = 1:mnmax_temp) bsupvmnc_vmec(mn,2:ns_vmec-1) = 0.5*(bsupvmnc_vmec(mn,2:ns_vmec-1) + bsupvmnc_vmec(mn,3:ns_vmec))
            FORALL(mn = 1:mnmax_temp) gmnc_vmec(mn,2:ns_vmec-1) = 0.5*(gmnc_vmec(mn,2:ns_vmec-1) + gmnc_vmec(mn,3:ns_vmec))
            FORALL(mn = 1:mnmax_temp) lmns_temp(mn,2:ns_vmec-1) = 0.5*(lmns_temp(mn,2:ns_vmec-1) + lmns_temp(mn,3:ns_vmec))
            bsupumnc_vmec(:,ns_vmec) = 2*bsupumnc_vmec(:,ns_vmec) - bsupumnc_vmec(:,ns_vmec-1)
            bsupvmnc_vmec(:,ns_vmec) = 2*bsupvmnc_vmec(:,ns_vmec) - bsupvmnc_vmec(:,ns_vmec-1)
            gmnc_vmec(:,ns_vmec)     = 2*gmnc_vmec(:,ns_vmec) - gmnc_vmec(:,ns_vmec-1)
            lmns_temp(:,ns_vmec) = 2*lmns_temp(:,ns_vmec) - lmns_temp(:,ns_vmec-1)
            ! Load STEL_TOOLS
            IF (lasym_vmec) THEN
               bsupumns_vmec(:,1) = 1.5*bsupumns_vmec(:,2) - 0.5*bsupumns_vmec(:,3)
               bsupvmns_vmec(:,1) = 1.5*bsupvmns_vmec(:,2) - 0.5*bsupvmns_vmec(:,3)
               gmns_vmec(:,1)     = (3*gmns_vmec(:,2) - gmns_vmec(:,3))*0.5D+00
               lmnc_temp(:,1) = (3*lmnc_temp(:,2) - lmnc_temp(:,3))*0.5D+00
               FORALL(mn = 1:mnmax_temp) bsupumns_vmec(mn,2:ns_vmec-1) = 0.5*(bsupumns_vmec(mn,2:ns_vmec-1) + bsupumns_vmec(mn,3:ns_vmec))
               FORALL(mn = 1:mnmax_temp) bsupvmns_vmec(mn,2:ns_vmec-1) = 0.5*(bsupvmns_vmec(mn,2:ns_vmec-1) + bsupvmns_vmec(mn,3:ns_vmec))
               FORALL(mn = 1:mnmax_temp) gmns_vmec(mn,2:ns_vmec-1) = 0.5*(gmns_vmec(mn,2:ns_vmec-1) + gmns_vmec(mn,3:ns_vmec))
               FORALL(mn = 1:mnmax_temp) lmnc_vmec(mn,2:ns_vmec-1) = 0.5*(lmnc_vmec(mn,2:ns_vmec-1) + lmnc_vmec(mn,3:ns_vmec))
               bsupumns_vmec(:,ns_vmec) = 2*bsupumns_vmec(:,ns_vmec) - bsupumns_vmec(:,ns_vmec-1)
               bsupvmns_vmec(:,ns_vmec) = 2*bsupvmns_vmec(:,ns_vmec) - bsupvmns_vmec(:,ns_vmec-1)
               gmns_vmec(:,ns_vmec)     = 2*gmns_vmec(:,ns_vmec) - gmns_vmec(:,ns_vmec-1)
               lmnc_temp(:,ns_vmec) = 2*lmnc_temp(:,ns_vmec) - lmnc_temp(:,ns_vmec-1)
               CALL load_fourier_geom(1,ns_vmec,mnmax_temp,nu,nv,INT(xm_temp),INT(-xn_temp),iflag,&
                           DBLE(rmnc_temp),DBLE(zmns_temp),RMNS=DBLE(rmns_temp),ZMNC=DBLE(zmnc_temp),&
                           BUMNC=DBLE(bsupumnc_vmec),BVMNC=DBLE(bsupvmnc_vmec),&
                           BUMNS=DBLE(bsupumns_vmec),BVMNS=DBLE(bsupvmns_vmec),&
                           LMNS=DBLE(lmns_temp),LMNC=DBLE(lmnc_temp),&
                           BMNC=DBLE(bmnc_vmec),BMNS=DBLE(bmns_vmec),&
                           GMNC=DBLE(gmnc_vmec),GMNS=DBLE(gmns_vmec))
               DEALLOCATE(rmns_temp,zmnc_temp,lmnc_temp)
            ELSE
               CALL load_fourier_geom(1,ns_vmec,mnmax_temp,nu,nv,INT(xm_temp),INT(-xn_temp),iflag,&
                           DBLE(rmnc_temp),DBLE(zmns_temp),&
                           BUMNC=DBLE(bsupumnc_vmec),BVMNC=DBLE(bsupvmnc_vmec),&
                           LMNS=DBLE(lmns_temp),&
                           BMNC=DBLE(bmnc_vmec),&
                           GMNC=DBLE(gmnc_vmec))
            END IF
            DEALLOCATE(xm_temp,xn_temp,rmnc_temp,zmns_temp,lmns_temp)
            ! Update R0,Z0
            CALL get_equil_RZ(DBLE(0),DBLE(0),DBLE(0),r0,z0,iflag)
            ! Update boundary coefficients (this only effects the input files as we use a reset anyway)
            ! SAL - 09/19/12  In free boundary runs this can cause problems.
            ! SAL - 04/25/17  In fixed boundary this can cause prolbems because of flip_theta in readin
            !IF (.not. lfreeb_vmec) THEN
            !   DO mn = 1, mnmax_vmec
            !      im = xm_vmec(mn)
            !      in = xn_vmec(mn)/nfp
            !      rbc_vmec(in,im)=rmnc_vmec(mn,ns_vmec)
            !      zbs_vmec(in,im)=zmns_vmec(mn,ns_vmec)
            !      IF (im == 0) THEN
            !         raxis_cc_vmec(abs(in)) = rmnc_vmec(mn,1)
            !         zaxis_cs_vmec(abs(in)) = zmns_vmec(mn,1)
            !      END IF
            !   END DO
            !   IF (lasym_vmec) THEN
            !      DO mn = 1, mnmax_vmec
            !         im = xm_vmec(mn)
            !         in = xn_vmec(mn)/nfp
            !         rbs_vmec(in,im)=rmns_vmec(mn,ns_vmec)
            !         zbc_vmec(in,im)=zmnc_vmec(mn,ns_vmec)
            !         IF (im == 0) THEN
            !            raxis_cs_vmec(abs(in)) = rmns_vmec(mn,1)
            !            zaxis_cc_vmec(abs(in)) = zmnc_vmec(mn,1)
            !         END IF
            !      END DO
            !   END IF
            !END IF
         CASE('spec')
         CASE('pies')
         CASE('siesta')
      END SELECT
      ! Setup the internal STELLOPT arrays
      IF (EZspline_allocated(phi_spl)) CALL EZspline_free(phi_spl,iflag)
      IF (EZspline_allocated(ne_spl)) CALL EZspline_free(ne_spl,iflag)
      IF (EZspline_allocated(te_spl)) CALL EZspline_free(te_spl,iflag)
      IF (EZspline_allocated(ti_spl)) CALL EZspline_free(ti_spl,iflag)
      IF (EZspline_allocated(th_spl)) CALL EZspline_free(th_spl,iflag)
      IF (EZspline_allocated(nustar_spl)) CALL EZspline_free(nustar_spl,iflag)
      IF (EZspline_allocated(zeff_spl)) CALL EZspline_free(zeff_spl,iflag)
      IF (EZspline_allocated(emis_xics_spl)) CALL EZspline_free(emis_xics_spl,iflag)
      dex = MINLOC(phi_aux_s(2:),DIM=1)
      IF (dex > 4) THEN
         CALL EZspline_init(phi_spl,dex,bcs0,iflag)
         phi_spl%x1 = phi_aux_s
         phi_spl%isHermite = 1
         CALL EZspline_setup(phi_spl,phi_aux_f,ier)
      END IF
      dex = MINLOC(ne_aux_s(2:),DIM=1)
      IF (dex > 4) THEN
         CALL EZspline_init(ne_spl,dex,bcs0,iflag)
         ne_spl%x1 = ne_aux_s(1:dex)
         ne_spl%isHermite = 1
         CALL EZspline_setup(ne_spl,ne_aux_f,ier)
      END IF
      dex = MINLOC(te_aux_s(2:),DIM=1)
      IF (dex > 4) THEN
         CALL EZspline_init(te_spl,dex,bcs0,iflag)
         te_spl%x1 = te_aux_s(1:dex)
         te_spl%isHermite = 1
         CALL EZspline_setup(te_spl,te_aux_f,ier)
      END IF
      dex = MINLOC(ti_aux_s(2:),DIM=1)
      IF (dex > 4) THEN
         CALL EZspline_init(ti_spl,dex,bcs0,iflag)
         ti_spl%x1 = ti_aux_s(1:dex)
         ti_spl%isHermite = 1
         CALL EZspline_setup(ti_spl,ti_aux_f,ier)
      END IF
      dex = MINLOC(th_aux_s(2:),DIM=1)
      IF (dex > 4) THEN
         CALL EZspline_init(th_spl,dex,bcs0,iflag)
         th_spl%x1 = th_aux_s(1:dex)
         th_spl%isHermite = 1
         CALL EZspline_setup(th_spl,phi_aux_f,ier)
      END IF
      dex = MINLOC(nustar_s(2:),DIM=1)
      IF (dex > 4) THEN
         CALL EZspline_init(nustar_spl,dex,bcs0,iflag)
         nustar_spl%x1 = nustar_s(1:dex)
         nustar_spl%isHermite = 1
         CALL EZspline_setup(nustar_spl,nustar_f,ier)
      END IF
      dex = MINLOC(zeff_aux_s(2:),DIM=1)
      IF (dex > 4) THEN
         CALL EZspline_init(zeff_spl,dex,bcs0,iflag)
         zeff_spl%x1 = zeff_aux_s(1:dex)
         zeff_spl%isHermite = 1
         CALL EZspline_setup(zeff_spl,zeff_aux_f,ier)
      END IF
      dex = MINLOC(emis_xics_s(2:),DIM=1)
      IF (dex > 4) THEN
         CALL EZspline_init(emis_xics_spl,dex,bcs0,iflag)
         emis_xics_spl%x1 = emis_xics_s(1:dex)
         emis_xics_spl%isHermite = 1
         CALL EZspline_setup(emis_xics_spl,emis_xics_f,ier)
      END IF
      IF (lscreen) THEN
         WRITE(6,'(A,F7.3)')   '     ASPECT RATIO:  ',aspect
         WRITE(6,'(A,F7.3,A)') '             BETA:  ',beta,'  (total)'
         WRITE(6,'(A,F7.3,A)') '                    ',betap,'  (poloidal)'
         WRITE(6,'(A,F7.3,A)') '                    ',betat,'  (toroidal)'
         WRITE(6,'(A,E20.12)') '  TORIDAL CURRENT:  ',curtor
         WRITE(6,'(A,F7.3)')   '     TORIDAL FLUX:  ',phiedge
         WRITE(6,'(A,F7.3)')   '           VOLUME:  ',volume    
         WRITE(6,'(A,F7.3)')   '     MAJOR RADIUS:  ',rmajor
         WRITE(6,'(A,F7.3)')   '     MINOR RADIUS:  ',aminor
         WRITE(6,'(A,F7.3)')   '       AXIS FIELD:  ',Baxis
         WRITE(6,'(A,E20.12)')   '    STORED ENERGY:  ',wp
         CALL FLUSH(6)
      END IF
      
      ! Testing output
      !DO u = 1, 20
      !      temp = 0.0+(0.05*u)
      !      phi_temp = pi2/3
      !      ier=0
      !      CALL get_equil_s(1.8_rprec,phi_temp,temp,s_temp,ier,u_temp)
      !      WRITE(77,*) temp,phi_temp,1.8,s_temp,u_temp,ier
      !      PRINT *,'r=',1.8,'phi=',phi_temp,'z=',temp,'s=',s_temp,ier, u_temp
      !END DO
      !stop
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_load_equil
