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
                               mnmax_nyq_vmec => mnmax_nyq,&
                               xm_nyq_vmec => xm_nyq, xn_nyq_vmec => xn_nyq, &
                               lasym_vmec => lasym, mpol_vmec => mpol,&
                               ntor_vmec => ntor, nfp_vmec => nfp, &
                               extcur_vmec => extcur, &
                               read_wout_file, read_wout_deallocate, &
                               jcurv_vmec => jcurv, rmax_vmec => rmax_surf, &
                               rmin_vmec => rmin_surf, zmax_vmec => zmax_surf, &
                               omega_vmec => omega, machsq_vmec => machsq, &
                               tpotb_vmec => tpotb, b0_vmec => b0, ierr_vmec
      USE vmec_input,    ONLY: rbc_vmec => rbc, rbs_vmec => rbs, &
                               zbc_vmec => zbc, zbs_vmec => zbs, &
                               raxis_cc_vmec => raxis_cc, raxis_cs_vmec => raxis_cs, &
                               zaxis_cc_vmec => zaxis_cc, zaxis_cs_vmec => zaxis_cs, &
                               lfreeb_vmec => lfreeb, &
                               ah_aux_s, ah_aux_f, at_aux_s, at_aux_f 
                               
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
      INTEGER :: im,in, k, ns_extrap
      INTEGER, ALLOCATABLE :: xm_temp(:), xn_temp(:)
      REAL(rprec) :: temp, s_temp, u_temp, phi_temp
      REAL(rprec), ALLOCATABLE :: xu(:), xv(:)
      REAL(rprec), ALLOCATABLE :: rmnc_temp(:,:), zmns_temp(:,:), lmns_temp(:,:)
      REAL(rprec), ALLOCATABLE :: rmns_temp(:,:), zmnc_temp(:,:), lmnc_temp(:,:)
      REAL(rprec), ALLOCATABLE :: gmnc_temp(:,:), bmnc_temp(:,:)
      REAL(rprec), ALLOCATABLE :: bsupumnc_temp(:,:), bsupvmnc_temp(:,:)
      REAL(rprec), ALLOCATABLE :: gmns_temp(:,:), bmns_temp(:,:)
      REAL(rprec), ALLOCATABLE :: bsupumns_temp(:,:), bsupvmns_temp(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: Vol(:)
      DOUBLE PRECISION, ALLOCATABLE :: mfact(:,:)

!      DOUBLE PRECISION, PARAMETER :: sextrap=2.25
      DOUBLE PRECISION, PARAMETER :: sextrap=1.00
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      ier = 0
      SELECT CASE (TRIM(equil_type))
         CASE('vmec2000','animec','flow','satire','parvmec','paravmec','vboot','vmec2000_oneeq')
            ! Read the VMEC output
            CALL read_wout_deallocate
            CALL read_wout_file(TRIM(proc_string),ier)
            IF (ier .ne. 0 .or. ierr_vmec .ne. 0) THEN
               iflag = -1
               RETURN
            END IF
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
            IF (ALLOCATED(shat)) DEALLOCATE(shat)
            ALLOCATE(rho(ns_vmec), shat(ns_vmec))
            FORALL(u=1:ns_vmec) shat(u) = REAL(u-1)/REAL(ns_vmec-1)
            rho = sqrt(shat); rho(1) = 0
            nfp     = nfp_vmec
            ! Get the external currents
            IF (ALLOCATED(extcur)) DEALLOCATE(extcur)
            ALLOCATE(extcur(nextcur))
            extcur(1:nextcur) = extcur_vmec(1:nextcur)
            ! Radial Splines Arrays
            CALL setup_prof_spline(pres_spl,  ns_vmec, shat, pres_vmec, iflag)
            CALL setup_prof_spline(iota_spl,  ns_vmec, shat, iota_vmec, iflag)
            !CALL setup_prof_spline(ip_spl,    ns_vmec, shat, ip_vmec,   iflag)
            CALL setup_prof_spline(jdotb_spl, ns_vmec, shat, jdotb_vmec, iflag)
            CALL setup_prof_spline(jcurv_spl, ns_vmec, shat, jcurv_vmec, iflag)
            ALLOCATE(Vol(ns_vmec))
            FORALL(u=1:ns_vmec) Vol(u) = SUM(vp_vmec(1:u))
            CALL setup_prof_spline(V_spl,     ns_vmec, shat, Vol, iflag)
            DEALLOCATE(Vol)

            ! Handle the FLOW arrays
            IF (TRIM(equil_type)=='flow' .or. TRIM(equil_type)=='satire') THEN
               mach0 = SQRT(machsq_vmec)
               CALL setup_prof_spline(omega_spl, ns_vmec, shat, omega_vmec, iflag)
            END IF

            ! Get the realspace R and Z and metric elements
            nu = 8 * mpol_vmec + 1
            nu = 2 ** CEILING(log(REAL(nu))/log(2.0_rprec))
            nv = 4 * ntor_vmec + 5                                      ! Use at least 5 toroidal points
            nv = 2 ** CEILING(log(REAL(nv))/log(2.0_rprec)) + 1  ! Odd so we get nfp/2 plane

            ! Place the non-nyquist sized arrays on nyquist sized arrays
            ! Allocate helper arrays for 
            ns_extrap = FLOOR(sextrap*(ns_vmec-1)+1)
            mnmax_temp = mnmax_nyq_vmec
            ALLOCATE(xm_temp(mnmax_temp),xn_temp(mnmax_temp))
            ALLOCATE(rmnc_temp(mnmax_temp,ns_extrap), zmns_temp(mnmax_temp,ns_extrap), lmns_temp(mnmax_temp,ns_extrap))
            rmnc_temp=0; zmns_temp=0; lmns_temp = 0
            IF (lasym_vmec) THEN
               ALLOCATE(rmns_temp(mnmax_temp,ns_extrap), zmnc_temp(mnmax_temp,ns_extrap), lmnc_temp(mnmax_temp,ns_extrap))
               rmns_temp=0; zmnc_temp=0; lmnc_temp = 0
            END IF
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

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Extrapolate Full grid quantities
            !    Photographic enlargement
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ALLOCATE(mfact(mnmax_temp,2))
            DO k = ns_vmec+1, ns_extrap
               mfact(:,1)= SQRT(DBLE(k-1)/DBLE(ns_vmec-1))
               rmnc_temp(:,k) = (rmnc_temp(:,ns_vmec)-rmnc_temp(:,1))*mfact(:,1)+rmnc_temp(:,1) 
               zmns_temp(:,k) = (zmns_temp(:,ns_vmec)-zmns_temp(:,1))*mfact(:,1)+zmns_temp(:,1)
               IF (lasym_vmec) THEN
                  rmns_temp(:,k) = (rmns_temp(:,ns_vmec)-rmns_temp(:,1))*mfact(:,1)+rmns_temp(:,1) 
                  zmnc_temp(:,k) = (zmnc_temp(:,ns_vmec)-zmnc_temp(:,1))*mfact(:,1)+zmnc_temp(:,1)
               END IF
            END DO

            ! Allocate Half grid helper quantities
            ALLOCATE(gmnc_temp(mnmax_temp,ns_extrap), bmnc_temp(mnmax_temp,ns_extrap))
            ALLOCATE(bsupumnc_temp(mnmax_temp,ns_extrap), bsupvmnc_temp(mnmax_temp,ns_extrap))
            gmnc_temp = 0; bmnc_temp = 0; bsupumnc_temp = 0; bsupvmnc_temp = 0
            IF (lasym_vmec) THEN
               ALLOCATE(gmns_temp(mnmax_temp,ns_extrap), bmns_temp(mnmax_temp,ns_extrap))
               ALLOCATE(bsupumns_temp(mnmax_temp,ns_extrap), bsupvmns_temp(mnmax_temp,ns_extrap))
               gmns_temp = 0; bmns_temp = 0; bsupumns_temp = 0; bsupvmns_temp = 0
            END IF

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Extrapolate Half grid quantities 
            !    x = [0,1] between bounding gridpoints
            !    whi_even = x
            !    wlo_even = 1-x
            !    whi_odd  = whi_even*sqrt(s/s_hi)
            !    wlo_odd  = wlo_even*sqrt(s/s_lo)
            !    f = f(jlo)*wlo + f(jhi)*whi  Interpolant
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !   First interpolate from half grid to full (respect overwrite indexing)
            DO k = 2, ns_vmec-1
               WHERE (MOD(NINT(REAL(xm_temp(:))),2) .eq. 0)
                  mfact(:,1)= 0.5
                  mfact(:,2)= 0.5
               ELSEWHERE
                  mfact(:,1)= 0.5*SQRT((k-1.0)/(k-1.5)) !rho/rholo
                  mfact(:,2)= 0.5*SQRT((k-1.0)/(k-0.5)) !rho/rhohi
               ENDWHERE
               lmns_temp(:,k) = mfact(:,1)*lmns_temp(:,k)+mfact(:,2)*lmns_temp(:,k+1)
               gmnc_temp(:,k) = mfact(:,1)*gmnc_vmec(:,k)+mfact(:,2)*gmnc_vmec(:,k+1)
               bmnc_temp(:,k) = mfact(:,1)*bmnc_vmec(:,k)+mfact(:,2)*bmnc_vmec(:,k+1)
               bsupumnc_temp(:,k) = mfact(:,1)*bsupumnc_vmec(:,k)+mfact(:,2)*bsupumnc_vmec(:,k+1)
               bsupvmnc_temp(:,k) = mfact(:,1)*bsupvmnc_vmec(:,k)+mfact(:,2)*bsupvmnc_vmec(:,k+1)
               IF (lasym_vmec) THEN
                  lmnc_temp(:,k) = mfact(:,1)*lmnc_temp(:,k)+mfact(:,2)*lmnc_temp(:,k+1)
                  gmns_temp(:,k) = mfact(:,1)*gmns_vmec(:,k)+mfact(:,2)*gmns_vmec(:,k+1)
                  bmns_temp(:,k) = mfact(:,1)*bmns_vmec(:,k)+mfact(:,2)*bmns_vmec(:,k+1)
                  bsupumns_temp(:,k) = mfact(:,1)*bsupumns_vmec(:,k)+mfact(:,2)*bsupumns_vmec(:,k+1)
                  bsupvmns_temp(:,k) = mfact(:,1)*bsupvmns_vmec(:,k)+mfact(:,2)*bsupvmns_vmec(:,k+1)
               END IF
            END DO

            !   Second, extrapolate to ns
            !       note that ns-1 is full grid but ns is on half grid
            k = ns_vmec
            WHERE (MOD(NINT(REAL(xm_temp(:))),2) .eq. 0)
               mfact(:,1)= 2.0 ! ns (half grid point)
               mfact(:,2)=-1.0 ! ns-1 (full grid point)
            ELSEWHERE
               mfact(:,1)= 2.0*SQRT((ns_vmec-1)/(k-1.5))
               mfact(:,2)=-1.0*SQRT((ns_vmec-1)/(k-2.0))
            ENDWHERE
            lmns_temp(:,k) = mfact(:,1)*lmns_temp(:,k)+mfact(:,2)*lmns_temp(:,k-1)
            gmnc_temp(:,k) = mfact(:,1)*gmnc_vmec(:,k)+mfact(:,2)*gmnc_vmec(:,k-1)
            bmnc_temp(:,k) = mfact(:,1)*bmnc_vmec(:,k)+mfact(:,2)*bmnc_vmec(:,k-1)
            bsupumnc_temp(:,k) = mfact(:,1)*bsupumnc_vmec(:,k)+mfact(:,2)*bsupumnc_vmec(:,k-1)
            bsupvmnc_temp(:,k) = mfact(:,1)*bsupvmnc_vmec(:,k)+mfact(:,2)*bsupvmnc_vmec(:,k-1)
            IF (lasym_vmec) THEN
               lmnc_temp(:,k) = mfact(:,1)*lmnc_temp(:,k)+mfact(:,2)*lmnc_temp(:,k-1)
               gmns_temp(:,k) = mfact(:,1)*gmns_vmec(:,k)+mfact(:,2)*gmns_vmec(:,k-1)
               bmns_temp(:,k) = mfact(:,1)*bmns_vmec(:,k)+mfact(:,2)*bmns_vmec(:,k-1)
               bsupumns_temp(:,k) = mfact(:,1)*bsupumns_vmec(:,k)+mfact(:,2)*bsupumns_vmec(:,k-1)
               bsupvmns_temp(:,k) = mfact(:,1)*bsupvmns_vmec(:,k)+mfact(:,2)*bsupvmns_vmec(:,k-1)
            
            !   Third extrapolate to axis This is the Samantha Lazerson bugfix 2023.01.29
            k = 1
            WHERE (MOD(NINT(REAL(xm_temp(:))),2) .eq. 0)
               mfact(:,1)= 2.0
               mfact(:,2)=-1.0
            ELSEWHERE
               mfact(:,1)= 0
               mfact(:,2)= 0
            ENDWHERE
            lmns_temp(:,k) = mfact(:,1)*lmns_temp(:,k+1)+mfact(:,2)*lmns_temp(:,k+2)
            gmnc_temp(:,k) = mfact(:,1)*gmnc_vmec(:,k+1)+mfact(:,2)*gmnc_vmec(:,k+2)
            bmnc_temp(:,k) = mfact(:,1)*bmnc_vmec(:,k+1)+mfact(:,2)*bmnc_vmec(:,k+2)
            bsupumnc_temp(:,k) = mfact(:,1)*bsupumnc_vmec(:,k+1)+mfact(:,2)*bsupumnc_vmec(:,k+2)
            bsupvmnc_temp(:,k) = mfact(:,1)*bsupvmnc_vmec(:,k+1)+mfact(:,2)*bsupvmnc_vmec(:,k+2)
            IF (lasym_vmec) THEN
               lmnc_temp(:,k) = mfact(:,1)*lmnc_temp(:,k+1)+mfact(:,2)*lmnc_temp(:,k+2)
               gmns_temp(:,k) = mfact(:,1)*gmns_vmec(:,k+1)+mfact(:,2)*gmns_vmec(:,k+2)
               bmns_temp(:,k) = mfact(:,1)*bmns_vmec(:,k+1)+mfact(:,2)*bmns_vmec(:,k+2)
               bsupumns_temp(:,k) = mfact(:,1)*bsupumns_vmec(:,k+1)+mfact(:,2)*bsupumns_vmec(:,k+2)
               bsupvmns_temp(:,k) = mfact(:,1)*bsupvmns_vmec(:,k+1)+mfact(:,2)*bsupvmns_vmec(:,k+2)
            END IF
            END IF
            DEALLOCATE(mfact)

            ! Load STEL_TOOLS
            IF (lasym_vmec) THEN
               CALL load_fourier_geom(1,ns_extrap,mnmax_temp,nu,nv,INT(xm_temp),INT(-xn_temp),iflag,&
                           DBLE(rmnc_temp),DBLE(zmns_temp),RMNS=DBLE(rmns_temp),ZMNC=DBLE(zmnc_temp),&
                           BUMNC=DBLE(bsupumnc_temp),BVMNC=DBLE(bsupvmnc_temp),&
                           BUMNS=DBLE(bsupumns_temp),BVMNS=DBLE(bsupvmns_temp),&
                           LMNS=DBLE(lmns_temp),LMNC=DBLE(lmnc_temp),&
                           BMNC=DBLE(bmnc_temp),BMNS=DBLE(bmns_temp),&
                           GMNC=DBLE(gmnc_temp),GMNS=DBLE(gmns_temp),SMAX=sextrap)
               DEALLOCATE(rmns_temp,zmnc_temp,lmnc_temp,&
                  bmns_temp,gmns_temp,bsupumns_temp,bsupvmns_temp)
            ELSE
               CALL load_fourier_geom(1,ns_extrap,mnmax_temp,nu,nv,INT(xm_temp),INT(-xn_temp),iflag,&
                           DBLE(rmnc_temp),DBLE(zmns_temp),&
                           BUMNC=DBLE(bsupumnc_temp),BVMNC=DBLE(bsupvmnc_temp),&
                           LMNS=DBLE(lmns_temp),&
                           BMNC=DBLE(bmnc_temp),&
                           GMNC=DBLE(gmnc_temp),SMAX=sextrap)
            END IF
            DEALLOCATE(xm_temp,xn_temp,rmnc_temp,zmns_temp,lmns_temp,&
               bmnc_temp,gmnc_temp,bsupumnc_temp,bsupvmnc_temp)
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
         CASE('test')
            ! Do nothing
      END SELECT
      ! Setup the internal STELLOPT arrays
      dex = MINLOC(phi_aux_s(2:),DIM=1)
      IF (dex > 2) CALL setup_prof_spline(phi_spl,dex,phi_aux_s(1:dex),phi_aux_f(1:dex),ier)
      dex = MINLOC(ne_aux_s(2:),DIM=1)
      IF (dex > 2) CALL setup_prof_spline(ne_spl,dex,ne_aux_s(1:dex),ne_aux_f(1:dex),ier)
      dex = MINLOC(te_aux_s(2:),DIM=1)
      IF (dex > 2) CALL setup_prof_spline(te_spl,dex,te_aux_s(1:dex),te_aux_f(1:dex),ier)
      dex = MINLOC(ti_aux_s(2:),DIM=1)
      IF (dex > 2) CALL setup_prof_spline(ti_spl,dex,ti_aux_s(1:dex),ti_aux_f(1:dex),ier)
      dex = MINLOC(th_aux_s(2:),DIM=1)
      IF (dex > 2) CALL setup_prof_spline(th_spl,dex,th_aux_s(1:dex),th_aux_f(1:dex),ier)
      dex = MINLOC(nustar_s(2:),DIM=1)
      IF (dex > 2) CALL setup_prof_spline(nustar_spl,dex,nustar_s(1:dex),nustar_f(1:dex),ier)
      dex = MINLOC(zeff_aux_s(2:),DIM=1)
      IF (dex > 2) CALL setup_prof_spline(zeff_spl,dex,zeff_aux_s(1:dex),zeff_aux_f(1:dex),ier)
      dex = MINLOC(ah_aux_s(2:),DIM=1)
      IF (dex > 2) CALL setup_prof_spline(ah_spl,dex,ah_aux_s(1:dex),ah_aux_f(1:dex),ier)
      dex = MINLOC(emis_xics_s(2:),DIM=1)
      IF (dex > 2) CALL setup_prof_spline(emis_xics_spl,dex,emis_xics_s(1:dex),emis_xics_f(1:dex),ier)

      ! Screen output
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

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_load_equil
