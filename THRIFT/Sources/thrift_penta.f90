!-----------------------------------------------------------------------
!     Subroutine:    thrift_penta
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          08/22/2024
!     Description:   This subroutine calculates the PENTA data.
!-----------------------------------------------------------------------
      SUBROUTINE thrift_penta(lscreen,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE thrift_runtime
      USE thrift_vars, nrho_thrift => nrho
      USE thrift_profiles_mod
      USE thrift_equil
      USE thrift_funcs
      USE phys_const, ONLY: p_mass
      USE penta_interface_mod
      USE mpi_params
      USE mpi_inc
      USE thrift_plasma_solver_mod, ONLY: Dn_NEO,cn_NEO,Dp_NEO,cp_NEO,&
      rho_plasma_grid,Nr_plasma_solver,mytimestep_plasma_solver,&
      G_NEO_complet,Q_NEO_complet,Nt_total_plasma_solver
!-----------------------------------------------------------------------
!     Subroutine Parameters
!        lscreen       Screen output
!        iflag         Error flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL, INTENT(in)    :: lscreen
      INTEGER, INTENT(inout) :: iflag
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      INTEGER :: ns_dkes, k, ier, j, i, ncstar, nestar, mystart, myend, &
                 mysurf, root_max_Er, jspecies
      REAL(rprec) :: s, rho, mytime
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: rho_k, iota, phip, chip, btheta, bzeta, bsq, vp, &
                        te, ne, dtedrho, dnedrho, EparB, JBS_PENTA, etapar_PENTA, Er_PENTA, rho_temp, J_temp, eta_temp, Er_temp
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: GNEO_PENTA, QNEO_PENTA, GNEO_temp, QNEO_temp
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: Dn_PENTA, cn_PENTA, Dn_temp, cn_temp
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: Dp_PENTA, cp_PENTA, Dp_temp, cp_temp
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: ni,ti, dtidrho, dnidrho
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: D11, D13, D33
      TYPE(EZspline1_r8) :: EparB_spl, J_spl, eta_spl, Er_spl, GNEO_spl, QNEO_spl
      TYPE(EZspline1_r8) :: Dn_spl, cn_spl, Dp_spl, cp_spl
      INTEGER :: bcs0(2)
      CHARACTER(LEN=32) :: temp_str, temp1_str
!-----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!-----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (lscreen) WRITE(6,'(a)') ' --------------------  NEOCLASSICAL BOOTSTRAP USING PENTA  -------------------'
      IF (lscreen) Write(*,*) " <r>/<a>","   Er root(s) (V/cm)"

      IF (lvmec) THEN
         ierr_mpi = 0
         ! PENTA is parallelized over radial surfaces in this routine.
         ns_dkes = 0
         DO k = 1, DKES_NS_MAX
            IF ((DKES_K(k) > 0)) ns_dkes = ns_dkes+1
         END DO

         ALLOCATE(rho_k(ns_dkes),iota(ns_dkes),phip(ns_dkes),chip(ns_dkes),btheta(ns_dkes),bzeta(ns_dkes),bsq(ns_dkes),vp(ns_dkes),EparB(ns_dkes))
         ALLOCATE(te(ns_dkes),ne(ns_dkes),dtedrho(ns_dkes),dnedrho(ns_dkes))
         ALLOCATE(ni(ns_dkes,nion_prof),ti(ns_dkes,nion_prof),dtidrho(ns_dkes,nion_prof),dnidrho(ns_dkes,nion_prof))
         
         IF (myworkid == master) THEN

            ALLOCATE(JBS_PENTA(ns_dkes),etapar_PENTA(ns_dkes),Er_PENTA(ns_dkes))
            ALLOCATE(GNEO_PENTA(nion_prof+1,ns_dkes),QNEO_PENTA(nion_prof+1,ns_dkes))
            ALLOCATE(Dn_PENTA(nion_prof+1,ns_dkes),cn_PENTA(nion_prof+1,ns_dkes))
            ALLOCATE(Dp_PENTA(nion_prof+1,ns_dkes),cp_PENTA(nion_prof+1,ns_dkes))

            JBS_PENTA = 0.0; etapar_PENTA = 0.0; Er_PENTA = 0.0
            GNEO_PENTA = 0.0; QNEO_PENTA = 0.0
            Dn_PENTA = 0.0; cn_PENTA = 0.0
            Dp_PENTA = 0.0; cp_PENTA = 0.0

            mytime = THRIFT_T(mytimestep)
            
            ! EparB Spline
            bcs1=(/ 0, 0/)
            CALL EZspline_init(EparB_spl,nsj,bcs1,ier)
            IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'thrift_penta: eparb',ier)
            EparB_spl%isHermite   = 0
            EparB_spl%x1 = SQRT(THRIFT_S)
            CALL EZspline_setup(EparB_spl,THRIFT_EPARB(:,mytimestep),ier,EXACT_DIM=.true.)
            IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'thrift_penta: eparb',ier)
            !

            DO k = 1, ns_dkes
                  mysurf = DKES_K(k)
                  s = DBLE(mysurf-1) / DBLE(ns_eq-1)
                  rho = SQRT(s)
                  rho_k(k) = rho

                  ier = 0
                  CALL EZSpline_interp(iota_spl, rho, iota(k), ier)
                  ier = 0
                  CALL EZSpline_interp(phip_spl, rho, phip(k), ier)
                  ! PENTA wants d/ds note that this won't work if rho=0
                  phip(k) = 0.5*phip(k)/rho
                  chip(k) = iota(k) * phip(k)

                  ier = 0
                  CALL EZSpline_interp(bu_spl, rho, btheta(k), ier)
                  ier = 0
                  CALL EZSpline_interp(bv_spl, rho, bzeta(k), ier)
                  ier = 0
                  CALL EZSpline_interp(bsq_spl, rho, bsq(k), ier)
                  ier = 0

                  CALL EZSpline_interp(vp_spl, rho, vp(k), ier)
                  ! Vp = dV/dPHI need dVds with VMEC normalization
                  vp(k) = vp(k) * phip(k) /(pi2*pi2)

                  ! Profiles
                  CALL get_prof_te(rho, THRIFT_T(mytimestep), te(k))
                  CALL get_prof_ne(rho, THRIFT_T(mytimestep), ne(k))
                  CALL get_prof_teprime(rho, THRIFT_T(mytimestep), dtedrho(k))
                  CALL get_prof_neprime(rho, THRIFT_T(mytimestep), dnedrho(k))
                  DO j = 1, nion_prof
                        CALL get_prof_ti(rho, THRIFT_T(mytimestep), j, ti(k,j))
                        CALL get_prof_ni(rho, THRIFT_T(mytimestep), j, ni(k,j))
                        CALL get_prof_tiprime(rho, THRIFT_T(mytimestep), j, dtidrho(k,j))
                        CALL get_prof_niprime(rho, THRIFT_T(mytimestep), j, dnidrho(k,j))
                  END DO
                  
                  ! EparB
                  ier = 0
                  CALL EZSpline_interp(EparB_spl, rho, EparB(k), ier)
            END DO
            CALL EZspline_free(EparB_spl,ier)
         END IF
         
#if defined(MPI_OPT)
         CALL MPI_BCAST(rho_k,ns_dkes,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(iota,ns_dkes,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(phip,ns_dkes,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(btheta,ns_dkes,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(bzeta,ns_dkes,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(bsq,ns_dkes,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(vp,ns_dkes,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(te,ns_dkes,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(ne,ns_dkes,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(dtedrho,ns_dkes,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(dnedrho,ns_dkes,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(ti,ns_dkes*nion_prof,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(ni,ns_dkes*nion_prof,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(dtidrho,ns_dkes*nion_prof,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(dnidrho,ns_dkes*nion_prof,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(EparB,ns_dkes,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
         ! VMEC quantities
         CALL MPI_BCAST(eq_Aminor,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(eq_Rmajor,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
         ! THRIFT quantities
         CALL MPI_BCAST(mytime,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(mytimestep,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
         
#endif
      
         ! DKES Data
         ncstar = COUNT(DKES_NUSTAR < 1E10)
         nestar = COUNT(DKES_ERSTAR < 1E10)

         DO k = 1,ns_dkes !mystart,myend
            ! PENTA
            CALL PENTA_SET_ION_PARAMS(nion_prof, DBLE(Zatom_prof), Matom_prof/p_mass)
            CALL PENTA_SET_COMMANDLINE(Er_min_Vcm,Er_max_Vcm,DKES_K(k),1,EparB(k),1,'','','')
            CALL PENTA_ALLOCATE_SPECIES
            CALL PENTA_SET_EQ_DATA(rho_k(k),eq_Aminor,eq_Rmajor,vp(k),chip(k),phip(k),iota(k),btheta(k),bzeta(k),bsq(k))
            CALL PENTA_SET_PPROF(ne(k),dnedrho(k)/eq_Aminor,te(k),dtedrho(k)/eq_Aminor,ni(k,:),dnidrho(k,:)/eq_Aminor,ti(k,:),dtidrho(k,:)/eq_Aminor)
            ! MAKE CORRECTIONS ON D31 AND D33 -- values coming from DKES2 miss Bsq factors (see J. Lore documentation)
            CALL PENTA_SET_DKES_STAR(ncstar,nestar,DKES_NUSTAR(1:ncstar),DKES_ERSTAR(1:nestar), &
               DKES_D11(k,:,:), DKES_D31(k,:,:)*SQRT(bsq(k)), DKES_D33(k,:,:)*bsq(k))
            CALL PENTA_SET_BEAM(0.0_rprec) ! Zero becasue we don't read
            CALL PENTA_SET_U2() ! Leave blank for default value
            CALL PENTA_READ_INPUT_FILES(.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.)
            CALL PENTA_SCREEN_INFO
            CALL PENTA_ALLOCATE_DKESCOEFF
            CALL PENTA_FIT_DXX_COEF

            WRITE(temp_str,'(i4.4)') k
            WRITE(temp1_str,'(i3.3)') mytimestep
            CALL PENTA_OPEN_OUTPUT(TRIM(temp1_str) // '_k' // TRIM(temp_str))
            CALL PENTA_FIT_RAD_TRANS

            CALL MPI_BARRIER(MPI_COMM_MYWORLD,ierr_mpi)
            CALL PENTA_RUN_2_EFIELD_3_FIND_ROOTS
            CALL MPI_BARRIER(MPI_COMM_MYWORLD,ierr_mpi)

            IF (myworkid == master) THEN
                  CALL PENTA_RUN_4_AMBIPOLAR

                  ! Save JBS corresponding to the root that has the largest Er
                  ! This because whenever there are 2 stable roots, a rule of thumb is to pick the one with largest Er
                  ! root_max_Er = MAXLOC(Er_roots(1:num_roots),1)
                  ! JBS_PENTA(k) = J_BS_ambi(root_max_Er)
                  ! etapar_PENTA(k) = 1.0_rprec / sigma_par_ambi(root_max_Er)
                  ! Er_PENTA(k) = MAXVAL(Er_roots(1:num_roots),1)
                  
                  ! The call to ROOT_ANALYSIS sets the array 'root_type' which decides which root will settle according to 
                  ! Maxwell construction criterium (see eg. Turkin et al. PoP 18, 022505, 2011)
                  ! This criterium substitutes the above (now commented) lines where the selected root corresponded to
                  ! the largest Er
                  CALL ROOT_ANALYSIS
                  ! Using root_type, pick the ambipolar root that will be saved by THRIFT
                  DO i=1,num_roots
                        IF(root_type(i)) THEN
                              JBS_PENTA(k) = J_BS_ambi(i)
                              etapar_PENTA(k) = 1.0_rprec / sigma_par_ambi(i)
                              Er_PENTA(k) = Er_roots(i)
                              ! NEO fluxes
                              GNEO_PENTA(:,k) = Gammas_ambi(:,i)
                              QNEO_PENTA(:,k) = QoTs_ambi(:,i) * Temps * e_charge
                              ! NEO particle transport coefficients
                              Dn_PENTA(:,k) = (/ (-L_n_ambi(i,j,j), j=1,nion_prof+1) /)
                              cn_PENTA(:,k) = MATMUL(L_T_ambi(i,:,:),e_charge*dTdrs) / dens + Er_PENTA(k)*SUM(L_Er_ambi(i,:,:),dim=2) / dens
                              ! NEO heat transport coefficients
                              Dp_PENTA(:,k) = (/ (-R_T_ambi(i,j,j), j=1,nion_prof+1) /) * e_charge*Temps / dens 
                              cp_PENTA(:,k) = MATMUL(R_n_ambi(i,:,:),dndrs)/dens - MATMUL(R_T_ambi(i,:,:),(e_charge*Temps/dens)*dndrs)/dens &
                                          + Er_PENTA(k)*SUM(R_Er_ambi(i,:,:),dim=2) / dens
                              EXIT
                        ENDIF
                  END DO
            END IF

            CALL MPI_BARRIER(MPI_COMM_MYWORLD,ierr_mpi)
            CALL PENTA_RUN_5_CLEANUP(lscreen)


         END DO
         
         IF (myworkid == master) THEN

            IF(save_all_ambipolar_roots) CALL PENTA_MERGE_AMBIPOLAR_FILES(ns_dkes,temp1_str,mytime)
            IF(save_fluxes_vs_Er) CALL PENTA_MERGE_FLUXES_VS_ER_FILES(ns_dkes,temp1_str,mytime)

            ! Interpolate JBS_PENTA, etapar_PENTA and Er_PENTA at rho=0 and rho=1
            ALLOCATE(J_temp(ns_dkes+2),eta_temp(ns_dkes+2),Er_temp(ns_dkes+2),rho_temp(ns_dkes+2))
            ALLOCATE(GNEO_temp(nion_prof+1,ns_dkes+2),QNEO_temp(nion_prof+1,ns_dkes+2))
            rho_temp(1)        = 0.0
            rho_temp(2:ns_dkes+1) = rho_k
            rho_temp(ns_dkes+2)   = 1.0
            !
            J_temp(2:ns_dkes+1)   = JBS_PENTA
            J_temp(1)             = J_temp(2) - (J_temp(3)-J_temp(2)) * rho_temp(2) / (rho_temp(3)-rho_temp(2))
            J_temp(ns_dkes+2)     = JBS_PENTA(ns_dkes-1) + (JBS_PENTA(ns_dkes)-JBS_PENTA(ns_dkes-1)) * (1-rho_k(ns_dkes-1)) / (rho_k(ns_dkes)-rho_k(ns_dkes-1))
            !
            eta_temp(2:ns_dkes+1)   = etapar_PENTA
            eta_temp(1)             = eta_temp(2) - (eta_temp(3)-eta_temp(2)) * rho_temp(2) / (rho_temp(3)-rho_temp(2))
            eta_temp(ns_dkes+2)     = etapar_PENTA(ns_dkes-1) + (etapar_PENTA(ns_dkes)-etapar_PENTA(ns_dkes-1)) * (1-rho_k(ns_dkes-1)) / (rho_k(ns_dkes)-rho_k(ns_dkes-1))
            !
            Er_temp(2:ns_dkes+1)   = Er_PENTA
            Er_temp(1)             = Er_temp(2) - (Er_temp(3)-Er_temp(2)) * rho_temp(2) / (rho_temp(3)-rho_temp(2))
            Er_temp(ns_dkes+2)     = Er_PENTA(ns_dkes-1) + (Er_PENTA(ns_dkes)-Er_PENTA(ns_dkes-1)) * (1-rho_k(ns_dkes-1)) / (rho_k(ns_dkes)-rho_k(ns_dkes-1))
            !
            GNEO_temp(:,2:ns_dkes+1) = GNEO_PENTA
            GNEO_temp(:,1)           = GNEO_temp(:,2) - (GNEO_temp(:,3)-GNEO_temp(:,2)) * rho_temp(2) / (rho_temp(3)-rho_temp(2))
            GNEO_temp(:,ns_dkes+2)   = GNEO_PENTA(:,ns_dkes-1) + (GNEO_PENTA(:,ns_dkes)-GNEO_PENTA(:,ns_dkes-1)) * (1-rho_k(ns_dkes-1)) / (rho_k(ns_dkes)-rho_k(ns_dkes-1))
            !
            QNEO_temp(:,2:ns_dkes+1) = QNEO_PENTA
            QNEO_temp(:,1)           = QNEO_temp(:,2) - (QNEO_temp(:,3)-QNEO_temp(:,2)) * rho_temp(2) / (rho_temp(3)-rho_temp(2))
            QNEO_temp(:,ns_dkes+2)   = QNEO_PENTA(:,ns_dkes-1) + (QNEO_PENTA(:,ns_dkes)-QNEO_PENTA(:,ns_dkes-1)) * (1-rho_k(ns_dkes-1)) / (rho_k(ns_dkes)-rho_k(ns_dkes-1))

            ! Splines
            bcs0=(/ 0, 0/)
            !JBS
            CALL EZspline_init(J_spl,ns_dkes+2,bcs0,ier)
            J_spl%x1        = rho_temp
            J_spl%isHermite = 0
            CALL EZspline_setup(J_spl,J_temp,ier,EXACT_DIM=.true.)
            !etapar
            CALL EZspline_init(eta_spl,ns_dkes+2,bcs0,ier)
            eta_spl%x1        = rho_temp
            eta_spl%isHermite = 0
            CALL EZspline_setup(eta_spl,eta_temp,ier,EXACT_DIM=.true.)
            !Er
            CALL EZspline_init(Er_spl,ns_dkes+2,bcs0,ier)
            Er_spl%x1        = rho_temp
            Er_spl%isHermite = 0
            CALL EZspline_setup(Er_spl,Er_temp,ier,EXACT_DIM=.true.)
            !
            DEALLOCATE(J_temp,eta_temp,Er_temp)

            ! Calculate J_BS, etapara and Er in THRFIT GRID
            CALL EZspline_interp(J_spl,nsj,SQRT(THRIFT_S),THRIFT_JBOOT(:,mytimestep),ier)
            CALL EZspline_interp(Er_spl,nsj,SQRT(THRIFT_S),THRIFT_ER(:,mytimestep),ier)
            IF( etapar_type == 'dkespenta') CALL EZspline_interp(eta_spl,nsj,SQRT(THRIFT_S),THRIFT_ETAPARA(:,mytimestep),ier)

            CALL EZspline_free(J_spl,ier)
            CALL EZspline_free(eta_spl,ier)
            CALL EZspline_free(Er_spl,ier)

            ! Spline of GNEO and QNEO; computation at THRIFT GRID
            DO jspecies=1,(nion_prof+1)
                  !GNEO
                  CALL EZspline_init(GNEO_spl,ns_dkes+2,bcs0,ier)
                  GNEO_spl%x1        = rho_temp
                  GNEO_spl%isHermite = 1
                  CALL EZspline_setup(GNEO_spl,GNEO_temp(jspecies,:),ier,EXACT_DIM=.true.)
                  !QNEO
                  CALL EZspline_init(QNEO_spl,ns_dkes+2,bcs0,ier)
                  QNEO_spl%x1        = rho_temp
                  QNEO_spl%isHermite = 1
                  CALL EZspline_setup(QNEO_spl,QNEO_temp(jspecies,:),ier,EXACT_DIM=.true.)
                  
                  ! Compute at THRIFT GRID
                  CALL EZspline_interp(GNEO_spl,nsj,SQRT(THRIFT_S),THRIFT_GNEO(jspecies,:,mytimestep),ier)
                  CALL EZspline_interp(QNEO_spl,nsj,SQRT(THRIFT_S),THRIFT_QNEO(jspecies,:,mytimestep),ier)

                  ! Save GNE0_complet and QNEO_complet at plasma solver grid
                  IF(solve_plasma_equations) THEN
                        CALL EZspline_interp(GNEO_spl,Nr_plasma_solver,rho_plasma_grid,G_NEO_complet(jspecies,mytimestep_plasma_solver,:),ier)
                        CALL EZspline_interp(QNEO_spl,Nr_plasma_solver,rho_plasma_grid,Q_NEO_complet(jspecies,mytimestep_plasma_solver,:),ier)
                  END IF

                  ! Deallocate splines
                  CALL EZspline_free(GNEO_spl,ier)
                  CALL EZspline_free(QNEO_spl,ier)
            END DO

            ! Compute NEO coefficients in plasma grid if transport equations being solved
            IF(solve_plasma_equations) THEN
                  ! Interpolate JBS_PENTA, etapar_PENTA and Er_PENTA at rho=0 and rho=1
                  ALLOCATE(Dn_temp(nion_prof+1,ns_dkes+2),cn_temp(nion_prof+1,ns_dkes+2))
                  ALLOCATE(Dp_temp(nion_prof+1,ns_dkes+2),cp_temp(nion_prof+1,ns_dkes+2))
                  !
                  Dn_temp(:,2:ns_dkes+1)   = Dn_PENTA
                  Dn_temp(:,1)             = 0.0_rprec
                  Dn_temp(:,ns_dkes+2)     = Dn_PENTA(:,ns_dkes-1) + (Dn_PENTA(:,ns_dkes)-Dn_PENTA(:,ns_dkes-1)) * (1-rho_k(ns_dkes-1)) / (rho_k(ns_dkes)-rho_k(ns_dkes-1))
                  !
                  cn_temp(:,2:ns_dkes+1)   = cn_PENTA
                  cn_temp(:,1)             = 0.0_rprec
                  cn_temp(:,ns_dkes+2)     = cn_PENTA(:,ns_dkes-1) + (cn_PENTA(:,ns_dkes)-cn_PENTA(:,ns_dkes-1)) * (1-rho_k(ns_dkes-1)) / (rho_k(ns_dkes)-rho_k(ns_dkes-1))
                  !
                  Dp_temp(:,2:ns_dkes+1)   = Dp_PENTA
                  Dp_temp(:,1)             = 0.0_rprec
                  Dp_temp(:,ns_dkes+2)     = Dp_PENTA(:,ns_dkes-1) + (Dp_PENTA(:,ns_dkes)-Dp_PENTA(:,ns_dkes-1)) * (1-rho_k(ns_dkes-1)) / (rho_k(ns_dkes)-rho_k(ns_dkes-1))
                  !
                  cp_temp(:,2:ns_dkes+1)   = cp_PENTA
                  cp_temp(:,1)             = 0.0_rprec
                  cp_temp(:,ns_dkes+2)     = cp_PENTA(:,ns_dkes-1) + (cp_PENTA(:,ns_dkes)-cp_PENTA(:,ns_dkes-1)) * (1-rho_k(ns_dkes-1)) / (rho_k(ns_dkes)-rho_k(ns_dkes-1))
                  !
                  ! Spline of Dn,cn,Dp,cp; computation at plasma grid
                  DO jspecies=1,(nion_prof+1)
                        !Dn
                        CALL EZspline_init(Dn_spl,ns_dkes+2,bcs0,ier)
                        Dn_spl%x1        = rho_temp
                        Dn_spl%isHermite = 1
                        CALL EZspline_setup(Dn_spl,Dn_temp(jspecies,:),ier,EXACT_DIM=.true.)
                        !cn
                        CALL EZspline_init(cn_spl,ns_dkes+2,bcs0,ier)
                        cn_spl%x1        = rho_temp
                        cn_spl%isHermite = 1
                        CALL EZspline_setup(cn_spl,cn_temp(jspecies,:),ier,EXACT_DIM=.true.)
                        !Dp
                        CALL EZspline_init(Dp_spl,ns_dkes+2,bcs0,ier)
                        Dp_spl%x1        = rho_temp
                        Dp_spl%isHermite = 1
                        CALL EZspline_setup(Dp_spl,Dp_temp(jspecies,:),ier,EXACT_DIM=.true.)
                        !cp
                        CALL EZspline_init(cp_spl,ns_dkes+2,bcs0,ier)
                        cp_spl%x1        = rho_temp
                        cp_spl%isHermite = 1
                        CALL EZspline_setup(cp_spl,cp_temp(jspecies,:),ier,EXACT_DIM=.true.)
                        
                        ! Compute at plasma solver grid
                        CALL EZspline_interp(Dn_spl,Nr_plasma_solver,rho_plasma_grid,Dn_NEO(jspecies,mytimestep_plasma_solver,:),ier)
                        CALL EZspline_interp(cn_spl,Nr_plasma_solver,rho_plasma_grid,cn_NEO(jspecies,mytimestep_plasma_solver,:),ier)
                        CALL EZspline_interp(Dp_spl,Nr_plasma_solver,rho_plasma_grid,Dp_NEO(jspecies,mytimestep_plasma_solver,:),ier)
                        CALL EZspline_interp(cp_spl,Nr_plasma_solver,rho_plasma_grid,cp_NEO(jspecies,mytimestep_plasma_solver,:),ier)

                        ! Deallocate splines
                        CALL EZspline_free(Dn_spl,ier)
                        CALL EZspline_free(cn_spl,ier)
                        CALL EZspline_free(Dp_spl,ier)
                        CALL EZspline_free(cp_spl,ier)
                  END DO
                  !
                  DEALLOCATE(Dn_temp,cn_temp,Dp_temp,cp_temp)
            END IF

            DEALLOCATE(GNEO_temp,QNEO_temp,rho_temp)
            DEALLOCATE(rho_k,iota,phip,chip,btheta,bzeta,bsq,vp,EparB)
            DEALLOCATE(te,ne,dtedrho,dnedrho)
            DEALLOCATE(ni,ti,dtidrho,dnidrho)
            DEALLOCATE(JBS_PENTA,etapar_PENTA,Er_PENTA,GNEO_PENTA,QNEO_PENTA)
            DEALLOCATE(Dn_PENTA,cn_PENTA,Dp_PENTA,cp_PENTA)

         ELSE !other threads
            DEALLOCATE(rho_k,iota,phip,chip,btheta,bzeta,bsq,vp,EparB)
            DEALLOCATE(te,ne,dtedrho,dnedrho)
            DEALLOCATE(ni,ti,dtidrho,dnidrho)
            RETURN
                        
         END IF

      ENDIF
      IF (lscreen) WRITE(6,'(a)') ' -------------------  NEOCLASSICAL BOOTSTRAP CALCULATION DONE  ---------------------'
      RETURN
!-----------------------------------------------------------------------
!     END SUBROUTINE
!-----------------------------------------------------------------------
      END SUBROUTINE thrift_penta