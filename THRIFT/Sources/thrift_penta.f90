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
      G_NEO_complet,Q_NEO_complet,Nt_total_plasma_solver,plasma_Er,Er_spline
      USE thrift_globals, ONLY: look_for_ambipolar,update_thrift_vars,update_transport_vars
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
                 mysurf, root_max_Er, jspecies, irho, it_prev
      REAL(rprec) :: s, rho, mytime
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: rho_k, iota, phip, chip, btheta, bzeta, bsq, vp, &
                        te, ne, dtedrho, dnedrho, EparB, JBS_PENTA, etapar_PENTA, Er_PENTA, rho_temp, J_temp, eta_temp, Er_temp, &
                        Er_k
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: GNEO_PENTA, QNEO_PENTA, GNEO_temp, QNEO_temp
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: Dn_PENTA, cn_PENTA, Dn_temp, cn_temp
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: Dp_PENTA, cp_PENTA, Dp_temp, cp_temp
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: ni,ti, dtidrho, dnidrho
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: D11, D13, D33
      TYPE(EZspline1_r8) :: EparB_spl, J_spl, eta_spl, Er_spl, GNEO_spl, QNEO_spl
      TYPE(EZspline1_r8) :: Dn_spl, cn_spl, Dp_spl, cp_spl
      INTEGER :: bcs0(2)=(/ 0, 0/)
      CHARACTER(LEN=32) :: temp_str, temp1_str
!-----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!-----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (lscreen) WRITE(6,'(a)') ' --------------------  NEOCLASSICAL BOOTSTRAP USING PENTA  -------------------'
      IF (lscreen) Write(*,*) " <r>/<a>","   Er root(s) (V/cm)"

      IF (.NOT. lvmec) RETURN
      

      ierr_mpi = 0
      ! PENTA is parallelized over radial surfaces in this routine.
      ns_dkes = 0
      DO k = 1, DKES_NS_MAX
      IF ((DKES_K(k) > 0)) ns_dkes = ns_dkes+1
      END DO
      ! Break up work
      CALL MPI_CALC_MYRANGE(MPI_COMM_MYWORLD,1,ns_dkes,mystart,myend)

      ALLOCATE(rho_k(ns_dkes),iota(ns_dkes),phip(ns_dkes),chip(ns_dkes),btheta(ns_dkes),bzeta(ns_dkes),bsq(ns_dkes),vp(ns_dkes),EparB(ns_dkes))
      ALLOCATE(te(ns_dkes),ne(ns_dkes),dtedrho(ns_dkes),dnedrho(ns_dkes))
      ALLOCATE(ni(ns_dkes,nion_prof),ti(ns_dkes,nion_prof),dtidrho(ns_dkes,nion_prof),dnidrho(ns_dkes,nion_prof))
      ALLOCATE(JBS_PENTA(ns_dkes),etapar_PENTA(ns_dkes),Er_PENTA(ns_dkes),Er_k(ns_dkes))
      ALLOCATE(GNEO_PENTA(nion_prof+1,ns_dkes),QNEO_PENTA(nion_prof+1,ns_dkes))
      ALLOCATE(Dn_PENTA(nion_prof+1,ns_dkes),cn_PENTA(nion_prof+1,ns_dkes))
      ALLOCATE(Dp_PENTA(nion_prof+1,ns_dkes),cp_PENTA(nion_prof+1,ns_dkes))

      JBS_PENTA = 0.0; etapar_PENTA = 0.0; Er_PENTA = 0.0; Er_k = 0.0
      GNEO_PENTA = 0.0; QNEO_PENTA = 0.0
      Dn_PENTA = 0.0; cn_PENTA = 0.0
      Dp_PENTA = 0.0; cp_PENTA = 0.0

      IF (myworkid == master) THEN
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
                  ! Er from previous plasma solver step
                  IF (solve_plasma_equations) THEN
                        ier = 0
                        CALL EZspline_interp(Er_spline, rho, Er_k(k), ier)
                  END IF
            END DO
            !
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
      CALL MPI_BCAST(Er_k,ns_dkes,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
      ! VMEC quantities
      CALL MPI_BCAST(eq_Aminor,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_BCAST(eq_Rmajor,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
      ! THRIFT quantities
      CALL MPI_BCAST(mytime,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_BCAST(mytimestep,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
      ! Plasma solver quantities
      CALL MPI_BCAST(solve_plasma_equations,1,MPI_LOGICAL,master,MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_BCAST(mytimestep_plasma_solver,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)       
      ! thrift_globals
      CALL MPI_BCAST(look_for_ambipolar,1,MPI_LOGICAL,master,MPI_COMM_MYWORLD,ierr_mpi)
#endif

      ! DKES Data
      ncstar = COUNT(DKES_NUSTAR < 1E10)
      nestar = COUNT(DKES_ERSTAR < 1E10)

      DO k = mystart,myend
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
            CALL PENTA_LMAT_MATRIX
            CALL PENTA_SET_INTEGRATION_ARRAYS

            WRITE(temp_str,'(i4.4)') k
            WRITE(temp1_str,'(i3.3)') mytimestep
            IF(look_for_ambipolar) CALL PENTA_OPEN_OUTPUT(TRIM(temp1_str) // '_k' // TRIM(temp_str))

            ! Only search new ambipolar Er solution if look_for_ambipolar = true
            ! If not, take the previous plasma_Er (stored in Er_spline)
            IF(solve_plasma_equations .AND. .NOT. look_for_ambipolar) THEN
                  ! Need to define num_roots and set Er_roots
                  num_roots = 1
                  Er_roots(1) = Er_k(k)
            ELSE
                  ! Now the basic steps
                  CALL PENTA_RUN_2_EFIELD
                  CALL PENTA_RUN_3_FIND_ROOTS
            END IF

            CALL PENTA_RUN_4_AMBIPOLAR
            
            ! The call to ROOT_ANALYSIS sets the array 'root_type' which decides which root to pick
            ! The criterium is to pick the 'ion_root'
            CALL ROOT_ANALYSIS(TRIM(Er_root_type))

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
            !       
            CALL PENTA_RUN_5_CLEANUP(lscreen)
      END DO

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_MYWORLD,ierr_mpi)
      IF (myworkid == master) THEN
            CALL MPI_REDUCE(MPI_IN_PLACE,JBS_PENTA,ns_dkes,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_REDUCE(MPI_IN_PLACE,etapar_PENTA,ns_dkes,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_REDUCE(MPI_IN_PLACE,Er_PENTA,ns_dkes,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_REDUCE(MPI_IN_PLACE,GNEO_PENTA,(nion_prof+1)*ns_dkes,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_REDUCE(MPI_IN_PLACE,QNEO_PENTA,(nion_prof+1)*ns_dkes,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_REDUCE(MPI_IN_PLACE,Dn_PENTA,(nion_prof+1)*ns_dkes,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_REDUCE(MPI_IN_PLACE,cn_PENTA,(nion_prof+1)*ns_dkes,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_REDUCE(MPI_IN_PLACE,Dp_PENTA,(nion_prof+1)*ns_dkes,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_REDUCE(MPI_IN_PLACE,cp_PENTA,(nion_prof+1)*ns_dkes,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
      ELSE 
            CALL MPI_REDUCE(JBS_PENTA,JBS_PENTA,ns_dkes,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_REDUCE(etapar_PENTA,etapar_PENTA,ns_dkes,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_REDUCE(Er_PENTA,Er_PENTA,ns_dkes,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_REDUCE(GNEO_PENTA,GNEO_PENTA,(nion_prof+1)*ns_dkes,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_REDUCE(QNEO_PENTA,QNEO_PENTA,(nion_prof+1)*ns_dkes,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_REDUCE(Dn_PENTA,Dn_PENTA,(nion_prof+1)*ns_dkes,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_REDUCE(cn_PENTA,cn_PENTA,(nion_prof+1)*ns_dkes,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_REDUCE(Dp_PENTA,Dp_PENTA,(nion_prof+1)*ns_dkes,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_REDUCE(cp_PENTA,cp_PENTA,(nion_prof+1)*ns_dkes,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL FLUSH(6)
            DEALLOCATE(rho_k,iota,phip,chip,btheta,bzeta,bsq,vp,EparB)
            DEALLOCATE(te,ne,dtedrho,dnedrho)
            DEALLOCATE(ni,ti,dtidrho,dnidrho)
            DEALLOCATE(JBS_PENTA,etapar_PENTA,Er_PENTA,Er_k,GNEO_PENTA,QNEO_PENTA,Dn_PENTA,cn_PENTA,Dp_PENTA,cp_PENTA)
            RETURN
      ENDIF
#endif
            
      IF (myworkid == master) THEN
            IF(look_for_ambipolar) THEN
                  IF(save_all_ambipolar_roots) CALL PENTA_MERGE_AMBIPOLAR_FILES(ns_dkes,temp1_str,mytime)
                  IF(save_fluxes_vs_Er) CALL PENTA_MERGE_FLUXES_VS_ER_FILES(ns_dkes,temp1_str,mytime)
            END IF

            IF(solve_plasma_equations) THEN
                  ! plasma_Er updates only when look_for_ambipolar is True
                  ! otherwise, can simply store that of previous iter
                  IF(look_for_ambipolar) THEN
                        CALL interpolate_from_PENTA(nrho_penta=ns_dkes, rho_penta=rho_k, y_penta=Er_PENTA,\
                                          nrho_out=Nr_plasma_solver, rho_out=rho_plasma_grid, y_out=plasma_Er(mytimestep_plasma_solver,:),\
                                          isHermite=1, useLog=.FALSE.)
                  ELSE
                        plasma_Er(mytimestep_plasma_solver,:) = plasma_Er(mytimestep_plasma_solver-1,:)
                  END IF
            END IF

            IF(update_thrift_vars) THEN
                  CALL interpolate_from_PENTA(nrho_penta=ns_dkes, rho_penta=rho_k, y_penta=JBS_PENTA,\
                                          nrho_out=nsj, rho_out=SQRT(THRIFT_S), y_out=THRIFT_JBOOT(:,mytimestep),\
                                          isHermite=1, useLog=.FALSE.)
                  !
                  CALL interpolate_from_PENTA(nrho_penta=ns_dkes, rho_penta=rho_k, y_penta=Er_PENTA,\
                                          nrho_out=nsj, rho_out=SQRT(THRIFT_S), y_out=THRIFT_ER(:,mytimestep),\
                                          isHermite=1, useLog=.FALSE.)
                  !
                  IF(etapar_type == 'dkespenta') THEN
                  CALL interpolate_from_PENTA(nrho_penta=ns_dkes, rho_penta=rho_k, y_penta=etapar_PENTA,\
                                          nrho_out=nsj, rho_out=SQRT(THRIFT_S), y_out=THRIFT_ETAPARA(:,mytimestep),\
                                          isHermite=0, useLog=.TRUE.) 
                  END IF
                  !
                  DO jspecies=1,(nion_prof+1)
                        CALL interpolate_from_PENTA(nrho_penta=ns_dkes, rho_penta=rho_k, y_penta=GNEO_PENTA(jspecies,:),\
                                          nrho_out=nsj, rho_out=SQRT(THRIFT_S), y_out=THRIFT_GNEO(jspecies,:,mytimestep),\
                                          isHermite=1, useLog=.FALSE.)
                        !
                        CALL interpolate_from_PENTA(nrho_penta=ns_dkes, rho_penta=rho_k, y_penta=QNEO_PENTA(jspecies,:),\
                                          nrho_out=nsj, rho_out=SQRT(THRIFT_S), y_out=THRIFT_QNEO(jspecies,:,mytimestep),\
                                          isHermite=1, useLog=.FALSE.)
                  END DO
            END IF
            
            IF(solve_plasma_equations .AND. update_transport_vars) THEN
                  DO jspecies=1,(nion_prof+1)
                        !
                        CALL interpolate_from_PENTA(nrho_penta=ns_dkes, rho_penta=rho_k, y_penta=GNEO_PENTA(jspecies,:),\
                                          nrho_out=Nr_plasma_solver, rho_out=rho_plasma_grid, y_out=G_NEO_complet(jspecies,mytimestep_plasma_solver,:),\
                                          isHermite=1, useLog=.FALSE.)
                        !
                        CALL interpolate_from_PENTA(nrho_penta=ns_dkes, rho_penta=rho_k, y_penta=QNEO_PENTA(jspecies,:),\
                                          nrho_out=Nr_plasma_solver, rho_out=rho_plasma_grid, y_out=Q_NEO_complet(jspecies,mytimestep_plasma_solver,:),\
                                          isHermite=1, useLog=.FALSE.)
                        !
                        CALL interpolate_from_PENTA(nrho_penta=ns_dkes, rho_penta=rho_k, y_penta=Dn_PENTA(jspecies,:),\
                                          nrho_out=Nr_plasma_solver, rho_out=rho_plasma_grid, y_out=Dn_NEO(jspecies,mytimestep_plasma_solver,:),\
                                          isHermite=1, useLog=.FALSE.)
                        !
                        CALL interpolate_from_PENTA(nrho_penta=ns_dkes, rho_penta=rho_k, y_penta=cn_PENTA(jspecies,:),\
                                          nrho_out=Nr_plasma_solver, rho_out=rho_plasma_grid, y_out=cn_NEO(jspecies,mytimestep_plasma_solver,:),\
                                          isHermite=1, useLog=.FALSE.)
                        !
                        CALL interpolate_from_PENTA(nrho_penta=ns_dkes, rho_penta=rho_k, y_penta=Dp_PENTA(jspecies,:),\
                                          nrho_out=Nr_plasma_solver, rho_out=rho_plasma_grid, y_out=Dp_NEO(jspecies,mytimestep_plasma_solver,:),\
                                          isHermite=1, useLog=.FALSE., preventNeg = .TRUE.)
                        !
                        CALL interpolate_from_PENTA(nrho_penta=ns_dkes, rho_penta=rho_k, y_penta=cp_PENTA(jspecies,:),\
                                          nrho_out=Nr_plasma_solver, rho_out=rho_plasma_grid, y_out=cp_NEO(jspecies,mytimestep_plasma_solver,:),\
                                          isHermite=1, useLog=.FALSE.)
                  END DO
            END IF
            DEALLOCATE(rho_k,iota,phip,chip,btheta,bzeta,bsq,vp,EparB)
            DEALLOCATE(te,ne,dtedrho,dnedrho)
            DEALLOCATE(ni,ti,dtidrho,dnidrho)
            DEALLOCATE(JBS_PENTA,etapar_PENTA,Er_PENTA,Er_k,GNEO_PENTA,QNEO_PENTA)
            DEALLOCATE(Dn_PENTA,cn_PENTA,Dp_PENTA,cp_PENTA)              
      END IF

      ! Restore booleans
      look_for_ambipolar = .FALSE.
      update_thrift_vars = .FALSE.
      update_transport_vars = .FALSE.

      IF (lscreen) WRITE(6,'(a)') ' -------------------  NEOCLASSICAL BOOTSTRAP CALCULATION DONE  ---------------------'
      RETURN
      !-----------------------------------------------------------------------
      !     END SUBROUTINE
      !-----------------------------------------------------------------------
      END SUBROUTINE thrift_penta