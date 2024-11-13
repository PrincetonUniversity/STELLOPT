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
                 mysurf, root_max_Er
      REAL(rprec) :: s, rho
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: rho_k, iota, phip, chip, btheta, bzeta, bsq, vp, &
                        te, ne, dtedrho, dnedrho, EparB, JBS_PENTA, etapar_PENTA, rho_temp, J_temp, eta_temp
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: ni,ti, dtidrho, dnidrho
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: D11, D13, D33
      TYPE(EZspline1_r8) :: J_spl, eta_spl
      INTEGER :: bcs0(2)
!-----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!-----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (lscreen) WRITE(6,'(a)') ' --------------------  NEOCLASSICAL BOOTSTRAP USING PENTA  -------------------'
      IF (lvmec) THEN
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
         ALLOCATE(JBS_PENTA(ns_dkes),etapar_PENTA(ns_dkes))

         JBS_PENTA = 0.0; etapar_PENTA = 0.0

         IF (myworkid == master) THEN
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

                  !EparB
                  EparB(k) = 0.0_rprec !! temporary... cannot do THRIFT_EPARB(mysurf,mytimestep) since this is defined in thrift grid, not vmec... interpolation?
            END DO


            ! print *, 'master_thrift_penta: dkes_d11', DKES_D11
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
         
#endif
      
         ! DKES Data
         ncstar = COUNT(DKES_NUSTAR < 1E10)
         nestar = COUNT(DKES_ERSTAR < 1E10)

         DO k = mystart,myend
            ! PENTA
            CALL PENTA_SET_ION_PARAMS(nion_prof, DBLE(Zatom_prof), Matom_prof/p_mass)

            !!! WE NEED TO CHANGE THE EPARB --- it needs to interact with THRIFT...!!!
            !!! ALSO: IS THIS -250:250 OK?? I PROPOSE TO MOVE THIS TO THE INPUT...
            CALL PENTA_SET_COMMANDLINE(-250.0_rprec,250.0_rprec,DKES_K(k),1,EparB(k),1,'','','')
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
            CALL PENTA_OPEN_OUTPUT(proc_string)
            CALL PENTA_FIT_RAD_TRANS
            ! Now the basic steps
            CALL PENTA_RUN_2_EFIELD
            CALL PENTA_RUN_3_AMBIPOLAR

            ! Save JBS corresponding to the root that has the largest Er
            ! This because whenever there are 2 stable roots, a rule of thumb is to pick the one with largest Er
            root_max_Er = MAXLOC(Er_roots(1:num_roots),1)
            JBS_PENTA(k) = J_BS_ambi(root_max_Er)
            etapar_PENTA(k) = sigma_par_ambi(root_max_Er)

            CALL PENTA_RUN_4_CLEANUP
         END DO

         !! Bootstrap interpolation onto THRIFT grid
#if defined(MPI_OPT)
         CALL MPI_BARRIER(MPI_COMM_MYWORLD,ierr_mpi)
         IF (myworkid == master) THEN
            CALL MPI_REDUCE(MPI_IN_PLACE,JBS_PENTA,ns_dkes,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_REDUCE(MPI_IN_PLACE,etapar_PENTA,ns_dkes,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
         ELSE 
            CALL MPI_REDUCE(JBS_PENTA,JBS_PENTA,ns_dkes,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_REDUCE(etapar_PENTA,etapar_PENTA,ns_dkes,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
            DEALLOCATE(rho_k,iota,phip,chip,btheta,bzeta,bsq,vp,EparB)
            DEALLOCATE(te,ne,dtedrho,dnedrho)
            DEALLOCATE(ni,ti,dtidrho,dnidrho)
            DEALLOCATE(JBS_PENTA,etapar_PENTA)
            RETURN
         ENDIF
#endif
         
         IF (myworkid == master) THEN

            ! Interpolate JBS_PENTA and etapar_PENTA at rho=0 and rho=1
            ALLOCATE(J_temp(ns_dkes+2),eta_temp(ns_dkes+2),rho_temp(ns_dkes+2))
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

            ! Splines
            bcs0=(/ 0, 0/)
            !JBS
            CALL EZspline_init(J_spl,ns_dkes+2,bcs0,ier)
            J_spl%x1        = rho_temp
            J_spl%isHermite = 1
            CALL EZspline_setup(J_spl,J_temp,ier,EXACT_DIM=.true.)
            !etapar
            CALL EZspline_init(eta_spl,ns_dkes+2,bcs0,ier)
            eta_spl%x1        = rho_temp
            eta_spl%isHermite = 1
            CALL EZspline_setup(eta_spl,eta_temp,ier,EXACT_DIM=.true.)
            !
            DEALLOCATE(J_temp,eta_temp,rho_temp)

            ! Calculate J_BS and etapara in THRFIT GRID
            DO i = 1, nsj
                  rho = SQRT( THRIFT_S(i) )
                  CALL EZspline_interp(J_spl,rho,THRIFT_JBOOT(i,mytimestep),ier)
                  IF( etapar_type == 'dkespenta') CALL EZspline_interp(eta_spl,rho,THRIFT_ETAPARA(i,mytimestep),ier)
            END DO
            CALL EZspline_free(J_spl,ier)
            CALL EZspline_free(eta_spl,ier)

            ! PRINT *, 'THRIFT_JBOOT=', THRIFT_JBOOT(:,mytimestep)
            ! PRINT *, 'THRIFT_ETAPARA=', THRIFT_ETAPARA(:,mytimestep)

            DEALLOCATE(rho_k,iota,phip,chip,btheta,bzeta,bsq,vp,EparB)
            DEALLOCATE(te,ne,dtedrho,dnedrho)
            DEALLOCATE(ni,ti,dtidrho,dnidrho)
            DEALLOCATE(JBS_PENTA,etapar_PENTA)
            
         END IF

      ENDIF
      IF (lscreen) WRITE(6,'(a)') ' -------------------  NEOCLASSICAL BOOTSTRAP CALCULATION DONE  ---------------------'
      RETURN
!-----------------------------------------------------------------------
!     END SUBROUTINE
!-----------------------------------------------------------------------
      END SUBROUTINE thrift_penta