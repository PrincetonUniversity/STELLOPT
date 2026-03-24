!-----------------------------------------------------------------------
!     Subroutine:    stellopt_penta
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!                    A. Coelho (antonio.coelho@gauss-fusion.com)
!     Date:          03/23/2026
!     Description:   This subroutine runs PENTA for a set of surfaces
!                    and Er/nu pairs. The code assumes DKES has already
!                    been run.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_penta(lscreen,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime, ONLY:  proc_string, bigno, rprec
      USE equil_utils, ONLY: get_equil_phi, nrad, shat, phi_type
      USE stellopt_targets, ONLY: nu_dkes, lbooz, nsd, &
                                  sigma_dkes_11, sigma_dkes_31, &
                                  sigma_dkes_33, sigma_dkes_boot, &
                                  E_dkes, nprof, nruns_dkes, &
                                  sigma_dkes_erdiff, Ep_DKES_Erdiff, &
                                  Em_DKES_Erdiff, Ep_DKES_alpha, &
                                  Em_DKES_alpha, sigma_dkes_alpha, &
                                  nu_dkes_Erdiff, num_dkes_alpha, &
                                  nup_dkes_alpha
      USE booz_params, mboz_xboozer => mboz, nboz_xboozer => nboz,&
                       lscreen_xboozer => lscreen, nfp_xboozer => nfp, &
                       ns_xboozer => ns, lasym_xboozer => lasym_b, &
                       mpol_b => mpol, ntor_b => ntor, bcast_boozer_vars
      ! PENTA library
      USE PENTA_INTERFACE_MOD
      
!-----------------------------------------------------------------------
!     Subroutine Parameters
!        iflag         Error flag
!----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL, INTENT(in)    :: lscreen
      INTEGER, INTENT(inout) :: iflag
!-----------------------------------------------------------------------
!     Local Variables
!        ik          Helper integer (index over penta surfaces)
!        ii          Helper integer (index over all surfaces)
!----------------------------------------------------------------------
      INTEGER :: ik, ii, mystart, myend
      INTEGER :: nsurf_penta
      INTEGER, DIMENSION(:), ALLOCATABLE :: ik_penta
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: Vp_local, chip_local, phip_local, iota_local,
         btheta_local, bzeta_local
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (lscreen) WRITE(6,'(a)') ' ---------------------------    PENTA CALCULATION     -------------------------'
!DEC$ IF DEFINED (MPI_OPT)
!      ierr_mpi = 0
!      CALL MPI_BCAST(nruns_dkes,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
!DEC$ ENDIF
      ! We need to initialize a few things.
      CALL bcast_boozer_vars(master, MPI_COMM_MYWORLD, ierr_mpi)
      ! First we need all the equilibrium/profile information (master)
      ! Second we need to count how many target surfaces we're calculating so 
      IF (myworkid == master) THEN
         nsurf_penta = COUNT(nsurf_penta)
         ALLOCATE(ik_penta(nsurf_penta))
         DO ik = 1, nsurf_penta
            IF (lneed_penta(ik)) ik_pent(ik) = ik
         END DO
      END IF
!DEC$ IF DEFINED (MPI_OPT)
      ierr_mpi = 0
      CALL MPI_BCAST(nsurf_penta,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_BCAST(Aminor,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_BCAST(Rmajor,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
!DEC$ ENDIF    
      IF (myworkid /= master) THEN
         ALLOCATE(ik_penta(nsurf_penta))
      END IF
!DEC$ IF DEFINED (MPI_OPT)
      ierr_mpi = 0
      CALL MPI_BCAST(ik_penta,nsurf_penta,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
!DEC$ ENDIF   
      CALL MPI_CALC_MYRANGE(MPI_COMM_MYWORLD,1,nsurf_penta,mystart,myend)
      ! Loop over radial surfaces
      DO ik = mystart,myend
         ! PENTA
         ! Not needed beacause we read indata namelist.
         !CALL PENTA_SET_ION_PARAMS(nion_prof, DBLE(Zatom_local), Matom_local)
         CALL PENTA_SET_COMMANDLINE(Er_min_Vcm,Er_max_Vcm,ik_penta(ik),1,EparB(k),1,'','','')
         CALL PENTA_ALLOCATE_SPECIES
         ! I'm passing actual rho here, so if you need s then use shat
         ii = ik_penta(ik)
         CALL PENTA_SET_EQ_DATA(rho(ii),Aminor,Rmajor,vp(k),chip(ii),phip_b(k),iota_b(ii),buco_b(ii),bvco_b(ii),bsq(k))
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
         CALL PENTA_FIT_RAD_TRANS

         WRITE(temp_str,'(i4.4)') k
         WRITE(temp1_str,'(i3.3)') mytimestep
         CALL PENTA_OPEN_OUTPUT(TRIM(temp1_str) // '_k' // TRIM(temp_str))

         ! Only search new ambipolar Er solution if look_for_ambipolar = true
         ! If not, take the current THRIFT_ER solution
         ! The boolean look_for_ambipolar is controled in thrift_plasma_solver_mod
         IF(solve_plasma_equations .AND. .NOT. look_for_ambipolar) THEN
               ! Get Er value from THRIFT
               irho =  MINLOC(ABS(SQRT(THRIFT_S) - rho_k(k)), dim=1)
               ! Need to define num_roots and set Er_roots
               num_roots = 1
               Er_roots(1) = THRIFT_ER(irho,mytimestep)
         ELSE
               ! Now the basic steps
               CALL PENTA_RUN_2_EFIELD
               CALL PENTA_RUN_3_FIND_ROOTS
         END IF

         CALL PENTA_RUN_4_AMBIPOLAR
         
         ! The call to ROOT_ANALYSIS sets the array 'root_type' which decides which root to pick
         ! The criterium is to pick the 'ion_root'
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
         !       
         CALL PENTA_RUN_5_CLEANUP(lscreen)
      END DO
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_MYWORLD,ierr_mpi)
      IF (myworkid == master) THEN
         CALL MPI_REDUCE(MPI_IN_PLACE, JBS_PENTA, nsurf_penta, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
         CALL MPI_REDUCE(MPI_IN_PLACE, etapar_PENTA, nsurf_penta, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
         CALL MPI_REDUCE(MPI_IN_PLACE, Er_PENTA, nsurf_penta, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
      ELSE
         CALL MPI_REDUCE(JBS_PENTA,    JBS_PENTA, nsurf_penta, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
         CALL MPI_REDUCE(etapar_PENTA, etapar_PENTA, nsurf_penta, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
         CALL MPI_REDUCE(Er_PENTA,     Er_PENTA, nsurf_penta, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
      END IF
!DEC$ ENDIF
      IF (lscreen) WRITE(6,'(a)') ' -----------------  PENTA CALCULATION (DONE) ----------------'
      !IF (myworkid .ne. master) CALL read_wout_deallocate
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_penta
