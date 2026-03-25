!-----------------------------------------------------------------------
!   Subroutine:    stellopt_penta
!   Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!                  A. Coelho (antonio.coelho@gauss-fusion.com)
!   Date:          03/23/2026
!   Description:   This subroutine runs PENTA for a set of surfaces
!                  and Er/nu pairs. The code assumes DKES has already
!                  been run.
!-----------------------------------------------------------------------
   SUBROUTINE stellopt_penta(lscreen,iflag)
!-----------------------------------------------------------------------
!       Libraries
!-----------------------------------------------------------------------
   USE stellopt_runtime, ONLY: proc_string, bigno, rprec
   USE stellopt_targets, ONLY: lneed_penta, lneed_dkes, nsd, nprof, &
      nu_dkes, E_dkes
   USE stellopt_vars, ONLY: ne_type, te_type, ti_type
   USE equil_utils, ONLY: get_equil_te, get_equil_ti, get_equil_ne, &
      get_equil_volume, get_equil_bdotb
   USE equil_vals, ONLY: Aminor, Rmajor, rho, shat, Er_PENTA, JBS_PENTA
   USE read_boozer_mod, ONLY: bcast_boozer_vars, phip_b, iota_b, buco_b, bvco_b
!DEC$ IF DEFINED (DKES_OPT)
   USE dkes_realspace, ONLY: DKES_L11m, DKES_L11p, DKES_L31m, &
      DKES_L31p, DKES_L33m, DKES_L33p
!DEC$ ENDIF
   ! PENTA library
   USE PENTA_INTERFACE_MOD
   USE mpi_params
   USE mpi_inc
      
!-----------------------------------------------------------------------
!  Subroutine Parameters
!     iflag         Error flag
!-----------------------------------------------------------------------
   IMPLICIT NONE
   LOGICAL, INTENT(in)    :: lscreen
   INTEGER, INTENT(inout) :: iflag

!-----------------------------------------------------------------------
!  Local Variables
!     ii,ij,ik,il       Helper index
!     ier               Error flag
!     mystart, myend    Helpers for parallelization over arrays
!     nsurf_penta       Number of PENTA surfaces
!     nion_prof         Number of ions
!     ncstar            Number of unique collisionalities
!     nestar            Number of unique Er
!     s_local           Normalized Toroidal Flux helper
!     s2_local          Normalized Toroidal Flux helper (for deriv.)
!     rho_local         Rho helper
!     dprof             Profile helper (derivative)
!     EparB             E.B ?????
!     ik_penta          Array of surface index values
!     XX_local          Arrays of size nsurf_penta to help run
!     temp_str          String helper for file names
!-----------------------------------------------------------------------
   INTEGER :: ii, ij, ik, il, im, ier, mystart, myend
   INTEGER :: nsurf_penta, nion_prof, ncstar, nestar, nsurf_dkes
   REAL(rprec) :: s_local, s2_local, rho_local, dprof, EparB, Er, Nu
   INTEGER, DIMENSION(:), ALLOCATABLE :: ik_penta
   REAL(rprec), DIMENSION(:), ALLOCATABLE :: te_local, ne_local, &
                              dtedrho_local, dnedrho_local, &
                              vp_local, bdotb_local, &
                              dkes_nustar, dkes_erstar, &
                              s_penta
   REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: ti_local, ni_local, &
                              dtidrho_local, dnidrho_local
   REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: DKES_D11, DKES_D31, DKES_D33

   CHARACTER(LEN=256) :: temp_str 

!-----------------------------------------------------------------------
!  PARAMETERS
!-----------------------------------------------------------------------
   LOGICAL, PARAMETER :: solve_plasma_equations = .FALSE.
   LOGICAL, PARAMETER :: look_for_ambipolar = .FALSE.

!-----------------------------------------------------------------------
!  BEGIN SUBROUTINE
!-----------------------------------------------------------------------
   IF (iflag < 0) RETURN
   IF (lscreen) WRITE(6,'(a)') ' ---------------------------    PENTA CALCULATION     -------------------------'
   ! This make sure everyone has boozer data
   CALL bcast_boozer_vars(master, MPI_COMM_MYWORLD, ierr_mpi)
   ! Master is the only thread who knows things
   ! We setup a bunch of helper arrays using lookup functions
   IF (myworkid == master) THEN
      nsurf_penta = COUNT(lneed_penta)
      nion_prof = num_ion_species
      ALLOCATE(ik_penta(nsurf_penta), s_penta(nsurf_penta), &
         ne_local(nsurf_penta), te_local(nsurf_penta), &
         ni_local(nsurf_penta,nion_prof), ti_local(nsurf_penta,nion_prof))
      ALLOCATE(dnedrho_local(nsurf_penta), dtedrho_local(nsurf_penta),&
         dnidrho_local(nsurf_penta,nion_prof), dtidrho_local(nsurf_penta,nion_prof))
      ALLOCATE(vp_local(nsurf_penta),bdotb_local(nsurf_penta))
      ii = 1
      DO ik = 1, nsd
         IF (lneed_penta(ik)) THEN
             ik_penta(ii) = ik
             s_penta(ii) = shat(ik)
             s_local = shat(ik)
             s2_local = shat(ik-1)
             rho_local = rho(ik)
             CALL get_equil_ne(s_local,TRIM(ne_type),ne_local(ii),ier)
             CALL get_equil_ne(s2_local,TRIM(ne_type),dprof,ier)
             dnedrho_local(ii) = 2.0*rho_local*(ne_local(ii)-dprof)/(s_local-s2_local)
             CALL get_equil_te(s_local,TRIM(te_type),te_local(ii),ier)
             CALL get_equil_te(s2_local,TRIM(te_type),dprof,ier)
             dtedrho_local(ii) = 2.0*rho_local*(te_local(ii)-dprof)/(s_local-s2_local)
             ! Ions 
             ni_local(ii,:) = ne_local(ii)/nion_prof ! Assume equal for now
             dnidrho_local(ii,:) = dnedrho_local(ii)/nion_prof ! Assume equal for now
             CALL get_equil_ti(s_local,TRIM(ti_type),ti_local(ii,1),ier)
             CALL get_equil_ti(s2_local,TRIM(ti_type),dprof,ier)
             dtidrho_local(ii,1) = 2.0*rho_local*(ti_local(ii,1)-dprof)/(s_local-s2_local)
             ti_local(ii,:) = ti_local(ii,1)
             dtidrho_local(ii,:) = dtidrho_local(ii,1)
             ! VP
             CALL get_equil_volume(s_local,dprof,ier,vp_local(ii))
             CALL get_equil_bdotb(s_local,bdotb_local(ii),ier)
             ii = ii +1
         END IF
      END DO
   END IF
   !!!!!!!!!!!!!!!!!!!!!!!!Sorting of NU_DKES and ER_DKES!!!!!!!!!!!!
   !!  This assumes that the arrays are in ascending order and ER 
   !!  varies faster than NU
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! First count unique values
   ncstar = 0; nestar = 0
   nu = -bigno; Er = -bigno
   DO ij = 1, nprof
      IF (nu_dkes(ij) > nu) THEN
         ncstar = ncstar + 1
         nu = nu_dkes(ij)
      END IF
      IF (E_dkes(ij)  > Er) THEN
         nestar = nestar + 1
         Er = E_dkes(ij)
      END IF
   END DO
   ! Now setup helper arrays
   ALLOCATE(dkes_nustar(ncstar), dkes_erstar(nestar))
   ! Now we send data to the rest of the threads
   ncstar = 0; nestar = 0
   nu = -bigno; Er = -bigno
   DO ij = 1, nprof
      IF (nu_dkes(ij) > nu) THEN
         ncstar = ncstar + 1
         nu = nu_dkes(ij)
         dkes_nustar(ncstar) = nu_dkes(ij)
      END IF
      IF (E_dkes(ij)  > Er) THEN
         nestar = nestar + 1
         Er = E_dkes(ij)
         dkes_erstar(nestar) = E_dkes(ij)
      END IF
   END DO
   ! Now do the DKES Coefficient arrays
   ! Note there can be more DKES surfaces than PENTA
   ALLOCATE(DKES_D11(nsurf_penta,ncstar,nestar),DKES_D31(nsurf_penta,ncstar,nestar),DKES_D33(nsurf_penta,ncstar,nestar))
   il = 1; im = 0
   nsurf_dkes = COUNT(lneed_dkes)
   DO ik = 1, nsurf_dkes
      IF (lneed_penta(ik)) im = im + 1
      DO ii = 1, ncstar
         DO ij = 1, nestar
            IF (lneed_penta(ik)) THEN
               DKES_D11(im,ii,ij) = DKES_L11p(il) + DKES_L11m(il)
               DKES_D31(im,ii,ij) = DKES_L31p(il) + DKES_L31m(il)
               DKES_D33(im,ii,ij) = DKES_L33p(il) + DKES_L33m(il)
            END IF
            il = il + 1
         END DO
      END DO
   END DO
   DKES_D11 = DKES_D11 * 0.5
   DKES_D31 = DKES_D31 * 0.5
   DKES_D33 = DKES_D33 * 0.5
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!DEC$ IF DEFINED (MPI_OPT)
   ierr_mpi = 0
   CALL MPI_BCAST(nsurf_penta,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
   CALL MPI_BCAST(nion_prof,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
   CALL MPI_BCAST(Aminor,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
   CALL MPI_BCAST(Rmajor,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
!DEC$ ENDIF
   ! Other threads allocate helpers  
   IF (myworkid /= master) THEN
      ALLOCATE(ik_penta(nsurf_penta), s_penta(nsurf_penta), &
         ne_local(nsurf_penta), te_local(nsurf_penta), &
         ni_local(nsurf_penta,nion_prof), ti_local(nsurf_penta,nion_prof))
      ALLOCATE(dnedrho_local(nsurf_penta), dtedrho_local(nsurf_penta),&
         dnidrho_local(nsurf_penta,nion_prof), dtidrho_local(nsurf_penta,nion_prof))
      ALLOCATE(vp_local(nsurf_penta),bdotb_local(nsurf_penta))
   END IF
   ! Now we broadcast the helpers to all threads
!DEC$ IF DEFINED (MPI_OPT)
   ierr_mpi = 0
   CALL MPI_BCAST(ik_penta,nsurf_penta,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
   CALL MPI_BCAST(s_penta,nsurf_penta,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
   CALL MPI_BCAST(vp_local,nsurf_penta,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
   CALL MPI_BCAST(bdotb_local,nsurf_penta,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
   CALL MPI_BCAST(ne_local,nsurf_penta,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
   CALL MPI_BCAST(te_local,nsurf_penta,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
   CALL MPI_BCAST(ti_local,nsurf_penta*nion_prof,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
   CALL MPI_BCAST(ni_local,nsurf_penta*nion_prof,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
   CALL MPI_BCAST(dnedrho_local,nsurf_penta,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
   CALL MPI_BCAST(dtedrho_local,nsurf_penta,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
   CALL MPI_BCAST(dtidrho_local,nsurf_penta*nion_prof,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
   CALL MPI_BCAST(dnidrho_local,nsurf_penta*nion_prof,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
   CALL MPI_BCAST(DKES_D11,nsurf_penta*ncstar*nestar,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
   CALL MPI_BCAST(DKES_D31,nsurf_penta*ncstar*nestar,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
   CALL MPI_BCAST(DKES_D33,nsurf_penta*ncstar*nestar,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
!DEC$ ENDIF   
   ! Everyone allocates the JBS and ER arrays
   IF (ALLOCATED(JBS_PENTA)) DEALLOCATE(JBS_PENTA)
   IF (ALLOCATED(Er_PENTA)) DEALLOCATE(Er_PENTA)
   ALLOCATE(JBS_PENTA(nsurf_penta), Er_PENTA(nsurf_penta))
   JBS_PENTA = 0.0; Er_PENTA = 0.0
   ! Break up the work
   CALL MPI_CALC_MYRANGE(MPI_COMM_MYWORLD,1,nsurf_penta,mystart,myend)
   ! Loop over radial surfaces
   DO ik = mystart,myend
      ! ii is the index in VMEC/Boozer grid, ik is over the PENTA surfaces
      ii = ik_penta(ik)
      s_local = s_penta(ik)
      rho_local = SQRT(s_local)
      ! Not needed beacause we read indata namelist (in chisq_penta_er, everyone does this)
      !CALL PENTA_SET_ION_PARAMS(nion_prof, DBLE(Zatom_local), Matom_local)
      EparB = 0.0 ! Ummm should this be zero for steady state?
      CALL PENTA_SET_COMMANDLINE(Er_min_Vcm,Er_max_Vcm,ii,1,EparB,1,'','','')
      CALL PENTA_ALLOCATE_SPECIES
      ! I'm passing actual rho here, so if you need s then use s_local
      CALL PENTA_SET_EQ_DATA(rho_local, Aminor, Rmajor, &
                         vp_local(ik), phip_b(ii)*iota_b(ii),&
                         phip_b(ii), iota_b(ii), &
                         buco_b(ii), bvco_b(ii), bdotb_local(ik))
      CALL PENTA_SET_PPROF(ne_local(ik),   dnedrho_local(ik)/Aminor,&
                           te_local(ik),   dtedrho_local(ik)/Aminor,&
                           ni_local(ik,:), dnidrho_local(ik,:)/Aminor,&
                           ti_local(ik,:), dtidrho_local(ik,:)/Aminor)
      ! MAKE CORRECTIONS ON D31 AND D33 -- values coming from DKES2 miss Bsq factors (see J. Lore documentation)
      CALL PENTA_SET_DKES_STAR(ncstar, nestar, DKES_NUSTAR(1:ncstar), DKES_ERSTAR(1:nestar), &
            DKES_D11(ik,:,:), DKES_D31(ik,:,:)*SQRT(bdotb_local(ik)), DKES_D33(ik,:,:)*bdotb_local(ik))
      CALL PENTA_SET_BEAM(0.0_rprec) ! Zero becasue we don't read
      CALL PENTA_SET_U2() ! Leave blank for default value
      CALL PENTA_READ_INPUT_FILES(.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.)
      CALL PENTA_SCREEN_INFO
      CALL PENTA_ALLOCATE_DKESCOEFF
      CALL PENTA_FIT_DXX_COEF
      CALL PENTA_FIT_RAD_TRANS

      ! Technically speaking proc_string contains the unique name of this equilibrium 
      WRITE(temp_str,'(A,A,i4.4)') TRIM(proc_string),'_k',ii
      CALL PENTA_OPEN_OUTPUT(TRIM(temp_str))

      ! Only search new ambipolar Er solution if look_for_ambipolar = true
      ! If not, take the current THRIFT_ER solution
      ! The boolean look_for_ambipolar is controled in thrift_plasma_solver_mod
      IF(solve_plasma_equations .AND. .NOT. look_for_ambipolar) THEN
      !  Need to think about if and how to implement this in STELLOPT
      !      ! Get Er value from THRIFT
      !      ij =  MINLOC(ABS(SQRT(THRIFT_S) - rho_local), dim=1)
      !      ! Need to define num_roots and set Er_roots
      !      num_roots = 1
      !      Er_roots(1) = THRIFT_ER(ij,mytimestep)
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
      DO ij=1,num_roots
            IF(root_type(ij)) THEN
                  JBS_PENTA(ik) = J_BS_ambi(ij)
                  !etapar_PENTA(ik) = 1.0_rprec / sigma_par_ambi(ij)
                  Er_PENTA(ik) = Er_roots(ij)
                  ! NEO fluxes
                  !GNEO_PENTA(:,k) = Gammas_ambi(:,i)
                  !QNEO_PENTA(:,k) = QoTs_ambi(:,i) * Temps * e_charge
                  ! NEO particle transport coefficients
                  !Dn_PENTA(:,k) = (/ (-L_n_ambi(i,j,j), j=1,nion_prof+1) /)
                  !cn_PENTA(:,k) = MATMUL(L_T_ambi(i,:,:),e_charge*dTdrs) / dens + Er_PENTA(k)*SUM(L_Er_ambi(i,:,:),dim=2) / dens
                  ! NEO heat transport coefficients
                  !Dp_PENTA(:,k) = (/ (-R_T_ambi(i,j,j), j=1,nion_prof+1) /) * e_charge*Temps / dens 
                  !cp_PENTA(:,k) = MATMUL(R_n_ambi(i,:,:),dndrs)/dens - MATMUL(R_T_ambi(i,:,:),(e_charge*Temps/dens)*dndrs)/dens &
                  !                  + Er_PENTA(k)*SUM(R_Er_ambi(i,:,:),dim=2) / dens
                  EXIT
            ENDIF
      END DO
      !       
      CALL PENTA_RUN_5_CLEANUP(lscreen)
   END DO
   ! Now deallocate all the helper arrays
   DEALLOCATE(ik_penta, s_penta, ne_local, te_local, ni_local, ti_local)
   DEALLOCATE(dnedrho_local, dtedrho_local, dnidrho_local, dtidrho_local)
   DEALLOCATE(vp_local, bdotb_local)
   DEALLOCATE(dkes_nustar, dkes_erstar)
   DEALLOCATE(DKES_D11, DKES_D31, DKES_D33)
   ! Now make sure the master thread has the full set of data.
!DEC$ IF DEFINED (MPI_OPT)
   CALL MPI_BARRIER(MPI_COMM_MYWORLD,ierr_mpi)
   IF (myworkid == master) THEN
      CALL MPI_REDUCE(MPI_IN_PLACE, JBS_PENTA, nsurf_penta, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
      !CALL MPI_REDUCE(MPI_IN_PLACE, etapar_PENTA, nsurf_penta, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_REDUCE(MPI_IN_PLACE, Er_PENTA, nsurf_penta, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
   ELSE
      CALL MPI_REDUCE(JBS_PENTA,    JBS_PENTA, nsurf_penta, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
      !CALL MPI_REDUCE(etapar_PENTA, etapar_PENTA, nsurf_penta, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_REDUCE(Er_PENTA,     Er_PENTA, nsurf_penta, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
   END IF
!DEC$ ENDIF
   IF (lscreen) WRITE(6,'(a)') ' -----------------  PENTA CALCULATION (DONE) ----------------'
   RETURN
!----------------------------------------------------------------------
!  END SUBROUTINE
!----------------------------------------------------------------------
   END SUBROUTINE stellopt_penta
