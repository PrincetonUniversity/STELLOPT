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
                 mysurf
      REAL(rprec) :: s, rho, iota, phip, chip, btheta, bzeta, bsq, vp, &
                        te, ne, dtedrho, dnedrho
      REAL(rprec), DIMENSION(num_ion_species) :: ni,ti, dtidrho, dnidrho
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: D11, D13, D33
!-----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!-----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (lscreen) WRITE(6,'(a)') ' --------------------  NEOCLASSICAL FLUX USING PENTA  -------------------'
      IF (lvmec) THEN
         ierr_mpi = 0
         ! PENTA is parallelized over radial surfaces in this routine.
         ns_dkes = 0
         DO k = 1, DKES_NS_MAX
            IF ((DKES_K(k) > 0)) ns_dkes = ns_dkes+1
         END DO
         ! Break up work
         CALL MPI_CALC_MYRANGE(MPI_COMM_MYWORLD,1,ns_dkes,mystart,myend)
         DO k = mystart,myend
            ! Calc some information
            mysurf = DKES_K(k)
            s = DBLE(mysurf)/DBLE(ns_eq)
            rho = sqrt(s)
            ier = 0
            CALL EZSpline_interp(iota_spl, rho, iota, ier)
            ier = 0
            CALL EZSpline_interp(phip_spl, rho, phip, ier)
            ! PENTA wants d/ds note that this won't work if rho=0
            phip = 0.5*phip/rho
            chip = iota * phip
            ier = 0
            CALL EZSpline_interp(bu_spl, rho, btheta, ier)
            ier = 0
            CALL EZSpline_interp(bv_spl, rho, bzeta, ier)
            ier = 0
            CALL EZSpline_interp(bsq_spl, rho, bsq, ier)
            ier = 0
            CALL EZSpline_interp(vp_spl, rho, vp, ier)
            ! Vp = dV/dPHI need dVds with VMEC normalization
            vp = vp * phip /(pi2*pi2)
            ! Profiles
            CALL get_prof_te(rho, THRIFT_T(mytimestep), te)
            CALL get_prof_ne(rho, THRIFT_T(mytimestep), ne)
            CALL get_prof_teprime(rho, THRIFT_T(mytimestep), dtedrho)
            CALL get_prof_neprime(rho, THRIFT_T(mytimestep), dnedrho)
            DO j = 1, num_ion_species
               CALL get_prof_ti(rho, THRIFT_T(mytimestep), j, ti(j))
               CALL get_prof_ni(rho, THRIFT_T(mytimestep), j, ni(j))
               CALL get_prof_tiprime(rho, THRIFT_T(mytimestep), j, dtidrho(j))
               CALL get_prof_niprime(rho, THRIFT_T(mytimestep), j, dnidrho(j))
            END DO
            ! DKES Data
            ncstar = COUNT(DKES_NUSTAR < 1E10)
            nestar = COUNT(DKES_ERSTAR < 1E10)

            ! PENTA
            CALL PENTA_SET_ION_PARAMS(nion_prof, DBLE(Zatom_prof), Matom_prof/AMU)
            CALL PENTA_SET_COMMANDLINE(-250.0_rprec,250.0_rprec,DKES_K(k),1,0.0_rprec,1,'','','')
            CALL PENTA_ALLOCATE_SPECIES
            CALL PENTA_SET_EQ_DATA(rho,eq_Aminor,eq_Rmajor,vp,chip,phip,iota,btheta,bzeta,bsq)
            CALL PENTA_SET_PPROF(ne,dnedrho,te,dtedrho,ni,dnidrho,ti,dtidrho)
            ! Check here we have D31 from DKES but need D13 in PENTA
            CALL PENTA_SET_DKES_STAR(ncstar,nestar,DKES_NUSTAR(1:ncstar),DKES_ERSTAR(1:nestar), &
               DKES_D11(mysurf,:,:),DKES_D31(mysurf,:,:),DKES_D33(mysurf,:,:))
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
            CALL PENTA_RUN_4_CLEANUP
         END DO
      ENDIF
      IF (lscreen) WRITE(6,'(a)') ' -------------------  NEOCLASSICAL FLUX CALCULATION DONE  ---------------------'
      RETURN
!-----------------------------------------------------------------------
!     END SUBROUTINE
!-----------------------------------------------------------------------
      END SUBROUTINE thrift_penta