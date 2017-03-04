!-----------------------------------------------------------------------
!     Subroutine:    stelltran_init
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          08/21/2015
!     Description:   This routine is for reading in the input files.
!-----------------------------------------------------------------------
      SUBROUTINE stelltran_init
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE stelltran_input_mod
      USE stelltran_runtime
!       USE stelltran_vars, ONLY: De, Xe, Xi, prof_length
      USE vmec_params, ONLY: norm_term_flag, bad_jacobian_flag,&
                             more_iter_flag, jac75_flag, input_error_flag,&
                             phiedge_error_flag, ns_error_flag, &
                             misc_error_flag, successful_term_flag, &
                             restart_flag, readin_flag, timestep_flag, &
                             output_flag, cleanup_flag, reset_jacdt_flag
      USE mpi_params
      USE bootsj_input
      USE stelltran_data
      USE safe_open_mod, ONLY: safe_open
      USE stelltran_vars, ONLY: prof_length, Xe, Xi, De ! Can't find a way around not needing these
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Local Error flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier, iunit,ik, n
      INTEGER ::  ictrl(5)
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! First read the STELLTRAN_INPUT namelist from the input file
      ier = 0
      CALL read_stelltran_input(id_string,ier,myid)
      

      ! Now read and initialize the database file
      SELECT CASE(TRIM(run_type))
         CASE('analysis')
            CALL stelltran_read_database
         CASE('predictive')
            ALLOCATE(Ip(ntimesteps),Vloop(ntimesteps))
            ALLOCATE(te_f(nprof,ntimesteps),ne_f(nprof,ntimesteps),&
                      ti_f(nprof,ntimesteps),zeff_f(nprof,ntimesteps),P_ecrh(3,ntimesteps))
             ALLOCATE(te_s(nprof),ne_s(nprof),ti_s(nprof),zeff_s(nprof))
             FORALL(ik = 1:nprof) ne_s(ik) = (ik-1)/REAL(nprof-1)
             FORALL(ik = 1:nprof) te_s(ik) = (ik-1)/REAL(nprof-1)
             FORALL(ik = 1:nprof) ti_s(ik) = (ik-1)/REAL(nprof-1)
             FORALL(ik = 1:nprof) zeff_s(ik) = (ik-1)/REAL(nprof-1)
             FORALL(ik = 1:nprof) ne_f(ik,1) = 5.0E18 * (1-ne_s(ik)**2)
             FORALL(ik = 1:nprof) te_f(ik,1) = 3.0E2  * (1-te_s(ik)**2)
             ti_f(:,1) = 0.0
             zeff_f(:,1) = 1.0
             P_ecrh(:,1) = (/ 5.0E-01, 7.0E-01, 6.0E-01 /) ! probably want to change this at some point
             Ip = 0.0
             Vloop = 0.0
!            CALL stelltran_read_database_predictive
         CASE DEFAULT
      END SELECT

      ! Now initalize the equilibrium calculations
      SELECT CASE(TRIM(equil_type))
         CASE('vmec2000')
            id_string = id_string(7:LEN(id_string))
            ictrl(1) = restart_flag + readin_flag + reset_jacdt_flag
            ictrl(2) = 0
            ictrl(3) = 50
            ictrl(4) = 0
            ictrl(5) = myid
            CALL runvmec(ictrl,id_string,.false.,0,'')
            ! Need to read BOOTSJ namelist as well
            CALL safe_open(iunit,ier,'input.'//TRIM(id_string),'old','formatted')
            !IF (ier < 0) RETURN
            CALL read_namelist (iunit, ier, 'bootin')
            CLOSE(iunit)
      END SELECT

      ! Allocate Global Helpers
      ALLOCATE(nu_star(prof_length,ntimesteps))
      ALLOCATE(johm_sav(prof_length,ntimesteps))
      ALLOCATE(jboot_sav(prof_length,ntimesteps))
      ALLOCATE(jecrh_sav(prof_length,ntimesteps))
      ALLOCATE(jbeam_sav(prof_length,ntimesteps))

      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE stelltran_init
