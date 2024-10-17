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
      USE thrift_input_mod
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
      INTEGER :: ns_dkes, k
      REAL(rprec) :: s
!-----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!-----------------------------------------------------------------------
      IMPLICIT NONE
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
            ! This part mimics penta_run_1 but without reading the namelists
            Er_min = -250.0
            Er_max =  250.0
            js     =  DKES_K(k)
            i_append = 0
            B_Eprl = 0.0
            Smax   = 1
            coeff_ext = ''
            run_ident = ''
            pprof_char = ''
            CALL PENTA_ALLOCATE_SPECIES
            CALL PENTA_READ_INPUT_FILES
            CALL PENTA_SCREEN_INFO
            CALL PENTA_ALLOCATE_DKESCOEFF
            CALL PENTA_FIT_DXX_COEF
            CALL PENTA_OPEN_OUTPUT
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