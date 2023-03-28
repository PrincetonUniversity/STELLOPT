!-----------------------------------------------------------------------
!     Subroutine:    thrift_run_bootstrap
!     Authors:       L. van Ham
!     Date:          11/XX/2022
!     Description:   This subroutine calls the relevant bootstrap
!                    current code.
!-----------------------------------------------------------------------
      SUBROUTINE thrift_run_bootstrap
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE thrift_runtime
      USE thrift_vars
      USE thrift_equil
      USE thrift_profiles_mod
      USE thrift_funcs
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: i, ier
      REAL(rprec) :: epsilon, pp1, pm1, ds
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: pprime, j_temp
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      ier = 0

      ! Check to make sure we're not zero beta
      IF (eq_beta == 0) THEN
         THRIFT_JBOOT(:,mytimestep) = 0
         RETURN
      END IF

      SELECT CASE(TRIM(bootstrap_type))
         CASE ('model','simple','test')
            ! Calculate dpds 
            ALLOCATE(pprime(nsj))
            ds = THRIFT_S(2)-THRIFT_S(1)
            pprime(1) = 0
            DO i = 2, nsj-1
               CALL get_prof_p( SQRT(THRIFT_S(i+1)), THRIFT_T(mytimestep), pp1)
               CALL get_prof_p( SQRT(THRIFT_S(i-1)), THRIFT_T(mytimestep), pm1)
               pprime(i) = (pp1-pm1)/(2*ds)
            END DO
            pprime(nsj) = 2*pprime(nsj-1) - pprime(nsj-2)

            ! j_BS = sqrt(epsilon) Rmajor *dp/dPhi
            ! epsilon = a/R (inverse aspect ratio)
            ! dp/dPhi : Pa/Wb (toroidal flux derivative)
            ! dp/dPhi = dp/ds * ds/dPhi
            !         = dp/ds / Phi_edge
            ! j_BS = s*sqrt(epsilon)*Rmajor/Phi_edge*dp/ds
            
            ALLOCATE(j_temp(nsj))
            epsilon = eq_Aminor/eq_Rmajor ! Inverse Aspect Ratio
            j_temp = SQRT(epsilon)*eq_Rmajor/eq_phiedge*THRIFT_S*pprime
            CALL Js_to_Jrho(j_temp, THRIFT_JBOOT(:,mytimestep))
            DEALLOCATE(pprime, j_temp)

         CASE ('bootsj')
            CALL thrift_paraexe('booz_xform',proc_string,lscreen_subcodes)
            CALL thrift_paraexe('bootsj',proc_string,lscreen_subcodes)
         CASE ('sfincs')
      END SELECT

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_run_bootstrap

