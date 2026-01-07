!-----------------------------------------------------------------------
!     Subroutine:    thrift_equil_p
!     Authors:       L. van Ham
!     Date:          11/XX/2022
!     Description:   This subroutine updated the equilbirium pressure
!-----------------------------------------------------------------------
      SUBROUTINE thrift_equil_p
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE thrift_runtime
      USE thrift_vars
      USE thrift_profiles_mod
      USE vmec_input, ONLY: am_aux_s, am_aux_f, ac_aux_s, ac_aux_f, &
                            ah_aux_s, ah_aux_f, at_aux_s, at_aux_f, &
                            pmass_type, pcurr_type, ph_type, pt_type, &
                            pres_scale
      USE thrift_globals, ONLY: solve_plasma_equations
      USE thrift_plasma_solver_mod, ONLY: evolve_plasma_equations, &
                                          get_fast_alphas_dens
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: i, j
      REAL(rprec) :: s_val, rho_val, p_val
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------


      IF(solve_plasma_equations) CALL evolve_plasma_equations
      ! The subroutine 'evolve_plasma_equations' evolves density and pressure from the previous THRIFT time step until
      ! the current THRIFT time step
      ! The routine MUST update the splines NE3D, NI3D, TI3D, P3D with the data from the new time step
      ! It is done like this because get_prof_ne, get_prof_ti, etc, are called EVERYWHERE throughout the code, so 
      ! I  THINK THE BEST AT THIS POINT IS TO CONSTRUCT THE SPLINES NE3D, NI3D, TE3D, TI3D and P3D at every THRIFT time step
      ! I'm just a bit afraid this might compromise speed of the code, but it is a choice between restructuring the whole code and performence...
      ! 
      ! 
      ! And then the following lines are run, independently if profiles are read from file or evolved:

      IF (lvmec) THEN
         PMASS_TYPE = 'akima_spline'
         PRES_SCALE = one
         DO i = 1, n_eq
           s_val = DBLE(i-1)/DBLE(n_eq-1)
           rho_val = sqrt(s_val)
           CALL get_prof_p(rho_val,THRIFT_T(mytimestep),p_val)
           !WRITE(6,*) rho_val,THRIFT_T(mytimestep),p_val
           !CALL FLUSH(6)
           AM_AUX_S(i) = s_val
           AM_AUX_F(i) = p_val
         END DO
      END IF

      ! Save profiles in the THRIFT_## arrays
      DO i = 1,nsj
            rho_val = SQRT( THRIFT_S(i) )
            CALL get_prof_te(rho_val, THRIFT_T(mytimestep), THRIFT_TEMP(1,i,mytimestep))
            CALL get_prof_ne(rho_val, THRIFT_T(mytimestep), THRIFT_DENS(1,i,mytimestep))
            THRIFT_PRESS(1,i,mytimestep) = THRIFT_DENS(1,i,mytimestep) * THRIFT_TEMP(1,i,mytimestep) * e_charge
            DO j = 1, nion_prof
                  CALL get_prof_ti(rho_val, THRIFT_T(mytimestep), j, THRIFT_TEMP(j+1,i,mytimestep))
                  CALL get_prof_ni(rho_val, THRIFT_T(mytimestep), j, THRIFT_DENS(j+1,i,mytimestep))
                  THRIFT_PRESS(j+1,i,mytimestep) = THRIFT_DENS(j+1,i,mytimestep) * THRIFT_TEMP(j+1,i,mytimestep) * e_charge
            END DO
      END DO
      ! Save N_FAST_ALPHAS
      IF(solve_plasma_equations) THEN 
            CALL get_fast_alphas_dens(SQRT(THRIFT_S),THRIFT_FAST_ALPHAS_DENS(:,mytimestep))
      ELSE 
            THRIFT_FAST_ALPHAS_DENS(:,mytimestep) = 0.0
      END IF

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_equil_p

