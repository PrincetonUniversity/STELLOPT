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
                            pmass_type, pcurr_type, ph_type, pt_type
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: i
      REAL(rprec) :: s_val, rho_val, p_val
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      IF (lvmec) THEN
         PMASS_TYPE = 'akima_spline'
         DO i = 1, 99
           s_val = DBLE(i-1)/DBLE(99-1)
           rho_val = sqrt(s_val)
           CALL get_prof_p(rho_val,THRIFT_T(mytimestep),p_val)
           AM_AUX_S(i) = s_val
           AM_AUX_F(i) = p_val
         END DO
      END IF

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_equil_p
