      SUBROUTINE add_resistive_E
      USE stel_kinds
      USE evolution, ONLY: nprecon, fsq_total, xc, l_update_state
      USE perturbation, ONLY: ftol, eta_factor
      USE descriptor_mod, ONLY: iam
      USE siesta_bfield      
#if defined(SKS)
      USE nscalingtools, ONLY: PARFUNCTISL, startglobrow, endglobrow
      USE siesta_state, ONLY: update_state_par, update_state
      USE siesta_displacement, ONLY: update_upperv_par, update_upperv
#else
      USE siesta_state, ONLY: update_state
      USE siesta_displacement, ONLY: update_upperv
#endif
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!      REAL(dp), PARAMETER  :: zero = 0
      INTEGER, PARAMETER   :: nits = 10
      REAL(dp)             :: save_eta
      INTEGER              :: itime
      LOGICAL, PARAMETER   :: lcurrent_only=.TRUE.
!-----------------------------------------------
      IF (fsq_total.LE.ftol) RETURN
!
!     Update resistive contribution to B with new currents
!       
      IF (iam .EQ. 0) THEN
      WRITE (6,'(/,a,i3)') ' UPDATING RESISTIVE E-FIELD: ITERATIONS=',nits
      ENDIF

!     No perturbation, vsupX = 0
      xc = 0
#if defined (SKS)
      IF (PARFUNCTISL) THEN
      CALL init_state_par(.FALSE.)
      CALL update_upperv_par

      save_eta = eta_factor
      eta_factor = eta_factor/MAX(1,nits)

!     Diffuse B-field (eta*j) from last update_state call
      DO itime=1,nits

         CALL init_state_par(lcurrent_only)                             !Use new B-field to get KsubXij currents
         CALL update_bfield_par(.TRUE.)                                 !Resistive update to B-field
         CALL update_state_par(.FALSE., zero, zero)                     !Updates B-field with resistive piece

      END DO
      ELSE
#endif
      CALL init_state(.FALSE.)
      CALL update_upperv

      save_eta = eta_factor
      eta_factor = eta_factor/MAX(1,nits)

!     Diffuse B-field (eta*j) from last update_state call
      DO itime=1,nits

         CALL init_state(lcurrent_only)                                 !Use new B-field to get KsubXij currents
         CALL update_bfield(.TRUE.)                                     !Resistive update to B-field
         CALL update_state(.FALSE., zero, zero)                         !Updates B-field with resistive piece

      END DO
#if defined (SKS)
      END IF
#endif

      eta_factor = save_eta

      END SUBROUTINE add_resistive_E
